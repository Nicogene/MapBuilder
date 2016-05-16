#include "MapBuilder.h"

//#include <yarp/os/all.h>
#include <yarp/sig/all.h>
#include <yarp/sig/Image.h>
#include <yarp/dev/PolyDriver.h>
#include <yarp/dev/IEncoders.h>
#include <yarp/os/Time.h>
#include <yarp/sig/Vector.h>
#include <yarp/math/Math.h>

#include <iCub/iKin/iKinFwd.h>
#include <iCub/iKin/iKinIpOpt.h>

#include <opencv2/core/eigen.hpp>


#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <fstream>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <iostream>
#include <stdint.h>
#include <iomanip>

#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <cvsba/cvsba.h>

#include "bal_problem.h"
#include <ceres/ceres.h>
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "my_reprojection_error.h"

#include <unordered_set>

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
//#include "g2o/math_groups/se3quat.h"
#include "g2o/solvers/structure_only/structure_only_solver.h"

using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;
using namespace iCub::iKin;





//MapBuilder Constructor

MapBuilder::MapBuilder(int numC, bool WithYDplayer, int mBA ){
    this->NumOfCouples=numC;
    this->WithYarpDataPlayer=WithYDplayer;
    this->methodBA=mBA;
    this->first=true;

}

MapBuilder::MapBuilder(){
    this->NumOfCouples=1;
    this->WithYarpDataPlayer=false;
    first=true;
    this->methodBA=-1;

}
//MapBuilder Destructor

MapBuilder::~MapBuilder(){

    this->images.clear();//deallocate memory
    this->buckets.clear();
    this->points.clear();
    this->visibility.clear();
    this->Rotations.clear();
    this->Translations.clear();
    this->ProjectionMatrices.clear();
    this->PointsXYZ.clear();
    std::cout<<"~MapBuilder"<<std::endl;
}
void MapBuilder:: allocateMem(){
    this->visibility.resize(this->NumOfCouples*2);
    this->Translations.resize(this->NumOfCouples*2);
    this->Rotations.resize(this->NumOfCouples*2);
    this->images.resize(this->NumOfCouples*2);//allocate memory
    this->points.resize(this->NumOfCouples*2);
    this->ProjectionMatrices.resize(this->NumOfCouples*2);
}
void MapBuilder::setNumOfCouples(int numC){
    this->NumOfCouples=numC;
}

int MapBuilder::getNumOfCouples(){ return this->NumOfCouples;}

void MapBuilder::setWithYarpDataPlayer(bool WYDP){this->WithYarpDataPlayer=WYDP;}

bool MapBuilder::getWithYarpDataPlayer(){return this->WithYarpDataPlayer;}

void MapBuilder::setCameraMatrices(double fxL, double fyL, double cxL, double cyL,
                                   double fxR, double fyR, double cxR, double cyR){
    this->CameraMatrixL = (cv::Mat_<double>(3,3) << fxL, 0, cxL, 0, fyL, cyL, 0, 0, 1);
    this->CameraMatrixR = (cv::Mat_<double>(3,3) << fxR, 0, cxR, 0, fyR, cyR, 0, 0, 1);

}

void MapBuilder::setDistorsionMatrices(double d0L,double d1L, double d2L, double d3L, double d4L,
                                       double d0R,double d1R, double d2R, double d3R, double d4R){

    this->distCoeffsL= (cv::Mat_<double>(5,1) << d0L, d1L, d2L, d3L, d4L);
    this->distCoeffsR= (cv::Mat_<double>(5,1) << d0R, d1R, d2R, d3R, d4R);

}


bool MapBuilder::process(){
    bool res=true;
    this->allocateMem();
    PointsXYZ.clear();
    for(int i=0;i<points.size();i++){
        points[i].clear();
    }
    this->collectImages();
    this->setCameraMatrices(407.386849650585816, 408.446478587298941,
                            335.142950296626907, 247.077012143219406,
                            411.389187422896157, 412.513836073373398,
                            320.913970488683390, 224.679749311050074);

    this->setDistorsionMatrices(-0.370087506134749, 0.122569815417904, -0.000557593423888,
                                -0.000650939992441, 0.000000000000000,
                                -0.381727559529518, 0.138669468725887, -0.000320490995377,
                                -0.000398273646263, 0.000000000000000
                                );
    if(this->findPoints())
    {

        //        this->getRotationsThroughEssential();
        this->getTransformationsToFirstLeft();
        this->checkPointIsVisible(); //after optimisation TODO
        this->removeRowsColsVisibility();
        this->initialize3DPoints();
        if(this->methodBA==CVSBA){
            if(this->cvSba())
                std::cout<<"BA riuscito!"<<std::endl;
            else
                std::cout<<"BA fallito!"<<std::endl;
        }
        else if(this->methodBA==CERESBA)
        {
            this->writeFileCeres();
            this->ceresBa();
        }
        this->visualize3DMap();
        return false;
    }
    else
        res=false;

    return res;
    //    return false;
}

bool MapBuilder::processCouples(){
    bool res=true;
    this->allocateMem();
    this->collectImages();
    this->setCameraMatrices(407.386849650585816, 408.446478587298941,
                            335.142950296626907, 247.077012143219406,
                            411.389187422896157, 412.513836073373398,
                            320.913970488683390, 224.679749311050074);

    this->setDistorsionMatrices(-0.370087506134749, 0.122569815417904, -0.000557593423888,
                                -0.000650939992441, 0.000000000000000,
                                -0.381727559529518, 0.138669468725887, -0.000320490995377,
                                -0.000398273646263, 0.000000000000000
                                );
    if(this->findPointsCouples())
    {
        this->getCouplesTransformationsThroughEssential();
        this->getTransformationsCouples();
        return false;
    }
    else
        res=false;

    return res;
    //    return false;
}

void MapBuilder::collectImages(){
    int i=0;
    bool res=false;
    Network yarp;
    BufferedPort<ImageOf<PixelRgb> > imagePortL,imagePortR;
    int startCountR=0, startCountL=0;
    imagePortL.open("/nico/MapBuilder/image/L");
    imagePortR.open("/nico/MapBuilder/image/R");
    yarp.connect("/icub/cam/left",imagePortL.getName());
    yarp.connect("/icub/cam/right",imagePortR.getName());
    while(i<2*this->NumOfCouples){

        ImageOf<PixelRgb> *yarpImageL=imagePortL.read();
        ImageOf<PixelRgb> *yarpImageR=imagePortR.read();
        Stamp sL,sR;
        if(imagePortL.getEnvelope(sL) && imagePortR.getEnvelope(sR)){
            if(first){
                start=sL.getTime();
                first=false;
                std::cout<<"Reinitialize"<<std::endl;
            }
            if(i==0){
                startCountL=sL.getCount();
                startCountR=sR.getCount();
            }
            if(abs(sL.getCount()-startCountL)>2 || abs(sR.getCount()-startCountR)>2 || fabs((sL.getTime())-(sR.getTime()))>0.03){//0.03 is the half delta t
                std::cout<<sL.getCount()-startCountL<<" !!ASYNCRO!! "<<sR.getCount()-startCountR<<std::endl;
                startCountL=sL.getCount();
                startCountR=sR.getCount();
                continue;
            }
            else {
                std::cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXX SYNCRO XXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<std::endl;
                std::cout<<"L-> seq : "<<sL.getCount()<<"  time:  "<<sL.getTime()-start<<std::endl;
                std::cout<<"R-> seq : "<<sR.getCount()<<"  time:  "<<sR.getTime()-start<<std::endl;

            }
            if (yarpImageL!=NULL && yarpImageR!=NULL){
                Image iL(*yarpImageL);
                IplImage* cvImageL=(IplImage*)iL.getIplImage();
                Image iR(*yarpImageR);
                IplImage* cvImageR=(IplImage*)iR.getIplImage();

                cvCvtColor(cvImageL, cvImageL, CV_RGB2BGR);
                cvCvtColor(cvImageR, cvImageR, CV_RGB2BGR);

                //Convert IplImage to cv::Mat
                cv::Mat LOld = cv::cvarrToMat((IplImage*) cvImageL);
                cv::Mat R= cv::cvarrToMat((IplImage*) cvImageR);
                cv::Mat L;
                //resize(L,L,R.size());
                //Mat L(R.rows,R.cols,R.type());

                if(this->WithYarpDataPlayer)
                    L = LOld(cv::Rect(3,3,640,480)).clone();//Le immagini left dell'acquisizione di Tariq sono 648x488
                else
                    L=LOld.clone();

                cv::medianBlur(R,R,3);
                cv::medianBlur(L,L,3);

                this->getTransformationsToRoot(yarp,i);
                //               std::cout<<std::endl<<std::endl;
                //               std::cout<<ProjectionMatrices[i]<<std::endl;
                //               std::cout<<ProjectionMatrices[i+1]<<std::endl;

                this->images[i]=L.clone();
                this->images[i+1]=R.clone();
                i=i+2;
                std::cout<<"Ho collezionato "<<i<<" immagini"<<std::endl;

            }
            else{
                std::cout<<"Non ho trovato nulla!"<<std::endl;
                continue;
            }
        }
        else{
            continue;
        }


    }
    //    for(int i=0;i<images.size();i++){
    //        std::stringstream imstream;
    ////        cv::imshow("Immagine di riferimento",Lkey);
    //        imstream<<"Image N"<<i;
    //        std::string str=imstream.str();
    //        cv::namedWindow(str,1);
    //        cv::moveWindow(str,100*(i+1),100);
    //        cv::circle(images[i],cv::Point(i*10,i*10),3,cv::Scalar( 0, 255, 0 ),4,8,0);
    //        cv::imshow(str,images[i]);
    //        if(cv::waitKey(30) >= 0) break;

    //    }

}

bool MapBuilder::findPoints(){
    bool res=true;
    cv::Mat l_gray,opt_gray;

    std::vector<cv::KeyPoint> keypointsL;
    cv::cvtColor(this->images[0],l_gray,CV_RGB2GRAY);
    cv::FAST(l_gray,keypointsL,15,false,cv::FastFeatureDetector::TYPE_9_16);
    if(keypointsL.size()>0){
        cv::KeyPoint::convert(keypointsL, /*pointsL*/points[0]);
        for(int i=0;i<visibility.size();i++){
            std::vector<int> vector1(points[0].size(), 1);//ok lo inizializza con una righe di uni.
            visibility[i]=vector1;
        }



        this->buckets=bucketing(points[0],4,4);
        //        std::cout<<CameraMatrixL<<std::endl; // CONTROLLATE SONO OK
        //        std::cout<<CameraMatrixR<<std::endl;
        //        std::cout<<distCoeffsL<<std::endl;
        //        std::cout<<distCoeffsR<<std::endl;

        for(int i=1;i<this->points.size();i++)
        {
            cv::Mat Lkey=images[0].clone();
            std::vector<uchar> status;
            cv::Mat error;
            cv::cvtColor(this->images[i],opt_gray,CV_RGB2GRAY);
            cv::TermCriteria termcrit(cv::TermCriteria::COUNT|cv::TermCriteria::EPS,20,0.03);
            cv::Size winSize(31,31);

            cv::calcOpticalFlowPyrLK(l_gray,opt_gray,points[0],points[i],status,error,winSize);
            //            std::cout<<error.rows<<"x"<<error.cols<<" channels: "<<error.channels()<<" type "<<error.type()<<std::endl;// Nx1 ch 1, type 5
            std::ofstream errorfile;
            errorfile.open("Error.txt");
            for(int j=0;j<error.rows;j++){
                errorfile<<error.at<float>(j,0)<<std::endl;
            }


            //            std::cout<<"error "<<error.rows<<"x"<<error.cols<<" channels "<<error.channels()<<" type "<<error.type()<<std::endl; //Nx1 ch 1 type 5

            if(points[i].size()>0){
                for (int j=0;j<points[i].size();j++){        // non posso piu' usare questo vincolo, due immagini possono anche essere lontane
                    // TODO check thresholds
                    if((std::isnan(error.at<float>(j,0))) || (error.at<float>(j,0))> 1000 || (int)status[j]==0 || (abs(points[0][j].x-points[i][j].x)<10 && abs(points[0][j].y-points[i][j].y)<10)){
                        visibility[i][j]=0;
                        //points[i].erase (points[i].begin()+j);//TODO NON SERVE VERO?
                    }
                    else{
                        cv::circle(Lkey,cv::Point(this->points[i][j].x,this->points[i][j].y),1,cv::Scalar( 255, 0, 0 ),1,8,0);
                        cv::line(Lkey, this->points[0][j], this->points[i][j], cv::Scalar( 255, 0, 0 ));
                    }
                }
                if(i%2==0)
                    cv::undistortPoints(points[i],points[i],this->CameraMatrixL,this->distCoeffsL);
                else
                    cv::undistortPoints(points[i],points[i],this->CameraMatrixR,this->distCoeffsR);

                //                std::cout<<"FIND POINTS  NUM points: "<<points.size()<<"x"<<points[i].size()<<std::endl;
            }
            else return false;

            //            std::stringstream filestream;
            //            std::ofstream myfile;
            //            filestream<<"Points"<<i<<".txt";
            //            myfile.open(filestream.str());
            //            for(int j=0;j<points.at(i).size();j++){
            //                myfile<<points[i][j].x<<" "<<points[i][j].y<<std::endl;
            //            }

//            std::stringstream imstream;
//            cv::imshow("Immagine di riferimento",Lkey);
//            imstream<<"Image N"<<i;
//            std::string str=imstream.str();
//            cv::namedWindow(str,1);
//            cv::moveWindow(str,200*(i+1),100);
//            cv::imshow(str,images[i]);
//            if(cv::waitKey(30) >= 0) return false;

        }
        //        std::vector<cv::Point2f> tmp;
        cv::undistortPoints(this->points[0],this->points[0],this->CameraMatrixL,this->distCoeffsL);
        //        points[0]=tmp;

    }
    else
        return false;

    //    std::cout<<"SIZE "<<visibility.size()<<std::endl;
    //    std::cout<<"Keypoints "<<keypointsL.size()<<std::endl;
    // //OK
    //    for(int i=0; i< visibility.size();i++){
    //        std::cout<<"Riga "<<i<<visibility[i].size()<<std::endl;
    //    }



    return res;

}

// TODO Remove points on the border of the size of the optical flow window

void MapBuilder::removeRowsColsVisibility(){
    //remove rows
    for(int i=0;i<this->visibility.size();i++){
        int sum=0;
        for (int n : visibility[i])
            sum += n;
        if(sum<10){
            ProjectionMatrices.erase (ProjectionMatrices.begin()+i);
            Translations.erase (Translations.begin()+i);
            Rotations.erase (Rotations.begin()+i);
            images.erase (images.begin()+i);
            visibility.erase (visibility.begin()+i);
            points.erase (points.begin()+i);
            //            std::cout<<"Ciao T "<<sum<<std::endl;
            std::cout<<"Frame erased, too few correspondences"<<std::endl;
            i--;
        }
        //        else
        //            std::cout<<"Ciao F "<<sum<<std::endl;
    }
    //remove cols
    for(int i=0;i<visibility[0].size();i++){
        int sum=0;
        for(int j=0;j<visibility.size();j++){
            sum+=visibility[j][i];
        }
        if(sum<2){
            for(int j=0;j<visibility.size();j++){
                visibility[j].erase (visibility[j].begin()+i); // TODO check that points are deleted in the same order
            }
//            std::cout<<"PRIMA "<<points[0].size()<<std::endl; //OK TESTED, FA quello che serve
            points[0].erase(points[0].begin()+i);
//            std::cout<<"DOPO "<<points[0].size()<<std::endl;
            i--;
            //std::cout<<"Column of visibility erased no correspondence for the point "<<i<<std::endl;
        }
    }
//    for(int j=0;j<visibility.size();j++) //OK TESTED
//        std::cout<<"PROVA CHE SIANO UGUALI"<<visibility[j].size()<<std::endl;

    std::ofstream myfile;
    myfile.open("Visibility.txt");
    for(int i=0;i<visibility.size();i++){
        for(int j=0;j<visibility[i].size();j++){
            myfile<<visibility[i][j]<<" ";
        }
        myfile<<std::endl;

    }


}

bool MapBuilder::findPointsCouples(){ //DEPRECIATED
    bool res=true;

    for(int i=0;i<this->points.size();i=i+2)
    {
        cv::Mat l_gray,opt_gray;
        std::vector<cv::Point2f> pointsL, pointsR;


        std::vector<cv::KeyPoint> keypointsL;
        cv::cvtColor(this->images[i],l_gray,CV_RGB2GRAY);
        FAST(l_gray,keypointsL,15,false,cv::FastFeatureDetector::TYPE_9_16);
        if(keypointsL.size()>0){

            cv::KeyPoint::convert(keypointsL, pointsL);
            //                std::ofstream controllopunti;
            //                controllopunti.open("ControlloPuntiOrigine.txt");
            //                for(int i=0;i<pointsL.size();i++){
            //                    controllopunti<<pointsL[i].x<<" "<<pointsL[i].y<<std::endl;
            //                }

            std::vector<uchar> status;
            cv::Mat error;
            cv::cvtColor(this->images[i+1],opt_gray,CV_RGB2GRAY);
            this->buckets=bucketing(pointsL,4,4);
            //            std::cout<<"SIZE BUCKETS"<<buckets.size()<<std::endl; OK!
            //            std::cout<<"POINTS "<<pointsL.size()<<std::endl;
            //            for(int k=0;k<buckets.size();k++){
            //                std::cout<<"SIZE "<<k<<" "<<buckets[k].size()<<std::endl;
            //            }
            calcOpticalFlowPyrLK(l_gray,opt_gray,pointsL,pointsR,status,error,cv::Size (15,15));
            //            controllopunti<<"XXXXXXXXXXXXXXXXXXXXXXXXXX R XXXXXXXXXXXXXXXXXXX"<<std::endl<<std::endl;
            //            for(int i=0;i<pointsR.size();i++){

            //                controllopunti<<pointsR[i].x<<" "<<pointsR[i].y<<std::endl;
            //            }

            if(pointsL.size()>0){
                for (int j=0;j<pointsL.size();j++){        // non posso piu' usare questo vincolo, due immagini possono anche essere lontane
                    if((std::isnan(error.at<float>(j,1))) || (norm(pointsL[j]-pointsR[j])>50) || (int)status[j]==0 || (abs(pointsL[j].x-pointsR[j].x)< 5 && abs(pointsL[j].y-pointsR[j].y)< 5 )){
                        continue;
                    }
                    else{
                        points[i].push_back(pointsL[j]);
                        points[i+1].push_back(pointsR[j]);
                    }

                }
                cv::undistortPoints(points[i],points[i],this->CameraMatrixL,this->distCoeffsL);
                cv::undistortPoints(points[i+1],points[i+1],this->CameraMatrixR,this->distCoeffsR);

                //                    std::ofstream controllopunti;
                //                    controllopunti.open("ControlloPuntiUndistort.txt");
                //                    for(int j=0;j<points[i].size();j++){
                //                        controllopunti<<points[i][j].x<<" "<<points[i][j].y<<std::endl;
                //                    }

                //                std::cout<<"FIND POINTS  NUM points: "<<points.size()<<"x"<<points[i].size()<<std::endl;
            }
            else return false;



        }




        else
            return false;
    }
    return res;

}

void MapBuilder::getCouplesTransformationsThroughEssential(){
    std::ofstream rotationFile;
    rotationFile.open("RodrAngles.txt");
    for(int i=0;i<points.size();i=i+2){
        cv::Mat Rot,t,E,rotvec,mask;
        double focal = 1.0;
        cv::Point2d pp(0.0, 0.0);
        //        cv::Mat F=cv::findFundamentalMat(points[i],points[i+1]); FUNZIONA PEGGIO SE FACCIAMO IL CORRECTMATCHES

        //        cv::correctMatches(F,points[i],points[i+1],points[i],points[i+1]);// correct the matches with the geometric error constraint


        //        E = cv::findEssentialMat(points[i], points[i+1], focal, pp, cv::RANSAC, 0.999, 1.0, mask); // mask: Nx1(ch 1)
        //        cv::recoverPose(E, points[i], points[i+1], Rot, t, focal, pp, mask);

        cv::Mat RotFiltered=featureSelection(100,points[i],points[i+1]);
        if(RotFiltered.empty()){
            return;
        }
        cv::Rodrigues(RotFiltered,rotvec);
        rotationFile<<rotvec.at<double>(0,0)*CTRL_RAD2DEG<<" "<<rotvec.at<double>(0,1)*CTRL_RAD2DEG<<" "<<rotvec.at<double>(0,2)*CTRL_RAD2DEG<<std::endl;
        //        std::cout<<rotvec<<std::endl; // Questo Funziona, l'altro accesso at<double> no
    }
}

void MapBuilder::getRotationsThroughEssential(){
    std::ofstream maskfile;
    maskfile.open("Masks.txt");

    for(int i=1;i<visibility.size();i++){
        std::vector<cv::Point2f> pointsL, pointsR;
        for(int j=0;j<visibility[0].size();j++)
        {
            if(visibility[i][j]==1 /*&& points[0][j].x>0 && points[0][j].y>0 && points[i][j].x>0 && points[i][j].y>0 */){
                pointsL.push_back(points[0][j]);
                pointsR.push_back(points[i][j]);
            }
        }
        if(pointsL.size()<8 || pointsR.size()<8)//8 per il feature selection
            continue;

        double focal = 1.0;
        cv::Point2d pp(0.0, 0.0);
        cv::Mat Rot,t,mask,E;

        //        E = cv::findEssentialMat(pointsL, pointsR, focal, pp, cv::RANSAC, 0.999, 0.001, mask); // mask: Nx1(ch 1)
        //        cv::recoverPose(E, pointsL, pointsR, Rot, t, focal, pp, mask);
        //        std::cout<<" PRIMA "<<ProjectionMatrices[i]<<std::endl; //OK
        Rot=featureSelection(25,pointsL,pointsR);
        //        ProjectionMatrices[i](cv::Rect(0,0,3,3))=Rot;
        Rot.copyTo(ProjectionMatrices[i](cv::Rect(0,0,3,3)));
        //        std::cout<<Rot<<std::endl;
        //        std::cout<<"DOPO "<<ProjectionMatrices[i]<<std::endl; //OK
        maskfile<<mask<<std::endl<<std::endl;
        //        std::cout<<"mask "<<mask.rows<<"x"<<mask.cols<<" channel "<<mask.channels()<<" type "<<mask.type()<<std::endl; Nx1 uchar
        for(int k=0;k<mask.rows;k++){
            if((int)mask.at<uchar>(k,1)==0){

                visibility[i][k]=0;

            }

        }


    }
}

std::vector<std::vector<int>> MapBuilder::bucketing(std::vector<cv::Point2f> pts, int width, int height)
{
                              std::vector<std::vector<int>> buckets;
                              buckets.resize(width*height);
                              int r,c;
                              for(int j=0; j<pts.size();j++){
    r=floor(pts[j].y/(480/height));
    c=floor(pts[j].x/(640/width));
    //            std::cout<<r<<" and "<<c<<std::endl;
    buckets[r*width+c].push_back(j);
}
return buckets; }

cv::Mat MapBuilder::myTriangulate(cv::Mat Rot,cv::Mat t, cv::Mat x1, cv::Mat x2 ){
    cv::Mat xf(3,x1.cols,CV_64FC1);//Verificare che effettivamente abbia N righe, possono essere anche cols.
    x1.convertTo(x1,CV_64FC1);
    x2.convertTo(x2,CV_64FC1);
    Rot.convertTo(Rot,CV_64FC1);
    t.convertTo(t,CV_64FC1);



    //calculate lines
    cv::Mat v1,v2;
    v1=x1;
    v2= Rot * x2;


    //distance formulas based on Schneider pp 409-412
    cv::Mat a,b,c,d,e,s,denom,num;
    cv::reduce(v1.mul(v1),a,0,CV_REDUCE_SUM);
    cv::reduce(v1.mul(v2),b,0,CV_REDUCE_SUM);
    cv::reduce(v2.mul(v2),c,0,CV_REDUCE_SUM);
    d = t.t() * v1;
    e = t.t() * v2;

    //Non ho implementato questa linea, a che serve?
    //denom(denom < eps) = 1; % accounts for parallel lines

    denom = a.mul(c) - b.mul(b);
    num = c.mul(d) - b.mul(e);
    cv::divide(num,denom,s);


    // compute lines length
    cv:: Mat r1,sqrta;
    cv::sqrt(a,sqrta);
    r1=sqrta.mul(s);

    //TODO mettere a 0 gli elementi della visibility con prof negativa.

    //3d points in P1 coordinates
    cv::Mat rim,z,v1tworows,v1pow,v1sum;
    v1tworows=v1(cv::Rect(0,0,x1.cols,2));
    cv::pow(v1tworows,2.0,v1pow);
    cv::reduce(v1pow,v1sum,0,CV_REDUCE_SUM);
    cv::sqrt(v1sum+1,rim);


    cv::divide(r1,rim,z);

    xf(cv::Rect(0,0,x1.cols,1)) = x1(cv::Rect(0,0,x1.cols,1)).mul(z);
    xf(cv::Rect(0,1,x1.cols,1)) = x1(cv::Rect(0,1,x1.cols,1)).mul(z);
    xf(cv::Rect(0,2,x1.cols,1)) = x1(cv::Rect(0,2,x1.cols,1)).mul(z);


    return xf;
}

cv::Mat MapBuilder::featureSelection(int numIter, std::vector<cv::Point2f> pointsL, std::vector<cv::Point2f> pointsR){
    if(numIter<=0)
        numIter=1;
    srand(time(0));
    for(int j=0; j<buckets.size();j++){
        if(buckets[j].empty()){
            cv::Mat m;
            std::cout<<"Pochi punti, c'e' un bucket vuoto!"<<std::endl;
            return m;
        }

    }
    std::vector<cv::Mat> Rots(numIter);
    std::vector<double> rpjerror(numIter);
    double focal = 1.0;
    cv::Point2d pp(0.0, 0.0);
    //    std::cout<<"Size dei points "<<pointsL.size()<<std::endl;
    std::vector<int> indecesIn(buckets.size()/2);
    for(int i=0;i<numIter;i++){
        std::vector<cv::Point2f> bucketsPointsL,bucketsPointsR;
        bucketsPointsL.resize(buckets.size()/2);bucketsPointsR.resize(buckets.size()/2);
        int it=0, start;
        if(i%2==0){
            start=0;
        }
        else{
            start=1;
        }
        for(int k=start;k<buckets.size();k=k+2,it++){
            int index= rand() % (*std::max_element(buckets[k].begin(),buckets[k].end())- *std::min_element(buckets[k].begin(),buckets[k].end()) + 1) + *std::min_element(buckets[k].begin(),buckets[k].end());
            //            std::cout<<"Index "<<index<<" "<<buckets[k].size()<<std::endl; //Inedex e size ok
            indecesIn[it]=index;
        }
        //        std::random_shuffle(indecesIn.begin(),indecesIn.end()); // l'ho commentato perche' gli do in pasto tutti gli 8 punti.
        for(int j=0;j<bucketsPointsL.size();j++){
            bucketsPointsL[j]=pointsL[indecesIn[j]];
            bucketsPointsR[j]=pointsR[indecesIn[j]];
            //            std::cout<<bucketsPointsL[j].x<<" "<<bucketsPointsL[j].y<<std::endl;
            //            std::cout<<bucketsPointsR[j].x<<" "<<bucketsPointsR[j].y<<std::endl;//Punti Ok
        }
        cv::Mat Proj2,out_h,out,rotvec,distcoeff,t, eye3x3, E , mask;
        eye3x3=eye3x3.eye(3,3,CV_32F);
        std::vector<cv::Point2f> pointreprojR;
        //        bucketsPointsL.push_back(points[0][10]);
        //        bucketsPointsR.push_back(points[1][10]);
        E = cv::findEssentialMat(bucketsPointsL, bucketsPointsR, focal, pp, cv::RANSAC, 0.999, 1.0, mask); // mask: Nx1(ch 1)
        //        std::cout<<bucketsPointsL.size()<<" "<<bucketsPointsR.size()<<std::endl; 5 e 5 ->ok!
        //        std::cout<<E.rows<<"x"<<E.cols<<" channels "<<E.channels()<<std::endl;
        //        std::cout<<E<<std::endl;
        cv::recoverPose(E, bucketsPointsL, bucketsPointsR, Rots[i], t, focal, pp, mask);
        cv::hconcat(Rots[i],t,Proj2);

        cv::Mat Proj0 = (cv::Mat_<float>(3,4) << 1, 0, 0, 0, // The 3D coordinates are referred to the left image system of reference
                         0, 1, 0, 0,
                         0, 0, 1, 0);
        cv::triangulatePoints(Proj0,Proj2,bucketsPointsL,bucketsPointsR,out_h);//Le dimensioni di TUTTE le matrici sono corrette.
        cv::convertPointsFromHomogeneous(cv::Mat(out_h.t()).reshape(4, 1),out);
        cv::Rodrigues(Rots[i],rotvec);
        cv::projectPoints(out,rotvec,t,eye3x3,distcoeff,pointreprojR);

        std::vector<double> ReprErrorvec(bucketsPointsL.size());
        for (int s=0;s<bucketsPointsL.size();s++){
            double ReprError=norm(pointreprojR[s]-bucketsPointsR[s]);
            ReprErrorvec[s]=ReprError;
        }
        cv::Scalar meanReprErr,stdErr;
        meanStdDev(ReprErrorvec,meanReprErr,stdErr);
        rpjerror[i]= (double)meanReprErr.val[0];
    }
    int best;
    best=std::fabs(std::distance(rpjerror.begin(),std::min_element(rpjerror.begin(),rpjerror.end())));
    //    std::cout<<"index "<<best<<" value "<<rpjerror[best]<<std::endl; //Va bene
    return Rots[best];

}

bool MapBuilder::getTransformationsToRoot(Network yarp,int i){
    bool res=true;

    double anglesHead[6], anglesTorso[3];
    Bottle  *headBottle, *torsoBottle;
    BufferedPort<Bottle> encHeadPort, encTorsoPort;
    Property optionsHead, optionsTorso;
    yarp::dev::PolyDriver ddHead(optionsHead), ddTorso(optionsTorso);
    yarp::dev::IEncoders *encHead, *encTorso;


    if(!this->WithYarpDataPlayer){



        //Head
        optionsHead.put("robot", "icub");
        optionsHead.put("device", "remote_controlboard");
        Value& robotname = optionsHead.find("robot");
        cv::String s("/");
        s += robotname.asString();
        s += "/head/control";
        optionsHead.put("local", s.c_str());
        s.clear();
        s += "/";
        s += robotname.asString();
        s += "/head";
        optionsHead.put("remote", s.c_str());

        s.clear();

        //Torso

        optionsTorso.put("robot", "icub");
        optionsTorso.put("device", "remote_controlboard");
        //        Value& robotname = optionsHead.find("robot");
        //        cv::String s("/");
        s="/";
        s += robotname.asString();
        s += "/torso/control";
        optionsHead.put("local", s.c_str());
        s.clear();
        s += "/";
        s += robotname.asString();
        s += "/torso";
        optionsHead.put("remote", s.c_str());


        if (!ddHead.isValid() || !ddTorso.isValid()) {
            std::printf("Device not available.  Here are the known devices:\n");
            std::printf("%s", yarp::dev::Drivers::factory().toString().c_str());
            Network::fini();
            return false;
        }
        else{
            ddHead.view(encHead);
            ddTorso.view(encTorso);
            encHead->getEncoders(anglesHead);
            encTorso->getEncoders(anglesTorso);
        }
    }
    else
    {

        //Head
        encHeadPort.open("/nico/MapBuilder/head/state:i");
        yarp.connect("/icub/head/state:o",encHeadPort.getName());
        headBottle=encHeadPort.read();
        anglesHead[0]=headBottle->get(0).asDouble();
        anglesHead[1]=headBottle->get(1).asDouble();
        anglesHead[2]=headBottle->get(2).asDouble();
        anglesHead[3]=headBottle->get(3).asDouble();
        anglesHead[4]=headBottle->get(4).asDouble();
        anglesHead[5]=headBottle->get(5).asDouble();

        //Torso
        encTorsoPort.open("/nico/MapBuilder/torso/state:i");
        yarp.connect("/icub/torso/state:o",encTorsoPort.getName());
        torsoBottle=encTorsoPort.read();
        anglesTorso[0]=torsoBottle->get(0).asDouble();
        anglesTorso[1]=torsoBottle->get(1).asDouble();
        anglesTorso[2]=torsoBottle->get(2).asDouble();
    }

    //    std::cout<<"head: "<<anglesHead[0]<<" "<<anglesHead[1]<<" "<<anglesHead[2]<<" "<<anglesHead[3]<<" "<<anglesHead[4]<<" "<<anglesHead[5]<<" "<<std::endl;
    //    std::cout<<"Torso: "<<anglesTorso[0]<<" "<<anglesTorso[1]<<" "<<anglesTorso[2]<<std::endl;

    //    double tilt=(anglesHead[3]/180.0)*M_PI;//Ce li da in gradi e mi servono in radianti per il coseno e seno
    //    double version=(anglesHead[4]/180.0)*M_PI;
    //    double vergence=(anglesHead[5]/180.0)*M_PI;
    //    std::cout<<"Tilt: "<<tilt<<" "<<"version: "<<version<<" "<<"vergence: "<<vergence<<std::endl;
    //    std::cout<<"Torso angles: "<<(anglesTorso[0]/180)*M_PI<<" "<<(anglesTorso[1]/180)*M_PI<<" "<<(anglesTorso[2]/180)*M_PI<<" "<<std::endl;

    //    double LangAbs=(version + vergence/2);
    //    double RangAbs=(version - vergence/2);
    //    double RangRel=0;
    ////    if(LangAbs>0 && RangAbs>0)
    ////        RangRel=LangAbs-RangAbs;
    ////    else if(LangAbs<0 && RangAbs<0)// tre casi di angolo relativo(vedi foglio quaderno)
    ////        RangRel=RangAbs-LangAbs;
    ////    else if(LangAbs>0 && RangAbs<0)
    //    RangRel=RangAbs-LangAbs;

    //    cv::Mat Proj1enc = (cv::Mat_<float>(3,4) << 1, 0, 0, 0, // The 3D coordinates are referred to the left image system of reference
    //                                        0, 1, 0, 0,
    //                                        0, 0, 1, 0);
    //    cv::Mat Proj2enc = (cv::Mat_<float>(3,4) << cos(RangRel), 0, sin(RangRel), 68, // The 3D coordinates are referred to the left image system of reference
    //                                        0, 1, 0, 0,
    //                                       -sin(RangRel), 0, cos(RangRel), 0);

    //    this->ProjectionMatrices[i]=Proj1enc.clone();
    //    this->ProjectionMatrices[i+1]=Proj2enc.clone();
    iCubEye eyeR("right"), eyeL("left");
    Vector q0,qf,qhat,xf,xhat;

    Matrix trasfR,trasfL;

    iKinChain *chain;

    chain=eyeR.asChain();

    chain->setAllConstraints(false);//Ugo suggest to drop the kinematics contraint because we are not connected to a robot for now

    q0=chain->getAng();
    std::cout << "Unblocking the torso joints... "<<std::endl;
    chain->releaseLink(0);
    chain->releaseLink(1);
    chain->releaseLink(2);

    std::cout << chain->getDOF() << " DOFs available" << std::endl;

    qf.resize(chain->getDOF());
    double version,vergence;

    version=anglesHead[4]*CTRL_DEG2RAD;
    vergence=anglesHead[5]*CTRL_DEG2RAD;

    qf[0]=anglesTorso[2]*CTRL_DEG2RAD;//the torso angles are inverted
    qf[1]=anglesTorso[1]*CTRL_DEG2RAD;
    qf[2]=anglesTorso[0]*CTRL_DEG2RAD;
    qf[3]=anglesHead[0]*CTRL_DEG2RAD;
    qf[4]=anglesHead[1]*CTRL_DEG2RAD;
    qf[5]=anglesHead[2]*CTRL_DEG2RAD;
    qf[6]=anglesHead[3]*CTRL_DEG2RAD;//tilt, the eye have two DOF, tilt and pan
    qf[7]=version - vergence/2;//pan



    //    qf=chain->setAng(qf);

    //    xf=chain->EndEffPose(false);//euler angles

    //    std::cout << "Current eye end-effector pose: " << xf.toString().c_str() << std::endl;
    trasfR=chain->getH(qf);


    //    for(int i=0;i<trasfR.rows();i++){
    //        std::cout<<std::endl;
    //        for(int j=0;j<trasfR.cols();j++)
    //            std::cout<<trasfR(i,j)<<" ";
    //    }
    chain->setAllConstraints(false);//elimino i vincoli perche' non so quali erano
    //quelli del robot
    chain->clear();

    chain=eyeL.asChain();


    chain->releaseLink(0);
    chain->releaseLink(1);
    chain->releaseLink(2);

    qf[7]=version + vergence/2;//pan

    trasfL=chain->getH(qf);

    //    for(int i=0;i<trasfL.rows();i++){
    //        std::cout<<std::endl;                 OK nel passaggio tra Matrix a cv::Mat
    //        for(int j=0;j<trasfL.cols();j++)
    //            std::cout<<trasfL(i,j)<<" ";
    //    }

    //    std::cout<<"DIMENSIONIL:"<<trasfL.rows()<<"x"<<trasfL.cols()<<std::endl;
    //    std::cout<<"DIMENSIONIR:"<<trasfR.rows()<<"x"<<trasfR.cols()<<std::endl;

    //    std::cout<<trasfL.toString()<<std::endl<<trasfR.toString()<<std::endl;
    //    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl; OK nel passaggio tra Matrix a cv::Mat

    cv::Mat ProjLenc = (cv::Mat_<double>(4,4) << trasfL(0,0), trasfL(0,1), trasfL(0,2), trasfL(0,3),
                        trasfL(1,0), trasfL(1,1), trasfL(1,2), trasfL(1,3),
                        trasfL(2,0), trasfL(2,1), trasfL(2,2), trasfL(2,3),
                        trasfL(3,0), trasfL(3,1), trasfL(3,2), trasfL(3,3));

    cv::Mat ProjRenc = (cv::Mat_<double>(4,4) << trasfR(0,0), trasfR(0,1), trasfR(0,2), trasfR(0,3),
                        trasfR(1,0), trasfR(1,1), trasfR(1,2), trasfR(1,3),
                        trasfR(2,0), trasfR(2,1), trasfR(2,2), trasfR(2,3),
                        trasfR(3,0), trasfR(3,1), trasfR(3,2), trasfR(3,3));

    //       std::cout<<ProjLenc<<std::endl;
    //       std::cout<<ProjRenc<<std::endl;
    //       cv::Mat Prova;
    //       cv::subtract(ProjLenc,ProjRenc,Prova);
    //       std::cout<<Prova<<std::endl;
    this->ProjectionMatrices[i]=ProjLenc.clone();
    this->ProjectionMatrices[i+1]=ProjRenc.clone();


    return res;

}

void MapBuilder::getTransformationsCouples(){
    std::ofstream rotationFile;
    rotationFile.open("RodrAnglesKin.txt");

    for(int i=0;i<ProjectionMatrices.size();i=i+2)
    {
        //        std::cout<<ProjectionMatrices[i]<<std::endl<<ProjectionMatrices[i+1]<<std::endl;
        Eigen::Matrix4f m0;
        Eigen::Matrix4f mi;
        cv::cv2eigen(ProjectionMatrices[i],m0);
        cv::cv2eigen(ProjectionMatrices[i+1],mi);
        //        std::cout<<mi<<std::endl<<std::endl<<m0.inverse()<<std::endl<<std::endl; OK inversa e prodotto funzionano
        mi=m0.inverse()*mi;
        //        std::cout<<mi<<std::endl<<std::endl;
        cv::eigen2cv(mi,ProjectionMatrices[i+1]);
        cv::Mat rotvec;
        Eigen::Vector3d roteig;
        cv::Rodrigues(ProjectionMatrices[i+1](cv::Rect(0,0,3,3)).clone(),rotvec);
        cv::cv2eigen(rotvec,roteig);// Faccio cio' perche' se stampo rotvec su std output ho 000, se stampo su file(usando .at<double>) ho dei numeri enormi, non capisco perche'
        rotationFile<<roteig(0,0)<<" "<<roteig(1,0)<<" "<<roteig(2,0)<<std::endl;// Ho controllato e' un 3x1
    }

}

void MapBuilder::getTransformationsToFirstLeft(){
    Eigen::Matrix4d m0;

    if(true){//Trasformations to the first left image

        cv::cv2eigen(ProjectionMatrices[0],m0);
        //        Eigen::Matrix4f rotx,roty;
        //        rotx=Eigen::Matrix4f::Identity(4,4);
        //        roty=Eigen::Matrix4f::Identity(4,4);
        //        rotx(1,1)=-1.0;//rotazione di 180 gradi attorno all'asse x
        //        rotx(2,2)=-1.0;
        //        roty(0,0)=0.0;roty(2,2)=0.0;
        //        roty(0,2)=1.0;roty(2,0)=1.0;
        //        std::cout<<m0<<std::endl<<m0.inverse()<<std::endl;// l'inversa e' ok, l ho controllata con matlab
        for(int i=0;i<ProjectionMatrices.size();i++)
        {
            Eigen::Matrix4d mi;
            cv::cv2eigen(ProjectionMatrices[i],mi); //Metodo B per odometry(vedi quaderno)
            //mi=m0.inverse()*mi*rotx*roty;//riporto le coordinate nel sistema di riferimento CV se vedi TODO.txt, questa trasformazione non serve a niente
            //            mi(0,3)=mi(0,3)*1000;
            //            mi(1,3)=mi(1,3)*1000;
            //            mi(2,3)=mi(2,3)*1000;
            mi=m0.inverse()*mi;//m0.inverse non cambia m0 e in ho controllato con octave la inversa l'ha fatta giusta.
            cv::eigen2cv(mi,ProjectionMatrices[i]);

            ProjectionMatrices[i](cv::Rect(3,0,1,3)).copyTo(Translations[i]);
            ProjectionMatrices[i](cv::Rect(0,0,3,3)).copyTo(Rotations[i]);


        }
        //    std::cout<<"DETERMINANTE "<<m0.determinant()<<std::endl;
//        m0=m0.inverse()*m0; // delete this lines start from 0

//        cv::eigen2cv(m0,ProjectionMatrices[0]);
//        ProjectionMatrices[0](cv::Rect(0,0,3,3)).copyTo(Rotations[0]);
//        Translations[0] = cv::Mat(3,1,CV_64FC1,cv::Scalar::all(0));
    }
    else
    {// Transformations useful only for the odometry, depreciated, not updated to float->double
        std::vector<cv::Mat> transfToRoot;
        transfToRoot.resize(NumOfCouples*2);
        transfToRoot[0]=ProjectionMatrices[0].clone();
        Eigen::Matrix4f prev=Eigen::Matrix4f::Identity(4,4);
        Eigen::Matrix4f next,prevTrasf,nextTrasf,trasf;
        cv::cv2eigen(ProjectionMatrices[0],m0);
        for(int i=1;i<ProjectionMatrices.size();i++){
            transfToRoot[i]=ProjectionMatrices[i].clone();
            cv::cv2eigen(transfToRoot[i-1],prevTrasf);
            cv::cv2eigen(ProjectionMatrices[i],nextTrasf); //Metodo A per odometry(vedi quaderno)
            trasf=prevTrasf.inverse()*nextTrasf;
            next=trasf*prev;
            cv::eigen2cv(next,ProjectionMatrices[i]);
            prev=next;

        }

    }

    std::ofstream file1,file2,file3,file4;
    file1.open("TraslationsL.txt");
    file2.open("TraslationsR.txt");
    file3.open("AnglesL.txt");
    file4.open("AnglesR.txt");
    for(int j=0;j<this->ProjectionMatrices.size();j++){
        std::ofstream fileprova;
        fileprova.open("PROVAACCESSO.txt");
        Eigen::Vector3f v(ProjectionMatrices[j].at<double>(0,3),ProjectionMatrices[j].at<double>(1,3),ProjectionMatrices[j].at<double>(2,3));
        //        std::cout<<ProjectionMatrices[j]<<std::endl;                 //OK ACCESSO .at<float>
        //        fileprova<<ProjectionMatrices[j].at<float>(0,0)<<" "<<ProjectionMatrices[j].at<float>(0,1)<<" "<<ProjectionMatrices[j].at<float>(0,2)<<std::endl;

        cv::Mat Rotvec;

        cv::Rodrigues(Rotations[j].clone(),Rotvec); //rotvec 3x1;

        if(j%2==0){
            file1<<v(0)<<" "<<v(1)<<" "<<v(2)<<std::endl;
            file3<<Rotvec.at<double>(0,0)*CTRL_RAD2DEG<<" "<<Rotvec.at<double>(1,0)*CTRL_RAD2DEG<<" "<<Rotvec.at<double>(2,0)*CTRL_RAD2DEG<<std::endl;
        }
        else{
            file2<<v(0)<<" "<<v(1)<<" "<<v(2)<<std::endl;
            file4<<Rotvec.at<double>(0,0)*CTRL_RAD2DEG<<" "<<Rotvec.at<double>(1,0)*CTRL_RAD2DEG<<" "<<Rotvec.at<double>(2,0)*CTRL_RAD2DEG<<std::endl;
        }


    }
}

std::vector<cv::Point2d> MapBuilder::getPointsR(std::vector<int> *indeces){
    std::vector<cv::Point2d> pointsR;
    for(int i=0;i<visibility[0].size();i++){
        std::vector<float> txs;//vettore che contiene le tx, vado a prendere le maggiori.
        for(int j=0; j<visibility.size();j++){

            cv::Mat t(Translations[j]);
            cv::Mat transformedt(3,1,CV_64FC1);
            transformedt=Rotations[j].t() * t * visibility[j][i];//lo zero nella visibility esclude che quel frame //TESTED
            txs.push_back(std::fabs(transformedt.at<double>(0,0)));//TESTED
        }
        int index;
        std::vector<float>::iterator it;
        it=std::max_element(txs.begin(),txs.end());
        index= std::distance(txs.begin(), it);
        indeces->push_back(index);
        if(index>images.size()-1)//OK TESTED
            std::cout<<"Problema, l'indice e'.."<<index<<std::endl;
        if(visibility[index][i]==0)//OK TESTED
        {
            std::cout<<"Problema, la visibility e' zero per quel punto in quel frame! "<<txs<<std::endl;}
        pointsR.push_back(points[index][i]);

    }

    return pointsR;

}

void MapBuilder::checkPointIsVisible(){
// remove the points that are seen from few cameras
    //std::cout<<ceil(visibility.size()/2)+1<<std::endl;
    for(int i=0; i<points[0].size();i++){
        int sum=0;
        for(int j=0;j<visibility.size();j++){
            sum+=visibility[j][i];
        }
        if(sum>10/*ceil(visibility.size()/2)+1*/){
        }
        else{
            for(int j=0;j<visibility.size();j++){
                visibility[j][i]==0;
            }
        }
    }
}

void MapBuilder::initialize3DPoints(){


    cv::Mat out_h,mL,mR;


    //std::cout<<ProjectionMatrices[0]<<std::endl;

    std::vector<cv::Point2d> pointsR;

//    for(int j=0; j<visibility[0].size();j++){
//        if(visibility[0][j]==1){
//            pointsL.push_back(points[0][j]); //LA PRIMA RIGA DOVREBBE ESSERE TUTTA DI UNI, CHE SENSO HA??
//        }
//        else
//            continue;
//    }

//    std::cout<<"QUI!!  "<<pointsL.size()<<" "<<points[0].size()<<std::endl;

    std::vector<int> indeces;

    pointsR=this->getPointsR(&indeces);

//    cv::Mat F=cv::findFundamentalMat(points[0],pointsR); //I'm ruining the points (maybe):Direi di si.


//    cv::correctMatches(F,points[0],pointsR,points[0],pointsR);// correct the matches with the geometric error constraint

    //    double focal = 1.0;
    //    cv::Mat mask;
    //    cv::Point2d pp(0.0, 0.0);

    //    cv::Mat   E = cv::findEssentialMat(pointsL, pointsR, focal, pp, cv::RANSAC, 0.999, 2.0, mask); // mask: Nx1(ch 1)
    //    std::ofstream filemask;
    //    filemask.open("Mask.txt");
    //    filemask<<mask;
    //    std::vector<cv::Point2d> pointsLfilt, pointsRfilt;
    //    // mask Nx1 1 canale, type 0, unsigned CV_8U
    ////    std::cout<<"MASK "<<mask.rows<<"x"<<mask.cols<<" channel "<<mask.channels()<<" type "<<mask.type()<<std::endl;
    //    for(int j=0;j<pointsL.size();j++){
    //        if((int)mask.at<uchar>(j, 0) ==1){
    //            pointsLfilt.push_back(pointsL[j]);
    //            pointsRfilt.push_back(pointsR[j]);

    //        }
    //        else
    //        {
    //            visibility[index][indeces[j]]=0;
    //        }
    //    }

    //    std::cout<<"Differenza di punti dopo filtraggio L "<<pointsL.size()-pointsLfilt.size()<<std::endl;
    //    std::cout<<"Differenza di punti dopo filtraggio R "<<pointsR.size()-pointsRfilt.size()<<std::endl;

    //    mL=cv::Mat(pointsLfilt).reshape(1,2);
    //    mR=cv::Mat(pointsRfilt).reshape(1,2);
    mL=cv::Mat(points[0]).reshape(1,2);
    mR=cv::Mat(pointsR).reshape(1,2);

    mL.convertTo(mL,CV_64F);
    mR.convertTo(mR,CV_64F);

    cv::Mat o;
    o=o.ones(1,mL.cols,CV_64F);

    cv::vconcat(mL,o,mL);
    cv::vconcat(mR,o,mR);

    //    std::cout<<pointsL.size()<<" "<<points[index].size()<<std::endl; //ok
    //    std::cout<<mL.rows<<"x"<<mL.cols<<" channels "<<mL.channels()<<std::endl; // 2xN 1 canale;

    //    cv::Mat mL(1,pointsL.size(),CV_64FC2); //stesso risultato che con l'altro reshape
    //    mL=cv::Mat(pointsL).reshape(2,1);
    //    cv::Mat mR(1,pointsR.size(),CV_64FC2);
    //    mR=cv::Mat(pointsR).reshape(2,1);


    //cv::Mat Proj0 = ProjectionMatrices[0](cv::Rect(0,0,4,3)).clone();
    //    cv::Mat ProjN = ProjectionMatrices[index](cv::Rect(0,0,4,3)).clone();//Fa quello che deve
    //    //std::cout<<ProjN.rows<<"x"<<ProjN.cols<<std::endl; //ok e' una 3x4

    cv::Mat Proj0 = (cv::Mat_<float>(3,4) << 1, 0, 0, 0, // The 3D coordinates are referred to the left image system of reference
                     0, 1, 0, 0,
                     0, 0, 1, 0);

    //    ProjN.at<float>(0,3)=68.0;
    //    ProjN.at<float>(1,3)=0.0;
    //    ProjN.at<float>(2,3)=0.0;

    //    cv::Mat Rot,t;
    //    Rot=Rotations[index];
    //    t=Translations[index]*1000;
    //    std::cout<<ProjN<<std::endl<<Rot<<std::endl<<t<<std::endl; OK fa quello che deve




    //    cv::triangulatePoints(Proj0,ProjN,mL,mR,out_h);//Le dimensioni di TUTTE le matrici sono corrette.

    //out_h 4xN
    //    std::cout<<out_h.rows<<"x"<<out_h.cols<<" Type "<<out_h.type()<<" channels "<<out_h.channels()<<std::endl;
    cv::Mat out,ProjN;
    //    std::ofstream indecesfile;
    //    indecesfile.open("ProvaIndeces.txt"); //ok sembra sia giusto ho fatto anche un controllo
    //    indecesfile<<indeces; // incrociato con la visibility

    //    std::cout<<"Indeces size "<<indeces.size()<<std::endl; //ok
    //    std::ofstream fileprova;

    //    fileprova.open("ProvaCutmat.txt");
    for(int i=0;i<points[0].size();i++){
        cv::Mat currOut;
        //        fileprova<<mL.at<double>(0,i)<<" "<<mL.at<double>(1,i)<<" "<<mL.at<double>(2,i)<<std::endl<<std::endl;
        //        fileprova<<mL(cv::Rect(i,0,1,3))<<std::endl;// ok ho fatto le prove, sono uguali.
        if(i==0){
            //            std::cout<<Translations[indeces[i]]*1000<<std::endl; //ok giusto

            out=myTriangulate(Rotations[indeces[i]],Translations[indeces[i]]*1000,mL(cv::Rect(i,0,1,3)),mR(cv::Rect(i,0,1,3)));
            //              cv::hconcat(Rotations[indeces[i]],Translations[indeces[i]]*1000,ProjN);
            //              cv::triangulatePoints(Proj0,ProjN,mL(cv::Rect(i,0,1,2)),mR(cv::Rect(i,0,1,2)),out_h);
        }
        else{
            currOut=myTriangulate(Rotations[indeces[i]],Translations[indeces[i]]*1000,mL(cv::Rect(i,0,1,3)),mR(cv::Rect(i,0,1,3)));
            hconcat(out,currOut,out);
            //              cv::hconcat(Rotations[indeces[i]],Translations[indeces[i]]*1000,ProjN);
            //              cv::triangulatePoints(Proj0,ProjN,mL(cv::Rect(i,0,1,2)),mR(cv::Rect(i,0,1,2)),currOut);
            ////              std::cout<<currOut.rows<<"x"<<currOut.cols<<std::endl;//ok
            //              hconcat(out_h,currOut,out_h);
        }
        //        std::cout<<out_h.rows<<"x"<<out_h.cols<<std::endl;//ok
    }
    //    cv::convertPointsFromHomogeneous(cv::Mat(out_h.t()).reshape(4, 1),out);// vuole una matrice multichannel come ingresso e uscita. OK TESTATO
    //    std::cout<<out.rows<<"x"<<out.cols<<std::endl; //ok
    std::ofstream output,out2;
    output.open("prv.txt");
    out2.open("aaa"); //OK come quello di Prova2
    output<<out;
    for(int j=0;j<out.cols;j++){
        out2<<out.at<double>(0,j)<<" "<<out.at<double>(1,j)<<" "<<out.at<double>(2,j)<<std::endl;
    }
    std::cout<<out.rows<<"x"<<out.cols<<" channels "<<out.channels()<<" type: "<<out.type()<<std::endl;

    for(int j=0;j<out.cols;j++){ //out.rows cv::triangulatePoints
        Eigen::Vector3d p(out.at<double>(0,j),out.at<double>(1,j),out.at<double>(2,j)); //invertire j e 2,0,1 con cv::triangualtePoints, e mettere float
        //    for(int j=0;j<out.rows;j++){ //out.rows cv::triangulatePoints
        //        Eigen::Vector3d p(out.at<cv::Vec3d>(j,1)[0],out.at<cv::Vec3d>(j,1)[1],out.at<cv::Vec3d>(j,1)[2]);
        PointsXYZ.push_back(p);
        //                true_points(j,0) = out.at<cv::Vec3f>(j,1)[0]; MODO DI TARIQ PROVARE
        //                true_points(j,1) = out.at<cv::Vec3f>(j,1)[1];
        //                true_points(j,2) = out.at<cv::Vec3f>(j,1)[2];
    }

    std::cout<<"PUNTI 3D INIZIALIZZATI! Size: "<<PointsXYZ.size()<<std::endl;

}

bool MapBuilder::cvSba(){
    cvsba::Sba sba;
    std::vector <cv::Point3f> points3BA;
    int NCAMS=images.size();
    std::vector< cv::Mat > cameraMatrix, distCoeffs, RBA, TBA;
    for(int j=0;j<this->PointsXYZ.size();j++){
        cv::Point3f p;
        p.x=this->PointsXYZ[j](0);
        p.y=this->PointsXYZ[j](1);
        p.z=this->PointsXYZ[j](2);
        points3BA.push_back(p);
    }



    // fill camera intrinsics (same intrinsics for all cameras)
    cameraMatrix.resize(NCAMS);
    for(int i=0; i<NCAMS; i++)
        cameraMatrix[i] = cv::Mat::eye(3,3,CV_32FC1);

    // fill distortion (assume no distortion)
    distCoeffs.resize(NCAMS);
    for(int i=0; i<NCAMS; i++) distCoeffs[i] = cv::Mat(5,1,CV_32FC1, cv::Scalar::all(0));

    // fill rotation and translation
    RBA.resize(NCAMS);
    RBA[0] = cv::Mat(3,1,CV_32FC1,cv::Scalar::all(0));
    for(int j=1;j<NCAMS;j++){

        cv::Rodrigues(Rotations[j](cv::Rect(0,0,3,3)),RBA[j]);
    }


    cvsba::Sba::Params parameters;
    parameters.type=cvsba::Sba::MOTIONSTRUCTURE;
    parameters.fixedDistortion=5;// 5 non cambio ne distorsione ne camera matrix, se metto 0 le cambia tutte
    parameters.fixedIntrinsics=5;
    parameters.iterations=100;//num massimo di iterazioni
    parameters.minError=0.001;//partiamo che il reprojection error deve venire 1.
    parameters.verbose=true;//Ti dice cose durante la ottimizzazione
    sba.setParams(parameters);

    std::ofstream visfile;
    visfile.open("VisibAnalysis");

    for(int i=1;i<visibility.size();i++){
        int count=0;
        for(int j=0;j<visibility[i].size();j++){
            if(visibility[i][j]==1)
                count++;
        }
        visfile<<"Riga "<<i<<" size "<< visibility[i].size()<<" number of ones  "<<count<<" number of points "<<points[i].size()<<std::endl;
    }

    /****** RUN BUNDLE ADJUSTMENT ******/
    try{
        sba.run(points3BA,  this->points, this->visibility,  cameraMatrix,  RBA,  this->Translations, distCoeffs);
    }
    catch(cv::Exception& e){
        cout << "CV exception: " << e.what() << endl;
        return false;
    }

    std::cout<<"Initial error="<<sba.getInitialReprjError()<<". Final error="<<sba.getFinalReprjError()<<std::endl;
    if(sba.getFinalReprjError()>10){
        std::cout<<"BA FAILED!"<<std::endl;
        return false;
    }

    for(int j=0;j<this->PointsXYZ.size();j++){
        PointsXYZ[j](0)=points3BA[j].x;
        PointsXYZ[j](1)=points3BA[j].y;
        PointsXYZ[j](2)=points3BA[j].z;

    }



    return true;

}

Eigen::Vector2d project2d(const Eigen::Vector3d& v){
    Eigen::Vector2d res;
    res(0) = v(0)/v(2);
    res(1) = v(1)/v(2);
    return res;
}

Eigen::Vector3d unproject2d(const Eigen::Vector2d& v){
    Eigen::Vector3d res;
    res(0) = v(0);
    res(1) = v(1);
    res(2) = 1;
    return res;
}

inline Eigen::Vector3d invert_depth(const Eigen::Vector3d & x){
    return unproject2d(x.head<2>())/x[2];
}

bool MapBuilder::g2oBa(){

    //g2o BA
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    std::cerr << "Using CHOLMOD" << std::endl;
    linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

    optimizer.setAlgorithm(solver);

    //truepoints e' gia' pointXYZ

    Eigen::Vector2d focal_length(1.0,1.0); // pixels default 500,500 Se i punti sono calibrati, sono giusti cosi f 1 e pp 0
    Eigen::Vector2d principal_point(0.0,0.0); // 640x480 image default:principal_point(320,240)






    return true;
}

void MapBuilder::writeFileCeres(){
    std::ofstream file, filepunti;
    filepunti.open("Prima.txt");
    file.open("Problem.txt");
    int num_obs=0;
    for(int i=0;i<visibility.size();i++){
        for(int j=0; j<visibility[i].size();j++){
            if(visibility[i][j]==1)
                num_obs++;
        }
    }
    file<<images.size()<<" "<<visibility[0].size()<<" "<<num_obs<<std::endl;//header
    for (int i=0;i<visibility[0].size();i++){
        for(int j=0;j<visibility.size();j++){
            if(visibility[j][i]==1){
                file<<j<<" "<<i<<" "<<points[j][i].x<<" "<<points[j][i].y<<std::endl;//camera vs point, aka visibility
            }
        }
    }
    for(int i=0;i<images.size();i++){ //camera parameters(ext and intr)
        cv::Mat rotvec(3,1,CV_32FC1);
        cv::Rodrigues(Rotations[i],rotvec);
        file<<rotvec.at<float>(0,0)<<std::endl<<rotvec.at<float>(1,0)<<std::endl<<rotvec.at<float>(2,0)<<std::endl;//Rot
        //        std::cout<<rotvec<<std::endl; //Ok accesso sia sopra che sotto.
        file<<Translations[i].at<float>(0,0)<<std::endl<<Translations[i].at<float>(1,0)<<std::endl<<Translations[i].at<float>(2,0)<<std::endl;
        //        std::cout<<Translations[i]<<std::endl; //ok accesso
        file<<0.0<<std::endl<<0.0<<std::endl<<0.0<<std::endl<<0.0<<std::endl;
    }
    for(int i=0;i<PointsXYZ.size();i++)
    {
        file<<PointsXYZ[i](0)<<std::endl<<PointsXYZ[i](1)<<std::endl<<PointsXYZ[i](2)<<std::endl; // 3D points
        filepunti<<PointsXYZ[i](0)<<" "<<PointsXYZ[i](1)<<" "<<PointsXYZ[i](2)<<std::endl;//Ok accesso sopra.
    }
}

bool MapBuilder::ceresBa(){

    ceres::examples::BALProblem bal_problem;
    if (!bal_problem.LoadFile("Problem.txt")) {
        std::cerr << "ERROR: unable to open file " << "Problem.txt" << "\n";
        return 1;
    }
    const double* observations = bal_problem.observations();
    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    ceres::Problem problem;
    for (int i = 0; i < bal_problem.num_observations(); ++i) {
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
                ceres::MyReprojectionError::Create(observations[2 * i + 0],
                observations[2 * i + 1]);
        problem.AddResidualBlock(cost_function,
                                 NULL /* squared loss */,
                                 bal_problem.mutable_camera_for_observation(i),
                                 bal_problem.mutable_point_for_observation(i));
    }
    // Make Ceres automatically detect the bundle structure. Note that the
    // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
    // for standard bundle adjustment problems.
    ceres::Solver::Options options;
    options.num_linear_solver_threads=4;
    options.update_state_every_iteration=true;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    options.minimizer_progress_to_stdout = true;
//    options.parameter_tolerance=1e-40; //questo e' il parametro per la convergenza, indica quanto deve essere piccola la differenza tra lo stato attuale
                                        // e quello successivo per considerare il raggiungimento di un minimo. Di default la soglia e' 1e-8
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";
    PointsXYZ.clear();
    PointsXYZ=bal_problem.write3Dpoints();

    //    std::cout<<"SIZE DOPO CERES "<<PointsXYZ.size()<<std::endl; // ok la dimensione e' ok
    std::ofstream filedopo;
    filedopo.open("Dopo.txt");
    for(int i=0;i<PointsXYZ.size();i++){

        filedopo<<PointsXYZ[i](0)<<" "<<PointsXYZ[i](1)<<" "<<PointsXYZ[i](2)<<std::endl;
    }



    return true;
}

int MapBuilder::getTheFarthestFrame(){//deprecated, use getPointsR
    std::vector<float> norms;
    for(int i=0; i<Translations.size();i++){
        cv::Vec3f t(Translations[i]);
        //        std::cout<<t<<std::endl; // ok
        norms.push_back(cv::norm(t));
    }
    std::vector<float>::iterator it;
    //    std::cout<<norms<<std::endl; // ok controllato  anche con octave
    it=std::max_element(norms.begin(),norms.end());
    return std::distance(norms.begin(), it);
}

void MapBuilder::visualize3DMap(){

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Tranformations"));
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    //    std::ofstream f2;

    //    f2.open("Prova2.txt");

    for(int j=0;j<this->PointsXYZ.size();j++)
    {
        //        f2<<PointsXYZ[j]<<std::endl<<std::endl; //OK
        //        f2<<PointsXYZ[j](0)<<" "<<PointsXYZ[j](1)<<" "<<PointsXYZ[j](2)<<" "<<std::endl<<std::endl<<std::endl;

        pcl::PointXYZRGB point;
        point.x = this->PointsXYZ[j](0);
        point.y = this->PointsXYZ[j](1);
        point.z = this->PointsXYZ[j](2);

        if(point.z<0 || point.z>2000){//2000 mm-> 2 m
            //            point.r = 255;
            //            point.g = 0;
            //            point.b = 0;
            continue;
        }
        else
        {
            point.r = 255;
            point.g = 255;
            point.b = 255;
            //            if(point.z<500){
            //                point.r = 255/**point.z/2000*/;
            //                point.g = 0;
            //                point.b = 0;
            //            }
            //            else if(point.z>500 && point.z<1500){
            //                point.r = 0/**point.z/2000*/;
            //                point.g = 255;
            //                point.b = 0;
            //            }
            //            else if(point.z>1500 && point.z<2000){
            //                point.r = 0/**point.z/2000*/;
            //                point.g = 0;
            //                point.b = 255;
            //            }
        }


        //if(/*point.z>0 &&*/ (ReprErrorvec[j]<(double)meanReprErr.val[0]+3*(double)stdErr.val[0] || ReprErrorvec[j]>(double)meanReprErr.val[0]-3*(double)stdErr.val[0]))
        cloud -> points.push_back(point); // I keep in the point cloud only points with a positive z.


    }

    cloud->width = (int)cloud->points.size();
    cloud->height = 1;

    if (cloud->width>0)
    {
        pcl::io::savePCDFileASCII("PointsXYZ.pcd", *cloud);
        //        viewer->addCoordinateSystem(1.0, 0);// Non funziona piu' non so perche'
        viewer->addPointCloud(cloud);
        while(!viewer->wasStopped()){
            viewer->spinOnce();
        }

    }


}




