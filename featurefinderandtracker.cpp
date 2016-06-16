#include "featurefinderandtracker.h"

extern "C" {
#include "fast.h"
}

FeatureFinderAndTracker::FeatureFinderAndTracker()
{
}

FeatureFinderAndTracker::FeatureFinderAndTracker(std::vector<std::vector<u_char>> &imgs, int r, int c){
    images=imgs;
    rows=r;
    cols=c;

}

bool FeatureFinderAndTracker::process(std::vector<std::vector<double>> &points, std::vector<std::vector<uchar>> &status,
                                      std::vector<std::vector<float>> &error){

    std::vector<cv::Mat> imgs;
    restoreCvMat(imgs);
//    std::cout<<"process "<<imgs.size()<<std::endl;
    findPoints(imgs,points,status,error);

}

bool FeatureFinderAndTracker::findPoints(std::vector<cv::Mat> &imgs, std::vector<std::vector<double> > &pts,
                                         std::vector<std::vector<uchar> > &status, std::vector<std::vector<float> > &error){
    bool res=true;
//    std::cout<<imgs.size()<<std::endl;
    std::vector<std::vector<cv::Point2f>> points(imgs.size());
    pts.resize(imgs.size());
    status.resize(imgs.size()-1);//non ha senso quello riferito al primo frame
    error.resize(imgs.size()-1);
    cv::Mat CameraMatrixL, CameraMatrixR, distCoeffsL, distCoeffsR;


    this->fillInstrinsics(CameraMatrixL, CameraMatrixR, distCoeffsL, distCoeffsR);

    std::vector<cv::KeyPoint> keypointsL;
    xy* features;//ciao
    int numcorn;
    features=fast9_detect(images[0].data(),640,480,640,15,&numcorn);// byte (aka unsigned char), se non ci sono 0 padding strides e' uguale a xsize
//    std::cout<<features[0].x<<" "<<features[0].y<<" numPoints "<< numcorn<<std::endl;//OK accesso.
    points[0].resize(numcorn);
    for(int i=0;i<points[0].size();i++){
        points[0][i].x=features[i].x;
        points[0][i].y=features[i].y;
    }
//    cv::FAST(imgs[0],keypointsL,15,false,cv::FastFeatureDetector::TYPE_9_16);
    if(/*keypointsL.size()>0*/numcorn>0){
//        cv::KeyPoint::convert(keypointsL, points[0]);
        for(int i=1;i<points.size();i++)
        {
            std::vector<uchar> sts;
            cv::Mat err;
            cv::TermCriteria termcrit(cv::TermCriteria::COUNT|cv::TermCriteria::EPS,20,0.03);
            cv::Size winSize(31,31);

            cv::calcOpticalFlowPyrLK(imgs[0],imgs[i],points[0],points[i],sts,err,winSize);


            if(points[i].size()>0){
                for(int j=0;j<sts.size();j++){
                    if(abs(points[0][j].x-points[i][j].x)<10 && abs(points[0][j].y-points[i][j].y)<10)
                        sts[j]=(uchar)0;
                }
//                std::cout<<CameraMatrixL<<std::endl<<CameraMatrixR<<std::endl<<distCoeffsL<<std::endl<<distCoeffsR<<std::endl; //OK
                if(i%2==0)
                    cv::undistortPoints(points[i],points[i],CameraMatrixL,distCoeffsL);
                else
                    cv::undistortPoints(points[i],points[i],CameraMatrixR,distCoeffsR);
                cvtPointsToVector(points[i],pts[i]);
                status[i-1]=sts;
                error[i-1].assign((float*)err.datastart, (float*)err.dataend);

            }
            else return false;


        }
        //        std::vector<cv::Point2f> tmp;
        cv::undistortPoints(points[0],points[0],CameraMatrixL,distCoeffsL);
        cvtPointsToVector(points[0],pts[0]);
    }
    else
        return false;



    return res;


}

void FeatureFinderAndTracker::fillInstrinsics(cv::Mat &CL, cv::Mat &CR, cv::Mat &dL, cv::Mat &dR){


        CL = (cv::Mat_<double>(3,3) << 407.386849650585816, 0, 335.142950296626907, 0, 408.446478587298941, 247.077012143219406, 0, 0, 1);
        CR = (cv::Mat_<double>(3,3) << 411.389187422896157, 0,  320.913970488683390, 0, 412.513836073373398, 224.679749311050074, 0, 0, 1);

        dL= (cv::Mat_<double>(5,1) << -0.370087506134749, 0.122569815417904, -0.000557593423888,
             -0.000650939992441, 0.000000000000000);
        dR= (cv::Mat_<double>(5,1) << -0.381727559529518, 0.138669468725887, -0.000320490995377,
             -0.000398273646263, 0.000000000000000);

}

void FeatureFinderAndTracker::restoreCvMat(std::vector<cv::Mat> &imgs){
    imgs.resize(images.size());
    for(int i=0;i<images.size();i++){
        cv::Mat m(rows,cols,CV_8UC1,images[i].data());
        imgs[i]=m;
//        cv::imshow("Prova",imgs[i]);
//        cv::waitKey(300);// OK!
    }
}

void FeatureFinderAndTracker::cvtPointsToVector(std::vector<cv::Point2f> &cvpts, std::vector<double> &pts){
    pts.resize(cvpts.size()*2);
    for(int i=0;i<cvpts.size();i++){
        pts[2*i]=cvpts[i].x;
        pts[2*i+1]=cvpts[i].y;
    }

}

