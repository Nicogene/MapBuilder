#include "MapBuilder.h"
#include "icubimgkingrabber.h"
#include "featurefinderandtracker.h"
#include "visibilitymatrix.h"
#include "triangulationclass.h"
#include "MatrixVecAlgebra.h"
#include <iostream>
#include <fstream>

#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

int main()
{
//    MapBuilder m(10,true,CVSBA);
//    while(true){

//        if(!m.process())
//            break;
//    }
    //while(true){
    iCubImgKinGrabber m(5,true,480,640);
    std::vector<std::vector<u_char>> images,status;
    std::vector<std::vector<double>> proj,points;
    std::vector<std::vector<float>> error;

    m.process(images,proj);
    FeatureFinderAndTracker f(images,480,640);
    f.process(points,status,error);
    VisibilityMatrix v(images.size());
    std::vector<std::vector<int>> vmat;
//    std::cout<<"Prima "<<std::endl;//OK TESTED
//    std::cout<<0<<" "<<status[0].size()<<" "<<points[0].size()<<std::endl;
//    for(int i=0;i<status.size();i++){
//        std::cout<<i+1<<" "<<status[i].size()<<" "<<points[i].size()<<std::endl;
//    }


    vmat=v.getVMat(points,status,error,proj,images);


//    std::vector<double> prova={1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16},R,T;//OK TESTED
//    getRandT(prova,R,T);
//    cv::Mat cvrot(3,3,CV_64F,R.data());
//    cv::Mat cvprova(4,4,CV_64FC1,prova.data());

//    std::cout<<"Proj "<<std::endl<<cvprova.t()<<std::endl<<"Rot "<<std::endl<<cvrot.t()<<std::endl<<"T "<<std::endl<<T<<std::endl;
    std::vector<std::vector<double>> map;
    TriangulationClass t(points[0].size()/2);
    map=t.get3DPoints(points,proj,vmat);
    std::cout<<map.size()<<"x"<<map[0].size()<<std::endl;
    map=t.optimizeWithCeres(points,proj,vmat);


    std::ofstream filePunti;
    filePunti.open("Punti.txt");
    for(int i=0;i<map.size();i++)
        filePunti<<map[i][0]<<" "<<map[i][1]<<" "<<map[i][2]<<std::endl;

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Tranformations"));
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);


    for(int j=0;j<map.size();j++)
    {

        pcl::PointXYZRGB point;
        point.x = map[j][0];
        point.y = map[j][1];
        point.z = map[j][2];

        if(point.z<0 || point.z>2000){//2000 mm-> 2 m
            continue;
        }
        else
        {
            point.r = 255;
            point.g = 255;
            point.b = 255;

        }



        cloud -> points.push_back(point); // I keep in the point cloud only points with a positive z.


    }

    cloud->width = (int)cloud->points.size();
    cloud->height = 1;

    if (cloud->width>0)
    {

        viewer->addPointCloud(cloud);
        while(!viewer->wasStopped()){
            viewer->spinOnce();
        }

    }
//}//e' la parentesi del while togli il commento se vuoi il loop




//    std::cout<<"Dopo "<<std::endl;//OK TESTED
//    for(int i=0;i<vmat.size();i++){
//        std::cout<<i<<" "<<vmat[i].size()<<" "<<points[i].size()<<std::endl;
//    }



//    for(int i=0;i<images.size()-1;i++){//OK points[i] e' il doppio di status e error
//        std::cout<<i<<" status "<<status[i].size()<<std::endl;
//        std::cout<<" error "<<error[i].size()<<std::endl;
//        std::cout<<" points "<<points[i].size()<<std::endl;
//        std::ofstream filePunti;
//        std::stringstream nome;
//        nome<<"punti_"<<i<<".txt";
//        filePunti.open(nome.str());
//        for(int j=0;j<error[i].size();j++){
//            filePunti<<points[i][2*j]<<" "<<points[i][2*j+1]<<std::endl; //OK sono diversi da frame a frame
//        }

//    }






//    for(int i=0;i<proj.size();i++){
                                                //OK le sinistre sono eye, le destre hanno l'offset giusto in t
//        cv::Mat m(4,4,CV_64F,proj[i].data());

//        std::cout<<i<<" "<<m.t()<<std::endl;//sempre per la cosa del major column

//    }
//    std::cout<<proj.size()<<std::endl;
//    for(int i=0;i<proj.size();i++){
//        std::cout<<i<<" "<<proj[i].size()<<std::endl; //OK TESTED
//    }
//    std::cout<<proj[3]<<std::endl;
//    cv::Mat A(4,4,CV_64FC1,proj[3].data());
//    std::cout<<A<<std::endl;
//    cv::Mat A= (cv::Mat_<double>(4,4) << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 0, 0, 0, 1.0);
//    std::cout<<A<<std::endl;
//    A=A.t();//column major
//    std::vector<double> v(16);

//    v.assign((double*)A.datastart, (double*)A.dataend);
//    std::vector<double> vt(16);
//    inverseProjMat(vt,v);

//    cv::Mat B(4,4,CV_64FC1,vt.data());
//    std::cout<<B.t()<<std::endl;//perche' il vector e' in major column, cv::Mat e' in major row


}
