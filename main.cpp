#include "MapBuilder.h"
#include "icubimgkingrabber.h"
#include "MatrixVecAlgebra.h"
#include <iostream>

int main()
{
//    MapBuilder m(10,true,CVSBA);
//    while(true){

//        if(!m.process())
//            break;
//    }
    iCubImgKinGrabber m(5,true,480,640);
    std::vector<std::vector<u_char>> images;
    std::vector<std::vector<double>> proj;
    m.process(images,proj);

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
