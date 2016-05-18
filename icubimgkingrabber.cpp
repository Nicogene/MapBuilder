#include "icubimgkingrabber.h"
#include "MatrixVecAlgebra.h"
#include <yarp/sig/Image.h>
#include <yarp/dev/PolyDriver.h>
#include <yarp/dev/IEncoders.h>
#include <yarp/os/Time.h>
#include <yarp/sig/Vector.h>
#include <yarp/math/Math.h>
#include <iCub/iKin/iKinFwd.h>
#include <iCub/iKin/iKinIpOpt.h>
#include <opencv2/opencv.hpp>

using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;
using namespace iCub::iKin;

iCubImgKinGrabber::iCubImgKinGrabber()
{
}

iCubImgKinGrabber::iCubImgKinGrabber(int numC, bool WYDP, int r, int c){
    NumOfCouples=numC;
    WithYarpDataPlayer=WYDP;
    rows=r;
    cols=c;
    first=true;
}

void iCubImgKinGrabber::setNumOfCouples(int numC){
    NumOfCouples=numC;
}

int iCubImgKinGrabber::getNumOfCouples(){ return NumOfCouples;}

void iCubImgKinGrabber::setWithYarpDataPlayer(bool WYDP){WithYarpDataPlayer=WYDP;}

bool iCubImgKinGrabber::getWithYarpDataPlayer(){return WithYarpDataPlayer;}

bool iCubImgKinGrabber::process(std::vector<std::vector<u_char>>   &images, std::vector<std::vector<double>> &ProjectionMatrices){

    images=collectImages(ProjectionMatrices);

    return true;

}



std::vector<std::vector<u_char>> iCubImgKinGrabber::collectImages(std::vector<std::vector<double>>& ProjectionMatrices)
{
    std::vector<std::vector<u_char>> imgs;
    imgs.resize(NumOfCouples*2);// TODO change OLD MEM management
    for(int i=0;i<imgs.size();i++){
        imgs[i].resize(rows*cols);

    }
    int i=0;
    bool res=false;
    Network yarp;
    BufferedPort<ImageOf<PixelRgb> > imagePortL,imagePortR;
    int startCountR=0, startCountL=0;
    imagePortL.open("/nico/MapBuilder/image/L");
    imagePortR.open("/nico/MapBuilder/image/R");
    yarp.connect("/icub/cam/left",imagePortL.getName());
    yarp.connect("/icub/cam/right",imagePortR.getName());
    while(i<2*NumOfCouples){

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

                if(this->WithYarpDataPlayer)
                    L = LOld(cv::Rect(3,3,cols,rows)).clone();//Le immagini left dell'acquisizione di Tariq sono 648x488
                else
                    L=LOld.clone();

                cv::medianBlur(R,R,3);
                cv::medianBlur(L,L,3);
                cv::Mat l_gray,r_gray;
                cv::cvtColor(L,l_gray,CV_RGB2GRAY);
                cv::cvtColor(R,r_gray,CV_RGB2GRAY);

                if (l_gray.isContinuous() && r_gray.isContinuous()) {//TODO per ora va bene cosi, magari c'e' da fare una funzione generica
                  imgs[i].assign(l_gray.datastart, l_gray.dataend);
                  imgs[i+1].assign(r_gray.datastart, r_gray.dataend);
                  std::cout<<"Continous "<<i<<std::endl;
                } else {
                  for (int j = 0; j < rows; ++j) {
                    imgs[i].insert(imgs[i].end(), l_gray.ptr<uchar>(j), l_gray.ptr<uchar>(j)+l_gray.cols);
                    imgs[i+1].insert(imgs[i+1].end(), r_gray.ptr<uchar>(j), r_gray.ptr<uchar>(j)+r_gray.cols);
                  }
                }
                getTransformationsToRoot(yarp, i, ProjectionMatrices);
//                std::cout<<imgs[i].size()<<" "<< imgs[i+1].size()<<std::endl;OKTested
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

        getTransformationsToFirstLeft(ProjectionMatrices);
        return imgs;
    }

bool iCubImgKinGrabber::getTransformationsToRoot(yarp::os::Network yarp, int i, std::vector<std::vector<double>> &ProjectionMatrices){
    bool res=true;
    std::cout<<"Ciao1"<<std::endl;
    ProjectionMatrices.resize(2*NumOfCouples);//TODO change OLD MEM management
    for(int j=0;j<ProjectionMatrices.size();j++){
        ProjectionMatrices[j].resize(16);
    }



    double anglesHead[6], anglesTorso[3];
    Bottle  *headBottle, *torsoBottle;
    BufferedPort<Bottle> encHeadPort, encTorsoPort;
    Property optionsHead, optionsTorso;
    yarp::dev::PolyDriver ddHead(optionsHead), ddTorso(optionsTorso);
    yarp::dev::IEncoders *encHead, *encTorso;


    if(!WithYarpDataPlayer){



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

    trasfR=chain->getH(qf);

    chain->setAllConstraints(false);//elimino i vincoli perche' non so quali erano
    //quelli del robot
    chain->clear();

    chain=eyeL.asChain();


    chain->releaseLink(0);
    chain->releaseLink(1);
    chain->releaseLink(2);

    qf[7]=version + vergence/2;//pan

    trasfL=chain->getH(qf);
//OK Tested con il column majors
//    std::cout<<i<<" "<<trasfL(0,0)<<" "<<trasfL(0,1)<<" "<<trasfL(0,2)<<" "<<trasfL(0,0)<<" "<<trasfL(0,3)<<std::endl<<" "<<trasfL(1,0)<<" "<<trasfL(1,1)<<" "<<trasfL(1,2)<<" "<<trasfL(1,3)<<std::endl<<" "<<trasfL(2,0)<<" "<<trasfL(2,1)<<" "<<trasfL(2,2)<<" "<<trasfL(2,3)<<std::endl<<" "<<trasfL(3,0)<<" "<<trasfL(3,1)<<" "<<trasfL(3,2)<<" "<<trasfL(3,3)<<std::endl;

    trasfL=trasfL.transposed();//data e' in row major, faccio la trasposta per metterlo in column
    trasfR=trasfR.transposed();

    std::vector<double> tempL(trasfL.data(),trasfL.data()+16);
    std::vector<double> tempR(trasfR.data(),trasfR.data()+16);

//    std::cout<<"Size projL "<<tempL.size()<<std::endl; //OK 16 elementi.
//    std::cout<<"Size projR "<<tempR.size()<<std::endl;

//    std::cout<<i<<" "<<tempL<<std::endl; // ok tested con il column major di sopra

    ProjectionMatrices[i]=tempL;
    ProjectionMatrices[i+1]=tempR;

//    std::cout<<"Size projL "<<ProjectionMatrices[i].size()<<std::endl;
//    std::cout<<"Size projR "<<ProjectionMatrices[i+1].size()<<std::endl; //OK 16 elementi


    return res;

}

void iCubImgKinGrabber::getTransformationsToFirstLeft(std::vector<std::vector<double>>& ProjectionMatrices){

//    //Trasformations to the first left image


    std::vector<double> m0(16);
    inverseProjMat(m0,ProjectionMatrices[0]);

//    cv::Mat mat0=cv::Mat(4,4,CV_64FC1,ProjectionMatrices[0].data());
//    std::cout<<"Original: "<<mat0.t()<<std::endl;

//    cv::Mat matinv0=cv::Mat(4,4,CV_64FC1,m0.data());
//    std::cout<<"Inverse: "<<matinv0.t()<<std::endl; OK, checked with octave
//    std::vector<double> prod(16);
//    multiply(prod,ProjectionMatrices[0],m0,4,4,4);

//    cv::Mat res=cv::Mat(4,4,CV_64FC1,prod.data());
//    std::cout<<"Res: "<<res.t()<<std::endl;

    for(int i=0;i<ProjectionMatrices.size();i++)
    {
        std::vector<double> tmp(16);
        multiply(tmp,m0,ProjectionMatrices[i],4,4,4);
        ProjectionMatrices[i]=tmp;


    }





}

