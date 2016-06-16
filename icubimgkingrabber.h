#ifndef ICUBIMGKINGRABBER_H
#define ICUBIMGKINGRABBER_H
#include <vector>
#include <yarp/os/all.h>


class iCubImgKinGrabber
{
private:
    int NumOfCouples;
    int rows;
    int cols;
    bool first;
    double start;
    bool WithYarpDataPlayer;

public:
    iCubImgKinGrabber();
    iCubImgKinGrabber(int numC, bool WYDP=false, int r=0, int c=0);
    bool process(std::vector<std::vector<u_char>>& images, std::vector<std::vector<double>>& ProjectionMatrices);
    void setNumOfCouples(int numC);
    int getNumOfCouples();
    void setWithYarpDataPlayer(bool WYDP);
    bool getWithYarpDataPlayer();
protected:
    std::vector<std::vector<u_char>> collectImages(std::vector<std::vector<double>>& ProjectionMatrices);//Anche le immagini srotolate
    bool getTransformationsToRoot(yarp::os::Network yarp, int i,std::vector<std::vector<double>>& ProjectionMatrices);//Matrici srotolate in vettori, major column
    void getTransformationsToFirstLeft(std::vector<std::vector<double>>& ProjectionMatrices);
    void checkZeroRotations(std::vector<double> &proj);
};

#endif // ICUBIMGKINGRABBER_H
