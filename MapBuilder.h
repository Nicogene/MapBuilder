#ifndef MAPBUILDER_H
#define MAPBUILDER_H
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <yarp/os/all.h>
#include <Eigen/StdVector>

#define G2O 0
#define CVSBA 1
#define CERESBA 2


class MapBuilder
{
private:
    int NumOfCouples;
    std::vector<cv::Mat> images;
    bool first;
    double start;
    int methodBA;
    bool WithYarpDataPlayer;
    cv::Mat CameraMatrixL = cv::Mat::eye(3,3,CV_64FC1);
    cv::Mat CameraMatrixR = cv::Mat::eye(3,3,CV_64FC1);
    cv::Mat distCoeffsL = cv::Mat(5,1,CV_64FC1, cv::Scalar::all(0));
    cv::Mat distCoeffsR = cv::Mat(5,1,CV_64FC1, cv::Scalar::all(0));
    std::vector<std::vector<int>> buckets;
    std::vector<std::vector<cv::Point2f>> points;
    std::vector<std::vector< int >> visibility;
    std::vector<cv::Mat> Rotations;
    std::vector<cv::Mat> Translations;
    std::vector<cv::Mat> ProjectionMatrices;
    std::vector<Eigen::Vector3d> PointsXYZ;


public:
    MapBuilder();
    MapBuilder(int numC, bool WithYDplayer, int mBA=-1);
    ~MapBuilder();
    void setNumOfCouples(int numC);
    int getNumOfCouples();
    void setWithYarpDataPlayer(bool WYDP);
    bool getWithYarpDataPlayer();
    void setCameraMatrices(double fxL, double fyL, double cxL, double cyL,
                           double fxR, double fyR, double cxR, double cyR);
    void setDistorsionMatrices(double d0L,double d1L, double d2L, double d3L, double d4L,
                               double d0R,double d1R, double d2R, double d3R, double d4R);
    //get cameraMatrix, get distorsion MISSING farli se servono
    bool process();
    bool processCouples();
protected:
    void allocateMem();
    void collectImages();
    bool findPoints();
    void removeRowsColsVisibility();
    bool findPointsCouples();
    std::vector<std::vector<int>> bucketing(std::vector<cv::Point2f> pts, int width, int height);
    cv::Mat myTriangulate(cv::Mat Rot,cv::Mat t, cv::Mat pts1, cv::Mat pts2 );
    void getCouplesTransformationsThroughEssential();
    void getRotationsThroughEssential();
    cv::Mat featureSelection(int numIter, std::vector<cv::Point2f> pointsL, std::vector<cv::Point2f> pointsR);
    bool getTransformationsToRoot(yarp::os::Network yarp, int i);
    std::vector<cv::Point2d> getPointsR(std::vector<int> *indeces);
    void checkPointIsVisible();
    void initialize3DPoints();
    int getTheFarthestFrame();
    bool cvSba();
    bool g2oBa();
    void writeFileCeres();
    bool ceresBa();
    void visualize3DMap();
    void getTransformationsToFirstLeft();
    void getTransformationsCouples();

};
 
#endif
