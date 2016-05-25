#ifndef TRIANGULATIONCLASS_H
#define TRIANGULATIONCLASS_H
#include <vector>
#include <opencv2/opencv.hpp>
class TriangulationClass
{
    std::vector<std::vector<double>> ptsXYZ;
public:
    TriangulationClass();
    TriangulationClass(int numP);
    std::vector<std::vector<double>> get3DPoints(std::vector<std::vector<double>> &pts, std::vector<std::vector<double>> &ProjectionMatrices,
                                                 std::vector<std::vector<int>> &vmat);
protected:
    void cvtVectorToPoints(std::vector<cv::Point2d> &cvpts, std::vector<double> &pts);
    cv::Mat myTriangulate(cv::Mat Rot,cv::Mat t, cv::Mat x1, cv::Mat x2 );
    std::vector<cv::Point2d> getPointsR(std::vector<std::vector<double>> &pts,std::vector<int> &indeces,std::vector<std::vector<double>> &ProjectionMatrices,std::vector<std::vector<int>> &vmat);

};

#endif // TRIANGULATIONCLASS_H
