#ifndef FEATUREFINDERANDTRACKER_H
#define FEATUREFINDERANDTRACKER_H

#include <vector>
#include <iostream>
#include <opencv2/opencv.hpp>

class FeatureFinderAndTracker
{
    std::vector<std::vector<u_char>> images;
    int rows;
    int cols;

public:
    FeatureFinderAndTracker();
    FeatureFinderAndTracker(std::vector<std::vector<u_char>> &images, int r, int c);
    bool process(std::vector<std::vector<double>> &points, std::vector<std::vector<uchar>> &status, std::vector<std::vector<float>> &error);

protected:
    bool findPoints(std::vector<cv::Mat> &imgs,std::vector<std::vector<double>> &pts, std::vector<std::vector<uchar>> &status, std::vector<std::vector<float>> &error);
    void restoreCvMat(std::vector<cv::Mat> &imgs);
    void fillInstrinsics(cv::Mat &CL,cv::Mat &CR,cv::Mat &dL,cv::Mat &dR);//TODO trovare un modo furbo con cui passare le matrici.
    void cvtPointsToVector(std::vector<cv::Point2f> &cvpts, std::vector<double> &pts);
};

#endif // FEATUREFINDERANDTRACKER_H
