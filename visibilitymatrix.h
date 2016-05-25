#ifndef VISIBILITYMATRIX_H
#define VISIBILITYMATRIX_H

#include <vector>
#include <iostream>

class VisibilityMatrix
{
    std::vector<std::vector<int>> vmat;
public:
    VisibilityMatrix();
    VisibilityMatrix(int numCam);

    std::vector<std::vector<int>> getVMat(std::vector<std::vector<double>> &pts, std::vector<std::vector<u_char>> &status,
                                          std::vector<std::vector<float>> &error, std::vector<std::vector<double>>& ProjectionMatrices,
                                          std::vector<std::vector<u_char>>& imgs);
protected:
    void buildVisibility(std::vector<std::vector<u_char>> &status, std::vector<std::vector<float>> &error);
    void removeRowsColsVisibility(std::vector<std::vector<double>> &pts, std::vector<std::vector<double>>& ProjectionMatrices,std::vector<std::vector<u_char>>& imgs);


};

#endif // VISIBILITYMATRIX_H
