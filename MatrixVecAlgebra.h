#ifndef MATRIXVECALGEBRA_H
#define MATRIXVECALGEBRA_H

#include"icubimgkingrabber.h"

void multiply(std::vector<double> &C, std::vector<double> &A, std::vector<double> &B, int nrows, int ncomm, int ncols);
std::vector<double> transpose(std::vector<double> &A,int r, int c);
void inverseProjMat(std::vector<double> &inv, std::vector<double> &A);

#endif // MATRIXVECALGEBRA_H
