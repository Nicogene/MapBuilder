#include "MatrixVecAlgebra.h"
#include <math.h>
#include <vector>

void multiply(std::vector<double> &C, std::vector<double> &A, std::vector<double> &B, int nrows, int ncomm, int ncols){
    // C = A*B
    // multiplies (nrows x ncomm) x (ncomm x ncols) = (nrows x ncols)
    for (int j = 0; j < ncols; j++){ // ncols of result
        for (int i = 0; i < nrows; i++){ // nrows of result
            for (int k = 0; k < ncomm; k++){ // common dimension
                C[nrows*j + i] = C[nrows*j + i]
                        + A[nrows*k + i] * B[ncomm*j + k];
            };
        };
    };
    for (int j = 0; j < nrows*ncols; j++){ // ncols of result
        if (sqrt(C[j]*C[j]) < 0.00000001)
            C[j] = 0;
    }
}

std::vector<double> transpose(std::vector<double> &A,int r, int c) //TESTED WORKS
{
    std::vector<double> res(r*c);

    for(int i=0;i<r;i++)
        for(int j=0;j<c;j++)
            res[i+c*j]=A[j+c*i];

    return res;

}

void inverseProjMat(std::vector<double> &inv, std::vector<double> &A){
    std::vector<double> R(9),Rtrans(9),t(3),Rt(3);

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            R[i+3*j]=A[i+4*j];
    Rtrans=transpose(R,3,3);
//    std::cout<<"R "<<R<<std::endl<<"Rtrans "<<Rtrans<<std::endl;//OK
    std::vector<double> Rtransneg(9);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Rtransneg[i+3*j]=-Rtrans[i+3*j];

    for(int i=0;i<3;i++)
        t[i]=A[i+3*4];//OK
//    std::cout<<"Neg "<<Rtransneg<<std::endl;//OK
//    std::cout<<"t "<<t<<std::endl;//OK


    multiply(Rt,Rtransneg,t,3,3,1);

//    std::cout<<"Rt "<<Rt<<std::endl; //OK

    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++){
            if(i<3 && j<3)
                inv[i+4*j]=Rtrans[i+3*j];//OK
            else if(i>2 && j<3){
                inv[j+4*i]=Rt[j+3*(i-3)];//OK
//                std::cout<<i<<" "<<j<<std::endl;
//                std::cout<<"index inv "<<i+4*j<<std::endl;
//                std::cout<<"index Rt "<<j+3*(i-3)<<std::endl;
            }
            else
                inv[i+4*j]=A[i+4*j];//OK
        }


}
