#include "visibilitymatrix.h"
#include <cmath>



VisibilityMatrix::VisibilityMatrix()
{
}

VisibilityMatrix::VisibilityMatrix(int numCam)
{
    vmat.resize(numCam);
}

std::vector<std::vector<int>> VisibilityMatrix::getVMat(std::vector<std::vector<double>> &pts, std::vector<std::vector<u_char>> &status,
                                      std::vector<std::vector<float>> &error, std::vector<std::vector<double>>& ProjectionMatrices,
                                        std::vector<std::vector<u_char>>& imgs ){
                              buildVisibility(status,error);
                              removeRowsColsVisibility(pts,ProjectionMatrices,imgs);
                              return vmat;
                              }
void VisibilityMatrix::buildVisibility(std::vector<std::vector<u_char> > &status,
                                       std::vector<std::vector<float> > &error){

                              for(int i=0;i<vmat.size();i++){
    std::vector<int> vector1(status[0].size(), 1);//Initialize as matrix of ones
    vmat[i]=vector1;//each pts[i] has the same size
}
    for(int i=1;i<vmat.size();i++){
        for (int j=0;j<status[0].size();j++){        // non posso piu' usare questo vincolo, due immagini possono anche essere lontane

            if((std::isnan(error[i-1][j])) || (error[i-1][j])> 1000 || (int)status[i-1][j]==0){//status and error have not the element of the first frame
                vmat[i][j]=0;
            }

        }
}
}
void VisibilityMatrix::removeRowsColsVisibility(std::vector<std::vector<double> > &pts, std::vector<std::vector<double> > &ProjectionMatrices,
                                                std::vector<std::vector<u_char> > &imgs){
    //remove rows
    for(int i=0;i<this->vmat.size();i++){
        int sum=0;
        for (int n : vmat[i])
            sum += n;
        if(sum<10){
            ProjectionMatrices.erase (ProjectionMatrices.begin()+i);
            imgs.erase (imgs.begin()+i);
            vmat.erase (vmat.begin()+i);
            pts.erase (pts.begin()+i);

            //            std::cout<<"Ciao T "<<sum<<std::endl;
            std::cout<<"Frame erased, too few correspondences"<<std::endl;
            i--;
        }
        //        else
        //            std::cout<<"Ciao F "<<sum<<std::endl;
    }
    for(int i=0;i<vmat[0].size();i++){
        int sum=0;
        for(int j=0;j<vmat.size();j++){
            sum+=vmat[j][i];
        }
        if(sum<2){
            for(int j=0;j<vmat.size();j++){
                vmat[j].erase (vmat[j].begin()+i);
                pts[j].erase(pts[j].begin()+2*i);//x
                pts[j].erase(pts[j].begin()+2*i);//y //Giusto senza il +1, il vettore si accorcia, e alla posizione 2*i ci va quello 2*i+1
            }
//            std::cout<<"PRIMA "<<points[0].size()<<std::endl; //OK TESTED, FA quello che serve
//            std::cout<<"DOPO "<<points[0].size()<<std::endl;
            i--;
            //std::cout<<"Column of visibility erased no correspondence for the point "<<i<<std::endl;
        }
    }

                              }

