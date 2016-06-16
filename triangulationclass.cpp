#include "triangulationclass.h"
#include "MatrixVecAlgebra.h"
#include "bal_problem.h"
#include <ceres/ceres.h>
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "my_reprojection_error.h"
#include <fstream>

TriangulationClass::TriangulationClass()
{
}

TriangulationClass::TriangulationClass(int numP){
    ptsXYZ.resize(numP);
    for(int i=0;i<ptsXYZ.size();i++){
        ptsXYZ[i].resize(3);
    }
}

void TriangulationClass::getPointsR(std::vector<cv::Point2d> &pointsR,std::vector<std::vector<double>> &pts,std::vector<int> &indeces,std::vector<std::vector<double>> &ProjectionMatrices,std::vector<std::vector<int>> &vmat){
    for(int i=0;i<vmat[0].size();i++){
        std::vector<double> txs;//vettore che contiene le tx, vado a prendere le maggiori.
        for(int j=0; j<vmat.size();j++){

            std::vector<double> R,Rtrans,t,transformedt(3);
            getRandT(ProjectionMatrices[j],R,t);
            Rtrans=transpose(R,3,3);//DA TESTARE

            multiply(transformedt,Rtrans,t,3,3,1);
            txs.push_back(std::fabs(transformedt[0]* vmat[j][i]));//DA TESTARE lo zero nella vmat esclude che quel frame
        }
        int index;
        std::vector<double>::iterator it;
        it=std::max_element(txs.begin(),txs.end());
        index= std::distance(txs.begin(), it);
        indeces.push_back(index);
        if(index>ProjectionMatrices.size()-1)//DA TESTARE
            std::cout<<"Problema, l'indice e'.."<<index<<std::endl;
        if(vmat[index][i]==0)//DA TESTARE
        {
            std::cout<<"Problema, la vmat e' zero per quel punto in quel frame! "<<txs<<std::endl;
        }
        cv::Point2d point;
        point.x=pts[index][2*i];
        point.y=pts[index][2*i+1];
        pointsR.push_back(point);


    }


}

void TriangulationClass::myTriangulate(cv::Mat &xf, cv::Mat Rot, cv::Mat t, cv::Mat x1, cv::Mat x2 ){
    x1.convertTo(x1,CV_64FC1);
    x2.convertTo(x2,CV_64FC1);
    Rot.convertTo(Rot,CV_64FC1);
    t.convertTo(t,CV_64FC1);



    //calculate lines
    cv::Mat v1,v2;
    v1=x1;
    v2= Rot * x2;


    //distance formulas based on Schneider pp 409-412
    cv::Mat a,b,c,d,e,s,denom,num;
    cv::reduce(v1.mul(v1),a,0,CV_REDUCE_SUM);
    cv::reduce(v1.mul(v2),b,0,CV_REDUCE_SUM);
    cv::reduce(v2.mul(v2),c,0,CV_REDUCE_SUM);
    d = t.t() * v1;
    e = t.t() * v2;

    //Non ho implementato questa linea, a che serve?
    //denom(denom < eps) = 1; % accounts for parallel lines

    denom = a.mul(c) - b.mul(b);
    num = c.mul(d) - b.mul(e);
    cv::divide(num,denom,s);


    // compute lines length
    cv:: Mat r1,sqrta;
    cv::sqrt(a,sqrta);
    r1=sqrta.mul(s);

    //TODO mettere a 0 gli elementi della visibility con prof negativa.

    //3d points in P1 coordinates
    cv::Mat rim,z,v1tworows,v1pow,v1sum;
    v1tworows=v1(cv::Rect(0,0,x1.cols,2));
    cv::pow(v1tworows,2.0,v1pow);
    cv::reduce(v1pow,v1sum,0,CV_REDUCE_SUM);
    cv::sqrt(v1sum+1,rim);


    cv::divide(r1,rim,z);

    xf(cv::Rect(0,0,x1.cols,1)) = x1(cv::Rect(0,0,x1.cols,1)).mul(z);
    xf(cv::Rect(0,1,x1.cols,1)) = x1(cv::Rect(0,1,x1.cols,1)).mul(z);
    xf(cv::Rect(0,2,x1.cols,1)) = x1(cv::Rect(0,2,x1.cols,1)).mul(z);

}

void TriangulationClass::cvtVectorToPoints(std::vector<cv::Point2d> &cvpts, std::vector<double> &pts){
    for(int i=0;i<pts.size()/2;i++){
        cv::Point2d point;
        point.x=pts[2*i];
        point.y=pts[2*i+1];
        cvpts.push_back(point);
    }

}

std::vector<std::vector<double>>  TriangulationClass::get3DPoints(std::vector<std::vector<double> > &pts, std::vector<std::vector<double> > &ProjectionMatrices,
                                             std::vector<std::vector<int> > &vmat){
                                  std::vector<cv::Point2d> pointsL,pointsR;
                                  std::vector<int> indeces;
                                  cvtVectorToPoints(pointsL,pts[0]);//OK
                                  getPointsR(pointsR,pts,indeces,ProjectionMatrices,vmat);//OK
                                  std::ofstream provaindici,provaPuntiLR;
                                  provaindici.open("Provaindici.txt");
                                  provaPuntiLR.open("ProvaPuntiLR.txt");
                                  for(int i=0;i<indeces.size();i++){
                                      provaindici<<indeces[i]<<std::endl;
                                      provaPuntiLR<<pointsL[i].x<<" "<<pointsL[i].y<<" "<<pointsR[i].x<<" "<<pointsR[i].y<<std::endl;
                                  }
                                  cv::Mat mL,mR;

                                  mL=cv::Mat(pointsL).reshape(1,2);
                                  mR=cv::Mat(pointsR).reshape(1,2);

                                  mL.convertTo(mL,CV_64F);
                                  mR.convertTo(mR,CV_64F);

                                  cv::Mat o;
                                  o=o.ones(1,mL.cols,CV_64F);

                                  cv::vconcat(mL,o,mL);
                                  cv::vconcat(mR,o,mR);

                                  cv::Mat out(3,1,CV_64FC1);
                                  for(int i=0;i<pointsL.size();i++){

                                      std::vector<double> proj(16);
                                      proj=transpose(ProjectionMatrices[indeces[i]],4,4);//because projectionMatrices is column major

                                      cv::Mat currOut(3,1,CV_64FC1),Proj(4,4,CV_64FC1,proj.data());

                                      if(i==0){


                                          myTriangulate(out,Proj(cv::Rect(0,0,3,3)),Proj(cv::Rect(3,0,1,3))*1000,mL(cv::Rect(i,0,1,3)),mR(cv::Rect(i,0,1,3)));
                                      }
                                      else{
                                          myTriangulate(currOut,Proj(cv::Rect(0,0,3,3)),Proj(cv::Rect(3,0,1,3))*1000,mL(cv::Rect(i,0,1,3)),mR(cv::Rect(i,0,1,3)));
                                          cv::hconcat(out,currOut,out);
                                      }
                                  }


                                  for(int j=0;j<out.cols;j++){

                                      ptsXYZ[j][0]=out.at<double>(0,j);
                                      ptsXYZ[j][1]=out.at<double>(1,j);
                                      ptsXYZ[j][2]=out.at<double>(2,j);

                                  }


                                  std::cout<<"PUNTI 3D INIZIALIZZATI! Size: "<<ptsXYZ.size()<<std::endl;

                                  return ptsXYZ;




                                  }

std::vector<std::vector<double>> TriangulationClass::optimizeWithCeres(std::vector<std::vector<double>> &pts, std::vector<std::vector<double>> &ProjectionMatrices,
                                           std::vector<std::vector<int>> &vmat){


                                 writeFileCeres(pts,ProjectionMatrices,vmat);


                                 ceres::examples::BALProblem bal_problem;
                                 if (!bal_problem.LoadFile("Problem.txt")) {
                                     std::vector<std::vector<double>> v;
                                     std::cerr << "ERROR: unable to open file " << "Problem.txt" << "\n";
                                     return v;
                                 }
                                 const double* observations = bal_problem.observations();
                                 // Create residuals for each observation in the bundle adjustment problem. The
                                 // parameters for cameras and points are added automatically.
                                 ceres::Problem problem;
                                 for (int i = 0; i < bal_problem.num_observations(); ++i) {
                                     // Each Residual block takes a point and a camera as input and outputs a 2
                                     // dimensional residual. Internally, the cost function stores the observed
                                     // image location and compares the reprojection against the observation.
                                     ceres::CostFunction* cost_function =
                                             ceres::MyReprojectionError::Create(observations[2 * i + 0],
                                             observations[2 * i + 1]);
                                     problem.AddResidualBlock(cost_function,
                                                              NULL /* squared loss */,
                                                              bal_problem.mutable_camera_for_observation(i),
                                                              bal_problem.mutable_point_for_observation(i));
                                 }
                                 // Make Ceres automatically detect the bundle structure. Note that the
                                 // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
                                 // for standard bundle adjustment problems.
                                 ceres::Solver::Options options;
                                 options.num_linear_solver_threads=4;
                                 options.update_state_every_iteration=true;
                                 options.linear_solver_type = ceres::SPARSE_SCHUR;
                                 options.minimizer_progress_to_stdout = true;
                             //    options.parameter_tolerance=1e-40; //questo e' il parametro per la convergenza, indica quanto deve essere piccola la differenza tra lo stato attuale
                                                                     // e quello successivo per considerare il raggiungimento di un minimo. Di default la soglia e' 1e-8
                                 ceres::Solver::Summary summary;
                                 ceres::Solve(options, &problem, &summary);
                                 std::cout << summary.FullReport() << "\n";
                                 bal_problem.write3Dpoints(ptsXYZ);

                                 //    std::cout<<"SIZE DOPO CERES "<<PointsXYZ.size()<<std::endl; // ok la dimensione e' ok
//                                 std::ofstream filedopo;
//                                 filedopo.open("Dopo.txt");
//                                 for(int i=0;i<PointsXYZ.size();i++){

//                                     filedopo<<PointsXYZ[i](0)<<" "<<PointsXYZ[i](1)<<" "<<PointsXYZ[i](2)<<std::endl;
//                                 }



                                 return ptsXYZ;

                                 }
void TriangulationClass::writeFileCeres(std::vector<std::vector<double> > &pts, std::vector<std::vector<double> > &ProjectionMatrices, std::vector<std::vector<int> > &vmat){
    std::ofstream file, filepunti;
//    filepunti.open("Prima.txt");
    file.open("Problem.txt");
    int num_obs=0;
    for(int i=0;i<vmat.size();i++){
        for(int j=0; j<vmat[i].size();j++){
            if(vmat[i][j]==1)
                num_obs++;
        }
    }
    file<<vmat.size()<<" "<<vmat[0].size()<<" "<<num_obs<<std::endl;//header
    for (int i=0;i<vmat[0].size();i++){
        for(int j=0;j<vmat.size();j++){
            if(vmat[j][i]==1){
                file<<j<<" "<<i<<" "<<pts[j][2*i]<<" "<<pts[j][2*i+1]<<std::endl;//camera vs point, aka visibility
            }
        }
    }
    for(int i=0;i<vmat.size();i++){ //camera parameters(ext and intr)
        std::vector<double> R(16),t(3);
        getRandT(ProjectionMatrices[i],R,t);
        transpose(R,4,4);//because cv::Mat need row major convention
        cv::Mat rotvec(3,1,CV_32FC1), Rotation(4,4,CV_64F,R.data());
        cv::Rodrigues(Rotation,rotvec);

        file<<rotvec.at<float>(0,0)<<std::endl<<rotvec.at<float>(1,0)<<std::endl<<rotvec.at<float>(2,0)<<std::endl;//Rot
        //        std::cout<<rotvec<<std::endl; //Ok accesso sia sopra che sotto.
        file<<t[0]<<std::endl<<t[1]<<std::endl<<t[2]<<std::endl;
        //        std::cout<<Translations[i]<<std::endl; //ok accesso
        file<<0.0<<std::endl<<0.0<<std::endl<<0.0<<std::endl<<0.0<<std::endl;
    }
    for(int i=0;i<ptsXYZ.size();i++)
    {
        file<<ptsXYZ[i][0]<<std::endl<<ptsXYZ[i][1]<<std::endl<<ptsXYZ[i][2]<<std::endl; // 3D points
//        filepunti<<PointsXYZ[i](0)<<" "<<PointsXYZ[i](1)<<" "<<PointsXYZ[i](2)<<std::endl;//Ok accesso sopra.
    }
}


