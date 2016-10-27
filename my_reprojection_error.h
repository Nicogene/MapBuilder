
#ifndef CERES_MY_REPROJECTION_ERROR_H_
#define CERES_MY_REPROJECTION_ERROR_H_

#include "ceres/rotation.h"


// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 10 parameters: 3 for rotation, 3 for translation, 2 for
// focal length and 2 for optical centre.
namespace ceres{
    struct MyReprojectionError {
      MyReprojectionError(double observed_x, double observed_y)
          : observed_x(observed_x), observed_y(observed_y) {}

      template <typename T>
      bool operator()(const T* const camera,
                      const T* const point,
                      T* residuals) const {
        // camera[0,1,2] are the angle-axis rotation.
        T p[3];
        ceres::AngleAxisRotatePoint(camera, point, p);

        // camera[3,4,5] are the translation.
        p[0] -= camera[3];
        p[1] -= camera[4];
        p[2] -= camera[5];

        // Compute the center of distortion. The sign change comes from
        // the camera model that Noah My's Bundler assumes, whereby
        // the camera coordinate system has a negative z axis.
        const T& focalx = camera[6];
        const T& focaly = camera[7];
        const T& ccx = camera[8];
        const T& ccy = camera[9];
        T xp =  p[0] / p[2];//ho tolto il meno, assumeva che la z fosse negativa
        T yp =  p[1] / p[2];

        // Compute final projected point position.
        T predicted_x = focalx * xp - ccx;
        T predicted_y = focaly * yp - ccy;

        // The error is the difference between the predicted and observed position.
        residuals[0] = T(observed_x) - predicted_x;
        residuals[1] = T(observed_y) - predicted_y  ;

        return true;
      }

      // Factory to hide the construction of the CostFunction object from
      // the client code.
      static ceres::CostFunction* Create(const double observed_x,
                                         const double observed_y) {
        return (new ceres::AutoDiffCostFunction<MyReprojectionError, 2/*dim residual*/, 10 /*dim of camera*/, 3/*dim of the point, 3D point*/>(
                    new MyReprojectionError(observed_x, observed_y)));
      }

      double observed_x;
      double observed_y;
    };

    // Templated pinhole camera model for used with Ceres.  The camera is
    // parameterized using 11 parameters. 4 for rotation, 3 for
    // translation, 2 for focal length and 2 for the optical center.
    struct MyReprojectionErrorWithQuaternions {
      // (u, v): the position of the observation with respect to the image
      // center point.
      MyReprojectionErrorWithQuaternions(double observed_x, double observed_y)
          : observed_x(observed_x), observed_y(observed_y) {}

      template <typename T>
      bool operator()(const T* const camera_rotation,
                      const T* const camera_translation_and_intrinsics,
                      const T* const point,
                      T* residuals) const {
        const T& focalx = camera_translation_and_intrinsics[3];
        const T& focaly = camera_translation_and_intrinsics[4];
        const T& ccx = camera_translation_and_intrinsics[5];
        const T& ccy = camera_translation_and_intrinsics[6];


        // Use a quaternion rotation that doesn't assume the quaternion is
        // normalized, since one of the ways to run the bundler is to let Ceres
        // optimize all 4 quaternion parameters unconstrained.
        T p[3];
        QuaternionRotatePoint(camera_rotation, point, p);

        p[0] += camera_translation_and_intrinsics[0];
        p[1] += camera_translation_and_intrinsics[1];
        p[2] += camera_translation_and_intrinsics[2];

        // Compute the center of distortion. The sign change comes from
        // the camera model that Noah My's Bundler assumes, whereby
        // the camera coordinate system has a negative z axis.
        T xp =  p[0] / p[2]; //ho tolto il meno, assumeva che la z fosse negativa
        T yp =  p[1] / p[2];

        // Compute final projected point position.
        T predicted_x = focalx * xp - ccx;
        T predicted_y = focaly * yp - ccy;

        // The error is the difference between the predicted and observed position.
        residuals[0] = predicted_x - T(observed_x);
        residuals[1] = predicted_y - T(observed_y);

        return true;
      }

      // Factory to hide the construction of the CostFunction object from
      // the client code.
      static ceres::CostFunction* Create(const double observed_x,
                                         const double observed_y) {
        return (new ceres::AutoDiffCostFunction<
                MyReprojectionErrorWithQuaternions, 2, 4, 7, 3>(
                    new MyReprojectionErrorWithQuaternions(observed_x,
                                                                observed_y)));
      }

      double observed_x;
      double observed_y;
    };
}


#endif  
