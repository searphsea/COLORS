#ifndef _TYPE_H
#define _TYPE_H
#include <Eigen/Core>
#include <Eigen/SVD>
//

using Eigen::MatrixXd;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> Array;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
#endif
