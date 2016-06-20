#ifndef _LRS_H
#define _LRS_H

#include <Eigen/Core>
#include <Eigen/SVD>
#include <vector>

#include "type.h"

using Eigen::MatrixXd;
using namespace std;

inline int count_larger_than(const Vector& v, double value);

void lrs(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& D,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& E,
        double lambda);

int cal_lrs(const vector<vector<double> >& M, double lambda, char* out_prefix);

#endif 
