#include "colors.h"
#include "matrix.h"
//
#include <cmath>
#include <iostream>

inline int count_larger_than(const Vector& v, double value) {
  int count = 0;
  for (int i = 0; i < v.size(); ++i) {
    if (v[i] > value) ++count;
  }
  return count;
}

void lrs(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& D,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& E, double lambda) {

    const int M = D.rows();
    const int N = D.cols();

    Matrix Y = D;
    A = Matrix::Zero(M, N);
    E = Matrix::Zero(M, N);
    Array zero = Matrix::Zero(M, N);

    const double rho = 1.5;

    Eigen::JacobiSVD<Matrix> svd_only_singlar_values(Y);
    const double norm_two =
        svd_only_singlar_values.singularValues()(0);
    const double norm_inf = Y.array().abs().maxCoeff() / lambda;
    const double dual_norm = std::max(norm_two, norm_inf);
    const double d_norm = D.norm();
    Y /= dual_norm;

    double mu = 1.25 / norm_two;
    const double mu_bar = mu * 1.0e+7;

    bool converged = false;
    int max_iter = 1000;
    double error_tolerance = 1.0e-7;
    int iter = 0;
    int total_svd = 0;
    int sv = 10;
    
    while (!converged) {
        Array temp_T = D - A + (1.0 / mu) * Y;
        E = (temp_T - lambda / mu).max(zero) + (temp_T + lambda / mu).min(zero);

        E = E.array().max(zero);

        Eigen::JacobiSVD<Matrix> svd(D - E + 1.0 / mu * Y,
                Eigen::ComputeFullU | Eigen::ComputeFullV);
        Matrix U = svd.matrixU();
        Matrix V = svd.matrixV();
        Vector singularValues = svd.singularValues();

        int svp = count_larger_than(singularValues, 1 / mu);
        if (svp < sv) {
            sv = std::min(svp + 1, N);
        } else {
            sv = std::min(svp + static_cast<int>(0.05 * N + 0.5), N);
        }

        Matrix S_th =
            (singularValues.head(svp).array() - 1.0 / mu).matrix().asDiagonal();
        A = U.leftCols(svp) * S_th * V.leftCols(svp).transpose();

        A = A.array().max(zero);

        total_svd += 1;
        Matrix Z = D - A - E;
        Y = Y + mu * Z;
        mu = std::min(mu * rho, mu_bar);

        double objective = Z.norm() / d_norm;

        if (objective < error_tolerance) {
            converged = true;
        }

        if (++iter >= max_iter) {
            break;
        }
    }
}

int cal_lrs(const vector<vector<double> >& m, double lambda, char*out_prefix){
    int len = m.size();
    
    MatrixXd M = MatrixXd::Zero(len, len);
    MatrixXd L = MatrixXd::Zero(len, len);
    MatrixXd S = MatrixXd::Zero(len, len);
    load_matrix(M, m);

    lrs(M, L, S, lambda);
    
    vector<vector<double > > sparse(len, vector<double>(len, 0.0));
    vector<vector<double > > apc(len, vector<double>(len, 0.0));
    load_matrix(sparse, S);
    cal_apc(m, apc);
    
    char out_path[2048];
    sprintf(out_path, "%s_lrs", out_prefix);
    output_matrix(out_path, sparse);
    sprintf(out_path, "%s_score", out_prefix);
    output_score(out_path, m, apc, sparse);
    
    return 0;
}

