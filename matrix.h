#ifndef _MATRIX_H
#define _MATRIX_H

#define ALPHA 21

#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/SVD>

#include "type.h"

using Eigen::MatrixXd;
using namespace std;

int cal_matrix(const vector<string>& msa, const vector<double>& weight,
        vector<vector<double> >& MI, vector<vector<double> >&OMES, vector<vector<double> >& COV);

double cal_mi(double pab[ALPHA][ALPHA], double pa[], double pb[]);
double cal_cov(double pab[ALPHA][ALPHA], double pa[], double pb[]);
double cal_omes(double pab[ALPHA][ALPHA], double pa[], double pb[]);

int cal_apc(const vector<vector<double> >&m, vector<vector<double> >&apc);

int load_matrix(MatrixXd & dir, const vector<vector<double> >& source);
int load_matrix(vector<vector<double> >& dir, const Matrix& source);
int output_matrix(char *out_path, const vector<vector<double> >& m);
int output_score(char *out_path, 
        const vector<vector<double> >& m,
        const vector<vector<double> >& apc, 
        const vector<vector<double> >& sparse);

#endif
