#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

#include "colors.h"
#include "sequence.h"
#include "matrix.h"

using Eigen::MatrixXd;
using namespace std;

int main(int argc, const char** argv) {
    if(argc != 3){
        cout << "usage:<1>msa-file <2>output-prefix\n";
        return -1;
    }
    const char *msa_path = argv[1];
    const char *out_prefix = argv[2];

    vector<string> msa;
    vector<double> weight;
    vector<vector<double> > MI;
    vector<vector<double> > OMES;
    vector<vector<double> > COV;

    double seq_id = 0.8;
    read_msa(msa_path, msa);
    cal_seq_weight(msa, weight, seq_id);
    cal_matrix(msa, weight, MI, OMES, COV);

    char out_mi_path[1024];
    char out_omes_path[1024];
    char out_cov_path[1024];
    sprintf(out_mi_path, "%s.mi",out_prefix);
    sprintf(out_omes_path, "%s.omes",out_prefix);
    sprintf(out_cov_path, "%s.cov",out_prefix);
    
    cal_lrs(MI, 0.13, out_mi_path);
    cal_lrs(OMES, 0.26, out_omes_path);
    cal_lrs(COV, 0.1, out_cov_path);

    return 0;
}


