#include "matrix.h"
//
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>

int cal_matrix(const vector<string>& msa, const vector<double>& weight,
        vector<vector<double> >& MI, vector<vector<double> >&OMES, vector<vector<double> >& COV){

//assert(msa.size() > 0 && msa[0].size() > 0);
    int M = msa.size();
    int N = msa[0].size();
    MI.assign(N, vector<double>(N, 0.0));
    OMES.assign(N, vector<double>(N, 0.0));
    COV.assign(N, vector<double>(N, 0.0));

    int pseudoc = 1;     

    double* pa = new double[N*ALPHA];
    double pab[ALPHA][ALPHA];

    double neff = 0.0;
    for(int k = 0; k < M; k++){
        neff += weight[k];
    }

    for (int i = 0; i < N; i++)
    {
        double *pi =  pa + i * ALPHA;
        for (int aa = 0; aa < ALPHA; aa++){
            pi[aa] = pseudoc;
        }
        
        for(int k = 0; k < M; k++){
            pi[msa[k][i]] += weight[k];
        }

        for (int aa = 0; aa < ALPHA; aa++){
            pi[aa] /= pseudoc * ALPHA * 1.0 + neff;
        }
    }

    double *pi, *pj;
    for(int i = 0; i < N; i++){
        pi = pa + i * ALPHA;
        for(int j = i; j < N; j++){
            pj = pa + j * ALPHA;
            for(int a = 0; a < ALPHA; a++){
                for(int b = 0; b < ALPHA; b++){
                    if (i == j){
                        pab[a][b] = (a == b) ? pi[a] : 0.0;
                    }
                    else{
                        pab[a][b] = pseudoc * 1.0 / ALPHA;
                    }
                }
            }
            if(i != j){
                for(int k = 0; k < M; k++){
                    int a = msa[k][i];
                    int b = msa[k][j];
                    pab[a][b] += weight[k];
                }
                for (int a = 0; a < ALPHA; a++){
                    for (int b = 0; b < ALPHA; b++){
                        pab[a][b] /= pseudoc * ALPHA * 1.0 + neff;
                    }
                }
            }
            double mi = cal_mi(pab, pi, pj);
            double omes = cal_omes(pab, pi, pj);
            double cov = cal_cov(pab, pi, pj);
            MI[i][j] = MI[j][i] = mi;
            OMES[i][j] = OMES[j][i] = omes;
            COV[i][j] = COV[j][i] = cov;
        }
    }
    delete pa;

    return 0;
}

double cal_mi(double pab[ALPHA][ALPHA], double pa[], double pb[]){
    double sum = 0.0;
    for(int i = 0; i < ALPHA; i++){
        for(int j = 0; j < ALPHA; j++){
            if (pab[i][j] >  0.0){
                sum += pab[i][j] * log(pab[i][j]/pa[i]/pb[j]);
            }
        }
    }
    return sum;
}

double cal_cov(double pab[ALPHA][ALPHA], double pa[], double pb[]){
    double sum = 0.0;
    for(int i = 0; i < ALPHA; i++){
        for(int j = 0; j < ALPHA; j++){
            if(pab[i][j] > 0.0){
                sum += fabs(pab[i][j] - pa[i]*pb[j]);
            }
        }
    }
    return sum;
}

double cal_omes(double pab[ALPHA][ALPHA], double pa[], double pb[]){
    double sum = 0.0;
    for(int i = 0; i < ALPHA; i++){
        for(int j = 0; j < ALPHA; j++){
            double expected = pa[i] * pb[j];
            if(expected == 0.0){
                continue;
            }
            sum += (pab[i][j] - expected)*(pab[i][j] - expected);
        }
    }
    return sum;
}

int cal_apc(const vector<vector<double> >&m, vector<vector<double> >&apc){
    vector<double> mean(m.size(), 0.0);
    double sum = 0.0;
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m.size(); j++){
            if(i == j){
                continue;
            }
            mean[i] += m[i][j];
            sum += m[i][j];
        }
        mean[i] /= m.size() - 1.0;
    }

    sum /= m.size() * (m.size() - 1.0);

    for(int i = 0; i < m.size(); i++){
        for(int j = i+1; j < m.size(); j++){
            apc[i][j] = apc[j][i] = m[i][j] - mean[i] * mean[j] / sum;
        }
    }

    return 0;
}

int load_matrix(vector<vector<double> >& dir, const Matrix& source){
    for(int i = 0; i < dir.size(); i++){
        for(int j = 0; j < dir[i].size(); j++){
            dir[i][j] = source(i, j);
        }
    }
    return 0;
}

int load_matrix(MatrixXd& dir, const vector<vector<double> >& source){
    for(int i = 0; i < source.size(); i++){
        for(int j = 0; j < source[i].size(); j++){
            dir(i,j) = source[i][j];
        }
    }
    return 0;
}

int output_matrix(char *out_path, const vector<vector<double> >& m){
    ofstream fout(out_path);
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m.size(); j++){
            fout << m[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
    return 0;
}

int output_score(char *out_path, 
        const vector<vector<double> >& m,
        const vector<vector<double> >& apc, 
        const vector<vector<double> >& sparse){
    ofstream fout(out_path);
    fout << "#pos_i pos_j original_score apc_score lrs_score\n";
    for(int i = 0; i < m.size(); i++){
        for(int j = i+1; j < m.size(); j++){
            fout << i+1 << " " << j+1 << " "  
                << m[i][j] <<  " "<< apc[i][j] << " " << sparse[i][j] << "\n";
        }
    }
    fout.close();
    return 0;
}
