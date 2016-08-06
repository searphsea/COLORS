#include "sequence.h"
#include <cstdlib>
#include <fstream>
#include <string> 
//
int read_msa(const char*file_path, vector<string>& msa){
    ifstream fin;
    fin.open(file_path);
    if(fin.fail()){
        printf("The path %s does not exist!\n", file_path);
        exit(-1);
    }
    string line;
    while(getline(fin, line, '\n')){
        if(line.size() == 0 || line[0] == '>'){
            continue;
        }
        string seq = line;
        for(int i = 0; i < line.size(); i++){
            seq[i] = aatoi(line[i]);
        }
        msa.push_back(seq);
    }
    fin.close();
}

int cal_seq_weight(const vector<string>& msa, vector<double>& weight, double seq_id){
    int M = msa.size();
    weight.assign(M, 1.0);
    
    for(int i = 0; i < M; i++){
        for(int j = i + 1; j < M; j++){
            double sim = cal_seq_sim(msa[i], msa[j]);
            if(sim >= seq_id){
                weight[i] += 1.0;
                weight[j] += 1.0;
            }
        }
    }

    for(int i = 0; i < M; i++){
        weight[i] = 1.0 / weight[i];
    }
}

double cal_seq_sim(const string& seq1, const string& seq2){
    double sim = 0.0;
    for(int i = 0; i < seq1.size(); i++){
        sim += seq1[i] == seq2[i];
    }
    return sim / seq1.size();
}

unsigned char aatoi(unsigned char aa) {
    char id;
    switch(aa){
        case '-':
            id = 0;
            break;
        case 'A':
            id = 1;
            break;
        case 'C':
            id = 2;
            break;
        case 'D':
            id = 3;
            break;
        case 'E':
            id = 4;
            break;
        case 'F':
            id = 5;
            break;
        case 'G':
            id = 6;
            break;
        case 'H':
            id = 7;
            break;
        case 'I':
            id = 8;
            break;
        case 'K':
            id = 9;
            break;
        case 'L':
            id = 10;
            break;
        case 'M':
            id = 11;
            break;
        case 'N':
            id = 12;
            break;
        case 'P':
            id = 13;
            break;
        case 'Q':
            id = 14;
            break;
        case 'R':
            id = 15;
            break;
        case 'S':
            id = 16;
            break;
        case 'T':
            id = 17;
            break;
        case 'V':
            id = 18;
            break;
        case 'W':
            id = 19;
            break;
        case 'Y':
            id = 20;
            break;
        default:
            id = 0;
    }
    return id;
}
