#ifndef _SEQUENCE_H
#define _SEQUENCE_H
//
#define SEQ_BUFFER_SIZE 10240

#include <cstdio>
#include <string>
#include <vector>

using namespace std;

int read_msa(const char*file_path, vector<string>& msa);
int cal_seq_weight(const vector<string>& msa, vector<double>& weight, double seq_id);
double cal_seq_sim(const string& seq1, const string& seq2);
unsigned char aatoi(unsigned char aa);


#endif
