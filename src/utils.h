
//
// Created by xinyu on 3/23/2022.
//

#ifndef RNA_DESIGN_UTILS_H
#define RNA_DESIGN_UTILS_H
#include <string>
#include <vector>
#include "Nussinov.h"
#include "default.h"
using namespace std;

int aa_index (char aa);

int to_int(char a);

void char2num(vector<int> &, string &);

int n_index(char n);

int n_index2(int n1, int n2);

bool complementary(int X, int Y);

int sigma(int a, int i);

//double CAI(int a, int x);

//double maxCAI(int p);

double getCAI(const vector<int> &, const vector<int> &);

double stand_getCAI(const vector<int> &, const vector<int> &);

int m(int n, int a, int b, int i, int j, int x=-1); //

void write_csv(string, const vector<pair<string, vector<double>>> &);

double getCAI_s(const vector<int> &, const vector<int> &);

double stand_getCAI_s(const vector<int> &, const vector<int> &);

bool compare(double x, double y, double epislon=0.00001);

bool greaterThan(double a, double b);

bool basepair(int X, int Y);

string num2String(vector<int> &);

void transform2num(vector<int>&, string);

int get_index(vector<string>& seqs, string seq);

//int interior_loop(int i, int j, int h, int k, int i1, int j1, int h1, int k1, int n1, int n2);

const int getxPos(int, vector<int> &);

//int auterminal(int a, int b);

bool add_auterminal(int a, int b);

bool add_ggmm(int a, int b);

bool add_uugamm(int a, int b);

int l(int a, int i, int b, int j);

vector<int> read_fasta(string, ostream&);

vector<int> read_rna(string input, ostream& fout);

double evaluate_CAI(string &,vector<int> &,int);

double evaluate_CAI(vector<int> &, vector<int> &);

double evaluate_CAI_N(string &,vector<int> &,int);

double evaluate_MFE(string &, vector<int> &);

double evaluate_MFE(vector<int> &, vector<int> &);

int evaluate_BP_N(string &, vector<int> &, int);

void usage(void);
void help(void);


#endif //RNA_DESIGN_UTILS_H
