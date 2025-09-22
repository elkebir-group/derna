
//
// Created by xinyu on 3/23/2022.
//

#ifndef RNA_DESIGN_UTILS_H
#define RNA_DESIGN_UTILS_H
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include "Nussinov.h"
#include "default.h"
using namespace std;

int aa_index (char aa);

int to_int(char a);

vector<string> filterCandidates(
        const vector<string>& candidates,
        const unordered_map<int, char>& letter_map);

void char2num(vector<int> &, string &);

int n_index(char n);

int n_index2(int n1, int n2);

bool complementary(int X, int Y);

int sigma(int a, int i);

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

int get_index(vector<string>& seqs, string & seq);

int getxPos(int, vector<int> &);


bool add_auterminal(int a, int b);

bool add_ggmm(int a, int b);

bool add_uugamm(int a, int b);

int l(int a, int i, int b, int j);

vector<int> read_fasta(string &, ostream&);

vector<int> read_rna(string & input);

pair<int, int> find_amino_acid_and_codon_index(const vector<int>& codon);

double evaluate_CAI(string &,vector<int> &,int);

double evaluate_CAI(string &, int type = 0);

double evaluate_CAI(vector<int> &, vector<int> &);

double evaluate_CAI_N(string &,vector<int> &,int);

double evaluate_MFE(string &);

double evaluate_MFE(vector<int> &, string & bp);

int evaluate_BP_N(string &, int);

void usage();
void help();

// Save and load helpers
template<typename T>
void save_vector_binary(const vector<T> &vec, const string &filename) {
    ofstream file(filename, ios::binary);
    size_t size = vec.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(T));
}

template<typename T>
void load_vector_binary(vector<T> &vec, const string &filename) {
    ifstream file(filename, ios::binary);
    size_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(T));
}

template<typename T>
void save_vector_vector_binary(const vector<vector<T>> &vecs, const string &filename) {
    ofstream file(filename, ios::binary);
    size_t outer_size = vecs.size();
    file.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));
    for (const auto& vec : vecs) {
        size_t inner_size = vec.size();
        file.write(reinterpret_cast<const char*>(&inner_size), sizeof(inner_size));
        file.write(reinterpret_cast<const char*>(vec.data()), inner_size * sizeof(T));
    }
}

template<typename T>
void load_vector_vector_binary(vector<vector<T>> &vecs, const string &filename) {
    ifstream file(filename, ios::binary);
    size_t outer_size;
    file.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
    vecs.resize(outer_size);
    for (auto& vec : vecs) {
        size_t inner_size;
        file.read(reinterpret_cast<char*>(&inner_size), sizeof(inner_size));
        vec.resize(inner_size);
        file.read(reinterpret_cast<char*>(vec.data()), inner_size * sizeof(T));
    }
}

#endif //RNA_DESIGN_UTILS_H
