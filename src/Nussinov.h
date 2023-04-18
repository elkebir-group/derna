//
// Created by xinyu on 4/10/2022.
//

#ifndef RNA_DESIGN_NUSSINOV_H
#define RNA_DESIGN_NUSSINOV_H
#include <vector>
#include <ostream>

using namespace std;

class Nussinov {
    vector<double> value;
    vector<int> protein;
//    vector<int> n_codon;
//    string rna;
    int n,g;
    vector<tuple<int,int,double>> codon_index;
    vector<tuple<int,int,double>> codon_index_x;
    vector<tuple<int,int,double>> codon_index_y;

public:
    Nussinov(vector<int> protein, int n=0, int g=0);
    Nussinov(const Nussinov &);
    ~Nussinov();
    tuple<double, string, string> nussinov (ostream& fout);
    string get_rna();
    string find_bp(int a, int b, int i, int j, int x, int y, int seq);
    void lambda_swipe(double incr, ostream& fout, string outfile);
    tuple<double, string> nussinov_CAI(double lambda, ostream& fout);
    string nussinov_CAI_tb(int a, int b, int i, int j, int x, int y, double lambda);
    inline double & Access(int a, int b, int i, int j, int x, int y);
private:
    void initialize();
};

inline double & Nussinov::Access(int a, int b, int i, int j, int x, int y) {
//    cout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << endl;
//    if (a == b) {
//        assert(x==y);
//        if (i > j) {
//            throw std::invalid_argument("i > j");
//        }
//        assert(j>=i);
//
//    }
//    cout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << ",index: " << 36*(9*(n-1)*(b-a)+9*(2*b-a)+(3*j-2*i))+6*x+y << ",size: " << value.size() << endl;
//    assert(36*(9*(n-1)*(b-a)+9*(2*b-a)+(3*j-2*i))+6*x+y < value.size());
    return value[36*(9*(n-1)*(b-a)+9*(2*b-a)+(6-3*i+j))+6*x+y];
}

#endif //RNA_DESIGN_NUSSINOV_H
