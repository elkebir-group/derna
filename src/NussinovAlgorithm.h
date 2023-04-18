//
// Created by xinyu on 4/24/2022.
//

#ifndef RNA_DESIGN_NUSSINOVALGORITHM_H
#define RNA_DESIGN_NUSSINOVALGORITHM_H
#include <string>
#include <vector>

using namespace std;

class NussinovAlgorithm {
    vector<int> rna;
    int n,g;
    vector<int> dp;

public:
    NussinovAlgorithm(vector<int> &rna, int n, int g);
    NussinovAlgorithm(const NussinovAlgorithm &);
    ~NussinovAlgorithm();
    int nussinov(int i, int j);
//    inline int nussinov();
    string get_bp(int i, int j);
    int get_nbp(int, int);

private:
    int index(int i, int j);

};
#endif //RNA_DESIGN_NUSSINOVALGORITHM_H
