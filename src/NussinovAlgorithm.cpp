//
// Created by xinyu on 4/24/2022.
//

#include "NussinovAlgorithm.h"
#include "utils.h"
#include <iostream>
using namespace std;

NussinovAlgorithm::NussinovAlgorithm(vector<int>& rna, int n, int g):rna(rna),n(n),g(g) {
    dp.resize(n*n, -1);
}

NussinovAlgorithm::~NussinovAlgorithm() {

}

int NussinovAlgorithm::nussinov(int i, int j) {
    if (dp[index(i,j)] != -1) return dp[index(i,j)];

    if (i+g >= j) {
        dp[index(i,j)] = 0;
        return 0;
    }

    int ret = 0;
    if (i+g < j) ret = max(ret, nussinov(i+1, j));
    if (i+g < j) ret = max(ret, nussinov(i, j-1));
//    cout << rna[i] << " " << rna[j] << endl;
    if (complementary(rna[i],rna[j])) {

        if (i+g < j) ret = max(ret, nussinov(i+1, j-1) + 1);

        if (i+g == j-1) ret = max(ret, 1);
    }
    if (i+g < j) {
        for (int k=i+g; k<j; k++) {
            ret = max(ret, nussinov(i,k) + nussinov(k+1, j));
        }
    }

    dp[index(i,j)] = ret;
    return ret;
}

string NussinovAlgorithm::get_bp(int i, int j) {
    if (i > j) return "";
    if (i == j) return ".";

    if (nussinov(i,j) == nussinov(i+1,j)) return "." + get_bp(i+1, j);
    if (nussinov(i,j) == nussinov(i,j-1)) return get_bp(i,j-1) + ".";

    if (complementary(rna[i], rna[j]) && nussinov(i,j) == nussinov(i+1, j-1) + 1) {
        return "(" + get_bp(i+1, j-1) + ")";
    }

    for (int k=i+1; k<j; k++) {
        if (nussinov(i,j) == nussinov(i,k) + nussinov(k+1,j)) {
            return get_bp(i,k) + get_bp(k+1, j);
        }
    }

    return "FAIL";

}

int NussinovAlgorithm::get_nbp(int i, int j) {
    string bp = get_bp(i,j);
    int count = 0;
    for (int k = 0; k < (int)bp.size(); k++) {
        if (bp[k] == '(') count++;
    }
    return count;
}

int NussinovAlgorithm::index(int i, int j) {
    return n*i+j;
}

NussinovAlgorithm::NussinovAlgorithm(const NussinovAlgorithm & Copy):rna(Copy.rna),n(Copy.n),dp(Copy.dp.size()) {
    copy(Copy.dp.begin(), Copy.dp.end(), dp.begin());
}