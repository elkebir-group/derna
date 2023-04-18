//
// Created by xinyu on 4/10/2022.
//

#include "Nussinov.h"
#include "NussinovAlgorithm.h"
#include "utils.h"
#include <iostream>
#include <tuple>
#include <chrono>
#include <cmath>
#include <vector>

using namespace std;


Nussinov::Nussinov(vector<int> protein,int n, int g):n(n),g(g),protein(protein) {
    value.resize(9*n*n*36, 0.0);
//    n_codon = {
//            4, 6, 2, 2, 2, 2, 2, 4, 2, 1, 3, 6, 2, 2, 4, 6, 4, 1, 2, 4
//    };
    codon_index.resize(9*n*n, make_tuple(0,0,0));
    codon_index_x.resize(9*6*n*n, make_tuple(0,0,0));
    codon_index_y.resize(9*6*n*n, make_tuple(0,0,0));
}

Nussinov::~Nussinov() = default;

//void Nussinov::initialize() {
//    codon_index.clear();
//    codon_index_x.clear();
//    codon_index_y.clear();
//    codon_index.resize(9*n*n, make_tuple(0,0,0));
//    codon_index_x.resize(9*6*n*n, make_tuple(0,0,0));
//    codon_index_y.resize(9*6*n*n, make_tuple(0,0,0));
//}

tuple<double, string, string> Nussinov::nussinov(ostream &fout) {
    auto start = chrono::high_resolution_clock::now();
    fout << "Nussinov" << endl;

    double ans = 0;
    double ret = 0;
    double local_max = 0;
//    vector<tuple<int,int,double>> codon_index(9*n*n, make_tuple(0,0,0));
//    vector<tuple<int,int,double>> codon_index_x(9*6*n*n, make_tuple(0,0,0));
//    vector<tuple<int,int,double>> codon_index_y(9*6*n*n, make_tuple(0,0,0));

    // when a == b
    for (int a = 0,b; a < n; a++) {
        b=a;
        const int n_codon_a = n_codon[protein[a]];
        for (int l = 0; l < 3; l++) {
            for (int i = 0 , j; i+l < 3; i++) {
                j = i + l;
                if (sigma(b, j) - sigma(a, i) <= g) continue;
                for (int x = 0; x < n_codon_a; x++) {
                    ret = 0;
                    if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][x][j])) {
                        // F [a][b][i + 1][j − 1][x][x] + 1
                        if (i <= 1 && j >= 1 && i < j - 1) {
                            ret = max(ret, Access(a, b,i + 1, j - 1, x, x) + 1);
                        }

                    }

                    //F [a][b][i + 1][j][x][x]
                    if (i <= 1 && i + 1 <= j) {
                        ret = max(ret, Access(a, b, i + 1, j, x, x));
                    }
                    //F [a][b][i][j − 1][x][x]
                    if (j >= 1 && i <= j - 1) {
                        ret = max(ret, Access(a, b, i, j - 1, x, x));
                    }

                    if (ret > get<2>(codon_index[m(n, a, b,i,j)])) {
                        codon_index[m(n,a,b,i,j)] = make_tuple(x,x,ret);
                    }
                    if (ret > get<2>(codon_index_x[m(n, a, b,i,j,x)])) {
                        codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                    }
                    if (ret > get<2>(codon_index_y[m(n, a, b,i,j,x)])) {
                        codon_index_y[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                    }
                    Access(a, b, i, j, x, x) = ret;
                }
//                cout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",index: " << (9*n*(b-a)+(3*j-2*i)) << endl;
            }
        }
    }


    for (int len = 1; len < n; len++) {
        for(int a = 0,b; a+len < n; a++){
            b=a+len;
            const int n_codon_a = n_codon[protein[a]];
            const int n_codon_b = n_codon[protein[b]];
            for (int i = 2; i >= 0; i--) {
                for (int j = 0; j < 3; j++) {
                    if (sigma(b,j) - sigma(a,i) <= g) continue;
                    for (int x = 0; x < n_codon_a; x++) {
                        for (int y = 0; y < n_codon_b; y++) {
                            ret = 0;
                            // i, j pair
                            if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][y][j])) {
                                //F [a][b][i + 1][j − 1][x][y] + 1
                                if (i <= 1 && j >= 1) {
                                    ret = max(ret, Access(a, b, i + 1,j - 1,x,y) + 1);
                                }

                                if (i == 2 && j >=1) {
                                    // 1 + maxx′ ∈S(va+1){1 + F [a + 1][b][1][j − 1][x′][y]}
                                    if (a + 1 < b) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            local_max = max(local_max, Access(a+1, b, 0, j - 1,p,y));
                                        }
                                        ret = max(ret, local_max + 1);
                                    }
                                    // 1 + F [a + 1][b][1][j − 1][y][y]
                                    if (a + 1 == b) {
                                        ret = max(ret, Access(a+1, b, 0, j-1, y, y) + 1);
                                    }
                                }

                                if (i <= 1 && j == 0) {
                                    // 1 + maxy′∈S(vb−1){F [a][b − 1][i + 1][3][x][y′]}
                                    if (a < b - 1) {
                                        const int n_codon_bp = n_codon[protein[b-1]];
                                        local_max = 0;
                                        for (int q = 0; q < n_codon_bp; q++) {
                                            local_max = max(local_max, Access(a, b-1, i + 1, 2, x, q));
                                        }
                                        ret = max(ret, local_max + 1);
                                    }
                                    // 1 + F [a][b − 1][i + 1][3][x][x]
                                    if (a == b - 1) {
                                        ret = max(ret, Access(a, b-1, i+1, 2, x, x) + 1);
                                    }
                                }

                                if (i == 2 && j == 0) {
                                    // 1 + maxx′ ∈S(va+1),y′ ∈S(vb−1){F [a + 1][b − 1][1][3][x′][y′]}
                                    if (a + 1 < b - 1) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        const int n_codon_bp = n_codon[protein[b-1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            for (int q = 0; q < n_codon_bp; q++) {
                                                local_max = max(local_max, Access(a+1, b-1, 0,2,p,q));
                                            }
                                        }
                                        ret = max(ret, local_max + 1);
                                    }
                                    // 1 + maxx′ ∈S(va+1){F [a + 1][b − 1][1][3][x′][x′]}
                                    if (a + 1 == b - 1) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            local_max = max(local_max, Access(a+1, b-1, 0,2,p,p));
                                        }
                                        ret = max(ret, local_max + 1);
                                    }
                                    // 1
                                    if (sigma(b,j) - 1 == sigma(a,i) + g) { //
                                        ret = max(ret, (double)1.0);
                                    }
                                }

                            }
                            // F [a][b][i + 1][j][x][y]
                            if (i <= 1) {
                                ret = max(ret, Access(a, b, i + 1,j,x,y));
                            }
                            // F [a][b][i][j − 1][x][y]
                            if (j >= 1) {
                                ret = max(ret, Access(a, b, i,j - 1,x,y));
                            }

                            if (i == 2) {
                                // maxx′ ∈S(va+1) F [a + 1][b][1][j][x′][y]
                                if (a + 1 < b) {
                                    const int n_codon_an = n_codon[protein[a+1]];
                                    for (int p = 0; p < n_codon_an; p++) {
                                        ret = max(ret, Access(a+1, b, 0, j,p,y));
                                    }
                                }
                                // F [a + 1][b][1][j][y][y]
                                if (a + 1 == b) {
                                    ret = max(ret, Access(a+1, b, 0, j, y, y));
                                }
                            }

                            if (j == 0) {
                                // maxy′ ∈S(vb−1 ) F [a][b − 1][i][3][x][y′]
                                if (a < b - 1) {
                                    const int n_codon_bp = n_codon[protein[b-1]];
                                    for (int q = 0; q < n_codon_bp; q++) {
                                        ret = max(ret, Access(a, b - 1, i, 2, x, q));
                                    }
                                }
                                // F [a][b − 1][i][3][x][x]
                                if (a == b - 1) {
                                    ret = max(ret, Access(a, b - 1, i, 2, x, x));
                                }
                            }

                            // bifurication

                            // maxa<c<b−1,y′ ∈S(vc ),x′∈S(vc+1){F [a][c][i][3][x][y′] + F [c + 1][b][1][j][x′][y]}
                            if (a < b - 2) {
                                for (int c = a + 1; c < b-1; c++) {
                                    const int n_codon_cn = n_codon[protein[c+1]];
                                    const int n_codon_c = n_codon[protein[c]];
                                    for (int p = 0; p < n_codon_cn; p++) {
                                        for (int q = 0; q < n_codon_c; q++) {
                                            ret = max(ret, Access(a,c, i,2,x,q) + Access(c + 1, b, 0,j,p,y));
                                        }
                                    }
                                }
                            }
                            // maxx′ ∈S(va+1){F [a][a][i][3][x][x] + F [a + 1][b][1][j][x′][y]}
                            if (a < b-1 && i <= 2) {
                                const int n_codon_an = n_codon[protein[a+1]];
                                for (int p = 0; p < n_codon_an; p++) {
                                    ret = max(ret, Access(a, a, i, 2, x, x) + Access(a + 1, b, 0, j, p, y));
                                }
                            }
                            // maxy′∈S(vb−1 ){F [a][b − 1][i][3][x][y′] + F [b][b][1][j][y][y]}
                            if (a < b - 1 && j >= 0) {
                                const int n_codon_bp = n_codon[protein[b-1]];
                                for (int q = 0; q < n_codon_bp; q++) {
                                    ret = max(ret, Access(a, b - 1, i, 2, x, q) + Access(b, b, 0, j, y, y));
                                }
                            }
                            // F [a][a][i][3][x][x] + F [b][b][1][j][y][y]
                            if (a == b - 1 && i <= 2 && j >= 0) {
                                ret = max(ret, Access(a, a, i, 2, x, x) + Access(b, b, 0, j, y, y));
                            }
                            // maxk<3,a<c<b,x′ ∈S(vc ){F [a][c][i][k][x][x′] + F [c][b][k + 1][j][x′][y]}
                            if (a < b - 1) {
                                for (int c = a+1; c < b; c++) {
                                    for (int k = 0; k < 2; k++) {
                                        const int n_codon_c = n_codon[protein[c]];
                                        for (int q = 0; q < n_codon_c; q++) {
                                            ret = max(ret, Access(a, c, i,k,x,q) + Access(c, b, k + 1,j,q,y));
                                        }
                                    }
                                }
                            }
                            // F [a][a][1][2][x][x] + F [a][b][3][j][x][y]
                            if (a < b && i <= 1) {
                                for (int k = i; k < 2; k++) {
                                    ret = max(ret, Access(a, a, i, k, x, x) + Access(a, b, k+1, j, x, y));
                                }
                            }
                            // F [a][b][i][1][x][y] + F [b][b][2][3][y][y]
                            if (a < b && j >= 1) {
                                for (int k = 0; k <= j - 1; k++) {
                                    ret = max(ret, Access(a, b, i, k, x, y) + Access(b, b, k+1, j, y, y));
                                }
                            }

                            if (ret > get<2>(codon_index[m(n, a, b,i,j)])) {
                                codon_index[m(n,a,b,i,j)] = make_tuple(x,y,ret);
                            }
                            if (ret > get<2>(codon_index_x[m(n, a, b,i,j,x)])) {
                                codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,y,ret);
                            }
                            if (ret > get<2>(codon_index_y[m(n, a, b,i,j,y)])) {
                                codon_index_y[m(n,a,b,i,j,y)] = make_tuple(x,y,ret);
                            }
                            Access(a, b, i, j, x, y) = ret;

                            ans = max(ans, ret);
                        }
                    }

//                    fout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",ret: " << ret << endl;
                }
            }

        }
    }


    // find rna and base pair
    string bp = find_bp(0, n-1, 0, 2, get<0>(codon_index[m(n,0,n-1,0,2)]), get<1>(codon_index[m(n,0,n-1,0,2)]), 0);
    string rna = find_bp(0, n-1, 0, 2, get<0>(codon_index[m(n,0,n-1,0,2)]), get<1>(codon_index[m(n,0,n-1,0,2)]), 1);


    fout << "corresponding base pairs: " << bp << ",length: " << bp.size() << endl;
//    fout << "corresponding base pairs (Nussinov Validation): " << nuss_bp << ",length: " << nuss_bp.size() << endl;
    fout << "corresponding rna: " << rna << ",size: " << rna.size() << endl;
//    fout << "corresponding rna2: " << rna2 << ",size: " << rna.size() << endl;
    fout << "max base pair: " << ans << endl;
//    fout << "max base pair: (Nussinov Validation): " << nuss_ans << endl;
    auto end = chrono::high_resolution_clock::now();
    long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
    fout << "Time taken by program is : " << time_take;
    fout << "sec" << endl;

    return make_tuple(Access(0, n-1, 0, 2, get<0>(codon_index[m(n,0,n-1,0,2)]), get<1>(codon_index[m(n,0,n-1,0,2)])), rna, bp); //,0,2
}

string Nussinov::find_bp(int a, int b, int i, int j, int x, int y, int seq) {
    if (a > b || (a == b && i > j)) return "";
    char left = to_char[nucleotides[protein[a]][x][i]];
    char right = to_char[nucleotides[protein[b]][y][j]];
    if (a == b && i == j) {
        switch (seq) {
            case 0:
                return ".";
            case 1:
                string s;
                return s + left;
        }
    }
//    if (sigma(b,j) - sigma(a,i) <= g) return ".";


//    cout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << endl;
    if (a == b) {
        if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][x][j])) {
            if (i <= 1 && j >= 1 && i < j - 1) { //sigma(b, j-1) - sigma(a, i+1) >= g,
                if (Access(a,b,i,j,x,y) == Access(a,b,i+1,j-1,x,y) + 1) {
//                    cout << 1 << endl;
                    switch(seq) {
                        case 0:
                            return "(" + find_bp(a,b,i+1,j-1,x,y,seq) + ")";
                        case 1:
                            return left + find_bp(a,b,i+1,j-1,x,y,seq) + right;
                    }

                }
            }
        }
        if (i <= 1 && i + 1 <= j) { //sigma(b,j) - sigma(a,i+1) >= g,
            if (Access(a,b,i,j,x,y) == Access(a,b,i+1,j,x,y)) {
//                cout << 2 << endl;
                switch(seq) {
                    case 0:
                        return "." + find_bp(a,b,i+1,j,x,y,seq);
                    case 1:
                        return left + find_bp(a,b,i+1,j,x,y,seq);
                }

            }
        }
        if (j >= 1 && i <= j-1) {
            if (Access(a,b,i,j,x,y) == Access(a,b,i,j-1,x,y)) {
//                cout << 3 << endl;
                switch(seq) {
                    case 0:
                        return find_bp(a,b,i,j-1,x,y,seq) + ".";
                    case 1:
                        return find_bp(a,b,i,j-1,x,y,seq) + right;
                }

            }
        }
    }

    if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][y][j])) {

        if (i <= 1 && j >= 1 && a < b) {
            if (Access(a,b,i,j,x,y) == Access(a,b,i+1,j-1,x,y) + 1) {
//                cout << 4 << endl;
                switch(seq) {
                    case 0:
                        return "(" + find_bp(a,b,i+1,j-1,x,y,seq) + ")";
                    case 1:
                        return left + find_bp(a,b,i+1,j-1,x,y,seq) + right;
                }

            }
        }
        if (i == 2 && j >=1 && a + 1 <= b) {
            int p = get<0>(codon_index_y[m(n,a+1,b,0,j-1,y)]);
            if (Access(a,b,i,j,x,y) == Access(a + 1,b,0,j-1,p,y) + 1) { //0,j-1
//                cout << 5 << endl;
                switch(seq) {
                    case 0:
                        return "(" + find_bp(a+1,b,0,j-1,p,y,seq) +")";
                    case 1:
                        return left + find_bp(a+1,b,0,j-1,p,y,seq) + right;
                }

            }
        }

        if (i <= 1 && j == 0 && a <= b-1) {
            int q = get<1>(codon_index_x[m(n,a,b-1,i+1,2,x)]);
            if (Access(a,b,i,j,x,y) == Access(a,b-1,i+1,2,x,q) + 1) { //i+1,2
//                cout << 6 << endl;
                switch(seq) {
                    case 0:
                        return "(" + find_bp(a,b-1,i+1,2,x,q,seq) + ")";
                    case 1:
                        return left + find_bp(a,b-1,i+1,2,x,q,seq) + right;
                }

            }
        }

        if (i == 2 && j == 0 && a+1 <= b-1) {
            int p = get<0>(codon_index[m(n,a+1,b-1,0,2)]);
            int q = get<1>(codon_index[m(n,a+1,b-1,0,2)]);
            if (Access(a,b,i,j,x,y) == Access(a+1,b-1,0,2,p,q) + 1) { //0,2
//                cout << 7 << endl;
                switch(seq) {
                    case 0:
                        return "(" + find_bp(a+1,b-1,0,2,p,q,seq) + ")";
                    case 1:
                        return left + find_bp(a+1,b-1,0,2,p,q,seq) + right;
                }

            }
        }

        if (i == 2 && j == 0 && sigma(b,j) - 1 == sigma(a,i) + g) {
//            cout << 8 << endl;
            if (Access(a,b,i,j,x,y) == 1) {
                switch (seq) {
                    case 0:
                        return "()";
                    case 1:
                        string s;
                        return s + left + right;
                }
            }


        }

    }

    if (i <= 1 && a < b) {
        if (Access(a,b,i,j,x,y) == Access(a,b,i+1,j,x,y)) {
//            cout << 9 << endl;
            switch (seq) {
                case 0:
                    return "." + find_bp(a,b,i+1,j,x,y,seq);
                case 1:
                    return left + find_bp(a,b,i+1,j,x,y,seq);
            }

        }
    }

    if (j >= 1 && a < b) {
        if (Access(a,b,i,j,x,y) == Access(a,b,i,j-1,x,y)) {
//            cout << 10 << endl;
            switch(seq) {
                case 0:
                    return find_bp(a,b,i,j-1,x,y,seq) + ".";
                case 1:
                    return find_bp(a,b,i,j-1,x,y,seq) + right;
            }

        }
    }

    if (i == 2 && a + 1 <= b) {
//        cout << 11 << endl;
        int p = get<0>(codon_index_y[m(n,a+1,b,0,j,y)]);
        if (Access(a,b,i,j,x,y) == Access(a+1,b,0,j,p,y)) { //0,j
//            cout << 11 << endl;
            switch (seq) {
                case 0:
                    return "." + find_bp(a + 1,b,0,j,p,y,seq);
                case 1:
                    return left + find_bp(a + 1,b,0,j,p,y,seq);
            }

        }
    }

    if (j == 0 && a <= b - 1) { //sigma(b-1,2) - sigma(a, i) >= g
        int q = get<1>(codon_index_x[m(n,a,b-1,i,2,x)]);
        if (Access(a,b,i,j,x,y) == Access(a,b-1,i,2,x,q)) { //i,2
//            cout << 12 << endl;
            switch (seq) {
                case 0:
                    return find_bp(a,b - 1,i,2,x,q,seq) + ".";
                case 1:
                    return find_bp(a,b - 1,i,2,x,q,seq) + right;
            }

        }
    }

    if (a < b - 2) { //sigma(b-2,j) - sigma(a,i) >= g
        for (int c = a + 1; c < b-1; c++) {
            int p = get<1>(codon_index_x[m(n,a,c,i,2,x)]);
            int q = get<0>(codon_index_y[m(n,c+1,b,0,j,y)]);
            if (Access(a,b,i,j,x,y) == Access(a,c,i,2,x,p) + Access(c+1,b,0,j,q,y)) { //,i,2 //,0,j
//                cout << 13 << endl;
                return find_bp(a,c,i,2,x,p,seq) + find_bp(c+1,b,0,j,q,y,seq);
            }
        }
    }

    if (a < b-1 && i <= 2) { //sigma(a,2) - sigma(a,i) >= g
        int p = get<0>(codon_index_y[m(n,a+1,b,0,j,y)]);
        if (Access(a,b,i,j,x,y) == Access(a,a,i,2,x,x) + Access(a+1,b,0,j,p,y)) { //,0,j
//            cout << 14 << endl;
            return find_bp(a,a,i,2,x,x,seq) + find_bp(a+1,b,0,j,p,y,seq);
        }
    }

    if (a < b - 1 && j >= 0) { //sigma(b,j) - sigma(b,0) >= g
        int q = get<1>(codon_index_x[m(n,a,b-1,i,2,x)]);
        if (Access(a,b,i,j,x,y) == Access(a,b-1,i,2,x,q) + Access(b,b,0,j,y,y)) { //,i,2
//            cout << 15 << endl;
            return find_bp(a,b-1,i,2,x,q,seq) + find_bp(b,b,0,j,y,y,seq);
        }
    }

    if (a == b - 1 && i <= 2 && j >= 0) { //sigma(a,i) >= g && sigma(b,j) - sigma(b,0) >= g
        if (Access(a,b,i,j,x,y) == Access(a,a,i,2,x,x) + Access(b,b,0,j,y,y)) {
//            cout << 16 << endl;
            return find_bp(a,a,i,2,x,x,seq) + find_bp(b,b,0,j,y,y,seq);
        }
    }

    if (a < b - 1) { //sigma(b-1,j) - sigma(a,i) >= g
        for (int c = a+1; c < b; c++) {
            int n_codon_c = n_codon[protein[c]];
            for (int k = 0; k < 2; k++) {
                for (int t = 0; t < n_codon_c; t++) {
                    if (Access(a,b,i,j,x,y) == Access(a,c,i,k,x,t) + Access(c,b,k+1,j,t,y)) { //,i,k //,k+1,j
//                    cout << 17 << endl;
                        return find_bp(a,c,i,k,x,t,seq) + find_bp(c,b,k+1,j,t,y,seq);
                    }
                }
            }
        }
    }

    if (a < b && i <= 1) {
        for (int k = i; k < 2; k++) {
            if (Access(a,b,i,j,x,y) == Access(a, a, i, k, x, x) + Access(a, b, k+1, j, x, y)) {
//                cout << 18 << endl;
                return find_bp(a,a,i,k,x,x,seq) + find_bp(a,b,k+1,j,x,y,seq);
            }
        }
    }

    if (a < b && j >= 1) {
        for (int k = 0; k <= j - 1; k++) {
            if (Access(a,b,i,j,x,y) == Access(a, b, i, k, x, y) + Access(b, b, k+1, j, y, y)) {
//                cout << 19 << endl;
                return find_bp(a,b,i,k,x,y,seq) + find_bp(b,b,k+1,j,y,y,seq);
            }
        }
    }

    cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << nucleotides[protein[a]][x][i] << " " << nucleotides[protein[b]][y][j] << endl;
    cout << Access(a,b,i,j,x,y) << endl;
    if (Access(a,b,i,j,x,y) == 0) return ".";
    return "FAIL";
}

void Nussinov::lambda_swipe(double incr, ostream &fout, string outfile) {
    int size = ceil(1/incr + 1);
    cout << "size: " << size << endl;
    vector<double> lambda_buffer(size, 0);
    vector<double> F_buffer(size, 0);
    vector<double> CAI_buffer(size, 0);
    vector<pair<string, vector<double>>> dataset(3);
    int index = 0;
    double O;
    double CAI = 0;
    tuple<double, string> temp;

    string rna;
    double lambda = 0;
    while (lambda < 1.001) {
        fout << "lambda: " << lambda << endl;
        temp = nussinov_CAI(lambda, fout);
        O = get<0>(temp);
        rna = get<1>(temp);
        vector<int> rna_n(rna.size());
        transform2num(rna_n, rna);
        CAI = getCAI_s(rna_n, protein);
        if (lambda == 1) CAI = 0;
        lambda_buffer[index] = lambda;
        NussinovAlgorithm F = NussinovAlgorithm(rna_n, rna_n.size(), g);
        cout << lambda << " " << O << " " << CAI << " " << F.nussinov(0, rna_n.size()-1) << " " << lambda*(F.nussinov(0, rna_n.size()-1))+(1-lambda)*CAI << endl;
        F_buffer[index] = F.nussinov(0, rna_n.size()-1);
        if (lambda == 0) F_buffer[index] = 0;
        CAI_buffer[index] = CAI;
        index++;
        lambda += incr;
    }
    dataset[0] = make_pair("lambda", lambda_buffer);
    dataset[1] = make_pair("F", F_buffer);
    dataset[2] = make_pair("CAI", CAI_buffer);
    write_csv(outfile, dataset);
    cout << "swipe done" << endl;
}

tuple<double, string> Nussinov::nussinov_CAI(double lambda, ostream &fout) {
    auto start = chrono::high_resolution_clock::now();
    fout << "Nussinov CAI" << endl;

    double ans = 0;
    double ret = 0;
    double local_max = 0;
//    vector<tuple<int,int,double>> codon_index(9*n*n, make_tuple(0,0,0));
//    vector<tuple<int,int,double>> codon_index_x(9*6*n*n, make_tuple(0,0,0));
//    vector<tuple<int,int,double>> codon_index_y(9*6*n*n, make_tuple(0,0,0));
//    initialize();
    // when a == b
    for (int a = 0,b; a < n; a++) {
        b=a;
        const int p_a = protein[a];
        const int n_codon_a = n_codon[protein[a]];
        for (int l = 0; l < 3; l++) {
            for (int i = 0 , j; i+l < 3; i++) {
                j = i + l;

                for (int x = 0; x < n_codon_a; x++) {
                    ret = 0;
                    if (sigma(b, j) - sigma(a, i) <= g) {
                        if (a == b && i == j) {
                            ret = max(ret, (1 - lambda) * codon_cai_s[p_a][x]);
                        }
                        if (i <= 1 && j >= 1 && i < j - 1) {
                            ret = max(ret, Access(a, b,i + 1, j - 1, x, x));
                        }
                        if (i <= 1 && i + 1 <= j) {
                            ret = max(ret, Access(a, b, i + 1, j, x, x));
                        }
                        //F [a][b][i][j − 1][x][x]
                        if (j >= 1 && i <= j - 1) {
                            ret = max(ret, Access(a, b, i, j - 1, x, x));
                        }

                        Access(a, b, i, j, x, x) = ret;
//                        debug << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << x << ",ret: " << ret << endl;

                        if (ret > get<2>(codon_index[m(n,a,b,i,j)])) codon_index[m(n,a,b,i,j)] = make_tuple(x,x,ret);
                        if (ret > get<2>(codon_index_x[m(n,a,b,i,j,x)])) codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                        if (ret > get<2>(codon_index_y[m(n,a,b,i,j,x)])) codon_index_y[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                        continue;
                    }
                    if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][x][j])) {
                        // F [a][b][i + 1][j − 1][x][x] + 1
                        if (i <= 1 && j >= 1 && i < j - 1) {
                            ret = max(ret, Access(a, b,i + 1, j - 1, x, x) + lambda);
                        }
                    }

                    //F [a][b][i + 1][j][x][x]
                    if (i <= 1 && i + 1 <= j) {
                        ret = max(ret, Access(a, b, i + 1, j, x, x));
                    }
                    //F [a][b][i][j − 1][x][x]
                    if (j >= 1 && i <= j - 1) {
                        ret = max(ret, Access(a, b, i, j - 1, x, x));
                    }

                    if (ret > get<2>(codon_index[m(n, a, b,i,j)])) {
                        codon_index[m(n,a,b,i,j)] = make_tuple(x,x,ret);
                    }
                    if (ret > get<2>(codon_index_x[m(n, a, b,i,j,x)])) {
                        codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                    }
                    if (ret > get<2>(codon_index_y[m(n, a, b,i,j,x)])) {
                        codon_index_y[m(n,a,b,i,j,x)] = make_tuple(x,x,ret);
                    }
                    Access(a, b, i, j, x, x) = ret;
//                    debug << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << x << ",ret: " << ret << endl;

                }

            }
        }
    }


    for (int len = 1; len < n; len++) {
        for(int a = 0,b; a+len < n; a++){
            b=a+len;
            const int p_a = protein[a];
            const int p_b = protein[b];
            const int n_codon_a = n_codon[protein[a]];
            const int n_codon_b = n_codon[protein[b]];
            for (int i = 2; i >= 0; i--) {
                for (int j = 0; j < 3; j++) {
                    for (int x = 0; x < n_codon_a; x++) {
                        for (int y = 0; y < n_codon_b; y++) {
                            ret = 0;
                            if (sigma(b,j) - sigma(a,i) <= g) {
                                const int n_codon_an = n_codon[protein[a+1]];
                                const int n_codon_bp = n_codon[protein[b-1]];
                                if (i <= 1 && j >= 1) {
                                    ret = max(ret, Access(a,b,i+1,j-1,x,y));
                                }
                                if (i <= 1) {
                                    ret = max(ret, Access(a,b,i+1,j,x,y));
                                }
                                if (j >= 1) {
                                    ret = max(ret, Access(a,b,i,j-1,x,y));
                                }
                                if (i == 2) {
                                    if (a+1 < b) {
                                        for (int p = 0; p < n_codon_an; p++) {
                                            ret = max(ret, Access(a+1, b, 0, j, p, y) + (1-lambda)* codon_cai_s[p_a][x]);
                                        }
                                    }
                                    if (a+1 == b) {
                                        ret = max(ret, Access(a+1, b, 0, j, y, y) + (1-lambda)* codon_cai_s[p_a][x]);
                                    }
                                }
                                if (j == 0) {
                                    if (a < b-1) {
                                        for (int q = 0; q < n_codon_bp; q++) {
                                            ret = max(ret, Access(a,b-1,i,2,x,q) + (1-lambda)* codon_cai_s[p_b][y]);
                                        }
                                    }
                                    if (a == b-1) {
                                        ret = max(ret, Access(a, b-1, i, 2, x, x) + (1-lambda)* codon_cai_s[p_b][y]);
                                    }
                                }
                                if (i == 2 && j == 0) {
                                    if (a+1 < b-1) {
                                        for (int p = 0; p < n_codon_an; p++) {
                                            for (int q = 0; q < n_codon_bp; q++) {
                                                ret = max(ret, Access(a+1,b-1,0,2,p,q) + (1-lambda)* (codon_cai_s[p_a][x] + codon_cai_s[p_b][y]));
                                            }
                                        }
                                    }
                                    if (a+1 == b-1) {
                                        for (int p = 0; p < n_codon_an; p++) {
                                            ret = max(ret, Access(a+1, b-1, 0,2,p,p) + (1-lambda)* (codon_cai_s[p_a][x] + codon_cai_s[p_b][y]));
                                        }
                                    }
                                }

                                if (ret > get<2>(codon_index[m(n,a,b,i,j)])) {
                                    codon_index[m(n,a,b,i,j)] = make_tuple(x,y,ret);
                                }
                                if (ret > get<2>(codon_index_x[m(n,a,b,i,j,x)])) {
                                    codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,y,ret);
                                }
                                if (ret > get<2>(codon_index_y[m(n,a,b,i,j,y)])) {
                                    codon_index_y[m(n,a,b,i,j,y)] = make_tuple(x,y,ret);
                                }
                                Access(a, b, i, j, x, y) = ret;
//                                debug << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << ",ret: " << ret << endl;
                                continue;
                            }

                            // i, j pair
                            if (complementary(nucleotides[p_a][x][i], nucleotides[p_b][y][j])) {
                                //F [a][b][i + 1][j − 1][x][y] + 1
                                if (i <= 1 && j >= 1) {
//                                    debug << 1 << " " << ret << endl;
                                    ret = max(ret, Access(a, b, i + 1,j - 1,x,y) + lambda);
                                }

                                if (i == 2 && j >=1) {
                                    // 1 + maxx′ ∈S(va+1){1 + F [a + 1][b][1][j − 1][x′][y]}
                                    if (a + 1 < b) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            local_max = max(local_max, Access(a+1, b, 0, j - 1,p,y));
                                        }
                                        ret = max(ret, local_max + lambda + (1-lambda) * codon_cai_s[p_a][x]);
//                                        debug << 2 << " " << ret << endl;
                                    }
                                    // 1 + F [a + 1][b][1][j − 1][y][y]
                                    if (a + 1 == b) {
                                        ret = max(ret, Access(a+1, b, 0, j-1, y, y) + lambda + (1 - lambda) * codon_cai_s[p_a][x]);
//                                        debug << 3 << " " << ret << endl;
                                    }
                                }

                                if (i <= 1 && j == 0) {
                                    // 1 + maxy′∈S(vb−1){F [a][b − 1][i + 1][3][x][y′]}
                                    if (a < b - 1) {
                                        const int n_codon_bp = n_codon[protein[b-1]];
//                                        debug << "n_codon_bp: " << n_codon_bp << endl;
                                        local_max = 0;
                                        for (int q = 0; q < n_codon_bp; q++) {
//                                            debug << "local_max: " << local_max << "x: " << x << "q: " << q << " " << O.Access(a, b-1, i + 1, 2, x, q) << endl;
                                            local_max = max(local_max, Access(a, b-1, i + 1, 2, x, q));
                                        }
                                        ret = max(ret, local_max + lambda + (1 - lambda) * codon_cai_s[p_b][y]);
//                                        debug << 4 << "local_max: " << local_max << " " << lambda << " " << (1 - lambda) * CAI(p_b, y) << " " <<  ret << endl;
//                                        debug << "b: " << b << "p_b: " << p_b << ",y: " << y << endl;
                                    }
                                    // 1 + F [a][b − 1][i + 1][3][x][x]
                                    if (a == b - 1) {
                                        ret = max(ret, Access(a, b-1, i+1, 2, x, x) + lambda + (1 - lambda) * codon_cai_s[p_b][y]);
//                                        debug << 5 << " " << ret << endl;
                                    }
                                }

                                if (i == 2 && j == 0) {
                                    // 1 + maxx′ ∈S(va+1),y′ ∈S(vb−1){F [a + 1][b − 1][1][3][x′][y′]}
                                    if (a + 1 < b - 1) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        const int n_codon_bp = n_codon[protein[b-1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            for (int q = 0; q < n_codon_bp; q++) {
                                                local_max = max(local_max, Access(a+1, b-1, 0,2,p,q));
                                            }
                                        }
                                        ret = max(ret, local_max + lambda + (1-lambda) * (codon_cai_s[p_a][x]+ codon_cai_s[p_b][y]));
//                                        debug << 6 << " " << ret << endl;
                                    }
                                    // 1 + maxx′ ∈S(va+1){F [a + 1][b − 1][1][3][x′][x′]}
                                    if (a + 1 == b - 1) {
                                        const int n_codon_an = n_codon[protein[a+1]];
                                        local_max = 0;
                                        for (int p = 0; p < n_codon_an; p++) {
                                            local_max = max(local_max, Access(a+1, b-1, 0,2,p,p));
                                        }
                                        ret = max(ret, local_max + lambda + (1 - lambda) * (codon_cai_s[p_a][x]+ codon_cai_s[p_b][y]));
//                                        debug << 7 << " " << ret << endl;
                                    }
                                    if (sigma(b,j) - 1 == sigma(a,i) + g) {
                                        ret = max(ret, lambda + (1-lambda) * (codon_cai_s[p_a][x] + codon_cai_s[p_b][y]));
//                                        debug << 8 << " " << ret << endl;
                                    }
                                }

                            }
                            // F [a][b][i + 1][j][x][y]
                            if (i <= 1) {
                                ret = max(ret, Access(a, b, i + 1,j,x,y));
//                                debug << 9 << " " << ret << endl;
                            }
                            // F [a][b][i][j − 1][x][y]
                            if (j >= 1) {
                                ret = max(ret, Access(a, b, i,j - 1,x,y));
//                                debug << 10 << " " << ret << endl;
                            }

                            if (i == 2) {
                                // maxx′ ∈S(va+1) F [a + 1][b][1][j][x′][y]
                                if (a + 1 < b) {
                                    const int n_codon_an = n_codon[protein[a+1]];
                                    local_max = 0;
                                    for (int p = 0; p < n_codon_an; p++) {
                                        local_max = max(local_max, Access(a+1, b, 0, j,p,y));
                                    }
                                    ret = max(ret, local_max + (1 - lambda) * codon_cai_s[p_a][x]);
//                                    debug << 11 << " " << ret << endl;
                                }
                                // F [a + 1][b][1][j][y][y]
                                if (a + 1 == b) {
                                    ret = max(ret, Access(a+1, b, 0, j, y, y) + (1-lambda) * codon_cai_s[p_a][x]);
//                                    debug << 12 << " " << ret << endl;
                                }

                            }

                            if (j == 0) { // sigma(b-1, 2) - sigma(a, i) >= g
                                // maxy′ ∈S(vb−1 ) F [a][b − 1][i][3][x][y′]
                                if (a < b - 1) {
                                    const int n_codon_bp = n_codon[protein[b-1]];
                                    local_max = 0;
                                    for (int q = 0; q < n_codon_bp; q++) {
                                        local_max = max(local_max, Access(a, b - 1, i, 2, x, q));
                                    }
                                    ret = max(ret, local_max + (1-lambda) * codon_cai_s[p_b][y]);
//                                    debug << 13 << " " << ret << endl;
                                }
                                // F [a][b − 1][i][3][x][x]
                                if (a == b - 1) {
                                    ret = max(ret, Access(a, b - 1, i, 2, x, x) + (1-lambda)*codon_cai_s[p_b][y]);
//                                    debug << 14 << " " << ret << endl;
                                }
                            }

                            // bifurication

                            // maxa<c<b−1,y′ ∈S(vc ),x′∈S(vc+1){F [a][c][i][3][x][y′] + F [c + 1][b][1][j][x′][y]}
                            if (a < b - 2) {
                                for (int c = a + 1; c < b-1; c++) {
                                    const int n_codon_cn = n_codon[protein[c+1]];
                                    const int n_codon_c = n_codon[protein[c]];
                                    for (int p = 0; p < n_codon_cn; p++) {
                                        for (int q = 0; q < n_codon_c; q++) {
                                            ret = max(ret, Access(a,c, i,2,x,q) + Access(c + 1, b, 0,j,p,y));
                                        }
                                    }
                                }
//                                debug << 15 << " " << ret << endl;
                            }
                            // maxx′ ∈S(va+1){F [a][a][i][3][x][x] + F [a + 1][b][1][j][x′][y]}
                            if (a < b-1 && i <= 2) {
                                const int n_codon_an = n_codon[protein[a+1]];
                                for (int p = 0; p < n_codon_an; p++) {
                                    ret = max(ret, Access(a, a, i, 2, x, x) + Access(a + 1, b, 0, j, p, y));
                                }
//                                debug << 16 << " " << ret << endl;
                            }
                            // maxy′∈S(vb−1 ){F [a][b − 1][i][3][x][y′] + F [b][b][1][j][y][y]}
                            if (a < b - 1 && j >= 0) {
                                const int n_codon_bp = n_codon[protein[b-1]];
                                for (int q = 0; q < n_codon_bp; q++) {
                                    ret = max(ret, Access(a, b - 1, i, 2, x, q) + Access(b, b, 0, j, y, y));
                                }
//                                debug << 17 << " " << ret << endl;
                            }
                            // F [a][a][i][3][x][x] + F [b][b][1][j][y][y]
                            if (a == b - 1 && i <= 2 && j >= 0) {
                                ret = max(ret, Access(a, a, i, 2, x, x) + Access(b, b, 0, j, y, y));
//                                debug << 18 << " " << ret << endl;
                            }
                            // maxk<3,a<c<b,x′ ∈S(vc ){F [a][c][i][k][x][x′] + F [c][b][k + 1][j][x′][y]}
                            if (a < b - 1) {
                                for (int c = a+1; c < b; c++) {
                                    const int n_codon_c = n_codon[protein[c]];
                                    for (int k = 0; k < 2; k++) {
                                        for (int q = 0; q < n_codon_c; q++) {
                                            ret = max(ret, Access(a, c, i,k,x,q) + Access(c, b, k + 1,j,q,y) - (1-lambda) * codon_cai_s[protein[c]][q]);
                                        }
                                    }
                                }
//                                debug << 19 << " " << ret << endl;
                            }
                            // F [a][a][1][2][x][x] + F [a][b][3][j][x][y]
                            if (a < b && i <= 1) {
                                for (int k = i; k < 2; k++) {
                                    ret = max(ret, Access(a, a, i, k, x, x) + Access(a, b, k+1, j, x, y) - (1-lambda) * codon_cai_s[p_a][x]);
                                }
//                                debug << 20 << " " << ret << endl;
                            }
                            // F [a][b][i][1][x][y] + F [b][b][2][3][y][y]
                            if (a < b && j >= 1) { //sigma(b,j) - sigma(a,i) >= g
                                for (int k = 0; k <= j - 1; k++) {
//                                    if (sigma(b,j) - sigma(b, k+1) < g) continue;
                                    ret = max(ret, Access(a, b, i, k, x, y) + Access(b, b, k+1, j, y, y)- (1-lambda) * codon_cai_s[p_b][y]);
                                }
//                                debug << 21 << " " << ret << endl;
                            }


                            if (ret > get<2>(codon_index[m(n, a, b,i,j)])) {
                                codon_index[m(n,a,b,i,j)] = make_tuple(x,y,ret);
                            }
                            if (ret > get<2>(codon_index_x[m(n, a, b,i,j,x)])) {
                                codon_index_x[m(n,a,b,i,j,x)] = make_tuple(x,y,ret);
                            }
                            if (ret > get<2>(codon_index_y[m(n, a, b,i,j,y)])) {
                                codon_index_y[m(n,a,b,i,j,y)] = make_tuple(x,y,ret);
                            }
                            Access(a, b, i, j, x, y) = ret;

                            ans = max(ans, ret);
//                            debug << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << ",ret: " << ret << endl;
                        }
                    }
                }
            }

        }
    }


    // find rna
    string rna = nussinov_CAI_tb(0, n-1, 0, 2, get<0>(codon_index[m(n,0,n-1,0,2)]), get<1>(codon_index[m(n,0,n-1,0,2)]),lambda);


    fout << "corresponding rna: " << rna << ",size: " << rna.size() << endl;
    fout << "max base pair + CAI: " << ans << endl;
    auto end = chrono::high_resolution_clock::now();
    long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
    fout << "Time taken by program is : " << time_take;
    fout << "sec" << endl;
    return make_tuple(Access(0,n-1,0,2,get<0>(codon_index[m(n,0,n-1,0,2)]), get<1>(codon_index[m(n,0,n-1,0,2)])), rna); //,0,2
}

string Nussinov::nussinov_CAI_tb(int a, int b, int i, int j, int x, int y, double lambda) {
    if (a > b || (a == b && i > j)) {
        return "";
    }

    char left = to_char[nucleotides[protein[a]][x][i]];
    char right = to_char[nucleotides[protein[b]][y][j]];
//    cout << "a: " << a << ",b: " << b << ",i: " << i << ",j: " << j << ",x: " << x << ",y: " << y << " " << left << " " << right << endl;
    if (a == b && i == j) {
        string s;
        return s + left;
    }
    int p_a = protein[a];
    int p_b = protein[b];

    if (a == b) {
        if (sigma(b, j) - sigma(a, i) <= g) {
            if (i <= 1 && j >= 1 && i < j - 1) {
                if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j-1,x,y))) {
                    return left + nussinov_CAI_tb(a,b,i+1,j-1,x,y,lambda) + right;
                }
            }
            if (i <= 1 && i + 1 <= j) {
                if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j,x,y))) {
                    return left + nussinov_CAI_tb(a,b,i+1,j,x,y,lambda);
                }
            }
            if (j >= 1 && i <= j-1) {
                if (compare(Access(a,b,i,j,x,y), Access(a,b,i,j-1,x,y))) {
                    return nussinov_CAI_tb(a,b,i,j-1,x,y,lambda) + right;
                }
            }
        }
        if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][x][j])) {
            if (i <= 1 && j >= 1 && i < j - 1) { //sigma(b, j-1) - sigma(a, i+1) >= g,
                if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j-1,x,y) + lambda)) {
                    return left + nussinov_CAI_tb(a,b,i+1,j-1,x,y,lambda) + right;
                }
            }
        }
        if (i <= 1 && i + 1 <= j) { //sigma(b,j) - sigma(a,i+1) >= g,
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j,x,y))) {
                return left + nussinov_CAI_tb(a,b,i+1,j,x,y,lambda);
            }
        }
        if (j >= 1 && i <= j-1) {
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i,j-1,x,y))) {
                return nussinov_CAI_tb(a,b,i,j-1,x,y,lambda) + right;
            }
        }
    }

    if (sigma(b,j) - sigma(a,i) <= g) {
        if (i <= 1 && j >= 1) {
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j-1,x,y))) {
                return left + nussinov_CAI_tb(a, b, i+1, j-1,x,y, lambda) + right;
            }
        }
        if (i <= 1) {
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j,x,y))) {
                return left + nussinov_CAI_tb(a, b, i+1, j,x,y,  lambda);
            }
        }
        if (j >= 1) {
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i,j-1,x,y))) {
                return nussinov_CAI_tb( a, b, i, j-1,x,y,  lambda) + right;
            }
        }
        if (i==2 && j==0 && a+1 <= b-1) {
            int p = get<0>(codon_index[m(n,a+1,b-1,0,2)]);
            int q = get<1>(codon_index[m(n,a+1,b-1,0,2)]);
            if (compare(Access(a, b, i, j, x, y), Access(a+1,b-1,0,2,p,q) + (1-lambda)* (codon_cai_s[p_a][x] + codon_cai_s[p_b][y]))) { //get<0>(codon_index[m(n,a+1,b-1,0,2)]) //get<1>(codon_index[m(n,a+1,b-1,0,2)])
                return left + nussinov_CAI_tb(a+1, b-1, i, j,p,q,  lambda) + right;
            }

        }
        if (i == 2 && a + 1 <= b) {
            int p = get<0>(codon_index_y[m(n,a+1,b,0,j,y)]);
            if (compare(Access(a, b, i, j, x, y), Access(a+1, b, 0, j, p, y) + (1-lambda)* codon_cai_s[p_a][x])) { //,0,j
                return left + nussinov_CAI_tb(a+1, b, 0, j,p,y, lambda);
            }
        }

        if (j == 0 && a <= b-1) {
            int q = get<1>(codon_index_x[m(n,a,b-1,i,2,x)]);
            if (compare(Access(a, b, i, j, x, y), Access(a,b-1,i,2,x,q) + (1-lambda)*codon_cai_s[p_b][y])) { //,i,2
                return nussinov_CAI_tb(a, b-1, i, 2, x,q, lambda) + right;
            }
        }
    }

    if (complementary(nucleotides[protein[a]][x][i], nucleotides[protein[b]][y][j])) {

        if (i <= 1 && j >= 1 && a < b) {
            if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j-1,x,y) + lambda)) {
                return left + nussinov_CAI_tb(a,b,i+1,j-1,x,y,lambda) + right;
            }
        }
        if (i == 2 && j >=1 && a + 1 <= b) {
            int p = get<0>(codon_index_y[m(n,a+1,b,0,j-1,y)]);
            if (compare(Access(a,b,i,j,x,y), Access(a + 1,b,0,j-1,p,y) + lambda + (1-lambda)*codon_cai_s[p_a][x])) { //get<0>(codon_index_y[m(n,a+1,b,0,j-1,y)])
                return left + nussinov_CAI_tb(a+1,b,0,j-1,p,y,lambda) + right;
            }

        }

        if (i <= 1 && j == 0 && a <= b-1) {
            int q  = get<1>(codon_index_x[m(n,a,b-1,i+1,2,x)]);
            if (compare(Access(a,b,i,j,x,y), Access(a,b-1,i+1,2,x,q) + lambda + (1-lambda)*codon_cai_s[p_b][y])) { //get<1>(codon_index_x[m(n,a,b-1,i+1,2,x)])
                return left + nussinov_CAI_tb(a,b-1,i+1,2,x,q,lambda) + right;
            }

        }

        if (i == 2 && j == 0 && a+1 <= b-1) {
            int p = get<0>(codon_index[m(n,a+1,b-1,0,2)]);
            int q = get<1>(codon_index[m(n,a+1,b-1,0,2)]);
            if (compare(Access(a,b,i,j,x,y), Access(a+1,b-1,0,2,p,q) + lambda + (1-lambda)*(codon_cai_s[p_a][x]+codon_cai_s[p_b][y]))) { //,0,2
                return left + nussinov_CAI_tb(a+1,b-1,0,2,p,q,lambda) + right;
            }
        }

        if (i == 2 && j == 0 && sigma(b,j) - 1 == sigma(a,i) + g) {
            if (compare(Access(a,b,i,j,x,y), lambda + (1-lambda)*(codon_cai_s[p_a][x]+ codon_cai_s[p_b][y]))) {
                string s;
                return s + left + right;
            }
        }

    }

    if (i <= 1 && a < b) {
        if (compare(Access(a,b,i,j,x,y), Access(a,b,i+1,j,x,y))) {
            return left + nussinov_CAI_tb(a,b,i+1,j,x,y,lambda);
        }
    }

    if (j >= 1 && a < b) {
        if (compare(Access(a,b,i,j,x,y), Access(a,b,i,j-1,x,y))) {
            return nussinov_CAI_tb(a,b,i,j-1,x,y,lambda) + right;
        }
    }

    if (i == 2 && a + 1 <= b) {
        int p = get<0>(codon_index_y[m(n,a+1,b,0,j,y)]);
        if (compare(Access(a,b,i,j,x,y), (1-lambda)*codon_cai_s[p_a][x] + Access(a+1,b,0,j,p,y))) { //get<0>(codon_index_y[m(n,a+1,b,0,j,y)])
            return left + nussinov_CAI_tb(a + 1,b,0,j,p,y,lambda);
        }
    }

    if (j == 0 && a <= b - 1) { //sigma(b-1,2) - sigma(a, i) >= g
        int q = get<1>(codon_index_x[m(n,a,b-1,i,2,x)]);
        if (compare(Access(a,b,i,j,x,y), (1-lambda)* codon_cai_s[p_b][y] + Access(a,b-1,i,2,x,q))) { //get<1>(codon_index_x[m(n,a,b-1,i,2,x)])
            return nussinov_CAI_tb(a,b - 1,i,2,x,q,lambda) + right;
        }
    }

    if (a < b - 2) { //sigma(b-2,j) - sigma(a,i) >= g
        for (int c = a + 1; c < b-1; c++) {
            int p = get<1>(codon_index_x[m(n,a,c,i,2,x)]);
            int q = get<0>(codon_index_y[m(n,c+1,b,0,j,y)]);
            if (compare(Access(a,b,i,j,x,y), Access(a,c,i,2,x,p) + Access(c+1,b,0,j,q,y))) { //get<1>(codon_index_x[m(n,a,c,i,2,x)]) //get<0>(codon_index_y[m(n,c+1,b,0,j,y)])
                return nussinov_CAI_tb(a,c,i,2,x,p,lambda) + nussinov_CAI_tb(c+1,b,0,j,q,y,lambda);
            }
        }
    }
    if (a < b-1 && i <= 2) { //sigma(a,2) - sigma(a,i) >= g
        int p = get<0>(codon_index_y[m(n,a+1,b,0,j,y)]);
        if (compare(Access(a,b,i,j,x,y), Access(a,a,i,2,x,x) + Access(a+1,b,0,j,p,y))) { //get<0>(codon_index_y[m(n,a+1,b,0,j,y)])
            return nussinov_CAI_tb(a,a,i,2,x,x,lambda) + nussinov_CAI_tb(a+1,b,0,j,p,y,lambda);
        }
    }

    if (a < b - 1 && j >= 0) { //sigma(b,j) - sigma(b,0) >= g
        int q = get<1>(codon_index_x[m(n,a,b-1,i,2,x)]);
        if (compare(Access(a,b,i,j,x,y), Access(a,b-1,i,2,x,q) + Access(b,b,0,j,y,y))) { //get<1>(codon_index_x[m(n,a,b-1,i,2,x)])
            return nussinov_CAI_tb(a,b-1,i,2,x,q,lambda) + nussinov_CAI_tb(b,b,0,j,y,y,lambda);
        }
    }

    if (a == b - 1 && i <= 2 && j >= 0) {
        if (compare(Access(a,b,i,j,x,y), Access(a,a,i,2,x,x) + Access(b,b,0,j,y,y))) {
            return nussinov_CAI_tb(a,a,i,2,x,x,lambda) + nussinov_CAI_tb(b,b,0,j,y,y,lambda);
        }
    }


    if (a < b - 1) { //sigma(b-1,j) - sigma(a,i) >= g
        for (int c = a+1; c < b; c++) {
            int n_codon_c = n_codon[protein[c]];
            for (int k = 0; k < 2; k++) {
                for (int t = 0; t < n_codon_c; t++) {
                    if (compare(Access(a,b,i,j,x,y), (lambda-1)*codon_cai_s[protein[c]][t] + Access(a,c,i,k,x,t) + Access(c,b,k+1,j,t,y))) {
                        return nussinov_CAI_tb(a,c,i,k,x,t,lambda) + nussinov_CAI_tb(c,b,k+1,j,t,y,lambda);
                    }
                }
            }
        }
    }

    if (a < b && i <= 1) {
        for (int k = i; k < 2; k++) {
            if (compare(Access(a,b,i,j,x,y), (lambda-1)* codon_cai_s[p_a][x] + Access(a, a, i, k, x, x) + Access(a, b, k+1, j, x, y))) {
                return nussinov_CAI_tb(a,a,i,k,x,x,lambda) + nussinov_CAI_tb(a,b,k+1,j,x,y,lambda);
            }
        }
    }

    if (a < b && j >= 1) {
        for (int k = 0; k <= j - 1; k++) {
            if (compare(Access(a,b,i,j,x,y), (lambda-1)* codon_cai_s[p_b][y] + Access(a, b, i, k, x, y) + Access(b, b, k+1, j, y, y))) {
                return nussinov_CAI_tb(a,b,i,k,x,y,lambda) + nussinov_CAI_tb(b,b,k+1,j,y,y,lambda);
            }
        }
    }

    cout << a << " " << b << " " << i << " " << j << " " << x << " " << get<1>(codon_index[m(n,a,b,i,j)]) << " " << nucleotides[protein[a]][x][i] << " " << nucleotides[protein[b]][y][j] << endl;
    cout << Access(a,b,i,j,x,y) << endl;
    if (Access(a,b,i,j,x,y) == 0) return ".";
    return "FAIL";
}


Nussinov::Nussinov(const Nussinov & Copy):n(Copy.n),g(Copy.g),value(Copy.value.size()) {
    cout << "is it here" << endl;
//    copy(Copy.value.begin(), Copy.value.end(), value.begin());
}