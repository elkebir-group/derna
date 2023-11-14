//
// Created by xinyu on 5/18/2022.
//

#include "Zuker.h"
#include "utils.h"
#include <chrono>
#include <tuple>
#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <queue>
#include "params/constants.h"

using namespace std;



Zuker::Zuker(int n, int mode, vector<int> & protein):protein(protein),n(n) {
    start_index.resize(n*n, 0);
    index_offset.resize(n*n, 0);
    int idx = 0;
    last_idx = 0;
//    cnt = 0;
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            idx = n*a+b;
            start_index[idx] = last_idx;
            index_offset[idx] = n_codon[protein[a]] * n_codon[protein[b]];
            last_idx += index_offset[idx] * 9;
        }
    }

    if (mode == 1) {
        Z.resize(last_idx, 0);
        E.resize(last_idx, inf);
        M.resize(last_idx, inf);
        TM0.resize(last_idx,inf);
        e0 = inf;
    } else {
        O.resize(last_idx, 0);
        E1.resize(last_idx, inf);
        M1.resize(last_idx, inf);
        TM.resize(last_idx, inf);

        Z2.resize(last_idx, 0);
        E2.resize(last_idx, inf);
        M2.resize(last_idx, inf);
        TM2.resize(last_idx,inf);

        Z_CAI.resize(last_idx, 0);
        E_CAI.resize(last_idx, 0);
        M_CAI.resize(last_idx, 0);
        TM_CAI.resize(last_idx,0);
        e = inf;
    }
    OB.resize(last_idx, 0);
    EB.resize(last_idx, 0);
    MB.resize(last_idx, 0);
    O_bt.resize(last_idx);
    E_bt.resize(last_idx);
    M_bt.resize(last_idx);
    TM_bt.resize(last_idx);
    codon_selection.resize(n,-1);
    minX = -1, minY = -1;


    init_values();
}

//((X == 0 && Y == 3) || (X == 3 && Y == 0) || (X == 1 && Y == 2) || (X == 2 && Y == 1) || (X == 2 && Y == 3) || (X == 3 && Y == 2))
void Zuker::init_values() {
    basepair.resize(16,0);
    basepair[(0<<2)+3] = 1;
    basepair[(3<<2)+0] = 1;
    basepair[(1<<2)+2] = 1;
    basepair[(2<<2)+1] = 1;
    basepair[(2<<2)+3] = 1;
    basepair[(3<<2)+2] = 1;
    nucle_seq.resize(3*n,-1);
    sector.resize(3*n);
    bp_bond.resize(3*n);
    ava_nucle_p.resize(3*6*n, 0);
    ava_nucle_m.resize(3*6*n, 0);
    for (int a = 0; a < n; ++a) {
        for (int x = 0; x < n_codon[protein[a]]; ++x) {
            for (int i = 0; i < 3; ++i) {
                int idx = index(a,x,i);
                ava_nucle_p[idx] = ava_nucleotides_int(a,x,i,1);
                ava_nucle_m[idx] = ava_nucleotides_int(a,x,i,-1);
            }
        }
    }
}

int Zuker::calculate_Z(ostream &fout) {

    auto start = chrono::high_resolution_clock::now();
    fout << "Zuker" << endl;

    int ret = inf, energy = inf;
    int t = 0;

    // fill vector E and M first before filling Z
    calculate_E();
    cout << "E done" << endl;


    int nuc_len = 3*n;
    int a = 0, i = 0;
    int la = sigma(a,i);
    int pa = protein[a];
    const int n_codon_a = n_codon[pa];


    for (int len = 1; len < nuc_len; ++len) {

        int j = (i+len) % 3;
        int b = (3*a+i+len) / 3;
        int lb = sigma(b,j);
        if (b == n) continue;

        std::cout << "\rZ -- Completed: " << (((n<<1)-b)*(b-1))/(n*(n-1.)) * 100 << "%" << std::flush;
        int pb = protein[b];
        const int n_codon_b = n_codon[pb];

        for (int x = 0; x < n_codon_a; x++) {
            for (int y = 0; y < n_codon_b; y++) {
                if (a == b && x != y) continue;
                ret = inf, energy = inf, t = 0;

                int idx = index(a, b, i, j, x, y);

                int yj = nucleotides[pb][y][j];
                int xi = nucleotides[pa][x][i];

                // check if current ends can be paired
                int type = BP_pair[xi+1][yj+1];

                // type > 0 denote paired
                if (type > 0) {
                    ret = min(ret,Access_E(idx) + AU[xi][yj]);
                    if (energy > ret) {
                        energy = ret;
                        t = -1;
                    }
                }

                if (j >= 1) {
                    ret = min(ret, Access_Z(a, b, i, j - 1, x, y)); // idx - 36
                    if (energy > ret) {
                        energy = ret;
                        t = -2;
                    }
                }

                if (j == 0) {
                    if (a < b - 1) {
                        int ppb = protein[b-1];
                        const int n_codon_bp = n_codon[ppb];
//                        int idx_t = index(a, b - 1, i, 2, x, 0);
                        for (int q = 0; q < n_codon_bp; q++) {

                            ret = min(ret, Access_Z(a, b - 1, i, 2, x, q)); //  ,idx_t + q
                        }
                    }
                    if (a == b - 1) {
                        ret = min(ret, Access_Z(a, b - 1, i, 2, x, x)); // a, b - 1, i, 2, x, x
                    }
                    if (energy > ret) {
                        energy = ret;
                        t = -3;
                    }
                }



                // bifurication
                for (int lc = la + 1; lc <= lb - 4; lc++) {
                    int c = lc / 3;
                    int i1 = lc % 3;
                    int pc = protein[c];
                    int n_codon_c = n_codon[pc];

                    for (int hx = 0; hx < n_codon_c; ++hx) {
                        if (c == a && x != hx) continue;
                         int hi = nucleotides[pc][hx][i1];
                         if (i1 >= 1) {
                             if (a == c && x != hx) continue;
                             ret = min(ret, Access_Z(a,c,i,i1-1,x,hx) + Access_E(c,b,i1,j,hx,y) + AU[
                                     hi][yj]); ////, idx_1+hx //, idx_2+6*hx
                         } else {
                             if (a == c-1) {
                                 ret = min(ret, Access_Z(a,c-1,i,2,x,x) + Access_E(c,b,i1,j,hx,y) + AU[hi][yj]);
                             } else {
                                 int n_codon_cp = n_codon[protein[c-1]];
                                 for (int ky = 0; ky < n_codon_cp; ++ky) {
                                     ret = min(ret, Access_Z(a,c-1,i,2,x,ky) + Access_E(c,b,i1,j,hx,y) + AU[hi][yj]); //, idx_3+ky //, idx_4+6*hx
                                 }
                             }

                         }
                    }

                }
                if (energy > ret) {
                    energy = ret;
                    t = -4;
                }

                // save minimum value with corresponding codon selection
                if (b == n-1 && ret < e0) {
                    e0 = ret;
                    minX = x, minY = y;
                }

                Access_Z(idx) = ret;
                Access_OB(idx) = t;

            }
        }

    }


    auto end = chrono::high_resolution_clock::now();
    long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
    fout << "Energy: " << Access_Z(0,n-1,0,2,minX, minY)/100.0 << endl;
    fout << "Time taken by DP is : " << time_take;
    fout << "sec" << endl;

    int res = Access_Z(0,n-1,0,2,minX, minY);
    return res;

}

void Zuker::calculate_E() {
    int min_energy = inf, energy = inf;
    int t;
    string an,bp,cp,dn;

    int nuc_len = 3*n;
    vector<int> T(1536,inf);
    vector<int> T1(1536,0);
    vector<int> H(96,inf);
    vector<int> MT(256, inf);

    static bool updatedT[1536];
    static bool updatedH[96];
    static bool updatedM[256];

    memset(updatedH, false, 96);
    memset(updatedM, false, 256);
    memset(updatedT, false, 1536);


    for (int len = 4; len < nuc_len; ++len) {
        std::cout << "\rE/M -- Completed: " << (len*1.0/nuc_len) * 100 << "%" << std::flush;

        int max_a = n - (int)floor(len/3);
        for (int a = 0; a < max_a; ++a) {
            for (int i = 0; i < 3; ++i) {

                int j = (i+len) % 3;
                int b = (3*a+i+len) / 3;

                if (b == n) continue;
                int pa = protein[a];
                int pb = protein[b];
                int pna = protein[a + 1];
                int ppb = protein[b - 1];
                const int n_codon_a = n_codon[pa];
                const int n_codon_b = n_codon[pb];
                const int n_codon_an = n_codon[pna];
                const int n_codon_bp = n_codon[ppb];
                int la = sigma(a,i);
                int lb = sigma(b,j);
                int l = lb - la;

                // set updatedT to false for each entry
                // set updatedH to false for each entry
                // memset
                memset(updatedH, false, 96);
                memset(updatedM, false, 256);
                memset(updatedT, false, 1536); //4096

                for (int x = 0; x < n_codon_a; x++) {
                    for (int y = 0; y < n_codon_b; y++) {
                        int xi = nucleotides[pa][x][i];
                        int yj = nucleotides[pb][y][j];

                        int type = BP_pair[xi+1][yj+1];

                        int xi_ = 0,xi2_ = 0, _2yj = 0, _yj = 0;
                        if (i == 0) {
                            xi_ = nucleotides[pa][x][i+1];
                            xi2_ = nucleotides[pa][x][i+2];
                        }
                        if (j == 2) {
                            _yj = nucleotides[pb][y][j-1];
                            _2yj = nucleotides[pb][y][j-2];
                        }
                        if (i == 1) {
                            xi_ = nucleotides[pa][x][i+1];
                        }
                        if (j == 1) {
                            _yj = nucleotides[pb][y][j-1];
                        }

                        // ends not paired
                        if (!type) {
                            calculate_M(a,b,i,j,x,y);
                            continue;
                        }

                        min_energy = inf, energy = inf, t = 0;


                        int an_int = ava_nucle_p[index(a,x,i)];
                        int bp_int = ava_nucle_m[index(b,y,j)];


                        int idx1 = T_index(type, xi_,xi2_,_2yj, _yj); //
                        if (!updatedT[idx1]) {
                            // hairpin
                            if (l + 1 > 8) {
                                int idx = H_index(type, xi_, _yj); //
                                if (!updatedH[idx]) {
                                    H[idx] = hairpin(la, lb, l, xi, yj, an_int, bp_int);
                                    updatedH[idx] = true;
                                }
                                min_energy = min(min_energy, H[idx]);
                            }
                            else {
                                min_energy = min(min_energy, hairpin(la, lb, l, pa, pb, pna, ppb, n_codon_an, n_codon_bp, xi, yj, i, j, x, y, an_int, bp_int));
                            }
                            if (energy > min_energy) {
                                energy = min_energy;
                                t = -1;
                            }
                            // interior: stacking + bulge + internal;
                            min_energy = min(min_energy, internal(a, b, i, j, x, y, la, lb, xi, yj, an_int, bp_int));
                            if (energy > min_energy) {
                                energy = min_energy;
                                t = EB[index(a,b,i,j,x,y)];
                            }
                            // multiloop
                            min_energy = min(min_energy, multi_loop(a, b, i, j, x, y, pa, pb, n_codon_an, n_codon_bp));
                            if (energy > min_energy) {
                                energy = min_energy;
                                t = -3;
                            }
                            T[idx1] = min_energy;
                            T1[idx1] = t;
                            updatedT[idx1] = true;
                        } else {
                            t = T1[idx1];
                        }
                        min_energy = T[idx1];

                        Access_E(a,b,i,j,x,y) =  min_energy;

                        Access_EB(a,b,i,j,x,y) = t;
                        calculate_M(a,b,i,j,x,y);

                    }
                }
            }
        }

    }

}

int Zuker::hairpin(int la, int lb, int l, int xi, int yj, int an_int, int bp_int) const {
    int hairpin_energy = inf;
    // xi, yj are the paired ends
    for (int xi_ = 0; xi_ < 4; ++xi_) { // xi_ is the nucleotide to the right of xi
        if (((1 << xi_) & an_int) == 0) continue;
        for (int _yj = 0; _yj < 4; ++_yj) { // _yj is the nucleotide to the left of yj
            if (((1 << _yj) & bp_int) == 0) continue;
            // check if the current values of xi_ and _yj are compatible in the same codon
            if (lb - la <= 4 &&  !rightCodon(la+1,lb-1,xi_,_yj)) continue; //(la+1)/3==(lb-1)/3 &&
            hairpin_energy = min(hairpin_energy, hairpin_loop(xi,yj,xi_,_yj,l-1));
        }
    }
    return hairpin_energy;
}

int Zuker::hairpin(int la, int lb, int l, int pa, int pb, int pna, int ppb, int n_codon_an, int n_codon_bp, int xi, int yj, int i, int j, int x, int y, int an_int, int bp_int) const {
    int hairpin_energy = inf;
    string s;
    switch (l) {
        int xi_, _yj, xi2_, _2yj, xi3_, _3yj;

        case 4:
            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                } else {
                    xi2_ = nucleotides[pb][y][0];
                }
                s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};

                if (hairpinE.count(s) > 0) {
                    if (hairpin_energy > hairpinE[s]) {
                        hairpin_energy = hairpinE[s];
                    }
                }

            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(s) > 0) {
                        if (hairpin_energy > hairpinE[s]) {
                            hairpin_energy = hairpinE[s];
                        }
                    }
                }
            }
            break;
        case 5:
            if (i == 0) {
                xi_ = nucleotides[pa][x][1];
                xi2_ = nucleotides[pa][x][2];
                _2yj = nucleotides[pb][y][0];
                _yj = nucleotides[pb][y][1];
                s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                if (hairpinE.count(s) > 0) {
                    if (hairpin_energy > hairpinE[s]) {
                        hairpin_energy = hairpinE[s];
                    }
                }

            } else if (i == 1) {
                xi_ = nucleotides[pa][x][2];

                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi2_ = nucleotides[pna][x1][0];
                    _2yj = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(s) > 0) {
                        if (hairpin_energy > hairpinE[s]) {
                            hairpin_energy = hairpinE[s];
                        }
                    }
                }
            } else {
                _yj = nucleotides[pb][y][0];
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _2yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(s) > 0) {
                        if (hairpin_energy > hairpinE[s]) {
                            hairpin_energy = hairpinE[s];
                        }
                    }
                }
            }
            break;
        case 7:
            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi3_ = nucleotides[pna][x1][0];
                        _3yj = nucleotides[pna][x1][1];
                        _2yj = nucleotides[pna][x1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};

                        if (hairpinE.count(s) > 0) {
                            if (hairpin_energy > hairpinE[s]) {
                                hairpin_energy = hairpinE[s];
                            }
                        }
                    }
                } else {
                    _2yj = nucleotides[pb][y][0];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi2_ = nucleotides[pna][x1][0];
                        xi3_ = nucleotides[pna][x1][1];
                        _3yj = nucleotides[pna][x1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
//
                        if (hairpinE.count(s) > 0) {
                            if (hairpin_energy > hairpinE[s]) {
                                hairpin_energy = hairpinE[s];
                            }
                        }
                    }
                }
            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    xi3_ = nucleotides[pna][x1][2];
                    for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                        _3yj = nucleotides[ppb][y1][0];
                        _2yj = nucleotides[ppb][y1][1];
                        _yj = nucleotides[ppb][y1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
//
                        if (hairpinE.count(s) > 0) {
                            if (hairpin_energy > hairpinE[s]) {
                                hairpin_energy = hairpinE[s];
                            }
                        }
                    }
                }
            }

            break;
        default:
            break;
    }
    int temp_e;
    // general case
    for (int xi_ = 0; xi_ < 4; ++xi_) {
        if (((1 << xi_) & an_int) == 0) continue;
        for (int _yj = 0; _yj < 4; ++_yj) {
            if (((1 << _yj) & bp_int) == 0) continue;

            // check if the current values of xi_ and _yj are compatible in the same codon
            if (lb - la <= 4 && !rightCodon(la+1,lb-1,xi_,_yj)) continue; //&& (la+1)/3==(lb-1)/3

            temp_e = hairpin_loop(xi,yj,xi_,_yj,l-1);
            if (hairpin_energy > temp_e) {
                hairpin_energy = temp_e;
            }

        }
    }

    return hairpin_energy;
}


double Zuker::hairpin_CAI(double lambda, int l,int a, int b, int pa, int pb, int pna, int ppb, int n_codon_an, int n_codon_bp, int xi, int yj, int i, int j, int x, int y) {
    double hairpin_energy = inf;
    double mfe, cai;
    string s;
    double temp_mfe, temp_cai;
    double temp_e;
    static vector<int> temp;
    int la = sigma(a,i), lb = sigma(b,j);
    // special case l = 4,5,7
    switch (l) {
        int xi_, _yj, xi2_, _2yj, xi3_, _3yj;

        case 4:
            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                } else {
                    xi2_ = nucleotides[pb][y][0];
                }
                s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};
//
                if (hairpinE.count(s) > 0) {
                    temp_mfe = lambda*hairpinE[s];
                    temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y);
                    temp_e = temp_mfe + temp_cai;
                    if (hairpin_energy > temp_e) {
                        hairpin_energy = temp_e;
                        temp = {l,xi_,xi2_,_yj,a,b,x,y};
                        mfe = temp_mfe, cai = temp_cai;

                    }
                }
            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(s) > 0) {
                        temp_mfe = lambda*hairpinE[s];
                        temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1);
                        temp_e = temp_mfe + temp_cai;
                        if (hairpin_energy > temp_e) {
                            hairpin_energy = temp_e;
                            temp = {l,xi_,xi2_,_yj,a,b,x,y,a+1,x1};
                            mfe = temp_mfe, cai = temp_cai;
                        }
                    }

                }
            }

            break;
        case 5:
            if (i == 0) {
                xi_ = nucleotides[pa][x][1];
                xi2_ = nucleotides[pa][x][2];
                _2yj = nucleotides[pb][y][0];
                _yj = nucleotides[pb][y][1];
                s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                if (hairpinE.count(s) > 0) {
                    temp_mfe = lambda*hairpinE[s];
                    temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y);
                    temp_e = temp_mfe + temp_cai;
                    if (hairpin_energy > temp_e) {
                        hairpin_energy = temp_e;
                        temp = {l,xi_,xi2_,_2yj,_yj,a,b,x,y};
                        mfe = temp_mfe, cai = temp_cai;

                    }
                }

            } else if (i == 1) {
                xi_ = nucleotides[pa][x][2];

                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi2_ = nucleotides[pna][x1][0];
                    _2yj = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};
                    if (hairpinE.count(s) > 0) {
                        temp_mfe = lambda*hairpinE[s];
                        temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1);
                        temp_e = temp_mfe + temp_cai;
                        if (hairpin_energy > temp_e) {
                            hairpin_energy = temp_e;
                            temp = {l,xi_,xi2_,_2yj,_yj,a,b,x,y,a+1,x1};
                            mfe = temp_mfe, cai = temp_cai;
//
                        }
                    }
                }
            } else {
                _yj = nucleotides[pb][y][0];
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _2yj = nucleotides[pna][x1][2];
                    s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(s) > 0) {
                        temp_mfe = lambda*hairpinE[s];
                        temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1);


                        temp_e = temp_mfe + temp_cai;
                        if (hairpin_energy > temp_e) {
                            hairpin_energy = temp_e;
                            temp = {l,xi_,xi2_,_2yj,_yj,a,b,x,y,a+1,x1};
                            mfe = temp_mfe, cai = temp_cai;
//
                        }
                    }
                }
            }
//
            break;

        case 7:

            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi3_ = nucleotides[pna][x1][0];
                        _3yj = nucleotides[pna][x1][1];
                        _2yj = nucleotides[pna][x1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};

                        if (hairpinE.count(s) > 0) {
                            temp_mfe = lambda*hairpinE[s];
                            temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1);
                            temp_e = temp_mfe + temp_cai;
                            if (hairpin_energy > temp_e) {
                                hairpin_energy = temp_e;
                                temp = {l,xi_,xi2_,xi3_,_3yj,_2yj,_yj,a,b,x,y,a+1,x1};
                                mfe = temp_mfe, cai = temp_cai;
                            }
                        }
                    }
                } else {
                    _2yj = nucleotides[pb][y][0];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi2_ = nucleotides[pna][x1][0];
                        xi3_ = nucleotides[pna][x1][1];
                        _3yj = nucleotides[pna][x1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};

                        if (hairpinE.count(s) > 0) {
                            temp_mfe = lambda*hairpinE[s];
                            temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1);
                            temp_e = temp_mfe + temp_cai;
                            if (hairpin_energy > temp_e) {
                                hairpin_energy = temp_e;
                                temp = {l,xi_,xi2_,xi3_,_3yj,_2yj,_yj,a,b,x,y,a+1,x1};
                                mfe = temp_mfe, cai = temp_cai;
                            }
                        }
                    }
                }
            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    xi3_ = nucleotides[pna][x1][2];
                    for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                        _3yj = nucleotides[ppb][y1][0];
                        _2yj = nucleotides[ppb][y1][1];
                        _yj = nucleotides[ppb][y1][2];
                        s = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
                        if (hairpinE.count(s) > 0) {
                            temp_mfe = lambda*hairpinE[s];
                            temp_cai = (lambda-1)*add_hairpin_CAI_8(a,b,x,y,a+1,x1,b-1,y1);
                            temp_e = temp_mfe + temp_cai;
                            if (hairpin_energy > temp_e) {
                                hairpin_energy = temp_e;
                                temp = {l,xi_,xi2_,xi3_,_3yj,_2yj,_yj,a,b,x,y,a+1,x1,b-1,y1};
                                mfe = temp_mfe, cai = temp_cai;
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;

    }


    int xi_, _yj;

    // general case

    if (i == 2 && j == 0) {
        for (int x1 = 0; x1 < n_codon_an; ++x1) {
            if (a + 1 == b && x1 != y) continue;
            xi_ = nucleotides[pna][x1][0];
            for (int y1 = 0; y1 < n_codon_bp; y1++) {
                if (a + 1 == b && y1 != x) continue;
                if (a+1 == b-1 && x1 != y1) continue;
                _yj = nucleotides[ppb][y1][2];
                if (lb - la <= 4 && (la+1)/3==(lb-1)/3 && !rightCodon(la+1,lb-1,xi_,_yj)) continue;
                temp_mfe = lambda*hairpin_loop(xi,yj,xi_,_yj,l-1);
                temp_cai = (lambda-1)*add_hairpin_CAI_2(a,b,x,y,a+1,x1,b-1,y1);
                temp_e = temp_mfe + temp_cai;
                if (hairpin_energy > temp_e) {
                    hairpin_energy = temp_e;
                    temp = {l,xi,xi_,_yj,yj,a+1,x1,b-1,y1,-1};
                    mfe = temp_mfe, cai = temp_cai;
                }

            }
        }

    }
    else if (i == 2) {
        _yj = nucleotides[pb][y][j-1];
        for (int x1 = 0; x1 < n_codon_an; ++x1) {
            if (a+1 == b && x1 != y) continue;
            xi_ = nucleotides[pna][x1][0];
            if (lb - la <= 4 && (la+1)/3==(lb-1)/3 && !rightCodon(la+1,lb-1,xi_,_yj)) continue;
            temp_mfe = lambda*hairpin_loop(xi,yj,xi_,_yj,l-1);
            temp_cai = (lambda-1)*add_hairpin_CAI_2(a,b,x,y,a+1,x1);
            temp_e = temp_mfe + temp_cai;
            if (hairpin_energy > temp_e) {
                hairpin_energy = temp_e;
                temp = {l,xi,xi_,_yj,yj,a+1,x1,b,y,-1};
                mfe = temp_mfe, cai = temp_cai;
            }
        }
    }
    else if (j == 0) {
        xi_ = nucleotides[pa][x][i+1];
        for (int y1 = 0; y1 < n_codon_bp; y1++) {
            if (a+1 == b && y1 != x) continue;
            _yj = nucleotides[ppb][y1][2];
            if (lb - la <= 4 && (la+1)/3==(lb-1)/3 && !rightCodon(la+1,lb-1,xi_,_yj)) continue;
            temp_mfe = lambda*hairpin_loop(xi,yj,xi_,_yj,l-1);
            temp_cai = (lambda-1)*add_hairpin_CAI_2(a,b,x,y,b-1,y1);
            temp_e = temp_mfe + temp_cai;
            if (hairpin_energy > temp_e) {
                hairpin_energy = temp_e;
                temp = {l,xi,xi_,_yj,yj,a,x,b-1,y1,-1};
                mfe = temp_mfe, cai = temp_cai;
            }
        }

    }
    else {
        xi_ = nucleotides[pa][x][i+1];
        _yj = nucleotides[pb][y][j-1];
        if (lb - la <= 4 && (la+1)/3==(lb-1)/3 && !rightCodon(la+1,lb-1,xi_,_yj)) {

        } else {
            temp_mfe = lambda*hairpin_loop(xi,yj,xi_,_yj,l-1);
            temp_cai = (lambda-1)*add_hairpin_CAI_2(a,b,x,y);
            temp_e = temp_mfe + temp_cai;
            if (hairpin_energy > temp_e) {
                hairpin_energy = temp_e;
                temp = {l,xi,xi_,_yj,yj,a,x,b,y,-1};
                mfe = temp_mfe, cai = temp_cai;
            }
        }

    }

    if (E1[index(a,b,i,j,x,y)] > hairpin_energy) {
        E1[index(a,b,i,j,x,y)] = hairpin_energy;
        E_bt[index(a,b,i,j,x,y)] = temp;
        E2[index(a,b,i,j,x,y)] = mfe;
        E_CAI[index(a,b,i,j,x,y)] = cai;

    }

    return hairpin_energy;
}

bool Zuker::rightCodon(int l1, int l2, int x, int y) const {

    int l = l1/3;
    if (l != l2/3) return true;
    int p = protein[l];
    int p1 = l1 - 3*l;
    int p2 = l2 - 3*l;
    return codonMatch[p][p1][p2][x][y];
}

int Zuker::internal(int a, int b, int i, int j, int x, int y, int la, int lb, int xi, int yj, int an_int, int bp_int) {
    int internal_energy = inf, energy = inf;
    int t = 0;
    static vector<int> H(96,inf);
    static vector<int> T(256,inf);

    static bool updatedH[96];
    memset(updatedH, false, 96);

    int idx_m = index(a,b,i,j,x,y);
    for (int lc = la + 1; lc <= min(lb-5, la+MAXLOOP+1); lc++) {
        int ll = lc - la;
        int c = lc / 3;
        int i1 = lc % 3;
        int min_ld = max(lc + 4, lb - la + lc - MAXLOOP - 2);

        for (int ld = lb-1; ld >= min_ld; ld--) {
            int lr = lb - ld;
            int d = ld / 3;
            int j1 = ld % 3;
            int pc = protein[c];
            int pd = protein[d];
            const int n_codon_c = n_codon[pc];
            const int n_codon_d = n_codon[pd];

            memset(updatedH, false, 96);



            for (int xh = 0; xh < n_codon_c; ++xh) {
                if (c == a && xh != x) continue;

                // hi and kj are the other paired ends
                int hi = nucleotides[pc][xh][i1];

                for (int xk = 0; xk < n_codon_d; ++xk) {
                    if (d == b && xk != y) continue;
                    int kj = nucleotides[pd][xk][j1];

                    int interior_energy = inf, en = inf;
                    t = 0;

                    int type2 = BP_pair[hi+1][kj+1];

                    if (type2 == 0) continue;

                    int _hi = 0, kj_ = 0;
                    if (i1 >= 1) {
                        _hi = nucleotides[pc][xh][i1-1];
                    }
                    if (j1 <= 1) {
                        kj_ = nucleotides[pd][xk][j1+1];
                    }

                    int idx = H_index(type2, _hi, kj_); // ((type2-1)<<4)+(_hi<<2)+kj_; // ;
                    if (!updatedH[idx]) {
                        updatedH[idx] = true;
                    }
                    else {
                        internal_energy = min(internal_energy, Access_E(c,d,i1,j1,xh,xk) + H[idx]);
                        if (energy > internal_energy) {
                            energy = internal_energy;
                            Access_EB(idx_m) = T[idx];
                        }
                        continue;
                    }


                    if (ll == 1 && lr == 1) {
                        interior_energy = min(interior_energy, stacking(xi,yj,hi,kj));
                        if (en > interior_energy) {
                            en = interior_energy;
                            t = -4;
                        }
                    }
                    else if (ll == 1 || lr == 1) {
                        interior_energy = min(interior_energy, bulge_loop(xi,yj,hi,kj,max(ll,lr)-1));
                        if (en > interior_energy) {
                            en = interior_energy;
                            t = -5;
                        }
                    }
                    else {
                        // internal loop structure
                        //   hi -- kj
                        // _hi      kj_
                        //     ....
                        // xi_      _yj
                        //   xi -- yj
                        int cp_int = ava_nucle_m[index(c,xh,i1)];
                        int dn_int = ava_nucle_p[index(d,xk,j1)];

                        for (int xi_ = 0; xi_ < 4; ++xi_) {
                            if (((1 << xi_) & an_int) == 0) continue;
                            if (lc -1 - la <= 2 && !rightCodon(la+1, lc,xi_,hi)) continue;
                            for (_hi = 0; _hi < 4; ++_hi) {
                                if (((1 << _hi) & cp_int) == 0) continue;

                                // check if the current values of xi_ and _hi are compatible in the same codon
                                if (lc - la <= 4  && !rightCodon(la+1,lc-1,xi_,_hi)) continue; //&& (la+1)/3==(lc-1)/3

                                for (kj_ = 0; kj_ < 4; ++kj_) {
                                    if (((1 << kj_) & dn_int) == 0) continue;
                                    if (lb -1 - ld <= 2 && !rightCodon(ld+1, lb,kj_,yj)) continue;
                                    for (int _yj = 0; _yj < 4; ++_yj) {
                                        if (((1 << _yj) & bp_int) == 0) continue;
                                        // check if the current values of kj_ and _yj are compatible in the same codon
                                        if (lb - ld <= 4 &&  !rightCodon(ld+1,lb-1,kj_,_yj)) continue; //(ld+1)/3==(lb-1)/3 &&

                                        interior_energy = min(interior_energy, interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1));

                                    }
                                }
                            }
                        }
                        if (en > interior_energy) {
                            en = interior_energy;
                            t = -6;
                        }
                    }

                    H[idx] = interior_energy;
                    T[idx] = t;
                    internal_energy = min(internal_energy, Access_E(c,d,i1,j1,xh,xk) + interior_energy);

                    if (energy > internal_energy) {
                        energy = internal_energy;
                        Access_EB(idx_m) = t;
                    }


                }
            }
        }
    }

    return internal_energy;
}

double Zuker::internal_CAI(double lambda, int a, int b,int i, int j, int x, int y, int la, int lb, int xi, int yj) {
    double internal_energy = inf, energy = inf;
    double mfe2, cai2, mfe, cai;
    int t = 0;

    // store back pointers
    static vector<int> tmp, temp, temp11;

    int idx_m = index(a,b,i,j,x,y);
    int cx = 0,cy = 0;
    int xi_ = 0, _yj = 0;
    if (i <= 1) {
        xi_ = nucleotides[protein[a]][x][i+1];
        cx = 1;
    }
    if (j >= 1) {
        _yj = nucleotides[protein[b]][y][j-1];
        cy = 2;
    }
    for (int lc = la + 1; lc <= min(lb-5, la+MAXLOOP+1); lc++) {
        int ll = lc - la;
        int c = lc / 3;
        int i1 = lc % 3;
        int min_ld = lb - la + lc - MAXLOOP - 2;
        if (min_ld < lc + 4) min_ld = lc + 4;
        for (int ld = lb-1; ld >= min_ld; ld--) {
            int lr = lb - ld;
            int d = ld / 3;
            int j1 = ld % 3;
            int pc = protein[c];
            int pd = protein[d];
            const int n_codon_c = n_codon[pc];
            const int n_codon_d = n_codon[pd];
            double temp_e, t_mfe, t_cai;


            for (int xh = 0; xh < n_codon_c; ++xh) {

                if (c == a && xh != x) continue;
                int hi = nucleotides[pc][xh][i1];
                for (int xk = 0; xk < n_codon_d; ++xk) {
                    if (d == b && xk != y) continue;
                    int kj = nucleotides[pd][xk][j1];

                    double interior_energy = inf, en = inf;
                    double temp_mfe, temp_cai;
                    mfe = inf, cai = inf;
                    t = 0;

                    int type2 = BP_pair[hi+1][kj+1];
                    if (type2 == 0) continue;


                    int _hi = 0, kj_ = 0;
                    int ch = 0, ck = 0;

                    if (i1 >= 1) {
                        _hi = nucleotides[pc][xh][i1-1];
                        ch = 5;
                    }
                    if (j1 <= 1) {
                        kj_ = nucleotides[pd][xk][j1+1];
                        ck = 9;
                    }

                    int tt = cx + cy + ch + ck;
                    int ppc = protein[c-1], pnd = protein[d+1];
                    int pna = protein[a+1], ppb = protein[b-1];
                    const int n_codon_pc = n_codon[ppc];
                    const int n_codon_dn = n_codon[pnd];
                    const int n_codon_an = n_codon[pna];
                    const int n_codon_pb = n_codon[ppb];


                    if (ll == 1 && lr == 1) {
                        temp_mfe = lambda*stacking(xi,yj,hi,kj);
                        temp_cai = (lambda-1)*(add_CAI(a,c,x) + add_CAI(b,d,y));
                        interior_energy = min(interior_energy, temp_mfe + temp_cai);
                        if (en > interior_energy) {
                            en = interior_energy;
                            temp = {c,d,i1,j1,xh,xk,hi,kj};

                            t = -4;
                            mfe = temp_mfe, cai = temp_cai;
                        }

                    }
                    else if (ll == 1 || lr == 1) {
                        temp_mfe = lambda*bulge_loop(xi,yj,hi,kj,max(ll,lr)-1);
                        temp_cai = (lambda-1)*(add_CAI(a,c,x) + add_CAI(b,d,y));
                        interior_energy = min(interior_energy, temp_mfe + temp_cai);
                        if (en > interior_energy) {
                            en = interior_energy;
                            temp = {c,d,i1,j1,xh,xk,hi,kj,max(ll,lr)-1};
                            t = -5;
                            mfe = temp_mfe, cai = temp_cai;
                        }
                    }
                    else {

                        switch (tt) {
                            case 0:
                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];
                                    for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                        if (a+1 == c && x != hx1) continue;
                                        if (a+1 == c-1 && x1 != hx1) continue;
                                        _hi = nucleotides[ppc][hx1][2];
                                        for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                            if (d+1 == b && y1 != xk) continue;
                                            _yj = nucleotides[ppb][y1][2];
                                            for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                                if (d+1 == b && y != ky1) continue;
                                                if (d+1 == b-1 && y1 != ky1) continue;
                                                kj_ = nucleotides[pnd][ky1][0];
                                                t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                                t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1,c-1,hx1) + add_interior_CAI_2(b,d,y,b-1,y1,d+1,ky1));
                                                temp_e = t_mfe + t_cai;
                                                if (interior_energy > temp_e) {
                                                    interior_energy = temp_e;
                                                    temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c-1,hx1,b-1,y1,d+1,ky1};
                                                    temp_mfe = t_mfe, temp_cai = t_cai;
                                                }

                                            }
                                        }
                                    }
                                }
                                break;
                            case 1:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }

                                for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                    if (d+1 == b && y1 != xk) continue;
                                    _yj = nucleotides[ppb][y1][2];
                                    for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                        if (a+1 == c && x != hx1) continue;
                                        _hi = nucleotides[ppc][hx1][2];
                                        for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                            if (d+1 == b && y != ky1) continue;
                                            if (d+1 == b-1 && y1 != ky1) continue;
                                            kj_ = nucleotides[pnd][ky1][0];
                                            t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                            t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,c-1,hx1) + add_interior_CAI_2(b,d,y,b-1,y1,d+1,ky1));
                                            temp_e = t_mfe + t_cai;
                                            if (interior_energy > temp_e) {
                                                interior_energy = temp_e;
                                                temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c-1,hx1,b-1,y1,d+1,ky1};
                                                temp_mfe = t_mfe, temp_cai = t_cai;

                                            }

                                        }
                                    }
                                }
                                break;
                            case 2:
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                }

                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];

                                    for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                        if (a+1 == c && x != hx1) continue;
                                        if (a+1 == c-1 && x1 != hx1) continue;
                                        _hi = nucleotides[ppc][hx1][2];
                                        for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                            if (d+1 == b && y != ky1) continue;
                                            kj_ = nucleotides[pnd][ky1][0];
                                            t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                            t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1,c-1,hx1) + add_interior_CAI_2(b,d,y,d+1,ky1));

                                            temp_e = t_mfe + t_cai;
                                            if (interior_energy > temp_e) {
                                                interior_energy = temp_e;
                                                temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c-1,hx1,b,y,d+1,ky1};
                                                temp_mfe = t_mfe, temp_cai = t_cai;
//
                                            }
                                        }
                                    }
                                }
                                break;
                            case 5:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                }

                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];
                                    for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                        if (d+1 == b && y1 != xk) continue;
                                        _yj = nucleotides[ppb][y1][2];
                                        for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                            if (d+1 == b && y != ky1) continue;
                                            if (d+1 == b-1 && y1 != ky1) continue;
                                            kj_ = nucleotides[pnd][ky1][0];
                                            t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                            t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1) + add_interior_CAI_2(b,d,y,b-1,y1,d+1,ky1));
                                            temp_e = t_mfe + t_cai;
                                            if (interior_energy > temp_e) {
                                                interior_energy = temp_e;
                                                temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c,xh,b-1,y1,d+1,ky1};
                                                temp_mfe = t_mfe, temp_cai = t_cai;
//
                                            }
                                        }
                                    }
                                }
                                break;
                            case 9:
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }

                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];
                                    for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                        if (d+1 == b && y1 != xk) continue;
                                        _yj = nucleotides[ppb][y1][2];
                                        for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                            if (a+1 == c && x != hx1) continue;
                                            if (a+1 == c-1 && x1 != hx1) continue;
                                            _hi = nucleotides[ppc][hx1][2];
                                            t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                            t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1,c-1,hx1) + add_interior_CAI_2(b,d,y,b-1,y1));
                                            temp_e = t_mfe + t_cai;
                                            if (interior_energy > temp_e) {
                                                interior_energy = temp_e;
                                                temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c-1,hx1,b-1,y1,d,xk};
                                                temp_mfe = t_mfe, temp_cai = t_cai;
//
                                            }
                                        }
                                    }
                                }
                                break;
                            case 3:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                }
                                for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                    if (a+1 == c && x != hx1) continue;
                                    _hi = nucleotides[ppc][hx1][2];
                                    for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                        if (d+1 == b && y != ky1) continue;
                                        kj_ = nucleotides[pnd][ky1][0];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,c-1,hx1) + add_interior_CAI_2(b,d,y,d+1,ky1));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c-1,hx1,b,y,d+1,ky1};
                                            temp_mfe = t_mfe, temp_cai = t_cai;
//
                                        }

                                    }
                                }
                                break;
                            case 6:
                                if (lc - la <= 4  && !rightCodon(la+1,lc-1,xi_,_hi)) continue; //&& (la+1)/3==(lc-1)/3
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                    if (d+1 == b && y1 != xk) continue;
                                    _yj = nucleotides[ppb][y1][2];
                                    for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                        if (d+1 == b && y != ky1) continue;
                                        if (d+1 == b-1 && y1 != ky1) continue;
                                        kj_ = nucleotides[pnd][ky1][0];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x) + add_interior_CAI_2(b,d,y,b-1,y1,d+1,ky1));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c,xh,b-1,y1,d+1,ky1};
                                            temp_mfe = t_mfe, temp_cai = t_cai;

                                        }
                                    }
                                }
                                break;
                            case 10:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                    if (d+1 == b && y1 != xk) continue;
                                    _yj = nucleotides[ppb][y1][2];
                                    for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                        if (a+1 == c && x != hx1) continue;
                                        _hi = nucleotides[ppc][hx1][2];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,c-1,hx1) + add_interior_CAI_2(b,d,y,b-1,y1));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c-1,hx1,b-1,y1,d,xk};
                                            temp_mfe = t_mfe, temp_cai = t_cai;
//
                                        }
                                    }
                                }
                                break;
                            case 7:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                }
                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];
                                    for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                        if (d+1 == b && y != ky1) continue;
                                        kj_ = nucleotides[pnd][ky1][0];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1) + add_interior_CAI_2(b,d,y,d+1,ky1));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c,xh,b,y,d+1,ky1};
                                            temp_mfe = t_mfe, temp_cai = t_cai;
//
                                        }
                                    }
                                }
                                break;
                            case 11:
                                if (lb - ld <= 4  && !rightCodon(ld+1,lb-1,kj_,_yj)) continue; //&& (ld+1)/3==(lb-1)/3
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }

                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];

                                    for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                        if (a+1 == c && x != hx1) continue;
                                        if (a+1 == c-1 && x1 != hx1) continue;
                                        _hi = nucleotides[ppc][hx1][2];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1,c-1,hx1) + add_interior_CAI_2(b,d,y));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c-1,hx1,b,y,d,xk};

                                            temp_mfe = t_mfe, temp_cai = t_cai;
//
                                        }

                                    }
                                }
                                break;
                            case 14:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];

                                    for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                        if (d+1 == b && y1 != xk) continue;
                                        _yj = nucleotides[ppb][y1][2];
                                        t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                        t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1) + add_interior_CAI_2(b,d,y,b-1,y1));
                                        temp_e = t_mfe + t_cai;
                                        if (interior_energy > temp_e) {
                                            interior_energy = temp_e;
                                            temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c,xh,b-1,y1,d,xk};

                                            temp_mfe = t_mfe, temp_cai = t_cai;
//

                                        }
                                    }
                                }
                                break;
                            case 8:
                                if (lc - la <= 4  && !rightCodon(la+1,lc-1,xi_,_hi)) continue; //&& (la+1)/3==(lc-1)/3
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                }
                                for (int ky1 = 0; ky1 < n_codon_dn; ++ky1) {
                                    if (d+1 == b && y != ky1) continue;
                                    kj_ = nucleotides[pnd][ky1][0];
                                    t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                    t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x) + add_interior_CAI_2(b,d,y,d+1,ky1));
                                    temp_e = t_mfe + t_cai;
                                    if (interior_energy > temp_e) {
                                        interior_energy = temp_e;
                                        temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c,xh,b,y,d+1,ky1};
                                        temp_mfe = t_mfe, temp_cai = t_cai;
//
                                    }
                                }
                                break;
                            case 12:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb - ld <= 4  && !rightCodon(ld+1,lb-1,kj_,_yj)) continue; //&& (ld+1)/3==(lb-1)/3
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                for (int hx1 = 0; hx1 < n_codon_pc; ++hx1) {
                                    if (a+1 == c && x != hx1) continue;
                                    _hi = nucleotides[ppc][hx1][2];
                                    t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                    t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,c-1,hx1) + add_interior_CAI_2(b,d,y));
                                    temp_e = t_mfe + t_cai;
                                    if (interior_energy > temp_e) {
                                        interior_energy = temp_e;
                                        temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c-1,hx1,b,y,d,xk};
                                        temp_mfe = t_mfe, temp_cai = t_cai;
//
                                    }
                                }
                                break;
                            case 15:
                                if (lc - la <= 4  && !rightCodon(la+1,lc-1,xi_,_hi)) continue; //&& (la+1)/3==(lc-1)/3
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                for (int y1 = 0; y1 < n_codon_pb; ++y1) {
                                    if (d+1 == b && y1 != xk) continue;
                                    _yj = nucleotides[ppb][y1][2];
                                    t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                    t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x) + add_interior_CAI_2(b,d,y,b-1,y1));
                                    temp_e = t_mfe + t_cai;
                                    if (interior_energy > temp_e) {
                                        interior_energy = temp_e;
                                        temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c,xh,b-1,y1,d,xk};
                                        temp_mfe = t_mfe, temp_cai = t_cai;
//
                                    }
                                }
                                break;
                            case 16:
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                }
                                if (lb - ld <= 4  && !rightCodon(ld+1,lb-1,kj_,_yj)) continue; //&& (ld+1)/3==(lb-1)/3
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                    if (a+1 == c && x1 != xh) continue;
                                    xi_ = nucleotides[pna][x1][0];
                                    t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                    t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x,a+1,x1) + add_interior_CAI_2(b,d,y));
                                    temp_e = t_mfe + t_cai;
                                    if (interior_energy > temp_e) {
                                        interior_energy = temp_e;
                                        temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a+1,x1,c,xh,b,y,d,xk};

                                        temp_mfe = t_mfe, temp_cai = t_cai;
                                    }
                                }
                                break;
                            default:
                                if (lc - la <= 4  && !rightCodon(la+1,lc-1,xi_,_hi)) continue; //&& (la+1)/3==(lc-1)/3
                                if (lc -1 - la <= 2) {
                                    // check if the current values of xi and _hi are compatible in the same codon
                                    if ( !rightCodon(la, lc-1,xi,_hi)) continue;
                                    // check if the current values of xi_ and hi are compatible in the same codon
                                    if ( !rightCodon(la+1, lc,xi_,hi)) continue;
                                }
                                if (lb - ld <= 4  && !rightCodon(ld+1,lb-1,kj_,_yj)) continue; //&& (ld+1)/3==(lb-1)/3
                                if (lb -1 - ld <= 2) {
                                    // check if the current values of kj and _yj are compatible in the same codon
                                    if (!rightCodon(ld, lb-1,kj,_yj)) continue; //ld/3==(lb-1)/3 &&
                                    // check if the current values of kj_ and yj are compatible in the same codon
                                    if ( !rightCodon(ld+1, lb,kj_,yj)) continue; //(ld+1)/3==lb/3 &&
                                }
                                t_mfe = lambda*interior_loop(xi,yj,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1);
                                t_cai = (lambda-1)*(add_interior_CAI_2(a,c,x) + add_interior_CAI_2(b,d,y));
                                temp_e = t_mfe + t_cai;
                                if (interior_energy > temp_e) {
                                    interior_energy = temp_e;
                                    temp11 = {c,d,i1,j1,xh,xk,hi,kj,xi_,_yj,_hi,kj_,ll-1,lr-1,a,x,c,xh,b,y,d,xk};
                                    temp_mfe = t_mfe, temp_cai = t_cai;
                                }
                                break;
                        }


                        if (en > interior_energy) {
                            en = interior_energy;
                            t = -6;
                            temp = temp11;
                            mfe = temp_mfe, cai = temp_cai;
                        }
                    }

                    internal_energy = min(internal_energy, Access_E1(c,d,i1,j1,xh,xk) + interior_energy);

                    if (energy > internal_energy) {
                        energy = internal_energy;
                        Access_EB(idx_m) = t;
                        E_bt[idx_m] = temp;
                        mfe2 = Access_E2(c,d,i1,j1,xh,xk) + mfe;
                        cai2 = E_CAI[index(c,d,i1,j1,xh,xk)] + cai;

                    }

                }
            }
        }
    }
    if (E1[idx_m] > internal_energy) {
        E1[idx_m] = internal_energy;
        E2[idx_m] = mfe2;
        E_CAI[idx_m] = cai2;
    }


    return internal_energy;
}

double Zuker::add_interior_CAI_2(int a, int c, int x, int na, int x1, int pc, int h1) const {
    if (a == c) return 0;
    double cai = codon_cai[protein[a]][x];
    if (na == -1 && x1 == -1 && pc == -1 && h1 == -1) {
        return cai;
    } else if (pc == -1 && h1 == -1) {
        if (na != a && na != c) {
            cai += codon_cai[protein[na]][x1];
        }

        return cai;
    } else {

        if (na != a && na != c) {
            cai += codon_cai[protein[na]][x1];
        }
        if (pc != c && pc != a && pc != na) {
            cai += codon_cai[protein[pc]][h1];
        }
        return cai;
    }

}

double Zuker::add_CAI(int p1, int p2, int x) const {
    if (p1 == p2) return 0;
    double cai = codon_cai[protein[p1]][x];
    return cai;
}

int Zuker::multi_loop(int a, int b, int i, int j, int x, int y, int pa, int pb, int n_codon_an, int n_codon_bp) {
    int multi_loop = inf;
    if (a < b-2 && i == 2 && j == 0) {
        for (int x1 = 0; x1 < n_codon_an; ++x1) {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                multi_loop = min(multi_loop, TM0[index(a+1,b-1,0,2,x1,y1)]); // , idx_t + 6*x1 + y1
            }
        }
    }

    if (a < b-1) {
        if (i == 2 && j >= 1) {
            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                multi_loop = min(multi_loop, TM0[index(a+1,b,0,j-1,x1,y)]);// , idx_t + 6*x1
            }
        }
        if (i <= 1 && j == 0) {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                multi_loop = min(multi_loop, TM0[index(a,b-1,i+1,2,x,y1)]); // , idx_t + y1
            }
        }

    }

    if (a < b && i <= 1 && j >= 1) {

        multi_loop = min(multi_loop, TM0[index(a,b,i+1,j-1,x,y)]); //index(a,b,i+1,j-1,x,y)
    }
    multi_loop = multi_loop + AU[nucleotides[pa][x][i]][nucleotides[pb][y][j]] + ML_closing + ML_intern;

    return multi_loop;
}

double Zuker::multi_loop_CAI(double lambda,int a, int b, int i, int j, int x, int y, int pa, int pb, int n_codon_an, int n_codon_bp) {
    double multi_loop = inf;
    double temp_e;
    double mfe, cai;
    static vector<int> temp;

    if (a < b-2 && i == 2 && j == 0) {
        for (int x1 = 0; x1 < n_codon_an; ++x1) {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                temp_e =  TM[index(a+1,b-1,0,2,x1,y1)] + (lambda-1)*(codon_cai[pa][x] + codon_cai[pb][y]); //idx_t + 6*x1 + y1
                if (multi_loop > temp_e) {
                    multi_loop = temp_e;
                    temp = {a+1,b-1,0,2,x1,y1};
                    mfe = TM2[index(a+1,b-1,0,2,x1,y1)];
                    cai = TM_CAI[index(a+1,b-1,0,2,x1,y1)] + (lambda-1)*(codon_cai[pa][x] + codon_cai[pb][y]);
                }
            }
        }
    }

    if (a < b-1) {
        if (i == 2 && j >= 1) {
            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                temp_e = TM[index(a+1,b,0,j-1,x1,y)] + (lambda-1)*codon_cai[pa][x]; //idx_t + 6*x1
                if (multi_loop > temp_e) {
                    multi_loop = temp_e;
                    temp = {a+1,b,0,j-1,x1,y};
                    mfe = TM2[index(a+1,b,0,j-1,x1,y)];
                    cai = TM_CAI[index(a+1,b,0,j-1,x1,y)] + (lambda-1)*codon_cai[pa][x];
                }
            }
        }
        if (i <= 1 && j == 0) {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                temp_e = TM[index(a,b-1,i+1,2,x,y1)] + (lambda-1)*codon_cai[pb][y]; //idx_t + y1
                if (multi_loop > temp_e) {
                    multi_loop = temp_e;
                    temp = {a,b-1,i+1,2,x,y1};
                    mfe = TM2[index(a,b-1,i+1,2,x,y1)];
                    cai = TM_CAI[index(a,b-1,i+1,2,x,y1)] + (lambda-1)*codon_cai[pb][y];
                }
            }
        }
    }

    if (a < b && i <= 1 && j >= 1) {
        temp_e = TM[index(a,b,i+1,j-1,x,y)];
        if (multi_loop > temp_e) {
            multi_loop = temp_e;
            temp = {a,b,i+1,j-1,x,y};
            mfe = TM2[index(a,b,i+1,j-1,x,y)];
            cai = TM_CAI[index(a,b,i+1,j-1,x,y)];
        }
    }

    multi_loop = multi_loop + lambda*(AU[nucleotides[pa][x][i]][nucleotides[pb][y][j]] + ML_closing + ML_intern);

    E_bt[index(a,b,i,j,x,y)] = temp;
    if (multi_loop < E1[index(a,b,i,j,x,y)]) {
        E2[index(a,b,i,j,x,y)] = mfe + lambda*(AU[nucleotides[pa][x][i]][nucleotides[pb][y][j]] + ML_closing + ML_intern);
        E_CAI[index(a,b,i,j,x,y)] = cai;
    }
    return multi_loop;
}

int Zuker::calculate_M(int a, int b, int i, int j, int x, int y) {

    int min_energy = inf, energy = inf;
    int t = 0;
    int la = sigma(a,i);
    int lb = sigma(b,j);
    int pna = protein[a+1];
    int ppb = protein[b-1];
    const int n_codon_an = n_codon[pna];
    const int n_codon_bp = n_codon[ppb];
    int pa = protein[a];
    int pb = protein[b];
    int ni = nucleotides[pa][x][i];
    int nj = nucleotides[pb][y][j];
    int type = BP_pair[ni+1][nj+1];

    int idx = index(a,b,i,j,x,y);

    if (type) {

        min_energy = min(min_energy, E[idx] + ML_intern + AU[ni][nj]); //Access_E(idx)
        if (energy > min_energy) {
            energy = min_energy;
            t = -1;
        }
    }

    if (i == 2) {

        {
            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                min_energy = min(min_energy, Access_M(a + 1, b, 0, j, x1, y) + ML_BASE); //, idx_t + 6*x1
            }
        }

        if (energy > min_energy) {
            energy = min_energy;
            t = -3;
        }
    }else {

        min_energy = min(min_energy, Access_M(a,b,i+1,j,x,y) + ML_BASE); //, idx-108
        if (energy > min_energy) {
            energy = min_energy;
            t = -2;
        }

    }

    if (j == 0) {

        {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                min_energy = min(min_energy, Access_M(a, b - 1, i, 2, x, y1) + ML_BASE); // idx_t + y1
            }
        }
        if (energy > min_energy) {
            energy = min_energy;
            t = -5;
        }

    } else {
            min_energy = min(min_energy, Access_M(a,b,i,j-1,x,y) + ML_BASE); //, idx-36
            if (energy > min_energy) {
                energy = min_energy;
                t = -4;
            }

    }


    // bifurication
    int bi_energy = inf;
    for (int lc = la + 5; lc <= lb-4; lc++) {
        int c = lc / 3;
        int i1 = lc % 3;
        int pc = protein[c];
        int n_codon_c = n_codon[pc];
        for (int hx = 0; hx < n_codon_c; ++hx) {
            if (i1 >= 1) {
                bi_energy = min(bi_energy, Access_M(a,c,i,i1-1,x,hx) + Access_M(c,b,i1,j,hx,y)); //, idx_1 + hx //, idx_2 + 6*hx
            }
            else {

                {
                    int n_codon_cp = n_codon[protein[c-1]];
                    for (int ky = 0; ky < n_codon_cp; ++ky) {
                        bi_energy = min(bi_energy, Access_M(a,c-1,i,2,x,ky) + Access_M(c,b,i1,j,hx,y)); //, idx_3 + ky //, idx_4 + 6*hx
                    }
                }
            }
        }

    }
    TM0[idx] = min(bi_energy, TM0[idx]);


    min_energy = min(min_energy, bi_energy);
    if (energy > min_energy) {
        energy = min_energy;
        t = -6;
    }

    Access_M(idx) = min_energy;
    Access_MB(idx) = t;
    return min_energy;
}

void Zuker::traceback() {
    int s = 0;
    int t = 0;
    sector[++s].a = 0;
    sector[s].b = n-1;
    sector[s].i = 0;
    sector[s].j = 2;
    sector[s].x = minX;
    sector[s].y = minY;
    sector[s].ml = 0;


    OUTLOOP:
    while (s > 0) {
        int zij, zi, eij; //zj
        int traced_ab, traced_ij, traced_xy;
        int c,d,i1,j1,hx,ky;
        int an_int, bp_int;
        int a = sector[s].a;
        int b = sector[s].b;
        int i  = sector[s].i;
        int j  = sector[s].j;
        int x  = sector[s].x;
        int y  = sector[s].y;
        int ml = sector[s--].ml;


        int pa = protein[a];
        int pb = protein[b];
        int li = sigma(a,i);
        int rj = sigma(b,j);

        int xi = nucleotides[pa][x][i];
        int yj = nucleotides[pb][y][j];

        int kj;

        if (a == b && i == j) break;

        int n_codon_bp = n_codon[protein[b-1]];
        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;

        if (ml==2) {
            bp_bond[++t].i = li;
            bp_bond[t].j   = rj;
            goto repeat1;
        }

        zij = (ml == 0) ? Access_Z(a,b,i,j,x,y) : Access_M(a,b,i,j,x,y);

        if (j >= 1) {
            zi = (ml == 0) ? Access_Z(a,b,i,j-1,x,y) : Access_M(a,b,i,j-1,x,y) + ML_BASE;
            if (zij == zi) kj = y;
        }
        else {
            for (ky = 0; ky < n_codon_bp; ++ky) {
                if (a == b-1 && x != ky) continue;
                zi = (ml == 0) ? Access_Z(a,b-1,i,2,x,ky) : Access_M(a,b-1,i,2,x,ky) + ML_BASE;
                if (zij == zi) {
                    kj = ky;
                    break;
                }
            }
        }


        if (zij == inf || zij == inf+ML_BASE || zij == inf+ML_intern) {
            goto repeat1;
        }

        if (zij == zi) {
            int tb, tj;
            if (j >= 1) {
                tb = b;
                tj = j-1;
            }
            else {
                tb = b-1;
                tj = 2;
            }

            sector[++s].a = a;
            sector[s].b = tb, sector[s].i = i, sector[s].j = tj, sector[s].x = x, sector[s].y = kj, sector[s].ml = ml;
            goto OUTLOOP;
        }

        if (ml == 0) {
            if (Access_basepair(xi, yj)==1) {
                int en_e = Access_E(a,b,i,j,x,y) + AU[xi][yj];
                int en_z = Access_Z(a,b,i,j,x,y);
                if (en_e == en_z) {
                    c = a;
                    i1 = i;
                    traced_ab = b;
                    traced_ij = j;
                    traced_xy = x;
                    goto LABEL1;
                }
            }
            for (c=b-1,traced_ab=0,traced_ij=0; c>=0; --c) {
                int pc = protein[c];
                int n_codon_c = n_codon[pc];
                for (i1 = 0; i1 < 3; ++i1) {
                    int lc = sigma(c,i1);
                    if (lc < 1 || rj - lc < 4) continue;
                    for (hx = 0; hx < n_codon_c; ++hx) {
                        
                        int hi = nucleotides[pc][hx][i1];
                        int en_e = Access_E(c,b,i1,j,hx,y) + AU[hi][yj];
                        int en_z;
                        if (i1 >= 1) {
                            if (a == c && x != hx) continue;
                            en_z = Access_Z(a,c,i,i1-1,x,hx);
                            if (zij == en_e + en_z) kj = hx;
                        } else {
                            int n_codon_cp = n_codon[protein[c-1]];
                            for (ky = 0; ky < n_codon_cp; ++ky) {
                                if (a == c-1 && x != ky) continue;
                                en_z = Access_Z(a,c-1,i,2,x,ky);
                                if (zij == en_e + en_z) {
                                    kj = ky;
                                    break;
                                }
                            }
                        }

                        if (zij == en_e + en_z) {
                            traced_ab = b;
                            traced_ij = j;
                            traced_xy = hx;

                            int tb, tj;
                            if (i1 >= 1) {
                                tb = c;
                                tj = i1-1;
                            }
                            else {
                                tb = c-1;
                                tj = 2;
                            }

                            sector[++s].a = a;
                            sector[s].i = i, sector[s].b = tb, sector[s].j = tj, sector[s].x = x, sector[s].y = kj, sector[s].ml = 0;

                            goto LABEL1;
                        }

                    }
                }
            }

            LABEL1:

            if (!traced_ab && !traced_ij) {
                exit(2);
            }
            a = c;
            i = i1;
            b = traced_ab;
            j = traced_ij;
            x = traced_xy;
            bp_bond[++t].i = sigma(a,i);
            bp_bond[t].j   = sigma(b,j);
            goto repeat1;
        }
        else {
            int hi = -1;
            int en_m;
            if (i <= 1) {
                en_m = Access_M(a, b, i + 1, j, x, y) + ML_BASE;
                if (en_m == zij) hi = x;
            }
            else {
                int n_codon_an = n_codon[protein[a + 1]];
                for (hx = 0; hx < n_codon_an; ++hx) {
                    if (a + 1 == b && hx != y) continue;
                    en_m = Access_M(a + 1, b, 0, j, hx, y) + ML_BASE;
                    if (en_m == zij) {
                        hi = hx;
                        break;
                    }
                }
            }
            if (hi != -1) {
                int ta, ti;
                if (i <= 1) {
                    ta = a;
                    ti = i + 1;
                }
                else {
                    ta = a + 1;
                    ti = 0;
                }

                sector[++s].a = ta;
                sector[s].i = ti, sector[s].b = b, sector[s].j = j, sector[s].x = hi, sector[s].y = y, sector[s].ml = ml;
                goto OUTLOOP;
            }

            if (zij == Access_E(a, b, i, j, x, y) + AU[xi][yj] + ML_intern) {
                bp_bond[++t].i = sigma(a, i);
                bp_bond[t].j = sigma(b, j);
                goto repeat1;
            }

            for (int lc = li + 5; lc <= rj-4; lc++) {
                c = lc / 3;
                i1 = lc % 3;
                int n_codon_c = n_codon[protein[c]];
                for (hx = 0; hx < n_codon_c; ++hx) {
                    int en_mb = inf;
                    if (i1 >= 1) {
                        en_mb = Access_M(a, c, i, i1 - 1, x, hx) + Access_M(c, b, i1, j, hx, y);
                        if (en_mb == zij) {
                            hi = hx;
                            kj = hx;
                        }
                    } else {
                        int n_codon_cp = n_codon[protein[c - 1]];
                        for (ky = 0; ky < n_codon_cp; ++ky) {
                            if (a == c - 1 && x != ky) continue;
                            en_mb = Access_M(a, c - 1, i, 2, x, ky) + Access_M(c, b, i1, j, hx, y);
                            if (en_mb == zij) {
                                hi = ky;
                                kj = hx;
                                break;
                            }
                        }
                    }

                    if (en_mb == zij) { //  compare(en_mb,zij)
                        int tb, tj;
                        if (i1 >= 1) {
                            tb = c;
                            tj = i1 - 1;
                        } else {
                            tb = c - 1;
                            tj = 2;
                        }

                        sector[++s].a = a;
                        sector[s].i = i, sector[s].b = tb, sector[s].j = tj, sector[s].x = x, sector[s].y = hi, sector[s].ml = ml;

                        sector[++s].a = c;
                        sector[s].i = i1, sector[s].b = b, sector[s].j = j, sector[s].x = kj, sector[s].y = y, sector[s].ml = ml;
                        goto OUTLOOP;
                    }
                }
//                }
            }
        }

        repeat1:
        eij = Access_E(a,b,i,j,x,y);
        rj = sigma(b,j);
        li = sigma(a,i);
        pa = protein[a];
        pb = protein[b];
        xi = nucleotides[pa][x][i];
        yj = nucleotides[pb][y][j];

        an_int = ava_nucle_p[index(a,x,i)];
        bp_int =  ava_nucle_m[index(b,y,j)];
        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;


        int pna = protein[a+1];
        int ppb = protein[b-1];
        int n_codon_an = n_codon[pna];
        n_codon_bp = n_codon[ppb];
        int l = rj-li;
        string h;

        if (l + 1 == 5) {
            int xi_, _yj, xi2_;
            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                } else {
                    xi2_ = nucleotides[pb][y][0];
                }
                h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};

                if (hairpinE.count(h) > 0) {
                    vector<int> values {xi_,xi2_,_yj};
                    if (eij == hairpinE[h]) {
                        assign(nucle_seq, values ,li+1);
                        goto OUTLOOP;
                    }
                }
            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};

                    if (hairpinE.count(h) > 0) {
                        vector<int> values {xi_,xi2_,_yj};
                        if (eij == hairpinE[h]) {
                            assign(nucle_seq, values ,li+1);
                            goto OUTLOOP;
                        }
                    }
                }
            }

            goto universal;
        }
        else if (l + 1 == 6) {
            int xi_, _yj, xi2_, _2yj;
            if (i == 0) {
                xi_ = nucleotides[pa][x][1];
                xi2_ = nucleotides[pa][x][2];
                _2yj = nucleotides[pb][y][0];
                _yj = nucleotides[pb][y][1];
                h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};
                vector<int> values {xi_,xi2_,_2yj,_yj};
                if (hairpinE.count(h) > 0) {
                    if (eij == hairpinE[h]) {
                        assign(nucle_seq,values,li+1);
                        goto OUTLOOP;
                    }
                }
            } else if (i == 1) {
                xi_ = nucleotides[pa][x][2];

                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi2_ = nucleotides[pna][x1][0];
                    _2yj = nucleotides[pna][x1][1];
                    _yj = nucleotides[pna][x1][2];
                    h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};
                    vector<int> values {xi_,xi2_,_2yj,_yj};
                    if (hairpinE.count(h) > 0) {
                        if (eij == hairpinE[h]) {
                            assign(nucle_seq,values,li+1);
                            goto OUTLOOP;
                        }
                    }
                }
            } else {
                _yj = nucleotides[pb][y][0];
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    _2yj = nucleotides[pna][x1][2];
                    h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};
                    vector<int> values {xi_,xi2_,_2yj,_yj};
                    if (hairpinE.count(h) > 0) {
                        if (eij == hairpinE[h]) {
                            assign(nucle_seq,values,li+1);
                            goto OUTLOOP;
                        }
                    }
                }
            }

            goto universal;
        }
        else if (l + 1 == 8) {
            int xi_, _yj, xi2_, _2yj, xi3_, _3yj;
            if (i <= 1) {
                xi_ = nucleotides[pa][x][i+1];
                _yj = nucleotides[pb][y][j-1];
                if (i == 0) {
                    xi2_ = nucleotides[pa][x][2];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi3_ = nucleotides[pna][x1][0];
                        _3yj = nucleotides[pna][x1][1];
                        _2yj = nucleotides[pna][x1][2];
                        h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
                        vector<int> values {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                        if (hairpinE.count(h) > 0) {
                            if (eij == hairpinE[h]) {
                                assign(nucle_seq,values,li+1);
                                goto OUTLOOP;
                            }
                        }
                    }
                } else {
                    _2yj = nucleotides[pb][y][0];
                    for (int x1 = 0; x1 < n_codon_an; ++x1) {
                        xi2_ = nucleotides[pna][x1][0];
                        xi3_ = nucleotides[pna][x1][1];
                        _3yj = nucleotides[pna][x1][2];
                        h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
                        vector<int> values {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                        if (hairpinE.count(h) > 0) {
                            if (eij == hairpinE[h]) {
                                assign(nucle_seq,values,li+1);
                                goto OUTLOOP;
                            }
                        }
                    }
                }
            } else {
                for (int x1 = 0; x1 < n_codon_an; ++x1) {
                    xi_ = nucleotides[pna][x1][0];
                    xi2_ = nucleotides[pna][x1][1];
                    xi3_ = nucleotides[pna][x1][2];
                    for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                        _3yj = nucleotides[ppb][y1][0];
                        _2yj = nucleotides[ppb][y1][1];
                        _yj = nucleotides[ppb][y1][2];
                        h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
                        vector<int> values {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                        if (hairpinE.count(h) > 0) {
                            if (eij == hairpinE[h]) {
                                assign(nucle_seq,values,li+1);
                                goto OUTLOOP;
                            }
                        }
                    }
                }
            }

            goto universal;
        }
        universal:
            for (int xi_ = 0; xi_ < 4; ++xi_) {
                if (((1 << xi_) & an_int) == 0) continue;
                for (int _yj = 0; _yj < 4; ++_yj) {
                    if (((1 << _yj) & bp_int) == 0) continue;
                    if (rj - li <= 4 && (li+1)/3==(rj-1)/3 && !rightCodon(li+1,rj-1,xi_,_yj)) continue;
                    if (eij == hairpin_loop(xi,yj,xi_,_yj,l-1)) {
                        nucle_seq[li+1] = xi_;
                        nucle_seq[rj-1] = _yj;
                        goto OUTLOOP;
                    }
                }
            }

        for (int lc = li + 1; lc <= min(rj-5, li+MAXLOOP+1); lc++) {
            int ll = lc - li;
            c = lc / 3;
            i1 = lc % 3;
            int min_ld = rj - li + lc - MAXLOOP - 2;
            if (min_ld < lc + 4) min_ld = lc + 4;
            for (int ld = rj-1; ld >= min_ld; ld--) {
                int lr = rj - ld;
                d = ld / 3;
                j1 = ld % 3;
                int pc = protein[c];
                int pd = protein[d];
                const int n_codon_c = n_codon[pc];
                const int n_codon_d = n_codon[pd];
                for (hx = 0; hx < n_codon_c; ++hx) {
                    if (c == a && hx != x) continue;
                    int hi = nucleotides[pc][hx][i1];
                    for (ky = 0; ky < n_codon_d; ++ky) {
                        if (d == b && ky != y) continue;
                        kj = nucleotides[pd][ky][j1];
                        int type = BP_pair[hi+1][kj+1];
                        if (type==0) continue;
                        int interior_energy = inf;
                        int internal_energy = eij - Access_E(c, d, i1, j1, hx, ky);
                        if (ll == 1 && lr == 1) {
                            interior_energy = stacking(xi,yj,hi,kj);
                            if (interior_energy == internal_energy) {
                                bp_bond[++t].i = lc;
                                bp_bond[t].j = ld;

                                nucle_seq[lc] = hi;
                                nucle_seq[ld] = kj;

                                a = c, b = d, i = i1, j = j1, x = hx, y = ky;
                                goto repeat1;
                            }

                        } else if (ll == 1 || lr == 1) {

                            if (ll == 1) {
                                interior_energy = bulge_loop(xi, yj, hi, kj, lr - 1);
                            } else {
                                interior_energy = bulge_loop(xi, yj, hi, kj, ll - 1);
                            }

                            if (interior_energy == internal_energy) {
                                bp_bond[++t].i = lc;
                                bp_bond[t].j = ld;

                                nucle_seq[lc] = hi;
                                nucle_seq[ld] = kj;

                                a = c, b = d, i = i1, j = j1, x = hx, y = ky;
                                goto repeat1;
                            }

                        } else {
                            int cp_int = ava_nucle_m[index(c,hx,i1)];
                            int dn_int = ava_nucle_p[index(d,ky,j1)];
                            for (int xi_ = 0; xi_ < 4; ++xi_) {
                                if (((1 << xi_) & an_int) == 0) continue;
                                for (int _hi = 0; _hi < 4; ++_hi) {
                                    if (((1 << _hi) & cp_int) == 0) continue;
                                    if (li + 1 == lc - 1 && xi_ != _hi) continue;
                                    if (lc - li <= 4 && (li+1)/3==(lc-1)/3 && !rightCodon(li+1,lc-1,xi_,_hi)) continue;
                                    if (lc -1 - li <= 2) {
                                        if (li/3==(lc-1)/3 && !rightCodon(li, lc-1,xi,_hi)) continue;
                                        if ((li+1)/3==lc/3 && !rightCodon(li+1, lc,xi_,hi)) continue;
                                    }
                                    for (int kj_ = 0; kj_ < 4; ++kj_) {
                                        if (((1 << kj_) & dn_int) == 0) continue;
                                        for (int _yj = 0; _yj < 4; ++_yj) {
                                            if (((1 << _yj) & bp_int) == 0) continue;
                                            if (ld + 1 == rj - 1 && kj_ != _yj) continue;
                                            if (rj - ld <= 4 && (ld+1)/3==(rj-1)/3 && !rightCodon(ld+1,rj-1,kj_,_yj)) continue;
                                            if (rj -1 - ld <= 2) {
                                                if (ld/3==(rj-1)/3 && !rightCodon(ld, rj-1,kj,_yj)) continue;
                                                if ((ld+1)/3==rj/3 && !rightCodon(ld+1, rj,kj_,yj)) continue;
                                            }

                                            interior_energy = interior_loop(xi, yj, hi, kj, xi_, _yj, _hi, kj_,
                                                                            ll - 1, lr - 1);
                                            if (interior_energy == internal_energy) {
                                                bp_bond[++t].i = lc;
                                                bp_bond[t].j = ld;

                                                nucle_seq[lc] = hi;
                                                nucle_seq[ld] = kj;

                                                nucle_seq[li + 1] = xi_;
                                                nucle_seq[rj - 1] = _yj;
                                                nucle_seq[lc - 1] = _hi;
                                                nucle_seq[ld + 1] = kj_;

                                                a = c, b = d, i = i1, j = j1, x = hx, y = ky;
                                                goto repeat1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        sector[s+1].ml = sector[s+2].ml = 1;
        int en = eij - AU[xi][yj] - ML_intern - ML_closing;
        int sa,sb,si,sj,sx,sy,sh,sk;

        for (c = a+1; c < b-1; ++c) {
            for (i1 = 0; i1 < 3; ++i1) {
                int lc = sigma(c,i1);
                if (lc - li < 4 || rj - lc < 4) continue;
                int n_codon_c = n_codon[protein[c]];
                for (hx = 0; hx < n_codon_c; ++hx) {
                    int en_mb;
                    if (i <= 1 && j >= 1) {
                        if (i1 >= 1) {
                            en_mb = Access_M(a,c,i+1,i1-1,x,hx) + Access_M(c,b,i1,j-1,hx,y);
                            if (en == en_mb) {
                                sx = x, sh = hx, sk = hx, sy = y;
                                goto LABEL2;
                            }
                        }
                        else {
                            int n_codon_cp = n_codon[protein[c-1]];
                            for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                if (a == c-1 && x != x1) continue;
                                en_mb = Access_M(a,c-1,i+1,2,x,x1) + Access_M(c,b,i1,j-1,hx,y);
                                if (en == en_mb) {
                                    sx = x, sh = x1, sk = hx, sy = y;
                                    goto LABEL2;
                                }
                            }
                        }
                    }
                    else if (i <= 1) {
                        for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                            if (i1 >= 1) {
                                if ((a == c && x != hx) || (c == b-1 && hx != y1)) continue;
                                en_mb = Access_M(a,c,i+1,i1-1,x,hx) + Access_M(c,b-1,i1,2,hx,y1);
                                if (en == en_mb) {
                                    sx = x, sh = hx, sk = hx, sy = y1;
                                    goto LABEL2;
                                }
                            }
                            else {
                                int n_codon_cp = n_codon[protein[c-1]];
                                for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                    if ((a == c-1 && x != x1) || (c == b-1 && hx != y1)) continue;
                                    en_mb = Access_M(a,c-1,i+1,2,x,x1) + Access_M(c,b-1,i1,2,hx,y1);
                                    if (en == en_mb) {
                                        sx = x, sh = x1, sk = hx, sy = y1;
                                        goto LABEL2;
                                    }
                                }
                            }
                        }
                    }
                    else if (j >= 1) {
                        for (int x2 = 0; x2 < n_codon_an; ++x2) {
                            if (i1 >= 1) {
                                en_mb = Access_M(a+1,c,0,i1-1,x2,hx) + Access_M(c,b,i1,j-1,hx,y);
                                if (en == en_mb) {
                                    sx = x2, sh = hx, sk = hx, sy = y;
                                    goto LABEL2;
                                }
                            }
                            else {
                                int n_codon_cp = n_codon[protein[c-1]];
                                for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                    en_mb = Access_M(a+1,c-1,0,2,x2,x1) + Access_M(c,b,i1,j-1,hx,y);
                                    if (en == en_mb) {
                                        sx = x2, sh = x1, sk = hx, sy = y;
                                        goto LABEL2;
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (int x2 = 0; x2 < n_codon_an; ++x2) {
                            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                                if (i1 >= 1) {
                                    en_mb = Access_M(a+1,c,0,i1-1,x2,hx) + Access_M(c,b-1,i1,2,hx,y1);
                                    if (en == en_mb) {
                                        sx = x2, sh = hx, sk = hx, sy = y1;
                                        goto LABEL2;
                                    }
                                }
                                else {
                                    int n_codon_cp = n_codon[protein[c-1]];
                                    for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                        en_mb = Access_M(a+1,c-1,0,2,x2,x1) + Access_M(c,b-1,i1,2,hx,y1);
                                        if (en == en_mb) {
                                            sx = x2, sh = x1, sk = hx, sy = y1;
                                            goto LABEL2;
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }

        LABEL2:

        if (rj - sigma(c,i1) > 4) {
            int tb, tj;
            if (i <= 1) {
                sa = a;
                si = i+1;
            }
            else {
                sa = a+1;
                si = 0;
            }
            if (j >= 1) {
                sb = b;
                sj = j-1;
            }
            else {
                sb = b-1;
                sj = 2;
            }
            if (i1 >= 1) {
                tb = c;
                tj = i1-1;
            } else {
                tb = c-1;
                tj = 2;
            }

            sector[++s].a = sa;
            sector[s].i = si, sector[s].b = tb, sector[s].j = tj, sector[s].x = sx, sector[s].y = sh;


            sector[++s].a = c;
            sector[s].i = i1, sector[s].b = sb, sector[s].j = sj, sector[s].x = sk, sector[s].y = sy;
        } else {
            cout << "backtracking failed in repeat: " << a << " " << b << " " << i << " " << j << " " << x << " " << y <<  endl;
        }
    }
    bp_bond[0].i = t;

}

// absolute minimal value
void Zuker::traceback_B() {
    int s = 0;
    int t = 0;
    int bt = 0;
    int oij, oi, eij; //oj
    int c,d,i1,j1,hx,ky;
    int traced_ab, traced_ij, traced_xy;
    sector[++s].a = 0;
    sector[s].b = n-1;
    sector[s].i = 0;
    sector[s].j = 2;
    sector[s].x = minX;
    sector[s].y = minY;
    sector[s].ml = 0;

    OUTLOOP:
    while (s > 0) {
        int an_int, bp_int;
        int a = sector[s].a;
        int b = sector[s].b;
        int i  = sector[s].i;
        int j  = sector[s].j;
        int x  = sector[s].x;
        int y  = sector[s].y;
        int ml = sector[s--].ml;
        int idx = index(a,b,i,j,x,y);
        int pa = protein[a];
        int pb = protein[b];
        int li = sigma(a,i);
        int rj = sigma(b,j);

        int xi = nucleotides[pa][x][i];
        int yj = nucleotides[pb][y][j];
        int n_codon_bp = n_codon[protein[b-1]];
        if (a == b && i == j) {
            break;
        }

        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;

        if (ml==2) {
            bp_bond[++t].i = li;
            bp_bond[t].j   = rj;
            goto repeat;
        }


        if (ml == 0) {
            int energy;
            bt = Access_OB(idx);
            oij = Access_Z(a,b,i,j,x,y);
            switch (bt) {
                int kj, c1, i1_;
                case -1:
                    oi = Access_E(a,b,i,j,x,y) + AU[xi][yj];
                    bp_bond[++t].i = sigma(a,i);
                    bp_bond[t].j   = sigma(b,j);
                    goto repeat;
                    break;
                case -2:
                    oi = Access_Z(a,b,i,j-1,x,y);
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i, sector[s].j = j-1, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -3:
                    energy = inf;
                    for (ky = 0; ky < n_codon_bp; ++ky) {
                        if (a == b-1 && x != ky) continue;
                        oi = Access_Z(a,b-1,i,2,x,ky);

                        if (energy > oi) {
                            energy = oi;
                            kj = ky;
                        }
                    }

                    if (compare(oij, energy)) {
                        sector[++s].a = a;
                        sector[s].b = b - 1, sector[s].i = i, sector[s].j = 2, sector[s].x = x, sector[s].y = kj, sector[s].ml = ml;
                        goto OUTLOOP;
                    }
                    break;
                case -4:
                    energy = inf;
                    for (c=b-1,traced_ab=0,traced_ij=0; c>=0; --c) {
                        int pc = protein[c];
                        int n_codon_c = n_codon[pc];
                        for (i1 = 0; i1 < 3; ++i1) {
                            int lc = sigma(c,i1);
                            if (lc < 1 || rj - lc < 4) continue;
                            for (hx = 0; hx < n_codon_c; ++hx) {

                                int hi = nucleotides[pc][hx][i1];
                                int en_e = Access_E(c,b,i1,j,hx,y) + AU[hi][yj];
                                int en_z;
                                if (i1 >= 1) {
                                    if (a == c && x != hx) continue;

                                    en_z = Access_Z(a,c,i,i1-1,x,hx) + en_e;
                                    if (energy > en_z) {
                                        energy = en_z;
                                        kj = hx;
                                        c1 = c, i1_ = i1;
                                    }
                                } else {
                                    int n_codon_cp = n_codon[protein[c-1]];
                                    for (ky = 0; ky < n_codon_cp; ++ky) {
                                        if (a == c-1 && x != ky) continue;

                                        en_z = Access_Z(a,c-1,i,2,x,ky) + en_e;
                                        if (energy > en_z) {
                                            energy = en_z;
                                            kj = ky;
                                            c1 = c, i1_ = i1;
                                        }
                                    }
                                }

//                                {
                                if (compare(oij, energy)) { //   ,  oij == en_e + en_z
                                    traced_ab = b;
                                    traced_ij = j;
                                    traced_xy = hx;

                                    int tb, tj;
                                    if (i1_ >= 1) {
                                        tb = c1;
                                        tj = i1_-1;
                                    }
                                    else {
                                        tb = c1-1;
                                        tj = 2;
                                    }

                                    sector[++s].a = a;
                                    sector[s].i = i, sector[s].b = tb, sector[s].j = tj, sector[s].x = x, sector[s].y = kj, sector[s].ml = 0;
                                    if (!traced_ab && !traced_ij) {
                                        cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << oij << endl;
                                        exit(2);
                                    }

                                    a = c1;
                                    i = i1_;
                                    b = traced_ab;
                                    j = traced_ij;
                                    x = traced_xy;
                                    bp_bond[++t].i = sigma(a,i);
                                    bp_bond[t].j   = sigma(b,j);
                                    goto repeat;
                                }

                            }
                        }
                    }
                default:
                    cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << oij << " " << bt << endl;
                    exit(1);
            }
        }
        else {
            bt = Access_MB(idx);
            oij = Access_M(a,b,i,j,x,y);
            int n_codon_an = n_codon[protein[a + 1]];
            switch (bt) {
                int hi, kj;
                int c1, i1_;
                int energy;
                case -1:
                    bp_bond[++t].i = sigma(a, i);
                    bp_bond[t].j = sigma(b, j);
                    goto repeat;
                    break;
                case -2:
                    oi = Access_M(a, b, i + 1, j, x, y) + ML_BASE;
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i+1, sector[s].j = j, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -3:
                    energy = inf;
                    for (hx = 0; hx < n_codon_an; ++hx) {
                        if (a + 1 == b && hx != y) continue;
                        oi = Access_M(a + 1, b, 0, j, hx, y)  + ML_BASE;
                        if (energy > oi) {
                            energy = oi;
                            hi = hx;
                        }
//
                    }
                    if (compare(energy, oij)) {
                        sector[++s].a = a + 1;
                        sector[s].b = b, sector[s].i = 0, sector[s].j = j, sector[s].x = hi, sector[s].y = y, sector[s].ml = ml;
                        goto OUTLOOP;
                    }
                    break;
                case -4:
                    oi = Access_M(a,b,i,j-1,x,y) +  ML_BASE;
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i, sector[s].j = j-1, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -5:
                    energy = inf;
                    for (ky = 0; ky < n_codon_bp; ++ky) {
                        if (a == b-1 && x != ky) continue;
                        oi = Access_M(a,b-1,i,2,x,ky) +  ML_BASE;
                        if (energy > oi) {
                            energy = oi;
                            kj = ky;
                        }
                    }
                    if (compare(energy, oij)) {
                        sector[++s].a = a;
                        sector[s].b = b - 1, sector[s].i = i, sector[s].j = 2, sector[s].x = x, sector[s].y = kj, sector[s].ml = ml;
                        goto OUTLOOP;
                    }
                    break;
                case -6:
                    energy = inf;
                    for (int lc = li + 5; lc <= rj-4; lc++) {
                        c = lc / 3;
                        i1 = lc % 3;
                        int pc = protein[c];
                        int n_codon_c = n_codon[pc];
                        for (hx = 0; hx < n_codon_c; ++hx) {
                            int en_mb;
                            if (i1 >= 1) {
                                en_mb = Access_M(a, c, i, i1 - 1, x, hx) + Access_M(c, b, i1, j, hx, y);
                                if (energy > en_mb) {
                                    energy = en_mb;
                                    hi = hx, kj = hx;
                                    c1 = c, i1_ = i1;
                                }
                            } else {
                                int n_codon_cp = n_codon[protein[c - 1]];
                                for (ky = 0; ky < n_codon_cp; ++ky) {
                                    if (a == c - 1 && x != ky) continue;

                                    en_mb = Access_M(a, c - 1, i, 2, x, ky) + Access_M(c, b, 0, j, hx, y);
                                    if (energy > en_mb) {
                                        energy = en_mb;
                                        hi = ky, kj = hx;
                                        c1 = c, i1_ = i1;
                                    }
                                }
                            }
//                            {
                            if (compare(energy, oij)) { //, en_mb == oij,
                                int tb, tj;
                                if (i1_ >= 1) {
                                    tb = c1;
                                    tj = i1_ - 1;
                                } else {
                                    tb = c1 - 1;
                                    tj = 2;
                                }

                                sector[++s].a = a;
                                sector[s].i = i, sector[s].b = tb, sector[s].j = tj, sector[s].x = x, sector[s].y = hi, sector[s].ml = ml;

                                sector[++s].a = c1;
                                sector[s].i = i1_, sector[s].b = b, sector[s].j = j, sector[s].x = kj, sector[s].y = y, sector[s].ml = ml;
                                goto OUTLOOP;
                            }
                        }
                    }
                    break;
                default:
                    cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << oij << endl;
                    exit(3);
            }
        }

        repeat:
        eij = Access_E(a,b,i,j,x,y);
        rj = sigma(b,j);
        li = sigma(a,i);
        pa = protein[a];
        pb = protein[b];
        xi = nucleotides[pa][x][i];
        yj = nucleotides[pb][y][j];

        if (a > b || (a == b && i >= j)) {
            break;
        }

        an_int = ava_nucle_p[index(a,x,i)];
        bp_int = ava_nucle_m[index(b,y,j)];
        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;
        int pna = protein[a+1];
        int ppb = protein[b-1];
        int n_codon_an = n_codon[pna];
        n_codon_bp = n_codon[ppb];
        int l = rj-li;
        idx = index(a,b,i,j,x,y);
        bt = Access_EB(idx);
        string h;
        vector<int> values;
        switch (bt) {
            int kj, lc, ld, pc, pd, n_codon_c, n_codon_d, min_ld;
            int lc1, ld1, hi1, kj1, hx1, ky1, c1, d1, i1_, j1_;
            int xi1_, _yj1;
            int interior_energy, energy; //internal_energy

            case 0:
                break;
            case -1:
                energy = inf;
                int uni;
                uni = 0;
                if (l + 1 == 5) {
                    int xi_, _yj, xi2_;
                    if (i <= 1) {
                        xi_ = nucleotides[pa][x][i+1];
                        _yj = nucleotides[pb][y][j-1];
                        if (i == 0) {
                            xi2_ = nucleotides[pa][x][2];
                        } else {
                            xi2_ = nucleotides[pb][y][0];
                        }
                        h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};
                        if (hairpinE.count(h) > 0) {

                            interior_energy = hairpinE[h];
                            if (energy > interior_energy) {
                                energy = interior_energy;
                                values =  {xi_,xi2_,_yj};
                            }
                        }
                    } else {
                        for (int x1 = 0; x1 < n_codon_an; ++x1) {
                            xi_ = nucleotides[pna][x1][0];
                            xi2_ = nucleotides[pna][x1][1];
                            _yj = nucleotides[pna][x1][2];
                            h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_yj],to_char[yj]};
                            if (hairpinE.count(h) > 0) {
                                interior_energy = hairpinE[h];
                                if (energy > interior_energy) {
                                    energy = interior_energy;
                                    values =  {xi_,xi2_,_yj};
                                }

                            }
                        }
                    }
                    goto universal;
                }
                else if (l + 1 == 6) {
                    int xi_, _yj, xi2_, _2yj;
                    if (i == 0) {
                        xi_ = nucleotides[pa][x][1];
                        xi2_ = nucleotides[pa][x][2];
                        _2yj = nucleotides[pb][y][0];
                        _yj = nucleotides[pb][y][1];
                        h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                        if (hairpinE.count(h) > 0) {
                            interior_energy = hairpinE[h];
                            if (energy > interior_energy) {
                                energy = interior_energy;
                                values = {xi_,xi2_,_2yj,_yj};
                            }
                        }
                    } else if (i == 1) {
                        xi_ = nucleotides[pa][x][2];
                        for (int x1 = 0; x1 < n_codon_an; ++x1) {
                            xi2_ = nucleotides[pna][x1][0];
                            _2yj = nucleotides[pna][x1][1];
                            _yj = nucleotides[pna][x1][2];
                            h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};

                            if (hairpinE.count(h) > 0) {
                                interior_energy = hairpinE[h] ;
                                if (energy > interior_energy) {
                                    energy = interior_energy;
                                    values = {xi_,xi2_,_2yj,_yj};
                                }
                            }
                        }
                    } else {
                        _yj = nucleotides[pb][y][0];
                        for (int x1 = 0; x1 < n_codon_an; ++x1) {
                            xi_ = nucleotides[pna][x1][0];
                            xi2_ = nucleotides[pna][x1][1];
                            _2yj = nucleotides[pna][x1][2];
                            h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[_2yj],to_char[_yj],to_char[yj]};
                            if (hairpinE.count(h) > 0) {
                                interior_energy = hairpinE[h];
                                if (energy > interior_energy) {
                                    energy = interior_energy;
                                    values = {xi_,xi2_,_2yj,_yj};
                                }
                            }
                        }
                    }

                    goto universal;
                }
                else if (l + 1 == 8) {
                    int xi_, _yj, xi2_, _2yj, xi3_, _3yj;
                    if (i <= 1) {
                        xi_ = nucleotides[pa][x][i+1];
                        _yj = nucleotides[pb][y][j-1];
                        if (i == 0) {
                            xi2_ = nucleotides[pa][x][2];
                            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                xi3_ = nucleotides[pna][x1][0];
                                _3yj = nucleotides[pna][x1][1];
                                _2yj = nucleotides[pna][x1][2];
                                h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};

                                if (hairpinE.count(h) > 0) {
                                    interior_energy = hairpinE[h] ;
                                    if (energy > interior_energy) {
                                        energy = interior_energy;
                                        values = {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                                    }
                                }
                            }
                        } else {
                            _2yj = nucleotides[pb][y][0];
                            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                                xi2_ = nucleotides[pna][x1][0];
                                xi3_ = nucleotides[pna][x1][1];
                                _3yj = nucleotides[pna][x1][2];
                                h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};
                                if (hairpinE.count(h) > 0) {
                                    interior_energy = hairpinE[h];
                                    if (energy > interior_energy) {
                                        energy = interior_energy;
                                        values = {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                                    }
                                }
                            }
                        }
                    } else {
                        for (int x1 = 0; x1 < n_codon_an; ++x1) {
                            xi_ = nucleotides[pna][x1][0];
                            xi2_ = nucleotides[pna][x1][1];
                            xi3_ = nucleotides[pna][x1][2];
                            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                                _3yj = nucleotides[ppb][y1][0];
                                _2yj = nucleotides[ppb][y1][1];
                                _yj = nucleotides[ppb][y1][2];
                                h = {to_char[xi],to_char[xi_],to_char[xi2_],to_char[xi3_],to_char[_3yj],to_char[_2yj],to_char[_yj],to_char[yj]};

                                if (hairpinE.count(h) > 0) {
                                    interior_energy = hairpinE[h];
                                    if (energy > interior_energy) {
                                        energy = interior_energy;
                                        values = {xi_,xi2_,xi3_,_3yj,_2yj,_yj};
                                    }
                                }
                            }
                        }
                    }
                    goto universal;
                }

            universal:
                for (int xi_ = 0; xi_ < 4; ++xi_) {
                    if (((1 << xi_) & an_int) == 0) continue;
                    for (int _yj = 0; _yj < 4; ++_yj) {
                        if (((1 << _yj) & bp_int) == 0) continue;
                        if (rj - li <= 4 && (li+1)/3==(rj-1)/3 && !rightCodon(li+1,rj-1,xi_,_yj)) continue;
                        interior_energy = hairpin_loop(xi,yj,xi_,_yj,l-1) ;
                        if (energy > interior_energy) {
                            energy = interior_energy;
                            uni = 1;
                            xi1_ = xi_, _yj1 = _yj;
                        }
                    }
                }

                if (compare(eij, energy)) {
                    if (uni == 0) {
                        assign(nucle_seq,values,li+1);
                        goto OUTLOOP;
                    } else {
                        nucle_seq[li+1] = xi1_;
                        nucle_seq[rj-1] = _yj1;
                        goto OUTLOOP;
                    }
                }
                break;
            case -3:
                energy = inf;
                sector[s+1].ml = sector[s+2].ml = 1;
                int en;
                en = eij - (AU[xi][yj] + ML_intern + ML_closing);

                int sa,sb,si,sj,sx,sy,sh,sk;

                for (c = a+1; c < b-1; ++c) {
                    for (i1 = 0; i1 < 3; ++i1) {
                        lc = sigma(c,i1);
                        if (lc - li < 4 || rj - lc < 4) continue;
                        pc = protein[c];
                        n_codon_c = n_codon[pc];
                        for (hx = 0; hx < n_codon_c; ++hx) {
                            int en_mb;
                            if (i <= 1 && j >= 1) {
                                if (i1 >= 1) {
                                    en_mb = Access_M(a,c,i+1,i1-1,x,hx) + Access_M(c,b,i1,j-1,hx,y) ;
                                    if (energy > en_mb) {
                                        energy = en_mb;
                                        sx = x, sh = hx, sk = hx, sy = y;
                                        c1 = c, i1_ = i1;
                                    }
                                }
                                else {
                                    int n_codon_cp = n_codon[protein[c-1]];
                                    for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                        if (a == c-1 && x != x1) continue;

                                        en_mb = Access_M(a,c-1,i+1,2,x,x1) + Access_M(c,b,i1,j-1,hx,y);
                                        if (energy > en_mb) {
                                            energy = en_mb;
                                            sx = x, sh = x1, sk = hx, sy = y;
                                            c1 = c, i1_ = i1;
                                        }
                                    }
                                }
                            }
                            else if (i <= 1) {
                                for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                                    if (i1 >= 1) {
                                        en_mb = Access_M(a,c,i+1,i1-1,x,hx) + Access_M(c,b-1,i1,2,hx,y1) ;
                                        if (energy > en_mb) {
                                            energy = en_mb;
                                            sx = x, sh = hx, sk = hx, sy = y1;
                                            c1 = c, i1_ = i1;
                                        }
                                    }
                                    else {
                                        int n_codon_cp = n_codon[protein[c-1]];
                                        for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                            if (a == c-1 && x != x1) continue;

                                            en_mb = Access_M(a,c-1,i+1,2,x,x1) + Access_M(c,b-1,i1,2,hx,y1);
                                            if (energy > en_mb) {
                                                energy = en_mb;
                                                sx = x, sh = x1, sk = hx, sy = y1;
                                                c1 = c, i1_ = i1;
                                            }
                                        }
                                    }
                                }
                            }
                            else if (j >= 1) {
                                for (int x2 = 0; x2 < n_codon_an; ++x2) {
                                    if (i1 >= 1) {
                                        if (a == c-1 && x2 != hx) continue;

                                        en_mb = Access_M(a+1,c,0,i1-1,x2,hx) + Access_M(c,b,i1,j-1,hx,y) ;
                                        if (energy > en_mb) {
                                            energy = en_mb;
                                            sx = x2, sh = hx, sk = hx, sy = y;
                                            c1 = c, i1_ = i1;
                                        }
                                    }
                                    else {
                                        int n_codon_cp = n_codon[protein[c-1]];
                                        for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                            if (a+1 == c-1 && x2 != x1) continue;

                                            en_mb = Access_M(a+1,c-1,0,2,x2,x1) + Access_M(c,b,i1,j-1,hx,y);
                                            if (energy > en_mb) {
                                                energy = en_mb;
                                                sx = x2, sh = x1, sk = hx, sy = y;
                                                c1 = c, i1_ = i1;
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                for (int x2 = 0; x2 < n_codon_an; ++x2) {
                                    for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                                        if (i1 >= 1) {
                                            if (a == c-1 && x2 != hx) continue;

                                            en_mb = Access_M(a+1,c,0,i1-1,x2,hx) + Access_M(c,b-1,i1,2,hx,y1); //+ (lambda-1)*(codon_cai[pa][x] + codon_cai[pb][y])
                                            if (energy > en_mb) {
                                                energy = en_mb;
                                                sx = x2, sh = hx, sk = hx, sy = y1;
                                                c1 = c, i1_ = i1;
                                            }
                                        }
                                        else {
                                            int n_codon_cp = n_codon[protein[c-1]];
                                            for (int x1 = 0; x1 < n_codon_cp; ++x1) {
                                                if (a+1 == c-1 && x2 != x1) continue;

                                                en_mb = Access_M(a+1,c-1,0,2,x2,x1) + Access_M(c,b-1,i1,2,hx,y1) ; //+ (lambda-1)*(codon_cai[pa][x] + codon_cai[pb][y])
                                                if (energy > en_mb) {
                                                    energy = en_mb;
                                                    sx = x2, sh = x1, sk = hx, sy = y1;
                                                    c1 = c, i1_ = i1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

//            LABEL2:
                if (compare(en, energy)) {

                    if (rj - sigma(c1, i1_) > 4) {
                        int tb, tj;
                        if (i <= 1) {
                            sa = a;
                            si = i + 1;
                        } else {
                            sa = a + 1;
                            si = 0;
                        }
                        if (j >= 1) {
                            sb = b;
                            sj = j - 1;
                        } else {
                            sb = b - 1;
                            sj = 2;
                        }
                        if (i1_ >= 1) {
                            tb = c1;
                            tj = i1_ - 1;
                        } else {
                            tb = c1 - 1;
                            tj = 2;
                        }

                        sector[++s].a = sa;
                        sector[s].i = si, sector[s].b = tb, sector[s].j = tj, sector[s].x = sx, sector[s].y = sh;


                        sector[++s].a = c1;
                        sector[s].i = i1_, sector[s].b = sb, sector[s].j = sj, sector[s].x = sk, sector[s].y = sy;
                        break;
                    }
                }
                break;
            case -4:
                energy = inf;
                lc = li + 1;
                c = lc / 3;
                i1 = lc % 3;
                ld = rj - 1;
                d = ld / 3;
                j1 = ld % 3;
                pc = protein[c];
                pd = protein[d];
                n_codon_c = n_codon[pc];
                n_codon_d = n_codon[pd];
                for (hx = 0; hx < n_codon_c; ++hx) {
                    if (c == a && hx != x) continue;
                    for (ky = 0; ky < n_codon_d; ++ky) {
                        if (d == b && ky != y) continue;
                        int hi = nucleotides[pc][hx][i1];
                        kj = nucleotides[pd][ky][j1];
                        int type2 = BP_pair[hi+1][kj+1];
                        if (type2==0) {
                            continue;
                        }

                        interior_energy = Access_E(c, d, i1, j1, hx, ky) + stacking(xi, yj, hi, kj) ;
                        if (energy > interior_energy) {
                            energy = interior_energy;
                            lc1 = lc, ld1 = ld, hi1 = hi, kj1 = kj, hx1 = hx, ky1 = ky, c1 = c, d1 = d, i1_ = i1, j1_ = j1;

                        }


                    }
                }
                if (compare(energy,eij)) { // ,  interior_energy == internal_energy
//                {
                    bp_bond[++t].i = lc1;
                    bp_bond[t].j = ld1;

                    nucle_seq[lc1] = hi1;
                    nucle_seq[ld1] = kj1;
                    a = c1, b = d1, i = i1_, j = j1_, x = hx1, y = ky1;
                    goto repeat;
                }
                break;
            case -5:

                energy = inf;

                for (lc = li + 1; lc <= min(rj-5, li+MAXLOOP+1); lc++) {
                    int ll = lc - li;
                    c = lc / 3;
                    i1 = lc % 3;
                    min_ld = rj - li + lc - MAXLOOP - 2;
                    if (min_ld < lc + 4) min_ld = lc + 4;
                    for (ld = rj-1; ld >= min_ld; ld--) {
                        int lr = rj - ld;
                        d = ld / 3;
                        j1 = ld % 3;
                        pc = protein[c];
                        pd = protein[d];
                        n_codon_c = n_codon[pc];
                        n_codon_d = n_codon[pd];
                        for (hx = 0; hx < n_codon_c; ++hx) {
                            if (c == a && hx != x) continue;
                            for (ky = 0; ky < n_codon_d; ++ky) {
                                if (d == b && ky != y) continue;
                                int hi = nucleotides[pc][hx][i1];
                                kj = nucleotides[pd][ky][j1];
                                int type2 = BP_pair[hi+1][kj+1];
                                if (type2==0) {
                                    continue;
                                }


                                if ((ll == 1 && lr > 1) || (lr == 1 && ll > 1)) {

                                    interior_energy = Access_E(c, d, i1, j1, hx, ky) + bulge_loop(xi, yj, hi, kj, max(ll,lr) - 1) ;
                                    if (energy > interior_energy) {
                                        energy = interior_energy;
                                        lc1 = lc, ld1 = ld, hi1 = hi, kj1 = kj, hx1 = hx, ky1 = ky, c1 = c, d1 = d, i1_ = i1, j1_ = j1;

                                    }


                                }
                            }
                        }
                    }
                }
                if (compare(energy,eij)) { //,  interior_energy == internal_energy
//                {
                    bp_bond[++t].i = lc1;
                    bp_bond[t].j = ld1;

                    nucle_seq[lc1] = hi1;
                    nucle_seq[ld1] = kj1;


                    a = c1, b = d1, i = i1_, j = j1_, x = hx1, y = ky1;
                    goto repeat;
                }
                break;
            case -6:
                int _hi1, kj1_;
                energy = inf;
                for (lc = li + 1; lc <= min(rj-5, li+MAXLOOP+1); lc++) {
                    int ll = lc - li;
                    c = lc / 3;
                    i1 = lc % 3;
                    min_ld = rj - li + lc - MAXLOOP - 2;
                    if (min_ld < lc + 4) min_ld = lc + 4;
                    for (ld = rj-1; ld >= min_ld; ld--) {
                        int lr = rj - ld;
                        d = ld / 3;
                        j1 = ld % 3;
                        pc = protein[c];
                        pd = protein[d];
                        n_codon_c = n_codon[pc];
                        n_codon_d = n_codon[pd];
                        for (hx = 0; hx < n_codon_c; ++hx) {
                            if (c == a && hx != x) continue;
                            for (ky = 0; ky < n_codon_d; ++ky) {
                                if (d == b && ky != y) continue;
                                int hi = nucleotides[pc][hx][i1];
                                kj = nucleotides[pd][ky][j1];
                                int type2 = BP_pair[hi+1][kj+1];
                                if (type2==0) {
                                    continue;
                                }

                                if (ll > 1 && lr > 1) {
                                    int cp_int = ava_nucle_m[index(c,hx,i1)];
                                    int dn_int = ava_nucle_p[index(d,ky,j1)];
                                    for (int xi_ = 0; xi_ < 4; ++xi_) {
                                        if (((1 << xi_) & an_int) == 0) continue;
                                        for (int _hi = 0; _hi < 4; ++_hi) {
                                            if (((1 << _hi) & cp_int) == 0) continue;
                                            if (li + 1 == lc - 1 && xi_ != _hi) continue;
                                            if (lc - li <= 4 && (li+1)/3==(lc-1)/3 && !rightCodon(li+1,lc-1,xi_,_hi)) continue;
                                            if (lc -1 - li <= 2) {
                                                if (li/3==(lc-1)/3 && !rightCodon(li, lc-1,xi,_hi)) continue;
                                                if ((li+1)/3==lc/3 && !rightCodon(li+1, lc,xi_,hi)) continue;
                                            }
                                            for (int kj_ = 0; kj_ < 4; ++kj_) {
                                                if (((1 << kj_) & dn_int) == 0) continue;
                                                for (int _yj = 0; _yj < 4; ++_yj) {
                                                    if (((1 << _yj) & bp_int) == 0) continue;
                                                    if (ld + 1 == rj - 1 && kj_ != _yj) continue;
                                                    if (rj - ld <= 4 && (ld+1)/3==(rj-1)/3 && !rightCodon(ld+1,rj-1,kj_,_yj)) continue;
                                                    if (rj -1 - ld <= 2) {
                                                        if (ld/3==(rj-1)/3 && !rightCodon(ld, rj-1,kj,_yj)) continue;
                                                        if ((ld+1)/3==rj/3 && !rightCodon(ld+1, rj,kj_,yj)) continue;
                                                    }

                                                    interior_energy = Access_E(c, d, i1, j1, hx, ky) + interior_loop(xi, yj, hi, kj, xi_, _yj, _hi, kj_,
                                                                                                                              ll - 1, lr - 1);
                                                    if (energy > interior_energy) {
                                                        energy = interior_energy;
                                                        lc1 = lc, ld1 = ld, hi1 = hi, kj1 = kj, hx1 = hx, ky1 = ky, c1 = c, d1 = d, i1_ = i1, j1_ = j1;
                                                        xi1_ = xi_, _yj1 = _yj, _hi1 = _hi, kj1_ = kj_;
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (compare(energy,eij)) { // ,  interior_energy == internal_energy
//                {
                    bp_bond[++t].i = lc1;
                    bp_bond[t].j = ld1;

                    nucle_seq[lc1] = hi1;
                    nucle_seq[ld1] = kj1;

                    nucle_seq[li + 1] = xi1_;
                    nucle_seq[rj - 1] = _yj1;
                    nucle_seq[lc1 - 1] = _hi1;
                    nucle_seq[ld1 + 1] = kj1_;


                    a = c1, b = d1, i = i1_, j = j1_, x = hx1, y = ky1;
                    goto repeat;
                }
                break;
            default:
                cout << "backtracking failed in repeat: " << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << en << " " << bt << " " << endl;
                exit(3);
        }
    }
    bp_bond[0].i = t;


}

void Zuker::traceback_B2(double lambda) {

    int s = 0;
    int t = 0;
    int bt = 0;
    double oij, oi, eij;
    int c,d,i1,j1,hx,ky;
    sector[++s].a = 0;
    sector[s].b = n-1;
    sector[s].i = 0;
    sector[s].j = 2;
    sector[s].x = minX;
    sector[s].y = minY;
    sector[s].ml = 0;

    codon_selection[0] = minX;
    codon_selection[n-1] = minY;

    OUTLOOP:
    while (s > 0) {
        int a = sector[s].a;
        int b = sector[s].b;
        int i  = sector[s].i;
        int j  = sector[s].j;
        int x  = sector[s].x;
        int y  = sector[s].y;
        int ml = sector[s--].ml;

        int idx = index(a,b,i,j,x,y);
        int pa = protein[a];
        int pb = protein[b];
        int li = sigma(a,i);
        int rj = sigma(b,j);



        int xi = nucleotides[pa][x][i];
        int yj = nucleotides[pb][y][j];
        if (a == b && i == j) {
            break;
        }

        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;

        codon_selection[a] = x;
        codon_selection[b] = y;



        if (ml==2) {
            bp_bond[++t].i = li;
            bp_bond[t].j   = rj;
            goto repeat;
        }

        if (ml == 0) {
            double energy;
            bt = Access_OB(idx);
            oij = Access_O(a,b,i,j,x,y);
            switch (bt) {
                int kj; //c1, i1_
                case -1:
                    oi = Access_E1(a,b,i,j,x,y) + lambda*AU[xi][yj];
                    bp_bond[++t].i = sigma(a,i);
                    bp_bond[t].j   = sigma(b,j);
                    goto repeat;
                    break;
                case -2:
                    oi = Access_O(a,b,i,j-1,x,y);
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i, sector[s].j = j-1, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -3:
                    kj = O_bt[idx][0];
                    energy = Access_O(a,b-1,i,2,x,kj) + (lambda-1)*codon_cai[pb][y];

                    sector[++s].a = a;
                    sector[s].b = b - 1, sector[s].i = i, sector[s].j = 2, sector[s].x = x, sector[s].y = kj, sector[s].ml = ml;

                    goto OUTLOOP;

                    break;
                case -4:
                    int al, il, xl, br, jr, yr, hi;
                    al = O_bt[idx][0], il = O_bt[idx][1], xl = O_bt[idx][2], br = O_bt[idx][3], jr = O_bt[idx][4],yr = O_bt[idx][5], hi = O_bt[idx][6];
                    if (jr >= 1) {
                        energy = Access_O(a,al,i,il,x,xl) + Access_E1(br,b,jr,j,yr,y) - (lambda-1)*codon_cai[protein[al]][xl] + lambda * AU[hi][yj];
                    } else {
                        energy = Access_O(a,al,i,il,x,xl) + Access_E1(br,b,jr,j,yr,y) + lambda * AU[hi][yj];
                    }

                    sector[++s].a = a;
                    sector[s].i = i, sector[s].b = al, sector[s].j = il, sector[s].x = x, sector[s].y = xl, sector[s].ml = 0;


                    a = br;
                    i = jr;
                    x = yr;
                    bp_bond[++t].i = sigma(a,i);
                    bp_bond[t].j   = sigma(b,j);
                    goto repeat;


                    break;
                default:
                    cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << oij << " " << bt << endl;
                    exit(1);
            }
        }
        else {
            bt = Access_MB(idx);
            oij = Access_M1(a,b,i,j,x,y);
            switch (bt) {
                int hi, kj;
                int c1, i1_;
                double energy;
                case -1:
                    bp_bond[++t].i = sigma(a, i);
                    bp_bond[t].j = sigma(b, j);
                    goto repeat;
                    break;
                case -2:
                    oi = Access_M1(a, b, i + 1, j, x, y) + lambda * ML_BASE;
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i+1, sector[s].j = j, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -3:
                    hi = M_bt[idx][0];
                    energy = Access_M1(a + 1, b, 0, j, hi, y) + (lambda-1) * codon_cai[pa][x] + lambda * ML_BASE;
                    sector[++s].a = a + 1;
                    sector[s].b = b, sector[s].i = 0, sector[s].j = j, sector[s].x = hi, sector[s].y = y, sector[s].ml = ml;
                    goto OUTLOOP;
                    break;
                case -4:
                    oi = Access_M1(a,b,i,j-1,x,y) + lambda * ML_BASE;
                    sector[++s].a = a;
                    sector[s].b = b, sector[s].i = i, sector[s].j = j-1, sector[s].x = x, sector[s].y = y, sector[s].ml = ml;

                    goto OUTLOOP;
                    break;
                case -5:
                    kj = M_bt[idx][0];
                    energy = Access_M1(a,b-1,i,2,x,kj) + (lambda-1)*codon_cai[pb][y] + lambda * ML_BASE;
                    sector[++s].a = a;
                    sector[s].b = b - 1, sector[s].i = i, sector[s].j = 2, sector[s].x = x, sector[s].y = kj, sector[s].ml = ml;

                    goto OUTLOOP;
//                    }
                    break;
                case -6:

                    c = M_bt[idx][0], i1 =M_bt[idx][1], hi = M_bt[idx][2], c1 = M_bt[idx][3], i1_ = M_bt[idx][4], kj = M_bt[idx][5];
                    if (i1_ >= 1) {
                        energy = Access_M1(a, c, i, i1, x, hi) + Access_M1(c1, b, i1_, j, kj, y) - (lambda-1) * codon_cai[protein[c]][kj];
                    } else {
                        energy = Access_M1(a, c, i, i1, x, hi) + Access_M1(c1, b, i1_, j, kj, y);
                    }

                    sector[++s].a = a;
                    sector[s].i = i, sector[s].b = c, sector[s].j = i1, sector[s].x = x, sector[s].y = hi, sector[s].ml = ml;

                    sector[++s].a = c1;
                    sector[s].i = i1_, sector[s].b = b, sector[s].j = j, sector[s].x = kj, sector[s].y = y, sector[s].ml = ml;

                    goto OUTLOOP;

                    break;
                default:
                    cout << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << oij << endl;
                    exit(3);
            }
        }

        repeat:
        eij = Access_E1(a,b,i,j,x,y);

        rj = sigma(b,j);
        li = sigma(a,i);
        pa = protein[a];
        pb = protein[b];
        xi = nucleotides[pa][x][i];
        yj = nucleotides[pb][y][j];

        if (a > b || (a == b && i >= j)) {
            break;
        }


        nucle_seq[li] = xi;
        nucle_seq[rj] = yj;
        codon_selection[a] = x;
        codon_selection[b] = y;

        idx = index(a,b,i,j,x,y);
        bt = Access_EB(idx);
        string h,str;
        vector<int> values, nuc;
        vector<int> ss, codon;
        switch (bt) {
            int hi, kj, lc, ld, ll; //pc, pd, n_codon_c, n_codon_d, min_ld
            int c1, i1_; //,d1, j1_
            int xi_, _yj; // cx, cy
            double energy; //interior_energy, internal_energy
            int l,a1,b1,x1,y1,a2,x2,b2,y2;
            char xi1, xi2,xi3,yj3,yj2,yj1;


            case -1:
                double hairpin, cai;

                values = E_bt[idx];

                if (find(values.begin(), values.end(), -1) == values.end()) {
                    l = values[0];
                    for (int pos = 1; pos < l; ++pos) {
                        switch (pos) {
                            case 1:
                                xi1 = to_char[values[pos]];
                                break;
                            case 2:
                                xi2 = to_char[values[pos]];
                                break;
                            case 3:
                                xi3 = to_char[values[pos]];
                                break;
                            case 4:
                                yj3 = to_char[values[pos]];
                                break;
                            case 5:
                                yj2 = to_char[values[pos]];
                                break;
                            case 6:
                                yj1 = to_char[values[pos]];
                                break;
                            default:
                                break;
                        }
                    }
                    for (int pos = l; pos < (int)values.size(); ++pos) {
                        switch (pos-l) {
                            case 0:
                                a1 = values[pos];
                                break;
                            case 1:
                                b1 = values[pos];
                                break;
                            case 2:
                                x1 = values[pos];
                                break;
                            case 3:
                                y1 = values[pos];
                                break;
                            case 4:
                                a2 = values[pos];
                                break;
                            case 5:
                                x2 = values[pos];
                                break;
                            case 6:
                                b2 = values[pos];
                                break;
                            case 7:
                                y2 = values[pos];
                                break;
                            default:
                                break;

                        }
                    }
                    switch (l) {
                        case 4:
                            nuc = {to_int(xi1), to_int(xi2),to_int(xi3)};
                            str = {to_char[xi],xi1,xi2,xi3,to_char[yj]};
                            break;
                        case 5:
                            nuc = {to_int(xi1),to_int(xi2),to_int(xi3),to_int(yj3)};
                            str = {to_char[xi],xi1,xi2,xi3,yj3,to_char[yj]};
                            break;
                        case 7:
                            nuc = {to_int(xi1),to_int(xi2),to_int(xi3),to_int(yj3),to_int(yj2),to_int(yj1)};
                            str = {to_char[xi],xi1,xi2,xi3,yj3,yj2,yj1,to_char[yj]};
                            break;
                        default:
                            cout << l << endl;
                            exit(4);
                            break;
                    }
                    hairpin = lambda * hairpinE[str];
                    switch (values.size() - l) {
                        case 4:
                            cai = (lambda - 1) * add_hairpin_CAI_8(a1,b1,x1,y1);
                            break;
                        case 6:
                            cai = (lambda - 1) * add_hairpin_CAI_8(a1,b1,x1,y1,a2,x2);
                            codon_selection[a2] = x2;
                            break;
                        case 8:
                            cai = (lambda - 1) * add_hairpin_CAI_8(a1,b1,x1,y1,a2,x2,b2,y2);
                            codon_selection[a2] = x2;
                            codon_selection[b2] = y2;
                            break;
                        default:
                            cout << values.size() - l - 1 << endl;
                            exit(6);
                            break;
                    }

                    assign(nucle_seq,nuc,li+1);
                    goto OUTLOOP;
                } else {

                    l = values[0], xi_ = values[2], _yj = values[3], a1 = values[5], x1 = values[6] , b1 = values[7], y1 = values[8];
                    energy = lambda*hairpin_loop(xi,yj,xi_,_yj,l-1) + (lambda-1)*add_hairpin_CAI_2(a,b,x,y,a1,x1,b1,y1);
                    nucle_seq[li+1] = xi_;
                    nucle_seq[rj-1] = _yj;
                    goto OUTLOOP;
                }
                break;
            case -3:

                sector[s+1].ml = sector[s+2].ml = 1;
                double en;
                en = eij - lambda*(AU[xi][yj] + ML_intern + ML_closing);
                int t_idx;
                int sa,sb,si,sj,sx,sy,sh,sk;


                sa = E_bt[idx][0], sb = E_bt[idx][1], si = E_bt[idx][2], sj = E_bt[idx][3], sx = E_bt[idx][4], sy = E_bt[idx][5];

                if (i == 2 && j == 0) {
                    en -= ((lambda-1)*(codon_cai[pa][x] + codon_cai[pb][y]));
                } else if (i == 2) {
                    en -= (lambda-1)*codon_cai[pa][x];
                } else if (j == 0) {
                    en -= (lambda-1)*codon_cai[pb][y];
                }
                t_idx = index(sa,sb,si,sj,sx,sy);


                c = TM_bt[t_idx][0], i1 =TM_bt[t_idx][1], sh = TM_bt[t_idx][2], c1 = TM_bt[t_idx][3], i1_ = TM_bt[t_idx][4], sk = TM_bt[t_idx][5];


                if (i1_ >= 1) {
                    energy = Access_M1(sa,c,si,i1,sx,sh) + Access_M1(c1,sb,i1_,sj,sk,sy) - (lambda-1)*codon_cai[protein[c1]][sk];
                } else {
                    energy = Access_M1(sa,c,si,i1,sx,sh) + Access_M1(c1,sb,i1_,sj,sk,sy);
                }


                sector[++s].a = sa;
                sector[s].i = si, sector[s].b = c, sector[s].j = i1, sector[s].x = sx, sector[s].y = sh;


                sector[++s].a = c1;
                sector[s].i = i1_, sector[s].b = sb, sector[s].j = sj, sector[s].x = sk, sector[s].y = sy;

                break;

            case -4:

                c = E_bt[idx][0], d = E_bt[idx][1], i1 = E_bt[idx][2], j1 = E_bt[idx][3], hx = E_bt[idx][4], ky = E_bt[idx][5], hi = E_bt[idx][6], kj = E_bt[idx][7];
                energy = Access_E1(c, d, i1, j1, hx, ky) + lambda * stacking(xi, yj, hi, kj) + (lambda-1)*(add_CAI(a,c,x) + add_CAI(b,d,y));
                lc = sigma(c,i1), ld = sigma(d,j1);

                bp_bond[++t].i = lc;
                bp_bond[t].j = ld;


                nucle_seq[lc] = hi;
                nucle_seq[ld] = kj;
                a = c, b = d, i = i1, j = j1, x = hx, y = ky;

                goto repeat;
                break;
            case -5:

                c = E_bt[idx][0], d = E_bt[idx][1], i1 = E_bt[idx][2], j1 = E_bt[idx][3], hx = E_bt[idx][4], ky = E_bt[idx][5], hi = E_bt[idx][6], kj = E_bt[idx][7], ll = E_bt[idx][8];
                energy = Access_E1(c, d, i1, j1, hx, ky) + lambda * bulge_loop(xi, yj, hi, kj, ll)  + (lambda-1)*(add_CAI(a,c,x) + add_CAI(b,d,y));
                lc = sigma(c,i1), ld = sigma(d,j1);

                bp_bond[++t].i = lc;
                bp_bond[t].j = ld;

                nucle_seq[lc] = hi;
                nucle_seq[ld] = kj;

                a = c, b = d, i = i1, j = j1, x = hx, y = ky;

                goto repeat;
//                }
                break;
            case -6:
                int _hi, kj_, lr;
                int na,xa,cp,xc,bp,xb,nd,xd;

                c = E_bt[idx][0], d = E_bt[idx][1], i1 = E_bt[idx][2], j1 = E_bt[idx][3], hx = E_bt[idx][4], ky = E_bt[idx][5], hi = E_bt[idx][6], kj = E_bt[idx][7];
                xi_ = E_bt[idx][8], _yj = E_bt[idx][9], _hi = E_bt[idx][10], kj_ = E_bt[idx][11],  ll = E_bt[idx][12], lr = E_bt[idx][13];
                na = E_bt[idx][14], xa = E_bt[idx][15], cp = E_bt[idx][16], xc = E_bt[idx][17], bp = E_bt[idx][18], xb = E_bt[idx][19], nd = E_bt[idx][20], xd = E_bt[idx][21];

                lc = sigma(c,i1), ld = sigma(d,j1);

                energy = Access_E1(c, d, i1, j1, hx, ky) + lambda* interior_loop(xi, yj, hi, kj, xi_, _yj, _hi, kj_,ll, lr)
                                                           + (lambda-1)*(add_interior_CAI_2(a,c,x,na,xa,cp,xc) + add_interior_CAI_2(b,d,y,bp,xb,nd,xd));


                bp_bond[++t].i = lc;
                bp_bond[t].j = ld;

                nucle_seq[lc] = hi;
                nucle_seq[ld] = kj;

                nucle_seq[li + 1] = xi_;
                nucle_seq[rj - 1] = _yj;
                nucle_seq[lc - 1] = _hi;
                nucle_seq[ld + 1] = kj_;

                a = c, b = d, i = i1, j = j1, x = hx, y = ky;

                goto repeat;

                break;
            default:
                cout << "backtracking failed in repeat: " << a << " " << b << " " << i << " " << j << " " << x << " " << y << " " << energy << " " << eij << " " << bt << " " << endl;
                exit(3);
        }
    }
    bp_bond[0].i = t;


}

void Zuker::assign_codon(vector<int> &s, int sp) {
    int m = (int)s.size()/3;
    for (int i = 0; i < m; ++i) {
        vector<int> codon(3);
        int p = protein[sp+i];
        for (int j = 0; j <= 2; j++) {

            codon[j] = s[3*i+j];

        }
        int x = getxPos(p, codon);
        codon_selection[p] = x;
    }
}

double Zuker::calculate_CAI_O(ostream & fout, double lambda) {
    auto start = chrono::high_resolution_clock::now();
    fout << "Zuker CAI" << endl;

    double ret = inf, energy = inf;
    double mfe, cai, t_mfe, t_cai;
    int t = 0;

    if (lambda != 1.0) {
        for (int a = 0; a < n; ++a) {
            int pa = protein[a];
            int n_codon_a = n_codon[pa];
            for (int i = 0; i < 3; ++i) {
                for (int j = i; j < 3; ++j) {
                    for (int x = 0; x < n_codon_a; ++x) {
                        Access_O(a, a, i, j, x, x) = (lambda - 1) * codon_cai[pa][x];
                        Z_CAI[index(a, a, i, j, x, x)] = (lambda - 1) * codon_cai[pa][x];
                        E_CAI[index(a, a, i, j, x, x)] = (lambda - 1) * codon_cai[pa][x];
                        M_CAI[index(a, a, i, j, x, x)] = (lambda - 1) * codon_cai[pa][x];
                    }
                }
            }
        }
    }

    ofstream debug("dd.txt");

    calculate_CAI_E(lambda);
    cout << "E done" << endl;


    int nuc_len = 3*n;
    int a = 0, i = 0;
    int la = sigma(a,i);
    int pa = protein[a];

    const int n_codon_a = n_codon[pa];
    for (int len = 1; len < nuc_len; ++len) {
        int j = (i+len) % 3;
        int b = (3*a+i+len) / 3;
        if (b == n) continue;

        int pb = protein[b];
        const int n_codon_b = n_codon[pb];
        int ppb = protein[b-1];
        for (int x = 0; x < n_codon_a; x++) {

            for (int y = 0; y < n_codon_b; y++) {
                if (a == b && x != y) continue;
                ret = inf, energy = inf, t = 0;
                // i, j pair
                int lb = sigma(b,j);
                int idx = index(a, b, i, j, x, y);
                int yj = nucleotides[pb][y][j];
                int xi = nucleotides[pa][x][i];
                int type = BP_pair[xi+1][yj+1];

                // paired ends
                if (type > 0) {
                    t_mfe = Access_E2(idx) + lambda * AU[xi][yj];
                    t_cai = E_CAI[idx];
                    ret = min(ret,t_mfe + t_cai);
                    if (energy > ret) {
                        energy = ret;
                        t = -1;
                        mfe = t_mfe;
                        cai = t_cai;
                    }
                }

                // [a][b][i][j − 1][x][y]
                if (j >= 1) {
                    t_mfe = Access_Z2(a,b,i,j-1,x,y);
                    t_cai = Z_CAI[index(a,b,i,j-1,x,y)];
                    ret = min(ret, t_mfe + t_cai); //idx - 36
                    if (energy > ret) {
                        energy = ret;
                        t = -2;
                        mfe = t_mfe;
                        cai = t_cai;
                    }
                }

                int my;
                if (j == 0) {
                    // miny′ ∈S(vb−1 )  [a][b − 1][i][3][x][y′]
                    if (a < b - 1) {
                        const int n_codon_bp = n_codon[ppb];

                        double temp_e = -inf;
                        for (int q = 0; q < n_codon_bp; q++) {
                            t_mfe = Access_Z2(a, b - 1, i, 2, x, q);
                            t_cai = Z_CAI[index(a, b - 1, i, 2, x, q)] + (lambda - 1) * codon_cai[pb][y];
                            temp_e = t_mfe + t_cai; //idx_t + q
                            if ((ret > temp_e)) {
                                ret = temp_e;
                                my = q;
                                mfe = t_mfe;
                                cai = t_cai;
                            }
                        }
                    }
                    //  [a][b − 1][i][3][x][x]
                    if (a == b - 1) {
                        t_mfe = Access_Z2(a, b - 1, i, 2, x, x);
                        t_cai = Z_CAI[index(a, b - 1, i, 2, x, x)] + (lambda - 1) * codon_cai[pb][y];
                        double temp_e = t_mfe + t_cai;
                        if (ret > temp_e) {
                            ret = temp_e;
                            my = x;
                            mfe = t_mfe;
                            cai = t_cai;
                        }

                    }
                    if (energy > ret) {
                        energy = ret;
                        t = -3;
                    }

                }


                // bifurication
                int al, il, xl, br, jr, yr, hi_;
                for (int lc = la + 1; lc <= lb - 4; lc++) {
                    int c = lc / 3;
                    int i1 = lc % 3;
                    int pc = protein[c];
                    int n_codon_c = n_codon[pc];

                    double temp_e;
                    for (int hx = 0; hx < n_codon_c; ++hx) {
                        if (c == a && x != hx) continue;
                        int hi = nucleotides[pc][hx][i1];
                        if (i1 >= 1) {
                            t_mfe = Access_Z2(a,c,i,i1-1,x,hx) + Access_E2(c,b,i1,j,hx,y) + lambda*AU[hi][yj];
                            t_cai = Z_CAI[index(a,c,i,i1-1,x,hx)] + E_CAI[index(c,b,i1,j,hx,y)] - (lambda-1) * codon_cai[pc][hx];
                            temp_e = t_mfe + t_cai; // , idx_1 + hx // , idx_2 + 6*hx
                            if (ret > temp_e) {
                                ret = temp_e;
                                al = c, il = i1-1, xl = hx, br = c, jr = i1, yr = hx, hi_=hi;
                                mfe = t_mfe;
                                cai = t_cai;
                            }

                        } else {
                            int n_codon_cp = n_codon[protein[c-1]];
                            for (int ky = 0; ky < n_codon_cp; ++ky) {
                                if (a == c-1 && x != ky) continue;
                                t_mfe = Access_Z2(a,c-1,i,2,x,ky) + Access_E2(c,b,i1,j,hx,y) + lambda*AU[hi][yj];
                                t_cai = Z_CAI[index(a,c-1,i,2,x,ky)] + E_CAI[index(c,b,i1,j,hx,y)];
                                temp_e = t_mfe + t_cai; //, idx_3 + ky //, idx_4 + 6*hx
                                if (ret > temp_e) {
                                    ret = temp_e;
                                    al = c-1, il = 2, xl = ky, br = c, jr = i1, yr = hx, hi_=hi;
                                    mfe = t_mfe;
                                    cai = t_cai;
                                }

                            }
                        }

                    }

                }

                if (energy > ret) {
                    energy = ret;
                    t = -4;
                }
                if (t == -4) {
                    O_bt[idx] = {al,il,xl,br,jr,yr,hi_};
                } else if (t == -3) {
                    O_bt[idx] = {my};
                }

                if (b == n-1 && ret < e) {
                    e = ret;
                    minX = x, minY = y;
                }

                Access_O(idx) = ret;
                Access_OB(idx) = t;


                Access_Z2(idx) = mfe;
                Z_CAI[idx] = cai;

            }
        }

    }


    auto end = chrono::high_resolution_clock::now();
    long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
    fout << "Energy: " << Access_O(0,n-1,0,2,minX, minY)/100 << endl;
    fout << "Time taken by DP is : " << time_take;
    fout << "sec" << endl;

    double res = Access_O(0,n-1,0,2,minX, minY);

    mfe = Access_Z2(0,n-1,0,2,minX, minY);
    cai = Z_CAI[index(0,n-1,0,2,minX, minY)];
    fout << "lambda: " << lambda << ",O: " << res << ",mfe: " << mfe/lambda << ",cai: " << cai/(lambda-1) << ",combined: " << mfe +  cai << endl;


    return res;
}

void Zuker::calculate_CAI_E(double lambda) {
    double min_energy = inf, energy = inf;
    int t;
    static vector<int> temp;

    int nuc_len = 3*n;
    for (int len = 4; len < nuc_len; ++len) {
        std::cout << "\rE/M -- Completed: " << (len*1.0/nuc_len) * 100 << "%" << std::flush;
        int max_a = n - (int)floor(len/3);
        for (int a = 0, b; a < max_a; ++a) {
            for (int i = 0; i < 3; ++i) {

                int j = (i+len) % 3;
                b = (3*a+i+len) / 3;

                if (b == n) continue;
                int pa = protein[a];
                int pb = protein[b];
                int pna = protein[a + 1];
                int ppb = protein[b - 1];
                const int n_codon_a = n_codon[pa];
                const int n_codon_b = n_codon[pb];
                const int n_codon_an = n_codon[pna];
                const int n_codon_bp = n_codon[ppb];
                int la = sigma(a,i);
                int lb = sigma(b,j);
                int l = lb - la;
                for (int x = 0; x < n_codon_a; x++) {
                    for (int y = 0; y < n_codon_b; y++) {
                        int xi = nucleotides[pa][x][i];
                        int yj = nucleotides[pb][y][j];
                        int type = BP_pair[xi+1][yj+1];


                        // the other end does not pair
                        if (type == 0) {
                            calculate_CAI_M(a,b,i,j,x,y,lambda);
                            continue;
                        }
                        min_energy = inf, energy = inf, t = 0;

                        // hairpin
                        min_energy = min(min_energy, hairpin_CAI(lambda,l,a,b, pa, pb, pna, ppb, n_codon_an, n_codon_bp, xi, yj, i, j, x, y));
                        if (greaterThan(energy, min_energy)) {
                            energy = min_energy;
                            temp = E_bt[index(a,b,i,j,x,y)];
                            t = -1;
                        }


                        // interior loops: stacking, bulge, or internal
                        min_energy = min(min_energy, internal_CAI(lambda,a,b,i,j, x, y, la, lb, xi, yj));


                        if (greaterThan(energy, min_energy)) {
                            energy = min_energy;
                            t = Access_EB(a,b,i,j,x,y);
                            temp = E_bt[index(a,b,i,j,x,y)];
                        }


                        // multiloop
                        min_energy = min(min_energy, multi_loop_CAI(lambda,a, b, i, j, x, y, pa, pb, n_codon_an, n_codon_bp));
                        if (greaterThan(energy, min_energy)) {
                            energy = min_energy;
                            t = -3;
                            temp = E_bt[index(a,b,i,j,x,y)];
                        }

                        Access_E1(a,b,i,j,x,y) = min_energy;

                        Access_EB(a,b,i,j,x,y) = t;
                        E_bt[index(a,b,i,j,x,y)] = temp;


                        calculate_CAI_M(a,b,i,j,x,y,lambda);
                    }
                }
            }
        }

    }
}

double Zuker::add_hairpin_CAI_2(int a, int b, int x, int y, int a1, int x1, int b1, int y1) const {
    double cai = 0;
    cai = codon_cai[protein[a]][x];
    if (b != a) cai += codon_cai[protein[b]][y];

    if (a1 == -1 && x1 == -1 && b1 == -1 && y1 == -1) {
        return cai;
    } else if (b1 == -1 && y1 == -1) {
        if (a1 != a && a1 != b) {
            cai += codon_cai[protein[a1]][x1];
        }
        return cai;
    } else {
        if (a1 != a && a1 != b) {
            cai += codon_cai[protein[a1]][x1];
        }
        if (b1 != a && b1 != b && b1 != a1) {
            cai += codon_cai[protein[b1]][y1];
        }
        return cai;
    }

}

double Zuker::add_hairpin_CAI_8(int a, int b, int x, int y, int a1, int x1, int b1, int y1) const {
    double cai = 0;
    cai = codon_cai[protein[a]][x];
    if (b != a) cai += codon_cai[protein[b]][y];
    if (a1 == -1 && x1 == -1 && b1 == -1 && y1 == -1) {
        return cai;
    }
    else if (b1 == -1 && y1 == -1) {
        if (a1 != a && a1 != b)
            cai += codon_cai[protein[a1]][x1];
        return cai;
    }
    else {

        if (a1 != a && a1 != b)
            cai += codon_cai[protein[a1]][x1];

        if (b1 != a && b1 != b && b1 != a1)
            cai += codon_cai[protein[b1]][y1];

        return cai;
    }


}

double Zuker::add_hairpin_CAI_3(vector<int> &s, int sp) const {
    int m = (int)s.size()/3;
//    cout << m << endl;
    double CAI_ans = 0;
    for (int i = 0; i < m; ++i) {
        vector<int> codon(3);
        int p = protein[sp+i];
        for (int j = 0; j <= 2; j++) {
            codon[j] = s[3*i+j];
        }
        int x = getxPos(p, codon);
        CAI_ans += codon_cai[p][x];

    }
    return CAI_ans;
}

void Zuker::calculate_CAI_M(int a, int b, int i, int j, int x, int y, double lambda) { //vector<double> & TM

    double mfe = inf, cai = 0;
    double min_energy = inf, energy = inf;
    int t = 0;
    int la = sigma(a,i);
    int lb = sigma(b,j);
    int pna = protein[a+1];
    int ppb = protein[b-1];
    const int n_codon_an = n_codon[pna];
    const int n_codon_bp = n_codon[ppb];
    int pa = protein[a];
    int pb = protein[b];
    int ni = nucleotides[pa][x][i];
    int nj = nucleotides[pb][y][j];
    int type = BP_pair[ni+1][nj+1];

    int idx = index(a,b,i,j,x,y);

    double temp_e;
    static vector<int> temp, temp2;

    if (type > 0) {
        min_energy = min(min_energy, Access_E1(idx) + lambda*ML_intern + lambda*AU[ni][nj]);
        if (greaterThan(energy, min_energy)) {
            energy = min_energy;
            t = -1;
            mfe = Access_E2(idx) + lambda*ML_intern + lambda*AU[ni][nj];
            cai = E_CAI[idx];
        }
    }

    if (i <= 1) {
        min_energy = min(min_energy, Access_M1(a,b,i+1,j,x,y) + lambda*ML_BASE); //, idx - 108
        if (greaterThan(energy, min_energy)) {
            energy = min_energy;
            t = -2;
            mfe = Access_M2(a,b,i+1,j,x,y) + lambda*ML_BASE;
            cai = M_CAI[index(a,b,i+1,j,x,y)];
        }
    }

    if (i == 2) {
        if (a + 1 < b) {
            for (int x1 = 0; x1 < n_codon_an; ++x1) {
                temp_e =  Access_M1(a+1,b,0,j,x1,y) + lambda*ML_BASE + (lambda-1)*codon_cai[pa][x]; //idx_t + 6*x1
                if (min_energy > temp_e) {
                    min_energy = temp_e;
                    temp2 = {x1};
                    mfe = Access_M2(a+1,b,0,j,x1,y) + lambda*ML_BASE;
                    cai = M_CAI[index(a+1,b,0,j,x1,y)] + (lambda-1)*codon_cai[pa][x];
                }
            }
        }

        if (greaterThan(energy, min_energy)) {
            energy = min_energy;
            temp = temp2;

            t = -3;
        }
    }


    if (j >= 1) {
        min_energy = min(min_energy, Access_M1(a,b,i,j-1,x,y) + lambda*ML_BASE); //, idx - 36
        if (greaterThan(energy, min_energy)) {
            energy = min_energy;
            t = -4;
            mfe = Access_M2(a,b,i,j-1,x,y) + lambda*ML_BASE;
            cai = M_CAI[index(a,b,i,j-1,x,y)];
        }
    }


    if (j == 0) {
        if (a < b - 1) {
            for (int y1 = 0; y1 < n_codon_bp; ++y1) {
                temp_e = Access_M1(a,b-1,i,2,x,y1) + lambda*ML_BASE + (lambda-1)*codon_cai[pb][y]; //idx_t + y1
                if (min_energy > temp_e) {
                    min_energy = temp_e;
                    temp2 = {y1};
                    mfe = Access_M2(a,b-1,i,2,x,y1) + lambda*ML_BASE;
                    cai = M_CAI[index(a,b-1,i,2,x,y1)] + (lambda-1)*codon_cai[pb][y];
                }
            }
        }
        if (greaterThan(energy, min_energy)) {
            energy = min_energy;
            temp = temp2;
            t = -5;
        }

    }


    // bifurication

    double bi_energy = inf;
    double bi_mfe, bi_cai;
    vector<int> bi_temp;
    for (int lc = la + 5; lc <= lb-4; lc++) {
        int c = lc / 3;
        int i1 = lc % 3;
        int pc = protein[c];
        int n_codon_c = n_codon[pc];
        for (int hx = 0; hx < n_codon_c; ++hx) {
            if (i1 >= 1) {
                temp_e = Access_M1(a,c,i,i1-1,x,hx) + Access_M1(c,b,i1,j,hx,y) - (lambda-1)*codon_cai[pc][hx]; //, idx_1 + hx // , idx_2 + 6*hx
                if (bi_energy > temp_e) {
                    bi_energy = temp_e;
                    bi_temp = {c,i1-1,hx,c,i1,hx};
                    bi_mfe = Access_M2(a,c,i,i1-1,x,hx) + Access_M2(c,b,i1,j,hx,y);
                    bi_cai = M_CAI[index(a,c,i,i1-1,x,hx)] + M_CAI[index(c,b,i1,j,hx,y)] - (lambda-1)*codon_cai[pc][hx];
                }

            }
            else {

                {
                    int n_codon_cp = n_codon[protein[c-1]];
                    for (int ky = 0; ky < n_codon_cp; ++ky) {

                        temp_e = Access_M1(a,c-1,i,2,x,ky) + Access_M1(c,b,i1,j,hx,y); //, idx_3 + ky //, idx_4 + 6*hx
                        if (bi_energy > temp_e) {
                            bi_energy = temp_e;
                            bi_temp = {c-1,2,ky,c,i1,hx};
                            bi_mfe = Access_M2(a,c-1,i,2,x,ky) + Access_M2(c,b,i1,j,hx,y);
                            bi_cai = M_CAI[index(a,c-1,i,2,x,ky)] + M_CAI[index(c,b,i1,j,hx,y)];
                        }
                    }
                }

            }

        }

    }

    min_energy = min(min_energy, bi_energy);

    if (greaterThan(energy, min_energy)) {
        energy = min_energy;
        t = -6;
        temp = bi_temp;
        mfe = bi_mfe, cai = bi_cai;
    }

    if (TM[idx] > bi_energy) {
        TM[idx] = bi_energy;
        TM_bt[idx] = bi_temp;
        TM2[idx] = bi_mfe;
        TM_CAI[idx] = bi_cai;
    }

    Access_M2(idx) = mfe;
    M_CAI[idx] = cai;


    Access_M1(idx) = min_energy;
    Access_MB(idx) = t;
    M_bt[idx] = temp;

}

void Zuker::lambda_swipe(double incr, ostream &fout, string & outfile) {
    int size = ceil(1/incr + 1);
    vector<double> O_buffer(size, 0);
    vector<double> lambda_buffer(size, 0);
    vector<double> F_buffer(size, 0);
    vector<double> CAI_buffer(size, 0);
    vector<double> stand_CAI(size, 0);
    vector<pair<string, vector<double>>> dataset(5);
    int index = 0;
    double lambda = 0;
    while (lambda < 1.001) {
        fout << "lambda: " << lambda << endl;
        reinit();
        string rna(3*n, '.'), bp(3*n, '.');
        double Ov = calculate_CAI_O(fout,lambda);
        O_buffer[index] = Ov;
        cout << "traceback" << endl;
        traceback_B2(lambda);
        get_rna_cai(rna);
        get_bp(bp);
        fout << "rna: " << rna << endl;
        fout << "bp: " << bp << endl;
        double CAI = evaluate_CAI(rna, protein, 1);
        double cai = evaluate_CAI(rna, protein, 0);
        lambda_buffer[index] = lambda;
        double mfe = evaluate_MFE(rna);
        cout << "lambda: " << lambda << ",O: " << Ov << ",cai: " << CAI << ",sCAI: " << cai << ",mfe: " << mfe << ",combined: " << lambda*mfe+(lambda-1)*CAI << endl;
        F_buffer[index] = mfe;
        CAI_buffer[index] = CAI;
        stand_CAI[index] = cai;
        index++;
        lambda += incr;
    }
    dataset[0] = make_pair("lambda", lambda_buffer);
    dataset[1] = make_pair("F", F_buffer);
    dataset[2] = make_pair("CAI", CAI_buffer);
    dataset[3] = make_pair("sCAI", stand_CAI);
    dataset[4] = make_pair('O', O_buffer);
    write_csv(outfile + ".csv", dataset);
    cout << "swipe done" << endl;
}

void Zuker::lambda_swipe_2(double threshold, double threshold2, ostream &fout, string & outfile) {
    double left_lambda;
    double right_lambda;
    vector<double> O_buffer,lambda_buffer,F_buffer,CAI_buffer,stand_CAI;
    vector<double> MFE_buffer, C_buffer;
    vector<pair<string, vector<double>>> dataset(7);
    queue<pair<double, double>> lambda;
    unordered_map<double, double> mfe_map;
    unordered_map<double, double> cai_map;
    double left_CAI = 0, left_cai = 0, left_mfe = 0, right_CAI = inf, right_cai = inf, right_mfe = inf;

    lambda.push({EPSILON,1-EPSILON});
    int idx;

    while (!lambda.empty()) {
        left_lambda=  lambda.front().first, right_lambda = lambda.front().second;
        lambda.pop();
        string left_rna(3*n, '.'), left_bp(3*n, '.'), right_rna(3*n, '.'), right_bp(3*n, '.');
        if (mfe_map.count(left_lambda) != 0) {
            left_cai = cai_map[left_lambda];
            left_mfe = mfe_map[left_lambda];
        } else {
            reinit();
            fout << "lambda: " << left_lambda << endl;
            double left_O = calculate_CAI_O(fout,left_lambda);
            traceback_B2(left_lambda);
            get_rna_cai(left_rna);
            get_bp(left_bp);
            left_CAI = evaluate_CAI(left_rna, protein, 1);
            left_cai = evaluate_CAI(left_rna, protein, 0);
            left_mfe = evaluate_MFE(left_rna);
            mfe_map[left_lambda] = left_mfe;
            cai_map[left_lambda] = left_cai;
            O_buffer.push_back(left_O);
            lambda_buffer.push_back(left_lambda);
            F_buffer.push_back(left_mfe);
            CAI_buffer.push_back(left_CAI);
            stand_CAI.push_back(left_cai);
            idx = index(0, n-1, 0, 2, minX, minY);
            MFE_buffer.push_back(Z2[idx] / left_lambda);
            C_buffer.push_back(Z_CAI[idx] / (left_lambda-1));
            fout << "rna: " << left_rna << endl;
            fout << "bp: " << left_bp << endl;
            cout << "lambda: " << left_lambda << ",O: " << left_O << ",cai: " << Z_CAI[idx] / (left_lambda-1) << ",mfe: " << Z2[idx] / left_lambda << ",integrated: " << Z2[idx] + Z_CAI[idx] << endl;
            cout << "lambda: " << left_lambda << ",O: " << left_O << ",cai: " << left_CAI << ",sCAI: " << left_cai << ",mfe: " << left_mfe << ",combined: " << left_lambda*left_mfe+(left_lambda-1)*left_CAI << endl;

        }

        if (mfe_map.count(right_lambda) != 0) {
            right_cai = cai_map[right_lambda];
            right_mfe = mfe_map[right_lambda];
        } else {
            reinit();
            fout << "lambda: " << right_lambda << endl;
            double right_O = calculate_CAI_O(fout,right_lambda);
            traceback_B2(right_lambda);
            get_rna_cai(right_rna);
            get_bp(right_bp);
            right_CAI = evaluate_CAI(right_rna, protein, 1);
            right_cai = evaluate_CAI(right_rna, protein, 0);
            right_mfe = evaluate_MFE(right_rna);
            mfe_map[right_lambda] = right_mfe;
            cai_map[right_lambda] = right_cai;
            O_buffer.push_back(right_O);
            lambda_buffer.push_back(right_lambda);
            F_buffer.push_back(right_mfe);
            CAI_buffer.push_back(right_CAI);
            stand_CAI.push_back(right_cai);
            idx = index(0, n-1, 0, 2, minX, minY);
            MFE_buffer.push_back(Z2[idx] / right_lambda);
            C_buffer.push_back(Z_CAI[idx] / (right_lambda-1));
            fout << "rna: " << right_rna << endl;
            fout << "bp: " << right_bp << endl;
            cout << "lambda: " << right_lambda << ",O: " << right_O << ",cai: " << Z_CAI[idx] / (right_lambda-1) << ",mfe: " << Z2[idx] / right_lambda << ",integrated: " << Z2[idx] + Z_CAI[idx] << endl;
            cout << "lambda: " << right_lambda << ",O: " << right_O << ",cai: " << right_CAI << ",sCAI: " << right_cai << ",mfe: " << right_mfe << ",combined: " << right_lambda*right_mfe+(right_lambda-1)*right_CAI << endl;

        }
        if (!compare(left_cai, right_cai) && !compare(left_mfe, right_mfe)) {
            cout << "lamda diff: " << right_lambda-left_lambda << endl;
            if (right_lambda < threshold) {
                if (!compare(left_lambda,right_lambda,threshold2)) {
                    double m = (left_lambda + right_lambda) / 2;
                    lambda.push({left_lambda, m});
                    lambda.push({m, right_lambda});
                }
            } else {
                if (!compare(left_lambda,right_lambda,threshold)) {
                    double m = (left_lambda + right_lambda) / 2;
                    lambda.push({left_lambda, m});
                    lambda.push({m, right_lambda});
                }
            }

        }
    }
    cout << "size : " << lambda.size() << endl;
    fout << "left lambda: " << left_lambda << ",right lambda: " << right_lambda << endl;
    dataset[0] = make_pair("lambda", lambda_buffer);
    dataset[1] = make_pair("rc_mfe", F_buffer);
    dataset[2] = make_pair("rc_CAI", CAI_buffer);
    dataset[3] = make_pair("rc_sCAI", stand_CAI);
    dataset[4] = make_pair('O', O_buffer);
    dataset[5] = make_pair("mfe", MFE_buffer);
    dataset[6] = make_pair("CAI", C_buffer);
    write_csv(outfile + ".csv", dataset);
    cout << "swipe done" << endl;
}


void Zuker::lambda_swipe_3(double threshold, double threshold2, ostream &fout, string & outfile) {
    double left_lambda;
    double right_lambda;
    vector<double> O_buffer,lambda_buffer,F_buffer,CAI_buffer,stand_CAI;
    vector<double> MFE_buffer, C_buffer;
    vector<pair<string, vector<double>>> dataset(7);
    queue<pair<double, double>> lambda;
    unordered_map<double, double> mfe_map;
    unordered_map<double, double> cai_map;
    double left_CAI = 0, left_cai = 0, left_mfe = 0, right_CAI = inf, right_cai = inf, right_mfe = inf;

    lambda.push({EPSILON,1-EPSILON});
    int idx;

    while (!lambda.empty()) {
        left_lambda=  lambda.front().first, right_lambda = lambda.front().second;
        lambda.pop();
        string left_rna(3*n, '.'), left_bp(3*n, '.'), right_rna(3*n, '.'), right_bp(3*n, '.');
        if (mfe_map.count(left_lambda) != 0) {
            left_cai = cai_map[left_lambda];
            left_mfe = mfe_map[left_lambda];
        } else {
            reinit();
            fout << "lambda: " << left_lambda << endl;
            double left_O = calculate_CAI_O(fout,left_lambda);
            traceback_B2(left_lambda);
            get_rna_cai(left_rna);
            get_bp(left_bp);
            left_CAI = evaluate_CAI(left_rna, protein, 1);
            left_cai = evaluate_CAI(left_rna, protein, 0);
            left_mfe = evaluate_MFE(left_rna);
            mfe_map[left_lambda] = left_mfe;
            cai_map[left_lambda] = left_cai;
            O_buffer.push_back(left_O);
            lambda_buffer.push_back(left_lambda);
            F_buffer.push_back(left_mfe);
            CAI_buffer.push_back(left_CAI);
            stand_CAI.push_back(left_cai);
            idx = index(0, n-1, 0, 2, minX, minY);
            MFE_buffer.push_back(Z2[idx] / left_lambda);
            C_buffer.push_back(Z_CAI[idx] / (left_lambda-1));
            fout << "rna: " << left_rna << endl;
            fout << "bp: " << left_bp << endl;
            cout << "lambda: " << left_lambda << ",O: " << left_O << ",cai: " << Z_CAI[idx] / (left_lambda-1) << ",mfe: " << Z2[idx] / left_lambda << ",integrated: " << Z2[idx] + Z_CAI[idx] << endl;
            cout << "lambda: " << left_lambda << ",O: " << left_O << ",cai: " << left_CAI << ",sCAI: " << left_cai << ",mfe: " << left_mfe << ",combined: " << left_lambda*left_mfe+(left_lambda-1)*left_CAI << endl;

        }

        if (mfe_map.count(right_lambda) != 0) {
            right_cai = cai_map[right_lambda];
            right_mfe = mfe_map[right_lambda];
        } else {
            reinit();
            fout << "lambda: " << right_lambda << endl;
            double right_O = calculate_CAI_O(fout,right_lambda);
            traceback_B2(right_lambda);
            get_rna_cai(right_rna);
            get_bp(right_bp);
            right_CAI = evaluate_CAI(right_rna, protein, 1);
            right_cai = evaluate_CAI(right_rna, protein, 0);
            right_mfe = evaluate_MFE(right_rna);
            mfe_map[right_lambda] = right_mfe;
            cai_map[right_lambda] = right_cai;
            O_buffer.push_back(right_O);
            lambda_buffer.push_back(right_lambda);
            F_buffer.push_back(right_mfe);
            CAI_buffer.push_back(right_CAI);
            stand_CAI.push_back(right_cai);
            idx = index(0, n-1, 0, 2, minX, minY);
            MFE_buffer.push_back(Z2[idx] / right_lambda);
            C_buffer.push_back(Z_CAI[idx] / (right_lambda-1));
            fout << "rna: " << right_rna << endl;
            fout << "bp: " << right_bp << endl;
            cout << "lambda: " << right_lambda << ",O: " << right_O << ",cai: " << Z_CAI[idx] / (right_lambda-1) << ",mfe: " << Z2[idx] / right_lambda << ",integrated: " << Z2[idx] + Z_CAI[idx] << endl;
            cout << "lambda: " << right_lambda << ",O: " << right_O << ",cai: " << right_CAI << ",sCAI: " << right_cai << ",mfe: " << right_mfe << ",combined: " << right_lambda*right_mfe+(right_lambda-1)*right_CAI << endl;

        }
        if (!compare(left_cai, right_cai) && !compare(left_mfe, right_mfe)) {
            cout << "lamda diff: " << right_lambda-left_lambda << endl;
            if (right_lambda < threshold) {
                if (!compare(left_lambda,right_lambda,threshold2)) {
                    double m = left_lambda * 0.4 + right_lambda * 0.6;
                    lambda.push({left_lambda, m});
                    lambda.push({m, right_lambda});
                }
            } else {
                if (!compare(left_lambda,right_lambda,threshold)) {
                    double m = left_lambda * 0.4 + right_lambda * 0.6;
                    lambda.push({left_lambda, m});
                    lambda.push({m, right_lambda});
                }
            }

        }
    }
    cout << "size : " << lambda.size() << endl;
    fout << "left lambda: " << left_lambda << ",right lambda: " << right_lambda << endl;
    dataset[0] = make_pair("lambda", lambda_buffer);
    dataset[1] = make_pair("rc_mfe", F_buffer);
    dataset[2] = make_pair("rc_CAI", CAI_buffer);
    dataset[3] = make_pair("rc_sCAI", stand_CAI);
    dataset[4] = make_pair('O', O_buffer);
    dataset[5] = make_pair("mfe", MFE_buffer);
    dataset[6] = make_pair("CAI", C_buffer);
    write_csv(outfile + ".csv", dataset);
    cout << "swipe done" << endl;
}

void Zuker::get_bp(string & bp) {
    for (int a = 1; a <= bp_bond[0].i; ++a) {
        bp[bp_bond[a].i] = '(';
        bp[bp_bond[a].j] = ')';
    }
}

void Zuker::get_rna_X(string & rna) {
    for (int i = 0; i < (int)rna.size(); ++i) {
        if (nucle_seq[i] == -1) {
            rna[i] = 'X';
        }
        else rna[i] = to_char[nucle_seq[i]];
    }

}

void Zuker::getAllRna(vector<string> & rna_array, int index, int i) {
    string & rna = rna_array[index];
    while (i < 3*n) {
        if (nucle_seq[i] == -1) {
            if (i % 3 == 0) {
                int n_codon_x = n_codon[protein[i/3]];

                if (nucle_seq[i+1] == -1 && nucle_seq[i+2] == -1) {
                    rna[i] = to_char[nucleotides[protein[i / 3]][0][0]];
                    rna[i + 1] = to_char[nucleotides[protein[i / 3]][0][1]];
                    rna[i + 2] = to_char[nucleotides[protein[i / 3]][0][2]];
                    for (int j = 1; j < n_codon_x; ++j) {
                        rna_array[++index] = rna;
                        rna_array[index][i] = to_char[nucleotides[protein[i / 3]][j][0]];
                        rna_array[index][i + 1] = to_char[nucleotides[protein[i / 3]][j][1]];
                        rna_array[index][i + 2] = to_char[nucleotides[protein[i / 3]][j][2]];
                        getAllRna(rna_array,index,i+3);
                    }
                    i = i+3;
                }
                else if (nucle_seq[i+1] == -1) {

                    rna[i+2] = to_char[nucle_seq[i+2]];
                    int b = n_codon_x;
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna[i+1] = to_char[nucleotides[protein[i/3]][j][1]];
                            b = j;
                            break;
                        }
                    }
                    for (int j = b + 1; j < n_codon_x; ++j) {
                        if (nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna_array[++index] = rna;
                            rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna_array[index][i+1] = to_char[nucleotides[protein[i/3]][j][1]];
                            getAllRna(rna_array,index,i+3);
                        }
                    }
                    i = i+3;
                }
                else if (nucle_seq[i+2] == -1) {
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    int b = n_codon_x;
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna[i+2] = to_char[nucleotides[protein[i/3]][j][2]];
                            b = j;
                            break;
                        }
                    }
                    for (int j = b + 1; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1]) {
                            rna_array[++index] = rna;
                            rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna_array[index][i+2] = to_char[nucleotides[protein[i/3]][j][2]];
                            getAllRna(rna_array,index,i+3);
                        }
                    }
                    i = i+3;
                }
                else {
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    rna[i+2] = to_char[nucle_seq[i+2]];
                    int b= n_codon_x;
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1] && nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            b = j;
                            break;
                        }
                    }
                    for (int j = b + 1; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1] && nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna_array[++index] = rna;
                            rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][0]];
                            getAllRna(rna_array,index,i+3);
                        }
                    }
                    i = i+3;
                }
            }
            if (i % 3 == 1) {
                int n_codon_x = n_codon[protein[i/3]];
                if (nucle_seq[i+1] == -1) {
                    int b = n_codon_x;
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][1]];
                            rna[i+1] = to_char[nucleotides[protein[i/3]][j][2]];
                            b = j;
                            break;
                        }
                    }
                    for (int j = b + 1; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0]) {
                            rna_array[++index] = rna;
                            rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][1]];
                            rna_array[index][i+1] = to_char[nucleotides[protein[i/3]][j][2]];
                            getAllRna(rna_array,index,i+2);
                        }
                    }
                    i = i+2;
                } else {
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    int b = n_codon_x;
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0] && nucle_seq[i+1] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][1]];
                            b = j;
                            break;
                        }
                    }
                    for (int j = b + 1; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0] && nucle_seq[i+1] == nucleotides[protein[i/3]][j][2]) {
                            rna_array[++index] = rna;
                            rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][1]];
                            getAllRna(rna_array,index,i+2);
                        }
                    }
                    i = i+2;
                }
            }
            if (i % 3 == 2) {
                int n_codon_x = n_codon[protein[i/3]];
                int b = n_codon_x;
                for (int j = 0; j < n_codon_x; ++j) {
                    if (nucle_seq[i-2] == nucleotides[protein[i/3]][j][0] && nucle_seq[i-1] == nucleotides[protein[i/3]][j][1]) {
                        rna[i] = to_char[nucleotides[protein[i/3]][j][2]];
                        b = j;
                        break;
                    }
                }
                for (int j = b + 1; j < n_codon_x; ++j) {
                    if (nucle_seq[i-2] == nucleotides[protein[i/3]][j][0] && nucle_seq[i-1] == nucleotides[protein[i/3]][j][1]) {
                        rna_array[++index] = rna;
                        rna_array[index][i] = to_char[nucleotides[protein[i/3]][j][2]];
                        getAllRna(rna_array,index,i+1);
                    }
                }
                i++;
            }
        }
        else {
            rna[i] = to_char[nucle_seq[i]];
            i++;
        }
    }
}

void Zuker::get_rna(string & rna) {
    int i = 0;
    while (i < 3*n) {
        if (nucle_seq[i] == -1) {
            if (i % 3 == 0) {
                if (nucle_seq[i+1] == -1 && nucle_seq[i+2] == -1) {

                    rna[i] = to_char[nucleotides[protein[i/3]][0][0]];
                    rna[i+1] = to_char[nucleotides[protein[i/3]][0][1]];
                    rna[i+2] = to_char[nucleotides[protein[i/3]][0][2]];
                    i = i+3;
                }
                else if (nucle_seq[i+1] == -1) {
                    int n_codon_x = n_codon[protein[i/3]];
                    rna[i+2] = to_char[nucle_seq[i+2]];
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna[i+1] = to_char[nucleotides[protein[i/3]][j][1]];
                            break;
                        }
                    }
                    i = i+3;
                }
                else if (nucle_seq[i+2] == -1) {
                    int n_codon_x = n_codon[protein[i/3]];
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            rna[i+2] = to_char[nucleotides[protein[i/3]][j][2]];
                            break;
                        }
                    }
                    i = i+3;
                }
                else {
                    int n_codon_x = n_codon[protein[i/3]];
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    rna[i+2] = to_char[nucle_seq[i+2]];
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i+1] == nucleotides[protein[i/3]][j][1] && nucle_seq[i+2] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][0]];
                            break;
                        }
                    }
                    i = i+3;
                }
            }
            if (i % 3 == 1) {
                if (nucle_seq[i+1] == -1) {
                    int n_codon_x = n_codon[protein[i/3]];
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][1]];
                            rna[i+1] = to_char[nucleotides[protein[i/3]][j][2]];
                            break;
                        }
                    }
                    i = i+2;
                } else {
                    int n_codon_x = n_codon[protein[i/3]];
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    for (int j = 0; j < n_codon_x; ++j) {
                        if (nucle_seq[i-1] == nucleotides[protein[i/3]][j][0] && nucle_seq[i+1] == nucleotides[protein[i/3]][j][2]) {
                            rna[i] = to_char[nucleotides[protein[i/3]][j][1]];
                            break;
                        }
                    }
                    i = i+2;
                }
            }
            if (i % 3 == 2) {
                int n_codon_x = n_codon[protein[i/3]];
                for (int j = 0; j < n_codon_x; ++j) {
                    if (nucle_seq[i-2] == nucleotides[protein[i/3]][j][0] && nucle_seq[i-1] == nucleotides[protein[i/3]][j][1]) {
                        rna[i] = to_char[nucleotides[protein[i/3]][j][2]];
                        break;
                    }
                }
                i++;
            }
        }
        else {
            rna[i] = to_char[nucle_seq[i]];
            i++;
        }
    }
}

void Zuker::reinit() {
    nucle_seq.clear(), sector.clear(), bp_bond.clear(), O.clear(), E1.clear(), M1.clear(), TM.clear();
    OB.clear(), EB.clear(), MB.clear();
    O_bt.clear(), E_bt.clear(), M_bt.clear(), TM_bt.clear();
    Z2.clear(), E2.clear(), M2.clear(), TM2.clear();
    Z_CAI.clear(), E_CAI.clear(), M_CAI.clear(), TM_CAI.clear();
    codon_selection.clear();
    sector.resize(3*n);
    bp_bond.resize(3*n);
    nucle_seq.resize(3*n,-1);
    O.resize(last_idx, 0);
    E1.resize(last_idx, inf);
    M1.resize(last_idx, inf);
    TM.resize(last_idx, inf);
    OB.resize(last_idx, 0);
    EB.resize(last_idx, 0);
    MB.resize(last_idx, 0);
    O_bt.resize(last_idx);
    E_bt.resize(last_idx);
    M_bt.resize(last_idx);
    TM_bt.resize(last_idx);
    Z2.resize(last_idx, 0);
    E2.resize(last_idx, inf);
    M2.resize(last_idx, inf);
    TM2.resize(last_idx,inf);

    Z_CAI.resize(last_idx, 0);
    E_CAI.resize(last_idx, 0);
    M_CAI.resize(last_idx, 0);
    TM_CAI.resize(last_idx,0);
    codon_selection.resize(n,-1);
    minX = -1, minY = -1;
    e = inf;
}

void Zuker::maxCAISeq() {
    for (int i = 0; i < n; ++i) {
        int p = protein[i];
        int idx = max_cai_pos[p];
        for (int j = 0; j < 3; ++j) {
            nucle_seq[3 * i + j] = nucleotides[p][idx][j];
        }
    }
}

void Zuker::get_rna_cai(string & rna) {


    int i = 0;
    while (i < 3*n) {
        if (nucle_seq[i] == -1) {
            if (i % 3 == 0) {
                if (nucle_seq[i+1] == -1 && nucle_seq[i+2] == -1) {
                    int p = protein[i/3];
                    int x = codon_selection[i/3];

                    if (x == -1)
                        x = max_cai_pos[p];
                    rna[i] = to_char[nucleotides[p][x][0]];
                    rna[i+1] = to_char[nucleotides[p][x][1]];
                    rna[i+2] = to_char[nucleotides[p][x][2]];
                    i = i+3;
                }
                else if (nucle_seq[i+1] == -1) {

                    int idx = codon_selection[i/3];
                    if (idx == -1)
                        idx = fill_rna(i);
                    rna[i+2] = to_char[nucle_seq[i+2]];
                    rna[i] = to_char[nucleotides[protein[i/3]][idx][0]];
                    rna[i+1] = to_char[nucleotides[protein[i/3]][idx][1]];

                    i = i+3;
                }
                else if (nucle_seq[i+2] == -1) {
                    int idx = codon_selection[i/3];
                    if (idx == -1)
                        idx = fill_rna(i);
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    rna[i] = to_char[nucleotides[protein[i/3]][idx][0]];
                    rna[i+2] = to_char[nucleotides[protein[i/3]][idx][2]];

                    i = i+3;
                }
                else {
                    int idx = codon_selection[i/3];
                    if (idx == -1)
                        idx = fill_rna(i);
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    rna[i+2] = to_char[nucle_seq[i+2]];
                    rna[i] = to_char[nucleotides[protein[i/3]][idx][0]];

                    i = i+3;
                }
            }
            if (i % 3 == 1) {
                if (nucle_seq[i+1] == -1) {

                    int idx = codon_selection[i/3];
                    if (idx == -1)
                        idx = fill_rna(i);
                    rna[i] = to_char[nucleotides[protein[i/3]][idx][1]];
                    rna[i+1] = to_char[nucleotides[protein[i/3]][idx][2]];


                    i = i+2;
                } else {
                    int idx = codon_selection[i/3];
                    if (idx == -1)
                        idx = fill_rna(i);
                    rna[i+1] = to_char[nucle_seq[i+1]];
                    rna[i] = to_char[nucleotides[protein[i/3]][idx][1]];

                    i = i+2;
                }
            }
            if (i % 3 == 2) {

                int idx = codon_selection[i/3];
                if (idx == -1)
                    idx = fill_rna(i);
                rna[i] = to_char[nucleotides[protein[i/3]][idx][2]];

                i++;
            }
        }
        else {
            rna[i] = to_char[nucle_seq[i]];
            i++;
        }
    }
}

int Zuker::fill_rna(int index) {
    int l0, l1, l2;
    if (index % 3 == 0) {
        l0 = nucle_seq[index];
        l1 = nucle_seq[index+1];
        l2 = nucle_seq[index+2];
    } else if (index % 3 == 1) {
        l0 = nucle_seq[index-1];
        l1 = nucle_seq[index];
        l2 = nucle_seq[index+1];
    } else {
        l0 = nucle_seq[index-2];
        l1 = nucle_seq[index-1];
        l2 = nucle_seq[index];
    }
    int pi = protein[index/3];
    int n_codon_pi = n_codon[pi];
    double max_cai = -INF;
    int idx = 0;
    // all three nucleotides in codon are specified
    if (l0 >=0 && l1 >= 0 && l2 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][0] == l0 && nucleotides[pi][i][1] == l1 && nucleotides[pi][i][2] == l2) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // the first two nucleotides in codon are specified
    else if (l0 >=0 && l1 >= 0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][0] == l0 && nucleotides[pi][i][1] == l1) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // the last two nucleotides in codon are specified
    else if (l1 >= 0 && l2 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][1] == l1 && nucleotides[pi][i][2] == l2) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // the first and the last nucleotides in codon are specified
    else if (l0 >= 0 && l2 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][0] == l0 && nucleotides[pi][i][2] == l2) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // only the first nucleotide in codon is specified
    else if (l0 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][0] == l0) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // only the second nucleotide in codon is specified
    else if (l1 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][1] == l1) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    }
        // only the third nucleotide in codon is specified
    else if (l2 >=0) {
        for (int i = 0; i < n_codon_pi; ++i) {
            if (nucleotides[pi][i][2] == l2) {
                if (max_cai < codon_cai[pi][i]) {
                    max_cai = codon_cai[pi][i];
                    idx = i;
                }
            }
        }
        return idx;
    } else {
        return max_cai_pos[pi];
    }
}

void Zuker::assign(vector<int> & target, vector<int> & values, int start) {
    int size = (int)values.size();
    int end = start + size;
    for (int i = start; i < end; ++i) {
        target[i] = values[i-start];
    }
}

Zuker::~Zuker() = default;

Zuker::Zuker(const Zuker & Copy):Z(Copy.Z.size()),n(Copy.n) {
    copy(Copy.Z.begin(), Copy.Z.end(), Z.begin());
}
