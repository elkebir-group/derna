//
// Created by xinyu on 4/26/2022.
//

#include <iostream>
#include <set>
#include <cassert>
#include <fstream>
#include <stack>
#include "ZukerAlgorithm.h"
#include "utils.h"

using namespace std;

//ofstream de("d.txt");

ZukerAlgorithm::ZukerAlgorithm(vector<int> & seq, int n):seq(seq),n(n){

    W.resize(n*n, 0);
    V.resize(n*n, inf);
    WM.resize(n*n, inf);
    TM.resize(n*n, inf);
    sector.resize(n);
    base_pair.resize(n);
    bp.resize(n, '.');

    to_char = {
            'A','C','G','U'
    };
}

ZukerAlgorithm::~ZukerAlgorithm() = default;

int ZukerAlgorithm::calculate_W() {
    int min_w = inf;

    calculate_V();
    int i = 0;
    for (int j = 1; j < n; ++j) {
        min_w = inf;
        if (basepair(seq[i],seq[j])) {
            min_w = min(min_w, V[index(i,j)] + AU[seq[i]][seq[j]]);
        }

        min_w = min(min_w, W[index(i,j-1)]);


        int temp = inf;
        for (int k = 1; k <= j-4; ++k) {
            temp = min(temp, W[index(i,k-1)] + V[index(k,j)] + AU[seq[k]][seq[j]]);
        }

        min_w = min(min_w, temp);
//        de << 3 << " " << min_w << endl;

        W[index(i,j)] = min_w;
//        de << "Z - l: " << i << ",r: " << j << ",L: " << seq[i] << ",R: " << seq[j] << ",ret: " << min_w << endl;
    }


    return W[index(0,n-1)];
}

void ZukerAlgorithm::calculate_V() {
    int min_v = inf;
    for (int l = 4; l < n; ++l) {
        for (int i = 0; i < n - l; ++i) {
            int j = i + l;
            if (!basepair(seq[i],seq[j])) {
                calculate_WM(i,j);
                continue;
            }
            min_v = inf;

            if (l+1 == 5 || l+1 == 6 || l+1 == 8) {
//                int index;
                string s(l + 1, '.');
                assign_seq2str(s, i);
//                cout << s << endl;

                switch (l) {
                    case 4:
                        if (hairpinE.count(s) > 0) {
                            min_v = min(min_v, hairpinE[s]);
                        }
                        break;
                    case 5:
                        if (hairpinE.count(s) > 0) {
                            min_v = min(min_v, hairpinE[s]);
                        }
                        break;
                    case 7:
                        if (hairpinE.count(s) > 0) {
                            min_v = min(min_v, hairpinE[s]);
                        }
                        break;
                    default:
                        cout << s << " " << l << endl;
                        exit(5);
                        break;
                }

            }

            min_v = min(min_v, hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1));


//            de << "1 " << min_v << endl;

            for (int p = i+1; p <= min(j-5, i+MAXLOOP+1); ++p) {
//                int minq = p + 4;
                int minq = j-i+p-MAXLOOP-2;
                if (minq < p+4) minq = p+4;
                int ll = p-i;
                for (int q = minq; q < j; q++) {
                    if (!basepair(seq[p],seq[q])) continue;
                    int lr = j-q;
                    int internal_min = inf;
                    if (ll == 1 && lr == 1) {
                        internal_min = min(internal_min, stacking(seq[i],seq[j],seq[p],seq[q]));
                    }
                    else if (ll == 1 || lr == 1) {
                        internal_min = min(internal_min, bulge_loop(seq[i],seq[j],seq[p],seq[q],max(ll,lr)-1));
                    }
                    else {
                        internal_min = min(internal_min, interior_loop(seq[i],seq[j],seq[p],seq[q],seq[i+1],seq[j-1],seq[p-1],seq[q+1],ll-1,lr-1));
                    }
                    min_v = min(min_v, internal_min + V[index(p,q)]);
                 }
            }

            min_v = min(min_v, TM[index(i+1,j-1)] + AU[seq[i]][seq[j]] + ML_closing + ML_intern);
//            de << "2 " << min_v << endl;
            V[index(i,j)] = min_v;

//            minV = max(min_v, minV);

//            de << "V - l: " << i << ",r: " << j << ",L: " << seq[i] << ",R: " << seq[j] << ",ret: " << min_v  << endl;


            calculate_WM(i,j);
        }
    }
}

void ZukerAlgorithm::assign_seq2str(string & s, int start) {
    int end = start + s.size();
    for (int i = start; i < end; ++i) {
        s[i-start] = to_char[seq[i]];
    }
}

void ZukerAlgorithm::calculate_WM(int i, int j) {
    int min_multi = inf;

    if (basepair(seq[i],seq[j])) {
        min_multi = min(min_multi, V[index(i,j)] + AU[seq[i]][seq[j]] + ML_intern);
    }

    if (i+1 < j) {
        min_multi = min(min_multi, WM[index(i+1,j)] + ML_BASE);

        min_multi = min(min_multi, WM[index(i,j-1)] + ML_BASE);
    }
//    de << 2 << " " << min_multi << endl;
    int energy_b = inf;
    for (int k = i + 4; k <= j-4; ++k) {
        energy_b = min(energy_b, WM[index(i,k-1)] + WM[index(k,j)]);
    }
    TM[index(i,j)] = min(energy_b, TM[index(i,j)]);
    min_multi = min(min_multi, energy_b);
//    de << 3 << " " << min_multi << endl;
//    de << "M - l: " << i << ",r: " << j << ",L: " << seq[i] << ",R: " << seq[j] << ",ret: " << min_multi  << endl;
    WM[index(i,j)] = min_multi;
}

int ZukerAlgorithm::index(int i, int j) {
    return n*i+j;
}

void ZukerAlgorithm::traceback() {
    int s = 0;
    int t = 0;
    sector[++s].i = 0;
    sector[s].j = n-1;
    sector[s].ml = 0;

    OUTLOOP:
    while (s>0) {
        int fij, fi, cij;
        int traced, i1, j1, k, p, q;

        int i = sector[s].i;
        int j = sector[s].j;
        int ml = sector[s--].ml;

        if (ml == 2) {
            base_pair[++t].i = i;
            if (t > n) exit(3);
            goto repeat1;
        }

        if (j == i) break;

//        cout << "outloop: " << "i: " << i << "j: " << j << endl;

        fij = (ml == 1) ? WM[index(i,j)] : W[index(i,j)];

        fi = (ml == 1) ? WM[index(i,j-1)] + ML_BASE : W[index(i,j-1)];

        if (fij == inf || fij == inf+ML_BASE || fij == inf+ML_intern) {
            goto repeat1;
        }

        if (fij == fi) {
            sector[++s].i = i;
            sector[s].j = j-1;
            sector[s].ml = ml;
//            cout << 3 << " " << "i: " << i << ",j: " << j-1 << " " << fij << endl;
            goto OUTLOOP;
        }

        if (ml == 0) {
//            cout << "w trace" << endl;
            if (basepair(seq[i],seq[j])) {
                int en_c = V[index(i,j)] + AU[seq[i]][seq[j]];
                int en_f = W[index(i,j)];
                if (en_c == en_f) {
                    k = i;
                    traced = j;
                    goto LABEL1;
                }
            }


//            cout << fij << endl;

            for (k = j-4,traced=0; k>=1; k--) {
                int en_c = V[index(k,j)] + AU[seq[k]][seq[j]];
                int en_f = W[index(i,k-1)];
//                cout << fij << " " << en_f << " " << en_c << " " << en_c + en_f << endl;
                if (fij == en_c + en_f) {
                    traced = j;
                    sector[++s].i = i;
                    sector[s].j = k-1;
                    sector[s].ml = 0;
//                    cout << 4 << " " << "i: " << i << ",j: " << k-1 << endl;
                    goto LABEL1;
                }
            }

            LABEL1:
            if (!traced) {
                cout << "backtrack failed in w, i: " << i << ",j: " << j << endl;
            }

            i = k;
            j = traced;
//            cout << 4.1 << "t: " << t << endl;
            base_pair[++t].i = i;
            base_pair[t].j = j;
//            cout << 4 << "t: " << t << endl;
            if (t > n) exit(3);
            goto repeat1;
        }
        else {
//            cout << fij << " " << WM[index(i+1,j)] << " " << WM[index(i+1,j)] + C << endl;
            if (fij == WM[index(i+1,j)] + ML_BASE) {
                sector[++s].i = i+1;
                sector[s].j = j;
                sector[s].ml = ml;
//                cout << 5 << " " << "i: " << i+1 << ",j: " << j << endl;
                goto OUTLOOP;
            }

//            cout << i << " " << j << " " << fij << " " << V[index(i,j)] << " " << V[index(i,j)] + B << endl;
            if (fij == V[index(i,j)] + AU[seq[i]][seq[j]] + ML_intern) {
//                cout << 5.1 << "t: " << t << endl;
                base_pair[++t].i = i;
                base_pair[t].j = j;
//                cout << 5 << "t: " << t << endl;
                if (t > n) exit(3);
                goto repeat1;
            }

            for (k = i+4; k <= j-4; k++){
                if (fij == WM[index(i,k-1)] + WM[index(k,j)]) {
                    sector[++s].i = i;
                    sector[s].j = k-1;
                    sector[s].ml = ml;
//                    cout << 6 << " " << "i: " << i << ",j: " << k-1 << endl;
                    sector[++s].i = k;
                    sector[s].j = j;
                    sector[s].ml = ml;
//                    cout << 7 << " " << "i: " << k << ",j: " << j << endl;
                    goto OUTLOOP;
                }
            }

            if (k>j-4) {
                cout << "backtrack failed in fML,i: " << i << ",j: " << j << endl;
                exit(1);
            }
        }

        repeat1:
        cij = V[index(i,j)];
//        cout << 14 << endl;
        int l = j-i;

//        cout << "repeat: " << "i: " << i << ",j: " << j << " " << cij << endl;

        if (l+1 == 5 || l+1 == 6 || l+1 == 8) {
            int index;
            string str(l+1,'.');
            assign_seq2str(str, i);
            if (l+1 == 5) {
                index = get_index(TriloopSeq, str);
                if (index != -1) {
                    if (cij == Triloop[index]) goto OUTLOOP;
                } else {
                    if (cij == hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1)) goto OUTLOOP;
                }
            }
            else if (l+1 == 6) {
                index = get_index(TetraloopSeq, str);
                if (index != -1) {
                   if (cij == Tetraloop[index]) goto OUTLOOP;
                } else {
                    if (cij == hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1)) goto OUTLOOP;
                }
            }
            else {
                index = get_index(HexaloopSeq, str);
                if (index != -1 ) {
                    if (cij == Hexaloop[index]) goto OUTLOOP;
                } else {
                    if (cij == hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1)) goto OUTLOOP;
                }
            }
        } else {
//            cout << hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1) << endl;
            if (cij == hairpin_loop(seq[i],seq[j],seq[i+1],seq[j-1],l-1)) goto OUTLOOP;
        }

        for (p = i+1; p <= min(j-5,i+MAXLOOP+1); ++p) {
            int minq = j-i+p-MAXLOOP-2;
            if (minq < p+4) minq = p+4;
            int ll = p-i;
            for (q = j-1; q >= minq; q--) {
                if (!basepair(seq[p],seq[q])) continue;
                int lr = j-q;
//                cout << "p: " << p << ",q: " << q << ",ll: " << ll << ",lr: " << lr << endl;
                int internal = cij - V[index(p,q)];
//                cout << 1 << endl;
                if (ll == 1 && lr == 1) {
//                    cout << 2 << endl;
                    if (internal == stacking(seq[i],seq[j],seq[p],seq[q])) {
//                        cout << 7.1 << "t: " << t << endl;
                        base_pair[++t].i = p;
                        base_pair[t].j = q;

//                        cout << 7 << "t: " << t << endl;
                        if (t > n) exit(3);
                        i = p, j = q;
//                        cout << 9 << " " << "i: " << i << ",j: " << j << endl;
                        goto repeat1;
                    }
                }
                else if (ll == 1 || lr == 1) {
                    int bulge_energy;
                    if (ll == 1) {
                        bulge_energy = bulge_loop(seq[i],seq[j],seq[p],seq[q],lr-1);
                    } else {
                        bulge_energy = bulge_loop(seq[i],seq[j],seq[p],seq[q],ll-1);
                    }
                    if (internal == bulge_energy) {
//                        cout << 1.1 << "t: " << t << endl;
                        base_pair[++t].i = p;
                        base_pair[t].j = q;

//                        cout << 1 << "t: " << t << endl;
                        if (t > n) exit(3);

                        i = p, j = q;
//                        cout << 10 << " " << "i: " << i << ",j: " << j << endl;
                        goto repeat1;
                    }
                }
                else {
                    int interior = interior_loop(seq[i],seq[j],seq[p],seq[q],seq[i+1],seq[j-1],seq[p-1],seq[q+1],ll-1,lr-1);
                    if (internal == interior) {
//                        cout << 2.1 << "t: " << t << endl;
                        base_pair[++t].i = p;
                        base_pair[t].j = q;
//                        cout << 2 << "t: " << t << endl;
                        if (t > n) exit(3);

                        i = p, j = q;
                        goto repeat1;
                    }
                }

            }
        }

        i1 = i+1;
        j1 = j-1;
        sector[s+1].ml = sector[s+2].ml = 1;

        int en = cij - AU[seq[i]][seq[j]] - ML_intern - ML_closing;



        for (k = i+5; k < j-4; k++) {
//            cout << en << " " << WM[index(i+1,k-1)] + WM[index(k,j-1)] << endl;
            if (en == WM[index(i+1,k-1)] + WM[index(k,j-1)]) {
                goto LABEL2;
            }
        }

        LABEL2:
        if (k <= j-5) {
            sector[++s].i = i1;
            sector[s].j = k-1;
//            cout << 1 << " " << "i: " << i1 << ",j: " << k-1 << endl;
            sector[++s].i = k;
            sector[s].j = j1;
//            cout << 2 << " " << "i: " << k << ",j: " << j1 << endl;

        }
        else {
            cout << "backtracking failed in repeat, i: " << i << ",j: " << j << ",ml: " << ml << endl;
            exit(2);
        }
    }

    base_pair[0].i = t;
}

void ZukerAlgorithm::get_bp(string & bp) {

//    cout << "base pair number: " << base_pair[0].i << endl;
    for (int a = 1; a <= base_pair[0].i; ++a) {
        bp[base_pair[a].i] = '(';
        bp[base_pair[a].j] = ')';
    }
//    cout << bp << endl;
}
