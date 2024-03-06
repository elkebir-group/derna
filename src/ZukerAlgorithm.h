//
// Created by xinyu on 4/26/2022.
//

#ifndef RNA_DESIGN_ZUKERALGORITHM_H
#define RNA_DESIGN_ZUKERALGORITHM_H

#include <string>
#include <vector>
#include <cmath>
#include "utils.h"
#include "default.h"


using namespace std;

class ZukerAlgorithm
{
    vector<int> & seq;
    string bp;
    int n;
    vector<int> W;
    vector<int> V;
    vector<int> WM;
    vector<int> TM;
    vector<int> WB;
    vector<int> VB;
    vector<int> MB;
    vector<st> sector;
    vector<bond> base_pair;
    vector<char> to_char;


public:
    ZukerAlgorithm(vector<int> & seq, int n);
    ZukerAlgorithm(const ZukerAlgorithm &);
    ~ZukerAlgorithm();
    void calculate_V();
    int calculate_W();
    void calculate_WM(int, int);
    inline int hairpin_loop(int, int, int,int,int) const;
    inline int stacking(int, int, int,int) const;
    inline int bulge_loop(int, int, int, int, int) const;
    inline int interior_loop(int, int, int, int,int,int,int,int,int,int) const;
    int index(int, int);
    void traceback();
    void traceback_2();
    void get_bp(string &);

private:
    void assign_seq2str(string &, int);
};

inline int ZukerAlgorithm::bulge_loop(int i, int j, int h, int k, int l) const {
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[h+1][k+1]];
    int bulge_energy;
    // add penalty based on size
    bulge_energy = (l <= MAXLOOP) ? bulge[l] : bulge[30]+(int)(lxc*Log[l]);

    // if len 1, include the delta G of intervening NN (SantaLucia 2004)
    if (l == 1) {
        bulge_energy += stackE[type][type2];//stacks[stack_index(i, k, j, h)];
    } else {
        bulge_energy += AU[i][j];
        bulge_energy += AU[k][h];
    }

    // add penalty for AU terminal
    return bulge_energy;
}

inline int ZukerAlgorithm::stacking(int i, int j, int p, int q) const {
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[p+1][q+1]];
    int stacking_energy = stackE[type][type2];//stacks[stack_index(i, j1, j, i1)];
    return stacking_energy;

}

inline int ZukerAlgorithm::interior_loop(int i, int j, int h, int k, int i1, int j1, int h1, int k1, int n1, int n2) const {
    int energy;
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[h+1][k+1]];
    int nl, ns;

    nl = max(n1, n2);
    ns = min(n1, n2);
    if (ns==1) {
        if (nl==1)                    /* 1x1 loop */
            return int11[type][type2][i1+1][j1+1];
        if (nl==2) {                  /* 2x1 loop */
            if (n1==1)
                energy = int21[type][type2][i1+1][k1+1][j1+1];
            else
                energy = int21[type2][type][k1+1][i1+1][h1+1];
            return energy;
        }
        else {  /* 1xn loop */
            energy = (nl+1<=MAXLOOP)?(internal_loop[nl+1]) : (internal_loop[30]+(int)(lxc*Log[nl+1]));
            energy += min(MAX_NINIO, (nl-ns)*ninio);
            energy += mismatch1nI[type][i1+1][j1+1] + mismatch1nI[type2][k1+1][h1+1];
            return energy;
        }
    }
    else if (ns==2) {
        if(nl==2)      {              /* 2x2 loop */
            return int22[type][type2][i1+1][h1+1][k1+1][j1+1];}
        else if (nl==3){              /* 2x3 loop */
            energy = internal_loop[5]+ninio;
            energy += mismatch23I[type][i1+1][j1+1] + mismatch23I[type2][k1+1][h1+1];
            return energy;
        }

    }
    { /* generic interior loop (no else here!)*/
        energy = (n1+n2<=MAXLOOP)?(internal_loop[n1+n2]) : (internal_loop[30]+(int)(lxc*Log[n1+n2]));
        energy += min(MAX_NINIO, (nl-ns)*ninio);

        energy += mismatchI[type][i1+1][j1+1] + mismatchI[type2][k1+1][h1+1];
    }
    return energy;
}

inline int ZukerAlgorithm::hairpin_loop(int xi, int yj, int xi_, int _yj, int l) const {
    int hairpin_energy;
    int type = BP_pair[xi+1][yj+1];

    // add penalty based on size
    hairpin_energy = (l <= 30) ? hairpins[l] : hairpins[30]+(int)(lxc*Log[l]);//hairpinLoops[l];

    if (l == 3) return hairpin_energy + AU[xi][yj];

    // add penalty for a terminal mismatch
    hairpin_energy += mismatchH[type][xi_+1][_yj+1];//T_mm[stack_index(xi, _yj, yj, xi_)];

//    hairpin_energy += add_auterminal(xi, yj,tempf);
    return hairpin_energy;

}

#endif //RNA_DESIGN_ZUKERALGORITHM_H
