//
// Created by Summer Gu on 6/4/25.
//

#ifndef DERNA_BEAMZUKER_H
#define DERNA_BEAMZUKER_H
#include "default.h"
#include "utils.h"

class BeamZuker {
    unordered_map<int, BeamEntry> O, E, M, TM;
//    unordered_map<int, BeamEntry> Z2,E2,M2,TM2;
    vector<int> protein, nucle_seq, basepair;
    vector<int> start_index, index_offset;
    vector<int> ava_nucle_p, ava_nucle_m, codon_selection;
    vector<bond> bp_bond;
    vector<stack_> sector;
    int n, k, minX, minY, last_idx;
    int e0;
    double e;

public:
    BeamZuker(int n,vector<int> &, int k = 10);
    void init_values();

    inline int index(int,int,int) const;
    inline int index(int,int,int,int,int,int) const;
    inline int ava_nucleotides_int(int a, int x, int i, int dir);


    /**
     * Fill the vector O where O[i] stores the minimum value at i
     * combining free energy and codon adaptation index
     *
     * @param ostream object for stdout
     * @param lambda is a weight between 0 and 1
     * */
    double calculate_CAI_O(ostream &, double lambda);

    /**
     * Fill the vector E1 where E1[i] stores the minimum value at i
     * if the two positions mapped to i are eligible base pairs
     *
     * @param ostream object for stdout
     * @param lambda is a weight between 0 and 1
     * */
    void calculate_CAI_E(double lambda);

    /**
     * Fill the vector M1 where M1[i] stores the minimum value of the
     * multiloop formed between the two positions mapped to i
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param lambda is a weight between 0 and 1
     * @param ostream object for stdout
     * */
    void calculate_CAI_M(int a, int b, int i, int j, int x, int y, double lambda); //vector<double> &

};

inline int BeamZuker::index(int a, int x, int i) const {
    return 18*a+3*x+i;
}

inline int BeamZuker::index(int a, int b, int i, int j, int x, int y) const {
    int idx = a*n+b;
    return start_index[idx] + index_offset[idx] * (3*i+j) + n_codon[protein[b]] * x + y;
}

inline int BeamZuker::ava_nucleotides_int(int a, int x, int i, int dir) {
    int s = 0;
    if (dir == 1) {
        if (i <= 1) {
            s |= (1 << nucleotides[protein[a]][x][i+1]);
        }
        else {
            int pna = protein[a+1];
            int an = n_codon[pna];
            for (int x1 = 0; x1 < an; ++x1) {
                s |= (1 << nucleotides[pna][x1][0]);
            }

        }
    }
    else {
        if (i >= 1) {
            s |= (1 << nucleotides[protein[a]][x][i-1]);
        }
        else {
            int ppa = protein[a-1];
            int ap = n_codon[ppa];
            for (int x1 = 0; x1 < ap; ++x1) {
                s |= (1 << nucleotides[ppa][x1][2]);
            }
        }
    }

    return s;
}
#endif //DERNA_BEAMZUKER_H
