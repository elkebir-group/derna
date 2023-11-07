//
// Created by xinyu on 5/18/2022.
//

#ifndef RNA_DESIGN_ZUKER_H
#define RNA_DESIGN_ZUKER_H
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <cmath>
#include "default.h"
#include "utils.h"

using namespace std;

class Zuker {

    vector<int> Z;
    vector<int> E;
    vector<int> M;
    vector<int> TM0;
    vector<double> Z2,E2,M2,TM2;
    vector<double> Z_CAI,E_CAI,M_CAI,TM_CAI;
    vector<double> O;
    vector<double> E1;
    vector<double> M1;
    vector<double> TM;
    vector<int> OB, EB, MB;
    vector<int> protein;
    vector<bond> bp_bond;
    vector<int> nucle_seq;
    vector<int> basepair;
    vector<stack_> sector;
    vector<int> ava_nucle_p;
    vector<int> ava_nucle_m;
    vector<int> codon_selection;
    vector<int> start_index, index_offset;
    int n, minX, minY, last_idx;
    int e0;
    double e;
    vector<vector<int>> O_bt, E_bt, M_bt, TM_bt;


public:
    Zuker(int n, int,vector<int> &);
    Zuker(const Zuker &);
    ~Zuker();
    void init_values();

    /**
     * Fill the vector Z where Z[i] stores the minimum free energy at i
     *
     * @param ostream object for stdout
     * */
    int calculate_Z(ostream &);


    /**
     * Fill the vector E where E[i] stores the minimum free energy at i
     * if the two positions mapped to i are eligible base pairs
     *
     * @param ostream object for stdout
     * */
    void calculate_E();



    /**
     * Fill the vector M where M[i] stores the minimum free energy of the
     * multiloop formed between the two positions mapped to i
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param ostream object for stdout
     * */
    int calculate_M(int a, int b, int i, int j, int x, int y); //, ostream &



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

    /**
     * Find RNA sequence and secondary structure
     * */
    void traceback();

    void traceback_B();

    /**
     * Find RNA sequence and secondary structure
     *
     * @param lambda is a weight between 0 and 1
     * */


    void traceback_B2(double lambda);

    void traceback_MFE(double lambda);

    /**
     * Returns RNA secondary structure
     *
     * @param bp is a string to store the secondary structure
     * */
    void get_bp(string & bp);

    /**
     * Returns RNA sequence
     *
     * @param rna is a string to store the sequence
     * */
    void get_rna(string & rna);

    /**
     * Returns RNA sequence when using both objectives (free energy and codon adaptation index)
     *
     * @param rna is a string to store the sequence
     * */
    void get_rna_cai(string & rna);

    int fill_rna(int i);

    /**
     * Returns RNA sequence with 'X' at positions where the nucleotides can be replaced
     *
     * @param rna is a string to store the sequence
     * */
    void get_rna_X(string & rna);

    /**
     * Returns all possible RNA sequences by replacing 'X' with eligible nucleotides
     *
     * @param rna_array is vector that stores all possible RNA sequences
     * @param index current index in array
     * @param i RNA nucleotide position
     * */
    void getAllRna(vector<string> & rna_array,int index,int i);

    /**
     * Perform lambda swipe on a given input.
     *
     * @param incr is the increment on lambda each loop [0,1]
     * @param fout ostream object for stdout
     * @param outfile output file to save swipe data, csv format
     */
    void lambda_swipe(double incr, ostream& fout, string & outfile);

    /**
     * Lambda swipe until finding points p,q where p.MFE = q.MFE and p.CAI = q.CAI
     * @param incr
     * @param fout
     * @param outfile
     */
    void lambda_swipe_2(double, double, ostream& fout, string & outfile);

    void lambda_swipe_3(double, double, ostream& fout, string & outfile);

    /**
     * Access vector Z given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline int &Access_Z(int a, int b, int i, int j, int x, int y);

    /**
     * Access vector E given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline int &Access_E(int a, int b, int i, int j, int x, int y);

    /**
     * Access vector M given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline int &Access_M(int a, int b, int i, int j, int x, int y);

    /**
     * Access vector O given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline double &Access_O(int a, int b, int i, int j, int x, int y);

    inline int &Access_OB(int a, int b, int i, int j, int x, int y);
    inline int &Access_EB(int a, int b, int i, int j, int x, int y);
    inline int &Access_MB(int a, int b, int i, int j, int x, int y);

    /**
     * Access vector E1 given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline double &Access_E1(int a, int b, int i, int j, int x, int y);

    /**
     * Access vector M1 given protein position, codon selection and codon index (a, b, i, j, x, y)
     *
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * */
    inline double &Access_M1(int a, int b, int i, int j, int x, int y);


    inline double &Access_Z2(int a, int b, int i, int j, int x, int y);
    inline double &Access_M2(int a, int b, int i, int j, int x, int y);
    inline double &Access_E2(int a, int b, int i, int j, int x, int y);
    inline double &Access_Z2(int index);
    inline double &Access_E2(int index);
    inline double &Access_M2(int index);



    /**
     * Access vector M given index
     *
     * @param index precomputed location of stored value in M
     * */
    inline int &Access_M(int index);

    /**
     * Access vector E given index
     *
     * @param index precomputed location of stored value in E
     * */
    inline int &Access_E(int index);

    /**
     * Access vector Z given index
     *
     * @param index precomputed location of stored value in Z
     * */
    inline int &Access_Z(int index);

    /**
     * Access vector O given index
     *
     * @param index precomputed location of stored value in O
     * */
    inline double &Access_O(int index);

    inline int &Access_OB(int index);
    inline int &Access_EB(int index);
    inline int &Access_MB(int index);

    /**
     * Access vector E1 given index
     *
     * @param index precomputed location of stored value in E1
     * */
    inline double &Access_E1(int index);

    /**
     * Access vector M1 given index
     *
     * @param index precomputed location of stored value in M1
     * */
    inline double &Access_M1(int index);

//    inline int ava_nucleotides_int(int, int, int, int, int);

private:
    /**
     * Return folding energy of hairpin loop structure
     *
     * @param xi  paired left end nucleotide of the hairpin
     * @param yj  paired right end nucleotide of the hairpin
     * @param xi_ the first unpaired left nucleotide in the loop
     * @param _yj the first unpaired right nucleotide in the loop
     * @param l   loop length
     * @param fout ostream object for stdout
     * @return folding energy of hairpin loop structure
     */
    static inline int hairpin_loop(int xi, int yj, int xi_, int _yj, int l);

    /**
     * Return folding energy of Coaxial Stacking
     *
     * @param i paired bottom left end nucleotide of the stacking
     * @param j paired bottom right end nucleotide of the stacking
     * @param i1 paired top left end nucleotide of the stacking
     * @param j1 paired top right end nucleotide of the stacking
     * @return folding energy of Coaxial Stacking
     */
    static inline int stacking(int i, int j, int i1, int j1);

    /**
     * Return folding energy of Bulge Loops
     *
     * @param i paired bottom left end nucleotide of the bulge loop
     * @param j paired bottom right end nucleotide of the bulge loop
     * @param h paired top left end nucleotide of the bulge loop
     * @param k paired top right end nucleotide of the bulge loop
     * @param l loop length
     * @return folding energy of Bulge Loops
     */
    static inline int bulge_loop(int i, int j, int h, int k, int l);

    /**
     * Return folding energy of Internal Loops
     * @param i paired bottom left end nucleotide of the bulge loop
     * @param j paired bottom right end nucleotide of the bulge loop
     * @param h paired top left end nucleotide of the bulge loop
     * @param k paired top right end nucleotide of the bulge loop
     * @param i1 the first unpaired bottom left nucleotide in the loop
     * @param j1 the first unpaired bottom right nucleotide in the loop
     * @param h1 the first unpaired top left nucleotide in the loop
     * @param k1 the first unpaired top right nucleotide in the loop
     * @param n1 left loop length
     * @param n2 right loop length
     * @param fout ostream object for stdout
     * @return folding energy of Internal Loops
     */
    static inline int interior_loop(int i, int j, int h, int k, int i1, int j1, int h1, int k1, int n1, int n2);


    inline int Access_basepair(int, int);

    /**
     * Return available nucleotides at relative position
     * @param a is the position index of the amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param dir -1 or 1: if -1 return available nucleotides at position before the given position,
     * else return available nucleotides at position after the given position
     * @return available nucleotides at relative position
     */
    inline int ava_nucleotides_int(int a, int x, int i, int dir);
//    int ava_nucleotides_int(int a, int x, int i, int dir);

    static inline int H_index(int, int, int, int);

    static inline int H_index(int, int, int);
    static inline int T_index(int, int, int, int, int);

    /**
     * Return minimum energy among stacking, bulge, and internal loops structures
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param la position index of the left nucleotide in the RNA seq
     * @param lb position index of the right nucleotide in the RNA seq
     * @param xi paired left end nucleotide of the structure
     * @param yj paired right end nucleotide of the structure
     * @param an_int available nucleotides at position after la
     * @param bp_int available nucleotides at position before lb
     * @param fout ostream object for stdout
     * @return minimum energy among stacking, bulge, and internal loops structures
     */
    int internal(int a, int b, int i, int j, int x, int y, int la, int lb, int xi, int yj, int an_int, int bp_int);


    /**
     * Return minimum energy among stacking, bulge, and internal loops structures
     * integrating with codon adaptation index
     *
     * @param lambda lambda value, 0 <= lambda <= 1
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param la position index of the left nucleotide in the RNA seq
     * @param lb position index of the right nucleotide in the RNA seq
     * @param xi paired left end nucleotide of the structure
     * @param yj position index of the right nucleotide in the RNA seq
     * @param an_int available nucleotides at position after la
     * @param bp_int available nucleotides at position before lb
     * @param fout ostream object for stdout
     * @return minimum energy among stacking, bulge, and internal loops structures
     * integrating with codon adaptation index
     */
    double internal_CAI(double lambda, int a, int b,int i, int j, int x, int y, int la, int lb, int xi, int yj);


    /**
     * Return minimum hairpin folding energy
     * @param la position index of the left nucleotide in the RNA seq
     * @param lb position index of the right nucleotide in the RNA seq
     * @param l  loop length
     * @param pa the amino acid at position a
     * @param pb the amino acid at position b
     * @param pna the amino acid to the right of position a
     * @param ppb the amino acid to the left of position b
     * @param n_codon_an number of codons available at position a+1
     * @param n_codon_bp number of codons available at position b-1
     * @param xi paired left end nucleotide of the hairpin
     * @param yj paired right end nucleotide of the hairpin
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param an_int available nucleotides at position after la
     * @param bp_int available nucleotides at position before lb
     * @param fout ostream object for stdout
     * @return minimum hairpin folding energy among all possible sequences between la and lb
     */
    int hairpin(int la, int lb, int l, int pa, int pb, int pna, int ppb, int n_codon_an, int n_codon_bp, int xi, int yj, int i, int j, int x, int y, int an_int, int bp_int) const;
    int hairpin(int,int,int,int,int,int,int) const;

    /**
     * Return minimum folding energy of multiloop structure
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param pa the amino acid at position a
     * @param pb the amino acid at position b
     * @param n_codon_an number of codons available at position a+1
     * @param n_codon_bp number of codons available at position b-1
     * @param fout ostream object for stdout
     * @param TM2 vector to store minimum energy of bifurication in multiloop structure
     * @return minimum folding energy of multiloop structure
     */
    int multi_loop(int a, int b, int i, int j, int x, int y, int pa, int pb, int n_codon_an, int n_codon_bp);
    int multi_loop_Mem(int a, int b, int i, int j, int x, int y, int pa, int pb, int n_codon_an, int n_codon_bp);
    double multi_loop_CAI(double,int,int,int,int,int,int,int,int,int,int);



    /**
     * Return minimum energy of hairpin structure
     * @param lambda lambda value, 0 <= lambda <= 1
     * @param l loop length
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param pa the amino acid at position a
     * @param pb the amino acid at position b
     * @param pna the amino acid to the right of position a
     * @param ppb the amino acid to the left of position b
     * @param n_codon_an number of codons available at position a+1
     * @param n_codon_bp number of codons available at position b-1
     * @param xi paired left end nucleotide of the hairpin
     * @param yj paired right end nucleotide of the hairpin
     * @param i is the position index in codon x, 0 <= i <= 2
     * @param j is the position index in codon y, 0 <= j <= 2
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @param an_int available nucleotides at position after la
     * @param bp_int available nucleotides at position before lb
     * @param fout ostream object for stdout
     * @return minimum energy of hairpin structure among all possible sequences between la and lb, considering CAI
     */
    double hairpin_CAI(double lambda, int l,int a, int b, int pa, int pb, int pna, int ppb, int n_codon_an, int n_codon_bp, int xi, int yj, int i, int j, int x, int y);

    /**
     * Return codon adptation index compensation for hairpin structure
     * @param a is the position index of the left amino acid in the protein seq
     * @param b is the position index of the right amino acid in the protein seq
     * @param x is a possible codon of the amino acid at position a
     * @param y is a possible codon of the amino acid at position b
     * @return codon adaptation index compensation for hairpin structure
     */

    double add_hairpin_CAI_2(int a, int b,  int x, int y, int a1=-1, int x1=-1,int b1 = -1, int y1 = -1) const;

    double add_hairpin_CAI_8(int a, int b,  int x, int y, int a1=-1, int x1=-1,int b1 = -1, int y1 = -1) const;

    double add_hairpin_CAI_3(vector<int> & s, int ) const;

    void assign_codon(vector<int> & s, int );

    double add_CAI(int,int,int) const;

    double add_interior_CAI_2(int a, int c, int x, int na=-1, int x1=-1, int pc=-1, int h1=-1) const;
    void static assign(vector<int> &, vector<int> &, int);

    /**
     * Return the max codon adaptation index at given position with given constraints
     * @param p is the position index of the amino acid in the protein seq
     * @param l0 is the first nucleotide in the amino acid, default -1
     * @param l1 is the second nucleotide in the amino acid, default -1
     * @param l2 is the third nucleotide in the amino acid, default -1
     * @return the max codon adaptation index at given position with given constraints
     */


    inline int index(int,int,int,int,int,int) const;

    inline int index(int,int,int) const;

    /**
     * Check if two nucleotides can be in the same codon
     * @param l1 position index of the left nucleotide in the RNA seq
     * @param l2 position index of the right nucleotide in the RNA seq
     * @param x left nucleotide
     * @param y right nucleotide
     * @return if two nucleotides can be in the same codon
     */
    bool rightCodon(int l1, int l2, int x, int y) const;
    void reinit();
    void maxCAISeq();
};

inline int Zuker::ava_nucleotides_int(int a, int x, int i, int dir) {
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


// nucleotide A = 0, C = 1, G = 2, U = 3
inline int Zuker::Access_basepair(int X, int Y)
{
    return basepair[(X<<2)+Y];
}

inline int Zuker::interior_loop(int i, int j, int h, int k, int i1, int j1, int h1, int k1, int n1, int n2) {
    int energy ;
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[h+1][k+1]];
    int nl, ns;

    nl = max(n1, n2);
    ns = min(n1, n2);

    if (ns==1) {
        if (nl==1)  {
            return int11[type][type2][i1+1][j1+1];
        }                   /* 1x1 loop */

        if (nl==2) {                  /* 2x1 loop */
            if (n1==1) {
                energy = int21[type][type2][i1+1][k1+1][j1+1];
            }
            else {
                energy = int21[type2][type][k1+1][i1+1][h1+1];
            }
            return energy;
        }
        else {  /* 1xn loop */
            energy = (nl+1<=MAXLOOP)?(internal_loop[nl+1]) : (internal_loop[30]+(int)(lxc*Log[nl+1]));//log((nl+1)/30.)
            energy += min(MAX_NINIO, (nl-ns)*ninio);
            energy += mismatch1nI[type][i1+1][j1+1] + mismatch1nI[type2][k1+1][h1+1];
            return energy;
        }
    }
    else if (ns==2) {
        if(nl==2)      {              /* 2x2 loop */
            return int22[type][type2][i1+1][h1+1][k1+1][j1+1];
        }
        else if (nl==3){              /* 2x3 loop */
            energy = internal_loop[5]+ninio;
            energy += mismatch23I[type][i1+1][j1+1] + mismatch23I[type2][k1+1][h1+1];
            return energy;
        }

    }
    { /* generic interior loop (no else here!)*/
        energy = (n1+n2<=MAXLOOP)?(internal_loop[n1+n2]) : (internal_loop[30]+(int)(lxc*Log[n1+n2]));//log((n1+n2)/30.)
        energy += min(MAX_NINIO, (nl-ns)*ninio);
        energy += mismatchI[type][i1+1][j1+1] + mismatchI[type2][k1+1][h1+1];


    }
    return energy;
}

inline int Zuker::bulge_loop(int i, int j, int h, int k, int l) {
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[h+1][k+1]];
    int bulge_energy;
    // add penalty based on size
    bulge_energy = (l <= MAXLOOP) ? bulge[l] : bulge[30]+(int)(lxc*Log[l]); //log((l)/30.)

    if (l == 1) {
        bulge_energy += stackE[type][type2];
    } else {
        // add penalty for AU terminal
        bulge_energy += AU[i][j];
        bulge_energy += AU[k][h];
    }

    return bulge_energy;
}

inline int Zuker::stacking(int i, int j, int i1, int j1) {
    int type = BP_pair[i+1][j+1];
    int type2 = rtype[BP_pair[i1+1][j1+1]];
    int stacking_energy = stackE[type][type2];//stacks[stack_index(i, j1, j, i1)];
    return stacking_energy;
}

inline int Zuker::hairpin_loop(int xi, int yj, int xi_, int _yj, int l) {
    int hairpin_energy;
    int type = BP_pair[xi+1][yj+1];

    // add penalty based on size
    hairpin_energy = (l <= 30) ? hairpins[l] : hairpins[30]+(int)(lxc*Log[l]);//hairpinLoops[l]; // log((l)/30.)

    if (l == 3) return hairpin_energy + AU[xi][yj];

    // add penalty for a terminal mismatch
    hairpin_energy += mismatchH[type][xi_+1][_yj+1];//T_mm[stack_index(xi, _yj, yj, xi_)];

//    hairpin_energy += add_auterminal(xi, yj,tempf);
    return hairpin_energy;
}

// TODO: precompute the indexes. offset in the loop rather than recomputing index.
inline int & Zuker::Access_Z(int a, int b, int i, int j, int x, int y) {
    return Z[index(a,b,i,j,x,y)];
}

inline int & Zuker::Access_E(int a, int b, int i, int j, int x, int y) {
    return E[index(a,b,i,j,x,y)];
}

inline int &Zuker::Access_M(int a, int b, int i, int j, int x, int y) {
    return M[index(a,b,i,j,x,y)];
}

inline int &Zuker::Access_M(int index) {

    return M[index];
}

inline int &Zuker::Access_Z(int index) {

    return Z[index];
}

inline int &Zuker::Access_E(int index) {

    return E[index];
}

inline double &Zuker::Access_O(int index) {

    return O[index];
}

inline int &Zuker::Access_OB(int index) {

    return OB[index];
}

inline int &Zuker::Access_EB(int index) {

    return EB[index];
}

inline int &Zuker::Access_MB(int index) {

    return MB[index];
}

inline double &Zuker::Access_E1(int index) {

    return E1[index];
}

inline double &Zuker::Access_M1(int index) {

    return M1[index];
}

inline double &Zuker::Access_Z2(int index) {

    return Z2[index];
}

inline double &Zuker::Access_E2(int index) {

    return E2[index];
}

inline double &Zuker::Access_M2(int index) {

    return M2[index];
}

inline int & Zuker::Access_OB(int a, int b, int i, int j, int x, int y) {
    return OB[index(a,b,i,j,x,y)];
}

inline int & Zuker::Access_EB(int a, int b, int i, int j, int x, int y) {
    return EB[index(a,b,i,j,x,y)];
}

inline int & Zuker::Access_MB(int a, int b, int i, int j, int x, int y) {
    return MB[index(a,b,i,j,x,y)];
}

inline double & Zuker::Access_O(int a, int b, int i, int j, int x, int y) {
    return O[index(a,b,i,j,x,y)];
}

inline double & Zuker::Access_E1(int a, int b, int i, int j, int x, int y) {
    return E1[index(a,b,i,j,x,y)];
}

inline double &Zuker::Access_M1(int a, int b, int i, int j, int x, int y) {
    return M1[index(a,b,i,j,x,y)];
}

inline double &Zuker::Access_Z2(int a, int b, int i, int j, int x, int y) {
    return Z2[index(a,b,i,j,x,y)];
}

inline double &Zuker::Access_E2(int a, int b, int i, int j, int x, int y) {
    return E2[index(a,b,i,j,x,y)];
}

inline double &Zuker::Access_M2(int a, int b, int i, int j, int x, int y) {
    return M2[index(a,b,i,j,x,y)];
}

inline int Zuker::index(int a, int b, int i, int j, int x, int y) const {
    int idx = a*n+b;
//    cout << start_index[idx] + index_offset[idx] * (3*i+j) + n_codon[protein[b]] * x + y << endl;
    return start_index[idx] + index_offset[idx] * (3*i+j) + n_codon[protein[b]] * x + y;
}

inline int Zuker::index(int a, int x, int i) const {
    return 18*a+3*x+i;
}

inline int Zuker::H_index(int n1, int n2, int n3, int n4) {
    return (n1<<6)+(n2<<4)+(n3<<2)+n4;
}

inline int Zuker::H_index(int type, int n1, int n2) {
    return ((type-1)<<4)+(n1<<2)+n2;
}


inline int Zuker::T_index(int type, int n1, int n2, int n3, int n4) {
    return ((type-1)<<8)+(n1<<6)+(n2<<4)+(n3<<2)+n4; //
}

#endif //RNA_DESIGN_ZUKER_H
