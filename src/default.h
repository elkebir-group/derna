//
// Created by Summer Gu on 8/16/22.
//

#ifndef RNA_DESIGN_DEFAULT_H
#define RNA_DESIGN_DEFAULT_H

#include <unordered_map>
#include <fstream>
using namespace std;

#include "params/constants.h"
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)
#define PUBLIC



typedef struct st {
    int i; // index i
    int j; // index j
    int ml; //?
} st;

typedef struct stack_ {
    int a;
    int b;
    int i; // index i
    int j; // index j
    int x; // nucleotide at i
    int y; // nucleotide at j
    int ml; //?
} stack_;

typedef struct bond {
    int i;
    int j;
} bond;

extern vector<double> Log;


extern double lxc37;   /* parameter for logarithmic loop
			  energy extrapolation            */

extern int BP_pair[5][5];
extern int rtype[7];
extern int nucleotides[20][6][3];
extern bool codonMatch[20][3][3][4][4];
extern double codon_cai[20][6];
extern double codon_cai_s[20][6];
extern double max_codon_cai[20][6][16];
extern int max_cai_pos[20];
extern int AU[4][4];
extern vector<char> to_char;
extern vector<int> n_codon;
extern int internal_loop_map[3][3];


extern double lxc;
extern int stackE[NBPAIRS+1][NBPAIRS+1];
extern int hairpins[31];
extern int bulge[31];
extern int internal_loop[31];
extern int mismatchI[NBPAIRS+1][5][5];
extern int mismatch1nI[NBPAIRS+1][5][5];
extern int mismatch23I[NBPAIRS+1][5][5];
extern int mismatchH[NBPAIRS+1][5][5];
extern int int11[NBPAIRS+1][NBPAIRS+1][5][5];
extern int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
extern int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
extern int ML_BASE;
extern int ML_closing;
extern int ML_intern;
extern int ninio;
extern int TerminalAU;

extern int stack37[NBPAIRS+1][NBPAIRS+1];
extern int stackdH[NBPAIRS+1][NBPAIRS+1]; /* stack enthalpies */

extern int hairpin37[31];
extern int hairpindH[31];
extern int bulge37[31];
extern int bulgedH[31];
extern int internal_loop37[31];
extern int internal_loopdH[31];
extern int mismatchI37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatchIdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatch1nI37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatch23I37[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatch1nIdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatch23IdH[NBPAIRS+1][5][5];  /* interior loop mismatches */
extern int mismatchH37[NBPAIRS+1][5][5];  /* same for hairpins */
extern int mismatchM37[NBPAIRS+1][5][5];  /* same for multiloops */
extern int mismatchHdH[NBPAIRS+1][5][5];  /* same for hairpins */
extern int mismatchMdH[NBPAIRS+1][5][5];  /* same for multiloops */
extern int mismatchExt37[NBPAIRS+1][5][5];
extern int mismatchExtdH[NBPAIRS+1][5][5];
extern int interiorLoop[4][4][4][4][4][4][4][4][MAXLOOP][MAXLOOP];

extern int dangle5_37[NBPAIRS+1][5];      /* 5' dangle exterior of pair */
extern int dangle3_37[NBPAIRS+1][5];      /* 3' dangle */
extern int dangle3_dH[NBPAIRS+1][5];       /* corresponding enthalpies */
extern int dangle5_dH[NBPAIRS+1][5];

extern int int11_37[NBPAIRS+1][NBPAIRS+1][5][5]; /* 1x1 interior loops */
extern int int11_dH[NBPAIRS+1][NBPAIRS+1][5][5];

extern int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5]; /* 2x1 interior loops */
extern int int21_dH[NBPAIRS+1][NBPAIRS+1][5][5][5];

extern int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5]; /* 2x2 interior loops */
extern int int22_dH[NBPAIRS+1][NBPAIRS+1][5][5][5][5];

/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */


extern int ML_BASE37;
extern int ML_BASEdH;

extern int ML_closing37;
extern int ML_closingdH;

extern int ML_intern37;
extern int ML_interndH;

//extern int TripleC37;
//extern int TripleCdH;
//extern int MultipleCA37;
//extern int MultipleCAdH;
//extern int MultipleCB37;
//extern int MultipleCBdH;

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
extern int  MAX_NINIO;                   /* maximum correction */

extern int ninio37;
extern int niniodH;
/* penalty for helices terminated by AU (actually not GC) */

extern int TerminalAU37;
extern int TerminalAUdH;
/* penalty for forming bi-molecular duplex */
//extern int DuplexInit37;
//extern int DuplexInitdH;
/* stabilizing contribution due to special hairpins of size 4 (tetraloops) */
extern vector<string> TetraloopSeq;  /* string containing the special tetraloops */
extern int  Tetraloop37[16];  /* Bonus energy for special tetraloops */
extern int  TetraloopdH[16];
extern vector<string> TriloopSeq;    /* string containing the special triloops */
extern int  Triloop37[2]; /* Bonus energy for special Triloops */
extern int  TriloopdH[2]; /* Bonus energy for special Triloops */
extern vector<string> HexaloopSeq;    /* string containing the special triloops */
extern int  Hexaloop37[4]; /* Bonus energy for special Triloops */
extern int  HexaloopdH[4]; /* Bonus energy for special Triloops */

extern int Tetraloop[16];
extern int Triloop[2];
extern int Hexaloop[4];
extern unordered_map<string, int> hairpinE;
extern unordered_map<string, int> max_cai_map;
extern unordered_map<string, int> pair2pos;


extern double Tmeasure;       /* temperature of param measurements */

void fill_stack(const string &filename, int data, char delimeter = ',');
void fill_mismatch(const string &filename, int data, char delimeter = ',');
void fill_intl11(const string &filename, int data, char delimeter = ',');
void fill_intl21(const string &filename, int data, char delimeter = ',');
void fill_intl22(const string &filename, int data, char delimeter = ',');
void fill_codon(const string &filename, char delimeter = ',');
void scale_params(const string &file = {}, const string &paramspath = {} ,double temp = 37);

inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}


#endif //RNA_DESIGN_DEFAULT_H
