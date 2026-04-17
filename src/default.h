//
// Created by Summer Gu on 8/16/22.
//

#ifndef RNA_DESIGN_DEFAULT_H
#define RNA_DESIGN_DEFAULT_H

#include <unordered_map>
#include <fstream>
#include <queue>
#include <initializer_list>
#include <cstdint>
#include <vector>

using namespace std;

#include "params/constants.h"
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)
#define PUBLIC

// Fixed-capacity inline backtrace buffer used in place of std::vector<int>
// for BeamEntry/XYVariant. Trivially copyable (memcpy-able), no heap allocs.
// Max sizes observed across call sites:
//   PositionBeamDP.cpp S_C_StoC: 9 ints (worst fixed-layout case)
//   PositionBeamDP.cpp SpecialHP: up to ~4 ints (codon choices for hairpin span)
//   Zuker.cpp bifurcation: 7 ints
//   CS-table tie append (update_cs): potentially unbounded but never read by traceback;
//       capped silently here, which is correct (only one tied entry is needed).
// CAP=16 gives comfortable headroom and nice alignment (16*4=64 bytes).
struct BtInfo {
    static constexpr int CAP = 16;
    int8_t len = 0;
    int data[CAP] = {};

    int size() const { return (int)len; }
    bool empty() const { return len == 0; }
    void clear() { len = 0; }

    int  operator[](int k) const { return data[k]; }
    int& operator[](int k)       { return data[k]; }
    int  operator[](size_t k) const { return data[k]; }
    int& operator[](size_t k)       { return data[k]; }

    void push_back(int v) { if (len < CAP) data[len++] = v; }

    void assign(std::initializer_list<int> il) {
        len = 0;
        for (int v : il) { if (len < CAP) data[len++] = v; }
    }
    BtInfo& operator=(std::initializer_list<int> il) { assign(il); return *this; }

    // Interop with legacy std::vector<int> call sites (Zuker.cpp, PositionBasedBeamZuker.cpp).
    BtInfo& operator=(const std::vector<int>& v) {
        len = 0;
        for (int x : v) { if (len < CAP) data[len++] = x; }
        return *this;
    }
    BtInfo(const std::vector<int>& v) {
        for (int x : v) { if (len < CAP) data[len++] = x; }
    }
    BtInfo(std::initializer_list<int> il) {
        for (int x : il) { if (len < CAP) data[len++] = x; }
    }
    BtInfo() = default;

    const int* begin() const { return data; }
    const int* end()   const { return data + len; }
    int* begin() { return data; }
    int* end()   { return data + len; }

    operator std::vector<int>() const { return std::vector<int>(data, data + len); }
};

// One codon-pair variant stored inside a merged structural BeamEntry.
// The BeamEntry's top-level (x,y,mfe,cai,score,bt_info,...) always reflect the BEST variant;
// the variants vector holds all (x,y) possibilities for extension.
struct XYVariant {
    int x = -1, y = -1;
    double mfe = 0.0, cai = 0.0, score = inf;
    int backtrace_type = 0;
    int last_closed_nuc = -1;
    BtInfo bt_info;
    // Per-variant CS / structural metadata. Every predecessor-identifying field
    // must live here, because a merged entry may hold variants whose source
    // transitions differ: the top-level mirror only reflects the BEST variant.
    int cs_inner_key = -1;
    int cs_right_s_key = -1;
    int cs_right_len = -1;
    int cs_single_start = -1;
    int cs_inner_left = -1;
    int cs_inner_right = -1;
    long long cs_pack_outer = -1LL;
    XYVariant() = default;
    XYVariant(int x_, int y_, double mfe_, double cai_, double sc_,
              int bt_type_, int last_closed_, const BtInfo& bt_)
        : x(x_), y(y_), mfe(mfe_), cai(cai_), score(sc_),
          backtrace_type(bt_type_), last_closed_nuc(last_closed_), bt_info(bt_) {}
};

struct BeamEntry {
    // --- 8-BYTE MEMBERS ---
    double score = inf; // Assuming inf is defined in your scope
    double mfe = 0.0;
    double cai = 0.0;

    // --- HEAVY MEMBERS (Pointers / Size_t) ---
    // Inline fixed-size backtrace (trivially copyable, no heap).
    BtInfo bt_info;

    // --- 4-BYTE MEMBERS (ints) ---
    int a = -1;
    int b = -1;
    int x = -1;
    int y = -1;

    int backtrace_type = 0;
    int last_closed_nuc = -1;

    int cs_inner_key = -1;
    int cs_right_s_key = -1;
    int cs_right_len = -1;
    int cs_single_start = -1;
    int cs_inner_left = -1;
    int cs_inner_right = -1;
    // cs_pack_outer: int CS key from cs_get_index (stored as long long; -1 = not from CS)
    long long cs_pack_outer = -1LL;
    long long cs_pack_inner_c = -1; // unused (kept for ABI compat)

    // --- 1-BYTE MEMBERS ---
    // i and j are guaranteed to be <= 2, so 8 bits is plenty.
    int8_t i = -1;
    int8_t j = -1;

    // double score;
    // int a, b, i, j, x, y;
    // int backtrace_type;
    // double mfe, cai;
    // vector<int> bt_info;  // backtrace keys only; do not use for last_closed_nuc
    // int last_closed_nuc = -1;  // dedicated: last closing position (for multi-loop SINGLE_MAX_LEN)
    // // Explicit CS-state payload so CS entries do not rely on geometry reconstruction.
    // int cs_inner_key = -1;      // predecessor C key used to build this CS state
    // int cs_right_s_key = -1;    // predecessor right-S key used to build this CS state
    // int cs_right_len = -1;      // length of the right single segment in this CS state
    // int cs_single_start = -1;   // starting nucleotide position of the right single segment
    // int cs_inner_left = -1;     // inner C left boundary position p
    // int cs_inner_right = -1;    // inner C right boundary position q

//    BeamEntry() : score(inf), a(-1), b(-1), i(-1), j(-1), x(-1), y(-1),
//                  backtrace_type(0), mfe(0.0), cai(0.0), bt_info() {}

    BeamEntry(double s = inf, int a_ = -1, int b_ = -1, int i_ = -1, int j_ = -1, int x_ = -1, int y_ = -1,
              int bt_type = 0, double mfe_ = 0, double cai_ = 0, const BtInfo& bt = BtInfo{}, int last_closed = -1)
            : score(s), mfe(mfe_), cai(cai_), bt_info(bt),
              a(a_), b(b_), x(x_), y(y_),
              backtrace_type(bt_type), last_closed_nuc(last_closed),
              i(i_), j(j_) {}

    bool operator<(const BeamEntry &other) const {
        return score < other.score;
    }

    // All (x,y) codon-pair variants for this structural state, for use during fill expansion.
    // The top-level x,y,mfe,cai,score,bt_info always mirror the best variant.
    std::vector<XYVariant> variants;
};

struct LambdaResult {
    double mfe_value;
    double cai_value;
    double CAI_value;
    double O_val;
};

typedef struct st {
    int i; // index i
    int j; // index j
    int ml; //?
} st;

// IDEALLY we would reorganize stack_ as follows:
// struct stack_ {
//     // --- 8-BYTE MEMBERS ---
//     double change = 0.0;
//
//     // --- 4-BYTE MEMBERS ---
//     int a = -1;
//     int b = -1;
//
//     // --- 1-BYTE MEMBERS ---
//     // Assuming i and j are the small 0-2 offsets from BeamEntry
//     int8_t i = -1;
//     int8_t j = -1;
//
//     // Codons
//     int8_t x = -1;
//     int8_t y = -1;
//
//     // Assuming ml is a multi-loop boolean flag or small state integer
//     bool ml = false;
//
//     // Total size: 24 bytes (down from 40 bytes)
//     // 8 (double) + 8 (two ints) + 5 (1-byte types) + 3 bytes padding = 24
// };

typedef struct stack_ {
    int a;
    int b;
    int i; // index i
    int j; // index j
    int x; // nucleotide at i
    int y; // nucleotide at j
    int ml; //?
    double change;
//    vector<int> bi;
} stack_;

typedef struct bond {
    int i;
    int j;
} bond;

struct Path {
    vector<stack_> sector_stack;
    vector<int> nucle_seq;
    vector<int> codon_selection;
    vector<bond> bp_bond;
    double change = 0;
    int t, s;


    Path(int n) {
        nucle_seq.resize(3*n, -1);
        codon_selection.resize(n, -1);
        bp_bond.resize(3*n);
//        sector_stack.resize(3*n);
        bp_bond[0].i = 0;
        t = 0, s = 0;
    }

    Path() = default;
};

struct PathHash {
    size_t operator()(const Path& p) const {
        size_t h = 0;

        // Hash nucle_seq
        for (int x : p.nucle_seq) {
            h ^= std::hash<int>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }

        // Hash change (rounded to avoid floating point noise)
        long long rounded = static_cast<long long>(p.change * 1e6);  // adjust precision as needed
        h ^= std::hash<long long>{}(rounded) + 0x9e3779b9 + (h << 6) + (h >> 2);

        return h;
    }
};

struct PathEqual {
    bool operator()(const Path& a, const Path& b) const {
        return a.nucle_seq == b.nucle_seq && std::abs(a.change - b.change) < 1e-6;
    }
};

struct OptionKey {
    int t;
    double e1, e2;

    bool operator<(const OptionKey& other) const {
        if (t != other.t) return t < other.t;
        if (std::abs(e1 - other.e1) > EPSILON) return e1 < other.e1;
        if (std::abs(e2 - other.e2) > EPSILON) return e2 < other.e2;
        return false;
    }
};

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

/** Per-pair-type multi-loop stem penalty (LCDSfold v_score_M1 / v_score_multi). tt = pair type 1..NBPAIRS. */
extern int E_MLstem(int tt);

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
void fill_miscellaneous(const string &filename, char delimeter = ',');
void scale_params(const string &file = {}, const string &paramspath = {} ,double temp = 37);

inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}


#endif //RNA_DESIGN_DEFAULT_H
