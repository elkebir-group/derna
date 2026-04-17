//
// Created by xinyu on 3/23/2022.
//

#include "utils.h"
#include "ZukerAlgorithm.h"
#include "NussinovAlgorithm.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

using namespace std;
int getxPos(int, const string &);

struct CodonTable {
    int aa;
    vector<int> nucleotides;
    vector<double> codon_usages;
};

typedef struct CodonTable   CodonTable;

// amino acid A = 0, R = 1, N = 2, D = 3, C = 4, Q = 5, E = 6, G = 7, H = 8, M = 9, I = 10, L = 11, K = 12, F = 13, P = 14
// S = 15, T = 16, W = 17, Y = 18, V = 19
// nucleotide A = 0, C = 1, G = 2, U = 3



// amino acid A = 0, R = 1, N = 2, D = 3, C = 4, Q = 5, E = 6, G = 7, H = 8, M = 9, I = 10, L = 11, K = 12, F = 13, P = 14
// S = 15, T = 16, W = 17, Y = 18, V = 19
int aa_index(char aa) {
    switch (aa) {
        case 'A': return 0;
        case 'R': return 1;
        case 'N': return 2;
        case 'D': return 3;
        case 'C': return 4;
        case 'Q': return 5;
        case 'E': return 6;
        case 'G': return 7;
        case 'H': return 8;
        case 'M': return 9;
        case 'I': return 10;
        case 'L': return 11;
        case 'K': return 12;
        case 'F': return 13;
        case 'P': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'W': return 17;
        case 'Y': return 18;
        case 'V': return 19;
        default:
            cout << (int)aa << " " << (char)aa << endl;
            throw invalid_argument("no char found");
    }
}

int n_index(char n) {
    switch (n) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'U': return 3;
        default: throw invalid_argument("invalid input");
    }
}

// AU, UA, GC, CG, GU, UG
// nucleotide A = 0, C = 1, G = 2, U = 3
int n_index2(int n1, int n2) {
    switch (n1-n2) {
        case -3: return 0;
        case 1:
            switch (n1) {
                case 2: return 2;
                default: return 5;
            }
        case -1:
            switch (n1) {
                case 2: return 4;
                default: return abs(n1-n2);
            }
        default:
            return abs(n1-n2);
    }
}


// nucleotide A = 0, C = 1, G = 2, U = 3
bool complementary(int X, int Y)
{
    return ((X == 0 && Y == 3) || (X == 3 && Y == 0) || (X == 1 && Y == 2) || (X == 2 && Y == 1));
}

// static constexpr int sigma(int a, int i) {
//     return 3*a+i;
// }

int m(int n, int a, int b, int i, int j, int x) { //
    switch (x) {
        case -1:
            return 9*(n-1)*(b-a)+9*(2*b-a)+(6-3*i+j);
        default:
            return 6*(9*(n-1)*(b-a)+9*(2*b-a)+(6-3*i+j))+x;
    }

}

void transform2num(vector<int> & target, string s) {
    for (int i = 0; i < (int)s.size(); i++) {
        target[i] = n_index(s[i]);
    }
}

string num2String(vector<int> & num) {
    int size = (int)num.size();
    string s(size,'.');
    for (int i = 0; i < size; i++) {
        s[i] = to_char[num[i]];
    }
    return s;
}

double getCAI(const vector<int> & rna, const vector<int> & protein) {
    int n = (int)protein.size();

    double CAI_ans = 0;
    for (int i = 0; i < n; ++i) {
        vector<int> codon(3);
        int p = protein[i];
        for (int j = 0; j <= 2; j++) {
            codon[j] = rna[3*i+j];
        }

        int x = getxPos(p, codon);

        CAI_ans += codon_cai[p][x];

    }
    return CAI_ans;
}

double stand_getCAI(const vector<int> & rna, const vector<int> & protein) {
    int n = (int)protein.size();
    double CAI_ans = 0;
    for (int i = 0; i < n; ++i) {
        vector<int> codon(3);
        int p = protein[i];
        for (int j = 0; j <= 2; j++) {
            codon[j] = rna[3*i+j];
        }
//        cout << i << endl;
        int x = getxPos(p, codon);
        CAI_ans += codon_cai[p][x];
    }
    cout << CAI_ans/n << endl;
    return exp(CAI_ans/n);
}

double getCAI_s(const vector<int> & rna, const vector<int> & protein) {
    int n = (int)protein.size();
    double CAI_ans = 0;
    for (int i = 0; i < n; ++i) {
        vector<int> codon(3);

        int p = protein[i];
        for (int j = 0; j <= 2; j++) {
            codon[j] = rna[3*i+j];
        }
//        cout << i << endl;
        int x = getxPos(p, codon);
        CAI_ans += codon_cai_s[p][x];
    }
    return CAI_ans;
}

double stand_getCAI_s(const vector<int> & rna, const vector<int> & protein) {
//    cout << "standard" << endl;
    int n = (int)protein.size();
    double CAI_ans = 0;
    for (int i = 0; i < n; ++i) {
        vector<int> codon(3);

        int p = protein[i];
        for (int j = 0; j <= 2; j++) {
            codon[j] = rna[3*i+j];
        }
//        cout << "index: " << i << ", protein: " << p << endl;
        int x = getxPos(p, codon);
        CAI_ans += codon_cai_s[p][x];
    }

    return exp(CAI_ans/n);
}

int getxPos(int p, vector<int> & codon) {
    for (int i = 0; i < 6; ++i) {
        vector<int> temp(begin(nucleotides[p][i]), end(nucleotides[p][i]));
        if (temp == codon) {
            return i;
        }
    }
    cout << p << endl;
    for (int i = 0; i < 3; ++i) {
        cout << codon[i];
    }
    cout << endl;
    throw invalid_argument("match not found");
}


// nucleotide A = 0, C = 1, G = 2, U = 3
int to_int(char a) {
    switch (a) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'U': return 3;
        default:
            cout << "char: " << char(a) << endl;
            throw invalid_argument("invalid argument to convert to number");
    }
}

void char2num(vector<int>& target, string & s) {
    for (int i = 0; i < (int)s.size(); ++i) {
        target[i] = to_int(s[i]);
    }
}


void write_csv(string filename, const vector<pair<string, vector<double>>> & dataset) {
    ofstream output(filename);
    for (int i = 0; i < (int)dataset.size(); i++) {
        output << dataset[i].first;
        if (i != (int)dataset.size() - 1) output << ",";
    }
    output << "\n";
    for (int i = 0; i < (int)dataset[0].second.size(); i++) {
        for (int j = 0; j < (int)dataset.size(); j++) {
            output << dataset[j].second[i];
            if(j != (int)dataset.size() - 1) output << ",";
        }
        output << "\n";
    }
    output.close();
}


bool compare(double x, double y, double epsilon) {
    if(fabs(x - y) < epsilon) return true;
    return false;
}

bool greaterThan(double a, double b)
{
    return a > b && !compare(a, b);
}


// nucleotide A = 0, C = 1, G = 2, U = 3
bool basepair(int X, int Y)
{
    return ((X == 0 && Y == 3) || (X == 3 && Y == 0) || (X == 1 && Y == 2) || (X == 2 && Y == 1) || (X == 2 && Y == 3) || (X == 3 && Y == 2));
}


bool add_auterminal(int a, int b) {
    return (a == 0 && b == 3) || (a == 3 && b == 0) || (a == 2 && b == 3) || (a == 3 && b == 2);
}


// nucleotide A = 0, C = 1, G = 2, U = 3
bool add_ggmm(int a, int b) {
    return (a == 2 && b == 2);
}


// nucleotide A = 0, C = 1, G = 2, U = 3
bool add_uugamm(int a, int b) {
    return (a == 3 && b == 3) || (a == 2 && b == 0);
}


int l(int a, int i, int b, int j) {
    return sigma(b,j) - sigma(a,i) - 1;
}


int get_index(vector<string> & seqs, string & seq) {
    auto index = find(seqs.begin(), seqs.end(), seq);

    if (index != seqs.end()) {
        return index - seqs.begin();
    }
    return -1;
}


void help()
{
    printf("Usage:\n");
    printf(" -i -- input file path\n");
    printf(" -o -- output file path\n");
    printf(" -m -- model <0,1,7,-1> , 0=Nussinov, 1=Zuker, 7=PositionBeamDP (fill+traceback), -1=Evaluation\n");
    printf(" -s -- mode <1,2,3>, 1 for MFE only, 2 for balancing MFE and CAI at fixed lambda, 3 for lambda sweep\n");
    printf(" -l -- lambda <[0,1]>\n");
    printf(" -a -- sweep increment <(0,1]>\n");
    printf(" -r -- input rna file path\n");
    printf(" -O -- sweep output csv file name\n");
    printf(" -g -- minimum gap allowed in Nussinov <[0,inf)>\n");
    printf(" -k -- beam width for PositionBeamDP (model 7) <[1,inf)>, default=10\n");
    printf(" -b -- beam start length <[4,inf)>, default=5 (length at which to start beam pruning)\n");
    printf(" -t -- threshold tau1 <(0,1)>\n");
    printf(" -p -- threshold tau2 <(0,1)>\n");
    printf(" -c -- codon usage table file path\n");
    printf(" -d -- directory to energy parameters\n");
    printf(" ...\n");
}

vector<int> read_rna(string & input) {
    ifstream fin(input);
    vector<int> rna;
    string line;
    char byte;
    if (!fin.is_open()) {
        cout << "Could not open the RNA file - '" << input << "'" << endl;
        exit(1);
    }
    if (fin.is_open()) {

        while (fin.get(byte)) {
            if (byte != '\n' && byte != '\r') {
                rna.push_back(to_int(byte));
            }
        }
    }
    return rna;
}


vector<int> read_fasta(string & input, ostream& fout) {
    ifstream fin(input);
    vector<int> protein;
    string line;
    char byte;
    if (!fin.is_open()) {
        fout << "Could not open the FASTA file - '" << input << "'" << endl;
        exit(1);
    }
    if (fin.is_open()) {
        getline (fin, line);
        fout << "protein sequence: ";
        while (fin.get(byte)) {
//            cout << byte << endl;
            if (byte != '\n' && byte != '\r') {
                fout << byte;
                protein.push_back(aa_index(byte));
            }
        }
        fout << endl;
    }
    return protein;
}

double evaluate_CAI(string & rna,vector<int> & protein,int type) {
    int l = int(rna.size());

//    cout << rna << endl;

    vector<int> seq(l);
    char2num(seq, rna);
    double CAI;
    if (type == 1) CAI = getCAI(seq, protein);
    else CAI = stand_getCAI_s(seq, protein);

    return CAI;

}

pair<int, int> find_amino_acid_and_codon_index(const vector<int>& codon) {
    if (codon.size() != 3) {
        throw std::invalid_argument("Codon must have exactly 3 elements");
    }

    for (int aa = 0; aa < 20; ++aa) {
        for (int k = 0; k < 6; ++k) {
            const int* triplet = nucleotides[aa][k];
            if (triplet[0] == -1) break;  // End of valid codons for this aa

            if (triplet[0] == codon[0] &&
                triplet[1] == codon[1] &&
                triplet[2] == codon[2]) {
                return {aa, k};  // Found: (amino acid index, codon index)
            }
        }
    }

    return {-1, -1};  // Not found
}

double evaluate_CAI(string & rna, int type) {
    int l = int(rna.size());
    vector<int> seq(l);
    char2num(seq, rna);
    int n = (int)(l/3);
    double CAI_ans = 0;
    int p, x;

    for (int i = 0; i < n; ++i) {
        vector<int> codon(3);

        for (int j = 0; j <= 2; j++) {
            codon[j] = seq[3*i+j];
        }
        tie(p, x) = find_amino_acid_and_codon_index(codon);
//        cout << "index: " << i << ", protein: " << p << endl;
        if (type == 0) {
            CAI_ans += codon_cai_s[p][x];
        } else {
            CAI_ans += codon_cai[p][x];
        }
    }

    if (type == 1) return CAI_ans;

    return exp(CAI_ans/n);

}


double evaluate_CAI(vector<int> & rna,vector<int> & protein) {

    double CAI = stand_getCAI_s(rna, protein);
    return CAI;

}

double evaluate_MFE(string & rna) {
    int l = int(rna.size());

//    cout << rna << endl;

    vector<int> seq(l);
    char2num(seq, rna);
    ZukerAlgorithm Zu = ZukerAlgorithm(seq,l);
    double mfe = Zu.calculate_W();
    return mfe;

}

vector<string> filterCandidates(
        const vector<string>& candidates,
        const unordered_map<int, char>& letter_map) {

    vector<string> filtered;

    for (const auto& s : candidates) {
        bool match = true;
        for (const auto& [idx, ch] : letter_map) {
            if ((size_t)idx >= s.size() || s[idx] != ch) {
                match = false;
                break;
            }
        }
        if (match) {
            filtered.push_back(s);
        }
    }

    return filtered;
}


double evaluate_MFE(vector<int> & rna, string & bp) {
    int l = int(rna.size());
    ZukerAlgorithm Zu = ZukerAlgorithm(rna,l);
    double mfe = Zu.calculate_W();
    if (!bp.empty()) {
        Zu.traceback_2();
        Zu.get_bp(bp);
    }
    return mfe;

}

// Evaluate the energy of a GIVEN (RNA, dot-bracket) pair using Turner parameters.
// Returns energy in cKcal (centikilo-calories), matching ZukerAlgorithm convention.
double evaluate_structure_energy(const vector<int>& seq, const string& bp_str) {
    int n = (int)seq.size();
    if (n == 0 || (int)bp_str.size() != n) return 0.0;

    // Build pair table: pair_table[i] = j if (i,j) paired, else -1.
    vector<int> pt(n, -1);
    vector<int> stk;
    for (int i = 0; i < n; i++) {
        if (bp_str[i] == '(') stk.push_back(i);
        else if (bp_str[i] == ')') {
            if (stk.empty()) return 1e18;  // malformed
            int j = stk.back(); stk.pop_back();
            pt[j] = i; pt[i] = j;
        }
    }
    if (!stk.empty()) return 1e18;  // malformed

    // Walk structure recursively using a stack-based approach.
    // For each base pair (i,j), find the enclosed loops and score them.
    // External loop: free energy = sum of AU penalties for stems + dangling ends (simplified).
    // We follow the Zuker convention: walk each closed pair and decompose.

    int total_energy = 0;

    // For each closing pair (i,j), find what's inside:
    // - If exactly one enclosed pair (p,q) with p=i+1,q=j-1: stacking
    // - If exactly one enclosed pair (p,q): hairpin if none, or internal/bulge loop
    // - If multiple enclosed pairs: multi-loop
    // External loop contributions (stems hanging off the exterior) are handled separately.

    // Process each closing pair (i,j) where i < j and pt[i] == j.
    for (int i = 0; i < n; i++) {
        int j = pt[i];
        if (j <= i) continue;  // only process i < j, and only where i is '('

        // Find enclosed pairs
        vector<pair<int,int>> enclosed;
        int k = i + 1;
        while (k < j) {
            if (pt[k] > k && pt[k] <= j) {
                enclosed.push_back({k, pt[k]});
                k = pt[k] + 1;
            } else {
                k++;
            }
        }

        if (enclosed.empty()) {
            // Hairpin loop: closing pair (i,j), loop length = j - i - 1
            int loop_len = j - i - 1;
            int xi = seq[i], yj = seq[j];
            int xi_ = (i + 1 < n) ? seq[i + 1] : -1;
            int _yj = (j - 1 >= 0) ? seq[j - 1] : -1;
            // Use ZukerAlgorithm's hairpin_loop formula inline
            int type = BP_pair[xi + 1][yj + 1];
            int hp_e;
            hp_e = (loop_len <= 30) ? hairpins[loop_len] : hairpins[30] + (int)(lxc * Log[loop_len]);

            if (loop_len == 3 || loop_len == 4 || loop_len == 6) {
                string s(loop_len + 2, '.');
                // Build the loop sequence string
                for (int p = 0; p < loop_len + 2; p++) {
                    int idx = i + p;
                    if (idx < n) {
                        switch (seq[idx]) {
                            case 0: s[p] = 'A'; break;
                            case 1: s[p] = 'C'; break;
                            case 2: s[p] = 'G'; break;
                            case 3: s[p] = 'U'; break;
                        }
                    }
                }
                if (loop_len == 3) {
                    if (hairpinE.count(s) > 0) { hp_e = hairpinE[s]; }
                    else { hp_e += AU[xi][yj]; }
                    total_energy += hp_e;
                    continue;
                }
                if (loop_len == 4 || loop_len == 6) {
                    if (hairpinE.count(s) > 0) { hp_e = hairpinE[s]; total_energy += hp_e; continue; }
                }
            }
            hp_e += mismatchH[type][xi_ + 1][_yj + 1];
            total_energy += hp_e;
        }
        else if (enclosed.size() == 1) {
            int p = enclosed[0].first, q = enclosed[0].second;
            int n1 = p - i - 1;  // unpaired on left
            int n2 = j - q - 1;  // unpaired on right

            if (n1 == 0 && n2 == 0) {
                // Stacking
                int type = BP_pair[seq[i] + 1][seq[j] + 1];
                int type2 = rtype[BP_pair[seq[p] + 1][seq[q] + 1]];
                total_energy += stackE[type][type2];
            }
            else if (n1 == 0 || n2 == 0) {
                // Bulge loop
                int bl = n1 + n2;  // bulge length
                int type = BP_pair[seq[i] + 1][seq[j] + 1];
                int type2 = rtype[BP_pair[seq[p] + 1][seq[q] + 1]];
                int bulge_e = (bl <= MAXLOOP) ? bulge[bl] : bulge[30] + (int)(lxc * Log[bl]);
                if (bl == 1) {
                    bulge_e += stackE[type][type2];
                } else {
                    bulge_e += AU[seq[i]][seq[j]];
                    bulge_e += AU[seq[q]][seq[p]];
                }
                total_energy += bulge_e;
            }
            else {
                // Internal loop
                int type = BP_pair[seq[i] + 1][seq[j] + 1];
                int type2 = rtype[BP_pair[seq[p] + 1][seq[q] + 1]];
                int i1 = seq[i + 1], j1 = seq[j - 1], h1 = seq[p - 1], k1 = seq[q + 1];
                int nl = max(n1, n2), ns = min(n1, n2);
                int energy;
                if (ns == 1) {
                    if (nl == 1) {
                        energy = int11[type][type2][i1 + 1][j1 + 1];
                    } else if (nl == 2) {
                        if (n1 == 1) energy = int21[type][type2][i1 + 1][k1 + 1][j1 + 1];
                        else energy = int21[type2][type][k1 + 1][i1 + 1][h1 + 1];
                    } else {
                        energy = (nl + 1 <= MAXLOOP) ? internal_loop[nl + 1] : internal_loop[30] + (int)(lxc * Log[nl + 1]);
                        energy += min(MAX_NINIO, (nl - ns) * ninio);
                        energy += mismatch1nI[type][i1 + 1][j1 + 1] + mismatch1nI[type2][k1 + 1][h1 + 1];
                    }
                } else if (ns == 2) {
                    if (nl == 2) {
                        energy = int22[type][type2][i1 + 1][h1 + 1][k1 + 1][j1 + 1];
                    } else if (nl == 3) {
                        energy = internal_loop[5] + ninio;
                        energy += mismatch23I[type][i1 + 1][j1 + 1] + mismatch23I[type2][k1 + 1][h1 + 1];
                    } else {
                        energy = (n1 + n2 <= MAXLOOP) ? internal_loop[n1 + n2] : internal_loop[30] + (int)(lxc * Log[n1 + n2]);
                        energy += min(MAX_NINIO, (nl - ns) * ninio);
                        energy += mismatchI[type][i1 + 1][j1 + 1] + mismatchI[type2][k1 + 1][h1 + 1];
                    }
                } else {
                    energy = (n1 + n2 <= MAXLOOP) ? internal_loop[n1 + n2] : internal_loop[30] + (int)(lxc * Log[n1 + n2]);
                    energy += min(MAX_NINIO, (nl - ns) * ninio);
                    energy += mismatchI[type][i1 + 1][j1 + 1] + mismatchI[type2][k1 + 1][h1 + 1];
                }
                total_energy += energy;
            }
        }
        else {
            // Multi-loop: closing pair (i,j) with multiple enclosed pairs
            int unpaired = 0;
            int prev = i + 1;
            for (auto& ep : enclosed) {
                unpaired += ep.first - prev;
                prev = ep.second + 1;
            }
            unpaired += j - prev;  // trailing unpaired before j

            int num_stems = (int)enclosed.size();
            int ml_e = ML_closing + ML_intern * (num_stems + 1) + ML_BASE * unpaired;
            // AU penalty for closing pair
            ml_e += AU[seq[i]][seq[j]];
            // AU penalty for each enclosed stem
            for (auto& ep : enclosed) {
                ml_e += AU[seq[ep.first]][seq[ep.second]];
            }
            total_energy += ml_e;
        }
    }

    // External loop: AU penalties for stems directly on the exterior
    // The Zuker algorithm adds AU penalties for external stems in W computation.
    // We add AU penalty for each stem at the external level.
    {
        int k = 0;
        while (k < n) {
            if (pt[k] > k) {
                total_energy += AU[seq[k]][seq[pt[k]]];
                k = pt[k] + 1;
            } else {
                k++;
            }
        }
    }

    return (double)total_energy;
}

double evaluate_CAI_N(string & rna,vector<int> & protein,int type) {
    int l = int(rna.size());
    vector<int> seq(l);
    transform2num(seq,rna);
    double CAI;
    if (type == 1) CAI = getCAI_s(seq, protein);
    else CAI = stand_getCAI_s(seq, protein);
    return CAI;
}

int evaluate_BP_N(string & rna, int g) {
    int l = int(rna.size());
    vector<int> seq(l);
    transform2num(seq,rna);
    NussinovAlgorithm F = NussinovAlgorithm(seq, l, g);
    int bp = F.nussinov(0, l-1);
    return bp;

}

bool is_complete_path(const Path& path) {
    return path.sector_stack.empty();
}
