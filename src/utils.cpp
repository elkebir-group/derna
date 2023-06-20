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

int sigma(int a, int i) {
    return 3*a+i;
}

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

//        cout << i << " " << x << endl;

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
//        CAI_ans += log(codon_cai[p][x]);
//        cout << codon_cai[p][x] << " " << log(codon_cai[p][x]) << " " << CAI_ans << endl;
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
//        cout << i << endl;
        int x = getxPos(p, codon);
//        CAI_ans *= codon_cai_s[p][x];
        CAI_ans += codon_cai_s[p][x];
//        cout << i << " " << codon_cai_s[p][x] << endl;
//        CAI_ans += log(codon_cai_s[p][x]);
    }
//    cout << CAI_ans << " " << n << endl;
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
            cout << char(a) << endl;
            throw invalid_argument("invalid argument to convert to number");
    }
}

void char2num(vector<int>& target, string & s) {
    for (int i = 0; i < (int)s.size(); ++i) {
//        cout << i << endl;
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
    printf(" -m -- model <0,1,-1> , 0 for Nussinov-based model, 1 for Zuker-based model, -1 for Evaluation\n");
    printf(" -s -- mode <1,2,3>, 1 for MFE only, 2 for balancing MFE and CAI at fixed lambda, 3 for lambda sweep\n");
    printf(" -l -- lambda <[0,1]>\n");
    printf(" -a -- sweep increment <(0,1]>\n");
    printf(" -r -- input rna file path\n");
    printf(" -O -- sweep output csv file name\n");
    printf(" -g -- minimum gap allowed in Nussinov <[0,inf)>\n");
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

double evaluate_MFE(vector<int> & rna, string & bp) {
    int l = int(rna.size());
    ZukerAlgorithm Zu = ZukerAlgorithm(rna,l);
    double mfe = Zu.calculate_W();
    if (!bp.empty()) {
        Zu.traceback();
        Zu.get_bp(bp);
    }
    return mfe;

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