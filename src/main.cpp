#include <iostream>

#include "Nussinov.h"
#include "Zuker.h"
#include "default.h"
#include "vector"
#include <string>
#include <tuple>
#include "utils.h"
#include "params/constants.h"

using namespace std;

int main(int argc, char *argv[]) {

    int n;//len of protein
    string input = "../data/uniprotSeq/P15421.fasta";
    string output = "output.txt";
    string rna_file,swipe_output;
    string codon_file = {};
    string param_path = {};
    int model = 1, mode = 1;
    double incr = inf, lambda = inf, threshold = 0.0025, threshold2 = 0.00075;
    int g = inf;

    if (argc < 2) {
        help();
    }

    try {
        size_t i = 1;
        while ((int)i+1 <= argc) {
            string param = argv[i];
            if (argv[i][0] == '-') {
                switch (argv[i][1]) {
                    case 'i':
                        input = argv[i+1];
                        break;
                    case 'o':
                        output = argv[i+1];
                        break;
                    case 'm':
                        model = std::stoi(argv[i+1]);
                        break;
                    case 's':
                        mode = std::stoi(argv[i+1]);
                        break;
                    case 'g':
                        g = std::stoi(argv[i+1]);
                        break;
                    case 'l':
                        lambda = std::stod(argv[i+1]);
                        break;
                    case 'a':
                        incr = std::stod(argv[i+1]);
                        break;
                    case 'r':
                        rna_file = argv[i+1];
                        break;
                    case 'O':
                        swipe_output = argv[i+1];
                        break;
                    case 'c':
                        codon_file = argv[i+1];
                        break;
                    case 'd':
                        param_path = argv[i+1];
                        break;
                    case 't':
                        threshold = stod(argv[i+1]);
                        break;
                    case 'p':
                        threshold2 = stod(argv[i+1]);
                        break;
                    default:
                        help();
                        return(0);
                }
            }
            i += 2;
        }

    } catch (const std::exception& e) {
        std::cout << "Exception!" << std::endl;
        help();
        return -1;
    }

    bool nussinov = false, zuker = false, test = false;
    bool subopt = false;
    switch (model) {
        case 0:
            nussinov = true;
            break;
        case 1:
            zuker = true;
            break;
        case 2:
            zuker = true;
            subopt = true;
            break;
        case -1:
            test = true;
            break;
        default:
            throw invalid_argument("Invalid Input for Model");
    }

    if (output.empty()) throw invalid_argument("Output File Needed");
    ofstream fout(output);
    scale_params(codon_file, param_path); //"../python/pfizer_codon_usage.csv"


    if (test) {
        if (rna_file.empty()) throw invalid_argument("RNA Input File Needed in Test Mode");
        if (input.empty()) throw invalid_argument("Protein Input File Needed in Test Mode");
        vector<int> protein = read_fasta(input, fout);

        vector<int> rna = read_rna(rna_file);
        string bp(rna.size(), '.');

        double cai = getCAI(rna, protein);
        double CAI = evaluate_CAI(rna, protein);
        double MFE = evaluate_MFE(rna, bp);

        fout << "secondary structure: " << bp << endl;
        fout << "eval MFE: " << MFE/100 << endl;
        fout << "eval CAI: " << cai << endl;
        fout << "eval standard CAI: " << CAI << endl;
        return 0;
    }

    bool mfe = false;
    bool mfe_cai = false;
    bool lambda_swipe = false;
    bool lambda_swipe2 = false;
    switch (mode) {
        case 1:
            mfe = true;
            break;
        case 2:
            mfe_cai = true;
            break;
        case 3:
            lambda_swipe = true;
            break;
        case 4:
            lambda_swipe2 = true;
            break;
        default:
            throw invalid_argument("Invalid Input for Mode");
    }



    if (input.empty()) throw invalid_argument("Input File Needed");

    vector<int> protein = read_fasta(input, fout);
    n = int(protein.size());
    fout << endl;
    double n_res = 0;
    string rna;


    if (nussinov && mfe) {
        if (g == inf) throw invalid_argument("Invalid Value of g");
        Nussinov N = Nussinov(protein, n, g);
        tuple<double, string, string> temp = N.nussinov(fout);
        n_res = get<0>(temp);
        rna = get<1>(temp);
        int bp = evaluate_BP_N(rna,g);
        fout << "nussinov bp count: " << bp << endl;
    }

    if (nussinov && mfe_cai) {
        if (g == inf) throw invalid_argument("Invalid Value of g");
        Nussinov N = Nussinov(protein, n, g);
        tuple<double, string> temp = N.nussinov_CAI(lambda, fout);
        n_res = get<0>(temp);
        rna = get<1>(temp);
        int bp = evaluate_BP_N(rna,g);
        int type = 0;
        double CAI = evaluate_CAI_N(rna,protein,type);
        fout << "lambda: " << lambda << endl;
        fout << "integrated energy: " << n_res << endl;
        fout << "CAI: " << CAI << endl;
        fout << "nussinov: " << bp << endl;
    }

    if (nussinov && lambda_swipe) {
        Nussinov N = Nussinov(protein, n, g);
        N.lambda_swipe(incr,fout, swipe_output);
    }


    if (zuker && mfe) {
        auto start = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        Z.calculate_Z(fout);
        Z.traceback_B();
        auto end = chrono::high_resolution_clock::now();
        long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
        fout << "Time taken : " << time_take;
        fout << "sec" << endl;
        string zuker_bp(3*n,'.'), zuker_rna(3*n,'.'), zuker_rna_X(3*n,'.');// zuker_bp2(3*n,'.'),zuker_bp1(3*n,'.');
        vector<string> rna_array(n, zuker_rna);
        Z.get_bp(zuker_bp);
        Z.get_rna(zuker_rna);
        Z.get_rna_X(zuker_rna_X);
        double cai = evaluate_CAI(zuker_rna, protein, 0);


        fout << "zuker bp:" << zuker_bp << ", size: " << zuker_bp.size() << endl;
        fout << "zuker rna:" << zuker_rna_X << ", size: " << zuker_rna.size() << endl;
        fout << "zuker rna:" << zuker_rna << ", size: " << zuker_rna.size() << endl;
        fout << "zuker cai: " << cai << endl;
        fout << "other rna: " << endl;
    }

    if (zuker && mfe_cai) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        Zuker Z = Zuker(n,mode,protein);
        fout << "lambda: " << lambda << endl;
        double energy_cai = Z.calculate_CAI_O(fout, lambda);

        Z.traceback_B2(lambda);
        string zuker_cai_rna(3*n,'.'), zuker_cai_bp(3*n,'.');
        string zuker_cai_rna_X(3*n, '.');
        Z.get_rna_X(zuker_cai_rna_X);
        Z.get_rna_cai(zuker_cai_rna);
        Z.get_bp(zuker_cai_bp);
        Z.save_all_vectors("./");

        int type = 0;

        double CAI_s = evaluate_CAI(zuker_cai_rna,protein,type);
        double CAI = evaluate_CAI(zuker_cai_rna,protein,1);
        double MFE = evaluate_MFE(zuker_cai_rna);

        cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << lambda*MFE+(lambda-1)*CAI << endl;
        fout << "zuker cai bp: " << zuker_cai_bp << ",size: " << zuker_cai_bp.size() << endl;
        fout << "zuker rna: " << zuker_cai_rna_X << ".size: " << zuker_cai_rna.size() << endl;
        fout << "zuker cai rna: " << zuker_cai_rna << ".size: " << zuker_cai_rna.size() << endl;

        fout << "Codon Adaptation Index: " << CAI_s << endl;
        fout << "Free Energy: " << MFE/100 << endl;

        if (subopt) {
            mt19937 rng(42);
            double i = 0.01;
            while (i < 1000) {
                cout << "temperature: " << i << endl;

                // temperature 100 missing the last codon
                // codon does not match amino acid
                // softmax_sample: options list is empty?
                //
                Z.traceback_suboptimal(lambda, i, rng);
                i *= 10;
                string subopt_cai_rna(3*n,'.'), subopt_cai_bp(3*n,'.');
                string subopt_cai_rna_X(3*n, '.');
//                cout << "traceback done" << endl;
                Z.get_rna_X(zuker_cai_rna_X);
//                cout << "getting rna X done" << endl;
                Z.get_rna_cai(zuker_cai_rna);
//                cout << "getting rna done" << endl;
                Z.get_bp(zuker_cai_bp);
//                cout << "getting bp done" << endl;

//                cout << zuker_cai_rna_X << endl;
//                cout << zuker_cai_rna << endl;
//                cout << zuker_cai_bp << endl;

                CAI_s = evaluate_CAI(zuker_cai_rna,protein,type);
                CAI = evaluate_CAI(zuker_cai_rna,protein,1);
//                CAI_s = evaluate_CAI(zuker_cai_rna);
//                CAI = evaluate_CAI(zuker_cai_rna,1);
                MFE = evaluate_MFE(zuker_cai_rna);

                cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << lambda*MFE+(lambda-1)*CAI << endl;
                fout << "zuker cai bp: " << zuker_cai_bp << ",size: " << zuker_cai_bp.size() << endl;
                fout << "zuker rna: " << zuker_cai_rna_X << ".size: " << zuker_cai_rna.size() << endl;
                fout << "zuker cai rna: " << zuker_cai_rna << ".size: " << zuker_cai_rna.size() << endl;

                fout << "Codon Adaptation Index: " << CAI_s << endl;
                fout << "Free Energy: " << MFE/100 << endl;
            }

        }
    }

    if (zuker && lambda_swipe) {
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_swipe_2(threshold,threshold2, fout,swipe_output);
    }

    if (zuker && lambda_swipe2) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        if (incr == inf) throw invalid_argument("Invalid increment");
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_swipe(incr,fout,swipe_output);
    }

    return 0;
}





