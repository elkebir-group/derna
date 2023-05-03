#include <iostream>
#include <fstream>
#include <chrono>

#include "Nussinov.h"
#include "NussinovAlgorithm.h"
#include "ZukerAlgorithm.h"
#include "Zuker.h"
#include "default.h"
#include "vector"
#include <string>
#include <tuple>
#include <regex>
#include "utils.h"
#include "params/constants.h"

using namespace std;
//namespace po = boost::program_options;
// TODO:
int main(int argc, char *argv[]) {

    int n;//len of protein
    string input,output,rna_file,swipe_output;
    string codon_file = {};
    int model, mode;
    double incr = inf, lambda = inf, threshold = 0.001;
    int g = inf;
    int const default_argc = 9;


//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/P15086.fasta", "-o", "res_zuker_P15086.txt", "-m", "1", "-s", "2", "-l", "0.06348505"};
//
//
//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/P0DTC2.fasta", "-o", "res_zuker_P0DTC2.txt", "-m", "1", "-s", "2", "-l", "0.00001"};

//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/Q14442.fasta", "-o", "res_zuker_Q14442.txt", "-m", "1", "-s", "2", "-l", "0.99999"};

//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/Q8NC38.fasta", "-o", "res_zuker_Q8NC38.txt", "-m", "1", "-s", "2", "-l", "0.00005"};
//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/Q92734.fasta", "-o", "res_zuker_Q92734.txt", "-m", "1", "-s", "2", "-l", "0.00005"};
//    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/P15421.fasta", "-o", "res_zuker_P15421.txt", "-m", "1", "-s", "2", "-l", "0.00005"};
    const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/Q8IW19.fasta", "-o", "res_zuker_Q8IW19.txt", "-m", "1", "-s", "1"};
//        const char* default_args[] = {"./RNA_Design","-i", "../data/uniprotSeq/P15421.fasta", "-o", "res_P15421_swipe.txt", "-m", "1", "-s", "3", "-O", "zuker_P15421_swipe.csv", "-l", "0.0", "-a", "0.04"};




    if (argc < 2) {
        argc = default_argc;
        argv = (char**)default_args;
    }

    try {
//        cout << argc << endl;
        size_t i = 1;
        while (i+1 < argc) {
            string param = argv[i];
//            cout << param << endl;
            if (argv[i][0] == '-') {
                switch (argv[i][1]) {
                    case 'i':
                        input = argv[i+1];
                        break;
                    case 'o':
                        output = argv[i+1];
                        break;
                    case 'm':
//                        cout << argv[i+1] << endl;
                        model = std::stoi(argv[i+1]);
                        break;
                    case 's':
//                        cout << argv[i+1] << endl;
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
                    case 't':
                        threshold = stod(argv[i+1]);
                        break;

                    case 'h':
                        help();
                    default:
                        usage();
                }
            }
            i += 2;
        }

    } catch (const std::exception& e) {
        std::cout << "Exception!" << std::endl;
        usage();
        return -1;
    }

//    cout << threshold << endl;

    bool nussinov = false;
    bool zuker = false;
    bool test = false;
    switch (model) {
        case 0:
            nussinov = true;
            break;
        case 1:
            zuker = true;
            break;
        case -1:
            test = true;
            break;
        default:
            cout << model << endl;
            throw invalid_argument("Invalid Input for Model");
    }

    if (output.empty()) throw invalid_argument("Output File Needed");
    ofstream fout(output);
    scale_params(codon_file); //"../python/pfizer_codon_usage.csv"


    if (test) {
        if (rna_file.empty()) throw invalid_argument("RNA Input File Needed in Test Mode");
        vector<int> protein = read_fasta(input, fout);

        vector<int> rna = read_rna(rna_file, fout);

        double cai = getCAI(rna, protein);
        double CAI = evaluate_CAI(rna, protein);
        double MFE = evaluate_MFE(rna, protein);

        fout << "eval MFE: " << MFE/100 << endl;
        fout << "eval CAI: " << cai << endl;
        fout << "eval standard CAI: " << CAI << endl;
        return 0;
    }

    bool mfe = false;
    bool mfe_cai = false;
    bool lambda_swipe = false;
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
        int bp = evaluate_BP_N(rna,protein,g);
        fout << "nussinov bp count: " << bp << endl;
    }

    if (nussinov && mfe_cai) {
        if (g == inf) throw invalid_argument("Invalid Value of g");
        Nussinov N = Nussinov(protein, n, g);
        tuple<double, string> temp = N.nussinov_CAI(lambda, fout);
        n_res = get<0>(temp);
        rna = get<1>(temp);
        int bp = evaluate_BP_N(rna,protein,g);
        int type = 0;
        double CAI = evaluate_CAI_N(rna,protein,type);
        fout << "lambda: " << lambda << endl;
        fout << "integrated energy: " << n_res << endl;
        fout << "CAI: " << CAI << endl;
        fout << "nussinov: " << bp << endl;
//        fout << "integrated validation: " << lambda*bp+(1-lambda)*CAI << endl;
    }

    if (nussinov && lambda_swipe) {
        Nussinov N = Nussinov(protein, n, g);
        N.lambda_swipe(incr,fout, swipe_output);
    }


    if (zuker && mfe) {
        auto start = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(input,fout,n,mode);
        Z.calculate_Z(fout);
//        Z.calculate_Z_Mem(fout);

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
        Zuker Z = Zuker(input,fout,n,mode);
        fout << "lambda: " << lambda << endl;
        double energy_cai = Z.calculate_CAI_O(fout, lambda);

        Z.traceback_B2(lambda);
//        cout << "traceback done" << endl;
        string zuker_cai_rna(3*n,'.'), zuker_cai_bp(3*n,'.');
        string zuker_cai_rna_X(3*n, '.');
        Z.get_rna_X(zuker_cai_rna_X);
        Z.get_rna_cai(zuker_cai_rna);
        Z.get_bp(zuker_cai_bp);


        int type = 0;

        double CAI_s = evaluate_CAI(zuker_cai_rna,protein,type);
//        cout << 1 << endl;
        double CAI = evaluate_CAI(zuker_cai_rna,protein,1);
//        cout << 2 << endl;
        double MFE = evaluate_MFE(zuker_cai_rna,protein);

        cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << lambda*MFE+(lambda-1)*CAI << endl;
        fout << "zuker cai bp: " << zuker_cai_bp << ",size: " << zuker_cai_bp.size() << endl;
        fout << "zuker rna: " << zuker_cai_rna_X << ".size: " << zuker_cai_rna.size() << endl;
        fout << "zuker cai rna: " << zuker_cai_rna << ".size: " << zuker_cai_rna.size() << endl;


        fout << "Codon Adaptation Index: " << CAI_s << endl;
        fout << "Minimum Free Energy: " << MFE/100 << endl;
    }

    if (zuker && lambda_swipe) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        Zuker Z = Zuker(input,fout,n,mode);
        Z.lambda_swipe_2(threshold, fout,swipe_output);
//        Z.lambda_swipe(incr,fout,swipe_output);

    }

    return 0;
}





