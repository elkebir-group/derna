#include <iostream>

#include "Nussinov.h"
#include "Zuker.h"
#include "PositionBeamDP.h"
#include "PositionBasedBeamZuker.h"
#include "LinearFoldCDS.h"
#include "default.h"
#include <string>
#include <tuple>
#include <vector>
#include <chrono>
#include "utils.h"
#include "params/constants.h"
#include <fstream>
#include <cstdlib>
#include <unistd.h>  // For _exit() system call
#ifdef _OPENMP
#include <omp.h>
#endif

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
    int g = inf, k = 10;  // k is beam width for PositionBeamDP (model 7)
    [[maybe_unused]] int beam_start_len = 5;  // Length to start beam pruning (default: 5)
    int num_threads = 0;  // 0 = let OpenMP choose (respects OMP_NUM_THREADS); set via -j

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
                    case 'k':
                        k = std::stoi(argv[i+1]);
                        break;
                    case 'b':
                        beam_start_len = std::stoi(argv[i+1]);
                        break;
                    case 'j':
                        num_threads = std::stoi(argv[i+1]);
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

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        std::cerr << "[threads] using " << num_threads << " OpenMP threads (-j)" << std::endl;
    }
#else
    if (num_threads > 0) {
        std::cerr << "[threads] -j " << num_threads
                  << " requested but derna was built without OpenMP; running single-threaded."
                  << std::endl;
    }
#endif

    // Tie-break fix: at strict lambda=1 the DP objective is 1*MFE + 0*CAI, so
    // among MFE-optima the codon pair is chosen by arbitrary (insertion-order)
    // tie-breaks that systematically land on low-frequency synonyms. Nudge
    // lambda by a tiny epsilon to make CAI the secondary objective while
    // keeping the MFE optimum exactly (eps * max|CAI sum| << MFE granularity).
    // Applies uniformly to -m 2 (Zuker) and -m 7 (PositionBeamDP).
    constexpr double CAI_TIEBREAK_EPS = 1e-6;
    if (lambda == 1.0) {
        std::cerr << "[tiebreak] lambda=1.0 -> using " << (1.0 - CAI_TIEBREAK_EPS)
                  << " internally so MFE ties are broken toward higher CAI." << std::endl;
        lambda = 1.0 - CAI_TIEBREAK_EPS;
    }

    bool nussinov = false, zuker = false, test = false, position_beam_dp = false, pos_based_beam_zuker = false;
    bool linear_fold_cds = false;
    bool subopt = false, subopt_all = false;
    switch (model) {
        case 0:
            nussinov = true;
            break;
        case 1:
            zuker = true;
            break;
        case 2:
            zuker = true;
            break;
        case 3:
            zuker = true;
            break;
        case 4:
            throw invalid_argument("Model 4 (BeamZuker) was removed; use model 7 (PositionBeamDP) for beam search");
        case 6:
            pos_based_beam_zuker = true;
            break;
        case 7:
            position_beam_dp = true;
            break;
        case 8:
            // LinearFoldCDS: left-to-right beam DP (Phase 3).
            // SKELETON ONLY — see scripts/PHASE3_M3_DESIGN.md. Recurrences not
            // yet implemented; exits after printing a not-implemented notice.
            linear_fold_cds = true;
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

    if (linear_fold_cds) {
        if (input.empty()) throw invalid_argument("Input File Needed");
        vector<int> protein = read_fasta(input, fout);
        int b_beam = (k > 0) ? k : 100;
        derna_lfcds::LinearFoldCDSResult res = derna_lfcds::run_linear_fold_cds(
            protein, lambda, b_beam, num_threads);
        fout << "LinearFoldCDS (-m 8) skeleton invoked. See "
                "scripts/PHASE3_M3_DESIGN.md for status." << endl;
        fout << "completed: " << (res.completed ? "yes" : "no") << endl;
        fout << "beam_b: " << b_beam << ", lambda: " << lambda << endl;
        return 0;
    }

    bool mfe = false;
    bool mfe_cai = false;
    bool lambda_sweep = false;
    bool lambda_sweep2 = false;
    switch (mode) {
        case 1:
            mfe = true;
            break;
        case 2:
            mfe_cai = true;
            break;
        case 3:
            lambda_sweep = true;
            break;
        case 4:
            lambda_sweep2 = true;
            break;
        default:
            throw invalid_argument("Invalid Input for Mode");
    }
    // Model 7 (PositionBeamDP) and model 6 (PositionBasedBeamZuker) are MFE+CAI
    if (position_beam_dp || pos_based_beam_zuker)
        mfe_cai = true;



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

    if (nussinov && lambda_sweep) {
        Nussinov N = Nussinov(protein, n, g);
        N.lambda_sweep(incr,fout, swipe_output);
    }


    if (zuker && mfe) {
        auto start = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        Z.calculate_Z(fout);
        Z.traceback_B();
        auto end = chrono::high_resolution_clock::now();
        double time_s = chrono::duration<double>(end - start).count();
        fout << "Time taken : " << time_s << " sec" << endl;
        cout << "Zuker (mfe) completed: time=" << time_s << " s" << endl;
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
        auto start_z = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        fout << "lambda: " << lambda << endl;
        double energy_cai = Z.calculate_CAI_O(fout, lambda);

        Z.traceback_B2(lambda);
        double time_z = chrono::duration<double>(chrono::high_resolution_clock::now() - start_z).count();
        cout << "Zuker (mfe_cai) completed: time=" << time_z << " s, O=" << energy_cai << endl;
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
        
        // VERIFICATION: Compare forward vs recomputed
        double recomputed_total = lambda * MFE + (lambda - 1) * CAI;
        double forward_total = energy_cai;
        double total_diff = abs(recomputed_total - forward_total);
        
        fout << "\n=== ZUKER BACKTRACE VERIFICATION ===" << endl;
        fout << "Forward Total: " << forward_total << endl;
        fout << "Recomputed Total: " << recomputed_total << endl;
        fout << "Difference: " << total_diff << endl;
        if (total_diff > 1.0) {
            fout << "WARNING: Zuker backtrace verification FAILED!" << endl;
        } else {
            fout << "Zuker backtrace verification PASSED." << endl;
        }
        fout << "=====================================\n" << endl;

        cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << lambda*MFE+(lambda-1)*CAI << endl;
        fout << "zuker cai bp: " << zuker_cai_bp << ",size: " << zuker_cai_bp.size() << endl;
        fout << "zuker rna: " << zuker_cai_rna_X << ".size: " << zuker_cai_rna.size() << endl;
        fout << "zuker cai rna: " << zuker_cai_rna << ".size: " << zuker_cai_rna.size() << endl;

        fout << "Codon Adaptation Index: " << CAI_s << endl;
        fout << "Free Energy: " << MFE/100 << endl;
// #ifdef DEBUG_ZUKER_LOGGING
//         // Dump O table for comparison (only when DEBUG_ZUKER_LOGGING is defined)
//         Z.dump_O_table(fout);
// #endif

        if (subopt) {
            mt19937 rng(60);
            double i = 1;

            // Open CSV file:
            ofstream csv_fout("zuker_subopt_sample.csv");
            csv_fout << "gamma,sequence,bp,MFE,CAI,combined,optimal,count\n";


            while (i >= 0.99) {
                cout << "gamma: " << i << endl;

                Z.traceback_suboptimal(lambda, i, rng);

                string subopt_cai_rna(3*n,'.'), subopt_cai_bp(3*n,'.');
                string subopt_cai_rna_X(3*n, '.');

                Z.get_rna_X(subopt_cai_rna_X);
                Z.get_rna_cai(subopt_cai_rna);
                Z.get_bp(subopt_cai_bp);

                CAI_s = evaluate_CAI(subopt_cai_rna,protein,type);
                CAI = evaluate_CAI(subopt_cai_rna,protein,1);
                MFE = evaluate_MFE(subopt_cai_rna);

                double combined = lambda*MFE + (lambda-1)*CAI;
                int is_optimal = (combined == energy_cai ? 1 : 0); // Always 0 for subopt paths in this block

//                vector<Path> all_paths;
//                Z.traceback_enumerate_dfs(lambda, i, 20000, all_paths);
                size_t count = Z.traceback_count_dfs(lambda, i);

                // Console output
                cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI
                     << ",cai_s: " << CAI_s << ",mfe: " << MFE
                     << ",combined: " << combined << ",count: " << count << endl;

                fout << "zuker cai bp: " << subopt_cai_bp << ",size: " << subopt_cai_bp.size() << endl;
                fout << "zuker rna: " << subopt_cai_rna_X << ".size: " << subopt_cai_rna.size() << endl;
                fout << "zuker cai rna: " << subopt_cai_rna << ".size: " << subopt_cai_rna.size() << endl;

                fout << "Codon Adaptation Index: " << CAI_s << endl;
                fout << "Free Energy: " << MFE/100 << endl;
                fout << "---------" << endl;

                // CSV output
                csv_fout << i << "," // gamma
                         << subopt_cai_rna << ","
                         << subopt_cai_bp << ","
                         << MFE/100 << ","
                         << CAI_s << ","
                         << combined << ","
                         << is_optimal << ","
                         << count << "\n";

                i -= 0.001;
            }



            csv_fout.close();
        }

        double gamma = 0.99;
        if (subopt_all) {
            vector<Path> all_paths;
            Z.traceback_enumerate_dfs(lambda, gamma, inf, all_paths); // e.g. 1000 paths

            cout << "Total paths: " << all_paths.size() << endl;

            ofstream csv_fout("zuker_subopt_paths_0.99.csv");
            csv_fout << "sequence,bp,MFE,CAI,combined,optimal,path\n";

            csv_fout << zuker_cai_rna << ","
                     << zuker_cai_bp << ","
                     << MFE/100 << ","      // as you output it in your console
                     << CAI_s << ","        // you said use CAI_s here
                     << energy_cai << ","
                     << 1 << ","
                     << energy_cai << "\n";

            for (size_t p = 0; p < all_paths.size(); ++p) {
                const Path& path = all_paths[p];

                // Assign path data into Z internal state:
                Z.load_path(path);

                // Prepare output strings:
                string subopt_cai_rna(3*n, '.');
                string subopt_cai_bp(3*n, '.');
                string subopt_cai_rna_X(3*n, '.');

                // Call your usual functions:
                Z.get_rna_X(subopt_cai_rna_X);
                Z.get_rna_cai(subopt_cai_rna);
                Z.get_bp(subopt_cai_bp);

                CAI_s = evaluate_CAI(subopt_cai_rna, protein, type);
                CAI   = evaluate_CAI(subopt_cai_rna, protein, 1);
                MFE   = evaluate_MFE(subopt_cai_rna);

                double combined = lambda*MFE + (lambda-1)*CAI;
                int is_optimal = (combined == energy_cai ? 1 : 0);

                // Console output
                cout << "Path " << p+1 << " / " << all_paths.size() << endl;
                cout << "RNA_X : " << subopt_cai_rna_X << endl;
                cout << "RNA   : " << subopt_cai_rna << endl;
                cout << "BP    : " << subopt_cai_bp << endl;
                cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << combined << ",path sum: " << path.change << endl;
                cout << "---------" << endl;

                // CSV output
                csv_fout << subopt_cai_rna << ","
                         << subopt_cai_bp << ","
                         << MFE/100 << ","      // as you output it in your console
                         << CAI_s << ","        // you said use CAI_s here
                         << combined << ","
                         << is_optimal << ","
                         << path.change << "\n";
            }
            csv_fout.close();
        }
    }

    if (zuker && lambda_sweep) {
        auto start_sw = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_sweep_2(threshold,threshold2, fout,swipe_output);
        cout << "Zuker (lambda_sweep_2) completed: time=" << chrono::duration<double>(chrono::high_resolution_clock::now() - start_sw).count() << " s" << endl;
    }

    if (zuker && lambda_sweep2) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        if (incr == inf) throw invalid_argument("Invalid increment");
        auto start_sw2 = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_sweep(incr,fout,swipe_output);
        cout << "Zuker (lambda_sweep) completed: time=" << chrono::duration<double>(chrono::high_resolution_clock::now() - start_sw2).count() << " s" << endl;
    }

    // PositionBasedBeamZuker (model 6): reference LCDSfold-style beam search
    if (pos_based_beam_zuker) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        auto start_pbbz = chrono::high_resolution_clock::now();
        PositionBasedBeamZuker pbbz(n, protein, k);
        double best_score_pbbz = pbbz.calculate_position_based(lambda);
        long pbbz_ms = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start_pbbz).count();
        cout << "PositionBasedBeamZuker completed: best score=" << best_score_pbbz << ", time=" << pbbz_ms << " ms" << endl;
        pbbz.cleanup_before_destruction();
    }

    // PositionBeamDP (model 7): fill_position_beam_tables + traceback_position_beam_tables
    if (position_beam_dp && mfe_cai) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        fout << "\n=== PositionBeamDP (fill + traceback) ===" << endl;
        fout << "Lambda: " << lambda << ", Beam width k: " << k << endl;
        AllTablesDerna tables;
        auto start_pdp = chrono::high_resolution_clock::now();
        fill_position_beam_tables(n, protein, lambda, k, tables);
        auto end_fill = chrono::high_resolution_clock::now();

        // Instrumentation: dump C[left,right] entries for a given (left,right) pair.
        // Env: DERNA_DUMP_C="left_lo,left_hi,right_lo,right_hi" — dump all C entries
        // where left_pos in [left_lo,left_hi] AND right_pos in [right_lo,right_hi].
        if (const char* env_dump = std::getenv("DERNA_DUMP_C")) {
            int ll, lh, rl, rh;
            if (sscanf(env_dump, "%d,%d,%d,%d", &ll, &lh, &rl, &rh) == 4) {
                int nuc_len = 3 * n;
                static const char to_ch[4] = {'A','C','G','U'};
                cerr << "[DERNA_DUMP_C] dumping C entries for left in [" << ll << "," << lh
                     << "] right in [" << rl << "," << rh << "]\n";
                for (int rp = rl; rp <= rh && rp < nuc_len; ++rp) {
                    if (rp < 0 || rp >= (int)tables.bestC.size()) continue;
                    for (auto& kv : tables.bestC[rp]) {
                        int key = kv.first;
                        int nuc_ro = key & 3;
                        int nuc_lo = (key >> 2) & 3;
                        int right_pos = (key >> 4) % nuc_len;
                        int left_pos = key / (nuc_len * 16);
                        if (left_pos < ll || left_pos > lh) continue;
                        const auto& e = kv.second;
                        cerr << "  C[" << left_pos << "," << right_pos << "]"
                             << " (nL=" << to_ch[nuc_lo] << ",nR=" << to_ch[nuc_ro] << ")"
                             << " score=" << e.score << " mfe=" << e.mfe << " cai=" << e.cai
                             << " manner=" << e.backtrace_type
                             << " nVar=" << e.variants.size() << "\n";
                        for (size_t vi = 0; vi < e.variants.size() && vi < 16; ++vi) {
                            const auto& v = e.variants[vi];
                            cerr << "    var[" << vi << "] x=" << v.x << " y=" << v.y
                                 << " score=" << v.score << " mfe=" << v.mfe << " cai=" << v.cai
                                 << " manner=" << v.backtrace_type << "\n";
                        }
                    }
                }
                cerr << "[DERNA_DUMP_C] done\n";
            }
        }

        // DERNA_SCAN_M1_LEFT="left" — scan all M1 right positions for entries with given left.
        if (const char* env_sm = std::getenv("DERNA_SCAN_M1_LEFT")) {
            int target_left;
            if (sscanf(env_sm, "%d", &target_left) == 1) {
                int nuc_len = 3*n;
                static const char to_ch[4] = {'A','C','G','U'};
                cerr << "[DERNA_SCAN_M1_LEFT] target_left=" << target_left << "\n";
                int seen = 0;
                for (int rp = 0; rp < (int)tables.bestM1.size() && rp < nuc_len; ++rp) {
                    for (auto& kv : tables.bestM1[rp]) {
                        int key = kv.first;
                        int left_pos = key >> 4;
                        if (left_pos != target_left) continue;
                        int nuc_R = key & 3, nuc_L = (key >> 2) & 3;
                        const auto& e = kv.second;
                        cerr << "  M1[" << left_pos << "," << rp << "] (nL=" << to_ch[nuc_L]
                             << ",nR=" << to_ch[nuc_R] << ") mfe=" << e.mfe
                             << " manner=" << e.backtrace_type
                             << " bt_size=" << e.bt_info.size() << "\n";
                        if (++seen >= 50) { cerr << "  ...truncated\n"; break; }
                    }
                    if (seen >= 50) break;
                }
                cerr << "[DERNA_SCAN_M1_LEFT] done (" << seen << ")\n";
            }
        }

        // DERNA_DUMP_MULTI="left_lo,left_hi,right_lo,right_hi" — dump Multi/M1/M2 entries.
        if (const char* env_m = std::getenv("DERNA_DUMP_MULTI")) {
            int ll, lh, rl, rh;
            if (sscanf(env_m, "%d,%d,%d,%d", &ll, &lh, &rl, &rh) == 4) {
                int nuc_len = 3 * n;
                static const char to_ch[4] = {'A','C','G','U'};
                auto dump_table = [&](const char* label, const std::vector<DernaBeamMap>& tab) {
                    cerr << "[DERNA_DUMP_" << label << "] left in [" << ll << "," << lh
                         << "] right in [" << rl << "," << rh << "]\n";
                    for (int rp = rl; rp <= rh && rp < (int)tab.size(); ++rp) {
                        for (auto& kv : tab[rp]) {
                            int key = kv.first;
                            int nuc_R = key & 3;
                            int nuc_L = (key >> 2) & 3;
                            int left_pos = key >> 4;
                            if (left_pos < ll || left_pos > lh) continue;
                            const auto& e = kv.second;
                            cerr << "  " << label << "[" << left_pos << "," << rp << "]"
                                 << " (nL=" << to_ch[nuc_L] << ",nR=" << to_ch[nuc_R] << ")"
                                 << " score=" << e.score << " mfe=" << e.mfe << " cai=" << e.cai
                                 << " manner=" << e.backtrace_type
                                 << " nVar=" << e.variants.size() << "\n";
                        }
                    }
                    (void)nuc_len;
                };
                dump_table("MULTI", tables.bestMulti);
                dump_table("M1",    tables.bestM1);
                dump_table("M2",    tables.bestM2);
                cerr << "[DERNA_DUMP_MULTI] table sizes:";
                for (int rp = rl; rp <= rh && rp < (int)tables.bestM1.size(); ++rp) {
                    cerr << " M1[" << rp << "]=" << tables.bestM1[rp].size()
                         << " M2[" << rp << "]=" << tables.bestM2[rp].size()
                         << " Multi[" << rp << "]=" << tables.bestMulti[rp].size();
                }
                cerr << "\n[DERNA_DUMP_MULTI] done\n";
            }
        }

        // DERNA_DUMP_F="pos1,pos2" — dump tab_f[pos1] and tab_f[pos2] entries.
        if (const char* env_f = std::getenv("DERNA_DUMP_F")) {
            int p1=-1,p2=-1;
            int n_parsed = sscanf(env_f, "%d,%d", &p1, &p2);
            (void)n_parsed;
            int nuc_len = 3*n;
            static const char to_ch[4] = {'A','C','G','U'};
            for (int pp : {p1, p2}) {
                if (pp < 0 || pp >= (int)tables.bestF.size()) continue;
                cerr << "[DERNA_DUMP_F] tab_f[" << pp << "] (" << tables.bestF[pp].size() << " entries)\n";
                int shown = 0;
                for (auto& kv : tables.bestF[pp]) {
                    int key = kv.first;
                    int nuc_R = key & 3;
                    int nuc_L = (key >> 2) & 3;
                    int left_pos = key >> 4;
                    const auto& e = kv.second;
                    cerr << "  F left=" << left_pos << " (nL=" << to_ch[nuc_L] << ",nR=" << to_ch[nuc_R] << ")"
                         << " score=" << e.score << " mfe=" << e.mfe << " cai=" << e.cai
                         << " manner=" << e.backtrace_type << " a=" << e.a << " b=" << e.b
                         << " i=" << (int)e.i << " j=" << (int)e.j << " nVar=" << e.variants.size() << "\n";
                    if (++shown >= 20) { cerr << "  ... (" << (tables.bestF[pp].size()-20) << " more)\n"; break; }
                }
                (void)nuc_len;
            }
            cerr << "[DERNA_DUMP_F] done\n";
        }

        vector<int> nucle_seq, codon_selection;
        vector<bond> bp_bond;
        traceback_position_beam_tables(tables, n, protein, nucle_seq, codon_selection, bp_bond);
        auto end_pdp = chrono::high_resolution_clock::now();
        long fill_ms = chrono::duration_cast<chrono::milliseconds>(end_fill - start_pdp).count();
        long trace_ms = chrono::duration_cast<chrono::milliseconds>(end_pdp - end_fill).count();

        // Fill any residue not visited during traceback with default codon 0 so the full sequence is defined
        for (int i = 0; i < n; ++i) {
            if (codon_selection[i] < 0 && protein[i] >= 0 && protein[i] < 20 && n_codon[protein[i]] > 0)
                codon_selection[i] = 0;
        }

        // Reconcile nucle_seq from codon_selection so each codon is consistent (avoids mixed single-base writes causing invalid codons in getxPos)
        for (int i = 0; i < n; ++i) {
            int px = codon_selection[i];
            if (px >= 0 && protein[i] >= 0 && protein[i] < 20 && px < n_codon[protein[i]]) {
                for (int j = 0; j < 3; ++j)
                    nucle_seq[3 * i + j] = nucleotides[protein[i]][px][j];
            }
        }

        int final_pos = 3 * n - 1;
        double best_score = 1e30;  // minimize (match Zuker)
        double best_dp_mfe = 0.0, best_dp_cai = 0.0;
        for (const auto& item : tables.bestF[final_pos]) {
            const auto& e = item.second;
            int left_bound = sigma(e.a, e.i), right_bound = sigma(e.b, e.j);
            if (left_bound == 0 && right_bound == final_pos && item.second.score < best_score) {
                best_score = item.second.score;
                best_dp_mfe = e.mfe;
                best_dp_cai = e.cai;
            }
        }
        if (best_score > 1e29) best_score = -1e30;  // no full-span F
        fout << "Fill time: " << fill_ms << " ms, Traceback time: " << trace_ms << " ms" << endl;

        string rna_str(3 * n, '.');
        string bp_str(3 * n, '.');
        for (int p = 0; p < 3 * n; ++p) {
            int aa_idx = p / 3;
            if (aa_idx < n && codon_selection[aa_idx] >= 0 && p < (int)nucle_seq.size() && nucle_seq[p] >= 0 && nucle_seq[p] <= 3)
                rna_str[p] = to_char[nucle_seq[p]];
        }
        for (const auto& b : bp_bond) {
            if (b.i >= 0 && b.i < 3 * n && b.j >= 0 && b.j < 3 * n) {
                bp_str[b.i] = '(';
                bp_str[b.j] = ')';
            }
        }

        double eval_CAI = (rna_str.find('.') == string::npos) ? evaluate_CAI(rna_str, protein, 0) : 0.0;
        double eval_CAI_raw = (rna_str.find('.') == string::npos) ? evaluate_CAI(rna_str, protein, 1) : 0.0;

        // Full Zuker refold: optimal MFE for this RNA (may use a different structure)
        string zuker_bp(3*n, '.');
        vector<int> rna_vec(3*n);
        for (int p = 0; p < 3*n; p++) rna_vec[p] = nucle_seq[p];
        double eval_MFE = evaluate_MFE(rna_vec, zuker_bp);

        // Evaluate the energy of the TRACEBACK structure specifically (not a refold)
        double tb_struct_energy = evaluate_structure_energy(rna_vec, bp_str);

        // Count structure differences between traceback and Zuker optimal
        int bp_diffs = 0;
        for (int p = 0; p < 3*n; p++) {
            if (bp_str[p] != zuker_bp[p]) bp_diffs++;
        }

        fout << "RNA (from traceback): " << rna_str << endl;
        fout << "Structure: " << bp_str << endl;
        fout << "Evaluated CAI (geom mean): " << eval_CAI << ", Evaluated CAI (raw sum): " << eval_CAI_raw << endl;
        fout << "Traceback structure MFE: " << tb_struct_energy << ", Zuker optimal MFE: " << eval_MFE
             << ", structure diffs: " << bp_diffs << endl;
        double recomputed_combined = lambda * tb_struct_energy + (lambda - 1) * eval_CAI_raw;
        fout << "Recomputed combined (traceback structure): " << recomputed_combined << endl;
        fout << "DP beam score: " << best_score << ", mfe: " << best_dp_mfe << ", cai: " << best_dp_cai << endl;
        if (abs(best_dp_mfe - tb_struct_energy) > 1.0)
            fout << "WARNING: DP accumulated MFE (" << best_dp_mfe << ") != traceback structure energy (" << tb_struct_energy << ")" << endl;

        cout << "PositionBeamDP completed: best score=" << best_score << ", fill=" << fill_ms << " ms, traceback=" << trace_ms << " ms" << endl;
    }

    return 0;
}





