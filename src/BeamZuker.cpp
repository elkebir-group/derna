//
// Created by Summer Gu on 6/4/25.
//

#include "Zuker.h"
#include "BeamZuker.h"
#include "utils.h"
#include <tuple>
#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <queue>
#include "params/constants.h"

using namespace std;

BeamZuker::BeamZuker(int n, vector<int> & protein, int k):protein(protein),n(n),k(k) {
    codon_selection.resize(n,-1);
    minX = -1, minY = -1;
    init_values();
    start_index.resize(n*n, 0);
    index_offset.resize(n*n, 0);
    int idx = 0;
    last_idx = 0;
//    cnt = 0;
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            idx = n*a+b;
            start_index[idx] = last_idx;
            index_offset[idx] = n_codon[protein[a]] * n_codon[protein[b]];
            last_idx += index_offset[idx] * 9;
        }
    }
}

//((X == 0 && Y == 3) || (X == 3 && Y == 0) || (X == 1 && Y == 2) || (X == 2 && Y == 1) || (X == 2 && Y == 3) || (X == 3 && Y == 2))
void BeamZuker::init_values() {
    basepair.resize(16,0);
    basepair[(0<<2)+3] = 1;
    basepair[(3<<2)+0] = 1;
    basepair[(1<<2)+2] = 1;
    basepair[(2<<2)+1] = 1;
    basepair[(2<<2)+3] = 1;
    basepair[(3<<2)+2] = 1;
    nucle_seq.resize(3*n,-1);
    sector.resize(3*n);
    bp_bond.resize(3*n);
    ava_nucle_p.resize(3*6*n, 0);
    ava_nucle_m.resize(3*6*n, 0);
    for (int a = 0; a < n; ++a) {
        for (int x = 0; x < n_codon[protein[a]]; ++x) {
            for (int i = 0; i < 3; ++i) {
                int idx = index(a,x,i);
                if (a < n-1) ava_nucle_p[idx] = ava_nucleotides_int(a,x,i,1);
                if (a > 0) ava_nucle_m[idx] = ava_nucleotides_int(a,x,i,-1);
            }
        }
    }
}

double BeamZuker::calculate_CAI_O(std::ostream & fout, double lambda) {
    auto start = std::chrono::high_resolution_clock::now();
    fout << "Zuker CAI" << std::endl;

    int idx;
    priority_queue<BeamEntry> prev_beam, curr_beam;
    unordered_map<int, priority_queue<BeamEntry>> pq;

    int a = 0;
    int pa = protein[a];
    const int n_codon_a = n_codon[pa];
    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
            for (int x = 0; x < n_codon_a; ++x) {
                double mfe = 0.0;
                double cai = (lambda - 1) * codon_cai[pa][x];
                double score = mfe + cai;
                BeamEntry en(score, a, a, i, j, x, x, 0, mfe, cai);
                prev_beam.push(en);
                if ((int)prev_beam.size() > k) prev_beam.pop();

                idx = index(a,a,i,j,x,x);
                O[idx] = en;
            }
        }
    }

    calculate_CAI_E(lambda);
    cout << "E done" << endl;

    int nuc_len = 3*n;
//    int la = sigma(a,i);

    for (int len = 1; len < nuc_len; ++len) {

        auto temp = prev_beam;

        for (const auto& prev : prev_beam) {
            a = prev.a;
            int b = prev.b, i = prev.i, j = prev.j;
            int x = prev.x, y = prev.y;
            idx = index(a,b,i,j,x,y);
            int pb = protein[b];

            {
                int xi = nucleotides[protein[a]][x][i];
                int yj = nucleotides[protein[b]][y][j];
                int type = BP_pair[xi+1][yj+1];

                if (type > 0) {
                    double mfe = E[idx].mfe + lambda * AU[xi][yj];
                    double cai = E[idx].cai;
                    double energy = mfe + cai;

                    BeamEntry prev_entry = O[idx];
                    BeamEntry en(energy, a, b, i, j, x, y, -1, mfe, cai,
                                 {});
                    next_pq.push(en);
                    if ((int)next_pq.size() > k) next_pq.pop();

                    if (O.count(idx) == 0 || en.score < O[idx].score) {
                        O[idx] = en;
                    }

                }

            }

            if (j < 2) {
                int new_j = j + 1;
//                int xi = nucleotides[protein[a]][x][i];
//                int yj = nucleotides[protein[b]][y][new_j];

                BeamEntry prev_entry = O[idx];
                BeamEntry en(prev_entry.score, a, b, i, new_j, x, y, -2, prev_entry.mfe, prev_entry.cai,
                             {});
                next_pq.push(en);
                if ((int)next_pq.size() > k) next_pq.pop();

                int new_idx = index(a,b,i,new_j,x,y);
                if (O.count(new_idx) == 0 || en.score < O[new_idx].score) {
                    O[new_idx] = en;
                }
            }

            if (j == 2 && b < n-1) {
                int new_b = b + 1;
                int new_j = 0;

                BeamEntry prev_entry = O[idx];
                double mfe = prev_entry.mfe;
                double cai = prev_entry.cai + (lambda - 1) * codon_cai[pb][y];
                for (int ky = 0; ky < n_codon[protein[b+1]]; ++ky) {
                    BeamEntry en(prev_entry.score, a, b, i, new_j, x, ky, -3, mfe, cai,
                                 {ky});
                    next_pq.push(en);
                    if ((int)next_pq.size() > k) next_pq.pop();

                    int new_idx = index(a,new_b,i,new_j,x,ky);
                    if (O.count(new_idx) == 0 || en.score < O[new_idx].score) {
                        O[new_idx] = en;
                    }
                }

            }

            {
                int la = sigma(a,i);
                int lb = sigma(b,j);

                for (int lc = la+1; lc <= lb-4; lc++) {
                    int c = lc / 3;
                    int i1 = lc % 3;

                    int pc = protein[c];
                    int n_codon_c = n_codon[pc];

                    for (int hx = 0; hx < n_codon_c; ++hx) {
                        if (c == a && x != hx) continue;
                        int left_idx;
                        if (i1 >= 1) {
                            left_idx = index(a, c, i, i1 - 1, x, hx);
                        } else {
                            int n_codon_cp = n_codon[protein[c-1]];
                            for (int ky = 0; ky < n_codon_cp; ++ky) {

                            }
                        }

                    }


                }

            }




        }

    }

}

void BeamZuker::calculate_CAI_E(double lambda) {
    unordered_map<int, vector<BeamEntry>> beams_by_len; // len -> BeamEntries
    unordered_map<int, BeamEntry> beam_map; // index -> best BeamEntry
    vector<BeamEntry> candidates;

    int nuc_len = 3 * n;

    // === Step 1: Initialization for len = 4 ===
    for (int a = 0; a < n - 1; ++a) {
        for (int i = 0; i < 3; ++i) {
            int b = a + 1;
            for (int j = 0; j < 3; ++j) {
                int pa = protein[a], pb = protein[b];
                for (int x = 0; x < n_codon[pa]; ++x) {
                    for (int y = 0; y < n_codon[pb]; ++y) {
                        int xi = nucleotides[pa][x][i];
                        int yj = nucleotides[pb][y][j];
                        if (BP_pair[xi+1][yj+1] == 0) continue;
                        int la = sigma(a,i), lb = sigma(b,j);
                        if (lb - la != 4) continue;

                        double mfe = hairpin_CAI(lambda, 4, a, b, pa, pb,
                                                 protein[a+1], protein[b-1],
                                                 n_codon[protein[a+1]],
                                                 n_codon[protein[b-1]],
                                                 xi, yj, i, j, x, y);
                        double cai = (lambda - 1) * (codon_cai[pa][x] + codon_cai[pb][y]);
                        double score = mfe + cai;
                        BeamEntry en(score, a, b, i, j, x, y, -1, mfe, cai);

                        int idx = index(a,b,i,j,x,y);
                        beams_by_len[4].push_back(en);
                        beam_map[idx] = en;
                        Access_E1(a,b,i,j,x,y) = mfe;
                        Access_EB(a,b,i,j,x,y) = -1;
                        E_bt[idx] = {};
                    }
                }
            }
        }
    }

    for (int len = 5; len < nuc_len; ++len) {
        cout << "\rE/M Beam -- Completed: " << (len * 100.0 / nuc_len) << "%" << flush;
        candidates.clear();

        for (auto &[prev_len, beam] : beams_by_len) {
            for (const auto &left : beam) {
                for (const auto &right : beam) {
                    if (left.b >= right.a) continue;

                    int la = sigma(left.a, left.i);
                    int lb = sigma(left.b, left.j);
                    int lc = sigma(right.a, right.i);
                    int ld = sigma(right.b, right.j);
                    if (lb < lc || ld - la != len) continue;

                    int a = left.a, b = right.b;
                    int i = left.i, j = right.j;
                    int x = left.x, y = right.y;
                    int pa = protein[a], pb = protein[b];
                    int xi = nucleotides[pa][x][i];
                    int yj = nucleotides[pb][y][j];

                    // try all 3: hairpin, internal, multi
                    double mfe1 = hairpin_CAI(lambda, len, a, b, pa, pb,
                                              protein[a+1], protein[b-1],
                                              n_codon[protein[a+1]], n_codon[protein[b-1]],
                                              xi, yj, i, j, x, y);
                    double mfe2 = internal_CAI(lambda, a, b, i, j, x, y, la, ld, xi, yj);
                    double mfe3 = multi_loop_CAI(lambda, a, b, i, j, x, y,
                                                 pa, pb, n_codon[protein[a+1]], n_codon[protein[b-1]]);
                    double mfe = min({mfe1, mfe2, mfe3});
                    double cai = (lambda - 1) * (codon_cai[pa][x] + codon_cai[pb][y]);
                    double score = mfe + cai;

                    BeamEntry en(score, a, b, i, j, x, y, -4, mfe, cai);
                    candidates.push_back(en);
                }
            }
        }

        // Top-k pruning
        sort(candidates.begin(), candidates.end());
        if ((int)candidates.size() > k) candidates.resize(k);

        for (auto &entry : candidates) {
            int idx = index(entry.a, entry.b, entry.i, entry.j, entry.x, entry.y);
            beams_by_len[len].push_back(entry);
            beam_map[idx] = entry;

            Access_E1(entry.a, entry.b, entry.i, entry.j, entry.x, entry.y) = entry.mfe;
            Access_EB(entry.a, entry.b, entry.i, entry.j, entry.x, entry.y) = entry.backtrace_type;
            E_bt[idx] = entry.bt_info;
        }
    }

    cout << "\nE/M Beam search complete with CAI" << endl;


}

void BeamZuker::calculate_CAI_M(int a, int b, int i, int j, int x, int y, double lambda) {

}