// LinearFoldCDS: left-to-right beam DP for mRNA design (Phase 3 / -m 3).
// Same (a, b, i, j, x, y) codon-pair state representation as -m 2 and -m 7,
// but iterates 5'→3' with per-step beam pruning and k-best bifurcation.
//
// Not yet functional — this is a scoping skeleton. See
// scripts/PHASE3_M3_DESIGN.md for the design and effort estimate.
#ifndef DERNA_LINEARFOLD_CDS_H
#define DERNA_LINEARFOLD_CDS_H

#include <cstddef>
#include <string>
#include <vector>

namespace derna_lfcds {

struct LinearFoldCDSResult {
    double beam_score = 0.0;     // lambda*mfe + (lambda-1)*cai of the best path
    double mfe = 0.0;            // in kcal/mol * 100 (raw Zuker units)
    double cai_raw = 0.0;        // sum of log(w(c)) across all codons
    double cai_geom = 0.0;       // exp(cai_raw / n)
    std::string rna;             // final mRNA (length = 3 * protein.size())
    std::string structure;       // dot-bracket (same length as rna)
    double fill_ms = 0.0;
    double traceback_ms = 0.0;
    bool completed = false;      // false if fallback or not-implemented
};

// Main entry point. protein[] = 0..19 amino-acid indices (as in default.cpp).
// beam_b = per-position beam width; 0 → LinearFold default (100).
// threads = OpenMP thread count; 0 = library default.
LinearFoldCDSResult run_linear_fold_cds(const std::vector<int>& protein,
                                         double lambda,
                                         int beam_b,
                                         int threads);

} // namespace derna_lfcds

#endif // DERNA_LINEARFOLD_CDS_H
