// LinearFoldCDS: left-to-right beam DP for mRNA design (Phase 3 / -m 3).
//
// This file intentionally ships NON-FUNCTIONAL. It provides the skeleton and
// recurrence outlines so the follow-up implementation has a clear starting
// point. See scripts/PHASE3_M3_DESIGN.md for the design.
//
// Expected implementation order (see design doc):
//   3a: F, C (pair closure), skip/push/pop        — ~2–3 days
//   3b: + S single-segment for bulges/internal    — ~3–5 days
//   3c: + Multi, M1, M2 with k-best heap          — ~4–7 days
//   3d: regression + tuning                       — ~2–3 days

#include "LinearFoldCDS.h"

#include <chrono>
#include <iostream>

namespace derna_lfcds {

namespace {

// Per-position chart item. Kept small for cache-friendliness during
// per-step beam sort.
struct Item {
    int a, b;           // aa indices at left / right boundaries of the item
    int8_t i, j;        // within-codon offsets {0,1,2} at boundaries
    int8_t x, y;        // codon indices at aa[a], aa[b]
    int8_t manner;      // transition kind (see enum below)
    int8_t flags;       // reserved (closed/open/hairpin markers)
    double mfe;
    double cai;
    double score;
    // Backpointer. For skeleton: a single 32-bit parent index into the
    // prior-position chart. Expand to richer struct when implementing
    // bifurcation (pair of parent indices for multi-loop combine).
    int32_t parent_idx;
};

enum class Manner : int8_t {
    Init = 0,
    Skip,
    Push,
    Pop,       // pair-close
    InternalClose,
    MultiCombine,
    MultiClose,
    External,
};

// Per-position beam. The actual implementation will prune to top-b by
// score after all extensions for a given position are produced.
using Beam = std::vector<Item>;

// k-best bifurcation helper. Given two sorted (ascending by score) beams
// A and B and an additive score function, return up to k items
// {a + b} with the smallest combined score. Standard lazy-expansion
// heap algorithm from Huang & Chiang (2005).
//
// NOTE: unimplemented. Left as a stub so the interface is visible.
template <typename Combine>
[[maybe_unused]] static void k_best_combine(const Beam& A, const Beam& B,
                                            int k,
                                            Combine combine,
                                            Beam& out) {
    (void)A; (void)B; (void)k; (void)combine; (void)out;
    // TODO: implement. Sketch:
    //   - priority_queue<tuple<double, int, int>> minheap
    //   - start with (A[0].score + B[0].score, 0, 0)
    //   - pop k times; each pop expands (ai+1, aj) and (ai, aj+1) if not seen
    //   - build out[] from popped combinations using `combine(A[ai], B[bj])`.
}

} // namespace

LinearFoldCDSResult run_linear_fold_cds(const std::vector<int>& protein,
                                         double lambda,
                                         int beam_b,
                                         int threads) {
    (void)lambda; (void)threads;
    LinearFoldCDSResult r;

    const int n = (int)protein.size();
    if (n <= 0) {
        std::cerr << "[LinearFoldCDS] empty protein; nothing to do.\n";
        return r;
    }

    const int b = (beam_b <= 0) ? 100 : beam_b;
    const int nuc_len = 3 * n;

    auto t0 = std::chrono::steady_clock::now();

    // Per-position chart. Skeleton reserves slots but never fills them.
    std::vector<Beam> chart(nuc_len + 1);
    chart[0].push_back(Item{});  // sentinel initial item

    // TODO: for pos = 0..nuc_len-1:
    //   1. For each item in chart[pos]:
    //      - push / skip / pop / multi-combine — emit zero or more
    //        successor items into chart[pos+1] (or farther if pop).
    //   2. Sort chart[pos+1] by score ascending; keep top-b.
    //   3. Maintain the (protein, codon x) dependency when position
    //      advances across a codon boundary.
    //
    // For now, skeleton just returns an unfilled result so the mode
    // wires through without producing garbage output.

    auto t1 = std::chrono::steady_clock::now();
    r.fill_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.completed = false;  // not yet implemented
    (void)nuc_len;
    (void)b;

    std::cerr << "[LinearFoldCDS] -m 8 not yet implemented; see "
                 "scripts/PHASE3_M3_DESIGN.md. Falling through to report only.\n";
    return r;
}

} // namespace derna_lfcds
