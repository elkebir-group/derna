//
// Position-based beam DP filling AllTablesDerna (BeamEntry and Zuker energy functions).
// If adding debug logging, guard with #ifdef DEBUG_POSITION_BEAM_DP_LOGGING (same strategy as Zuker.cpp).
//

#include "PositionBeamDP.h"
#include "Zuker.h"
#include "utils.h"
#include "default.h"
#include "codon_agg.h"
#include "params/constants.h"

// USE_CODON_LOOKUP_CLOSURE: compile-time flag to enable Step 3 closure rewrite.
// 0 (default): existing Block 1 K_S × K_C × K_S iteration.
// 1         : new closure that iterates 36-cell codon-pair agg grids (~5 populated
//             cells per grid in practice). Both paths produce structurally-identical
//             curr_c entries (differences allowed only via tied alternatives at equal
//             score).
#ifndef USE_CODON_LOOKUP_CLOSURE
#define USE_CODON_LOOKUP_CLOSURE 0
#endif
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace std::chrono;

// LCDSfold v_score_external_paired: terminal AU penalty only (type > 2 => GU/UG/AU/UA).
// PositionBeam nucleotides 0..3 (A,C,G,U); BP_pair is 1-indexed.
static inline int v_score_external_paired_pb(int nuci0, int nucj0) {
    const int type = BP_pair[nuci0 + 1][nucj0 + 1];  // 1..6 for valid pairs
    return (type > 2) ? TerminalAU37 : 0;
}

// LCDSfold v_score_M1 / v_score_multi equivalents for PositionBeam nucleotides 0..3.
// Return **positive** Turner/Vienna penalties (0.1 kcal/mol units), consistent with constants.h tables.
// LCDSfold E_MLstem adds terminal-AU for type>2; if this codebase's E_MLstem does not, we add it here
// so C->M1 and M1+C->M2 match LCDSfold (fixes pos=5 A-U M1 score gap of 25 = 50 MFE units).
static inline int v_score_M1_pb(int nuci0, int nucj0) {
    const int type = BP_pair[nuci0 + 1][nucj0 + 1];
    if (type == 0) return (int)inf;
    return E_MLstem(type) + ((type > 2) ? TerminalAU37 : 0);
}

static inline int v_score_multi_pb(int nuci0, int nucj0) {
    const int tt = BP_pair[nucj0 + 1][nuci0 + 1];  // reversed
    if (tt == 0) return (int)inf;
    // Match LCDSfold: E_MLstem(tt) + ML_closing37.  Since this codebase's E_MLstem
    // returns only ML_intern (no TerminalAU), add it explicitly for AU/GU closing pairs.
    return E_MLstem(tt) + ((tt > 2) ? TerminalAU37 : 0) + ML_closing37;
}

static inline bool same_or_next_codon(int from_b, int from_j, int bb, int slot, bool& same_codon) {
    same_codon = (from_b == bb && slot == from_j + 1);
    const bool next_codon = (from_b == bb - 1 && from_j == 2 && slot == 0);
    return same_codon || next_codon;
}

// Combined score: lambda*mfe + (lambda-1)*cai. mfe and cai are stored unweighted in BeamEntry.
static inline double combined_score(double lambda, double mfe, double cai) {
    return lambda * mfe + (lambda - 1) * cai;
}

// update_derna: merged-state variant update.
//
// Key is index_struct(a,b,i,j,n) — structural position only, no (x,y).
// Each entry carries a vector<XYVariant> with one slot per distinct (x,y) codon pair.
// The entry's top-level fields (score, x, y, mfe, cai, bt_info, ...) always mirror the best variant.
//
// This collapses all codon variants of the same structural state into one beam slot,
// reducing the C/N/Multi beam sizes from ~8400 to ~234 entries per position — matching
// LCDSfold's effective beam cardinality while retaining exact CAI tracking per variant.
static inline void fill_variant_cs_fields(XYVariant& v, const BeamEntry& en) {
    v.cs_inner_key    = en.cs_inner_key;
    v.cs_right_s_key  = en.cs_right_s_key;
    v.cs_right_len    = en.cs_right_len;
    v.cs_single_start = en.cs_single_start;
    v.cs_inner_left   = en.cs_inner_left;
    v.cs_inner_right  = en.cs_inner_right;
    v.cs_pack_outer   = en.cs_pack_outer;
}

// Canonical total order on XYVariant: (score, x, y, backtrace_type, bt_info, mfe, cai).
// Independent of insertion order so downstream per-variant iteration is deterministic
// regardless of map iteration order.
static inline bool xyv_lt(const XYVariant& a, const XYVariant& b) {
    if (a.score != b.score) return a.score < b.score;
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.backtrace_type != b.backtrace_type) return a.backtrace_type < b.backtrace_type;
    if (a.bt_info.len != b.bt_info.len) return a.bt_info.len < b.bt_info.len;
    for (int k = 0; k < a.bt_info.len; ++k)
        if (a.bt_info.data[k] != b.bt_info.data[k]) return a.bt_info.data[k] < b.bt_info.data[k];
    if (a.mfe != b.mfe) return a.mfe < b.mfe;
    return a.cai < b.cai;
}

// Insertion sort for small arrays. variants vector is capped at 36 and typical
// size is 4-12, where insertion sort beats std::sort (introsort) by 2-3x due to
// lower overhead. Stable and cache-friendly.
template <typename It, typename Comp>
static inline void insertion_sort_variants(It first, It last, Comp lt) {
    for (It i = first + 1; i < last; ++i) {
        auto key = std::move(*i);
        It j = i;
        while (j > first && lt(key, *(j - 1))) {
            *j = std::move(*(j - 1));
            --j;
        }
        *j = std::move(key);
    }
}

// Seed a freshly-inserted map entry with its first (x, y) variant. Extracted
// so batch callers that already hold an iterator can skip the hash lookup.
static inline void _update_derna_seed(BeamEntry& entry, double score, const BeamEntry& en) {
    for (auto& row : entry.variant_idx) for (auto& v2 : row) v2 = -1;
    entry.variants.clear();
    entry.variants.emplace_back(en.x, en.y, en.mfe, en.cai, score,
                                en.backtrace_type, en.last_closed_nuc, en.bt_info);
    fill_variant_cs_fields(entry.variants.back(), en);
    if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6)
        entry.variant_idx[en.x][en.y] = 0;
}

// Body of update_derna for an existing map entry (no try_emplace). Extracted so
// batch callers can reuse after hoisting the hash lookup outside an inner loop.
static void _update_derna_existing(BeamEntry& ex, double score, const BeamEntry& en);

// Returns true iff a new slot was created in the map (caller can append idx to a
// task_order log without a second find()).
static bool update_derna(DernaBeamMap& m, int idx, double score, const BeamEntry& en) {
    auto [it, inserted] = m.try_emplace(idx, en);
    if (inserted) {
        _update_derna_seed(it->second, score, en);
        return true;
    }
    _update_derna_existing(it->second, score, en);
    return false;
}

// OPT2: variant of update_derna for hot-path callers (Block 1 emit_rb / emit_il).
// Avoids the implicit copy of `en` at try_emplace and the redundant variant_idx
// loop reset in _update_derna_seed (the just-default-constructed entry already
// has variant_idx = -1 from BeamEntry's class default initializer). On insert we
// only write the fields the emit lambdas actually set; everything else stays at
// the default ctor's -1 / -1LL, which is what the prior path also produced.
//
// IMPORTANT: do not use this from callers whose `en` carries non-default
// last_closed_nuc / cs_inner_key / cs_right_s_key / cs_pack_outer — emit_rb /
// emit_il leave those at -1 / -1LL, which matches the default ctor.
static bool update_derna_fresh_scratch(DernaBeamMap& m, int idx, double score, const BeamEntry& en) {
    auto [it, inserted] = m.try_emplace(idx);
    BeamEntry& ex = it->second;
    if (inserted) {
        // Default-constructed in place: variant_idx is already all -1. Write
        // only the fields that emit_rb / emit_il populate.
        ex.score = score;
        ex.mfe = en.mfe; ex.cai = en.cai;
        ex.a = en.a; ex.b = en.b;
        ex.i = en.i; ex.j = en.j;
        ex.x = en.x; ex.y = en.y;
        ex.backtrace_type = en.backtrace_type;
        ex.bt_info = en.bt_info;
        ex.cs_inner_left   = en.cs_inner_left;
        ex.cs_right_len    = en.cs_right_len;
        ex.cs_single_start = en.cs_single_start;
        ex.cs_inner_right  = en.cs_inner_right;
        // Seed variant 0; skip the redundant variant_idx-loop reset.
        ex.variants.emplace_back(en.x, en.y, en.mfe, en.cai, score,
                                 en.backtrace_type, en.last_closed_nuc, en.bt_info);
        fill_variant_cs_fields(ex.variants.back(), en);
        if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6)
            ex.variant_idx[en.x][en.y] = 0;
        return true;
    }
    _update_derna_existing(ex, score, en);
    return false;
}

// VARDIST diagnostic counters (Phase 0a). Only updated when DERNA_VARDIST_DIAG=1.
// linear_find_calls = enters the (x<0 || x>=6) fallback branch in _update_derna_existing.
// fast_path_calls   = uses variant_idx[x][y] O(1) lookup.
struct VarDistDiag {
    uint64_t fast_path_calls = 0;
    uint64_t linear_find_calls = 0;
};
static thread_local VarDistDiag g_vardist_diag;
static std::mutex g_vardist_diag_mu;
static VarDistDiag g_vardist_diag_total;
static bool vardist_diag_enabled() {
    static const bool v = []() {
        const char* e = std::getenv("DERNA_VARDIST_DIAG");
        return e && e[0] == '1';
    }();
    return v;
}

// VARDIST: TSV output stream for per-pos variants distribution stats.
// Opened lazily on first write; closed at end of fill_position_beam_tables.
static std::ofstream s_vardist_tsv;
static void vardist_open_tsv_if_needed() {
    if (!vardist_diag_enabled()) return;
    if (s_vardist_tsv.is_open()) return;
    s_vardist_tsv.open("/tmp/derna_variants_dist.tsv", std::ios::out);
    if (s_vardist_tsv.is_open()) {
        s_vardist_tsv << "pos\ttable\ttotal_states\ttotal_variants\tavg_variants\tmax_variants"
                         "\tavg_buckets_used\tavg_variants_per_bucket\n";
    }
}

// Measure variants[] distribution across all entries of a DernaBeamMap, treating
// each variant's (nuc_lo, nuc_ro) bucket as an integer 0..15 derived from the
// variant's own (x,y) and the entry's structural protein/codon-slot indices.
//
// Reports (pos, table, total_states, total_variants, avg_variants, max_variants,
// avg_buckets_used, avg_variants_per_bucket).
//
// avg_buckets_used: for each state, count distinct (nuc_lo, nuc_ro) bins among its
// variants, then average across populated states. <= 1 means partitioning by
// boundary nucs would NOT split any state — i.e., the state is already
// effectively partitioned.
//
// avg_variants_per_bucket: under the proposed boundary-nuc partitioning, average
// variants[] size per resulting sub-state.
static void vardist_measure_map(int pos, const char* table_name,
                                const DernaBeamMap& m,
                                const std::vector<int>& protein) {
    if (!vardist_diag_enabled()) return;
    if (m.empty()) return;
    vardist_open_tsv_if_needed();
    if (!s_vardist_tsv.is_open()) return;

    uint64_t total_states = 0;
    uint64_t total_variants = 0;
    uint64_t sum_buckets_used = 0;
    uint64_t sum_subbuckets = 0;  // sum over (state, bucket) of variants_in_bucket
    uint64_t sum_subbucket_count = 0;  // # of nonempty (state, bucket) pairs
    uint64_t max_variants = 0;

    for (const auto& kv : m) {
        const BeamEntry& e = kv.second;
        if (e.variants.empty()) continue;
        total_states++;
        total_variants += e.variants.size();
        if (e.variants.size() > max_variants) max_variants = e.variants.size();

        // Bucket variants by (nuc_lo, nuc_ro).
        uint8_t bucket_count[16] = {0};
        for (const auto& v : e.variants) {
            // Skip uninitialized/special (-1) variants defensively.
            if (v.x < 0 || v.x >= 6 || v.y < 0 || v.y >= 6) continue;
            // protein[a] and protein[b] indices, with i,j slots.
            int paa = protein[e.a];
            int pbb = protein[e.b];
            int nlo = nucleotides[paa][v.x][e.i];
            int nro = nucleotides[pbb][v.y][e.j];
            int b = nlo * 4 + nro;
            if (b >= 0 && b < 16) bucket_count[b]++;
        }
        int buckets_used = 0;
        for (int b = 0; b < 16; ++b) {
            if (bucket_count[b] > 0) {
                buckets_used++;
                sum_subbuckets += bucket_count[b];
                sum_subbucket_count++;
            }
        }
        sum_buckets_used += buckets_used;
    }

    if (total_states == 0) return;
    double avg_variants = (double)total_variants / (double)total_states;
    double avg_buckets_used = (double)sum_buckets_used / (double)total_states;
    double avg_variants_per_bucket = (sum_subbucket_count > 0)
        ? (double)sum_subbuckets / (double)sum_subbucket_count : 0.0;

    s_vardist_tsv << pos << "\t" << table_name
                  << "\t" << total_states
                  << "\t" << total_variants
                  << "\t" << avg_variants
                  << "\t" << max_variants
                  << "\t" << avg_buckets_used
                  << "\t" << avg_variants_per_bucket << "\n";
}

// Forward decl: defined after DenseBeamTable struct below.
struct DenseBeamTable;
static void vardist_measure_dense(int pos, const char* table_name,
                                  const DenseBeamTable& td,
                                  const std::vector<int>& key_order,
                                  const std::vector<int>& protein);

// Implementation of the existing-entry variant update (post-try_emplace logic).
static void _update_derna_existing(BeamEntry& ex, double score, const BeamEntry& en) {
    // Helper to propagate ALL top-level fields from `en` when `en` becomes the new best.
    // CRITICAL: cs_* metadata fields must be copied along with bt_info, otherwise
    // the new bt_info points to predecessors described by the OLD cs_* fields,
    // causing traceback to look in wrong tables.
    auto promote_to_best = [&]() {
        ex.score = score; ex.x = en.x; ex.y = en.y;
        ex.mfe = en.mfe; ex.cai = en.cai;
        ex.backtrace_type = en.backtrace_type;
        ex.bt_info = en.bt_info;
        ex.last_closed_nuc = en.last_closed_nuc;
        ex.cs_inner_left   = en.cs_inner_left;
        ex.cs_right_len    = en.cs_right_len;
        ex.cs_single_start = en.cs_single_start;
        ex.cs_inner_right  = en.cs_inner_right;
        ex.cs_inner_key    = en.cs_inner_key;
        ex.cs_right_s_key  = en.cs_right_s_key;
        ex.cs_pack_outer   = en.cs_pack_outer;
    };
    // Lex compare BtInfo: (len, data...) for deterministic tie-break.
    auto bt_lt = [](const BtInfo& a, const BtInfo& b) {
        if (a.len != b.len) return a.len < b.len;
        for (int k = 0; k < a.len; ++k) if (a.data[k] != b.data[k]) return a.data[k] < b.data[k];
        return false;
    };
    // Direct (x,y) -> variants[] index lookup: O(1).
    int vi = -1;
    if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6) {
        vi = ex.variant_idx[en.x][en.y];
        if (vardist_diag_enabled()) g_vardist_diag.fast_path_calls++;
    } else {
        if (vardist_diag_enabled()) g_vardist_diag.linear_find_calls++;
        for (size_t k = 0; k < ex.variants.size(); ++k) {
            if (ex.variants[k].x == en.x && ex.variants[k].y == en.y) { vi = (int)k; break; }
        }
    }
    if (vi >= 0) {
        XYVariant& v = ex.variants[vi];
        bool replace = false;
        if (score < v.score) replace = true;
        else if (score == v.score) {
            if (en.backtrace_type < v.backtrace_type) replace = true;
            else if (en.backtrace_type == v.backtrace_type && bt_lt(en.bt_info, v.bt_info)) replace = true;
        }
        if (replace) {
            v = XYVariant(en.x, en.y, en.mfe, en.cai, score,
                          en.backtrace_type, en.last_closed_nuc, en.bt_info);
            fill_variant_cs_fields(v, en);
            if (score < ex.score ||
                (score == ex.score &&
                 (en.backtrace_type < ex.backtrace_type ||
                  (en.backtrace_type == ex.backtrace_type && bt_lt(en.bt_info, ex.bt_info)))))
                promote_to_best();
            // OPT-VARIANT-INCR: in-place replacement keeps variants[vi] at the same
            // index, so variant_idx[x][y]=vi is already correct. The previous
            // insertion_sort_variants + variant_idx rebuild was for deterministic
            // iteration order across runs, not correctness — variants are read
            // index-agnostically (via variant_idx O(1) lookups or for-each),
            // and VARIANT_CAP=36=6×6 means truncation never actually fires.
        }
        return;
    }
    // New (x,y) for this structural state. Append + record its index incrementally;
    // no sort, no full variant_idx rebuild. Sizes can never exceed 36 (=6×6) since
    // variants are uniquely keyed by (x,y).
    int8_t new_idx = (int8_t)ex.variants.size();
    ex.variants.emplace_back(en.x, en.y, en.mfe, en.cai, score,
                             en.backtrace_type, en.last_closed_nuc, en.bt_info);
    fill_variant_cs_fields(ex.variants.back(), en);
    if (score < ex.score) promote_to_best();
    if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6)
        ex.variant_idx[en.x][en.y] = new_idx;
}

// Batch-friendly variant: upsert map entry and apply update, given the caller has
// already placed a template_entry. If the caller's iterator is provided via `hint_it`,
// reuse that to avoid re-hashing.
static inline DernaBeamMap::iterator update_derna_batch_upsert(
    DernaBeamMap& m, int idx, double score, const BeamEntry& en,
    std::vector<int>& order) {
    auto [it, inserted] = m.try_emplace(idx, en);
    if (inserted) {
        _update_derna_seed(it->second, score, en);
        order.push_back(idx);
    } else {
        _update_derna_existing(it->second, score, en);
    }
    return it;
}

// After the first call with a key, subsequent updates with the SAME key can
// skip the hash lookup by using the iterator returned above.
static inline void update_derna_batch_same_key(
    DernaBeamMap::iterator it, double score, const BeamEntry& en) {
    _update_derna_existing(it->second, score, en);
}

// BLK1_AGG_DIAG instrumentation: globals for measuring per-(idx,x,y) variant
// slot redundancy on Block 1's emit path. Default OFF (zero overhead unless
// the env var is set). When BLK1_AGG_DIAG=1, _update_derna_existing /
// DenseBeamTable::update increment these counters; end-of-run logs report
// avg/max redundancy and a wall-clock estimate of per-update cost.
struct Blk1AggDiag {
    uint64_t seed_calls = 0;            // new compact_idx (first hit)
    uint64_t new_variant_calls = 0;     // existing idx, new (x,y)
    uint64_t existing_variant_calls = 0;// existing idx AND existing (x,y) — redundant
    uint64_t total_update_ns = 0;       // wall ns spent inside update()
    uint64_t per_pos_max_redundancy = 0;// max hits-on-single-(idx,x,y)-slot seen
};
static thread_local Blk1AggDiag g_blk1_agg_diag;
static std::vector<Blk1AggDiag> g_blk1_agg_diag_threads;
static std::mutex g_blk1_agg_diag_mu;
static bool blk1_agg_diag_enabled() {
    static const bool v = []() {
        const char* e = std::getenv("BLK1_AGG_DIAG");
        return e && e[0] == '1';
    }();
    return v;
}

// Per-pos thread-local hit-count vector for the (compact_idx, x, y) slot, used
// only when BLK1_AGG_DIAG=1. Sized to dense_cap * 36 and reset at pos start.
// The (x,y) range is [0,6) each so per slot we use 36 sub-slots.
static thread_local std::vector<uint8_t> g_blk1_agg_hits_per_xy;

// DenseBeamTable: dense-vector replacement for DernaBeamMap on Block 1's
// per-thread task_map hot path. The Block 1 emit lambdas key candidates by
// nuc_key_c_nuc(left_pos, pos, nuc_lo, nuc_ro), which decomposes into a
// per-pos compact index (left_pos*16 + nuc_lo*4 + nuc_ro). Per-pos peak
// per-thread density is ~1000 entries (avg 233) and the dense address space
// is at most nuc_len*16 (≈ 5712 for 1000-aa proteins) — much smaller than
// the raw key span (~2M), so a dense vector with sentinel-score reset is
// far cheaper than the unordered_dense::map's hash + insertion path.
//
// Reset semantics: only slots inserted since the last reset are zeroed
// (via populated_keys_), so reset is O(populated) not O(capacity).
//
// NOT YET WIRED IN. Callers continue to use DernaBeamMap. This struct is
// added in isolation to confirm compile correctness; Phase 2 will replace
// the hot-path task_map with this and benchmark.
struct DenseBeamTable {
    std::vector<BeamEntry> data_;          // sized to capacity once; never shrinks
    std::vector<int> populated_keys_;      // dense indices touched since last reset
    int capacity_ = 0;

    // Reset only the slots populated since the last reset. Optionally grow
    // capacity if the new pos demands more. Slot reset uses `score = inf` as
    // the sentinel for "unpopulated" to match BeamEntry's default ctor.
    void reset(int new_capacity) {
        if (new_capacity > capacity_) {
            data_.resize(new_capacity);
            capacity_ = new_capacity;
        }
        for (int k : populated_keys_) {
            BeamEntry& slot = data_[k];
            slot.variants.clear();
            slot = BeamEntry{};  // reset to default-constructed (score = inf)
        }
        populated_keys_.clear();
    }

    bool is_populated(int idx) const {
        return idx >= 0 && idx < capacity_ && data_[idx].score < inf;
    }

    BeamEntry& slot(int idx) { return data_[idx]; }
    const BeamEntry& slot(int idx) const { return data_[idx]; }

    int size() const { return (int)populated_keys_.size(); }

    // Mirror of free-function update_derna for the dense-table path.
    // Returns true iff a new slot was populated (caller can append idx to a
    // task_order log without a second is_populated() call).
    bool update(int idx, double score, const BeamEntry& en) {
        BeamEntry& ex = data_[idx];
        if (ex.score >= inf) {
            // Unpopulated: copy en into the slot, then run the seed routine.
            if (blk1_agg_diag_enabled()) {
                auto t0 = std::chrono::steady_clock::now();
                ex = en;
                _update_derna_seed(ex, score, en);
                populated_keys_.push_back(idx);
                g_blk1_agg_diag.seed_calls++;
                if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6 &&
                    !g_blk1_agg_hits_per_xy.empty()) {
                    const size_t k = (size_t)idx * 36 + (size_t)en.x * 6 + (size_t)en.y;
                    if (k < g_blk1_agg_hits_per_xy.size()) {
                        uint8_t& h = g_blk1_agg_hits_per_xy[k];
                        if (h < 255) h++;
                        if (h > g_blk1_agg_diag.per_pos_max_redundancy)
                            g_blk1_agg_diag.per_pos_max_redundancy = h;
                    }
                }
                auto t1 = std::chrono::steady_clock::now();
                g_blk1_agg_diag.total_update_ns +=
                    (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
                return true;
            }
            ex = en;
            _update_derna_seed(ex, score, en);
            populated_keys_.push_back(idx);
            return true;
        }
        if (blk1_agg_diag_enabled()) {
            auto t0 = std::chrono::steady_clock::now();
            // Detect whether (x,y) is a NEW variant within this entry.
            bool variant_existed = false;
            if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6) {
                variant_existed = (ex.variant_idx[en.x][en.y] >= 0);
            }
            _update_derna_existing(ex, score, en);
            if (variant_existed) g_blk1_agg_diag.existing_variant_calls++;
            else                  g_blk1_agg_diag.new_variant_calls++;
            if (en.x >= 0 && en.x < 6 && en.y >= 0 && en.y < 6 &&
                !g_blk1_agg_hits_per_xy.empty()) {
                const size_t k = (size_t)idx * 36 + (size_t)en.x * 6 + (size_t)en.y;
                if (k < g_blk1_agg_hits_per_xy.size()) {
                    uint8_t& h = g_blk1_agg_hits_per_xy[k];
                    if (h < 255) h++;
                    if (h > g_blk1_agg_diag.per_pos_max_redundancy)
                        g_blk1_agg_diag.per_pos_max_redundancy = h;
                }
            }
            auto t1 = std::chrono::steady_clock::now();
            g_blk1_agg_diag.total_update_ns +=
                (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
            return false;
        }
        _update_derna_existing(ex, score, en);
        return false;
    }
};

// VARDIST: defined here because it uses DenseBeamTable.
static void vardist_measure_dense(int pos, const char* table_name,
                                  const DenseBeamTable& td,
                                  const std::vector<int>& key_order,
                                  const std::vector<int>& protein) {
    if (!vardist_diag_enabled()) return;
    if (key_order.empty()) return;
    vardist_open_tsv_if_needed();
    if (!s_vardist_tsv.is_open()) return;

    uint64_t total_states = 0;
    uint64_t total_variants = 0;
    uint64_t sum_buckets_used = 0;
    uint64_t sum_subbuckets = 0;
    uint64_t sum_subbucket_count = 0;
    uint64_t max_variants = 0;

    for (int compact_idx : key_order) {
        const BeamEntry& e = td.slot(compact_idx);
        if (e.score >= inf) continue;
        if (e.variants.empty()) continue;
        total_states++;
        total_variants += e.variants.size();
        if (e.variants.size() > max_variants) max_variants = e.variants.size();

        uint8_t bucket_count[16] = {0};
        for (const auto& v : e.variants) {
            if (v.x < 0 || v.x >= 6 || v.y < 0 || v.y >= 6) continue;
            int paa = protein[e.a];
            int pbb = protein[e.b];
            int nlo = nucleotides[paa][v.x][e.i];
            int nro = nucleotides[pbb][v.y][e.j];
            int b = nlo * 4 + nro;
            if (b >= 0 && b < 16) bucket_count[b]++;
        }
        int buckets_used = 0;
        for (int b = 0; b < 16; ++b) {
            if (bucket_count[b] > 0) {
                buckets_used++;
                sum_subbuckets += bucket_count[b];
                sum_subbucket_count++;
            }
        }
        sum_buckets_used += buckets_used;
    }

    if (total_states == 0) return;
    double avg_variants = (double)total_variants / (double)total_states;
    double avg_buckets_used = (double)sum_buckets_used / (double)total_states;
    double avg_variants_per_bucket = (sum_subbucket_count > 0)
        ? (double)sum_subbuckets / (double)sum_subbucket_count : 0.0;

    s_vardist_tsv << pos << "\t" << table_name
                  << "\t" << total_states
                  << "\t" << total_variants
                  << "\t" << avg_variants
                  << "\t" << max_variants
                  << "\t" << avg_buckets_used
                  << "\t" << avg_variants_per_bucket << "\n";
}

// LCDS-style CS key: encodes (cs_inner_left, nuc_inner_L, nuc_outer_R, nuc_single_start, nuc_jm1, single_len).
// Packing layout (low → high bits):
//   [bit 0-1]  nuc_single_start  (×1)         — first nuc of right_S
//   [bit 2-3]  nuc_outer_R       (×4)         — inner-pair right
//   [bit 4-5]  nuc_inner_L       (×16)        — inner-pair left
//   [bit 6-7]  nuc_jm1           (×64)        — last nuc of right_S (= CS's last nuc)
//   [bits 8+]  len               (×256)       — right_S length
//   [highest]  inner_left        (×256*(SINGLE_MAX_LEN+1))
// nuc_jm1 was added (2026-04-28) to prevent two CS entries with same boundary
// nucleotides but different right-S last-nuc from collapsing — those two cases
// produce different downstream loop_e energies in the CS_to_C closure (via
// cs_best_by_y[cy].nuc_jm1, line ~1637) and must be kept separately.
// The low-6-bit decoder at line ~1505 (key % 64) is unaffected by this addition.
static inline int cs_get_index(int inner_left, int nuc_inner_L, int nuc_outer_R,
                               int nuc_single_start, int nuc_jm1, int len) {
    return inner_left * 256 * (SINGLE_MAX_LEN + 1)
         + len        * 256
         + nuc_jm1    * 64
         + nuc_inner_L * 16
         + nuc_outer_R *  4
         + nuc_single_start;
}


// Update CS table: keep best score; on tie, append (c_key,s_key) to bt_info; on worse, skip.
[[maybe_unused]] static void update_cs(DernaBeamMap& m, int key, double score,
                      int c_key, int s_key, const BeamEntry& tmpl) {
    const double eps = 1e-9;
    auto [it, inserted] = m.try_emplace(key, tmpl);

    if (inserted) {
        // New entry: set bt_info
        it->second.bt_info = {c_key, s_key};
        return;
    }

    // Key exists: decide what to do based on score comparison
    if (score < it->second.score - eps) {
        // Better: replace
        it->second = tmpl;
        it->second.bt_info = {c_key, s_key};
    } else if (score < it->second.score + eps) {
        // Tied: append if not duplicate
        bool dup = false;
        for (int k = 0; k + 1 < (int)it->second.bt_info.size(); k += 2) {
            if (it->second.bt_info[k] == c_key && it->second.bt_info[k+1] == s_key) {
                dup = true;
                break;
            }
        }
        if (!dup) {
            it->second.bt_info.push_back(c_key);
            it->second.bt_info.push_back(s_key);
        }
    }
    // else: worse, skip
}

// CAI values like 2.65e-314 are finite but denormal/uninitialized garbage; valid CAI is in a sane range.
constexpr double CAI_SUSPICIOUS_THRESHOLD = 1e-200;
static inline bool cai_looks_garbage(double cai) {
    return (cai != 0.0 && std::abs(cai) < CAI_SUSPICIOUS_THRESHOLD);
}

// Same invariant for C entries: catch which path writes garbage cai into curr_c.
// In release builds (NDEBUG) the assertions become no-ops; skip all work since the
// arithmetic itself is called from hot inner loops (Block 1 internal: ~660M/run).
static inline void check_c_entry_invariant(double lambda, const BeamEntry& ent, const char* path_name, int pos) {
#ifndef NDEBUG
    if (!std::isfinite(ent.cai) || !std::isfinite(ent.mfe)) {
        std::cerr << "[C invariant FAIL] path=" << path_name << " pos=" << pos << " cai=" << ent.cai << " mfe=" << ent.mfe << " a=" << ent.a << " b=" << ent.b << "\n";
        assert(false && "C entry cai/mfe must be finite");
    }
    if (cai_looks_garbage(ent.cai)) {
        std::cerr << "[C invariant FAIL] path=" << path_name << " pos=" << pos << " cai=" << ent.cai << " (garbage) mfe=" << ent.mfe << " a=" << ent.a << " b=" << ent.b << " i=" << static_cast<int>(ent.i) << " j=" << static_cast<int>(ent.j) << " x=" << ent.x << " y=" << ent.y << "\n";
        assert(false && "C entry cai looks like uninitialized/denormal garbage");
    }
    double expected = combined_score(lambda, ent.mfe, ent.cai);
    if (std::abs(ent.score - expected) >= 1e-8) {
        std::cerr << "[C invariant FAIL] path=" << path_name << " pos=" << pos << " score=" << ent.score << " expected=" << expected << "\n";
        assert(false && "C entry score must equal combined_score(lambda, mfe, cai)");
    }
#else
    (void)lambda; (void)ent; (void)path_name; (void)pos;
#endif
}

// Invariant check before writing to M1: catch bad CAI propagation (e.g. uninitialized/garbage cai
// producing artificially low score that overwrites correct state). Call before every update_derna into curr_m1.
static inline void check_m1_entry_invariant(double lambda, const BeamEntry& ent, const char* path_name) {
    (void)path_name;  // for documentation at call site (C->M1, M1->M1, M2->M1)
#ifndef NDEBUG
    assert(std::isfinite(ent.cai) && "M1 entry cai must be finite (bad propagation or uninitialized)");
    assert(std::isfinite(ent.mfe) && "M1 entry mfe must be finite");
    assert(!cai_looks_garbage(ent.cai) && "M1 entry cai looks like uninitialized/denormal garbage");
    double expected = combined_score(lambda, ent.mfe, ent.cai);
    assert(std::abs(ent.score - expected) < 1e-8 && "M1 entry score must equal combined_score(lambda, mfe, cai)");
#else
    (void)lambda; (void)ent;
#endif
}

// Prune to top-k by score (keep smallest; minimize combined score), with legality checks.
// LCDSfold's BeamPrune() does more than top-k; at minimum we enforce that each entry is
// position-consistent (i.e., its right boundary matches the table position) and internally sane.
// This prevents many "locally-good but globally-invalid" states from surviving.
enum class DernaTableKind {
    N,
    S,
    F,
    C,
    CS,
    M1,
    M2,
    Multi
};

static inline bool entry_has_basic_position_consistency(const BeamEntry& e, int pos, int n) {
    // 1. The most restrictive check first. If it doesn't end at 'pos', fail immediately.
    if (sigma(e.b, e.j) != pos) return false;

    // 2. The Unsigned Cast Trick: combines < 0 and >= bounds checks into one instruction.
    // If e.i is negative, casting to unsigned makes it huge, which naturally fails > 2u.
    if (static_cast<unsigned>(e.i) > 2u || static_cast<unsigned>(e.j) > 2u) return false;

    if (static_cast<unsigned>(e.a) >= static_cast<unsigned>(n) ||
        static_cast<unsigned>(e.b) >= static_cast<unsigned>(n)) return false;

    // 3. Because we just proved above that e.a >= 0 and e.i >= 0,
    // sigma(e.a, e.i) physically cannot be negative.
    // We can safely drop the 'left < 0' check and only test the upper bound.
    if (sigma(e.a, e.i) > pos) return false;

    return true;
}

static inline bool entry_is_position_consistent(const BeamEntry& e,
                                                int pos,
                                                int n,
                                                DernaTableKind kind,
                                                int aux_len = -1) {
    if (!entry_has_basic_position_consistency(e, pos, n)) return false;

    switch (kind) {
        case DernaTableKind::S: {
            // For tab_s[pos][l], the state must start exactly at pos-l+1.
            if (aux_len <= 0 || aux_len > SINGLE_MAX_LEN) return false;
            const int left = sigma(e.a, e.i);
            return left == pos - aux_len + 1;
        }
        case DernaTableKind::CS: {
            // CS explicitly stores the right single segment and the inner closed pair.
            if (e.cs_right_len <= 0 || e.cs_right_len > SINGLE_MAX_LEN) return false;
            if (e.cs_inner_key < 0 || e.cs_right_s_key < 0) return false;
            if (e.cs_single_start < 0 || e.cs_single_start > pos) return false;
            if (e.cs_inner_left < 0 || e.cs_inner_right < 0) return false;
            if (e.cs_inner_left > e.cs_inner_right) return false;

            // The BeamEntry geometry for CS is the right S segment itself.
            const int left = sigma(e.a, e.i);
            if (left != e.cs_single_start) return false;
            if (e.cs_single_start != pos - e.cs_right_len + 1) return false;
            if (e.cs_inner_right != e.cs_single_start - 1) return false;
            // MEK: commented out because same as line 158
            // if (e.cs_inner_left > e.cs_inner_right) return false;

            return true;
        }
        case DernaTableKind::N:
        case DernaTableKind::F:
        case DernaTableKind::C:
        case DernaTableKind::M1:
        case DernaTableKind::M2:
        case DernaTableKind::Multi:
        default:
            return true;
    }
}

// Prune log file: opened at start of fill_position_beam_tables, closed at end.
static std::ofstream s_prune_log;
// Current protein (for logging nuc choices); set during fill.
static const vector<int>* s_prune_protein = nullptr;

// Dual-beam diagnostic CSV: full pre-prune candidate dump with rank under
// combined-score, mfe-only, cai-only orderings. Opened only if env var
// DERNA_DUMP_PRUNES is set (path).
static std::ofstream s_dualbeam_csv;
static int s_dualbeam_step_id = 0;

static const char nuc_char[4] = {'A', 'C', 'G', 'U'};

// Forward decl: helper to dump full ranking before pruning.
static void dualbeam_dump_full(const DernaBeamMap& states,
                               int beamsize, int pos, int n,
                               const char* kind_name, int aux_len,
                               const vector<vector<double>>* best_f_prefix);

// Human-readable name for transition/manner (which transition produced this state).
static const char* manner_name(int bt) {
    switch (bt) {
        case (int)DernaManner::MANNER_NONEtoN: return "NONEtoN";
        case (int)DernaManner::MANNER_NONEtoF: return "NONEtoF";
        case (int)DernaManner::MANNER_NONEtoS: return "NONEtoS";
        case (int)DernaManner::MANNER_N_EtoN: return "N_EtoN";
        case (int)DernaManner::MANNER_NtoC: return "NtoC";
        case (int)DernaManner::MANNER_S_EtoS: return "S_EtoS";
        case (int)DernaManner::MANNER_CStoC: return "CStoC";
        case (int)DernaManner::MANNER_S_CStoC: return "S_CStoC";
        case (int)DernaManner::MANNER_C_StoCS: return "C_StoCS";
        case (int)DernaManner::MANNER_C_StoC: return "C_StoC";
        case (int)DernaManner::MANNER_S_C_StoC: return "S_C_StoC";
        case (int)DernaManner::MANNER_CtoC: return "CtoC";
        case (int)DernaManner::MANNER_S_CtoC: return "S_CtoC";
        case (int)DernaManner::MANNER_Multi_EtoMulti: return "Multi_EtoMulti";
        case (int)DernaManner::MANNER_MultitoC: return "MultitoC";
        case (int)DernaManner::MANNER_CtoM1: return "CtoM1";
        case (int)DernaManner::MANNER_M1_CtoM2: return "M1_CtoM2";
        case (int)DernaManner::MANNER_M1_EtoM1: return "M1_EtoM1";
        case (int)DernaManner::MANNER_M2toM1: return "M2toM1";
        case (int)DernaManner::MANNER_M2toMulti: return "M2toMulti";
        case (int)DernaManner::MANNER_S_M2toMulti: return "S_M2toMulti";
        case (int)DernaManner::MANNER_F_EtoF: return "F_EtoF";
        case (int)DernaManner::MANNER_CtoF: return "CtoF";
        case (int)DernaManner::MANNER_F_CtoF: return "F_CtoF";
        default: return "?";
    }
}

// Dual-beam diagnostic dump: emits one CSV row per candidate in `states`
// (BEFORE pruning), with rank under combined-score / MFE / CAI orderings.
// Combined-score = states[].score (lambda*mfe + (lambda-1)*cai), already what
// the DP uses. MFE-only ranks by .mfe ascending. CAI-only ranks by .cai
// ascending (i.e., lower CAI sum = "more negative" = closer to optimal CAI;
// remember internally CAI is stored as a negative log-prob sum).
static void dualbeam_dump_full(const DernaBeamMap& states,
                               int beamsize, int pos, int n,
                               const char* kind_name, int aux_len,
                               const vector<vector<double>>* best_f_prefix) {
    if (!s_dualbeam_csv.is_open() || states.empty()) return;
    const int N = (int)states.size();
    s_dualbeam_step_id++;
    const int step_id = s_dualbeam_step_id;

    // Materialize candidates so we can compute three rankings.
    struct Cand {
        int key;
        double score;   // combined
        double mfe;
        double cai;
        int a, b;
        int8_t i_, j_;
        int x, y;
        int left_pos, right_pos;
        int manner;
    };
    std::vector<Cand> cands;
    cands.reserve(N);
    for (const auto& kv : states) {
        const auto& e = kv.second;
        Cand c;
        c.key = kv.first;
        c.score = e.score;
        c.mfe = e.mfe;
        c.cai = e.cai;
        c.a = e.a; c.b = e.b;
        c.i_ = e.i; c.j_ = e.j;
        c.x = e.x; c.y = e.y;
        c.left_pos = (e.a >= 0) ? sigma(e.a, e.i) : -1;
        c.right_pos = (e.b >= 0) ? sigma(e.b, e.j) : pos;
        c.manner = e.backtrace_type;
        cands.push_back(c);
    }
    // Compute ranks: lower is better for all three (combined, mfe, cai).
    auto rank_under = [&](auto field) {
        std::vector<int> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) {
                      double va = field(cands[a]);
                      double vb = field(cands[b]);
                      if (va != vb) return va < vb;
                      return cands[a].key < cands[b].key;
                  });
        std::vector<int> rank(N);
        for (int r = 0; r < N; ++r) rank[idx[r]] = r + 1;
        return rank;
    };
    auto rank_combined = rank_under([](const Cand& c) { return c.score; });
    auto rank_mfe = rank_under([](const Cand& c) { return c.mfe; });
    auto rank_cai = rank_under([](const Cand& c) { return c.cai; });

    // Optional cumulative-score (for C/M1/M2/Multi where prune_cumulative is used).
    std::vector<int> rank_cum;
    if (best_f_prefix) {
        const int np = (int)best_f_prefix->size();
        auto cum_of = [&](const Cand& c) -> double {
            int prev = c.left_pos - 1;
            if (prev < 0 || prev >= np) return c.score;
            double best = inf;
            for (int y = 0; y < 4; ++y)
                if (std::isfinite((*best_f_prefix)[prev][y]))
                    best = std::min(best, (*best_f_prefix)[prev][y]);
            return (best < inf) ? best + c.score : c.score;
        };
        std::vector<int> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) {
                      double va = cum_of(cands[a]);
                      double vb = cum_of(cands[b]);
                      if (va != vb) return va < vb;
                      return cands[a].key < cands[b].key;
                  });
        rank_cum.assign(N, 0);
        for (int r = 0; r < N; ++r) rank_cum[idx[r]] = r + 1;
    }

    for (int idx = 0; idx < N; ++idx) {
        const Cand& c = cands[idx];
        s_dualbeam_csv
            << step_id << ',' << kind_name << ',' << pos << ',' << aux_len
            << ',' << beamsize << ',' << N
            << ',' << c.left_pos << ',' << c.right_pos
            << ',' << c.a << ',' << c.b
            << ',' << (int)c.i_ << ',' << (int)c.j_
            << ',' << c.x << ',' << c.y
            << ',' << c.key << ',' << c.manner
            << ',' << c.score << ',' << c.mfe << ',' << c.cai
            << ',' << rank_combined[idx]
            << ',' << rank_mfe[idx]
            << ',' << rank_cai[idx]
            << ',' << (rank_cum.empty() ? -1 : rank_cum[idx])
            << ',' << ((rank_combined[idx] <= beamsize) ? 1 : 0)
            << '\n';
    }
    (void)n;
}

// Cumulative-score pruning for the C table: rank by best_f_prefix[left_pos-1] + local_score.
// This matches LCDSfold's BeamPrune which uses bestF[i-1] + cand.score as the effective score,
// so that globally-promising C entries survive even if their local score is not the best.
// best_f_prefix[p] = best (min) F score at position p (0 if no F entry exists at p).
static void prune_cumulative(DernaBeamMap& states,
                             int beamsize,
                             int pos,
                             int n,
                             DernaTableKind kind,
                             const vector<vector<double>>& best_f_prefix,
                             [[maybe_unused]] const vector<int>& protein,
                             int aux_len = -1) {
    if (states.empty()) return;
    for (auto it = states.begin(); it != states.end();) {
        if (!std::isfinite(it->second.score) || !entry_is_position_consistent(it->second, pos, n, kind, aux_len))
            it = states.erase(it);
        else ++it;
    }
    // Diagnostic: dump full pre-prune ranking under combined / mfe / cai / cum.
    if (s_dualbeam_csv.is_open()) {
        const char* kn = "?";
        switch (kind) {
            case DernaTableKind::N: kn = "N"; break;
            case DernaTableKind::S: kn = "S"; break;
            case DernaTableKind::F: kn = "F"; break;
            case DernaTableKind::C: kn = "C"; break;
            case DernaTableKind::CS: kn = "CS"; break;
            case DernaTableKind::M1: kn = "M1"; break;
            case DernaTableKind::M2: kn = "M2"; break;
            case DernaTableKind::Multi: kn = "Multi"; break;
            default: break;
        }
        dualbeam_dump_full(states, beamsize, pos, n, kn, aux_len, &best_f_prefix);
    }
    if (beamsize <= 0 || states.size() <= (size_t)beamsize) return;

    const int np = (int)best_f_prefix.size();
    // Codon-compatible F-prefix lookup: when prev and left_pos are in the same
    // codon, F's RIGHT codon must equal candidate's LEFT codon (same amino = same
    // codon commitment). This is stricter than LCDSfold's nucleotide equality,
    // but appropriate for our codon-level state space where each entry commits
    // to specific codons for all 3 slots of each amino acid it touches.
    auto cum = [&](const BeamEntry& e) -> double {
        int left_pos = sigma(e.a, e.i);
        int prev = left_pos - 1;
        if (prev < 0 || prev >= np) return e.score;
        // LCDSfold-style: take min over all F nucleotide buckets at prev. Codon-consistency
        // between F's right codon and the candidate's left codon is ignored here (as in LCDSfold);
        // the pruning is an upper-bound heuristic for ranking.
        double best = inf;
        for (int y = 0; y < 4; ++y)
            if (std::isfinite(best_f_prefix[prev][y]))
                best = min(best, best_f_prefix[prev][y]);
        return (best < inf) ? best + e.score : e.score;
    };

    vector<pair<double, int>> vals;
    vals.reserve(states.size());
    for (const auto& kv : states) vals.push_back({cum(kv.second), kv.first});
    auto nth = vals.begin() + (beamsize - 1);
    nth_element(vals.begin(), nth, vals.end(),
                [](const auto& a, const auto& b) {
                    if (a.first != b.first) return a.first < b.first;
                    return a.second < b.second;
                });
    double threshold_score = nth->first;
    int threshold_key = nth->second;
    for (auto it = states.begin(); it != states.end();) {
        double sv = cum(it->second);
        if (sv > threshold_score || (sv == threshold_score && it->first > threshold_key))
            it = states.erase(it);
        else ++it;
    }
}

// Pack (closed_right, seg_len, s_right_key, c_inner_key) into a uint64_t for
// the per-pos "Block 1 fallback needed" set. closed_right uses 12 bits;
// seg_len uses 5 bits; s_right_key 24 bits; c_inner_key 23 bits = 64 bits.
// Sufficient for current SINGLE_MAX_LEN=30 and protein lengths up to ~1300.
static inline uint64_t pack_cs_blk1_unit(int closed_right, int seg_len,
                                         int s_right_key, int c_inner_key) {
    return ((uint64_t)(closed_right & 0xFFF) << 52) |
           ((uint64_t)(seg_len & 0x1F) << 47) |
           ((uint64_t)(s_right_key & 0xFFFFFF) << 23) |
           ((uint64_t)(c_inner_key & 0x7FFFFF));
}

// Hybrid prune for CS: keep the union of top-K by raw-score AND top-K by
// cumulative (bestF[cs_il-1] + score). The pure cumulative prune drops CS
// entries with high cs_il (large inner_left), because their bestF lookahead
// is from far back in the sequence and can't compete with low-cs_il entries
// whose score has already absorbed huge negative MFE from prior closures.
// Verified failure on P01707: an internal-loop CS entry at cs_il=30 with
// score=-339.659 was pruned (cum=-1938.05, threshold@K=-3386.6) even though
// its closure produces a globally optimal C entry. Block 1 historically
// recovered this combination by enumerating raw-score top entries directly.
//
// Hybrid prune restores raw-score top-K, lossless against Block 1's path on
// the validated 19-case sweep (see Phase 4 in feedback_cs_path_coverage.md).
//
// If dropped_blk1_units is non-null, the prune logs the packed
// (closed_right, seg_len, s_right_key, c_inner_key) tuple of every CS entry
// that was dropped — the consumer can use this to restrict Block 1's work to
// only the (closed_right, seg_len, s_key, c_key) units whose CS coverage was
// lost (selective Block 1).
static void prune_cs_hybrid(DernaBeamMap& states,
                            int beamsize,
                            int pos,
                            int n,
                            const vector<vector<double>>& best_f_prefix,
                            const vector<int>& protein,
                            std::unordered_set<uint64_t>* dropped_blk1_units = nullptr) {
    if (states.empty()) return;
    // Validity: drop non-finite scores or position-inconsistent entries.
    for (auto it = states.begin(); it != states.end();) {
        if (!std::isfinite(it->second.score) ||
            !entry_is_position_consistent(it->second, pos, n, DernaTableKind::CS, -1))
            it = states.erase(it);
        else ++it;
    }
    if (beamsize <= 0 || states.size() <= (size_t)(2 * beamsize)) return;

    const int np = (int)best_f_prefix.size();
    auto cum = [&](const BeamEntry& e) -> double {
        int left_pos = sigma(e.a, e.i);
        int prev = left_pos - 1;
        if (prev < 0 || prev >= np) return e.score;
        double best = inf;
        for (int y = 0; y < 4; ++y)
            if (std::isfinite(best_f_prefix[prev][y]))
                best = std::min(best, best_f_prefix[prev][y]);
        return (best < inf) ? best + e.score : e.score;
    };
    (void)protein;

    // Compute (raw, cum, key) for every entry and sort once for each metric.
    std::vector<int> all_keys; all_keys.reserve(states.size());
    std::vector<double> raw_scores; raw_scores.reserve(states.size());
    std::vector<double> cum_scores; cum_scores.reserve(states.size());
    for (const auto& kv : states) {
        all_keys.push_back(kv.first);
        raw_scores.push_back(kv.second.score);
        cum_scores.push_back(cum(kv.second));
    }
    // Find top-K by raw via partial sort of indices.
    std::vector<int> idx_raw(all_keys.size()), idx_cum(all_keys.size());
    std::iota(idx_raw.begin(), idx_raw.end(), 0);
    std::iota(idx_cum.begin(), idx_cum.end(), 0);
    auto raw_lt = [&](int a, int b) {
        if (raw_scores[a] != raw_scores[b]) return raw_scores[a] < raw_scores[b];
        return all_keys[a] < all_keys[b];
    };
    auto cum_lt = [&](int a, int b) {
        if (cum_scores[a] != cum_scores[b]) return cum_scores[a] < cum_scores[b];
        return all_keys[a] < all_keys[b];
    };
    std::nth_element(idx_raw.begin(), idx_raw.begin() + beamsize, idx_raw.end(), raw_lt);
    std::nth_element(idx_cum.begin(), idx_cum.begin() + beamsize, idx_cum.end(), cum_lt);
    std::unordered_set<int> keep_keys;
    keep_keys.reserve(2 * beamsize);
    for (int i = 0; i < beamsize; ++i) keep_keys.insert(all_keys[idx_raw[i]]);
    for (int i = 0; i < beamsize; ++i) keep_keys.insert(all_keys[idx_cum[i]]);
    for (auto it = states.begin(); it != states.end();) {
        if (keep_keys.count(it->first) == 0) {
            // Log this tuple as needing Block 1 fallback at the next pos.
            if (dropped_blk1_units) {
                const BeamEntry& e = it->second;
                if (e.cs_inner_key >= 0 && e.cs_right_s_key >= 0 &&
                    e.cs_right_len > 0 && e.cs_inner_right >= 0) {
                    dropped_blk1_units->insert(pack_cs_blk1_unit(
                        e.cs_inner_right, e.cs_right_len,
                        e.cs_right_s_key, e.cs_inner_key));
                }
            }
            it = states.erase(it);
        }
        else ++it;
    }
}

// Lookahead-aware CS prune (Approach a): rank CS entries by an estimate of
// their post-closure score, not raw or cumulative.
//
// Lookahead score:
//   L(cs) = bestF[cs.cs_inner_left - 1] + cs.score + min_sl_at[cs.cs_inner_left - 1]
//
// where min_sl_at[p] = min over (seg_len, key) of bestS[p][seg_len][key].score
// — i.e. the best left-S contribution that would be available if we eventually
// close cs into C via S_CS_to_C with S_left ending at cs_il-1. Since S scores
// are non-positive, min_sl_at[p] is a stricter (better) lookahead than 0.
//
// This is a stricter (more pessimistic, i.e. larger) ranking than raw score
// since it adds bestF[cs_il-1] and min_sl_at[cs_il-1] (both non-positive),
// hopefully selecting CS entries that produce better C closures.
//
// Used in conjunction with raw and cumulative metrics to form a 3-way hybrid.
// If dropped_blk1_units is non-null, dropped tuples are logged as in
// prune_cs_hybrid.
static void prune_cs_hybrid_lookahead(
        DernaBeamMap& states,
        int beamsize,
        int pos,
        int n,
        const vector<vector<double>>& best_f_prefix,
        const vector<double>& min_sl_at,
        const vector<int>& protein,
        std::unordered_set<uint64_t>* dropped_blk1_units = nullptr) {
    if (states.empty()) return;
    for (auto it = states.begin(); it != states.end();) {
        if (!std::isfinite(it->second.score) ||
            !entry_is_position_consistent(it->second, pos, n, DernaTableKind::CS, -1))
            it = states.erase(it);
        else ++it;
    }
    if (beamsize <= 0 || states.size() <= (size_t)(3 * beamsize)) return;
    (void)protein;

    const int np = (int)best_f_prefix.size();
    auto cum = [&](const BeamEntry& e) -> double {
        int left_pos = sigma(e.a, e.i);
        int prev = left_pos - 1;
        if (prev < 0 || prev >= np) return e.score;
        double best = inf;
        for (int y = 0; y < 4; ++y)
            if (std::isfinite(best_f_prefix[prev][y]))
                best = std::min(best, best_f_prefix[prev][y]);
        return (best < inf) ? best + e.score : e.score;
    };
    auto lookahead = [&](const BeamEntry& e) -> double {
        // cs_il = sigma(e.a, e.i); seam = cs_il - 1 = where S_left ends.
        // F prefix is at cs_il - 1 (same as cum, since S_left "covers" that
        // position when present). Use the cum baseline + min_sl bonus.
        double base = cum(e);
        if (!std::isfinite(base)) return base;
        int cs_il = sigma(e.a, e.i);
        int seam = cs_il - 1;
        if (seam >= 0 && seam < (int)min_sl_at.size() && std::isfinite(min_sl_at[seam])) {
            // S_left is non-positive; min_sl_at[seam] tightens the lookahead.
            base += std::min(0.0, min_sl_at[seam]);
        }
        return base;
    };

    std::vector<int> all_keys; all_keys.reserve(states.size());
    std::vector<double> raw_scores; raw_scores.reserve(states.size());
    std::vector<double> cum_scores; cum_scores.reserve(states.size());
    std::vector<double> la_scores; la_scores.reserve(states.size());
    for (const auto& kv : states) {
        all_keys.push_back(kv.first);
        raw_scores.push_back(kv.second.score);
        cum_scores.push_back(cum(kv.second));
        la_scores.push_back(lookahead(kv.second));
    }
    auto take_topk = [&](const std::vector<double>& scs, std::unordered_set<int>& keep) {
        std::vector<int> idx(all_keys.size());
        std::iota(idx.begin(), idx.end(), 0);
        auto lt = [&](int a, int b) {
            if (scs[a] != scs[b]) return scs[a] < scs[b];
            return all_keys[a] < all_keys[b];
        };
        std::nth_element(idx.begin(), idx.begin() + beamsize, idx.end(), lt);
        for (int i = 0; i < beamsize; ++i) keep.insert(all_keys[idx[i]]);
    };
    std::unordered_set<int> keep_keys;
    keep_keys.reserve(3 * beamsize);
    take_topk(raw_scores, keep_keys);
    take_topk(cum_scores, keep_keys);
    take_topk(la_scores,  keep_keys);

    for (auto it = states.begin(); it != states.end();) {
        if (keep_keys.count(it->first) == 0) {
            if (dropped_blk1_units) {
                const BeamEntry& e = it->second;
                if (e.cs_inner_key >= 0 && e.cs_right_s_key >= 0 &&
                    e.cs_right_len > 0 && e.cs_inner_right >= 0) {
                    dropped_blk1_units->insert(pack_cs_blk1_unit(
                        e.cs_inner_right, e.cs_right_len,
                        e.cs_right_s_key, e.cs_inner_key));
                }
            }
            it = states.erase(it);
        }
        else ++it;
    }
}

// Hybrid prune for SCLeft: keep the union of top-(M*K) by raw-score and
// top-(M*K) by cumulative-score (raw + best_f_prefix[outer_left_nuc - 1]).
//
// Why hybrid + multiplier: SCLeft entries have unique (negative) keys so they
// don't merge; without a cap, |SCLeft| blows up. Pure raw-score top-K dropped
// the optimum's predecessor for some instances (P01707). Pure cumulative-score
// top-K biases against entries with outer pair near position 0 (tiny F-prefix
// → looks worse than deep-F-prefix entries). Hybrid at 1*K still missed P01707
// because the optimum's SCLeft was below rank K by BOTH metrics. Multiplier
// M relaxes the cap to admit more candidates while still bounding runtime.
constexpr int SC_LEFT_PRUNE_MULTIPLIER = 3;
static void prune_sc_left_hybrid(DernaBeamMap& states,
                                 int beamsize,
                                 const vector<vector<double>>& best_f_prefix) {
    if (states.empty()) return;
    // Validity: drop non-finite scores.
    for (auto it = states.begin(); it != states.end();) {
        if (!std::isfinite(it->second.score)) it = states.erase(it);
        else ++it;
    }
    if (beamsize <= 0) return;
    const int eff_k = beamsize * SC_LEFT_PRUNE_MULTIPLIER;
    if (states.size() <= (size_t)(2 * eff_k)) return;

    const int np = (int)best_f_prefix.size();
    auto cum = [&](const BeamEntry& e) -> double {
        // sigma(e.a, e.i) = s_left_start; eventual outer_left_nuc = s_left_start - 1.
        // F-prefix lookahead is at outer_left_nuc - 1 = sigma(e.a, e.i) - 2.
        int prev = sigma(e.a, e.i) - 2;
        if (prev < 0 || prev >= np) return e.score;
        double best = inf;
        for (int y = 0; y < 4; ++y)
            if (std::isfinite(best_f_prefix[prev][y]))
                best = min(best, best_f_prefix[prev][y]);
        return (best < inf) ? best + e.score : e.score;
    };

    // Build keep-set: top-(M*K) by raw, plus top-(M*K) by cumulative.
    std::unordered_set<int> keep;
    keep.reserve((size_t)(2 * eff_k));
    auto take_top_k = [&](auto score_fn) {
        vector<pair<double, int>> vals;
        vals.reserve(states.size());
        for (const auto& kv : states) vals.push_back({score_fn(kv.second), kv.first});
        auto nth = vals.begin() + (eff_k - 1);
        std::nth_element(vals.begin(), nth, vals.end(),
                         [](const auto& a, const auto& b) {
                             if (a.first != b.first) return a.first < b.first;
                             return a.second < b.second;
                         });
        for (auto it = vals.begin(); it <= nth; ++it) keep.insert(it->second);
    };
    take_top_k([](const BeamEntry& e) { return e.score; });
    take_top_k([&](const BeamEntry& e) { return cum(e); });
    for (auto it = states.begin(); it != states.end();) {
        if (keep.count(it->first) == 0) it = states.erase(it);
        else ++it;
    }
}

static void prune_beam_derna_checked(DernaBeamMap& states,
                                    int beamsize,
                                    int pos,
                                    int n,
                                    DernaTableKind kind,
                                    int aux_len = -1,
                                    bool do_log = true) {
    if (states.empty()) return;

    // 1) Drop illegal / inconsistent entries first.
    for (auto it = states.begin(); it != states.end();) {
        if (!entry_is_position_consistent(it->second, pos, n, kind, aux_len)) it = states.erase(it);
        else ++it;
    }

    // Diagnostic: dump full pre-prune ranking.
    if (s_dualbeam_csv.is_open()) {
        const char* kn = "?";
        switch (kind) {
            case DernaTableKind::N: kn = "N"; break;
            case DernaTableKind::S: kn = "S"; break;
            case DernaTableKind::F: kn = "F"; break;
            case DernaTableKind::C: kn = "C"; break;
            case DernaTableKind::CS: kn = "CS"; break;
            case DernaTableKind::M1: kn = "M1"; break;
            case DernaTableKind::M2: kn = "M2"; break;
            case DernaTableKind::Multi: kn = "Multi"; break;
            default: break;
        }
        dualbeam_dump_full(states, beamsize, pos, n, kn, aux_len, nullptr);
    }

    // 2) Top-k by score (minimize) [MEK: and break ties using the key].
    if (beamsize > 0 && states.size() > (size_t)beamsize) {
        vector<pair<double, int>> vals;
        vals.reserve(states.size());
        for (const auto& kv : states) vals.push_back({kv.second.score, kv.first});

        auto nth = vals.begin() + (beamsize - 1);
        nth_element(vals.begin(), nth, vals.end(),
                    [](const auto& a, const auto& b) {
                        if (a.first != b.first) return a.first < b.first;
                        return a.second < b.second;
                    });
        double threshold = nth->first;
        int threshold_key = nth->second;

        for (auto it = states.begin(); it != states.end();)
        {
            double sv = it->second.score;
            if (sv > threshold || (sv == threshold && it->first > threshold_key))
                it = states.erase(it);
            else ++it;
        }
    }

    // Print top 20 kept entries (to prune log when open) for every table with entries, except S.
    if (do_log && s_prune_log.is_open() && !states.empty() && kind != DernaTableKind::S) {
        const char* kind_name = "?";
        switch (kind) {
            case DernaTableKind::N: kind_name = "N"; break;
            case DernaTableKind::S: kind_name = "S"; break;
            case DernaTableKind::F: kind_name = "F"; break;
            case DernaTableKind::C: kind_name = "C"; break;
            case DernaTableKind::CS: kind_name = "CS"; break;
            case DernaTableKind::M1: kind_name = "M1"; break;
            case DernaTableKind::M2: kind_name = "M2"; break;
            case DernaTableKind::Multi: kind_name = "Multi"; break;
            default: break;
        }
        vector<pair<double, int>> by_score;
        by_score.reserve(states.size());
        for (const auto& kv : states) by_score.push_back({kv.second.score, kv.first});
        sort(by_score.begin(), by_score.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        const int nprint = (int)min(size_t(20), by_score.size());
        s_prune_log << "[prune pos=" << pos << " " << kind_name;
        if (aux_len >= 0) s_prune_log << " len=" << aux_len;
        s_prune_log << "] top " << nprint << " (score, mfe, cai, a, b, i, j, x, y, n, nucL, nucR, manner, from_keys):\n";
        for (int i = 0; i < nprint; ++i) {
            auto it = states.find(by_score[i].second);
            if (it != states.end()) {
                const BeamEntry& e = it->second;
                int key = it->first;
                s_prune_log << "  " << (i + 1) << ": score=" << e.score << " mfe=" << e.mfe << " cai=" << e.cai
                            << " a=" << e.a << " b=" << e.b << " i=" << static_cast<int>(e.i) << " j=" << static_cast<int>(e.j) << " x=" << e.x << " y=" << e.y << " n=" << n;
                if (s_prune_protein && e.a >= 0 && e.a < n && e.b >= 0 && e.b < n && (size_t)e.a < s_prune_protein->size() && (size_t)e.b < s_prune_protein->size()) {
                    int pa = (*s_prune_protein)[e.a], pb = (*s_prune_protein)[e.b];
                    int nucL = nucleotides[pa][e.x][e.i], nucR = nucleotides[pb][e.y][e.j];
                    s_prune_log << " nucL=" << (nucL >= 0 && nucL < 4 ? nuc_char[nucL] : '?')
                               << " nucR=" << (nucR >= 0 && nucR < 4 ? nuc_char[nucR] : '?');
                } else {
                    s_prune_log << " nucL=? nucR=?";
                }
                s_prune_log << " manner=" << manner_name(e.backtrace_type) << " key=" << key;
                if (!e.bt_info.empty()) {
                    s_prune_log << " from_keys=[";
                    for (int k = 0; k < e.bt_info.size(); ++k) s_prune_log << (k ? "," : "") << e.bt_info[k];
                    s_prune_log << "]";
                } else {
                    s_prune_log << " from_keys=[]";
                }
                if (kind == DernaTableKind::CS) {
                    s_prune_log << " cs_right_len=" << e.cs_right_len << " cs_single_start=" << e.cs_single_start
                        << " cs_inner_left=" << e.cs_inner_left << " cs_inner_right=" << e.cs_inner_right
                        << " n_from_keys=" << (e.bt_info.size() / 2);
                }
                s_prune_log << "\n";
            }
        }
    }
}


// Helper: is position the last nucleotide of a codon?
static inline bool is_last_nuc(int pos) { return pos % 3 == 2; }

// LCDSfold special-hairpin lengths (tri/tetra/hexa loops: 5/6/8 nucleotide sub-strings
// including both closing bases). These are the full set of keys in hairpinE.
static const int SP_HP_LENS[] = {5, 6, 8};

// Enumerate codon assignments for amino acids [aa_lo .. aa_hi] (inclusive) that produce
// the target nucleotide sequence `hp` (length L) spanning nucleotide positions
// [first_idx .. last_idx] (last_idx = first_idx + L - 1). For each matching assignment,
// call `on_match(codon_choices, cai_interior)` with:
//   codon_choices : vector<int> of length (aa_hi - aa_lo + 1) holding the chosen codon index per aa
//   cai_interior  : sum of codon_cai[protein[aa]][codon_choices[k]] for every aa fully inside the loop,
//                   which LCDSfold accumulates via the `isLastNuc(x)` branch in initialize_Special_HP_LD.
// NOTE: In this codebase's scoring convention, the closing-pair codon CAIs (aa_lo and aa_hi) are
// added separately by the caller via add_hairpin_CAI_2 so that we do not double-count. Here we
// only emit CAIs for aa's strictly between aa_lo and aa_hi (the loop-interior aa's).
template <typename F>
static void enum_sp_hp_codons(int aa_lo, int aa_hi, int first_idx, int last_idx,
                              const std::string& hp, const std::vector<int>& protein,
                              F&& on_match) {
    const int n_aa = aa_hi - aa_lo + 1;
    std::vector<int> choices(n_aa, 0);
    // Per-aa: list of codon indices consistent with the target substring at the aa's positions.
    std::vector<std::vector<int>> compatible(n_aa);
    static const char to_char_local[4] = {'A','C','G','U'};
    for (int k = 0; k < n_aa; ++k) {
        int aa = aa_lo + k;
        int pa = protein[aa];
        int ncod = n_codon[pa];
        int aa_start = 3 * aa;              // first nt index of this aa
        for (int c = 0; c < ncod; ++c) {
            bool ok = true;
            for (int s = 0; s < 3; ++s) {
                int nt_idx = aa_start + s;
                if (nt_idx < first_idx || nt_idx > last_idx) continue; // not constrained by hp
                int want = hp[nt_idx - first_idx];
                int have = nucleotides[pa][c][s];
                if (to_char_local[have] != (char)want) { ok = false; break; }
            }
            if (ok) compatible[k].push_back(c);
        }
        if (compatible[k].empty()) return;
    }
    // Cartesian product over compatible codon choices.
    std::function<void(int)> rec = [&](int k) {
        if (k == n_aa) {
            // Compute cai_interior: aa's strictly between aa_lo and aa_hi (interior) - but only aa's whose
            // is_last_nuc (nt index 3*aa+2) lies in [first_idx, last_idx] — otherwise the codon CAI is
            // NOT charged by LCDSfold initialize_Special_HP_LD (which only adds when x is inside loop span).
            //
            // Actually LCDSfold adds CAI for every aa whose last nuc (3*aa+2) is in [first_idx, last_idx].
            // For boundary aa's (aa_lo, aa_hi), that is the code handled separately by add_hairpin_CAI_2
            // at the caller. Here we emit the sum only for aa's whose last nuc lies strictly inside the loop
            // and is NOT the boundary aa (to avoid double count with caller's add_hairpin_CAI_2).
            double cai_int = 0.0;
            for (int kk = 0; kk < n_aa; ++kk) {
                int aa = aa_lo + kk;
                int last_nt = 3 * aa + 2;
                if (last_nt < first_idx || last_nt > last_idx) continue;
                if (aa == aa_lo || aa == aa_hi) continue; // boundary aa CAI handled by caller
                cai_int += codon_cai[protein[aa]][choices[kk]];
            }
            on_match(choices, cai_int);
            return;
        }
        for (int c : compatible[k]) {
            choices[k] = c;
            rec(k + 1);
        }
    };
    rec(0);
}

// ---- Bulge / internal loop scoring ----
// Map PositionBeam's (outer pair, inner pair, #unpaired left/right) to Turner/Vienna single-loop scoring.
//
// In LCDSfold, `v_score_single(i,j,p,q, ...)` depends on mismatch nucleotides adjacent to the pairs.
// PositionBeam does not explicitly track those adjacent unpaired identities (it only tracks paired bases
// and segment lengths). To remain drop-in and consistent with PositionBeam's state representation, we
// use a conservative approximation for mismatch-dependent terms by setting adjacent nucleotides on the
// loop side to the paired nucleotides (this matches stacking exactly and preserves bulge sizing; for
// general internal loops it matches loop-size + asymmetry penalties but approximates mismatch terms).
//
// IMPORTANT: This change fixes the *model mismatch* (stack/bulge/internal cases, size/asymmetry penalties,
// special 1x1/2x1/2x2 cases) compared to the prior bulge-only approximation. If you later enrich the
// state to track adjacent unpaired bases, you should wire those in below.

// Precomputed size/asymmetry term for generic internal loops (ns>=1, n1>0, n2>0).
// sil_size_term[n1][n2] = internal_loop37[u or capped] + min(MAX_NINIO, (nl-ns)*ninio37).
// Built once at program start; avoids repeated log() calls on hot Block 1 path.
static const int SIL_TABLE_DIM = SINGLE_MAX_LEN + 2;
static double sil_size_term[SIL_TABLE_DIM * SIL_TABLE_DIM];
static bool   sil_size_inited = false;
static void init_sil_size_table() {
    for (int n1 = 0; n1 < SIL_TABLE_DIM; ++n1) {
        for (int n2 = 0; n2 < SIL_TABLE_DIM; ++n2) {
            int u  = n1 + n2;
            int nl = (n1 > n2) ? n1 : n2;
            int ns = (n1 > n2) ? n2 : n1;
            int energy;
            if (ns == 1 && nl >= 3) {
                energy = (nl + 1 <= MAXLOOP) ? internal_loop37[nl + 1]
                                             : (internal_loop37[30] + (int)(lxc37 * log((nl + 1) / 30.0)));
                energy += std::min((int)MAX_NINIO, (nl - ns) * ninio37);
            } else if (ns >= 2) {
                energy = (u <= MAXLOOP) ? internal_loop37[u]
                                        : (internal_loop37[30] + (int)(lxc37 * log(u / 30.0)));
                energy += std::min((int)MAX_NINIO, (nl - ns) * ninio37);
            } else {
                energy = 0;
            }
            sil_size_term[n1 * SIL_TABLE_DIM + n2] = (double)energy;
        }
    }
    sil_size_inited = true;
}
namespace { struct SilSizeTableInit { SilSizeTableInit() { init_sil_size_table(); } };
static SilSizeTableInit g_sil_size_init; }

[[maybe_unused]] static inline double score_bulge_or_internal_loop([[maybe_unused]] double lambda,
                                                  int nuc_lo, int nuc_ro,
                                                  int nuc_li, int nuc_ri,
                                                  int left_unpaired, int right_unpaired) {
    // Normalize to the Vienna-style n1/n2 convention used by v_score_single:
    // n1 = # unpaired on 5' side (between outer-left i and inner-left p)
    // n2 = # unpaired on 3' side (between inner-right q and outer-right j)
    const int n1 = std::max(0, left_unpaired);
    const int n2 = std::max(0, right_unpaired);

    // Pair "types" using the same encoding already used throughout DERNA/Zuker (BP_pair).
    // Outer pair is (nuc_lo, nuc_ro). Inner pair in v_score_single is (nucq, nucp) = (nuc_ri, nuc_li).
    const int type  = BP_pair[nuc_lo + 1][nuc_ro + 1];
    const int type2 = BP_pair[nuc_ri + 1][nuc_li + 1];

    // If either cannot pair, treat as infeasible.
    if (type == 0 || type2 == 0) return inf;

    // Stacking (no unpaired on either side): exactly match LCDSfold behavior. Return raw energy (caller applies lambda).
    if (n1 == 0 && n2 == 0) {
        return (double)Zuker::stacking(nuc_lo, nuc_ro, nuc_li, nuc_ri);
    }

    // Bulge (one side has 0 unpaired): use the same bulge_loop routine. Return raw energy.
    if (n1 == 0 || n2 == 0) {
        const int nl = (n1 > n2) ? n1 : n2;
        return (double)Zuker::bulge_loop(nuc_lo, nuc_ro, nuc_li, nuc_ri, nl);
    }

    // Internal loop (both sides have unpaired): implement Turner/Vienna-style single-loop energy.
    // Because PositionBeam does not track the adjacent unpaired bases, we approximate mismatch-dependent
    // terms by setting adjacent nucleotides on the loop side to the paired nucleotides.
    //
    // The tables in params/constants.h are indexed in the same style as the LinearCDSfold utility:
    //   - base indices are (1=A,2=C,3=G,4=U), so we use +1.
    const int nuci  = nuc_lo + 1;
    const int nucj  = nuc_ro + 1;
    const int nucp  = nuc_li + 1;
    const int nucq  = nuc_ri + 1;

    // Adjacent bases on loop side (approximated).
    const int nuci1  = nuci;
    const int nucj_1 = nucj;
    const int nucp_1 = nucp;
    const int nucq1  = nucq;

    int nl = (n1 > n2) ? n1 : n2;
    int ns = (n1 > n2) ? n2 : n1;

    int energy = 0;

    // Special cases: 1x1, 2x1, 2x2, 1xn, 2x3 follow the same structure as v_score_single. Return raw energy.
    if (ns == 1) {
        if (nl == 1) {
            energy = int11_37[type][type2][nuci1][nucj_1];
            return (double)energy;
        }
        if (nl == 2) {
            if (n1 == 1)
                energy = int21_37[type][type2][nuci1][nucq1][nucj_1];
            else
                energy = int21_37[type2][type][nucq1][nuci1][nucp_1];
            return (double)energy;
        }
        // Size+asymmetry term precomputed in sil_size_term[n1][n2].
        energy = (int)sil_size_term[n1 * SIL_TABLE_DIM + n2];
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type2][nucq1][nucp_1];
        return (double)energy;
    }

    if (ns == 2) {
        if (nl == 2) {
            energy = int22_37[type][type2][nuci1][nucp_1][nucq1][nucj_1];
            return (double)energy;
        }
        if (nl == 3) {
            energy = internal_loop37[5] + ninio37;
            energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type2][nucq1][nucp_1];
            return (double)energy;
        }
    }

    {
        energy = (int)sil_size_term[n1 * SIL_TABLE_DIM + n2];
        energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type2][nucq1][nucp_1];
    }

    return (double)energy;
}

// Full Turner/Vienna single-loop scoring with explicit mismatch-adjacent bases.
// All nucleotide arguments are in PositionBeam encoding 0..3 (A,C,G,U).
// If any of nuc_i1/nuc_jm1/nuc_pm1/nuc_qp1 are negative, we fall back to using the paired bases.
static inline double score_single_loop_with_mismatch([[maybe_unused]] double lambda,
                                                    int nuc_lo, int nuc_ro,
                                                    int nuc_li, int nuc_ri,
                                                    int left_unpaired, int right_unpaired,
                                                    int nuc_i1, int nuc_jm1,
                                                    int nuc_pm1, int nuc_qp1) {
    const int n1 = std::max(0, left_unpaired);
    const int n2 = std::max(0, right_unpaired);

    const int type  = BP_pair[nuc_lo + 1][nuc_ro + 1];
    const int type2 = BP_pair[nuc_ri + 1][nuc_li + 1];
    if (type == 0 || type2 == 0) return inf;

    if (n1 == 0 && n2 == 0) {
        return (double)Zuker::stacking(nuc_lo, nuc_ro, nuc_li, nuc_ri);
    }

    if (n1 == 0 || n2 == 0) {
        const int nl = (n1 > n2) ? n1 : n2;
        return (double)Zuker::bulge_loop(nuc_lo, nuc_ro, nuc_li, nuc_ri, nl);
    }

    // Internal loop: use Turner/Vienna v_score_single logic with 37 tables.
    // Convert to 1..4 indices for mismatch tables.
    const int nuci  = nuc_lo + 1;
    const int nucj  = nuc_ro + 1;
    const int nucp  = nuc_li + 1;  // inner-left base at p
    const int nucq  = nuc_ri + 1;  // inner-right base at q

    // Adjacent bases on loop side (use provided if available; otherwise fall back to paired bases).
    const int nuci1  = (nuc_i1  >= 0) ? (nuc_i1  + 1) : nuci;   // i+1
    const int nucj_1 = (nuc_jm1 >= 0) ? (nuc_jm1 + 1) : nucj;   // j-1
    const int nucp_1 = (nuc_pm1 >= 0) ? (nuc_pm1 + 1) : nucp;   // p-1
    const int nucq1  = (nuc_qp1 >= 0) ? (nuc_qp1 + 1) : nucq;   // q+1

    int nl = (n1 > n2) ? n1 : n2;
    int ns = (n1 > n2) ? n2 : n1;

    int energy = 0;

    if (ns == 1) {
        if (nl == 1) {
            energy = int11_37[type][type2][nuci1][nucj_1];
            return (double)energy;
        }
        if (nl == 2) {
            if (n1 == 1)
                energy = int21_37[type][type2][nuci1][nucq1][nucj_1];
            else
                energy = int21_37[type2][type][nucq1][nuci1][nucp_1];
            return (double)energy;
        }
        // Size+asymmetry term precomputed in sil_size_term[n1][n2].
        energy = (int)sil_size_term[n1 * SIL_TABLE_DIM + n2];
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type2][nucq1][nucp_1];
        return (double)energy;
    }

    if (ns == 2) {
        if (nl == 2) {
            energy = int22_37[type][type2][nuci1][nucp_1][nucq1][nucj_1];
            return (double)energy;
        }
        if (nl == 3) {
            energy = internal_loop37[5] + ninio37;
            energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type2][nucq1][nucp_1];
            return (double)energy;
        }
    }

    {
        energy = (int)sil_size_term[n1 * SIL_TABLE_DIM + n2];
        energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type2][nucq1][nucp_1];
    }

    return (double)energy;
}

// Unique key for (a,b,i,j,x,y) so table entries don't collide.
// Encodes a,b in [0,n), i,j in [0,3), x,y in [0, 6) (max codons per aa).
//
// NOTE (vs LCDSfold): LCDSfold keys by (i, nuci, nucj) = position + boundary nucleotides,
// so different (x,y) that yield the same nucL,nucR collapse into one entry. Here we key by
// (a,b,i,j,x,y), so multiple entries can have the same (a,b,i,j) and same nucL,nucR but
// different (x,y). They consume multiple beam slots and can cause fewer distinct cases
// MERGED-STATE KEY: structural position (a,b,i,j) only — drops (x,y) so all codon-pair variants
// of the same structural state share one beam slot. Variants are stored in BeamEntry.variants.
[[maybe_unused]] static inline int index_struct(int a, int b, int i, int j, int n) {
    return (a * n + b) * 9 + (3 * i + j);
}

// Nucleotide-level beam key (matching LCDSfold's GetIndex).
// Codons producing same boundary nucleotides collapse to one beam slot.
// Codon info preserved in variants. 16 slots per position.
static inline int codon_beam_key(int left_pos, int nuc_left, int nuc_right) {
    return left_pos * 16 + nuc_left * 4 + nuc_right;
}
static inline int codon_beam_key_left_pos(int key) { return key / 16; }

// Legacy codon-augmented key (kept for decode helpers and CS-special path).
[[maybe_unused]] static inline int index_derna(int a, int b, int i, int j, int x, int y, int n) {
    return ((a * n + b) * 9 + (3 * i + j)) * 36 + x * 6 + y;
}

// Decode (a, b) from index_derna key (for traceback fallback when key was pruned).
[[maybe_unused]] static inline void index_derna_decode_ab(int key, int n, int& a, int& b) {
    int K = key / 36;
    int ab = K / 9;
    a = ab / n;
    b = ab % n;
}

// Right-bound position of the C state encoded by key: sigma(b, j) = 3*b + j.
[[maybe_unused]] static inline int index_derna_key_to_close(int key, int n) {
    int K = key / 36;
    int ab = K / 9;
    int b = ab % n;
    int slot = K % 9;  // (3*i+j)
    int j = slot % 3;
    return 3 * b + j;
}

// Right-bound position for an index_struct key (no x*6+y factor).
// index_struct = (a*n+b)*9 + (3*i+j), so b = (key/9) % n, j = key % 3.
[[maybe_unused]] static inline int index_struct_to_close(int key, int n) {
    int ab = key / 9;
    int b = ab % n;
    int j = key % 3;
    return 3 * b + j;
}

// C-table key by boundary nucleotide positions only — matches LCDSfold's GetIndex(left_pos, nuc_lo, nuc_ro).
// left_pos  = sigma(a, i) = 3*a+i  (left boundary nucleotide index in sequence)
// right_pos = sigma(b, j) = 3*b+j  (right boundary nucleotide index; == pos during fill)
// nuc_lo    = nucleotides[protein[a]][x][i]  (left boundary nucleotide, 0..3)
// nuc_ro    = nucleotides[protein[b]][y][j]  (right boundary nucleotide, 0..3)
//
// Formula: (left_pos * nuc_len + right_pos) * 16 + nuc_lo * 4 + nuc_ro
// Max entries per tab_c[pos] = nuc_len * 16 (right_pos == pos for all entries in that table).
// For P15421 (nuc_len=234): 234*16 = 3744 — matches LCDSfold's key space exactly.
// Codon variants (x1,x2) producing the same nuc_lo map to the same key and compete via update_derna.
static inline int nuc_key_c_nuc(int left_pos, int right_pos, int nuc_lo, int nuc_ro, int nuc_len) {
    return (left_pos * nuc_len + right_pos) * 16 + nuc_lo * 4 + nuc_ro;
}

// Structural+nucleotide key for SCLeft entries. Two SCLeft entries with the same
// (s_left_start, seg_len_left, nuc_li, nuc_ri, nuc_i1, nuc_pm1) can merge: they
// have identical closure energy at MANNER_SCLefttoC, so per (key, x=svar.x, y=cvar.y)
// we keep the best (mfe, cai). c_inner_left is derived (= s_left_start + seg_len_left).
// Layout (msb→lsb): [13: s_left_start][5: seg_len_left][2: nuc_li][2: nuc_ri][2: nuc_i1][2: nuc_pm1]
// = 26 bits, fits in non-negative int (so it does not collide with negative-key
// schemes elsewhere).
static inline int sc_left_key_encode(int s_left_start, int seg_len_left,
                                     int nuc_li, int nuc_ri, int nuc_i1, int nuc_pm1) {
    return  ((s_left_start & 0x1FFF) << 13)
          | ((seg_len_left & 0x1F) << 8)
          | ((nuc_li & 0x3) << 6)
          | ((nuc_ri & 0x3) << 4)
          | ((nuc_i1 & 0x3) << 2)
          | (nuc_pm1 & 0x3);
}

// Right-bound nucleotide position (sigma(b,j) = pos at fill time) decoded from a nuc_key_c_nuc key.
static inline int nuc_key_c_to_close(int key, int nuc_len) {
    return (key / 16) % nuc_len;
}
// Left-bound nucleotide position decoded from a nuc_key_c_nuc key.
static inline int nuc_key_c_to_left(int key, int nuc_len) {
    return (key / 16) / nuc_len;
}


// #define DEBUG

void fill_position_beam_tables(int n, vector<int>& protein, double lambda,
                              int beamsize, AllTablesDerna& tables) {
    const int nuc_len = 3 * n;
    Zuker z(n, 1, protein);  // for hairpin CAI and energy tables


#ifdef DEBUG
    s_prune_log.open("position_beam_prune.log", std::ios::out);
    s_prune_protein = &protein;
    if (s_prune_log.is_open())
        s_prune_log << "--- position_beam_prune.log (top 20 per prune: score, mfe, cai, a, b, i, j, x, y, n, nucL, nucR, manner, from_keys) ---\n";
#endif
    // Dual-beam diagnostic CSV (env DERNA_DUMP_PRUNES=path).
    if (const char* env_path = std::getenv("DERNA_DUMP_PRUNES")) {
        s_dualbeam_csv.open(env_path, std::ios::out);
        s_prune_protein = &protein;
        s_dualbeam_step_id = 0;
        if (s_dualbeam_csv.is_open()) {
            s_dualbeam_csv << "step,table,pos,aux_len,beamsize,n_cands,"
                              "left_pos,right_pos,a,b,i,j,x,y,key,manner,"
                              "score,mfe,cai,"
                              "rank_combined,rank_mfe,rank_cai,rank_cum,kept_in_topk\n";
            std::cerr << "[dualbeam-diag] dumping prune candidates to " << env_path << "\n";
        }
    }

    auto& tab_n = tables.bestN;
    auto& tab_s = tables.bestS;
    auto& tab_f = tables.bestF;
    auto& tab_c = tables.bestC;
    auto& tab_cs = tables.bestCS;
    auto& tab_sc_left = tables.bestSCLeft;
    auto& tab_m1 = tables.bestM1;
    auto& tab_m2 = tables.bestM2;
    auto& tab_multi = tables.bestMulti;

    tab_n.resize(nuc_len);
    tab_s.resize(nuc_len);
    for (int p = 0; p < nuc_len; ++p) tab_s[p].resize(SINGLE_MAX_LEN + 1);
    tab_f.resize(nuc_len);
    tab_c.resize(nuc_len);
    tab_cs.resize(nuc_len);
    tab_sc_left.resize(nuc_len);
    tab_m1.resize(nuc_len);
    tab_m2.resize(nuc_len);
    tab_multi.resize(nuc_len);

    // ---- Codon-pair aggregation infrastructure (Step 1-2 of closure refactor) ----
    // s_agg[pos][seg_len][key] / c_agg[pos][key] are populated AFTER final pruning
    // at each position. Currently BUILT but not yet CONSUMED — closure logic is
    // unchanged. See codon_agg.h for design.
    CodonAggTables agg_tables;
    agg_tables.resize(nuc_len, SINGLE_MAX_LEN + 1);
    const bool agg_diag_on = agg_diag_enabled();

    // ---------- Base case at position 0 ----------
    int aa0 = 0, ii0 = 0, slot0 = 0;
    int paa0 = protein[aa0];
    const int ncod0 = n_codon[paa0];
    for (int xx0 = 0; xx0 < ncod0; ++xx0) {
        int nuc0 = nucleotides[paa0][xx0][0];
        int key0 = codon_beam_key(0, nuc0, nuc0);
        double mfe0 = 0.0, cai0 = (is_last_nuc(0)) ? codon_cai[paa0][xx0] : 0.0;
        double sc0 = combined_score(lambda, mfe0, cai0);
        BeamEntry ent0(sc0, aa0, aa0, ii0, slot0, xx0, xx0, (int)DernaManner::MANNER_NONEtoN, mfe0, cai0);
        update_derna(tab_n[0], key0, sc0, ent0);
        ent0.backtrace_type = (int)DernaManner::MANNER_NONEtoF;
        update_derna(tab_f[0], key0, sc0, ent0);
        ent0.backtrace_type = (int)DernaManner::MANNER_NONEtoS;
        if (tab_s[0].size() > 1) update_derna(tab_s[0][1], key0, sc0, ent0);
    }

    // best_f_prefix[p][y] = best (min) F score at position p among entries with right codon y.
    // Used by prune_cumulative so entries are ranked by codon-compatible F prefix potential.
    // Matches LCDSfold BeamPrune which checks IsLegal at the F-candidate boundary.
    const double inf_f = std::numeric_limits<double>::infinity();
    const int max_nucs = 4;
    vector<vector<double>> best_f_prefix(nuc_len, vector<double>(max_nucs, inf_f));
    // Lookahead-aware CS prune (Approach a): per-pos min S_left score ending
    // at p. Updated incrementally after the S table at pos p is built each
    // iteration, before the CS prune at pos p reads it.
    vector<double> min_sl_at_pos(nuc_len, std::numeric_limits<double>::infinity());
    // Seed pos 0 from base-case S entries.
    if (nuc_len > 0 && tab_s[0].size() > 1) {
        for (int seg_len = 1; seg_len < (int)tab_s[0].size(); ++seg_len) {
            for (const auto& kv : tab_s[0][seg_len]) {
                if (std::isfinite(kv.second.score) && kv.second.score < min_sl_at_pos[0])
                    min_sl_at_pos[0] = kv.second.score;
            }
        }
    }

    // Selective Block 1 (Approach c): per-pos dropped_tuples set.
    // CS entries dropped by the CS prune at pos N → packed
    // (closed_right, seg_len, s_right_key, c_inner_key) tuple. Block 1 at
    // pos N+1 only enumerates (c_kv, s_rkv) tuples that ARE in this set,
    // since those are the only tuples whose CS coverage was lost.
    //
    // Tuples NOT in this set fall into two categories — both safe to skip:
    //   (a) covered: their CS entry survived, so CS_to_C / S_CS_to_C will close
    //       them at pos+1 (with iv-wrap completing the inner-C variant span).
    //   (b) non-producer: their (c_key, s_key) failed geometry checks in
    //       C+S→CS, so no CS entry was ever produced — Block 1 has no useful
    //       work to do for them either (its own geometry filter would skip).
    //
    // tab_cs_selective_active[N] = whether selective filtering is valid at N.
    vector<std::unordered_set<uint64_t>> tab_cs_dropped_units(nuc_len);
    vector<bool> tab_cs_selective_active(nuc_len, false);
    // Seed position 0 F score from base-case entries.
    for (const auto& kv : tab_f[0]) {
        int nuc_r = nucleotides[protein[kv.second.b]][kv.second.y][kv.second.j];
        best_f_prefix[0][nuc_r] = min(best_f_prefix[0][nuc_r], kv.second.score);
    }

    // ---------- Main loop: position 1 .. nuc_len-1 ----------
    auto t_fill_start = high_resolution_clock::now();
    double ms_N_S = 0, ms_C = 0, ms_Multi = 0, ms_F = 0;
    double ms_CS_to_C = 0, ms_C_S_to_CS = 0, ms_C_to_C = 0, ms_S_C_to_C = 0;  // C sub-phases
    // Transition case counts (candidates examined per state transition) for runtime profiling.
    uint64_t cnt_N_EtoN = 0, cnt_NtoC = 0, cnt_S_EtoS = 0, cnt_CS_to_C = 0, cnt_S_CS_to_C = 0;
    uint64_t cnt_C_S_to_CS = 0, cnt_C_to_C = 0, cnt_S_C_to_C = 0, cnt_Multi_EtoMulti = 0, cnt_MultitoC = 0;
    uint64_t cnt_C_to_M1 = 0, cnt_M1_C_to_M2 = 0, cnt_M1_EtoM1 = 0, cnt_M2_to_M1 = 0;
    uint64_t cnt_F_EtoF = 0, cnt_C_to_F = 0, cnt_F_C_to_F = 0;
    uint64_t cnt_SC_left_build = 0, cnt_SCLeft_to_C = 0;
    uint64_t cnt_blk1_rightbulge = 0, cnt_blk1_internal = 0;
    uint64_t blk1_c_kv_kept_dbg = 0, blk1_c_kv_filtered_dbg = 0;

    // Block 1 dense per-thread task tables. Hoisted to function scope so
    // their data_ vectors are allocated once and reused across positions
    // (reset() zeroes only the slots populated by the previous pos).
#ifdef _OPENMP
    const int num_threads_blk1_outer = omp_get_max_threads();
#else
    const int num_threads_blk1_outer = 1;
#endif
    std::vector<DenseBeamTable> blk1_thread_dense_maps(num_threads_blk1_outer);
    std::vector<std::vector<int>> blk1_thread_key_order(num_threads_blk1_outer);

    // Whether to use the dense-vector path for Block 1's per-thread task_map.
    // Default ON. Set BLK1_DENSE=0 to fall back to the prior DernaBeamMap path
    // (kept in-place for direct A/B comparison). Lossless: both paths produce
    // identical curr_c contributions.
    static const bool blk1_use_dense = []() {
        const char* e = std::getenv("BLK1_DENSE");
        return !(e && e[0] == '0');
    }();

    // BLK1_AGG: per-(idx,x,y) early-rejection in emit_il/emit_rb.
    // Phase 0 instrumentation showed 99.83% of update_derna calls are redundant
    // (existing variant whose score is not improved); skipping the scratch_il
    // build and the update call entirely on those redundant emits saves the
    // bt_info copy + function-call work. Lossless: equivalent to running the
    // update path which would have been a no-op anyway.
    // Default ON (set BLK1_AGG=0 to disable for A/B comparison). Validated
    // gap-neutral on the 19-case K=1000 sweep; ~6% wall-time win on P01707.
    static const bool blk1_agg_enabled = []() {
        const char* e = std::getenv("BLK1_AGG");
        return !(e && e[0] == '0');  // default ON; set BLK1_AGG=0 to disable
    }();

    for (int pos = 1; pos < nuc_len; ++pos) {
        auto t0 = high_resolution_clock::now();
        int bb = pos / 3, slot = pos % 3;
        if (bb >= n) continue;
        int pbb = protein[bb];
        const int ncod_bb = n_codon[pbb];
        int pos_prev = pos - 1;

        auto& curr_n = tab_n[pos];
        auto& curr_s = tab_s[pos];
        auto& curr_f = tab_f[pos];
        auto& curr_c = tab_c[pos];
        auto& curr_cs = tab_cs[pos];
        auto& curr_sc_left = tab_sc_left[pos];
        auto& curr_m1 = tab_m1[pos];
        auto& curr_m2 = tab_m2[pos];
        auto& curr_multi = tab_multi[pos];

        // ---------- NONE -> S, NONE -> N (start new unpaired/single at pos) ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            int nuc_yy = nucleotides[pbb][yy][slot];
            int key = codon_beam_key(pos, nuc_yy, nuc_yy);
            double mfe0 = 0.0, cai0 = (is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0;
            double sc_cand = combined_score(lambda, mfe0, cai0);
            BeamEntry ent_cand(sc_cand, bb, bb, slot, slot, yy, yy, (int)DernaManner::MANNER_NONEtoS, mfe0, cai0);
            if (curr_s.size() > 1) update_derna(curr_s[1], key, sc_cand, ent_cand);
            ent_cand.backtrace_type = (int)DernaManner::MANNER_NONEtoN;
            update_derna(curr_n, key, sc_cand, ent_cand);
        }

        // ---------- N -> N (extend unpaired) and N -> C (hairpin) ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            DernaBeamMap& prev_n = tab_n[pos_prev];
            for (auto& kv : prev_n) {
                int from_key = kv.first;
                BeamEntry& from_ent = kv.second;
                int from_a = from_ent.a, from_b = from_ent.b, from_i = from_ent.i, from_j = from_ent.j;

                bool same_codon = (from_b == bb && slot == from_j + 1);
                bool next_codon = (from_b == bb - 1 && from_j == 2 && slot == 0);
                if (!same_codon && !next_codon) continue;

                for (auto& var : from_ent.variants) {
                    if (same_codon && yy != var.y) continue;
                    cnt_N_EtoN++;

                    // N -> N
                    double mfe_ext = var.mfe, cai_ext = var.cai + ((is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0);
                    double sc_ext = combined_score(lambda, mfe_ext, cai_ext);
                    int nuc_L_n = nucleotides[protein[from_a]][var.x][from_i];
                    int nuc_R_n = nucleotides[pbb][yy][slot];
                    int key_ext = codon_beam_key(sigma(from_a, from_i), nuc_L_n, nuc_R_n);

                    BeamEntry ent_ext(sc_ext, from_a, bb, from_i, slot, var.x, yy, (int)DernaManner::MANNER_N_EtoN, mfe_ext, cai_ext);
                    ent_ext.bt_info = {from_key, var.x, var.y};
                    update_derna(curr_n, key_ext, sc_ext, ent_ext);

                    // N -> C: hairpin with pair (left_bound-1, pos); use Zuker::hairpin_loop + z.add_hairpin_CAI_2
                    int left_bound = sigma(from_a, from_i);
                    int left_nuc = left_bound - 1;
                    if (left_nuc >= 0 && (pos - left_nuc) > HAIRPIN_GAP) {
                        const int loop_len = pos - left_nuc - 1;  // unpaired length
                        if (loop_len < 0) continue;
                        if (loop_len > SINGLE_MAX_LEN) continue;
                        int aa_loop = left_nuc / 3, ii_loop = left_nuc % 3;
                        int paa_loop = protein[aa_loop];
                        const int ncod_loop = n_codon[paa_loop];
                        int nuc_right = nucleotides[pbb][yy][slot];
                        int pna = (aa_loop + 1 < n) ? protein[aa_loop + 1] : -1;
                        int ppb_prev = (bb > 0) ? protein[bb - 1] : -1;
                        int ncod_an = (pna >= 0) ? n_codon[pna] : 0;
                        int ncod_bp = (ppb_prev >= 0) ? n_codon[ppb_prev] : 0;
                        for (int xx_loop = 0; xx_loop < ncod_loop; ++xx_loop) {
                            if (aa_loop == from_a && xx_loop != var.x) continue;
                            int nuc_left = nucleotides[paa_loop][xx_loop][ii_loop];
                            if (BP_pair[nuc_left + 1][nuc_right + 1] == 0) continue;
                            // FIX: include var.cai (loop-interior codon cai already accumulated
                            // by N extensions) plus pair boundary codons only.  The a1/b1 intermediate
                            // codons used by the old add_hairpin_CAI_2 call were already in var.cai,
                            // so passing a1/b1=-1 avoids double-counting.
                            double pair_cai_only = z.add_hairpin_CAI_2(aa_loop, bb, xx_loop, yy,
                                                                -1, -1, -1, -1, ii_loop, slot);
                            double cai_hp = var.cai + pair_cai_only;
                            double best_total = inf, best_mfe = 0;
                            double best_cai = cai_hp;  // standard path: constant across x1_cur/y1_cur
                            const int n_left = (ii_loop == 2 && ncod_an > 0) ? ncod_an : 1;
                            const int n_right = (slot == 0 && ncod_bp > 0) ? ncod_bp : 1;
                            // Loop-length in Zuker's l convention = pos - left_nuc = loop_len + 1.
                            const int l_zuker = loop_len + 1;
                            const bool try_special = (l_zuker == 4 || l_zuker == 5 || l_zuker == 7 || l_zuker == 9);
                            for (int x1_cur = 0; x1_cur < n_left; ++x1_cur) {
                                // If aa_loop+1's codon is already determined by the N
                                // entry's variant, x1_cur must match that codon.
                                if (ii_loop == 2 && pna >= 0) {
                                    int aa_next = aa_loop + 1;
                                    if (aa_next == from_a && x1_cur != var.x) continue;
                                    if (aa_next == from_b && x1_cur != var.y) continue;
                                }
                                for (int y1_cur = 0; y1_cur < n_right; ++y1_cur) {
                                    // If (bb-1)'s codon is already determined by the N
                                    // entry's variant, y1_cur must match that codon.
                                    if (slot == 0 && ppb_prev >= 0) {
                                        int aa_prev = bb - 1;
                                        if (aa_prev == from_a && y1_cur != var.x) continue;
                                        if (aa_prev == from_b && y1_cur != var.y) continue;
                                    }
                                    cnt_NtoC++;
                                    int xi_ = (ii_loop < 2) ? nucleotides[paa_loop][xx_loop][ii_loop + 1]
                                            : nucleotides[pna][x1_cur][0];
                                    int _yj = (slot > 0) ? nucleotides[pbb][yy][slot - 1]
                                            : nucleotides[ppb_prev][y1_cur][2];
                                    // Special-hairpin (triloop/tetraloop/hexaloop) path: Turner model mandates
                                    // special energy whenever the loop sequence matches a special key — it is
                                    // NOT a competing alternative to the regular formula. If z.hairpin_special_CAI_pub
                                    // returns he_sp < inf, the loop string matched a special entry; use that
                                    // energy directly. Otherwise fall through to the regular formula.
                                    // (Prior code took min(regular, special), which produced a phantom energy
                                    // when the regular formula happened to be more favorable than the special
                                    // — evaluate_structure_energy always uses the special for matching sequences,
                                    // so the beam's NtoC.mfe disagreed with the traceback structure energy.)
                                    double mfe_hp;
                                    double cai_cand = cai_hp;
                                    bool matched_special = false;
                                    if (try_special) {
                                        // When pna == ppb_prev (the middle-loop aa is a SINGLE
                                        // amino acid shared by both x1 and y1 parameters), force
                                        // x1 to match the iterated y1_cur so hairpin_special_CAI's
                                        // tetraloop / hexaloop codon lookup sees a consistent
                                        // codon choice. Without this, n_left==1 pins x1 to 0
                                        // and the special test misses matches such as CUACGG
                                        // whenever aa=15's codon is not codon 0.
                                        int x1_eff = x1_cur;
                                        if (pna >= 0 && pna == ppb_prev) x1_eff = y1_cur;
                                        double he_sp, mfe_sp_l, cai_sp_l;
                                        std::vector<int> tmp_v;
                                        std::tie(he_sp, mfe_sp_l, cai_sp_l, tmp_v) = z.hairpin_special_CAI_pub(
                                            lambda, l_zuker, aa_loop, bb, paa_loop,
                                            pbb, pna, ppb_prev, x1_eff, y1_cur,
                                            nuc_left, xi_, _yj, nuc_right, xx_loop, yy,
                                            ii_loop, slot);
                                        if (he_sp < inf) {
                                            // Convert back: hairpin_special_CAI multiplies mfe by lambda and
                                            // cai by (lambda-1).
                                            double mfe_sp_raw = (lambda != 0.0) ? (mfe_sp_l / lambda) : mfe_sp_l;
                                            double pair_cai_only_sp = z.add_hairpin_CAI_2(aa_loop, bb, xx_loop, yy,
                                                                                    -1, -1, -1, -1, ii_loop, slot);
                                            cai_cand = var.cai + pair_cai_only_sp;
                                            mfe_hp = mfe_sp_raw;
                                            matched_special = true;
                                        }
                                    }
                                    if (!matched_special) {
                                        mfe_hp = (double)Zuker::hairpin_loop(nuc_left, nuc_right, xi_, _yj, loop_len);
                                    }
                                    double total_e = combined_score(lambda, mfe_hp, cai_cand);
                                    if (total_e < best_total) {
                                        best_total = total_e;
                                        best_mfe = mfe_hp;
                                        best_cai = cai_cand;
                                    }
                                }
                            }
                            if (best_total >= inf) continue;
                            int key_cl = nuc_key_c_nuc(sigma(aa_loop, ii_loop), pos, nuc_left, nuc_right, nuc_len);
                            BeamEntry ent_cl(best_total, aa_loop, bb, ii_loop, slot, xx_loop, yy, (int)DernaManner::MANNER_NtoC, best_mfe, best_cai);
                            // Record the N-chain variant (var.x, var.y) whose (cai, mfe) chain produced
                            // this hairpin's (best_mfe, best_cai) so traceback selects the matching variant.
                            ent_cl.bt_info = {from_key, var.x, var.y};
                            check_c_entry_invariant(lambda, ent_cl, "N->C", pos);
                            update_derna(curr_c, key_cl, best_total, ent_cl);
                        }
                    }
                }
            }
        }
        // N prune: cumulative (bestF[left-1] + local), mirrors LCDSfold BeamPrune(false).
        prune_cumulative(curr_n, beamsize, pos, n, DernaTableKind::N, best_f_prefix, protein);

        // ---------- S -> S (extend single-stranded) ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            for (int seg_len = 1; seg_len <= min(SINGLE_MAX_LEN - 1, pos); ++seg_len) {
                if (seg_len + 1 >= (int)curr_s.size()) continue;
                DernaBeamMap& prev_s = tab_s[pos_prev][seg_len];
                for (auto& kv : prev_s) {
                    BeamEntry& from_ent = kv.second;
                    int from_a = from_ent.a, from_b = from_ent.b, from_i = from_ent.i, from_j = from_ent.j;
                    bool same_codon = (from_b == bb && slot == from_j + 1);
                    bool next_codon = (from_b == bb - 1 && from_j == 2 && slot == 0);
                    if (!same_codon && !next_codon) continue;

                    for (auto& var : from_ent.variants) {
                        if (same_codon && yy != var.y) continue;
                        cnt_S_EtoS++;
                        double mfe_ext = var.mfe, cai_ext = var.cai + ((is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0);
                        double sc_ext = combined_score(lambda, mfe_ext, cai_ext);
                        int nuc_L_s = nucleotides[protein[from_a]][var.x][from_i];
                        int nuc_R_s = nucleotides[pbb][yy][slot];
                        int key_ext = codon_beam_key(sigma(from_a, from_i), nuc_L_s, nuc_R_s);
                        BeamEntry ent_ext(sc_ext, from_a, bb, from_i, slot, var.x, yy, (int)DernaManner::MANNER_S_EtoS, mfe_ext, cai_ext);
                        ent_ext.bt_info = {kv.first, var.x, var.y};
                        update_derna(curr_s[seg_len + 1], key_ext, sc_ext, ent_ext);
                    }
                }
            }
        }
        for (int seg_len = 1; seg_len < (int)curr_s.size(); ++seg_len)
            prune_beam_derna_checked(curr_s[seg_len], beamsize, pos, n, DernaTableKind::S, seg_len);

        // ---- Build codon-pair aggregation for S table at this pos (post-prune) ----
        // S table is finalized for pos after the prune above. Block 1 / closure
        // logic still uses tab_s directly; agg table is BUILT but not CONSUMED.
        {
            CodonAggTables::PerPosDiag* s_diag_ptr = agg_diag_on ? &agg_tables.diag[pos] : nullptr;
            for (int seg_len = 1; seg_len < (int)curr_s.size(); ++seg_len) {
                build_s_codon_agg(curr_s[seg_len],
                                  agg_tables.s_agg[pos][seg_len],
                                  s_diag_ptr);
            }
        }

        // Update min_sl_at_pos[pos]: best (min) S_left score ending at this pos
        // across all seg_lens and all keys. Used by Approach (a) lookahead prune.
        {
            double best_sl = std::numeric_limits<double>::infinity();
            for (int seg_len = 1; seg_len < (int)curr_s.size(); ++seg_len) {
                for (const auto& kv : curr_s[seg_len]) {
                    if (std::isfinite(kv.second.score) && kv.second.score < best_sl)
                        best_sl = kv.second.score;
                }
            }
            min_sl_at_pos[pos] = best_sl;
        }

        // ---------- Local flattened caches for hot C/CS recurrences ----------
        // Keep semantics identical; only avoid repeatedly walking unordered_map buckets
        // inside the hottest nested loops.
        vector<vector<pair<int, const BeamEntry*>>> curr_s_flat(curr_s.size());
        for (int seg_len = 1; seg_len < (int)curr_s.size(); ++seg_len) {
            auto& flat = curr_s_flat[seg_len];
            flat.reserve(curr_s[seg_len].size());
            for (auto& kv : curr_s[seg_len]) flat.push_back({kv.first, &kv.second});
        }

        vector<pair<int, const BeamEntry*>> prev_cs_flat;
        prev_cs_flat.reserve(tab_cs[pos_prev].size());
        for (auto& kv : tab_cs[pos_prev]) prev_cs_flat.push_back({kv.first, &kv.second});

        vector<pair<int, const BeamEntry*>> prev_sc_left_flat;
        prev_sc_left_flat.reserve(tab_sc_left[pos_prev].size());
        for (auto& kv : tab_sc_left[pos_prev]) prev_sc_left_flat.push_back({kv.first, &kv.second});

        unordered_map<int, vector<pair<int, const BeamEntry*>>> c_flat_cache;
        auto get_c_flat = [&](int right_pos) -> const vector<pair<int, const BeamEntry*>>& {
            auto it = c_flat_cache.find(right_pos);
            if (it != c_flat_cache.end()) return it->second;
            auto& flat = c_flat_cache[right_pos];
            if (right_pos >= 0 && right_pos < nuc_len) {
                flat.reserve(tab_c[right_pos].size());
                for (auto& kv : tab_c[right_pos]) flat.push_back({kv.first, &kv.second});
                std::sort(flat.begin(), flat.end(),
                          [](const auto& a, const auto& b){ return a.first < b.first; });
            }
            return flat;
        };

        unordered_map<int, vector<pair<int, const BeamEntry*>>> s_flat_cache;
        auto get_s_flat = [&](int end_pos, int seg_len) -> const vector<pair<int, const BeamEntry*>>& {
            const int cache_key = end_pos * (SINGLE_MAX_LEN + 1) + seg_len;
            auto it = s_flat_cache.find(cache_key);
            if (it != s_flat_cache.end()) return it->second;
            auto& flat = s_flat_cache[cache_key];
            if (end_pos >= 0 && end_pos < nuc_len && seg_len >= 0 && seg_len < (int)tab_s[end_pos].size()) {
                flat.reserve(tab_s[end_pos][seg_len].size());
                for (auto& kv : tab_s[end_pos][seg_len]) flat.push_back({kv.first, &kv.second});
                std::sort(flat.begin(), flat.end(),
                          [](const auto& a, const auto& b){ return a.first < b.first; });
            }
            return flat;
        };

        ms_N_S += duration<double, milli>(high_resolution_clock::now() - t0).count();

        t0 = high_resolution_clock::now();
        // ---------- CS -> C (close CS with an outer pair; allow bulge/internal) ----------
        {
        auto t_c0 = high_resolution_clock::now();
        // This mirrors LCDSfoldCAI_DN_beam:
        //   (1) CS -> C : wrap an outer pair around an existing (C + right-single) construction
        //   (2) S + CS -> C : prepend an additional left-single segment before wrapping
        // NOTE: We approximate the general single-loop energy using Zuker::bulge_loop when one side is a bulge.
        // If you have a full internal-loop scorer, swap it in here.
        // Option 4 hoist: per-position constants invariant across all loops here.
        {
            const bool pos_is_last_cs = is_last_nuc(pos);
            // Precompute yy array: nuc_ro per yy, and codon_cai per yy.
            int cs_yy_nuc_ro[6];
            double cs_yy_cai[6];
            for (int yy = 0; yy < ncod_bb; ++yy) {
                cs_yy_nuc_ro[yy] = nucleotides[pbb][yy][slot];
                cs_yy_cai[yy] = pos_is_last_cs ? codon_cai[pbb][yy] : 0.0;
            }

            // ---- Parallelization setup ----
            // Each entry in prev_cs_flat defines an independent writer into curr_c.
            // We buffer writes per-thread and deterministically merge afterwards, mirroring
            // the pattern used in Block 1 below.
            const int num_cs_tasks = (int)prev_cs_flat.size();
#ifdef _OPENMP
            const int num_threads_cs = omp_get_max_threads();
#else
            const int num_threads_cs = 1;
#endif
            // Pre-warm s_flat_cache for every (seam_pos, seg_len_left) that the CS->C body
            // might touch. The cache is a std::unordered_map that lazily inserts on access
            // via get_s_flat(); inside the parallel region that would be a data race.
            for (const auto& kv_pw : prev_cs_flat) {
                const BeamEntry& cs_ent_pw = *kv_pw.second;
                int seam_pos_pw = cs_ent_pw.cs_inner_left - 1;
                if (seam_pos_pw < 0 || seam_pos_pw >= nuc_len) continue;
                int max_sl_pw = min(SINGLE_MAX_LEN, seam_pos_pw + 1);
                for (int sl = 1; sl <= max_sl_pw; ++sl) {
                    (void)get_s_flat(seam_pos_pw, sl);
                }
            }

            std::vector<DernaBeamMap> cs_thread_maps(num_threads_cs);
            std::vector<std::vector<int>> cs_thread_key_order(num_threads_cs);
            std::vector<uint64_t> cnt_CS_to_C_local(num_threads_cs, 0);
            std::vector<uint64_t> cnt_S_CS_to_C_local(num_threads_cs, 0);

            // Swapped outer loops: iterate prev_cs_flat outer, yy inner. Per-kv setup runs once.
            #pragma omp parallel for schedule(dynamic, 1)
            for (int cs_task_idx = 0; cs_task_idx < num_cs_tasks; ++cs_task_idx) {
                const auto& kv = prev_cs_flat[cs_task_idx];
                int cs_tab_k = kv.first;
                const BeamEntry& cs_ent = *kv.second;
#ifdef _OPENMP
                const int tid_cs = omp_get_thread_num();
#else
                const int tid_cs = 0;
#endif
                DernaBeamMap& task_map_cs = cs_thread_maps[tid_cs];
                auto& task_order_cs = cs_thread_key_order[tid_cs];

                bool same_codon_cs = false;
                if (!same_or_next_codon(cs_ent.b, cs_ent.j, bb, slot, same_codon_cs)) continue;

                int single_start = cs_ent.cs_single_start;
                int closed_right = single_start - 1;
                if (single_start < 0 || closed_right < 0 || closed_right >= nuc_len) continue;

                int seg_len_right = cs_ent.cs_right_len;
                if (seg_len_right <= 0 || seg_len_right >= (int)tab_s[pos_prev].size()) continue;

                int inner_left = cs_ent.cs_inner_left;
                int outer_left_nuc = inner_left - 1;
                if (outer_left_nuc < 0) continue;

                int aa_o = outer_left_nuc / 3;
                int ii_o = outer_left_nuc % 3;
                int paa_o = protein[aa_o];
                const int ncod_o = n_codon[paa_o];

                int inner_key = cs_ent.cs_inner_key;
                int s_right_key = cs_ent.cs_right_s_key;
                if (inner_key < 0 || s_right_key < 0) continue;
                if (tab_c[closed_right].count(inner_key) == 0) continue;
                if (tab_s[pos_prev][seg_len_right].count(s_right_key) == 0) continue;

                const BeamEntry& inner_ent = tab_c[closed_right].at(inner_key);
                const BeamEntry& s_right_ent = tab_s[pos_prev][seg_len_right].at(s_right_key);

                // Decode CS key to recover the committed nuc_inner_L and nuc_single_start.
                // cs_get_index packs the low 6 bits as
                //   nuc_inner_L * 16 + nuc_outer_R * 4 + nuc_single_start,
                // so triple = key % 64 (NOT key / 31). The previous decoder divided
                // by SINGLE_MAX_LEN+1 first, which yielded near-random values in
                // [0,3] — Block 1 was masking the resulting CS-path errors.
                int inner_nuc_li = -1;
                int nuc_qp1_committed = -1;
                {
                    long long csk = (long long)cs_ent.cs_pack_outer;
                    int triple = (int)(csk % 64);
                    int nuc_inner_L = triple / 16;
                    int nuc_single_start = triple % 4;
                    inner_nuc_li = nuc_inner_L;
                    nuc_qp1_committed = nuc_single_start;
                }

                // Coverage-gap fix: collect inner-C variants matching the CS's committed
                // nuc_inner_L. The yy-loop body below is wrapped in a per-iv loop so each
                // inner codon assignment gets its own seam-compatible xx_o (and slv.y in
                // S_CS path).
                //
                // Critical correctness point (2026-04-28): the previous fixed-cap of 6
                // dropped (x,y) variants when inner_ent had >6 matching variants — the
                // iteration order placed all-y=Y0 variants first, evicting better-scoring
                // (x, Y1) tuples. For c states with multiple x values per nuc_inner_L AND
                // multiple y values per nuc_outer_R (e.g. Leucine inner pair, where ALL
                // 6 codons produce U at slot 1), 12 valid variants existed but only 6
                // worse-scored ones were sampled.
                //
                // Selection rule:
                //   * For each iv.x, pick the BEST-score variant (any iv.y producing the
                //     same nuc_ri — which is invariant since nuc_ri is encoded in the
                //     c_ent key).
                //   * When inner_ent.b == bb (closure's outer-right shares a codon with
                //     inner-C's right), iv.y also locks yy via line 1756, so we must keep
                //     ALL (iv.x, iv.y) ivs (max 36) to expose every valid yy choice.
                std::array<std::array<const XYVariant*, 6>, 6> best_iv_by_xy;
                std::array<std::array<double, 6>, 6> best_iv_score_by_xy;
                for (auto& row : best_iv_by_xy) for (auto& p : row) p = nullptr;
                for (auto& row : best_iv_score_by_xy) for (auto& s : row) s = 1e30;
                for (const auto& iv : inner_ent.variants) {
                    if (iv.x < 0 || iv.x >= 6 || iv.y < 0 || iv.y >= 6) continue;
                    int nl = nucleotides[protein[inner_ent.a]][iv.x][inner_ent.i];
                    if (nl != inner_nuc_li) continue;
                    if (iv.score < best_iv_score_by_xy[iv.x][iv.y]) {
                        best_iv_score_by_xy[iv.x][iv.y] = iv.score;
                        best_iv_by_xy[iv.x][iv.y] = &iv;
                    }
                }
                std::array<const XYVariant*, 36> matching_ivs;
                int n_matching_ivs = 0;
                const bool inner_b_eq_bb = (inner_ent.b == bb);
                if (inner_b_eq_bb) {
                    // Keep every (x,y) — line 1756 locks yy to inner_y_v.
                    for (int xv = 0; xv < 6; ++xv) for (int yv = 0; yv < 6; ++yv) {
                        if (best_iv_by_xy[xv][yv] && n_matching_ivs < 36)
                            matching_ivs[n_matching_ivs++] = best_iv_by_xy[xv][yv];
                    }
                } else {
                    // Best y per x — yy is free of iv.y, so picking best-score iv per x is lossless.
                    for (int xv = 0; xv < 6; ++xv) {
                        const XYVariant* best = nullptr;
                        double best_sc = 1e30;
                        for (int yv = 0; yv < 6; ++yv) {
                            if (best_iv_by_xy[xv][yv] && best_iv_score_by_xy[xv][yv] < best_sc) {
                                best = best_iv_by_xy[xv][yv];
                                best_sc = best_iv_score_by_xy[xv][yv];
                            }
                        }
                        if (best && n_matching_ivs < 36) matching_ivs[n_matching_ivs++] = best;
                    }
                }
                if (n_matching_ivs == 0) continue;

                // --- Variant-elimination optimization for CS->C ---
                // All variants of a merged entry share the same boundary nucleotides (encoded
                // in the key). Therefore loop_e is constant across variant combinations for a
                // given entry combination, and variant scores are separable.

                // (A) Pre-build best-variant-by-codon lookup tables for S_left entries.
                struct BestVar { double score; const XYVariant* ptr; };
                static constexpr double BV_INF_CS = 1e30;
                struct SLCSTask {
                    const BeamEntry* s_left_ent;
                    int s_left_key;
                    int seg_len_left;
                    int outer_left2;
                    int aa_o2, ii_o2, paa_o2, ncod_o2;
                    bool xx_o2_locked;
                    bool same_codon_sl;
                    int nuc_i1_sl, nuc_pm1_sl;  // from key, constant across variants
                    std::array<BestVar, 6> sl_by_x;
                    std::array<BestVar, 6> sl_by_y;
                    std::array<std::array<BestVar, 6>, 6> sl_by_xy;
                    BestVar sl_all;
                };
                std::vector<SLCSTask> sl_cs_tasks;
                {
                    int seam_pos = cs_ent.cs_inner_left - 1;
                    if (seam_pos >= 0 && seam_pos < nuc_len) {
                        int max_sl = min(SINGLE_MAX_LEN, seam_pos + 1);
                        for (int seg_len_left = 1; seg_len_left <= max_sl; ++seg_len_left) {
                            const vector<pair<int, const BeamEntry*>>* s_left_scan = &get_s_flat(seam_pos, seg_len_left);
                            if (s_left_scan->empty()) continue;
                            for (const auto& lkv : *s_left_scan) {
                                const BeamEntry& s_left_ent = *lkv.second;
                                if (sigma(s_left_ent.b, s_left_ent.j) != seam_pos) continue;
                                int left_start = sigma(s_left_ent.a, s_left_ent.i);
                                int outer_left2 = left_start - 1;
                                if (outer_left2 < 0) continue;
                                int aa_o2 = outer_left2 / 3;
                                int ii_o2 = outer_left2 % 3;
                                int paa_o2 = protein[aa_o2];
                                const int ncod_o2 = n_codon[paa_o2];
                                bool same_codon_sl = false;
                                (void)same_or_next_codon(s_left_ent.b, s_left_ent.j, inner_ent.a, inner_ent.i, same_codon_sl);
                                const bool xx_o2_locked = (aa_o2 == s_left_ent.a);
                                const int sl_key = lkv.first;
                                SLCSTask slt_new;
                                slt_new.s_left_ent = &s_left_ent;
                                slt_new.s_left_key = sl_key;
                                slt_new.seg_len_left = seg_len_left;
                                slt_new.outer_left2 = outer_left2;
                                slt_new.aa_o2 = aa_o2; slt_new.ii_o2 = ii_o2;
                                slt_new.paa_o2 = paa_o2; slt_new.ncod_o2 = ncod_o2;
                                slt_new.xx_o2_locked = xx_o2_locked;
                                slt_new.same_codon_sl = same_codon_sl;
                                slt_new.nuc_i1_sl = (sl_key & 0xF) >> 2;
                                slt_new.nuc_pm1_sl = sl_key & 3;
                                for (auto& bv : slt_new.sl_by_x) bv = {BV_INF_CS, nullptr};
                                for (auto& bv : slt_new.sl_by_y) bv = {BV_INF_CS, nullptr};
                                for (auto& row : slt_new.sl_by_xy) for (auto& bv : row) bv = {BV_INF_CS, nullptr};
                                slt_new.sl_all = {BV_INF_CS, nullptr};
                                for (const auto& sv : s_left_ent.variants) {
                                    if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                    if (sv.x >= 0 && sv.x < 6 && sv.score < slt_new.sl_by_x[sv.x].score)
                                        slt_new.sl_by_x[sv.x] = {sv.score, &sv};
                                    if (sv.y >= 0 && sv.y < 6 && sv.score < slt_new.sl_by_y[sv.y].score)
                                        slt_new.sl_by_y[sv.y] = {sv.score, &sv};
                                    if (sv.x >= 0 && sv.x < 6 && sv.y >= 0 && sv.y < 6 &&
                                        sv.score < slt_new.sl_by_xy[sv.x][sv.y].score)
                                        slt_new.sl_by_xy[sv.x][sv.y] = {sv.score, &sv};
                                    if (sv.score < slt_new.sl_all.score)
                                        slt_new.sl_all = {sv.score, &sv};
                                }
                                sl_cs_tasks.push_back(std::move(slt_new));
                            }
                        }
                    }
                }

                // (B) Pre-build best csvar+sr lookup by csvar.y.
                // For each valid csvar.y, find the best S_right variant matching (nuc_qp1_committed, csvar.y),
                // and store the combined sr score + nuc_jm1.
                struct CSBestByY {
                    double sr_mfe, sr_cai, sr_combined;
                    int nuc_jm1;
                    bool valid;
                };
                std::array<CSBestByY, 6> cs_best_by_y;
                for (auto& e : cs_best_by_y) e = {0, 0, BV_INF_CS, -1, false};
                {
                    // First find best srvar for each csvar.y value.
                    // A csvar.y is valid if (a) some csvar has that y, and (b) some srvar matches.
                    std::array<bool, 6> csvar_y_exists{};
                    for (const auto& csvar : cs_ent.variants) {
                        if (!std::isfinite(csvar.cai) || !std::isfinite(csvar.mfe) || cai_looks_garbage(csvar.cai)) continue;
                        if (csvar.y >= 0 && csvar.y < 6) csvar_y_exists[csvar.y] = true;
                    }
                    for (int cy = 0; cy < 6; ++cy) {
                        if (!csvar_y_exists[cy]) continue;
                        const XYVariant* best_sr = nullptr;
                        for (const auto& sv : s_right_ent.variants) {
                            int nx = nucleotides[protein[s_right_ent.a]][sv.x][s_right_ent.i];
                            if (nx != nuc_qp1_committed) continue;
                            if (sv.y != cy) continue;
                            if (!best_sr || sv.score < best_sr->score) best_sr = &sv;
                        }
                        if (!best_sr) continue;
                        int njm1 = nucleotides[protein[cs_ent.b]][cy][cs_ent.j];
                        cs_best_by_y[cy] = {best_sr->mfe, best_sr->cai,
                                            combined_score(lambda, best_sr->mfe, best_sr->cai),
                                            njm1, true};
                    }
                }

                // (C) Also build best csvar+sr by nuc_jm1 (for !same_codon_cs path where
                // csvar.y is not locked to yy).
                struct CSBestByJm1 {
                    double sr_mfe, sr_cai, sr_combined;
                    int csvar_y;
                    bool valid;
                };
                std::array<CSBestByJm1, 4> cs_best_by_jm1;
                for (auto& e : cs_best_by_jm1) e = {0, 0, BV_INF_CS, -1, false};
                for (int cy = 0; cy < 6; ++cy) {
                    if (!cs_best_by_y[cy].valid) continue;
                    int jm1 = cs_best_by_y[cy].nuc_jm1;
                    if (cs_best_by_y[cy].sr_combined < cs_best_by_jm1[jm1].sr_combined) {
                        cs_best_by_jm1[jm1] = {cs_best_by_y[cy].sr_mfe, cs_best_by_y[cy].sr_cai,
                                                cs_best_by_y[cy].sr_combined, cy, true};
                    }
                }

                // (D) Iterate (yy, xx_o) for CStoC and (yy, slt, xx_o2) for S_CStoC.
                // loop_e depends only on key-derived boundary nucs, so it's computed once
                // per (outer pair, nuc_jm1) combination. Variant scores are looked up via tables.
                //
                // iv-wrap: per matching inner-C variant, run the yy iteration. Precomputations
                // (sl_cs_tasks, cs_best_by_y, cs_best_by_jm1) are iv-independent.
                //
                // CSPATH_LOSS_DIAG instrumentation (env-gated): replace the single-best
                // aggregations cs_best_by_y / cs_best_by_jm1 / sl_by_* with full-enumeration
                // closures, to test whether they cause coverage loss when Block 1 is stubbed.
                // Default behavior preserved when the env vars are unset.
                static const bool cspath_cy_full = []() {
                    const char* e = std::getenv("CSPATH_CY_FULL");
                    return e && e[0] == '1';
                }();
                static const bool cspath_jm1_full = []() {
                    const char* e = std::getenv("CSPATH_JM1_FULL");
                    return e && e[0] == '1';
                }();
                static const bool cspath_sl_full = []() {
                    const char* e = std::getenv("CSPATH_SL_FULL");
                    return e && e[0] == '1';
                }();
                for (int iv_idx = 0; iv_idx < n_matching_ivs; ++iv_idx) {
                    const XYVariant& cur_iv = *matching_ivs[iv_idx];
                    const int    inner_x_v   = cur_iv.x;
                    const int    inner_y_v   = cur_iv.y;
                    const double inner_mfe_v = cur_iv.mfe;
                    const double inner_cai_v = cur_iv.cai;
                    const int    nuc_li = nucleotides[protein[inner_ent.a]][inner_x_v][inner_ent.i];
                    const int    nuc_ri = nucleotides[protein[inner_ent.b]][inner_y_v][inner_ent.j];

                for (int yy = 0; yy < ncod_bb; ++yy) {
                    int nuc_ro = cs_yy_nuc_ro[yy];
                    double cai_yy_pbb = cs_yy_cai[yy];

                    // Determine yy_eff and effective values based on same_codon_cs.
                    // When same_codon_cs: yy is locked to csvar.y, so we only proceed if
                    // cs_best_by_y[yy] is valid.
                    // When same_codon_cs && inner_ent.b == bb: yy_eff is locked to inner_y_v,
                    // and csvar.y must == inner_y_v, so only yy == inner_y_v works.
                    if (same_codon_cs) {
                        if (inner_ent.b == bb && yy != inner_y_v) continue;
                        if (!cs_best_by_y[yy].valid) continue;
                    }

                    const int yy_eff = (same_codon_cs && inner_ent.b == bb) ? inner_y_v : yy;
                    const int nuc_ro_eff = same_codon_cs ? nucleotides[pbb][yy_eff][slot] : nuc_ro;
                    const double cai_yy_eff = same_codon_cs
                        ? (pos_is_last_cs ? codon_cai[pbb][yy_eff] : 0.0)
                        : cai_yy_pbb;

                    // -------------------- (1) CS -> C --------------------
                    for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                        if (aa_o == inner_ent.a && xx_o != inner_x_v) continue;
                        int nuc_lo = nucleotides[paa_o][xx_o][ii_o];
                        if (BP_pair[nuc_lo + 1][nuc_ro_eff + 1] == 0) continue;

                        // Find best csvar+sr combination. loop_e depends only on
                        // nuc_jm1 (from csvar.y), which takes at most 4 distinct values.
                        double best_mfe_cand = 0, best_cai_cand = 0;
                        double best_sc = BV_INF_CS;
                        if (same_codon_cs) {
                            if (cspath_cy_full) {
                                // FULL ENUMERATION: iterate all srs matching (nuc_qp1, sv.y==yy).
                                // For each, compute loop_e (with nuc_jm1 derived from yy),
                                // and track best.
                                int njm1_locked = nucleotides[protein[cs_ent.b]][yy][cs_ent.j];
                                double loop_e = score_single_loop_with_mismatch(lambda,
                                    nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                    0, seg_len_right,
                                    nuc_li, njm1_locked, nuc_lo, nuc_qp1_committed);
                                for (const auto& sv : s_right_ent.variants) {
                                    if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                    int nx = nucleotides[protein[s_right_ent.a]][sv.x][s_right_ent.i];
                                    if (nx != nuc_qp1_committed) continue;
                                    if (sv.y != yy) continue;
                                    double mfe_c = inner_mfe_v + sv.mfe + loop_e;
                                    double cai_c = inner_cai_v + sv.cai + cai_yy_eff
                                        + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                                    double sc = combined_score(lambda, mfe_c, cai_c);
                                    if (sc < best_sc) {
                                        best_sc = sc; best_mfe_cand = mfe_c; best_cai_cand = cai_c;
                                    }
                                }
                            } else {
                                // yy is locked to csvar.y → single sr+loop_e computation.
                                if (cs_best_by_y[yy].valid) {
                                    double loop_e = score_single_loop_with_mismatch(lambda,
                                        nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                        0, seg_len_right,
                                        nuc_li, cs_best_by_y[yy].nuc_jm1, nuc_lo, nuc_qp1_committed);
                                    best_mfe_cand = inner_mfe_v + cs_best_by_y[yy].sr_mfe + loop_e;
                                    best_cai_cand = inner_cai_v + cs_best_by_y[yy].sr_cai + cai_yy_eff
                                        + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                                    best_sc = combined_score(lambda, best_mfe_cand, best_cai_cand);
                                }
                            }
                        } else {
                            if (cspath_jm1_full) {
                                // FULL ENUMERATION: iterate all (cy, sr) pairs matching constraints.
                                // For each, jm1 is derived from cy. nuc_jm1 takes 4 values total
                                // but multiple (cy, sr) pairs can map to each.
                                std::array<bool, 6> csvar_y_exists{};
                                for (const auto& csvar : cs_ent.variants) {
                                    if (!std::isfinite(csvar.cai) || !std::isfinite(csvar.mfe) || cai_looks_garbage(csvar.cai)) continue;
                                    if (csvar.y >= 0 && csvar.y < 6) csvar_y_exists[csvar.y] = true;
                                }
                                for (int cy = 0; cy < 6; ++cy) {
                                    if (!csvar_y_exists[cy]) continue;
                                    int njm1 = nucleotides[protein[cs_ent.b]][cy][cs_ent.j];
                                    double loop_e = score_single_loop_with_mismatch(lambda,
                                        nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                        0, seg_len_right,
                                        nuc_li, njm1, nuc_lo, nuc_qp1_committed);
                                    for (const auto& sv : s_right_ent.variants) {
                                        if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                        int nx = nucleotides[protein[s_right_ent.a]][sv.x][s_right_ent.i];
                                        if (nx != nuc_qp1_committed) continue;
                                        if (sv.y != cy) continue;
                                        double mfe_c = inner_mfe_v + sv.mfe + loop_e;
                                        double cai_c = inner_cai_v + sv.cai + cai_yy_eff
                                            + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                                        double sc = combined_score(lambda, mfe_c, cai_c);
                                        if (sc < best_sc) {
                                            best_sc = sc; best_mfe_cand = mfe_c; best_cai_cand = cai_c;
                                        }
                                    }
                                }
                            } else {
                                // Try all distinct nuc_jm1 values, keep best total.
                                for (int jm1 = 0; jm1 < 4; ++jm1) {
                                    if (!cs_best_by_jm1[jm1].valid) continue;
                                    double loop_e = score_single_loop_with_mismatch(lambda,
                                        nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                        0, seg_len_right,
                                        nuc_li, jm1, nuc_lo, nuc_qp1_committed);
                                    double mfe_c = inner_mfe_v + cs_best_by_jm1[jm1].sr_mfe + loop_e;
                                    double cai_c = inner_cai_v + cs_best_by_jm1[jm1].sr_cai + cai_yy_eff
                                        + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                                    double sc = combined_score(lambda, mfe_c, cai_c);
                                    if (sc < best_sc) {
                                        best_sc = sc; best_mfe_cand = mfe_c; best_cai_cand = cai_c;
                                    }
                                }
                            }
                        }
                        if (best_sc >= BV_INF_CS) continue;
                        cnt_CS_to_C_local[tid_cs]++;

                        int key_cl = nuc_key_c_nuc(sigma(aa_o, ii_o), pos, nuc_lo, nuc_ro_eff, nuc_len);

                        BeamEntry ent_cl(best_sc, aa_o, bb, ii_o, slot, xx_o, yy_eff,
                                         (int)DernaManner::MANNER_CStoC, best_mfe_cand, best_cai_cand);
                        ent_cl.bt_info.clear();
                        ent_cl.cs_pack_outer = (long long)cs_tab_k;
                        check_c_entry_invariant(lambda, ent_cl, "CS->C", pos);
                        if (update_derna(task_map_cs, key_cl, best_sc, ent_cl)) task_order_cs.push_back(key_cl);
                    }

                    // -------------------- (2) S + CS -> C --------------------
                    for (const SLCSTask& slt : sl_cs_tasks) {
                        for (int xx_o2 = 0; xx_o2 < slt.ncod_o2; ++xx_o2) {
                            int nuc_lo2 = nucleotides[slt.paa_o2][xx_o2][slt.ii_o2];
                            if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;

                            // For each nuc_jm1 value, compute loop_e once, then find
                            // the best (csvar+sr, slvar) combination.
                            double best_total = BV_INF_CS;
                            double best_mfe2 = 0, best_cai2 = 0;
                            int best_slv_x = -1, best_slv_y = -1;
                            // FIX (phantom-energy bug): when we pick the best `jm1` across
                            // the 4 options, we also implicitly pick a specific CS variant
                            // whose right-S-boundary nucleotide equals `jm1`. Without
                            // recording which, traceback can select a DIFFERENT CS variant,
                            // and the final rna's nucleotide at pos-1 won't match the one
                            // the fill assumed → phantom over-credit of loop_e2. Record
                            // best_csvar_y so traceback can constrain the CS lookup.
                            int best_csvar_y = -1;
                            int best_jm1 = -1;

                            // Helper lambda: try a single (sr_mfe, sr_cai, jm1, cand_csvar_y) +
                            // a single sl_pick combination, update best_*.
                            auto try_emit = [&](double sr_mfe_j, double sr_cai_j,
                                                int jm1, int cand_csvar_y,
                                                const XYVariant* sl_pick, double loop_e2) {
                                if (!sl_pick) return;
                                double mfe2 = inner_mfe_v + sr_mfe_j + sl_pick->mfe + loop_e2;
                                double cai2 = inner_cai_v + sr_cai_j + cai_yy_eff + sl_pick->cai
                                    + (slt.ii_o2 == 2 ? codon_cai[slt.paa_o2][xx_o2] : 0.0);
                                double total = combined_score(lambda, mfe2, cai2);
                                if (total < best_total) {
                                    best_total = total;
                                    best_mfe2 = mfe2;
                                    best_cai2 = cai2;
                                    best_slv_x = sl_pick->x;
                                    best_slv_y = sl_pick->y;
                                    best_csvar_y = cand_csvar_y;
                                    best_jm1 = jm1;
                                }
                            };

                            // Helper lambda: given (sr_mfe, sr_cai, jm1, cand_csvar_y, loop_e2),
                            // dispatch to single-best sl pick (default) or full sl iteration
                            // (CSPATH_SL_FULL=1).
                            auto handle_sl_pick = [&](double sr_mfe_j, double sr_cai_j,
                                                       int jm1, int cand_csvar_y, double loop_e2) {
                                if (cspath_sl_full) {
                                    // FULL ENUMERATION: iterate all slvars matching seam
                                    // constraints (same_codon_sl ? sv.y==inner_x_v : (any))
                                    // and (xx_o2_locked ? sv.x==xx_o2 : (any)).
                                    for (const auto& sv : slt.s_left_ent->variants) {
                                        if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                        if (slt.same_codon_sl && sv.y != inner_x_v) continue;
                                        if (slt.xx_o2_locked && sv.x != xx_o2) continue;
                                        try_emit(sr_mfe_j, sr_cai_j, jm1, cand_csvar_y, &sv, loop_e2);
                                    }
                                } else {
                                    const XYVariant* sl_pick = nullptr;
                                    if (slt.same_codon_sl) {
                                        if (inner_x_v >= 0 && inner_x_v < 6) {
                                            sl_pick = slt.xx_o2_locked
                                                ? slt.sl_by_xy[xx_o2][inner_x_v].ptr
                                                : slt.sl_by_y[inner_x_v].ptr;
                                        }
                                    } else {
                                        sl_pick = slt.xx_o2_locked
                                            ? slt.sl_by_x[xx_o2].ptr
                                            : slt.sl_all.ptr;
                                    }
                                    try_emit(sr_mfe_j, sr_cai_j, jm1, cand_csvar_y, sl_pick, loop_e2);
                                }
                            };

                            if (same_codon_cs) {
                                if (cspath_cy_full) {
                                    // FULL ENUMERATION over srs matching (nuc_qp1, sv.y==yy).
                                    int njm1_locked = nucleotides[protein[cs_ent.b]][yy][cs_ent.j];
                                    double loop_e2 = score_single_loop_with_mismatch(lambda,
                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                        slt.seg_len_left, seg_len_right,
                                        slt.nuc_i1_sl, njm1_locked, slt.nuc_pm1_sl, nuc_qp1_committed);
                                    for (const auto& sv : s_right_ent.variants) {
                                        if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                        int nx = nucleotides[protein[s_right_ent.a]][sv.x][s_right_ent.i];
                                        if (nx != nuc_qp1_committed) continue;
                                        if (sv.y != yy) continue;
                                        handle_sl_pick(sv.mfe, sv.cai, njm1_locked, yy, loop_e2);
                                    }
                                } else if (cs_best_by_y[yy].valid) {
                                    int jm1_locked = cs_best_by_y[yy].nuc_jm1;
                                    double loop_e2 = score_single_loop_with_mismatch(lambda,
                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                        slt.seg_len_left, seg_len_right,
                                        slt.nuc_i1_sl, jm1_locked, slt.nuc_pm1_sl, nuc_qp1_committed);
                                    handle_sl_pick(cs_best_by_y[yy].sr_mfe, cs_best_by_y[yy].sr_cai,
                                                    jm1_locked, yy, loop_e2);
                                }
                            } else {
                                if (cspath_jm1_full) {
                                    // FULL ENUMERATION over (cy, sr) pairs.
                                    std::array<bool, 6> csvar_y_exists{};
                                    for (const auto& csvar : cs_ent.variants) {
                                        if (!std::isfinite(csvar.cai) || !std::isfinite(csvar.mfe) || cai_looks_garbage(csvar.cai)) continue;
                                        if (csvar.y >= 0 && csvar.y < 6) csvar_y_exists[csvar.y] = true;
                                    }
                                    for (int cy = 0; cy < 6; ++cy) {
                                        if (!csvar_y_exists[cy]) continue;
                                        int njm1 = nucleotides[protein[cs_ent.b]][cy][cs_ent.j];
                                        double loop_e2 = score_single_loop_with_mismatch(lambda,
                                            nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                            slt.seg_len_left, seg_len_right,
                                            slt.nuc_i1_sl, njm1, slt.nuc_pm1_sl, nuc_qp1_committed);
                                        for (const auto& sv : s_right_ent.variants) {
                                            if (!std::isfinite(sv.cai) || !std::isfinite(sv.mfe) || cai_looks_garbage(sv.cai)) continue;
                                            int nx = nucleotides[protein[s_right_ent.a]][sv.x][s_right_ent.i];
                                            if (nx != nuc_qp1_committed) continue;
                                            if (sv.y != cy) continue;
                                            handle_sl_pick(sv.mfe, sv.cai, njm1, cy, loop_e2);
                                        }
                                    }
                                } else {
                                    for (int jm1 = 0; jm1 < 4; ++jm1) {
                                        if (!cs_best_by_jm1[jm1].valid) continue;
                                        double loop_e2 = score_single_loop_with_mismatch(lambda,
                                            nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                            slt.seg_len_left, seg_len_right,
                                            slt.nuc_i1_sl, jm1, slt.nuc_pm1_sl, nuc_qp1_committed);
                                        handle_sl_pick(cs_best_by_jm1[jm1].sr_mfe, cs_best_by_jm1[jm1].sr_cai,
                                                       jm1, cs_best_by_jm1[jm1].csvar_y, loop_e2);
                                    }
                                }
                            }

                            if (best_total >= BV_INF_CS) continue;
                            if (slt.xx_o2_locked && best_slv_x != xx_o2) continue;
                            cnt_S_CS_to_C_local[tid_cs]++;
                            (void)best_jm1;

                            int key2 = nuc_key_c_nuc(sigma(slt.aa_o2, slt.ii_o2), pos, nuc_lo2, nuc_ro_eff, nuc_len);

                            BeamEntry ent2(best_total, slt.aa_o2, bb, slt.ii_o2, slot, xx_o2, yy_eff,
                                          (int)DernaManner::MANNER_S_CStoC, best_mfe2, best_cai2);
                            // bt_info layout (FIX): {s_left_key, 0, slv.x, slv.y,
                            //                        best_csvar_y, seg_len_left, seam_pos_left}
                            // bt[5] = seg_len_left, bt[6] = seam_pos_left = inner_left-1.
                            // The S_CStoC traceback uses these to look up the left-S
                            // directly in tab_s[seam_pos_left][seg_len_left]. The prior
                            // code iterated seg from 1 on the WRONG tab_s row (indexed by
                            // s_left_key's left_pos instead of the seam/right position) and
                            // picked the first-matching key, which could be a DIFFERENT
                            // S_left entry with different codons. That was the variant-
                            // inconsistency source of the phantom-energy bug.
                            int seam_pos_left = cs_ent.cs_inner_left - 1;
                            ent2.bt_info = {slt.s_left_key, 0, best_slv_x, best_slv_y,
                                            best_csvar_y, slt.seg_len_left, seam_pos_left};
                            ent2.cs_pack_outer = (long long)cs_tab_k;
                            check_c_entry_invariant(lambda, ent2, "S+CS->C", pos);
                            if (update_derna(task_map_cs, key2, best_total, ent2)) task_order_cs.push_back(key2);
                        }
                    }
                }  // end yy
                }  // end iv_idx (matching inner-C variants)
            }  // end parallel for over prev_cs_flat tasks

            // -------- Merge phase (mirrors Block 1) --------
            // For each thread, replay its entries (and every variant within) into curr_c
            // via update_derna. Iterate threads in ascending tid, and within each thread in
            // insertion order. Tie-break in update_derna is a total order so the final state
            // is order-independent, but the ordered iteration yields reproducible merges.
            {
                BeamEntry scratch_cs;
                for (int t = 0; t < num_threads_cs; ++t) {
                    DernaBeamMap& tm = cs_thread_maps[t];
                    const auto& order = cs_thread_key_order[t];
                    for (int key : order) {
                        auto it = tm.find(key);
                        if (it == tm.end()) continue;
                        const BeamEntry& src = it->second;
                        for (const auto& v : src.variants) {
                            scratch_cs = src;
                            scratch_cs.x = v.x; scratch_cs.y = v.y;
                            scratch_cs.mfe = v.mfe; scratch_cs.cai = v.cai; scratch_cs.score = v.score;
                            scratch_cs.backtrace_type = v.backtrace_type;
                            scratch_cs.last_closed_nuc = v.last_closed_nuc;
                            scratch_cs.bt_info = v.bt_info;
                            scratch_cs.cs_inner_key    = v.cs_inner_key;
                            scratch_cs.cs_right_s_key  = v.cs_right_s_key;
                            scratch_cs.cs_right_len    = v.cs_right_len;
                            scratch_cs.cs_single_start = v.cs_single_start;
                            scratch_cs.cs_inner_left   = v.cs_inner_left;
                            scratch_cs.cs_inner_right  = v.cs_inner_right;
                            scratch_cs.cs_pack_outer   = v.cs_pack_outer;
                            scratch_cs.variants.clear();
                            update_derna(curr_c, key, v.score, scratch_cs);
                        }
                    }
                }
                for (int t = 0; t < num_threads_cs; ++t) {
                    cnt_CS_to_C   += cnt_CS_to_C_local[t];
                    cnt_S_CS_to_C += cnt_S_CS_to_C_local[t];
                }
            }
        }
        // ---------- Block 1 (LCDSfold): Direct right-bulge S[j-1]+C[q]→C[j] and internal S+C+S→C ----------
        // Mirrors LCDSfold Block 1 (~lines 760-820 in LCDSfold.h):
        //   For each right S ending at pos_prev (length seg_len_right), look up C at closed_right=seg_start-1.
        //   Right-bulge (MANNER_C_StoC): outer_left = inner_left-1, n1=0, n2=seg_len_right
        //   Internal loop (MANNER_S_C_StoC): additionally prepend left S ending at inner_left-1
        // This avoids the intermediate CS composite pruning that limits CS→C to k candidates.
#if USE_CODON_LOOKUP_CLOSURE
        if (pos >= 5) {
            // ---- New codon-pair-iterating closure (Step 3 of agg refactor) ----
            // PERF-OPTIMIZED VERSION (Perf Iteration). Mirrors baseline Block 1's
            // emit cadence: ~1 update_derna call per output (xx_o, yy_eff) slot,
            // not per inner (ci, sl, sr) tuple. Adds:
            //   A) Inverted inner loop: per output slot, find best inner combination
            //      via per-grid best_by_x/y/xy summaries (mirrors baseline's
            //      best_csr_by_xy precomputation).
            //   B) loop_e2 caching (4 x SLL_MAX x 4 x 4 table per inner block).
            //   C) OpenMP parallelization over (seg_len_right, sr_left_pos) tasks
            //      with per-thread dense beam tables and a deterministic merge.
            //
            // The per-grid summaries are partitioned by nucleotide quadrant
            //   c-grid:   (nuc_li, nuc_ri)   — c.x → nuc_li, c.y → nuc_ri
            //   sr-grid:  (nuc_qp1, nuc_jm1) — sr.x → nuc_qp1, sr.y → nuc_jm1
            //   sl-grid:  (nuc_i1, nuc_pm1)  — sl.x → nuc_i1, sl.y → nuc_pm1
            // Each (4x4) quadrant pretends to be one "key" and is treated like a
            // baseline c_ent/sr_ent. That gives correctness across the codon-pair
            // collapse: cells in the same quadrant share boundary nucleotides
            // and so share the same loop_e contribution.
            const auto& s_right_agg_outer = agg_tables.s_agg[pos_prev];

            // ---- Helper: build per-grid quadrant summaries. ----
            // For a given grid, partition cells by (boundary_x_nuc, boundary_y_nuc)
            // quadrant. Within each populated quadrant, build:
            //   best_by_x[6]   : min score per (cell.x), -1 if absent
            //   best_by_y[6]   : min score per (cell.y), -1 if absent
            //   best_by_xy[6][6] : the cell at (x,y) score, -1 if absent
            //   best_all       : overall best score in quadrant
            // Each chosen "best" remembers the cell idx for backtrace.
            struct AggBest {
                double score;
                float mfe, cai;
                int8_t cell_idx;     // index into grid.cells (= x*NCOD+y), -1 if invalid
                bool valid;
            };
            static constexpr double AGG_INF = 1e30;
            struct QuadrantSummary {
                std::array<AggBest, 6> by_x;
                std::array<AggBest, 6> by_y;
                std::array<std::array<AggBest, 6>, 6> by_xy;
                AggBest all;
                bool active = false;
            };
            // Quadrant indexed by (nuc_a * 4 + nuc_b) where nuc_a is x-side nuc, nuc_b is y-side nuc.
            using GridQuadrants = std::array<QuadrantSummary, 16>;

            // Lambda to build per-quadrant summaries for an agg grid.
            //   xy_to_nuc_a[xx] = nucleotide for cell.x = xx (boundary on x-side)
            //   xy_to_nuc_b[yy] = nucleotide for cell.y = yy (boundary on y-side)
            // Cells with x or y outside valid ncod range are skipped.
            // Sentinel "invalid" used for AggBest fields. Building lazily: we only
            // initialize quadrants that get touched. Caller must inspect `active_qids`
            // (a list of quadrant ids that became active) to drive iteration.
            auto build_grid_quadrants = [&](const CodonAggGrid& grid,
                                            int paa_a, int slot_a, int ncod_a,
                                            int paa_b, int slot_b, int ncod_b,
                                            GridQuadrants& out,
                                            std::array<int8_t, 16>& active_qids,
                                            int& num_active) {
                num_active = 0;
                for (uint8_t ci_idx : grid.populated_indices) {
                    const CodonAggEntry& cell = grid.cells[ci_idx];
                    if (!cell.valid) continue;
                    const int xx = ci_idx / NCOD;
                    const int yy = ci_idx % NCOD;
                    if (xx >= ncod_a) continue;
                    if (yy >= ncod_b) continue;
                    const int nuc_a = nucleotides[paa_a][xx][slot_a];
                    const int nuc_b = nucleotides[paa_b][yy][slot_b];
                    const int qid = nuc_a * 4 + nuc_b;
                    QuadrantSummary& q = out[qid];
                    if (!q.active) {
                        // Lazy init this quadrant's tables.
                        q.active = true;
                        q.all = AggBest{AGG_INF, 0.0f, 0.0f, -1, false};
                        for (auto& bv : q.by_x) bv = AggBest{AGG_INF, 0.0f, 0.0f, -1, false};
                        for (auto& bv : q.by_y) bv = AggBest{AGG_INF, 0.0f, 0.0f, -1, false};
                        for (auto& row : q.by_xy)
                            for (auto& bv : row) bv = AggBest{AGG_INF, 0.0f, 0.0f, -1, false};
                        if (num_active < 16) active_qids[num_active++] = (int8_t)qid;
                    }
                    AggBest cand{(double)cell.score, cell.mfe, cell.cai, (int8_t)ci_idx, true};
                    if (cand.score < q.by_x[xx].score) q.by_x[xx] = cand;
                    if (cand.score < q.by_y[yy].score) q.by_y[yy] = cand;
                    if (cand.score < q.by_xy[xx][yy].score) q.by_xy[xx][yy] = cand;
                    if (cand.score < q.all.score) q.all = cand;
                }
            };

            // Build SR grid quadrants per (seg_len_right). Indexed by sr-grid pointer.
            // We do this per task within the OpenMP loop below.

            // ---- Build a list of (seg_len_right, sr-grid-pos seg_start) tasks. ----
            // For determinism in parallel merge, tasks are ordered by (seg_len_right, seg_start).
            struct Blk1AggTask {
                int seg_len_right;
                int seg_start;
                int inner_right_pos;
            };
            std::vector<Blk1AggTask> agg_tasks;
            agg_tasks.reserve(64);
            for (int seg_len_right = 1; seg_len_right <= min(SINGLE_MAX_LEN, pos_prev); ++seg_len_right) {
                if (seg_len_right >= (int)s_right_agg_outer.size()) continue;
                const auto& s_right_agg_map = s_right_agg_outer[seg_len_right];
                if (s_right_agg_map.empty()) continue;
                const int seg_start = pos_prev - seg_len_right + 1;
                const int inner_right_pos = seg_start - 1;
                if (inner_right_pos < 0 || inner_right_pos >= nuc_len) continue;
                if (s_right_agg_map.find(seg_start) == s_right_agg_map.end()) continue;
                if (inner_right_pos >= (int)agg_tables.c_agg.size()) continue;
                if (agg_tables.c_agg[inner_right_pos].empty()) continue;
                agg_tasks.push_back({seg_len_right, seg_start, inner_right_pos});
            }
            const int num_agg_tasks = (int)agg_tasks.size();

            // ---- Parallelization setup (mirrors baseline Block 1). ----
#ifdef _OPENMP
            const int num_threads_blk1c = omp_get_max_threads();
#else
            const int num_threads_blk1c = 1;
#endif
            // Reuse function-scope blk1_thread_dense_maps from baseline. Reset
            // populated slots from prior position(s).
            const int dense_cap_c = nuc_len * 16;
            for (int t = 0; t < num_threads_blk1_outer; ++t) {
                blk1_thread_dense_maps[t].reset(dense_cap_c);
                blk1_thread_key_order[t].clear();
            }
            std::vector<uint64_t> cnt_rb_local_c(num_threads_blk1c, 0);
            std::vector<uint64_t> cnt_in_local_c(num_threads_blk1c, 0);

            const bool pos_is_last_g = is_last_nuc(pos);

            #pragma omp parallel for schedule(dynamic, 1) if(num_agg_tasks > 1)
            for (int task_idx = 0; task_idx < num_agg_tasks; ++task_idx) {
                const Blk1AggTask& tk = agg_tasks[task_idx];
                const int seg_len_right = tk.seg_len_right;
                const int seg_start = tk.seg_start;
                const int inner_right_pos = tk.inner_right_pos;

#ifdef _OPENMP
                const int tid_c = omp_get_thread_num();
#else
                const int tid_c = 0;
#endif
                DenseBeamTable& task_dense = blk1_thread_dense_maps[tid_c];
                std::vector<int>& task_order = blk1_thread_key_order[tid_c];

                const auto& s_right_agg_map = s_right_agg_outer[seg_len_right];
                auto sr_grid_it = s_right_agg_map.find(seg_start);
                if (sr_grid_it == s_right_agg_map.end()) continue;
                const CodonAggGrid& sr_grid = sr_grid_it->second;
                if (sr_grid.populated_indices.empty()) continue;

                // S_right boundary aa/slot.
                const int s_right_b_aa = pos_prev / 3;
                const int s_right_b_slot = pos_prev % 3;
                const int p_s_right_b = protein[s_right_b_aa];
                const int ncod_sr_b = n_codon[p_s_right_b];
                const int s_right_a_aa = seg_start / 3;
                const int s_right_a_slot = seg_start % 3;
                const int p_s_right_a = protein[s_right_a_aa];
                const int ncod_sr_a = n_codon[p_s_right_a];

                bool same_codon_sr = false;
                if (!same_or_next_codon(s_right_b_aa, s_right_b_slot, bb, slot, same_codon_sr)) continue;

                const bool seam_sr_inner_same_codon = ((inner_right_pos % 3) != 2);

                // Build SR grid quadrant summaries. Indexed by (nuc_qp1, nuc_jm1).
                GridQuadrants sr_quadrants{};
                std::array<int8_t, 16> sr_active_qids{};
                int sr_num_active = 0;
                build_grid_quadrants(sr_grid,
                    p_s_right_a, s_right_a_slot, ncod_sr_a,
                    p_s_right_b, s_right_b_slot, ncod_sr_b,
                    sr_quadrants, sr_active_qids, sr_num_active);

                // Inner C agg map.
                const auto& c_inner_agg_map = agg_tables.c_agg[inner_right_pos];
                if (c_inner_agg_map.empty()) continue;

                const int inner_b_aa = inner_right_pos / 3;
                const int inner_b_slot = inner_right_pos % 3;
                const int p_inner_b = protein[inner_b_aa];
                const int ncod_inner_b = n_codon[p_inner_b];

                // Scratch buffers reused across emits within this task.
                BeamEntry scratch_rb;
                BeamEntry scratch_il;

                for (const auto& c_kv : c_inner_agg_map) {
                    const int inner_left = c_kv.first;
                    const CodonAggGrid& c_inner_grid = c_kv.second;
                    if (c_inner_grid.populated_indices.empty()) continue;
                    const int outer_left_nuc = inner_left - 1;
                    if (outer_left_nuc < 0) continue;
                    if ((pos - outer_left_nuc) <= (HAIRPIN_GAP + 2)) continue;

                    const int inner_a_aa = inner_left / 3;
                    const int inner_a_slot = inner_left % 3;
                    const int p_inner_a = protein[inner_a_aa];
                    const int ncod_inner_a = n_codon[p_inner_a];

                    const int aa_o = outer_left_nuc / 3;
                    const int ii_o = outer_left_nuc % 3;
                    const int paa_o = protein[aa_o];
                    const int ncod_o = n_codon[paa_o];
                    const bool xx_o_locked = (aa_o == inner_a_aa);

                    // Build C grid quadrant summaries. Indexed by (nuc_li, nuc_ri).
                    // Reuse outer-scope storage by zeroing only previously-active
                    // quadrants.
                    GridQuadrants c_quadrants{};
                    std::array<int8_t, 16> c_active_qids{};
                    int c_num_active = 0;
                    build_grid_quadrants(c_inner_grid,
                        p_inner_a, inner_a_slot, ncod_inner_a,
                        p_inner_b, inner_b_slot, ncod_inner_b,
                        c_quadrants, c_active_qids, c_num_active);

                    // Precompute best xx_o per nuc_lo (when !xx_o_locked).
                    struct NucBest { int xx; double cai; bool valid; };
                    std::array<NucBest, 4> best_xxo_per_nuc;
                    for (auto& nb : best_xxo_per_nuc) nb = {-1, 0.0, false};
                    if (!xx_o_locked) {
                        for (int xo = 0; xo < ncod_o; ++xo) {
                            int nlo = nucleotides[paa_o][xo][ii_o];
                            double xcai = (ii_o == 2) ? codon_cai[paa_o][xo] : 0.0;
                            if (!best_xxo_per_nuc[nlo].valid || xcai > best_xxo_per_nuc[nlo].cai) {
                                best_xxo_per_nuc[nlo] = {xo, xcai, true};
                            }
                        }
                    }

                    // ============================================================
                    // Pre-build SL tasks once per (c_grid). Mirrors baseline sl_tasks.
                    // ============================================================
                    struct SLTaskAgg {
                        int sl_left_pos;
                        int seg_len_left;
                        int outer_left2;
                        int aa_o2, ii_o2, paa_o2, ncod_o2;
                        bool xx_o2_locked;
                        bool same_codon_sl;
                        int p_sl_a, p_sl_b;
                        int sl_a_slot, sl_b_slot;
                        int ncod_sl_a, ncod_sl_b;
                        const CodonAggGrid* grid;
                        GridQuadrants quadrants;  // indexed by (nuc_i1, nuc_pm1)
                        std::array<int8_t, 16> active_qids;
                        int num_active;
                        std::array<NucBest, 4> best_xx_per_nuc;
                    };
                    std::vector<SLTaskAgg> sl_tasks_agg;
                    sl_tasks_agg.reserve(SINGLE_MAX_LEN);
                    {
                        const int seam_pos_left = outer_left_nuc;
                        if (seam_pos_left >= 1 && seam_pos_left < nuc_len &&
                            seam_pos_left < (int)agg_tables.s_agg.size()) {
                            const auto& sl_agg_outer = agg_tables.s_agg[seam_pos_left];
                            const int max_sl = min(SINGLE_MAX_LEN, seam_pos_left + 1);
                            for (int seg_len_left = 1; seg_len_left <= max_sl; ++seg_len_left) {
                                if (seg_len_left >= (int)sl_agg_outer.size()) continue;
                                const auto& sl_agg_map = sl_agg_outer[seg_len_left];
                                if (sl_agg_map.empty()) continue;
                                const int outer_left2 = inner_left - 1 - seg_len_left;
                                if (outer_left2 < 0) continue;
                                if ((pos - outer_left2) <= (HAIRPIN_GAP + 2)) continue;
                                const int sl_left_pos = inner_left - seg_len_left;
                                auto sl_it = sl_agg_map.find(sl_left_pos);
                                if (sl_it == sl_agg_map.end()) continue;
                                const CodonAggGrid& sl_grid = sl_it->second;
                                if (sl_grid.populated_indices.empty()) continue;

                                SLTaskAgg slt;
                                slt.sl_left_pos = sl_left_pos;
                                slt.seg_len_left = seg_len_left;
                                slt.outer_left2 = outer_left2;
                                slt.aa_o2 = outer_left2 / 3;
                                slt.ii_o2 = outer_left2 % 3;
                                slt.paa_o2 = protein[slt.aa_o2];
                                slt.ncod_o2 = n_codon[slt.paa_o2];
                                slt.same_codon_sl = ((inner_left % 3) != 0);
                                slt.sl_a_slot = sl_left_pos % 3;
                                slt.sl_b_slot = seam_pos_left % 3;
                                slt.p_sl_a = protein[sl_left_pos / 3];
                                slt.p_sl_b = protein[seam_pos_left / 3];
                                slt.ncod_sl_a = n_codon[slt.p_sl_a];
                                slt.ncod_sl_b = n_codon[slt.p_sl_b];
                                slt.xx_o2_locked = (slt.aa_o2 == sl_left_pos / 3);
                                slt.grid = &sl_grid;
                                slt.num_active = 0;
                                // SLTaskAgg's quadrants is fresh (default-initialized
                                // from struct ctor: active=false for all 16).
                                build_grid_quadrants(sl_grid,
                                    slt.p_sl_a, slt.sl_a_slot, slt.ncod_sl_a,
                                    slt.p_sl_b, slt.sl_b_slot, slt.ncod_sl_b,
                                    slt.quadrants, slt.active_qids, slt.num_active);
                                for (auto& nb : slt.best_xx_per_nuc) nb = {-1, 0.0, false};
                                for (int xo = 0; xo < slt.ncod_o2; ++xo) {
                                    int nlo = nucleotides[slt.paa_o2][xo][slt.ii_o2];
                                    double xcai = (slt.ii_o2 == 2) ? codon_cai[slt.paa_o2][xo] : 0.0;
                                    if (!slt.best_xx_per_nuc[nlo].valid ||
                                        xcai > slt.best_xx_per_nuc[nlo].cai) {
                                        slt.best_xx_per_nuc[nlo] = {xo, xcai, true};
                                    }
                                }
                                sl_tasks_agg.push_back(std::move(slt));
                            }
                        }
                    }

                    // ============================================================
                    // Lambda emit_rb_agg: emit one (xx_o, yy_eff) right-bulge candidate.
                    // ============================================================
                    auto emit_rb_agg = [&](const AggBest& c_best, const AggBest& sr_best,
                                           int xx_inner, int yy_inner,
                                           int xx_sr, int yy_sr,
                                           int nuc_li, int nuc_ri,
                                           int nuc_qp1, int nuc_jm1,
                                           int xx_o, int yy_eff,
                                           int nuc_lo, int nuc_ro_eff,
                                           double loop_e, double cai_yy) {
                        const double mfe_cand = (double)c_best.mfe + (double)sr_best.mfe + loop_e;
                        const double cai_cand = (double)c_best.cai + (double)sr_best.cai + cai_yy
                            + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                        const double sc_cand = combined_score(lambda, mfe_cand, cai_cand);
                        const int compact_idx = outer_left_nuc * 16 + nuc_lo * 4 + nuc_ro_eff;
                        // Early-rejection like baseline BLK1_AGG.
                        const BeamEntry& ex_pre = task_dense.slot(compact_idx);
                        if (ex_pre.score < inf &&
                            xx_o >= 0 && xx_o < 6 && yy_eff >= 0 && yy_eff < 6) {
                            const int8_t vi_pre = ex_pre.variant_idx[xx_o][yy_eff];
                            if (vi_pre >= 0 && sc_cand > ex_pre.variants[vi_pre].score) return;
                        }
                        const int c_inner_key = nuc_key_c_nuc(inner_left, inner_right_pos,
                                                              nuc_li, nuc_ri, nuc_len);
                        const int s_right_key = codon_beam_key(seg_start, nuc_qp1, nuc_jm1);
                        scratch_rb.score = sc_cand;
                        scratch_rb.mfe = mfe_cand;
                        scratch_rb.cai = cai_cand;
                        scratch_rb.a = aa_o; scratch_rb.b = bb;
                        scratch_rb.i = ii_o; scratch_rb.j = slot;
                        scratch_rb.x = xx_o; scratch_rb.y = yy_eff;
                        scratch_rb.backtrace_type = (int)DernaManner::MANNER_C_StoC;
                        scratch_rb.bt_info.len = 6;
                        scratch_rb.bt_info.data[0] = c_inner_key;
                        scratch_rb.bt_info.data[1] = s_right_key;
                        scratch_rb.bt_info.data[2] = xx_inner;
                        scratch_rb.bt_info.data[3] = yy_inner;
                        scratch_rb.bt_info.data[4] = xx_sr;
                        scratch_rb.bt_info.data[5] = yy_sr;
                        scratch_rb.cs_inner_left = inner_right_pos;
                        scratch_rb.cs_right_len = seg_len_right;
                        scratch_rb.last_closed_nuc = -1;
                        scratch_rb.cs_inner_key = -1;
                        scratch_rb.cs_right_s_key = -1;
                        scratch_rb.cs_single_start = -1;
                        scratch_rb.cs_inner_right = -1;
                        scratch_rb.cs_pack_outer = -1LL;
                        if (task_dense.update(compact_idx, sc_cand, scratch_rb))
                            task_order.push_back(compact_idx);
                    };

                    auto emit_il_agg = [&](const AggBest& c_best, const AggBest& sr_best,
                                           const AggBest& sl_best,
                                           int xx_inner, int yy_inner,
                                           int xx_sr, int yy_sr,
                                           int xx_sl, int yy_sl,
                                           int nuc_li, int nuc_ri,
                                           int nuc_qp1, int nuc_jm1,
                                           int nuc_i1, int nuc_pm1,
                                           int xx_o2, int yy_eff,
                                           int nuc_lo2, int nuc_ro_eff,
                                           double loop_e2, double cai_yy,
                                           int sl_left_pos, int outer_left2,
                                           int aa_o2, int ii_o2, int paa_o2,
                                           int seg_len_left) {
                        const double mfe2 = (double)c_best.mfe + (double)sr_best.mfe
                                          + (double)sl_best.mfe + loop_e2;
                        const double cai2 = (double)c_best.cai + (double)sr_best.cai
                                          + (double)sl_best.cai + cai_yy
                                          + (ii_o2 == 2 ? codon_cai[paa_o2][xx_o2] : 0.0);
                        const double sc2 = combined_score(lambda, mfe2, cai2);
                        const int compact_idx2 = outer_left2 * 16 + nuc_lo2 * 4 + nuc_ro_eff;
                        const BeamEntry& ex_pre = task_dense.slot(compact_idx2);
                        if (ex_pre.score < inf &&
                            xx_o2 >= 0 && xx_o2 < 6 && yy_eff >= 0 && yy_eff < 6) {
                            const int8_t vi_pre = ex_pre.variant_idx[xx_o2][yy_eff];
                            if (vi_pre >= 0 && sc2 > ex_pre.variants[vi_pre].score) return;
                        }
                        const int c_inner_key = nuc_key_c_nuc(inner_left, inner_right_pos,
                                                              nuc_li, nuc_ri, nuc_len);
                        const int s_right_key = codon_beam_key(seg_start, nuc_qp1, nuc_jm1);
                        const int s_left_key = codon_beam_key(sl_left_pos, nuc_i1, nuc_pm1);
                        scratch_il.score = sc2;
                        scratch_il.mfe = mfe2;
                        scratch_il.cai = cai2;
                        scratch_il.a = aa_o2; scratch_il.b = bb;
                        scratch_il.i = ii_o2; scratch_il.j = slot;
                        scratch_il.x = xx_o2; scratch_il.y = yy_eff;
                        scratch_il.backtrace_type = (int)DernaManner::MANNER_S_C_StoC;
                        scratch_il.bt_info.len = 9;
                        scratch_il.bt_info.data[0] = s_left_key;
                        scratch_il.bt_info.data[1] = c_inner_key;
                        scratch_il.bt_info.data[2] = s_right_key;
                        scratch_il.bt_info.data[3] = xx_sl;
                        scratch_il.bt_info.data[4] = yy_sl;
                        scratch_il.bt_info.data[5] = xx_inner;
                        scratch_il.bt_info.data[6] = yy_inner;
                        scratch_il.bt_info.data[7] = xx_sr;
                        scratch_il.bt_info.data[8] = yy_sr;
                        scratch_il.cs_inner_left = inner_right_pos;
                        scratch_il.cs_right_len = seg_len_right;
                        scratch_il.cs_single_start = outer_left_nuc;
                        scratch_il.cs_inner_right = seg_len_left;
                        scratch_il.last_closed_nuc = -1;
                        scratch_il.cs_inner_key = -1;
                        scratch_il.cs_right_s_key = -1;
                        scratch_il.cs_pack_outer = -1LL;
                        if (task_dense.update(compact_idx2, sc2, scratch_il))
                            task_order.push_back(compact_idx2);
                    };

                    // ============================================================
                    // Iterate (c-quadrant, sr-quadrant) pairs. Mirrors baseline's
                    // (c_ent, sr_ent) outer iteration: each quadrant pair fixes
                    // (nuc_li, nuc_ri, nuc_qp1, nuc_jm1).
                    // ============================================================
                    for (int c_ai = 0; c_ai < c_num_active; ++c_ai) {
                        const int c_qid = c_active_qids[c_ai];
                        const QuadrantSummary& cq = c_quadrants[c_qid];
                        const int nuc_li  = c_qid >> 2;
                        const int nuc_ri  = c_qid & 3;

                        for (int sr_ai = 0; sr_ai < sr_num_active; ++sr_ai) {
                            const int sr_qid = sr_active_qids[sr_ai];
                            const QuadrantSummary& srq = sr_quadrants[sr_qid];
                            const int nuc_qp1 = sr_qid >> 2;
                            const int nuc_jm1 = sr_qid & 3;

                            // Iterate yy_eff (outer right codon).
                            for (int yy_eff = 0; yy_eff < ncod_bb; ++yy_eff) {
                                if (same_codon_sr) {
                                    // sr.y must equal yy_eff.
                                    if (!srq.by_y[yy_eff].valid) continue;
                                }
                                const int nuc_ro_eff = nucleotides[pbb][yy_eff][slot];
                                const double cai_yy = pos_is_last_g ? codon_cai[pbb][yy_eff] : 0.0;

                                // Compute best (c + sr) per c.x  (= xx_inner) — call it best_csr_by_cx.
                                // Mirrors baseline's separable logic; the "inner_key" relation is:
                                //   seam_sr_inner_same_codon → sr.x == c.y
                                //   same_codon_sr             → sr.y == yy_eff (already filtered above)
                                struct BestCSR {
                                    double score;
                                    AggBest c_best;
                                    AggBest sr_best;
                                    bool valid;
                                };
                                std::array<BestCSR, 6> best_csr_by_cx;
                                BestCSR best_csr_all{AGG_INF, AggBest{AGG_INF,0,0,-1,false}, AggBest{AGG_INF,0,0,-1,false}, false};
                                for (auto& b : best_csr_by_cx)
                                    b = {AGG_INF, AggBest{AGG_INF,0,0,-1,false}, AggBest{AGG_INF,0,0,-1,false}, false};

                                // SR pick: depends on (seam_sr_inner_same_codon, same_codon_sr).
                                if (!seam_sr_inner_same_codon && !same_codon_sr) {
                                    // SR fully free → pick srq.all.
                                    if (!srq.all.valid) continue;
                                    for (int cx = 0; cx < 6; ++cx) {
                                        if (!cq.by_x[cx].valid) continue;
                                        double sc = cq.by_x[cx].score + srq.all.score;
                                        best_csr_by_cx[cx] = {sc, cq.by_x[cx], srq.all, true};
                                        if (sc < best_csr_all.score) best_csr_all = best_csr_by_cx[cx];
                                    }
                                } else if (!seam_sr_inner_same_codon /*  same_codon_sr */) {
                                    if (!srq.by_y[yy_eff].valid) continue;
                                    const auto& sr_pick = srq.by_y[yy_eff];
                                    for (int cx = 0; cx < 6; ++cx) {
                                        if (!cq.by_x[cx].valid) continue;
                                        double sc = cq.by_x[cx].score + sr_pick.score;
                                        best_csr_by_cx[cx] = {sc, cq.by_x[cx], sr_pick, true};
                                        if (sc < best_csr_all.score) best_csr_all = best_csr_by_cx[cx];
                                    }
                                } else if (!same_codon_sr /* seam_sr_inner_same_codon */) {
                                    // sr.x == c.y. Use cq.by_xy[cx][cy] paired with srq.by_x[cy].
                                    for (int cx = 0; cx < 6; ++cx) {
                                        for (int cy = 0; cy < 6; ++cy) {
                                            const auto& bc = cq.by_xy[cx][cy];
                                            if (!bc.valid) continue;
                                            const auto& bsr = srq.by_x[cy];
                                            if (!bsr.valid) continue;
                                            double sc = bc.score + bsr.score;
                                            if (sc < best_csr_by_cx[cx].score)
                                                best_csr_by_cx[cx] = {sc, bc, bsr, true};
                                        }
                                        if (best_csr_by_cx[cx].valid &&
                                            best_csr_by_cx[cx].score < best_csr_all.score)
                                            best_csr_all = best_csr_by_cx[cx];
                                    }
                                } else {
                                    // Both: sr.x == c.y AND sr.y == yy_eff.
                                    for (int cx = 0; cx < 6; ++cx) {
                                        for (int cy = 0; cy < 6; ++cy) {
                                            const auto& bc = cq.by_xy[cx][cy];
                                            if (!bc.valid) continue;
                                            const auto& bsr = srq.by_xy[cy][yy_eff];
                                            if (!bsr.valid) continue;
                                            double sc = bc.score + bsr.score;
                                            if (sc < best_csr_by_cx[cx].score)
                                                best_csr_by_cx[cx] = {sc, bc, bsr, true};
                                        }
                                        if (best_csr_by_cx[cx].valid &&
                                            best_csr_by_cx[cx].score < best_csr_all.score)
                                            best_csr_all = best_csr_by_cx[cx];
                                    }
                                }
                                if (!best_csr_all.valid) continue;

                                // ---- (1) Right-bulge: emit per nuc_lo, then per xx_o. ----
                                for (int nuc_lo = 0; nuc_lo < 4; ++nuc_lo) {
                                    if (!best_xxo_per_nuc[nuc_lo].valid && !xx_o_locked) continue;
                                    if (BP_pair[nuc_lo + 1][nuc_ro_eff + 1] == 0) continue;
                                    const double loop_e = score_single_loop_with_mismatch(lambda,
                                        nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                        0, seg_len_right,
                                        nuc_li, nuc_jm1, nuc_lo, nuc_qp1);
                                    if (loop_e >= inf) continue;
                                    for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                                        if (nucleotides[paa_o][xx_o][ii_o] != nuc_lo) continue;
                                        // Pick CSR: locked → best_csr_by_cx[xx_o]; free → all.
                                        const BestCSR* csr_pick = xx_o_locked
                                            ? &best_csr_by_cx[xx_o] : &best_csr_all;
                                        if (!csr_pick->valid) continue;
                                        cnt_rb_local_c[tid_c]++;
                                        const AggBest& c_b = csr_pick->c_best;
                                        const AggBest& sr_b = csr_pick->sr_best;
                                        const int xx_inner = c_b.cell_idx / NCOD;
                                        const int yy_inner = c_b.cell_idx % NCOD;
                                        const int xx_sr = sr_b.cell_idx / NCOD;
                                        const int yy_sr = sr_b.cell_idx % NCOD;
                                        emit_rb_agg(c_b, sr_b, xx_inner, yy_inner, xx_sr, yy_sr,
                                                    nuc_li, nuc_ri, nuc_qp1, nuc_jm1,
                                                    xx_o, yy_eff, nuc_lo, nuc_ro_eff,
                                                    loop_e, cai_yy);
                                    }
                                }

                                // ---- (2) Internal loop: iterate sl_tasks. ----
                                // loop_e2 cache: (nuc_lo2, sll_idx, nuc_i1, nuc_pm1) → energy.
                                // Reset on each yy_eff iteration since nuc_ro_eff changes
                                // (nuc_ro_eff is part of loop_e2 inputs but not part of the
                                // 4-d cache key — caller provides ro). Cache keyed by remaining
                                // varying inputs. nuc_ro_eff is fixed for this inner block.
                                static constexpr int SLL_MAX_C = SINGLE_MAX_LEN;
                                const bool use_loop_table_c = ((int)sl_tasks_agg.size() >= 2);
                                double loop_e2_table_c[4][SLL_MAX_C][4][4];
                                int8_t loop_e2_table_state_c[4][SLL_MAX_C][4][4];
                                if (use_loop_table_c) {
                                    memset(loop_e2_table_state_c, 0, sizeof(loop_e2_table_state_c));
                                }

                                for (const SLTaskAgg& slt : sl_tasks_agg) {
                                    const int sll_idx = slt.seg_len_left - 1;
                                    if (sll_idx < 0 || sll_idx >= SLL_MAX_C) continue;

                                    // Iterate sl quadrants.
                                    for (int sl_ai = 0; sl_ai < slt.num_active; ++sl_ai) {
                                        const int sl_qid = slt.active_qids[sl_ai];
                                        const QuadrantSummary& slq = slt.quadrants[sl_qid];
                                        const int nuc_i1  = sl_qid >> 2;
                                        const int nuc_pm1 = sl_qid & 3;

                                        // Combine best_csr with sl_pick. Constraints:
                                        //   slt.same_codon_sl → sl.y == c.x (== xx_inner from c_best)
                                        //   slt.xx_o2_locked  → xx_o2 == sl.x
                                        // Need to pick the winning (c, sr, sl) combination per xx_o2.

                                        // Fast path A: !same_codon_sl && !xx_o2_locked → pick slq.all and best_csr_all.
                                        if (!slt.same_codon_sl && !slt.xx_o2_locked) {
                                            if (!best_csr_all.valid) continue;
                                            if (!slq.all.valid) continue;
                                            for (int nuc_lo2 = 0; nuc_lo2 < 4; ++nuc_lo2) {
                                                const auto& nb = slt.best_xx_per_nuc[nuc_lo2];
                                                if (!nb.valid) continue;
                                                if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                                double loop_e2;
                                                if (use_loop_table_c) {
                                                    auto& st = loop_e2_table_state_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                                    if (st == 0) {
                                                        double le = score_single_loop_with_mismatch(lambda,
                                                            nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                            slt.seg_len_left, seg_len_right,
                                                            nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                        if (le < inf) {
                                                            loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1] = le;
                                                            st = 1;
                                                        } else { st = -1; }
                                                    }
                                                    if (st < 0) continue;
                                                    loop_e2 = loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                                } else {
                                                    loop_e2 = score_single_loop_with_mismatch(lambda,
                                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                        slt.seg_len_left, seg_len_right,
                                                        nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                    if (loop_e2 >= inf) continue;
                                                }
                                                const AggBest& sl_b = slq.all;
                                                const AggBest& c_b = best_csr_all.c_best;
                                                const AggBest& sr_b = best_csr_all.sr_best;
                                                const int xx_inner = c_b.cell_idx / NCOD;
                                                const int yy_inner = c_b.cell_idx % NCOD;
                                                const int xx_sr = sr_b.cell_idx / NCOD;
                                                const int yy_sr = sr_b.cell_idx % NCOD;
                                                const int xx_sl = sl_b.cell_idx / NCOD;
                                                const int yy_sl = sl_b.cell_idx % NCOD;
                                                for (int xo = 0; xo < slt.ncod_o2; ++xo) {
                                                    if (nucleotides[slt.paa_o2][xo][slt.ii_o2] != nuc_lo2) continue;
                                                    cnt_in_local_c[tid_c]++;
                                                    emit_il_agg(c_b, sr_b, sl_b,
                                                                xx_inner, yy_inner,
                                                                xx_sr, yy_sr,
                                                                xx_sl, yy_sl,
                                                                nuc_li, nuc_ri, nuc_qp1, nuc_jm1,
                                                                nuc_i1, nuc_pm1,
                                                                xo, yy_eff, nuc_lo2, nuc_ro_eff,
                                                                loop_e2, cai_yy,
                                                                slt.sl_left_pos, slt.outer_left2,
                                                                slt.aa_o2, slt.ii_o2, slt.paa_o2,
                                                                slt.seg_len_left);
                                                }
                                            }
                                            continue;
                                        }

                                        // Fast path B: same_codon_sl && !xx_o2_locked.
                                        // sl.y == c.x = xx_inner (= cx). Pick best across cx of (csr_by_cx + slq.by_y[cx]).
                                        if (slt.same_codon_sl && !slt.xx_o2_locked) {
                                            BestCSR sf_pick{AGG_INF, AggBest{AGG_INF,0,0,-1,false}, AggBest{AGG_INF,0,0,-1,false}, false};
                                            AggBest sf_sl{AGG_INF, 0, 0, -1, false};
                                            for (int cx = 0; cx < 6; ++cx) {
                                                if (!best_csr_by_cx[cx].valid) continue;
                                                const auto& sl_pick = slq.by_y[cx];
                                                if (!sl_pick.valid) continue;
                                                double total = best_csr_by_cx[cx].score + sl_pick.score;
                                                if (total < sf_pick.score) {
                                                    sf_pick = {total, best_csr_by_cx[cx].c_best,
                                                               best_csr_by_cx[cx].sr_best, true};
                                                    sf_sl = sl_pick;
                                                }
                                            }
                                            if (!sf_pick.valid) continue;
                                            for (int nuc_lo2 = 0; nuc_lo2 < 4; ++nuc_lo2) {
                                                const auto& nb = slt.best_xx_per_nuc[nuc_lo2];
                                                if (!nb.valid) continue;
                                                if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                                double loop_e2;
                                                if (use_loop_table_c) {
                                                    auto& st = loop_e2_table_state_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                                    if (st == 0) {
                                                        double le = score_single_loop_with_mismatch(lambda,
                                                            nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                            slt.seg_len_left, seg_len_right,
                                                            nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                        if (le < inf) {
                                                            loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1] = le;
                                                            st = 1;
                                                        } else { st = -1; }
                                                    }
                                                    if (st < 0) continue;
                                                    loop_e2 = loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                                } else {
                                                    loop_e2 = score_single_loop_with_mismatch(lambda,
                                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                        slt.seg_len_left, seg_len_right,
                                                        nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                    if (loop_e2 >= inf) continue;
                                                }
                                                const AggBest& c_b = sf_pick.c_best;
                                                const AggBest& sr_b = sf_pick.sr_best;
                                                const int xx_inner = c_b.cell_idx / NCOD;
                                                const int yy_inner = c_b.cell_idx % NCOD;
                                                const int xx_sr = sr_b.cell_idx / NCOD;
                                                const int yy_sr = sr_b.cell_idx % NCOD;
                                                const int xx_sl = sf_sl.cell_idx / NCOD;
                                                const int yy_sl = sf_sl.cell_idx % NCOD;
                                                for (int xo = 0; xo < slt.ncod_o2; ++xo) {
                                                    if (nucleotides[slt.paa_o2][xo][slt.ii_o2] != nuc_lo2) continue;
                                                    cnt_in_local_c[tid_c]++;
                                                    emit_il_agg(c_b, sr_b, sf_sl,
                                                                xx_inner, yy_inner,
                                                                xx_sr, yy_sr,
                                                                xx_sl, yy_sl,
                                                                nuc_li, nuc_ri, nuc_qp1, nuc_jm1,
                                                                nuc_i1, nuc_pm1,
                                                                xo, yy_eff, nuc_lo2, nuc_ro_eff,
                                                                loop_e2, cai_yy,
                                                                slt.sl_left_pos, slt.outer_left2,
                                                                slt.aa_o2, slt.ii_o2, slt.paa_o2,
                                                                slt.seg_len_left);
                                                }
                                            }
                                            continue;
                                        }

                                        // Slow path: variant lookup depends on xx_o2.
                                        for (int xx_o2 = 0; xx_o2 < slt.ncod_o2; ++xx_o2) {
                                            const int nuc_lo2 = nucleotides[slt.paa_o2][xx_o2][slt.ii_o2];
                                            if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                            double loop_e2;
                                            if (use_loop_table_c) {
                                                auto& st = loop_e2_table_state_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                                if (st == 0) {
                                                    double le = score_single_loop_with_mismatch(lambda,
                                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                        slt.seg_len_left, seg_len_right,
                                                        nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                    if (le < inf) {
                                                        loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1] = le;
                                                        st = 1;
                                                    } else { st = -1; }
                                                }
                                                if (st < 0) continue;
                                                loop_e2 = loop_e2_table_c[nuc_lo2][sll_idx][nuc_i1][nuc_pm1];
                                            } else {
                                                loop_e2 = score_single_loop_with_mismatch(lambda,
                                                    nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                    slt.seg_len_left, seg_len_right,
                                                    nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                                                if (loop_e2 >= inf) continue;
                                            }
                                            cnt_in_local_c[tid_c]++;

                                            // Pick best (c, sr, sl) for this xx_o2.
                                            BestCSR il_pick{AGG_INF, AggBest{AGG_INF,0,0,-1,false}, AggBest{AGG_INF,0,0,-1,false}, false};
                                            AggBest il_sl{AGG_INF, 0, 0, -1, false};
                                            if (slt.same_codon_sl) {
                                                // sl.y == cx; if xx_o2_locked, sl.x == xx_o2 → use slq.by_xy[xx_o2][cx].
                                                for (int cx = 0; cx < 6; ++cx) {
                                                    if (!best_csr_by_cx[cx].valid) continue;
                                                    const AggBest sl_pick = slt.xx_o2_locked
                                                        ? slq.by_xy[xx_o2][cx]
                                                        : slq.by_y[cx];
                                                    if (!sl_pick.valid) continue;
                                                    double total = best_csr_by_cx[cx].score + sl_pick.score;
                                                    if (total < il_pick.score) {
                                                        il_pick = {total, best_csr_by_cx[cx].c_best,
                                                                   best_csr_by_cx[cx].sr_best, true};
                                                        il_sl = sl_pick;
                                                    }
                                                }
                                            } else {
                                                // !same_codon_sl: SL fully free (mod xx_o2_locked).
                                                AggBest sl_pick = slt.xx_o2_locked
                                                    ? slq.by_x[xx_o2]
                                                    : slq.all;
                                                if (sl_pick.valid && best_csr_all.valid) {
                                                    il_pick = {best_csr_all.score + sl_pick.score,
                                                               best_csr_all.c_best,
                                                               best_csr_all.sr_best, true};
                                                    il_sl = sl_pick;
                                                }
                                            }
                                            if (!il_pick.valid) continue;
                                            const AggBest& c_b = il_pick.c_best;
                                            const AggBest& sr_b = il_pick.sr_best;
                                            const int xx_inner = c_b.cell_idx / NCOD;
                                            const int yy_inner = c_b.cell_idx % NCOD;
                                            const int xx_sr = sr_b.cell_idx / NCOD;
                                            const int yy_sr = sr_b.cell_idx % NCOD;
                                            const int xx_sl = il_sl.cell_idx / NCOD;
                                            const int yy_sl = il_sl.cell_idx % NCOD;
                                            emit_il_agg(c_b, sr_b, il_sl,
                                                        xx_inner, yy_inner,
                                                        xx_sr, yy_sr,
                                                        xx_sl, yy_sl,
                                                        nuc_li, nuc_ri, nuc_qp1, nuc_jm1,
                                                        nuc_i1, nuc_pm1,
                                                        xx_o2, yy_eff, nuc_lo2, nuc_ro_eff,
                                                        loop_e2, cai_yy,
                                                        slt.sl_left_pos, slt.outer_left2,
                                                        slt.aa_o2, slt.ii_o2, slt.paa_o2,
                                                        slt.seg_len_left);
                                        }  // end xx_o2
                                    }  // end sl quadrants
                                }  // end sl_tasks_agg
                            }  // end yy_eff
                        }  // end sr_qid
                    }  // end c_qid
                }  // end inner_left (c_inner_agg_map)
            }  // end OpenMP for over agg_tasks

            // -------- Merge phase (mirrors baseline Block 1). --------
            {
                BeamEntry scratch;
                auto replay_entry = [&](int orig_key, const BeamEntry& src) {
                    for (const auto& v : src.variants) {
                        scratch = src;
                        scratch.x = v.x; scratch.y = v.y;
                        scratch.mfe = v.mfe; scratch.cai = v.cai; scratch.score = v.score;
                        scratch.backtrace_type = v.backtrace_type;
                        scratch.last_closed_nuc = v.last_closed_nuc;
                        scratch.bt_info = v.bt_info;
                        scratch.cs_inner_key    = v.cs_inner_key;
                        scratch.cs_right_s_key  = v.cs_right_s_key;
                        scratch.cs_right_len    = v.cs_right_len;
                        scratch.cs_single_start = v.cs_single_start;
                        scratch.cs_inner_left   = v.cs_inner_left;
                        scratch.cs_inner_right  = v.cs_inner_right;
                        scratch.cs_pack_outer   = v.cs_pack_outer;
                        scratch.variants.clear();
                        update_derna(curr_c, orig_key, v.score, scratch);
                    }
                };
                for (int t = 0; t < num_threads_blk1c; ++t) {
                    DenseBeamTable& td = blk1_thread_dense_maps[t];
                    const auto& order = blk1_thread_key_order[t];
                    for (int compact_idx : order) {
                        const BeamEntry& src = td.slot(compact_idx);
                        if (src.score >= inf) continue;
                        const int left_pos = compact_idx >> 4;
                        const int low      = compact_idx & 0xF;
                        const int orig_key = (left_pos * nuc_len + pos) * 16 + low;
                        replay_entry(orig_key, src);
                    }
                }
                for (int t = 0; t < num_threads_blk1c; ++t) {
                    cnt_blk1_rightbulge += cnt_rb_local_c[t];
                    cnt_blk1_internal  += cnt_in_local_c[t];
                }
            }
        }  // end new closure
#else  // !USE_CODON_LOOKUP_CLOSURE — original Block 1 path
        if (pos >= 5) {  // Block 1 direct right-bulge + left-bulge (matches LCDSfold S_CtoC + CStoC behavior)
            // --- Opt A: Precompute beam_threshold for early-exit score bound ---
            // Use the beamsize-th best score (the cutoff for pruning) as threshold.
            // Any candidate worse than this cannot survive the post-Block 1 prune.
            double blk1_beam_threshold = inf;
            if ((int)curr_c.size() >= beamsize && beamsize > 0) {
                std::vector<double> c_scores;
                c_scores.reserve(curr_c.size());
                for (const auto& kv : curr_c) c_scores.push_back(kv.second.score);
                auto nth = c_scores.begin() + (beamsize - 1);
                std::nth_element(c_scores.begin(), nth, c_scores.end());
                blk1_beam_threshold = *nth;
            }

            // Opt A: S_left score lower bound (0 is a safe loose bound).
            static constexpr double blk1_global_best_sl = 0.0;

            // --- Opt A: min_loop_energy (most negative possible loop_e from score_single_loop_with_mismatch) ---
            // Stacking minimum ≈ -340 cKcal (only relevant for n1=n2=0, which Block 1 doesn't
            // hit). For Block 1's geometry: int11_37=-230, int22_37=-260, plus mismatch slack.
            // Tightened tests (-310, -260) gave NO measurable runtime improvement — the bound
            // prune doesn't bind often. Kept loose at -340 for safety margin.
            static const double blk1_min_loop_e = -340.0;

            // ---- Parallelization setup ----
            // Build a flat list of (seg_len_right, s_right_key, s_right_ent*) tasks so we can
            // parallelize over all of them with good load balance. Iterations across tasks are
            // independent writers into curr_c; we buffer writes per task and merge them in task
            // order after the parallel region — this preserves bit-identity with the serial
            // execution because update_derna's merge is deterministic under fixed insertion order.
            struct Blk1Task {
                int seg_len_right;
                int closed_right;
                int seg_start;
                const BeamEntry* s_right_ent;
                int s_right_key;
                const std::vector<std::pair<int, const BeamEntry*>>* c_flat_ptr;
            };
            std::vector<Blk1Task> blk1_tasks;
            blk1_tasks.reserve(64);
            for (int seg_len_right = 1; seg_len_right <= min(SINGLE_MAX_LEN, pos_prev); ++seg_len_right) {
                if (seg_len_right >= (int)tab_s[pos_prev].size()) continue;
                const auto& s_right_map_outer = tab_s[pos_prev][seg_len_right];
                if (s_right_map_outer.empty()) continue;
                int seg_start_o = pos_prev - seg_len_right + 1;
                int closed_right_o = seg_start_o - 1;
                if (closed_right_o < 0 || closed_right_o >= nuc_len) continue;
                const auto& c_flat_prewarm = get_c_flat(closed_right_o);
                if (c_flat_prewarm.empty()) continue;
                for (const auto& s_rkv : s_right_map_outer) {
                    blk1_tasks.push_back({seg_len_right, closed_right_o, seg_start_o,
                                          &s_rkv.second, s_rkv.first, &c_flat_prewarm});
                }
            }

            // Selective Block 1 (Approach c): only run Block 1 for tuples whose
            // CS entry was DROPPED by the prune at pos_prev. Tuples whose CS
            // entry survived (or which never produced a CS entry due to geometry)
            // are skipped. Gated by env var DERNA_SELECTIVE_BLK1=1.
            static const bool selective_blk1 = []() {
                const char* e = std::getenv("DERNA_SELECTIVE_BLK1");
                return e && e[0] == '1';
            }();
            // DERNA_BLK1_STUB=1 disables Block 1 entirely (for testing whether
            // the CS path alone is lossless).
            static const bool blk1_stub = []() {
                const char* e = std::getenv("DERNA_BLK1_STUB");
                return e && e[0] == '1';
            }();
            const std::unordered_set<uint64_t>* dropped_units_ptr = nullptr;
            if (selective_blk1 && pos_prev >= 0 && pos_prev < (int)tab_cs_dropped_units.size()
                && tab_cs_selective_active[pos_prev]) {
                dropped_units_ptr = &tab_cs_dropped_units[pos_prev];
            }
            if (blk1_stub) blk1_tasks.clear();

            if (!blk1_tasks.empty()) {
                const int num_tasks = (int)blk1_tasks.size();
#ifdef _OPENMP
                const int num_threads_blk1 = omp_get_max_threads();
#else
                const int num_threads_blk1 = 1;
#endif
                // Per-thread local DernaBeamMap + insertion-order log. The thread performs
                // its dedup (variant merge) locally via update_derna — avoiding the cost
                // of pushing every candidate into a raw buffer. After the parallel region
                // we replay each thread's entries (in its insertion order, in thread-id
                // order) into curr_c so the final result is deterministic and bit-identical.
                //
                // Default path uses blk1_thread_dense_maps (function-scope hoisted dense
                // vectors, reset per pos). BLK1_DENSE=0 falls back to the prior
                // DernaBeamMap path for A/B comparison.
                std::vector<DernaBeamMap> thread_maps(blk1_use_dense ? 0 : num_threads_blk1);
                std::vector<std::vector<int>> thread_key_order_legacy(blk1_use_dense ? 0 : num_threads_blk1);
                if (blk1_use_dense) {
                    // Reset only the populated slots from the previous pos. data_ stays
                    // allocated across positions; capacity grows on demand.
                    const int dense_cap = nuc_len * 16;
                    for (int t = 0; t < num_threads_blk1_outer; ++t) {
                        blk1_thread_dense_maps[t].reset(dense_cap);
                        blk1_thread_key_order[t].clear();
                    }
                }
                // BLK1_AGG_DIAG: per-thread per-(idx,x,y) hit-count vector, sized
                // and zeroed at the start of every pos. This vector is thread-local
                // and accessed from inside DenseBeamTable::update; we use a small
                // parallel section here to reset on each thread.
                if (blk1_agg_diag_enabled() && blk1_use_dense) {
                    const size_t need = (size_t)nuc_len * 16ull * 36ull;
                    #pragma omp parallel
                    {
                        if (g_blk1_agg_hits_per_xy.size() < need) {
                            g_blk1_agg_hits_per_xy.assign(need, 0);
                        } else {
                            std::fill(g_blk1_agg_hits_per_xy.begin(),
                                      g_blk1_agg_hits_per_xy.begin() + need, 0);
                        }
                    }
                }
                std::vector<uint64_t> cnt_rb_local(num_threads_blk1, 0);
                std::vector<uint64_t> cnt_in_local(num_threads_blk1, 0);

                // Per-thread SLTask struct definition with pre-computed variant lookup tables.
                struct BestVar { double score; const XYVariant* ptr; };
                static constexpr double BV_INF = 1e30;
                struct SLTask {
                    const BeamEntry* s_left_ent;
                    int s_left_key;
                    int seg_len_left;
                    int outer_left2;
                    int aa_o2, ii_o2, paa_o2, ncod_o2;
                    bool xx_o2_locked;
                    bool same_codon_sl;
                    int nuc_i1, nuc_pm1;
                    std::array<BestVar, 6> sl_by_x;
                    std::array<BestVar, 6> sl_by_y;
                    std::array<std::array<BestVar, 6>, 6> sl_by_xy;
                    BestVar sl_all;
                    // Precomputed: for each nuc_lo2 value (0..3), the best xx_o2 (by cai).
                    // best_xx_per_nuc[nuc].xx = best codon index, .cai = its cai contribution.
                    struct NucBest { int xx; double cai; bool valid; };
                    std::array<NucBest, 4> best_xx_per_nuc;
                };
                std::vector<std::vector<SLTask>> thread_sl_tasks(num_threads_blk1);

                #pragma omp parallel for schedule(dynamic, 1)
                for (int task_idx = 0; task_idx < num_tasks; ++task_idx) {
                    const Blk1Task& tk = blk1_tasks[task_idx];
                    const int seg_len_right = tk.seg_len_right;
                    const int closed_right = tk.closed_right;
                    const int seg_start = tk.seg_start;
                    const int s_right_key = tk.s_right_key;
                    const BeamEntry& s_right_ent = *tk.s_right_ent;
#ifdef _OPENMP
                    const int tid_blk1 = omp_get_thread_num();
#else
                    const int tid_blk1 = 0;
#endif
                    DernaBeamMap* task_map_ptr = blk1_use_dense ? nullptr : &thread_maps[tid_blk1];
                    DenseBeamTable* task_dense_ptr = blk1_use_dense ? &blk1_thread_dense_maps[tid_blk1] : nullptr;
                    std::vector<int>& task_order = blk1_use_dense ? blk1_thread_key_order[tid_blk1] : thread_key_order_legacy[tid_blk1];
                    const auto& c_flat_blk1 = *tk.c_flat_ptr;
                    if (sigma(s_right_ent.b, s_right_ent.j) != pos_prev) continue;
                    if (sigma(s_right_ent.a, s_right_ent.i) != seg_start) continue;
                    bool same_codon_sr = false;
                    if (!same_or_next_codon(s_right_ent.b, s_right_ent.j, bb, slot, same_codon_sr)) continue;

                    const bool seam_sr_applies_blk1 = ((closed_right % 3) != 2);

                    // Pre-compute best S_right variant score by x, by y, by (x,y),
                    // and overall best. Key insight: all variants of the same entry share
                    // identical boundary nucleotides (encoded in the key), so loop_e is
                    // constant across variants — only the variant scores differ.
                    std::array<BestVar, 6> best_sr_by_x;
                    std::array<BestVar, 6> best_sr_by_y;
                    std::array<std::array<BestVar, 6>, 6> best_sr_by_xy;
                    for (auto& bv : best_sr_by_x) bv = {BV_INF, nullptr};
                    for (auto& bv : best_sr_by_y) bv = {BV_INF, nullptr};
                    for (auto& row : best_sr_by_xy) for (auto& bv : row) bv = {BV_INF, nullptr};
                    BestVar best_sr_all = {BV_INF, nullptr};
                    for (const auto& sv : s_right_ent.variants) {
                        if (sv.x >= 0 && sv.x < 6 && sv.score < best_sr_by_x[sv.x].score)
                            best_sr_by_x[sv.x] = {sv.score, &sv};
                        if (sv.y >= 0 && sv.y < 6 && sv.score < best_sr_by_y[sv.y].score)
                            best_sr_by_y[sv.y] = {sv.score, &sv};
                        if (sv.x >= 0 && sv.x < 6 && sv.y >= 0 && sv.y < 6 &&
                            sv.score < best_sr_by_xy[sv.x][sv.y].score)
                            best_sr_by_xy[sv.x][sv.y] = {sv.score, &sv};
                        if (sv.score < best_sr_all.score)
                            best_sr_all = {sv.score, &sv};
                    }

                    // Extract S_right boundary nucs from key (constant across all srvars).
                    const int nuc_qp1_blk1 = (s_right_key & 0xF) >> 2;  // left boundary nuc
                    const int nuc_jm1_blk1 = s_right_key & 3;            // right boundary nuc

                    for (const auto& c_kv : c_flat_blk1) {
                        const int c_key = c_kv.first;
                        const BeamEntry& c_ent = *c_kv.second;
                        if (sigma(c_ent.b, c_ent.j) != closed_right) continue;

                        // Selective Block 1: skip (c_key, s_right_key) tuples whose
                        // CS entry was NOT dropped by the prune at pos_prev.
                        // Surviving CS tuples are covered by CS_to_C / S_CS_to_C
                        // closures (with iv-wrap completing the inner-C variant
                        // span). Non-producer tuples (no CS entry at all) are also
                        // safely skipped — they have no Block 1 work.
                        if (dropped_units_ptr) {
                            const uint64_t unit_pack = pack_cs_blk1_unit(
                                closed_right, seg_len_right, s_right_key, c_key);
                            if (dropped_units_ptr->count(unit_pack) == 0) continue;
                        }

                        // Seam consistency between C right end and S right start.
                        if ((closed_right % 3) != 2) {
                            if (s_right_ent.a != c_ent.b) continue;
                            if (s_right_ent.i != c_ent.j + 1) continue;
                        } else {
                            if (s_right_ent.a != c_ent.b + 1) continue;
                            if (s_right_ent.i != 0) continue;
                        }

                        int inner_left = sigma(c_ent.a, c_ent.i);
                        int outer_left_nuc = inner_left - 1;  // right-bulge: outer pair left
                        if (outer_left_nuc < 0) continue;
                        if ((pos - outer_left_nuc) <= (HAIRPIN_GAP + 2)) continue;

                        // Opt A: Early-exit score bound.
                        // Right-bulge: total >= c + sr + min_loop_e
                        // Internal:    total >= c + sr + best_sl + min_loop_e
                        // If both bounds exceed beam threshold, skip this pair entirely.
                        // Opt A: Early-exit score bound.
                        const double blk1_base_bound = c_ent.score + s_right_ent.score + blk1_min_loop_e;
                        if (blk1_base_bound + blk1_global_best_sl > blk1_beam_threshold) continue;

                        int aa_o = outer_left_nuc / 3, ii_o = outer_left_nuc % 3;
                        int paa_o = protein[aa_o];
                        const int ncod_o = n_codon[paa_o];
                        const bool pos_is_last = is_last_nuc(pos);
                        const bool xx_o_locked = (aa_o == c_ent.a);

                        // Extract inner C boundary nucs from key (constant across all cvars).
                        const int nuc_li = (c_key & 0xF) >> 2;
                        const int nuc_ri = c_key & 3;

                        // Pre-compute best C variant by x, by y, and by (x,y) for separable lookup.
                        std::array<BestVar, 6> best_c_by_x;
                        std::array<BestVar, 6> best_c_by_y;
                        std::array<std::array<BestVar, 6>, 6> best_c_by_xy;
                        for (auto& bv : best_c_by_x) bv = {BV_INF, nullptr};
                        for (auto& bv : best_c_by_y) bv = {BV_INF, nullptr};
                        for (auto& row : best_c_by_xy) for (auto& bv : row) bv = {BV_INF, nullptr};
                        BestVar best_c_all = {BV_INF, nullptr};
                        for (const auto& cv : c_ent.variants) {
                            if (cv.x >= 0 && cv.x < 6 && cv.score < best_c_by_x[cv.x].score)
                                best_c_by_x[cv.x] = {cv.score, &cv};
                            if (cv.y >= 0 && cv.y < 6 && cv.score < best_c_by_y[cv.y].score)
                                best_c_by_y[cv.y] = {cv.score, &cv};
                            if (cv.x >= 0 && cv.x < 6 && cv.y >= 0 && cv.y < 6 &&
                                cv.score < best_c_by_xy[cv.x][cv.y].score)
                                best_c_by_xy[cv.x][cv.y] = {cv.score, &cv};
                            if (cv.score < best_c_all.score)
                                best_c_all = {cv.score, &cv};
                        }

                        // Opt B: Precompute best xx_o per nucleotide for right-bulge outer left.
                        // When free (!xx_o_locked): pick the best-CAI codon per nucleotide,
                        // since loop_e depends only on nuc and variant lookup is xx_o-independent.
                        // When locked (xx_o_locked): iterate all valid codons since each has
                        // a different csr score (best_csr_by_cx[xx_o]). Still iterate by nuc
                        // to avoid redundant loop_e computations.
                        typedef SLTask::NucBest NucBest;
                        std::array<NucBest, 4> best_xxo_per_nuc;
                        for (auto& nb : best_xxo_per_nuc) nb = {-1, 0.0, false};
                        if (!xx_o_locked) {
                            for (int xo = 0; xo < ncod_o; ++xo) {
                                int nlo = nucleotides[paa_o][xo][ii_o];
                                double xcai = (ii_o == 2) ? codon_cai[paa_o][xo] : 0.0;
                                if (!best_xxo_per_nuc[nlo].valid || xcai > best_xxo_per_nuc[nlo].cai) {
                                    best_xxo_per_nuc[nlo] = {xo, xcai, true};
                                }
                            }
                        }

                        // Pre-build sl_tasks ONCE per c_ent (hoisted out of variant loops).
                        auto& sl_tasks = thread_sl_tasks[tid_blk1];
                        sl_tasks.clear();
                        {
                            int seam_pos_left = outer_left_nuc;
                            if (seam_pos_left >= 1 && seam_pos_left < nuc_len) {
                                int max_sl = min(SINGLE_MAX_LEN, seam_pos_left + 1);
                                for (int seg_len_left = 1; seg_len_left <= max_sl; ++seg_len_left) {
                                    if (seg_len_left >= (int)tab_s[seam_pos_left].size()) continue;
                                    const auto& s_left_map = tab_s[seam_pos_left][seg_len_left];
                                    if (s_left_map.empty()) continue;
                                    for (const auto& s_lkv : s_left_map) {
                                        const BeamEntry& s_left_ent = s_lkv.second;
                                        if (sigma(s_left_ent.b, s_left_ent.j) != seam_pos_left) continue;
                                        bool same_codon_sl = false;
                                        if (!same_or_next_codon(s_left_ent.b, s_left_ent.j, c_ent.a, c_ent.i, same_codon_sl)) continue;
                                        int outer_left2 = sigma(s_left_ent.a, s_left_ent.i) - 1;
                                        if (outer_left2 < 0) continue;
                                        if ((pos - outer_left2) <= (HAIRPIN_GAP + 2)) continue;
                                        int aa_o2 = outer_left2 / 3, ii_o2 = outer_left2 % 3;
                                        int paa_o2 = protein[aa_o2];
                                        const int ncod_o2 = n_codon[paa_o2];
                                        const bool xx_o2_locked = (aa_o2 == s_left_ent.a);
                                        const int sl_key = s_lkv.first;
                                        SLTask slt_new;
                                        slt_new.s_left_ent = &s_left_ent;
                                        slt_new.s_left_key = sl_key;
                                        slt_new.seg_len_left = seg_len_left;
                                        slt_new.outer_left2 = outer_left2;
                                        slt_new.aa_o2 = aa_o2; slt_new.ii_o2 = ii_o2;
                                        slt_new.paa_o2 = paa_o2; slt_new.ncod_o2 = ncod_o2;
                                        slt_new.xx_o2_locked = xx_o2_locked;
                                        slt_new.same_codon_sl = same_codon_sl;
                                        slt_new.nuc_i1 = (sl_key & 0xF) >> 2;
                                        slt_new.nuc_pm1 = sl_key & 3;
                                        for (auto& bv : slt_new.sl_by_x) bv = {BV_INF, nullptr};
                                        for (auto& bv : slt_new.sl_by_y) bv = {BV_INF, nullptr};
                                        for (auto& row : slt_new.sl_by_xy) for (auto& bv : row) bv = {BV_INF, nullptr};
                                        slt_new.sl_all = {BV_INF, nullptr};
                                        for (const auto& sv : s_left_ent.variants) {
                                            if (sv.x >= 0 && sv.x < 6 && sv.score < slt_new.sl_by_x[sv.x].score)
                                                slt_new.sl_by_x[sv.x] = {sv.score, &sv};
                                            if (sv.y >= 0 && sv.y < 6 && sv.score < slt_new.sl_by_y[sv.y].score)
                                                slt_new.sl_by_y[sv.y] = {sv.score, &sv};
                                            if (sv.x >= 0 && sv.x < 6 && sv.y >= 0 && sv.y < 6 &&
                                                sv.score < slt_new.sl_by_xy[sv.x][sv.y].score)
                                                slt_new.sl_by_xy[sv.x][sv.y] = {sv.score, &sv};
                                            if (sv.score < slt_new.sl_all.score)
                                                slt_new.sl_all = {sv.score, &sv};
                                        }
                                        // Precompute best xx_o2 per nuc_lo2 value for fast
                                        // inner loop when variant lookup is xx_o2-independent.
                                        for (auto& nb : slt_new.best_xx_per_nuc) nb = {-1, 0.0, false};
                                        for (int xo = 0; xo < ncod_o2; ++xo) {
                                            int nlo = nucleotides[paa_o2][xo][ii_o2];
                                            double xcai = (ii_o2 == 2) ? codon_cai[paa_o2][xo] : 0.0;
                                            if (!slt_new.best_xx_per_nuc[nlo].valid || xcai > slt_new.best_xx_per_nuc[nlo].cai) {
                                                slt_new.best_xx_per_nuc[nlo] = {xo, xcai, true};
                                            }
                                        }
                                        sl_tasks.push_back(std::move(slt_new));
                                    }
                                }
                            }
                        }

                        // --- Separable best-variant lookup for Block 1 ---
                        // loop_e depends only on key-derived boundary nucs (nuc_li, nuc_ri,
                        // nuc_qp1, nuc_jm1 for right-bulge; + nuc_i1, nuc_pm1 for internal)
                        // plus the outer pair nucs (nuc_lo, nuc_ro from xx_o/yy).
                        // score = cvar.score + srvar.score [+ slvar.score] + lambda*loop_e + const
                        // → find best seam-compatible variants in O(8) instead of O(512).

                        // Lambda to emit a right-bulge candidate for a given (cvar, srvar, xx_o, yy).
                        // OPT1: thread-local scratch BeamEntry reused across emit calls.
                        // Default-constructed once per Block 1 task; the per-call fields are
                        // overwritten in the lambdas. Saves ~120 B of default-init per call
                        // (variant_idx[6][6] + cs_* fields + last_closed_nuc), and every
                        // emit_rb/emit_il invocation now skips the BeamEntry constructor.
                        // Fields not touched by these lambdas (last_closed_nuc, cs_inner_key,
                        // cs_right_s_key, cs_pack_outer) match the constructor defaults
                        // (-1 / -1LL), which matches the prior behavior.
                        BeamEntry scratch_rb;
                        BeamEntry scratch_il;

                        auto emit_rb = [&](const XYVariant& cvar, const XYVariant& srvar,
                                           int xx_o, int yy_eff, int nuc_lo, int nuc_ro_eff,
                                           double loop_e, double cai_yy) {
                            double mfe_cand = cvar.mfe + srvar.mfe + loop_e;
                            double cai_cand = cvar.cai + srvar.cai + cai_yy
                                + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                            double sc_cand = combined_score(lambda, mfe_cand, cai_cand);
                            int key_cl = nuc_key_c_nuc(outer_left_nuc, pos, nuc_lo, nuc_ro_eff, nuc_len);
                            // BLK1_AGG: early-rejection (same logic as emit_il above).
                            if (blk1_agg_enabled && task_dense_ptr) {
                                const int compact_idx_pre = outer_left_nuc * 16 + nuc_lo * 4 + nuc_ro_eff;
                                const BeamEntry& ex_pre = task_dense_ptr->slot(compact_idx_pre);
                                if (ex_pre.score < inf &&
                                    xx_o >= 0 && xx_o < 6 && yy_eff >= 0 && yy_eff < 6) {
                                    const int8_t vi_pre = ex_pre.variant_idx[xx_o][yy_eff];
                                    if (vi_pre >= 0 && sc_cand > ex_pre.variants[vi_pre].score) {
                                        return;
                                    }
                                }
                            }
                            scratch_rb.score = sc_cand;
                            scratch_rb.mfe = mfe_cand;
                            scratch_rb.cai = cai_cand;
                            scratch_rb.a = aa_o; scratch_rb.b = bb;
                            scratch_rb.i = ii_o; scratch_rb.j = slot;
                            scratch_rb.x = xx_o; scratch_rb.y = yy_eff;
                            scratch_rb.backtrace_type = (int)DernaManner::MANNER_C_StoC;
                            scratch_rb.bt_info.len = 6;
                            scratch_rb.bt_info.data[0] = c_key;
                            scratch_rb.bt_info.data[1] = s_right_key;
                            scratch_rb.bt_info.data[2] = cvar.x;
                            scratch_rb.bt_info.data[3] = cvar.y;
                            scratch_rb.bt_info.data[4] = srvar.x;
                            scratch_rb.bt_info.data[5] = srvar.y;
                            scratch_rb.cs_inner_left = closed_right;
                            scratch_rb.cs_right_len = seg_len_right;
                            check_c_entry_invariant(lambda, scratch_rb, "Blk1_C_StoC", pos);
                            if (task_dense_ptr) {
                                const int compact_idx = outer_left_nuc * 16 + nuc_lo * 4 + nuc_ro_eff;
                                if (task_dense_ptr->update(compact_idx, sc_cand, scratch_rb))
                                    task_order.push_back(compact_idx);
                            } else {
                                if (update_derna(*task_map_ptr, key_cl, sc_cand, scratch_rb)) task_order.push_back(key_cl);
                            }
                        };

                        // Lambda to emit an internal loop candidate.
                        auto emit_il = [&](const XYVariant& slvar, const XYVariant& cvar,
                                           const XYVariant& srvar, int s_left_key,
                                           int xx_o2, int yy_eff, int nuc_lo2, int nuc_ro_eff,
                                           double loop_e2, double cai_yy,
                                           int outer_left2, int aa_o2, int ii_o2, int paa_o2,
                                           int seg_len_left) {
                            double mfe2 = cvar.mfe + srvar.mfe + slvar.mfe + loop_e2;
                            double cai2 = cvar.cai + srvar.cai + slvar.cai + cai_yy
                                + (ii_o2 == 2 ? codon_cai[paa_o2][xx_o2] : 0.0);
                            double sc2 = combined_score(lambda, mfe2, cai2);
                            int key2 = nuc_key_c_nuc(outer_left2, pos, nuc_lo2, nuc_ro_eff, nuc_len);
                            // BLK1_AGG: early-rejection. If the dense slot already holds
                            // a (xx_o2, yy_eff) variant with strictly-better score, skip
                            // the scratch_il build + update_derna call. The variant-merge
                            // path would have early-returned on this anyway (99.83% of
                            // calls per BLK1_AGG_DIAG measurement on P01707).
                            if (blk1_agg_enabled && task_dense_ptr) {
                                const int compact_idx2_pre = outer_left2 * 16 + nuc_lo2 * 4 + nuc_ro_eff;
                                const BeamEntry& ex_pre = task_dense_ptr->slot(compact_idx2_pre);
                                if (ex_pre.score < inf &&
                                    xx_o2 >= 0 && xx_o2 < 6 && yy_eff >= 0 && yy_eff < 6) {
                                    const int8_t vi_pre = ex_pre.variant_idx[xx_o2][yy_eff];
                                    if (vi_pre >= 0 && sc2 > ex_pre.variants[vi_pre].score) {
                                        return;  // skip — full update path would have done so too
                                    }
                                }
                            }
                            scratch_il.score = sc2;
                            scratch_il.mfe = mfe2;
                            scratch_il.cai = cai2;
                            scratch_il.a = aa_o2; scratch_il.b = bb;
                            scratch_il.i = ii_o2; scratch_il.j = slot;
                            scratch_il.x = xx_o2; scratch_il.y = yy_eff;
                            scratch_il.backtrace_type = (int)DernaManner::MANNER_S_C_StoC;
                            scratch_il.bt_info.len = 9;
                            scratch_il.bt_info.data[0] = s_left_key;
                            scratch_il.bt_info.data[1] = c_key;
                            scratch_il.bt_info.data[2] = s_right_key;
                            scratch_il.bt_info.data[3] = slvar.x;
                            scratch_il.bt_info.data[4] = slvar.y;
                            scratch_il.bt_info.data[5] = cvar.x;
                            scratch_il.bt_info.data[6] = cvar.y;
                            scratch_il.bt_info.data[7] = srvar.x;
                            scratch_il.bt_info.data[8] = srvar.y;
                            scratch_il.cs_inner_left = closed_right;
                            scratch_il.cs_right_len = seg_len_right;
                            scratch_il.cs_single_start = outer_left_nuc;
                            scratch_il.cs_inner_right = seg_len_left;
                            check_c_entry_invariant(lambda, scratch_il, "Blk1_S_C_StoC", pos);
                            if (task_dense_ptr) {
                                const int compact_idx2 = outer_left2 * 16 + nuc_lo2 * 4 + nuc_ro_eff;
                                if (task_dense_ptr->update(compact_idx2, sc2, scratch_il))
                                    task_order.push_back(compact_idx2);
                            } else {
                                if (update_derna(*task_map_ptr, key2, sc2, scratch_il)) task_order.push_back(key2);
                            }
                        };

                        // Iterate yy values, then (xx_o) for right-bulge and
                        // (sl_task, xx_o2) for internal loops. Variant selection is always
                        // via O(1) lookup tables.
                        //
                        // When same_codon_sr: yy must equal srvar.y, so we iterate distinct
                        // yy values (0..ncod_bb-1) and require best_sr_by_y[yy] to exist.
                        // When !same_codon_sr: yy iterates freely (0..ncod_bb-1).
                        {
                            for (int yy_idx = 0; yy_idx < ncod_bb; ++yy_idx) {
                                const int yy_eff = yy_idx;
                                if (same_codon_sr && !best_sr_by_y[yy_eff].ptr) continue;
                                const int nuc_ro_eff = nucleotides[pbb][yy_eff][slot];
                                const double cai_yy = pos_is_last ? codon_cai[pbb][yy_eff] : 0.0;

                                // Pre-compute combined (c + sr) best by cx for this yy_eff.
                                // Used for both right-bulge and internal loop lookups.
                                struct BestCSR { double score; const XYVariant* c_ptr; const XYVariant* sr_ptr; };
                                std::array<BestCSR, 6> best_csr_by_cx;
                                BestCSR best_csr_all = {BV_INF, nullptr, nullptr};
                                for (auto& b : best_csr_by_cx) b = {BV_INF, nullptr, nullptr};

                                if (!seam_sr_applies_blk1 && !same_codon_sr) {
                                    if (best_sr_all.ptr) {
                                        for (int cx = 0; cx < 6; ++cx) {
                                            if (best_c_by_x[cx].ptr) {
                                                double sc = best_c_by_x[cx].score + best_sr_all.score;
                                                best_csr_by_cx[cx] = {sc, best_c_by_x[cx].ptr, best_sr_all.ptr};
                                                if (sc < best_csr_all.score) best_csr_all = best_csr_by_cx[cx];
                                            }
                                        }
                                    }
                                } else if (!seam_sr_applies_blk1) {
                                    if (best_sr_by_y[yy_eff].ptr) {
                                        double sr_sc = best_sr_by_y[yy_eff].score;
                                        for (int cx = 0; cx < 6; ++cx) {
                                            if (best_c_by_x[cx].ptr) {
                                                double sc = best_c_by_x[cx].score + sr_sc;
                                                best_csr_by_cx[cx] = {sc, best_c_by_x[cx].ptr, best_sr_by_y[yy_eff].ptr};
                                                if (sc < best_csr_all.score) best_csr_all = best_csr_by_cx[cx];
                                            }
                                        }
                                    }
                                } else if (!same_codon_sr) {
                                    for (int cx = 0; cx < 6; ++cx) {
                                        for (int cy = 0; cy < 6; ++cy) {
                                            const auto& bc = best_c_by_xy[cx][cy];
                                            if (!bc.ptr) continue;
                                            const auto& bsr = best_sr_by_x[cy];
                                            if (!bsr.ptr) continue;
                                            double sc = bc.score + bsr.score;
                                            if (sc < best_csr_by_cx[cx].score)
                                                best_csr_by_cx[cx] = {sc, bc.ptr, bsr.ptr};
                                        }
                                        if (best_csr_by_cx[cx].score < best_csr_all.score)
                                            best_csr_all = best_csr_by_cx[cx];
                                    }
                                } else {
                                    for (int cx = 0; cx < 6; ++cx) {
                                        for (int cy = 0; cy < 6; ++cy) {
                                            const auto& bc = best_c_by_xy[cx][cy];
                                            if (!bc.ptr) continue;
                                            const auto& bsr = best_sr_by_xy[cy][yy_eff];
                                            if (!bsr.ptr) continue;
                                            double sc = bc.score + bsr.score;
                                            if (sc < best_csr_by_cx[cx].score)
                                                best_csr_by_cx[cx] = {sc, bc.ptr, bsr.ptr};
                                        }
                                        if (best_csr_by_cx[cx].score < best_csr_all.score)
                                            best_csr_all = best_csr_by_cx[cx];
                                    }
                                }

                                // ---- (1) Right-bulge (Opt C: group by nuc_lo, compute loop_e once) ----
                                // For each nucleotide nuc_lo, compute loop_e once, then emit
                                // all xx_o codons producing that nucleotide.
                                for (int nuc_lo = 0; nuc_lo < 4; ++nuc_lo) {
                                    if (!best_xxo_per_nuc[nuc_lo].valid && !xx_o_locked) continue;
                                    if (BP_pair[nuc_lo + 1][nuc_ro_eff + 1] == 0) continue;
                                    double loop_e = score_single_loop_with_mismatch(lambda,
                                        nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                        0, seg_len_right,
                                        nuc_li, nuc_jm1_blk1, nuc_lo, nuc_qp1_blk1);
                                    if (loop_e >= inf) continue;
                                    // Emit all codons producing this nucleotide
                                    for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                                        if (nucleotides[paa_o][xx_o][ii_o] != nuc_lo) continue;
                                        cnt_rb_local[tid_blk1]++;
                                        const auto& csr = xx_o_locked ? best_csr_by_cx[xx_o] : best_csr_all;
                                        if (!csr.c_ptr) continue;
                                        emit_rb(*csr.c_ptr, *csr.sr_ptr, xx_o, yy_eff,
                                                nuc_lo, nuc_ro_eff, loop_e, cai_yy);
                                    }
                                }

                                // ---- (2) Internal loop ----
                                // Precompute loop_e table indexed by (nuc_lo2, seg_len_left, nuc_i1, nuc_pm1).
                                // This hoists score_single_loop_with_mismatch out of the sl_tasks loop:
                                // at most 1280 table entries (4*20*4*4) replace millions of redundant calls.
                                // Only built when there are enough sl_tasks to amortize the setup cost.
                                static constexpr int SLL_MAX = SINGLE_MAX_LEN;
                                const bool use_loop_table = ((int)sl_tasks.size() >= 4);
                                double loop_e2_table[4][SLL_MAX][4][4];
                                int8_t loop_e2_table_state[4][SLL_MAX][4][4]; // 0=uncomputed, 1=valid, -1=invalid
                                if (use_loop_table) {
                                    memset(loop_e2_table_state, 0, sizeof(loop_e2_table_state));
                                }

                                for (const SLTask& slt : sl_tasks) {
                                    const int sll_idx = slt.seg_len_left - 1;

                                    // Fast path: variant lookup is completely independent of xx_o2.
                                    if (!slt.same_codon_sl && !slt.xx_o2_locked) {
                                        if (!best_csr_all.c_ptr || !slt.sl_all.ptr) continue;
                                        for (int nuc_lo2 = 0; nuc_lo2 < 4; ++nuc_lo2) {
                                            const auto& nb = slt.best_xx_per_nuc[nuc_lo2];
                                            if (!nb.valid) continue;
                                            if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                            double loop_e2;
                                            if (use_loop_table) {
                                                auto& st = loop_e2_table_state[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                                if (st == 0) {
                                                    double le = score_single_loop_with_mismatch(lambda,
                                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                        slt.seg_len_left, seg_len_right,
                                                        slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                                    if (le < inf) { loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1] = le; st = 1; }
                                                    else { st = -1; }
                                                }
                                                if (st < 0) continue;
                                                loop_e2 = loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                            } else {
                                                loop_e2 = score_single_loop_with_mismatch(lambda,
                                                    nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                    slt.seg_len_left, seg_len_right,
                                                    slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                                if (loop_e2 >= inf) continue;
                                            }
                                            // Emit all codons producing this nucleotide (lossless).
                                            for (int xo = 0; xo < slt.ncod_o2; ++xo) {
                                                if (nucleotides[slt.paa_o2][xo][slt.ii_o2] != nuc_lo2) continue;
                                                cnt_in_local[tid_blk1]++;
                                                emit_il(*slt.sl_all.ptr, *best_csr_all.c_ptr, *best_csr_all.sr_ptr,
                                                        slt.s_left_key, xo, yy_eff, nuc_lo2, nuc_ro_eff,
                                                        loop_e2, cai_yy, slt.outer_left2, slt.aa_o2, slt.ii_o2,
                                                        slt.paa_o2, slt.seg_len_left);
                                            }
                                        }
                                        continue;
                                    }

                                    // Semi-fast path: same_codon_sl && !xx_o2_locked.
                                    // Variant lookup is sl_by_y[cx], independent of xx_o2.
                                    // Use nuc_lo2 iteration (Opt B).
                                    if (slt.same_codon_sl && !slt.xx_o2_locked) {
                                        // Find best (csr, sl) combination across cx values
                                        const XYVariant* best_cv_sf = nullptr;
                                        const XYVariant* best_sr_sf = nullptr;
                                        const XYVariant* best_sl_sf = nullptr;
                                        double best_sf_total = BV_INF;
                                        for (int cx = 0; cx < 6; ++cx) {
                                            const auto& csr = best_csr_by_cx[cx];
                                            if (!csr.c_ptr) continue;
                                            const XYVariant* sl_pick = slt.sl_by_y[cx].ptr;
                                            if (!sl_pick) continue;
                                            double total = csr.score + sl_pick->score;
                                            if (total < best_sf_total) {
                                                best_sf_total = total;
                                                best_cv_sf = csr.c_ptr;
                                                best_sr_sf = csr.sr_ptr;
                                                best_sl_sf = sl_pick;
                                            }
                                        }
                                        if (!best_cv_sf) continue;
                                        for (int nuc_lo2 = 0; nuc_lo2 < 4; ++nuc_lo2) {
                                            const auto& nb = slt.best_xx_per_nuc[nuc_lo2];
                                            if (!nb.valid) continue;
                                            if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                            double loop_e2;
                                            if (use_loop_table) {
                                                auto& st = loop_e2_table_state[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                                if (st == 0) {
                                                    double le = score_single_loop_with_mismatch(lambda,
                                                        nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                        slt.seg_len_left, seg_len_right,
                                                        slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                                    if (le < inf) { loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1] = le; st = 1; }
                                                    else { st = -1; }
                                                }
                                                if (st < 0) continue;
                                                loop_e2 = loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                            } else {
                                                loop_e2 = score_single_loop_with_mismatch(lambda,
                                                    nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                    slt.seg_len_left, seg_len_right,
                                                    slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                                if (loop_e2 >= inf) continue;
                                            }
                                            // Emit all codons producing this nucleotide (lossless).
                                            for (int xo = 0; xo < slt.ncod_o2; ++xo) {
                                                if (nucleotides[slt.paa_o2][xo][slt.ii_o2] != nuc_lo2) continue;
                                                cnt_in_local[tid_blk1]++;
                                                emit_il(*best_sl_sf, *best_cv_sf, *best_sr_sf,
                                                        slt.s_left_key, xo, yy_eff, nuc_lo2, nuc_ro_eff,
                                                        loop_e2, cai_yy, slt.outer_left2, slt.aa_o2, slt.ii_o2,
                                                        slt.paa_o2, slt.seg_len_left);
                                            }
                                        }
                                        continue;
                                    }

                                    // Remaining slow path: variant lookup depends on xx_o2.
                                    if (!slt.same_codon_sl && !best_csr_all.c_ptr) continue;

                                    for (int xx_o2 = 0; xx_o2 < slt.ncod_o2; ++xx_o2) {
                                        int nuc_lo2 = nucleotides[slt.paa_o2][xx_o2][slt.ii_o2];
                                        if (BP_pair[nuc_lo2 + 1][nuc_ro_eff + 1] == 0) continue;
                                        double loop_e2;
                                        if (use_loop_table) {
                                            auto& st = loop_e2_table_state[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                            if (st == 0) {
                                                double le = score_single_loop_with_mismatch(lambda,
                                                    nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                    slt.seg_len_left, seg_len_right,
                                                    slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                                if (le < inf) { loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1] = le; st = 1; }
                                                else { st = -1; }
                                            }
                                            if (st < 0) continue;
                                            loop_e2 = loop_e2_table[nuc_lo2][sll_idx][slt.nuc_i1][slt.nuc_pm1];
                                        } else {
                                            loop_e2 = score_single_loop_with_mismatch(lambda,
                                                nuc_lo2, nuc_ro_eff, nuc_li, nuc_ri,
                                                slt.seg_len_left, seg_len_right,
                                                slt.nuc_i1, nuc_jm1_blk1, slt.nuc_pm1, nuc_qp1_blk1);
                                            if (loop_e2 >= inf) continue;
                                        }
                                        cnt_in_local[tid_blk1]++;

                                        const XYVariant* best_cv_il = nullptr;
                                        const XYVariant* best_sr_il = nullptr;
                                        const XYVariant* best_sl_il = nullptr;
                                        double best_il_total = BV_INF;

                                        if (slt.same_codon_sl) {
                                            for (int cx = 0; cx < 6; ++cx) {
                                                const auto& csr = best_csr_by_cx[cx];
                                                if (!csr.c_ptr) continue;
                                                const XYVariant* sl_pick = slt.xx_o2_locked
                                                    ? slt.sl_by_xy[xx_o2][cx].ptr
                                                    : slt.sl_by_y[cx].ptr;
                                                if (!sl_pick) continue;
                                                double total = csr.score + sl_pick->score;
                                                if (total < best_il_total) {
                                                    best_il_total = total;
                                                    best_cv_il = csr.c_ptr;
                                                    best_sr_il = csr.sr_ptr;
                                                    best_sl_il = sl_pick;
                                                }
                                            }
                                        } else {
                                            const XYVariant* sl_pick = slt.sl_by_x[xx_o2].ptr;
                                            if (sl_pick) {
                                                best_cv_il = best_csr_all.c_ptr;
                                                best_sr_il = best_csr_all.sr_ptr;
                                                best_sl_il = sl_pick;
                                                best_il_total = best_csr_all.score + sl_pick->score;
                                            }
                                        }

                                        if (!best_cv_il || !best_sr_il || !best_sl_il) continue;
                                        if (slt.xx_o2_locked && best_sl_il->x != xx_o2) continue;

                                        emit_il(*best_sl_il, *best_cv_il, *best_sr_il, slt.s_left_key,
                                                xx_o2, yy_eff, nuc_lo2, nuc_ro_eff, loop_e2, cai_yy,
                                                slt.outer_left2, slt.aa_o2, slt.ii_o2, slt.paa_o2, slt.seg_len_left);
                                    }  // end xx_o2 loop
                                }  // end sl_tasks
                            }  // end yy loop
                        }  // end unified same_codon_sr / free-yy branch
                    }  // end c_flat_blk1 loop
                }  // end per-task body
                // -------- Merge phase --------
                // For each thread, replay its entries (and every variant within) into curr_c
                // via update_derna. Iterate threads in ascending thread-id, and within each
                // thread in insertion order (thread_key_order) — but the tie-break in
                // update_derna is a total order on (score, backtrace_type, bt_info), so the
                // resulting curr_c state is actually order-independent. The ordered iteration
                // just makes the merge reproducible across runs.
                BeamEntry scratch;
                auto replay_entry = [&](int orig_key, const BeamEntry& src) {
                    // Inject each variant as a candidate into curr_c so variant-level
                    // dedup/tiebreak runs against whatever curr_c already holds.
                    for (const auto& v : src.variants) {
                        scratch = src;            // copy static fields (a,b,i,j,...)
                        scratch.x = v.x; scratch.y = v.y;
                        scratch.mfe = v.mfe; scratch.cai = v.cai; scratch.score = v.score;
                        scratch.backtrace_type = v.backtrace_type;
                        scratch.last_closed_nuc = v.last_closed_nuc;
                        scratch.bt_info = v.bt_info;
                        scratch.cs_inner_key    = v.cs_inner_key;
                        scratch.cs_right_s_key  = v.cs_right_s_key;
                        scratch.cs_right_len    = v.cs_right_len;
                        scratch.cs_single_start = v.cs_single_start;
                        scratch.cs_inner_left   = v.cs_inner_left;
                        scratch.cs_inner_right  = v.cs_inner_right;
                        scratch.cs_pack_outer   = v.cs_pack_outer;
                        scratch.variants.clear();  // update_derna reads only top-level.
                        update_derna(curr_c, orig_key, v.score, scratch);
                    }
                };
                for (int t = 0; t < num_threads_blk1; ++t) {
                    if (blk1_use_dense) {
                        DenseBeamTable& td = blk1_thread_dense_maps[t];
                        const auto& order = blk1_thread_key_order[t];
                        for (int compact_idx : order) {
                            const BeamEntry& src = td.slot(compact_idx);
                            if (src.score >= inf) continue;  // shouldn't happen but safe
                            // Reconstruct original key:
                            // compact_idx = left_pos*16 + (nuc_lo*4 + nuc_ro)
                            // orig_key   = (left_pos*nuc_len + pos)*16 + (nuc_lo*4 + nuc_ro)
                            const int left_pos = compact_idx >> 4;
                            const int low      = compact_idx & 0xF;
                            const int orig_key = (left_pos * nuc_len + pos) * 16 + low;
                            replay_entry(orig_key, src);
                        }
                    } else {
                        DernaBeamMap& tm = thread_maps[t];
                        const auto& order = thread_key_order_legacy[t];
                        for (int key : order) {
                            auto it = tm.find(key);
                            if (it == tm.end()) continue;
                            replay_entry(key, it->second);
                        }
                    }
                }
                // Fold per-thread counters into the global diagnostic counts.
                for (int t = 0; t < num_threads_blk1; ++t) {
                    cnt_blk1_rightbulge += cnt_rb_local[t];
                    cnt_blk1_internal  += cnt_in_local[t];
                }
                // VARDIST diagnostic: Block 1's per-thread dense task_maps. Each
                // populated dense slot is by-construction a single (left_pos,
                // nuc_lo, nuc_ro) bucket, so avg_buckets_used should be 1.0
                // and avg_variants_per_bucket should equal avg_variants.
                if (vardist_diag_enabled() && blk1_use_dense) {
                    for (int t = 0; t < num_threads_blk1; ++t) {
                        char tag[32];
                        snprintf(tag, sizeof(tag), "Blk1_t%d", t);
                        vardist_measure_dense(pos, tag,
                            blk1_thread_dense_maps[t],
                            blk1_thread_key_order[t],
                            protein);
                    }
                }
                // BLK1_DIAG=1: per-pos task_map density log. Default off (lossless).
                {
                    static const bool blk1_diag = []() {
                        const char* e = std::getenv("BLK1_DIAG");
                        return e && e[0] == '1';
                    }();
                    if (blk1_diag) {
                        for (int t = 0; t < num_threads_blk1; ++t) {
                            int min_key = INT_MAX, max_key = INT_MIN;
                            std::set<int> distinct_lefts;
                            size_t sz = 0;
                            if (blk1_use_dense) {
                                const auto& order = blk1_thread_key_order[t];
                                if (order.empty()) continue;
                                sz = order.size();
                                for (int compact_idx : order) {
                                    const int left_pos = compact_idx >> 4;
                                    const int low      = compact_idx & 0xF;
                                    const int orig_k = (left_pos * nuc_len + pos) * 16 + low;
                                    if (orig_k < min_key) min_key = orig_k;
                                    if (orig_k > max_key) max_key = orig_k;
                                    distinct_lefts.insert(left_pos);
                                }
                            } else {
                                const DernaBeamMap& tm = thread_maps[t];
                                if (tm.empty()) continue;
                                sz = tm.size();
                                for (const auto& kv : tm) {
                                    int k = kv.first;
                                    if (k < min_key) min_key = k;
                                    if (k > max_key) max_key = k;
                                    distinct_lefts.insert(nuc_key_c_to_left(k, nuc_len));
                                }
                            }
                            fprintf(stderr,
                                "BLK1_DIAG pos=%d tid=%d size=%zu min_key=%d max_key=%d "
                                "distinct_left=%zu nuc_len=%d span=%d\n",
                                pos, t, sz, min_key, max_key,
                                distinct_lefts.size(), nuc_len, max_key - min_key + 1);
                        }
                    }
                }
            }  // end if (!blk1_tasks.empty())
        }  // end Block 1
#endif  // USE_CODON_LOOKUP_CLOSURE

        // ---------- Special-hairpin seeding (LCDSfold initialize_Special_HP_LD port) ----------
        // For each sp_loops key whose length L ends at `pos`, enumerate codon choices realizing
        // the target substring at positions [first_idx .. pos]. Seed curr_c with MANNER_SpecialHP.
        // This bypasses the N-chain prune: even if the N-extension path to (first_idx+1..pos-1)
        // has been pruned, the special-hairpin closure is still seeded into C directly.
        {
            for (int L : SP_HP_LENS) {
                int first_idx = pos - L + 1;
                if (first_idx < 0) continue;
                // Check all sp_loops entries of this length.
                for (const auto& kv : hairpinE) {
                    const std::string& hp = kv.first;
                    if ((int)hp.size() != L) continue;
                    // Filter to only special hairpin keys (those matching LCDSfold sp_loops).
                    // hairpinE contains Triloop/Tetraloop/Hexaloop entries (len 5/6/8 inclusive of pair).
                    int aa_lo = first_idx / 3;
                    int aa_hi = pos / 3;
                    int i_left = first_idx % 3;
                    int j_right = pos % 3;
                    int pa = protein[aa_lo];
                    int pb = protein[aa_hi];
                    // hairpin_loop distance check: pos - first_idx + 1 == L, which is L loop-nts incl closing.
                    // Loop length (unpaired) = L - 2. HAIRPIN_GAP requires (pos - first_idx) > 3 → L > 4 always.
                    if (L <= HAIRPIN_GAP + 1) continue;
                    enum_sp_hp_codons(aa_lo, aa_hi, first_idx, pos, hp, protein,
                        [&](const std::vector<int>& choices, double cai_int) {
                            int x_lo = choices.front();
                            int y_hi = choices.back();
                            // Pair nucleotides from boundary codons.
                            int nuc_lo = nucleotides[pa][x_lo][i_left];
                            int nuc_ro = nucleotides[pb][y_hi][j_right];
                            if (BP_pair[nuc_lo + 1][nuc_ro + 1] == 0) return;
                            // MFE from hairpinE (already-rescaled closed-hairpin special energy).
                            int mfe_raw = hairpinE.at(hp);
                            // Boundary-codon CAI via add_hairpin_CAI_2 (matches N->C convention).
                            double pair_cai = z.add_hairpin_CAI_2(aa_lo, aa_hi, x_lo, y_hi,
                                                                  -1, -1, -1, -1, i_left, j_right);
                            double cai_total = pair_cai + cai_int;
                            double mfe_d = (double)mfe_raw;
                            double sc = combined_score(lambda, mfe_d, cai_total);

                            int key_cl = nuc_key_c_nuc(first_idx, pos, nuc_lo, nuc_ro, nuc_len);
                            BeamEntry ent(sc, aa_lo, aa_hi, i_left, j_right, x_lo, y_hi,
                                          (int)DernaManner::MANNER_SpecialHP, mfe_d, cai_total);
                            // Store all codon choices for aa's in [aa_lo..aa_hi] so traceback can
                            // reconstruct codon_selection without relying on N-chain entries.
                            ent.bt_info = choices;
                            update_derna(curr_c, key_cl, sc, ent);
                        });
                }
            }
        }

        // First C prune: LCDSfold prunes C after CStoC/S_CStoC, before C_StoCS/CtoC.
        prune_cumulative(curr_c, beamsize, pos, n, DernaTableKind::C, best_f_prefix, protein);
        ms_CS_to_C += duration<double, milli>(high_resolution_clock::now() - t_c0).count();

        // ---------- SC_left -> C (close left-bulge: outer pair (outer_left, pos), inner C ends at pos-1) ----------
        // SCLeft[pos_prev] = S_left + C_inner composites where C_inner ends at pos_prev.
        // The outer pair right end is pos; left end is S_left_start - 1.
        // Loop shape: left_unpaired = S_left_len, right_unpaired = 0.
        // NOTE: this path is NOT redundant with MANNER_S_CtoC at line 2815+ — that path
        // only uses the BEST S variant when codons differ (line 2854-2865), while SCLeft
        // enumerates all S variants. Stubbing this drops P01707 to gap=49.9.
        if (!prev_sc_left_flat.empty()) {
            for (int yy = 0; yy < ncod_bb; ++yy) {
                for (const auto& kv : prev_sc_left_flat) {
                    const BeamEntry& scl_ent = *kv.second;

                    // SCLeft must be immediately to the left of pos (same_or_next_codon).
                    // The seam check uses the entry's geometry (a, b, i, j), which is
                    // shared across all variants in the merged entry.
                    bool same_codon_scl = false;
                    if (!same_or_next_codon(scl_ent.b, scl_ent.j, bb, slot, same_codon_scl)) continue;

                    // Recover structural geometry from cs_* fields (shared across variants).
                    int s_left_start   = scl_ent.cs_single_start;
                    int c_inner_left   = scl_ent.cs_inner_left;
                    int seg_len_left   = scl_ent.cs_right_len;
                    if (s_left_start < 0 || c_inner_left < 0 || seg_len_left <= 0) continue;

                    int outer_left_nuc = s_left_start - 1;
                    if (outer_left_nuc < 0) continue;
                    if ((pos - outer_left_nuc) <= (HAIRPIN_GAP + 2)) continue;

                    int aa_o  = outer_left_nuc / 3, ii_o  = outer_left_nuc % 3;
                    int paa_o = protein[aa_o];
                    const int ncod_o = n_codon[paa_o];
                    const int seam_pos_sl = c_inner_left - 1;

                    // Right outer nuc (used when !same_codon_scl)
                    int nuc_ro = nucleotides[pbb][yy][slot];

                    // Iterate variants of the merged SCLeft entry. With dual keying, each
                    // variant carries its own (svar.x = entry.x, cvar.y = entry.y) and
                    // bt_info = {c_key, s_key, cvar.x, svar.y} (the internal codons).
                    for (const auto& scl_var : scl_ent.variants) {
                        // Same-codon seam at the outer-right side: the variant's y (=cvar.y)
                        // must match yy at pos+1 when they fall in one codon.
                        if (same_codon_scl && yy != scl_var.y) continue;
                        const auto& var_bt = scl_var.bt_info;
                        if (var_bt.size() < 4) continue;
                        const int c_key_v   = var_bt[0];
                        const int s_key_v   = var_bt[1];
                        const int inner_x_v = var_bt[2];      // cvar.x
                        const int sleft_y_v = var_bt[3];      // svar.y
                        const int sleft_x_v = scl_var.x;      // svar.x = variant's x
                        const int inner_y_v = scl_var.y;      // cvar.y = variant's y

                        // Look up the inner C and S_left entries (still needed for protein/codon
                        // indices used to compute outer-pair nucleotides).
                        if (c_key_v < 0 || tab_c[pos_prev].count(c_key_v) == 0) continue;
                        const BeamEntry& inner_ent = tab_c[pos_prev].at(c_key_v);
                        if (s_key_v < 0 || seam_pos_sl < 0 || seam_pos_sl >= nuc_len) continue;
                        if (seg_len_left >= (int)tab_s[seam_pos_sl].size()) continue;
                        if (tab_s[seam_pos_sl][seg_len_left].count(s_key_v) == 0) continue;
                        const BeamEntry& s_left_ent = tab_s[seam_pos_sl][seg_len_left].at(s_key_v);

                        const int yy_eff = (same_codon_scl && inner_ent.b == bb) ? inner_y_v : yy;
                        const int nuc_ro_eff = same_codon_scl ? nucleotides[pbb][yy_eff][slot] : nuc_ro;

                        // Inner-pair and S_left mismatch nucs are determined by the variant.
                        int nuc_li  = nucleotides[protein[inner_ent.a]][inner_x_v][inner_ent.i];
                        int nuc_ri  = nucleotides[protein[inner_ent.b]][inner_y_v][inner_ent.j];
                        int nuc_i1  = nucleotides[protein[s_left_ent.a]][sleft_x_v][s_left_ent.i];
                        int nuc_pm1 = nucleotides[protein[s_left_ent.b]][sleft_y_v][s_left_ent.j];
                        int nuc_jm1 = nuc_ri;

                        for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                            if (aa_o == s_left_ent.a && xx_o != sleft_x_v) continue;
                            int nuc_lo = nucleotides[paa_o][xx_o][ii_o];
                            if (BP_pair[nuc_lo + 1][nuc_ro_eff + 1] == 0) continue;

                            double loop_e = score_single_loop_with_mismatch(lambda,
                                                                            nuc_lo, nuc_ro_eff,
                                                                            nuc_li, nuc_ri,
                                                                            seg_len_left, 0,
                                                                            nuc_i1, nuc_jm1,
                                                                            nuc_pm1, -1);
                            if (loop_e >= inf) continue;

                            double mfe_cand = scl_var.mfe + loop_e;
                            double cai_cand = scl_var.cai
                                + (is_last_nuc(pos) ? codon_cai[pbb][yy_eff] : 0.0)
                                + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                            double sc_cand = combined_score(lambda, mfe_cand, cai_cand);

                            int key_cl = nuc_key_c_nuc(outer_left_nuc, pos, nuc_lo, nuc_ro_eff, nuc_len);
                            BeamEntry ent_cl(sc_cand, aa_o, bb, ii_o, slot, xx_o, yy_eff,
                                             (int)DernaManner::MANNER_SCLefttoC, mfe_cand, cai_cand);
                            // Encode the SCLeft entry key + variant (x, y) used, so traceback
                            // can recover the correct merged-table variant.
                            ent_cl.cs_pack_outer = (long long)kv.first;
                            ent_cl.bt_info = {sleft_x_v, inner_y_v};
                            check_c_entry_invariant(lambda, ent_cl, "SCLeft->C", pos);
                            update_derna(curr_c, key_cl, sc_cand, ent_cl);
                            cnt_SCLeft_to_C++;
                        }
                    }
                }
            }
        }
        }

        auto t_c1 = high_resolution_clock::now();
        // ---------- C + S -> CS ----------
        for (int seg_len = 1; seg_len <= min(SINGLE_MAX_LEN, pos - 4); ++seg_len) {
            if (seg_len >= (int)curr_s_flat.size()) continue;
            const auto& s_vec = curr_s_flat[seg_len];
            if (s_vec.empty()) continue;
            int seg_start = pos - seg_len + 1;
            if (seg_start <= 0) continue;
            int closed_right = seg_start - 1;
            if (closed_right < 0 || closed_right >= nuc_len) continue;

            const auto& c_vec_all = get_c_flat(closed_right);
            if (c_vec_all.empty()) continue;

            for (const auto& skv : s_vec) {
                const int s_key = skv.first;
                const BeamEntry& s_ent = *skv.second;
                if (sigma(s_ent.a, s_ent.i) != seg_start || s_ent.b != bb || s_ent.j != slot) continue;

                // With nuc_key_c, the C table merges codon variants by nucleotides.
                // Use the full C vector and check nucleotide compatibility in the inner loop
                // instead of the codon-index bucket optimization (which would miss merged variants).
                const vector<pair<int, const BeamEntry*>>* c_scan = &c_vec_all;

                for (const auto& ckv : *c_scan) {
                    const int c_key = ckv.first;
                    const BeamEntry& c_ent = *ckv.second;
                    if (sigma(c_ent.b, c_ent.j) != closed_right) continue;

                    // Geometry checks first (codon-agnostic).
                    if ((closed_right % 3) != 2) {
                        if (s_ent.a != c_ent.b) continue;
                        if (s_ent.i != c_ent.j + 1) continue;
                    } else {
                        if (s_ent.a != c_ent.b + 1) continue;
                        if (s_ent.i != 0) continue;
                    }

                    // Iterate c_ent variants for codon compatibility with s_ent (and s_ent variants too).
                    for (const auto& cvar : c_ent.variants) {
                        if ((closed_right % 3) != 2) {
                            // FIX (point 3): iterate ALL s_ent variants with sv.x == cvar.y,
                            // not just the first. Multiple svars matching differ in sv.y
                            // (which affects nuc_jm1 = right_S's last nuc → different
                            // downstream closure energies). The previous `break;` after
                            // first match was a coverage gap that forced Block 1 to
                            // brute-force enumerate (S_left × C × S_right) directly.
                            for (const auto& sv : s_ent.variants) {
                                if (sv.x != cvar.y) continue;
                                cnt_C_S_to_CS++;

                                double c_mfe_v = cvar.mfe, c_cai_v = cvar.cai;
                                double s_mfe_v = sv.mfe, s_cai_v = sv.cai;
                                double mfe_cs = c_mfe_v + s_mfe_v;
                                double cai_cs = c_cai_v + s_cai_v;
                                double sc_cs = combined_score(lambda, mfe_cs, cai_cs);

                                int cs_il = sigma(c_ent.a, c_ent.i);
                                int nuc_inner_L = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                                int nuc_closed_R = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];
                                int nuc_single_start = nucleotides[protein[s_ent.a]][sv.x][s_ent.i];
                                int nuc_jm1 = nucleotides[protein[s_ent.b]][sv.y][s_ent.j];
                                int key_cs = cs_get_index(cs_il, nuc_inner_L, nuc_closed_R, nuc_single_start, nuc_jm1, seg_len);
                                BeamEntry ent_cs(sc_cs, s_ent.a, s_ent.b, s_ent.i, s_ent.j, sv.x, sv.y,
                                                 (int)DernaManner::MANNER_C_StoCS, mfe_cs, cai_cs);
                                ent_cs.cs_pack_outer = (long long)key_cs;
                                ent_cs.cs_inner_key = c_key;
                                ent_cs.cs_right_s_key = s_key;
                                ent_cs.cs_right_len = seg_len;
                                ent_cs.cs_single_start = seg_start;
                                ent_cs.cs_inner_left = cs_il;
                                ent_cs.cs_inner_right = sigma(c_ent.b, c_ent.j);
                                ent_cs.bt_info = {c_key, s_key, cvar.x, cvar.y, sv.x, sv.y};
                                update_derna(curr_cs, key_cs, sc_cs, ent_cs);
                            }
                        } else {
                            // Seam at codon boundary: c and s are in different amino acids; iterate all s_ent variants.
                            for (const auto& svar : s_ent.variants) {
                                cnt_C_S_to_CS++;
                                double mfe_cs = cvar.mfe + svar.mfe;
                                double cai_cs = cvar.cai + svar.cai;
                                double sc_cs = combined_score(lambda, mfe_cs, cai_cs);
                                int cs_il = sigma(c_ent.a, c_ent.i);
                                int nuc_inner_L = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                                int nuc_closed_R = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];
                                int nuc_single_start = nucleotides[protein[s_ent.a]][svar.x][s_ent.i];
                                int nuc_jm1 = nucleotides[protein[s_ent.b]][svar.y][s_ent.j];
                                int key_cs = cs_get_index(cs_il, nuc_inner_L, nuc_closed_R, nuc_single_start, nuc_jm1, seg_len);
                                BeamEntry ent_cs(sc_cs, s_ent.a, s_ent.b, s_ent.i, s_ent.j, svar.x, svar.y,
                                                 (int)DernaManner::MANNER_C_StoCS, mfe_cs, cai_cs);
                                ent_cs.cs_pack_outer = (long long)key_cs;
                                ent_cs.cs_inner_key = c_key;
                                ent_cs.cs_right_s_key = s_key;
                                ent_cs.cs_right_len = seg_len;
                                ent_cs.cs_single_start = seg_start;
                                ent_cs.cs_inner_left = cs_il;
                                ent_cs.cs_inner_right = sigma(c_ent.b, c_ent.j);
                                ent_cs.bt_info = {c_key, s_key, cvar.x, cvar.y, svar.x, svar.y};
                                update_derna(curr_cs, key_cs, sc_cs, ent_cs);
                            }
                        }
                    }
                }
            }
        }
        // CS prune: hybrid (top-K raw ∪ top-K cumulative [∪ top-K lookahead if
        // DERNA_LOOKAHEAD_PRUNE=1]). Pure cumulative can drop CS entries whose
        // raw score is locally great but whose cs_il is far from the F-prefix's
        // most-folded prefix — the hybrid keeps both the globally promising
        // (cumulative) and the locally great (raw) entries. The optional
        // lookahead metric tightens by adding min_sl_at[cs_il-1].
        //
        // Selective Block 1 (Approach c): log dropped CS tuples so Block 1 at
        // pos+1 only runs for tuples whose CS coverage was lost.
        static const bool lookahead_prune = []() {
            const char* e = std::getenv("DERNA_LOOKAHEAD_PRUNE");
            return e && e[0] == '1';
        }();
        // Pass dropped-units pointer only when selective Block 1 is active, to
        // keep prune_cs_hybrid's behavior bit-identical to the baseline path
        // when no selective filtering is requested.
        static const bool sel_blk1_active_test = []() {
            const char* e = std::getenv("DERNA_SELECTIVE_BLK1");
            return e && e[0] == '1';
        }();
        std::unordered_set<uint64_t>* drop_arg = sel_blk1_active_test
            ? &tab_cs_dropped_units[pos] : nullptr;
        // Approach (b): two-tier CS storage. DERNA_CS_PRUNE_MULT=N caps CS at
        // N*beamsize entries (default behavior at MULT=2 is hybrid prune).
        // Setting MULT=0 or unset keeps the original 2K cap; MULT=20 ≈ unbounded.
        static const int cs_prune_mult = []() {
            const char* e = std::getenv("DERNA_CS_PRUNE_MULT");
            if (!e) return 0;
            int v = atoi(e);
            return v > 0 ? v : 0;
        }();
        // DERNA_CS_NO_PRUNE=1 disables CS prune entirely (= MULT=infinity).
        static const bool cs_no_prune = []() {
            const char* e = std::getenv("DERNA_CS_NO_PRUNE");
            return e && e[0] == '1';
        }();
        if (cs_no_prune) {
            // No prune — drop only invalid entries.
            for (auto it = curr_cs.begin(); it != curr_cs.end();) {
                if (!std::isfinite(it->second.score) ||
                    !entry_is_position_consistent(it->second, pos, n,
                                                  DernaTableKind::CS, -1))
                    it = curr_cs.erase(it);
                else ++it;
            }
        } else if (cs_prune_mult > 0) {
            // Two-tier prune: keep top-(MULT*K) raw, drop rest. Coarser than
            // hybrid but cheaper than unlimited.
            const int eff_k = beamsize * cs_prune_mult;
            for (auto it = curr_cs.begin(); it != curr_cs.end();) {
                if (!std::isfinite(it->second.score) ||
                    !entry_is_position_consistent(it->second, pos, n,
                                                  DernaTableKind::CS, -1))
                    it = curr_cs.erase(it);
                else ++it;
            }
            if (beamsize > 0 && (int)curr_cs.size() > eff_k) {
                std::vector<std::pair<double, int>> vals;
                vals.reserve(curr_cs.size());
                for (const auto& kv : curr_cs)
                    vals.push_back({kv.second.score, kv.first});
                auto nth = vals.begin() + (eff_k - 1);
                std::nth_element(vals.begin(), nth, vals.end(),
                                 [](const auto& a, const auto& b) {
                                     if (a.first != b.first) return a.first < b.first;
                                     return a.second < b.second;
                                 });
                double thr = nth->first;
                int thr_key = nth->second;
                for (auto it = curr_cs.begin(); it != curr_cs.end();) {
                    double sv = it->second.score;
                    if (sv > thr || (sv == thr && it->first > thr_key)) {
                        if (drop_arg) {
                            const BeamEntry& e = it->second;
                            if (e.cs_inner_key >= 0 && e.cs_right_s_key >= 0 &&
                                e.cs_right_len > 0 && e.cs_inner_right >= 0) {
                                drop_arg->insert(pack_cs_blk1_unit(
                                    e.cs_inner_right, e.cs_right_len,
                                    e.cs_right_s_key, e.cs_inner_key));
                            }
                        }
                        it = curr_cs.erase(it);
                    } else ++it;
                }
            }
        } else if (lookahead_prune) {
            prune_cs_hybrid_lookahead(curr_cs, beamsize, pos, n, best_f_prefix,
                                      min_sl_at_pos, protein, drop_arg);
        } else {
            prune_cs_hybrid(curr_cs, beamsize, pos, n, best_f_prefix, protein, drop_arg);
        }
        tab_cs_selective_active[pos] = (beamsize > 0);
        ms_C_S_to_CS += duration<double, milli>(high_resolution_clock::now() - t_c1).count();

        // [DIAGNOSTIC] CS_VAR_PRED_DIAG: compare per-variant bt_info.data[1] (the
        // s_right_key recorded at C+S→CS construction) against the top-level
        // cs_right_s_key (the winner's s_right_key). Tests whether a per-variant
        // s_right pointer would expose distinct s_right_ents inaccessible to the
        // top-level pointer alone. Default behavior unchanged when env unset.
        if (std::getenv("CS_VAR_PRED_DIAG")) {
            static uint64_t total_variants = 0;
            static uint64_t mismatched_variants = 0;
            static uint64_t total_entries = 0;
            static uint64_t cs_entries_with_any_mismatch = 0;
            uint64_t pos_entries = 0;
            uint64_t pos_total_variants = 0;
            uint64_t pos_mismatched = 0;
            for (const auto& kv : curr_cs) {
                const BeamEntry& ent = kv.second;
                ++total_entries;
                ++pos_entries;
                bool any_mismatch = false;
                for (const auto& v : ent.variants) {
                    ++total_variants;
                    ++pos_total_variants;
                    int variant_s_key = (v.bt_info.len > 1) ? v.bt_info.data[1] : -2;
                    if (variant_s_key != ent.cs_right_s_key) {
                        ++mismatched_variants;
                        ++pos_mismatched;
                        any_mismatch = true;
                    }
                }
                if (any_mismatch) ++cs_entries_with_any_mismatch;
            }
            std::cerr << "[CS_VAR_PRED_DIAG] pos=" << pos
                      << " pos_entries=" << pos_entries
                      << " pos_variants=" << pos_total_variants
                      << " pos_mismatched=" << pos_mismatched
                      << " cum_entries=" << total_entries
                      << " cum_variants=" << total_variants
                      << " cum_mismatched=" << mismatched_variants
                      << " cum_mismatch_pct=" << (total_variants ? 100.0 * mismatched_variants / total_variants : 0.0)
                      << " cum_entries_with_any_mismatch=" << cs_entries_with_any_mismatch << "\n";
        }

        auto t_c2 = high_resolution_clock::now();
        // ---------- C -> C (stacking) ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            DernaBeamMap& prev_c = tab_c[pos_prev];
            for (auto& kv : prev_c) {
                BeamEntry& from_ent = kv.second;
                int from_a = from_ent.a, from_b = from_ent.b, from_i = from_ent.i, from_j = from_ent.j;
                bool same_codon_c = false;
                if (!same_or_next_codon(from_b, from_j, bb, slot, same_codon_c)) continue;

                for (auto& var : from_ent.variants) {
                    if (same_codon_c && yy != var.y) continue;
                    int outer_left_nuc = sigma(from_a, from_i) - 1;
                    if (outer_left_nuc < 0 || (pos - outer_left_nuc) <= (HAIRPIN_GAP + 2)) continue;
                    int aa_o = outer_left_nuc / 3, ii_o = outer_left_nuc % 3;
                    int paa_o = protein[aa_o];
                    const int ncod_o = n_codon[paa_o];
                    int nuc_li_prev = nucleotides[protein[from_a]][var.x][from_i];
                    // nucleotide at pos_prev (right end of previous C). If we crossed a codon boundary, from_b may be bb-1.
                    int nuc_ri_prev = nucleotides[protein[from_b]][var.y][from_j];
                    int nuc_r = nucleotides[pbb][yy][slot];
                    for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                        cnt_C_to_C++;
                        if (aa_o == from_a && xx_o != var.x) continue;
                        int nuc_lo = nucleotides[paa_o][xx_o][ii_o];
                        if (BP_pair[nuc_lo + 1][nuc_r + 1] == 0) continue;
                        int st = Zuker::stacking(nuc_lo, nuc_r, nuc_li_prev, nuc_ri_prev);
                        double mfe_cand = var.mfe + st;
                        double cai_cand = var.cai
                            + (is_last_nuc(pos) ? codon_cai[pbb][yy] : 0.0)
                            + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                        double sc_cand = combined_score(lambda, mfe_cand, cai_cand);
                        int key_cl = nuc_key_c_nuc(outer_left_nuc, pos, nuc_lo, nuc_r, nuc_len);
                        BeamEntry ent_cl(sc_cand, aa_o, bb, ii_o, slot, xx_o, yy, (int)DernaManner::MANNER_CtoC, mfe_cand, cai_cand);
                        ent_cl.bt_info = {kv.first, var.x, var.y};
                        check_c_entry_invariant(lambda, ent_cl, "C->C", pos);
                        update_derna(curr_c, key_cl, sc_cand, ent_cl);
                    }
                }
            }
        }
        // Second C prune: cumulative (bestF[left-1] + local). LCDSfold prune 1 (line 1389).
        prune_cumulative(curr_c, beamsize, pos, n, DernaTableKind::C, best_f_prefix, protein);
        ms_C_to_C += duration<double, milli>(high_resolution_clock::now() - t_c2).count();

        auto t_c3 = high_resolution_clock::now();
        // ---------- S + C -> C (left bulge, right_unpaired=0) ----------
        // Matches LCDSfold MANNER_S_CtoC geometry: inner C at pos-1 (inner_right=pos-1,
        // directly adjacent to outer_right=pos), left S segment of seg_len ending at inner_left-1.
        // Outer pair: (outer_left, pos); inner pair: (inner_left, pos-1). n1=seg_len, n2=0 (left bulge).
        if (pos >= 3) {  // S+C->C left-bulge (matching LCDSfold S_CtoC)
            for (auto& ckv : tab_c[pos - 1]) {
                const int c_key = ckv.first;
                const BeamEntry& c_ent = ckv.second;
                if (sigma(c_ent.b, c_ent.j) != pos - 1) continue;

                const int inner_left  = sigma(c_ent.a, c_ent.i);
                const int seam_pos    = inner_left - 1;  // S must end here
                if (seam_pos < 1) continue;              // need >=1 S nuc + room for outer_left
                if (seam_pos >= (int)tab_s.size()) continue;

                // Left bulge (n2=0): inner_right=pos-1 is directly adjacent to outer_right=pos.
                // Require same_or_next_codon between (c_ent.b, c_ent.j) and (bb, slot).
                bool same_codon_cr = false;
                if (!same_or_next_codon(c_ent.b, c_ent.j, bb, slot, same_codon_cr)) continue;

                // Iterate over C variants so each (x,y) codon pair is considered individually.
                for (const auto& cvar : c_ent.variants) {
                    const int nuc_li = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                    const int nuc_ri = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];

                for (int seg_len = 1; seg_len <= min(SINGLE_MAX_LEN, seam_pos); ++seg_len) {
                    if (seg_len >= (int)tab_s[seam_pos].size()) continue;
                    const auto& s_map = tab_s[seam_pos][seg_len];
                    if (s_map.empty()) continue;

                    for (auto& skv : s_map) {
                        const int s_key = skv.first;
                        const BeamEntry& s_ent = skv.second;
                        // Validate S right end is at seam_pos
                        if (sigma(s_ent.b, s_ent.j) != seam_pos) continue;

                        // Seam consistency: S right end must be immediately before C left end
                        bool same_codon_sc = false;
                        if (!same_or_next_codon(s_ent.b, s_ent.j, c_ent.a, c_ent.i, same_codon_sc)) continue;
                        // Pick S variant with y matching cvar.x (if same codon); else top-level.
                        const XYVariant* svar = nullptr;
                        if (same_codon_sc) {
                            for (const auto& sv : s_ent.variants) {
                                if (sv.y == cvar.x) { svar = &sv; break; }
                            }
                            if (!svar) continue;
                        }
                        int s_x = svar ? svar->x : s_ent.x;
                        int s_y = svar ? svar->y : s_ent.y;
                        double s_mfe_v = svar ? svar->mfe : s_ent.mfe;
                        double s_cai_v = svar ? svar->cai : s_ent.cai;

                        const int outer_left = sigma(s_ent.a, s_ent.i) - 1;
                        if (outer_left < 0) continue;
                        if ((pos - outer_left) <= (HAIRPIN_GAP + 2)) continue;

                        const int aa_o   = outer_left / 3;
                        const int ii_o   = outer_left % 3;
                        const int paa_o  = protein[aa_o];
                        const int ncod_o = n_codon[paa_o];

                        // nuc_i1: nucleotide at outer_left+1 (first nuc of left S = i+1)
                        const int nuc_i1 = nucleotides[protein[s_ent.a]][s_x][s_ent.i];
                        // nuc_pm1: nucleotide at inner_left-1 = seam_pos (last nuc of left S = p-1)
                        const int nuc_pm1 = nucleotides[protein[s_ent.b]][s_y][s_ent.j];

                        for (int yy = 0; yy < ncod_bb; ++yy) {
                            // Use per-variant c codon, not top-level c_ent.y
                            if (same_codon_cr && yy != cvar.y) continue;

                            const int nuc_ro = nucleotides[pbb][yy][slot];

                            for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                                if (aa_o == s_ent.a && xx_o != s_x) continue;
                                const int nuc_lo = nucleotides[paa_o][xx_o][ii_o];
                                if (BP_pair[nuc_lo + 1][nuc_ro + 1] == 0) continue;
                                cnt_S_C_to_C++;

                                // Left bulge (n2=0): score_single_loop_with_mismatch calls Zuker::bulge_loop.
                                // nuc_i1/nuc_pm1 are passed for completeness (not used by bulge scorer).
                                double loop_e = score_single_loop_with_mismatch(
                                    lambda,
                                    nuc_lo, nuc_ro,
                                    nuc_li, nuc_ri,
                                    seg_len, 0,        // left_unpaired=seg_len, right_unpaired=0 (left bulge)
                                    nuc_i1, -1,
                                    nuc_pm1, -1
                                );
                                if (loop_e >= inf) continue;

                                double mfe_cand = cvar.mfe + s_mfe_v + loop_e;
                                double cai_cand = cvar.cai + s_cai_v
                                                + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0)
                                                + (is_last_nuc(pos) ? codon_cai[pbb][yy] : 0.0);
                                double sc_cand = combined_score(lambda, mfe_cand, cai_cand);

                                int key_cl = nuc_key_c_nuc(outer_left, pos, nuc_lo, nuc_ro, nuc_len);
                                BeamEntry ent_cl(
                                    sc_cand, aa_o, bb, ii_o, slot, xx_o, yy,
                                    (int)DernaManner::MANNER_S_CtoC, mfe_cand, cai_cand
                                );
                                ent_cl.bt_info       = {s_key, c_key, cvar.x, cvar.y, s_x, s_y};
                                ent_cl.cs_inner_left = inner_left;
                                ent_cl.cs_right_len  = seg_len;
                                check_c_entry_invariant(lambda, ent_cl, "S+C->C", pos);
                                update_derna(curr_c, key_cl, sc_cand, ent_cl);
                            }
                        }
                    }
                }
                }  // end c_ent.variants loop
            }
        }
        ms_S_C_to_C += duration<double, milli>(high_resolution_clock::now() - t_c3).count();
        ms_C += duration<double, milli>(high_resolution_clock::now() - t0).count();

        t0 = high_resolution_clock::now();
        // ---------- Multi -> Multi (extend) and Multi -> C (close multi) ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            DernaBeamMap& prev_multi = tab_multi[pos_prev];
            for (auto& kv : prev_multi) {
                BeamEntry& from_ent = kv.second;
                int from_a = from_ent.a, from_i = from_ent.i;
                bool same_codon_m = false;
                if (!same_or_next_codon(from_ent.b, from_ent.j, bb, slot, same_codon_m)) continue;

                for (auto& var : from_ent.variants) {
                    // Strict filter for Multi: Multi uses index_struct keys (no yy in key),
                    // so allowing multiple yy values inflates the beam with redundant entries.
                    // Correctness: all yy values are still reachable via different Multi entries.
                    if (same_codon_m && yy != var.y) continue;
                    int last_closed = (from_ent.last_closed_nuc >= 0) ? from_ent.last_closed_nuc : pos_prev;
                    if (pos - last_closed > SINGLE_MAX_LEN) continue;
                    cnt_Multi_EtoMulti++;
                    double mfe_ext = var.mfe, cai_ext = var.cai + ((is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0);
                    double sc_ext = combined_score(lambda, mfe_ext, cai_ext);
                    int nuc_L_mu = nucleotides[protein[from_a]][var.x][from_i];
                    int nuc_R_mu = nucleotides[pbb][yy][slot];
                    int key_ext = codon_beam_key(sigma(from_a, from_i), nuc_L_mu, nuc_R_mu);

                    BeamEntry ent_ext(sc_ext, from_a, bb, from_i, slot, var.x, yy, (int)DernaManner::MANNER_Multi_EtoMulti, mfe_ext, cai_ext, {kv.first, var.x, var.y}, last_closed);
                    // Propagate chain_start so Multi beam span checks work correctly.
                    ent_ext.cs_inner_left = from_ent.cs_inner_left;
                    update_derna(curr_multi, key_ext, sc_ext, ent_ext);

                    int left_nuc = sigma(from_a, from_i) - 1;
                    if (left_nuc < 0) continue;
                    int aa_loop = left_nuc / 3, ii_loop = left_nuc % 3;
                    int paa_loop = protein[aa_loop];
                    const int ncod_loop = n_codon[paa_loop];
                    int nuc_r = nucleotides[pbb][yy][slot];
                    for (int xx_loop = 0; xx_loop < ncod_loop; ++xx_loop) {
                        if (aa_loop == from_a && xx_loop != var.x) continue;
                        int nuc_lo = nucleotides[paa_loop][xx_loop][ii_loop];
                        if (BP_pair[nuc_lo + 1][nuc_r + 1] == 0) continue;
                        cnt_MultitoC++;
                        // LCDSfold-equivalent multi closing penalty: v_score_multi = E_MLstem(pair_type(reversed)) + ML_closing37.
                        double close_pen = (double)v_score_multi_pb(nuc_lo, nuc_r);
                        double cai_cl = var.cai
                            + (is_last_nuc(pos) ? codon_cai[pbb][yy] : 0.0)
                            + (is_last_nuc(left_nuc) ? codon_cai[paa_loop][xx_loop] : 0.0);

                        double sc_cl = combined_score(lambda, var.mfe + close_pen, cai_cl);

                        int key_cl = nuc_key_c_nuc(left_nuc, pos, nuc_lo, nuc_r, nuc_len);
                        BeamEntry ent_cl(sc_cl, aa_loop, bb, ii_loop, slot, xx_loop, yy, (int)DernaManner::MANNER_MultitoC,
                            var.mfe + close_pen, cai_cl);
                        ent_cl.bt_info = {kv.first, var.x, var.y};
                        check_c_entry_invariant(lambda, ent_cl, "Multi->C", pos);
                        update_derna(curr_c, key_cl, sc_cl, ent_cl);
                    }
                }
            }
        }
        // Third C prune: cumulative (bestF[left-1] + local), LCDSfold prune 2 (line 1434).
        prune_cumulative(curr_c, beamsize, pos, n, DernaTableKind::C, best_f_prefix, protein);
        // VARDIST diagnostic: per-pos C-table variants[] distribution.
        vardist_measure_map(pos, "C", curr_c, protein);

        // ---- Build codon-pair aggregation for C table at this pos (post-prune) ----
        // The Block 1 supplemental pass below is gated off (`if (false && ...)`), so
        // this prune is the FINAL C prune for this position. If that ever changes,
        // move this build call after the fourth prune at line ~4292.
        {
            CodonAggTables::PerPosDiag* c_diag_ptr = agg_diag_on ? &agg_tables.diag[pos] : nullptr;
            build_c_codon_agg(curr_c, agg_tables.c_agg[pos], c_diag_ptr);
        }

        // Per-pos diagnostic dump: write S+C row to output/agg_diagnostics.log.
        if (agg_diag_on) {
            agg_diag_write_pos(pos, agg_tables.diag[pos]);
        }

        // ---------- Block 1 (LCDSfold) supplemental pass: runs AFTER third prune ----------
        // By running here, Block 1 supplements the top-k C entries with combinations missed
        // by the CS composite intermediate pruning, without polluting the CS→C prune pool.
        // After Block 1 adds new entries, a fourth prune limits curr_c back to beamsize.
        if (false && pos >= 5) {  // Block 1 supplemental pass disabled (hurts small k)
            for (int seg_len_right = 1; seg_len_right <= min(SINGLE_MAX_LEN, pos_prev); ++seg_len_right) {
                if (seg_len_right >= (int)tab_s[pos_prev].size()) continue;
                const auto& s_right_map = tab_s[pos_prev][seg_len_right];
                if (s_right_map.empty()) continue;

                int seg_start = pos_prev - seg_len_right + 1;
                int closed_right = seg_start - 1;
                if (closed_right < 0 || closed_right >= nuc_len) continue;

                const auto& c_flat_blk1 = get_c_flat(closed_right);
                if (c_flat_blk1.empty()) continue;

                for (const auto& s_rkv : s_right_map) {
                    const int s_right_key = s_rkv.first;
                    const BeamEntry& s_right_ent = s_rkv.second;
                    if (sigma(s_right_ent.b, s_right_ent.j) != pos_prev) continue;
                    if (sigma(s_right_ent.a, s_right_ent.i) != seg_start) continue;
                    bool same_codon_sr = false;
                    if (!same_or_next_codon(s_right_ent.b, s_right_ent.j, bb, slot, same_codon_sr)) continue;

                    for (const auto& c_kv : c_flat_blk1) {
                        const int c_key = c_kv.first;
                        const BeamEntry& c_ent = *c_kv.second;
                        if (sigma(c_ent.b, c_ent.j) != closed_right) continue;

                        if ((closed_right % 3) != 2) {
                            if (s_right_ent.a != c_ent.b) continue;
                            if (s_right_ent.i != c_ent.j + 1) continue;
                            if (s_right_ent.x != c_ent.y) continue;
                        } else {
                            if (s_right_ent.a != c_ent.b + 1) continue;
                            if (s_right_ent.i != 0) continue;
                        }

                        int inner_left = sigma(c_ent.a, c_ent.i);
                        int outer_left_nuc = inner_left - 1;
                        if (outer_left_nuc < 0) continue;
                        if ((pos - outer_left_nuc) <= (HAIRPIN_GAP + 2)) continue;

                        int aa_o = outer_left_nuc / 3, ii_o = outer_left_nuc % 3;
                        int paa_o = protein[aa_o];
                        const int ncod_o = n_codon[paa_o];
                        int nuc_li = nucleotides[protein[c_ent.a]][c_ent.x][c_ent.i];
                        int nuc_ri = nucleotides[protein[c_ent.b]][c_ent.y][c_ent.j];
                        const int nuc_jm1_blk1 = nucleotides[protein[s_right_ent.b]][s_right_ent.y][s_right_ent.j];
                        const int nuc_qp1_blk1 = nucleotides[protein[s_right_ent.a]][s_right_ent.x][s_right_ent.i];

                        for (int yy = 0; yy < ncod_bb; ++yy) {
                            if (same_codon_sr && yy != s_right_ent.y) continue;
                            const int yy_eff = (same_codon_sr && c_ent.b == bb) ? c_ent.y : yy;
                            const int nuc_ro_eff = same_codon_sr ? nucleotides[pbb][yy_eff][slot] : nucleotides[pbb][yy][slot];

                            for (int xx_o = 0; xx_o < ncod_o; ++xx_o) {
                                if (aa_o == c_ent.a && xx_o != c_ent.x) continue;
                                int nuc_lo = nucleotides[paa_o][xx_o][ii_o];
                                if (BP_pair[nuc_lo + 1][nuc_ro_eff + 1] == 0) continue;
                                cnt_blk1_rightbulge++;

                                double loop_e = score_single_loop_with_mismatch(lambda,
                                    nuc_lo, nuc_ro_eff, nuc_li, nuc_ri,
                                    0, seg_len_right,
                                    nuc_li, nuc_jm1_blk1, nuc_lo, nuc_qp1_blk1);
                                if (loop_e >= inf) continue;

                                double mfe_cand = c_ent.mfe + s_right_ent.mfe + loop_e;
                                double cai_cand = c_ent.cai + s_right_ent.cai
                                    + (is_last_nuc(pos) ? codon_cai[pbb][yy_eff] : 0.0)
                                    + (ii_o == 2 ? codon_cai[paa_o][xx_o] : 0.0);
                                double sc_cand = combined_score(lambda, mfe_cand, cai_cand);

                                int key_cl = nuc_key_c_nuc(outer_left_nuc, pos, nuc_lo, nuc_ro_eff, nuc_len);
                                BeamEntry ent_cl(sc_cand, aa_o, bb, ii_o, slot, xx_o, yy_eff,
                                                 (int)DernaManner::MANNER_C_StoC, mfe_cand, cai_cand);
                                ent_cl.bt_info = {c_key, s_right_key};
                                ent_cl.cs_inner_left = closed_right;
                                ent_cl.cs_right_len = seg_len_right;
                                check_c_entry_invariant(lambda, ent_cl, "Blk1_C_StoC", pos);
                                update_derna(curr_c, key_cl, sc_cand, ent_cl);
                            }
                        }  // end yy loop
                    }  // end c_flat_blk1 loop
                }  // end s_right_map loop
            }  // end seg_len_right loop
            // Fourth C prune: limit back to beamsize after Block 1 supplemental pass.
            prune_beam_derna_checked(curr_c, beamsize, pos, n, DernaTableKind::C);
        }  // end Block 1 supplemental pass

        // ---------- Build SC_left = S_left + C composites (for left-bulge closing at pos+1) ----------
        // Mirror of C+S->CS: for each C in curr_c (pruned), find S_left segments ending at
        // sigma(C.a, C.i) - 1. The SC_left composite is used at pos+1 to close a left-bulge.
        // NOT REDUNDANT with MANNER_S_CtoC at line 2815+: that path uses only the BEST
        // (top-level) S variant in the !same_codon case, while SCLeft enumerates all S
        // variants. Verified empirically: stubbing SCLeft drops P01707 to gap=49.9.
        if (pos + 1 < nuc_len) {
            for (auto& ckv : curr_c) {
                const int c_key = ckv.first;
                const BeamEntry& c_ent = ckv.second;
                if (!std::isfinite(c_ent.mfe) || !std::isfinite(c_ent.cai)) continue;

                const int c_inner_left = sigma(c_ent.a, c_ent.i);
                const int seam_pos_sl  = c_inner_left - 1;
                if (seam_pos_sl < 0 || seam_pos_sl >= nuc_len) continue;

                for (int seg_len = 1; seg_len <= min(SINGLE_MAX_LEN, seam_pos_sl + 1); ++seg_len) {
                    if (seg_len >= (int)tab_s[seam_pos_sl].size()) continue;
                    const auto& s_map = tab_s[seam_pos_sl][seg_len];
                    if (s_map.empty()) continue;

                    for (auto& skv : s_map) {
                        const int s_key = skv.first;
                        const BeamEntry& s_ent = skv.second;
                        if (sigma(s_ent.b, s_ent.j) != seam_pos_sl) continue;
                        // S_left must end immediately before C's left boundary
                        bool same_codon_sl = false;
                        if (!same_or_next_codon(s_ent.b, s_ent.j, c_ent.a, c_ent.i, same_codon_sl)) continue;

                        const int s_left_start = sigma(s_ent.a, s_ent.i);
                        const int outer_left_nuc = s_left_start - 1;
                        if (outer_left_nuc < 0) continue;

                        // Iterate both C and S variants; seam (s.y == c.x) is per-variant.
                        // Dual-key merge: structural+nucleotide key collapses (cvar, svar) tuples
                        // that produce identical closure-relevant nucleotides; per (key, x=svar.x,
                        // y=cvar.y) the best score wins. bt_info holds the *internal* codons
                        // (cvar.x, svar.y) needed for traceback recovery — the boundary codons
                        // (svar.x, cvar.y) are already in (entry.x, entry.y).
                        for (const auto& cvar : c_ent.variants) {
                            const int nuc_li_v = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                            const int nuc_ri_v = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];
                            for (const auto& svar : s_ent.variants) {
                                if (same_codon_sl && svar.y != cvar.x) continue;
                                const int nuc_i1_v  = nucleotides[protein[s_ent.a]][svar.x][s_ent.i];
                                const int nuc_pm1_v = nucleotides[protein[s_ent.b]][svar.y][s_ent.j];

                                double mfe_scl = cvar.mfe + svar.mfe;
                                double cai_scl = cvar.cai + svar.cai;
                                double sc_scl  = combined_score(lambda, mfe_scl, cai_scl);

                                int key_scl = sc_left_key_encode(s_left_start, seg_len,
                                                                 nuc_li_v, nuc_ri_v, nuc_i1_v, nuc_pm1_v);

                                BeamEntry ent_scl(sc_scl, s_ent.a, c_ent.b, s_ent.i, c_ent.j,
                                                  svar.x, cvar.y, (int)DernaManner::MANNER_C_StoSCLeft,
                                                  mfe_scl, cai_scl);
                                ent_scl.cs_inner_key    = c_key;
                                ent_scl.cs_right_s_key  = s_key;
                                ent_scl.cs_right_len    = seg_len;
                                ent_scl.cs_single_start = s_left_start;
                                ent_scl.cs_inner_left   = c_inner_left;
                                ent_scl.bt_info = {c_key, s_key, cvar.x, svar.y};
                                update_derna(curr_sc_left, key_scl, sc_scl, ent_scl);
                                cnt_SC_left_build++;
                            }
                        }
                    }
                }
            }
            // With dual-keyed SCLeft, entries merge by structural+nucleotide signature
            // and per (key, x, y) the best (mfe, cai) is kept. |SCLeft| is bounded by
            // the number of distinct (s_left_start, seg_len_left, 4 nucs) tuples present
            // at this pos, well under beamsize in practice. We still apply a cumulative
            // prune for safety on extreme cases.
            prune_sc_left_hybrid(curr_sc_left, beamsize, best_f_prefix);
        }

        // ---------- C -> M1 and M1 + C -> M2 ----------
        for (auto& ckv : curr_c) {
            int key_c = ckv.first;
            BeamEntry& c_ent = ckv.second;
            if (!std::isfinite(c_ent.cai) || !std::isfinite(c_ent.mfe) || cai_looks_garbage(c_ent.cai)) {
                continue;  // do not propagate bad C entry
            }
            int seam_pos = sigma(c_ent.a, c_ent.i) - 1;
            int req_aa = (c_ent.i == 0) ? (c_ent.a - 1) : c_ent.a;
            int req_ii = (c_ent.i == 0) ? 2 : (c_ent.i - 1);
            // Iterate C variants: each (cvar.x, cvar.y) gives different M1 keys/penalties.
            for (const auto& cvar : c_ent.variants) {
                if (!std::isfinite(cvar.cai) || !std::isfinite(cvar.mfe) || cai_looks_garbage(cvar.cai)) continue;
                int nuc_left  = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                int nuc_right = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];
                double pen_ml = (double)v_score_M1_pb(nuc_left, nuc_right);
                cnt_C_to_M1++;
                double sc_ml1 = combined_score(lambda, cvar.mfe + pen_ml, cvar.cai);
                BeamEntry ent_ml1(sc_ml1, c_ent.a, c_ent.b, c_ent.i, c_ent.j, cvar.x, cvar.y, (int)DernaManner::MANNER_CtoM1,
                    cvar.mfe + pen_ml, cvar.cai, {key_c, cvar.x, cvar.y}, pos);
                ent_ml1.cs_inner_left = sigma(c_ent.a, c_ent.i);
                check_m1_entry_invariant(lambda, ent_ml1, "C->M1");
                update_derna(curr_m1, codon_beam_key(sigma(c_ent.a, c_ent.i), nuc_left, nuc_right), sc_ml1, ent_ml1);

                if (seam_pos >= HAIRPIN_GAP) {
                    for (auto& m1kv : tab_m1[seam_pos]) {
                        BeamEntry& m1_ent = m1kv.second;
                        // Filter removed (was: m1_ent.bt_info.size() < 2 continue) to allow
                        // CtoM1-derived M1 entries into M1+C→M2. This is required for the
                        // LCDS-optimal multi-loop assembly on O14880/Q14442/Q15669/P01707 etc.
                        // bt_info enrichment + per-variant manner dispatch (below) handle traceback.
                        if (!std::isfinite(m1_ent.cai) || !std::isfinite(m1_ent.mfe) || cai_looks_garbage(m1_ent.cai)) continue;
                        if (sigma(m1_ent.b, m1_ent.j) != seam_pos) continue;
                        if (m1_ent.b != req_aa || m1_ent.j != req_ii) continue;
                        // Iterate M1 variants; seam compatibility (m1.y == cvar.x) is per-variant.
                        for (const auto& m1var : m1_ent.variants) {
                            if (!std::isfinite(m1var.cai) || !std::isfinite(m1var.mfe) || cai_looks_garbage(m1var.cai)) continue;
                            if (c_ent.i > 0 && (m1_ent.b != c_ent.a || m1var.y != cvar.x)) continue;
                            cnt_M1_C_to_M2++;
                            double pen_ml2 = pen_ml;  // invariant in m1var loop
                            double mfe_ml2 = m1var.mfe + cvar.mfe + pen_ml2;
                            int chain_start_cs = (m1_ent.cs_inner_left >= 0)
                                                  ? m1_ent.cs_inner_left
                                                  : sigma(m1_ent.a, m1_ent.i);
                            double cai_ml2 = m1var.cai + cvar.cai;
                            double sc_ml2 = combined_score(lambda, mfe_ml2, cai_ml2);
                            int nuc_L_m2 = nucleotides[protein[m1_ent.a]][m1var.x][m1_ent.i];
                            int nuc_R_m2 = nucleotides[protein[c_ent.b]][cvar.y][c_ent.j];
                            int key_ml2 = codon_beam_key(sigma(m1_ent.a, m1_ent.i), nuc_L_m2, nuc_R_m2);
                            BeamEntry ent_ml2(sc_ml2, m1_ent.a, c_ent.b, m1_ent.i, c_ent.j, m1var.x, cvar.y, (int)DernaManner::MANNER_M1_CtoM2,
                                mfe_ml2, cai_ml2, {m1kv.first, key_c, m1var.x, m1var.y, cvar.x, cvar.y}, pos);
                            ent_ml2.cs_inner_left = chain_start_cs;
                            update_derna(curr_m2, key_ml2, sc_ml2, ent_ml2);
                        }
                    }
                }
            }
        }
        // M2 prune: cumulative (bestF[left-1] + local), mirrors LCDSfold BeamPrune(false).
        prune_cumulative(curr_m2, beamsize, pos, n, DernaTableKind::M2, best_f_prefix, protein);

        // ---------- M1 -> M1 (extend unpaired in multi) ----------
        // Predecessor: right end at pos_prev. At pos we advance right by one nuc.
        // Same codon: from_ent.b == bb, from_ent.j == slot - 1. Codon boundary: from_ent.b == bb-1, from_ent.j == 2, slot == 0.
        for (int yy = 0; yy < ncod_bb; ++yy) {
            for (auto& kv : tab_m1[pos_prev]) {
                BeamEntry& from_ent = kv.second;
                int last_closed = (from_ent.last_closed_nuc >= 0) ? from_ent.last_closed_nuc : pos_prev;
                if (pos - last_closed > SINGLE_MAX_LEN) continue;
                bool right_adjacent = (slot > 0 && from_ent.b == bb && from_ent.j == slot - 1)
                    || (slot == 0 && from_ent.b == bb - 1 && from_ent.j == 2);
                if (!right_adjacent) continue;
                if (!std::isfinite(from_ent.cai) || !std::isfinite(from_ent.mfe) || cai_looks_garbage(from_ent.cai)) {
                    continue;  // do not propagate bad M1
                }
                // Two-pass cross-level blocking: if the M1's chain_start is to the left of a
                // C stacking chain that covers the current position, block the extension.

                for (auto& var : from_ent.variants) {
                    // Same-codon consistency: if extending within the same codon,
                    // the codon choice must match (same check as N->N, S->S, F->F, Multi->Multi).
                    bool same_codon_m1 = (slot > 0 && from_ent.b == bb && from_ent.j == slot - 1);
                    if (same_codon_m1 && yy != var.y) continue;
                    cnt_M1_EtoM1++;
                    double mfe_ext = var.mfe, cai_ext = var.cai + ((is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0);
                    double sc_ext = combined_score(lambda, mfe_ext, cai_ext);
                    int nuc_L_m1 = nucleotides[protein[from_ent.a]][var.x][from_ent.i];
                    int nuc_R_m1 = nucleotides[pbb][yy][slot];
                    int key_ext = codon_beam_key(sigma(from_ent.a, from_ent.i), nuc_L_m1, nuc_R_m1);
                    BeamEntry ent_ext(sc_ext, from_ent.a, bb, from_ent.i, slot, var.x, yy, (int)DernaManner::MANNER_M1_EtoM1, mfe_ext, cai_ext, {kv.first, var.x, var.y}, last_closed);
                    // Propagate chain_start through M1 extension.
                    ent_ext.cs_inner_left = from_ent.cs_inner_left;
                    check_m1_entry_invariant(lambda, ent_ext, "M1->M1");
                    update_derna(curr_m1, key_ext, sc_ext, ent_ext);
                }
            }
        }

        // ---------- M2 -> M1, M2 -> Multi, S + M2 -> Multi ----------
        for (auto& m2kv : curr_m2) {
            const int m2_key = m2kv.first;
            const BeamEntry& m2_ent = m2kv.second;
            int last_closed = (m2_ent.last_closed_nuc >= 0) ? m2_ent.last_closed_nuc : pos;
        
            if (!std::isfinite(m2_ent.cai) || !std::isfinite(m2_ent.mfe) || cai_looks_garbage(m2_ent.cai)) {
                continue;  // do not propagate bad M2
            }
            // Pass-through: propagate all M2 variants into M1 and Multi.
            for (const auto& mvar : m2_ent.variants) {
                if (!std::isfinite(mvar.cai) || !std::isfinite(mvar.mfe) || cai_looks_garbage(mvar.cai)) continue;
                BeamEntry ent_m2_to_m1(mvar.score, m2_ent.a, m2_ent.b, m2_ent.i, m2_ent.j, mvar.x, mvar.y,
                                       (int)DernaManner::MANNER_M2toM1, mvar.mfe, mvar.cai, {m2_key}, last_closed);
                ent_m2_to_m1.cs_inner_left = m2_ent.cs_inner_left;
                BeamEntry ent_m2_to_multi(mvar.score, m2_ent.a, m2_ent.b, m2_ent.i, m2_ent.j, mvar.x, mvar.y,
                                          (int)DernaManner::MANNER_M2toMulti, mvar.mfe, mvar.cai, {m2_key}, last_closed);
                ent_m2_to_multi.cs_inner_left = m2_ent.cs_inner_left;
                cnt_M2_to_M1++;
                check_m1_entry_invariant(lambda, ent_m2_to_m1, "M2->M1");
                update_derna(curr_m1, m2_key, mvar.score, ent_m2_to_m1);
                update_derna(curr_multi, m2_key, mvar.score, ent_m2_to_multi);
            }

            int seam_pos = sigma(m2_ent.a, m2_ent.i) - 1;
            if (seam_pos >= 0) {
                for (int seg_len = 1; seg_len <= min(SINGLE_MAX_LEN, seam_pos); ++seg_len) {
                    if (seam_pos - seg_len + 1 < 0) continue;
                    int seg_start = seam_pos - seg_len + 1;
                    for (auto& skv : tab_s[seam_pos][seg_len]) {
                        BeamEntry& s_ent = skv.second;
                        if (sigma(s_ent.a, s_ent.i) != seg_start) continue;
                        if (sigma(s_ent.b, s_ent.j) != seam_pos) continue;

                        // Iterate both M2 and S variants; seam compat (s.y == m2.x) is per-variant.
                        for (const auto& m2var : m2_ent.variants) {
                            if (!std::isfinite(m2var.cai) || !std::isfinite(m2var.mfe) || cai_looks_garbage(m2var.cai)) continue;
                            for (auto& s_var : s_ent.variants) {
                                if (m2_ent.i != 0 && s_var.y != m2var.x) continue;
                                double mfe_comb = m2var.mfe + s_var.mfe, cai_comb = m2var.cai + s_var.cai;
                                double sc_comb = combined_score(lambda, mfe_comb, cai_comb);
                                int nuc_L_sm2 = nucleotides[protein[s_ent.a]][s_var.x][s_ent.i];
                                int nuc_R_sm2 = nucleotides[protein[m2_ent.b]][m2var.y][m2_ent.j];
                                int key_comb = codon_beam_key(sigma(s_ent.a, s_ent.i), nuc_L_sm2, nuc_R_sm2);
                                BeamEntry ent_comb(sc_comb, s_ent.a, m2_ent.b, s_ent.i, m2_ent.j, s_var.x, m2var.y, (int)DernaManner::MANNER_S_M2toMulti,
                                    mfe_comb, cai_comb);

                                // bt[0]=s_key, bt[1]=m2_key, bt[2]=s_var.y,
                                // bt[3]=seg_len, bt[4]=seam_pos
                                ent_comb.bt_info = {skv.first, m2_key, s_var.y, seg_len, seam_pos};
                                ent_comb.last_closed_nuc = last_closed;
                                ent_comb.cs_inner_left = m2_ent.cs_inner_left;
                                update_derna(curr_multi, key_comb, sc_comb, ent_comb);
                            }
                        }
                    }
                }
            }
        }
        // M1/Multi prune: cumulative (bestF[left-1] + local), mirrors LCDSfold BeamPrune(false).
        prune_cumulative(curr_m1, beamsize, pos, n, DernaTableKind::M1, best_f_prefix, protein);
        prune_cumulative(curr_multi, beamsize, pos, n, DernaTableKind::Multi, best_f_prefix, protein);
        ms_Multi += duration<double, milli>(high_resolution_clock::now() - t0).count();

        t0 = high_resolution_clock::now();
        // ---------- F: extend F and C -> F, F + C -> F ----------
        for (int yy = 0; yy < ncod_bb; ++yy) {
            for (auto& kv : tab_f[pos_prev]) {
                BeamEntry& from_ent = kv.second;
                bool right_adjacent = (slot > 0 && from_ent.b == bb && from_ent.j == slot - 1)
                    || (slot == 0 && from_ent.b == bb - 1 && from_ent.j == 2);
                if (!right_adjacent) continue;

                for (auto& var : from_ent.variants) {
                    // Fix 11: enforce same-codon consistency (same as N_EtoN and S_EtoS)
                    bool same_codon_f = (slot > 0 && from_ent.b == bb && from_ent.j == slot - 1);
                    if (same_codon_f && yy != var.y) continue;
                    cnt_F_EtoF++;
                    double mfe_ext = var.mfe, cai_ext = var.cai + ((is_last_nuc(pos)) ? codon_cai[pbb][yy] : 0.0);
                    double sc_ext = combined_score(lambda, mfe_ext, cai_ext);
                    int nuc_L_f = nucleotides[protein[from_ent.a]][var.x][from_ent.i];
                    int nuc_R_f = nucleotides[pbb][yy][slot];
                    int key_ext = codon_beam_key(sigma(from_ent.a, from_ent.i), nuc_L_f, nuc_R_f);
                    BeamEntry ent_ext(sc_ext, from_ent.a, bb, from_ent.i, slot, var.x, yy, (int)DernaManner::MANNER_F_EtoF, mfe_ext, cai_ext);
                    ent_ext.bt_info = {kv.first, var.x, var.y};
                    update_derna(curr_f, key_ext, sc_ext, ent_ext);
                }
            }
        }
        for (auto& ckv : curr_c) {
            BeamEntry& c_ent = ckv.second;
            int left_bound = sigma(c_ent.a, c_ent.i);
            int seam_pos = left_bound - 1;
            // Iterate C variants — each gives different nuc_left/nuc_right and F keys.
            for (const auto& cvar : c_ent.variants) {
                if (!std::isfinite(cvar.cai) || !std::isfinite(cvar.mfe) || cai_looks_garbage(cvar.cai)) continue;
                int nuc_left  = nucleotides[protein[c_ent.a]][cvar.x][c_ent.i];
                int nuc_right = nucleotides[pbb][cvar.y][c_ent.j];
                double ext_pen_raw = v_score_external_paired_pb(nuc_left, nuc_right);

                if (left_bound == 0) {
                    cnt_C_to_F++;
                    double mfe_f = cvar.mfe + ext_pen_raw;
                    double sc_f = combined_score(lambda, mfe_f, cvar.cai);
                    int key_f = codon_beam_key(sigma(c_ent.a, c_ent.i), nuc_left, nuc_right);
                    BeamEntry ent_f_c(sc_f, c_ent.a, c_ent.b, c_ent.i, c_ent.j, cvar.x, cvar.y, (int)DernaManner::MANNER_CtoF, mfe_f, cvar.cai);
                    update_derna(curr_f, key_f, sc_f, ent_f_c);
                }
                if (seam_pos >= 0) {
                    for (auto& fkv : tab_f[seam_pos]) {
                        BeamEntry& f_ent = fkv.second;
                        int f_left = sigma(f_ent.a, f_ent.i);
                        if (f_left != 0 && pos != nuc_len - 1) continue;
                        if (sigma(f_ent.b, f_ent.j) != seam_pos) continue;
                        bool same_codon_fc = (f_ent.b == c_ent.a && f_ent.j + 1 == c_ent.i);
                        bool cross_codon_fc = (c_ent.i == 0 && f_ent.b + 1 == c_ent.a && f_ent.j == 2);
                        if (!same_codon_fc && !cross_codon_fc) continue;

                        for (auto& f_var : f_ent.variants) {
                            // Per-variant seam check against this C variant.
                            if (same_codon_fc && f_var.y != cvar.x) continue;
                            cnt_F_C_to_F++;
                            double mfe_f = f_var.mfe + cvar.mfe + ext_pen_raw;
                            double cai_f = f_var.cai + cvar.cai;
                            double sc_f = combined_score(lambda, mfe_f, cai_f);
                            int nuc_L_fc = nucleotides[protein[f_ent.a]][f_var.x][f_ent.i];
                            int nuc_R_fc = nuc_right;
                            int key_f = codon_beam_key(sigma(f_ent.a, f_ent.i), nuc_L_fc, nuc_R_fc);
                            BeamEntry ent_f(sc_f, f_ent.a, c_ent.b, f_ent.i, c_ent.j, f_var.x, cvar.y, (int)DernaManner::MANNER_F_CtoF,
                                mfe_f, cai_f);
                            ent_f.bt_info = {fkv.first, ckv.first, seam_pos, f_var.y, cvar.x, cvar.y};
                            update_derna(curr_f, key_f, sc_f, ent_f);
                        }
                    }
                }
            }
        }
        // LCDSfold does NOT prune F — keep all F entries for accurate best_f_prefix.
        // prune_beam_derna_checked(curr_f, beamsize, pos, n, DernaTableKind::F);
        // Update best_f_prefix per codon for codon-compatible cumulative pruning.
        for (const auto& kv : curr_f) {
            for (const auto& v : kv.second.variants) {
                int nuc_r = nucleotides[protein[kv.second.b]][v.y][kv.second.j];
                best_f_prefix[pos][nuc_r] = min(best_f_prefix[pos][nuc_r], v.score);
            }
        }
        ms_F += duration<double, milli>(high_resolution_clock::now() - t0).count();
    }
    s_prune_protein = nullptr;
    if (s_prune_log.is_open()) {
        s_prune_log.close();
    }
    if (s_dualbeam_csv.is_open()) {
        s_dualbeam_csv.close();
    }
    // VARDIST: aggregate per-thread linear-find counters and emit summary.
    if (vardist_diag_enabled()) {
        VarDistDiag totals;
        #pragma omp parallel
        {
            std::lock_guard<std::mutex> lk(g_vardist_diag_mu);
            totals.fast_path_calls   += g_vardist_diag.fast_path_calls;
            totals.linear_find_calls += g_vardist_diag.linear_find_calls;
            g_vardist_diag = VarDistDiag{};
        }
        const uint64_t total = totals.fast_path_calls + totals.linear_find_calls;
        const double pct_linear = (total > 0)
            ? (100.0 * (double)totals.linear_find_calls / (double)total) : 0.0;
        cerr << "[VARDIST_DIAG] update_derna_existing fast_path=" << totals.fast_path_calls
             << " linear_find=" << totals.linear_find_calls
             << " total=" << total
             << " linear_find_pct=" << pct_linear << "%" << endl;
        if (s_vardist_tsv.is_open()) s_vardist_tsv.close();
    }
    // BLK1_AGG_DIAG: gather per-thread diagnostic counters and emit summary.
    if (blk1_agg_diag_enabled()) {
        // Collect each OpenMP thread's thread-local g_blk1_agg_diag into a shared vector.
        Blk1AggDiag totals;
        #pragma omp parallel
        {
            std::lock_guard<std::mutex> lk(g_blk1_agg_diag_mu);
            totals.seed_calls             += g_blk1_agg_diag.seed_calls;
            totals.new_variant_calls      += g_blk1_agg_diag.new_variant_calls;
            totals.existing_variant_calls += g_blk1_agg_diag.existing_variant_calls;
            totals.total_update_ns        += g_blk1_agg_diag.total_update_ns;
            if (g_blk1_agg_diag.per_pos_max_redundancy > totals.per_pos_max_redundancy)
                totals.per_pos_max_redundancy = g_blk1_agg_diag.per_pos_max_redundancy;
            // Reset thread-local for any subsequent runs.
            g_blk1_agg_diag = Blk1AggDiag{};
        }
        const uint64_t total_calls    = totals.seed_calls + totals.new_variant_calls + totals.existing_variant_calls;
        const uint64_t distinct_slots = totals.seed_calls + totals.new_variant_calls;
        const double avg_redundancy = (distinct_slots > 0)
            ? (double)total_calls / (double)distinct_slots : 0.0;
        const double ns_per_call = (total_calls > 0)
            ? (double)totals.total_update_ns / (double)total_calls : 0.0;
        // Wall-clock total update time across all threads.
        cerr << "[BLK1_AGG_DIAG] total_calls=" << total_calls
             << " seed=" << totals.seed_calls
             << " new_variant=" << totals.new_variant_calls
             << " existing_variant=" << totals.existing_variant_calls
             << " distinct_slots=" << distinct_slots
             << " avg_redundancy=" << avg_redundancy
             << " max_redundancy=" << totals.per_pos_max_redundancy
             << " sum_thread_ns=" << totals.total_update_ns
             << " ns_per_call=" << ns_per_call
             << endl;
    }
    // Write transition case counts to file for runtime profiling.
    {
        ofstream count_log("position_beam_transition_counts.txt", ios::out);
        if (count_log.is_open()) {
            count_log << "# PositionBeamDP transition case counts (candidates examined per state transition)\n";
            count_log << "N_EtoN\t" << cnt_N_EtoN << "\n";
            count_log << "NtoC\t" << cnt_NtoC << "\n";
            count_log << "S_EtoS\t" << cnt_S_EtoS << "\n";
            count_log << "CS_to_C\t" << cnt_CS_to_C << "\n";
            count_log << "S_CS_to_C\t" << cnt_S_CS_to_C << "\n";
            count_log << "C_S_to_CS\t" << cnt_C_S_to_CS << "\n";
            count_log << "C_to_C\t" << cnt_C_to_C << "\n";
            count_log << "S_C_to_C\t" << cnt_S_C_to_C << "\n";
            count_log << "Multi_EtoMulti\t" << cnt_Multi_EtoMulti << "\n";
            count_log << "MultitoC\t" << cnt_MultitoC << "\n";
            count_log << "C_to_M1\t" << cnt_C_to_M1 << "\n";
            count_log << "M1_C_to_M2\t" << cnt_M1_C_to_M2 << "\n";
            count_log << "M1_EtoM1\t" << cnt_M1_EtoM1 << "\n";
            count_log << "M2_to_M1\t" << cnt_M2_to_M1 << "\n";
            count_log << "F_EtoF\t" << cnt_F_EtoF << "\n";
            count_log << "C_to_F\t" << cnt_C_to_F << "\n";
            count_log << "F_C_to_F\t" << cnt_F_C_to_F << "\n";
            count_log << "SC_left_build\t" << cnt_SC_left_build << "\n";
            count_log << "SCLeft_to_C\t" << cnt_SCLeft_to_C << "\n";
            count_log << "blk1_rightbulge\t" << cnt_blk1_rightbulge << "\n";
            count_log << "blk1_internal\t" << cnt_blk1_internal << "\n";
            count_log.close();
        }
    }
    double fill_total_ms = duration<double, milli>(high_resolution_clock::now() - t_fill_start).count();
    cerr << "[PositionBeamDP profile] fill_total_ms=" << fill_total_ms
         << " N_S=" << ms_N_S << " C=" << ms_C << " Multi=" << ms_Multi << " F=" << ms_F << endl;
    cerr << "[PositionBeamDP profile] C breakdown: CS_to_C=" << ms_CS_to_C << " C_S_to_CS=" << ms_C_S_to_CS
         << " C_to_C=" << ms_C_to_C << " S_C_to_C=" << ms_S_C_to_C << " ms" << endl;
    if (std::getenv("DERNA_BLK1_DEBUG")) {
        cerr << "[PositionBeamDP blk1] kept=" << blk1_c_kv_kept_dbg
             << " filtered=" << blk1_c_kv_filtered_dbg
             << " (filtered fraction = "
             << (blk1_c_kv_kept_dbg + blk1_c_kv_filtered_dbg ?
                 100.0 * blk1_c_kv_filtered_dbg / (blk1_c_kv_kept_dbg + blk1_c_kv_filtered_dbg) : 0.0)
             << "%)" << endl;
    }
    // Close codon-agg diagnostics log if open.
    agg_diag_close();
}

// -----------------------------------------------------------------------------
// Traceback: reconstruct sequence and structure from filled tables
// -----------------------------------------------------------------------------
void traceback_position_beam_tables(const AllTablesDerna& tables, int n, const vector<int>& protein,
                                    vector<int>& nucle_seq, vector<int>& codon_selection,
                                    vector<bond>& bp_bond) {
    auto t_tb_start = high_resolution_clock::now();
    const int nuc_len = 3 * n;
    nucle_seq.assign(nuc_len, -1);
    codon_selection.assign(n, -1);
    bp_bond.clear();
    // Track which positions are already paired to avoid cross-level duplicate base pairs.
    vector<bool> paired(nuc_len, false);

    const auto& tab_n = tables.bestN;
    const auto& tab_s = tables.bestS;
    const auto& tab_f = tables.bestF;
    const auto& tab_c = tables.bestC;
    const auto& tab_cs = tables.bestCS;
    const auto& tab_m1 = tables.bestM1;
    const auto& tab_m2 = tables.bestM2;
    const auto& tab_multi = tables.bestMulti;
    const auto& tab_sc_left = tables.bestSCLeft;

    if (nuc_len == 0 || tab_f.size() <= (size_t)(nuc_len - 1)) return;

    int final_pos = nuc_len - 1;
    double best_sc = 1e30;  // minimize (match Zuker)
    int best_key = -1;
    for (const auto& kv : tab_f[final_pos]) {
        const BeamEntry& e = kv.second;
        int left_bound = sigma(e.a, e.i);
        int right_bound = sigma(e.b, e.j);
        if (left_bound != 0 || right_bound != final_pos) continue;
        if (e.score < best_sc) { best_sc = e.score; best_key = kv.first; }
    }
    if (best_key < 0 || tab_f[final_pos].count(best_key) == 0) return;

    int start_manner = tab_f[final_pos].at(best_key).backtrace_type;
    struct Item {
        int pos = -1, key = 0, manner = 0, seg_len = -1;
        bool cs_key_mode = false;
        long long cs_o = 0;
        long long cs_ic = -1;
        // Merged-state traceback: if >= 0, prefer the variant with this (x,y) codon pair.
        // This ensures same-codon extension chains follow a consistent codon throughout.
        int want_x = -1, want_y = -1;
    };
    vector<Item> stack;
    stack.push_back({final_pos, best_key, start_manner, -1, false, 0, -1LL, -1, -1});

    int max_steps = nuc_len * 20;
    int steps = 0;
    int skip_key_not_found = 0;

    // Collected chain states for 2-pass bug_check. Populated during the main
    // traceback loop; recomputed at end using fully-resolved nucle_seq.
    struct RecheckStep {
        int step_idx;
        int manner_val;
        int pos;
        double e_mfe;
        // Parent C state
        double pc_mfe;
        int inner_a, inner_b;
        int inner_i, inner_j;
        int inner_x, inner_y;
        // Outer pair and loop geometry
        int outer_a, outer_b;
        int outer_i, outer_j;
        int outer_x, outer_y;
        int seg_len_left;
        int seg_len_right;
    };
    std::vector<RecheckStep> recheck;

    while (!stack.empty() && steps++ < max_steps) {
        Item it = stack.back();
        stack.pop_back();
        int pos = it.pos, key = it.key, manner = it.manner, seg_len = it.seg_len;

        const BeamEntry* e_ptr = nullptr;
        int use_pos = pos;
        int use_seg = seg_len;

        auto get_map = [&](int ptry, int segtry, int manner_q) -> const DernaBeamMap* {
            if (ptry < 0 || (size_t)ptry >= tab_f.size()) return nullptr;
            if (manner_q == (int)DernaManner::MANNER_C_StoCS)
                return nullptr;
            if (manner_q == (int)DernaManner::MANNER_F_EtoF || manner_q == (int)DernaManner::MANNER_F_CtoF || manner_q == (int)DernaManner::MANNER_CtoF)
                return &tab_f[ptry];
            if (manner_q == (int)DernaManner::MANNER_NtoC || manner_q == (int)DernaManner::MANNER_N_EtoN || manner_q == (int)DernaManner::MANNER_NONEtoN)
                return &tab_n[ptry];
            if (manner_q == (int)DernaManner::MANNER_S_EtoS || manner_q == (int)DernaManner::MANNER_NONEtoS)
                return (segtry > 0 && segtry < (int)tab_s[ptry].size()) ? &tab_s[ptry][segtry] : nullptr;
            if (manner_q == (int)DernaManner::MANNER_CStoC || manner_q == (int)DernaManner::MANNER_CtoC ||
                manner_q == (int)DernaManner::MANNER_S_CtoC || manner_q == (int)DernaManner::MANNER_S_CStoC ||
                manner_q == (int)DernaManner::MANNER_MultitoC ||
                manner_q == (int)DernaManner::MANNER_C_StoC || manner_q == (int)DernaManner::MANNER_S_C_StoC ||
                manner_q == (int)DernaManner::MANNER_SpecialHP)
                return &tab_c[ptry];
            if (manner_q == (int)DernaManner::MANNER_CtoM1 ||
                manner_q == (int)DernaManner::MANNER_M1_EtoM1 ||
                manner_q == (int)DernaManner::MANNER_M2toM1)
                return &tab_m1[ptry];
            if (manner_q == (int)DernaManner::MANNER_M1_CtoM2)
                return &tab_m2[ptry];
            if (manner_q == (int)DernaManner::MANNER_M2toMulti ||
                manner_q == (int)DernaManner::MANNER_S_M2toMulti ||
                manner_q == (int)DernaManner::MANNER_Multi_EtoMulti)
                return &tab_multi[ptry];
            return nullptr;
        };

        auto try_find = [&](int ptry, int segtry) -> bool {
            if (manner == (int)DernaManner::MANNER_C_StoCS && it.cs_key_mode) {
                if (ptry < 0 || (size_t)ptry >= tab_cs.size()) return false;
                int cs_k = (int)it.cs_o;  // cs_o stores the int CS key
                auto jt = tab_cs[ptry].find(cs_k);
                if (jt != tab_cs[ptry].end()) { e_ptr = &jt->second; use_pos = ptry; use_seg = segtry; return true; }
                return false;
            }
            const DernaBeamMap* mp = get_map(ptry, segtry, manner);
            if (!mp) return false;
            auto it2 = mp->find(key);
            if (it2 != mp->end()) { e_ptr = &it2->second; use_pos = ptry; use_seg = segtry; return true; }
            return false;
        };

        bool found = try_find(pos, seg_len);
        if (!found) found = try_find(pos - 1, seg_len);

        if (!found) {
            if (!it.cs_key_mode) {
                // C table uses nuc_key_c keys; all other tables use index_struct keys.
                bool is_c_manner = (manner == (int)DernaManner::MANNER_NtoC ||
                                    manner == (int)DernaManner::MANNER_CStoC ||
                                    manner == (int)DernaManner::MANNER_CtoC ||
                                    manner == (int)DernaManner::MANNER_S_CtoC ||
                                    manner == (int)DernaManner::MANNER_S_CStoC ||
                                    manner == (int)DernaManner::MANNER_C_StoC ||
                                    manner == (int)DernaManner::MANNER_S_C_StoC ||
                                    manner == (int)DernaManner::MANNER_SCLefttoC ||
                                    manner == (int)DernaManner::MANNER_MultitoC ||
                                    manner == (int)DernaManner::MANNER_SpecialHP);
                int pkey = is_c_manner ? nuc_key_c_to_close(key, nuc_len) : codon_beam_key_left_pos(key);
                if (pkey >= 0 && pkey < nuc_len) found = try_find(pkey, seg_len);

                if (!found && (manner == (int)DernaManner::MANNER_S_EtoS || manner == (int)DernaManner::MANNER_NONEtoS) &&
                    pkey >= 0 && pkey < nuc_len) {
                    for (int s = 1; s < (int)tab_s[pkey].size(); ++s) {
                        const DernaBeamMap* mp = get_map(pkey, s, manner);
                        if (!mp) continue;
                        auto it3 = mp->find(key);
                        if (it3 != mp->end()) { e_ptr = &it3->second; use_pos = pkey; use_seg = s; found = true; break; }
                    }
                }
            }
        }

        if (!found || !e_ptr) {
            skip_key_not_found++;
            continue;
        }

        pos = use_pos;
        seg_len = use_seg;
        const BeamEntry& e = *e_ptr;
        int a = e.a, b = e.b, i = e.i, j = e.j;
        if (a < 0 || a >= n || b < 0 || b >= n) continue;

        // Optional step-by-step trace: emit (step, push_manner, created_manner,
        // pos, a, b, i, j, state.mfe, state.cai) so we can see the real chain
        // of states producing the final DP score.  Enable with
        // environment variable DERNA_TB_TRACE=1.  Note: `manner` is how this
        // state was pushed onto the traceback stack (always CtoC for C
        // descendents), while `e.backtrace_type` is the ACTUAL manner that
        // created this state — the one we care about for audit.
        {
            static const bool tb_trace_enabled = (std::getenv("DERNA_TB_TRACE") != nullptr);
            if (tb_trace_enabled) {
                std::cerr << "[tb_trace] step=" << steps
                          << " pos=" << pos
                          << " push_m=" << manner_name(manner)
                          << " created_m=" << manner_name(e.backtrace_type)
                          << " a=" << a << " b=" << b
                          << " i=" << i << " j=" << j
                          << " x=" << e.x << " y=" << e.y
                          << " e.mfe=" << e.mfe << " e.cai=" << e.cai
                          << " e.score=" << e.score
                          << std::endl;
            }
        }

        // Merged-state traceback: resolve (x, y, bt_info, backtrace_type) from the variant
        // that matches the expected codon pair from the caller (want_x / want_y).
        // Falls back to the top-level best variant if no match is found.
        int x = e.x, y = e.y;
        int mt = e.backtrace_type;
        const BtInfo* bt_ptr = &e.bt_info;
        // Per-variant cs_* metadata (see XYVariant): the top-level cs_* fields only
        // reflect the best variant, so a non-best variant selected via want_x/want_y
        // must use its own predecessor metadata.
        long long var_cs_pack_outer = e.cs_pack_outer;
        int var_cs_right_len = e.cs_right_len;
        int var_cs_single_start = e.cs_single_start;
        int var_cs_inner_left = e.cs_inner_left;
        int var_cs_inner_right = e.cs_inner_right;
        double var_mfe = e.mfe;  // trace-only: follow variant's mfe (not top-level) for bug_check2
        bool variant_matched = (it.want_x < 0 && it.want_y < 0);
        if (it.want_x >= 0 || it.want_y >= 0) {
            for (const auto& v : e.variants) {
                if ((it.want_x < 0 || v.x == it.want_x) &&
                    (it.want_y < 0 || v.y == it.want_y)) {
                    x = v.x; y = v.y;
                    mt = v.backtrace_type;
                    bt_ptr = &v.bt_info;
                    var_cs_pack_outer   = v.cs_pack_outer;
                    var_cs_right_len    = v.cs_right_len;
                    var_cs_single_start = v.cs_single_start;
                    var_cs_inner_left   = v.cs_inner_left;
                    var_cs_inner_right  = v.cs_inner_right;
                    var_mfe             = v.mfe;
                    variant_matched = true;
                    break;
                }
            }
            if (!variant_matched) {
                // Caller demanded a specific (x,y) variant that was evicted from
                // this state. Falling through to the top-level best variant would
                // overwrite codon_selection with the wrong codon and corrupt the
                // reconstructed structure. Skip this branch instead — better an
                // incomplete traceback than non-WC base pairs.
                skip_key_not_found++;
                continue;
            }
        }
        const BtInfo& bt = *bt_ptr;

        // === DERNA_TB_TRACE 2-pass bug_check (pass 1: collect) ===
        // Save chain states on C_StoC / S_CtoC / S_C_StoC so we can re-check
        // with fully-resolved nucle_seq after the traceback finishes.
        //   C_StoC  (right bulge):       bt[0]=c_key, bt[1]=s_right_key
        //   S_CtoC  (left bulge):        bt[0]=s_key, bt[1]=c_key
        //   S_C_StoC (internal loop):    bt[0]=s_left_key, bt[1]=c_key, bt[2]=s_right_key
        {
            static const bool tb_trace_enabled = (std::getenv("DERNA_TB_TRACE") != nullptr);
            if (tb_trace_enabled) {
                int ct = e.backtrace_type;
                bool is_rb = (ct == (int)DernaManner::MANNER_C_StoC);
                bool is_lb = (ct == (int)DernaManner::MANNER_S_CtoC);
                bool is_il = (ct == (int)DernaManner::MANNER_S_C_StoC);
                bool is_hp = (ct == (int)DernaManner::MANNER_NtoC);
                bool is_stk = (ct == (int)DernaManner::MANNER_CtoC);
                bool is_ilcs = (ct == (int)DernaManner::MANNER_S_CStoC);
                if (is_ilcs) {
                    // S_CStoC: closes an internal loop via CS composite.
                    // Inner C lives inside the CS at tab_cs[pos-1][cs_pack_outer].
                    // Walk through to get inner C's mfe. bt_info = {s_left_key, 0, slv.x, slv.y}.
                    int cs_pos = pos - 1;
                    long long cs_key = var_cs_pack_outer;
                    if (cs_pos >= 0 && (size_t)cs_pos < tab_cs.size()
                            && cs_key >= 0
                            && tab_cs[cs_pos].count((int)cs_key)) {
                        const BeamEntry& cs = tab_cs[cs_pos].at((int)cs_key);
                        // CS composite contains inner_C at tab_c[cs.cs_inner_right]
                        // Wait: in C_StoCS, cs_inner_right = sigma(c_ent.b, c_ent.j),
                        // i.e., where the inner C's right boundary lives.
                        int inner_c_pos = cs.cs_inner_right;
                        int inner_c_key = cs.cs_inner_key;
                        int slen_right = cs.cs_right_len;
                        int slen_left  = var_cs_right_len;  // for S_CStoC, cs_right_len of the outer state stores left seg
                        if (inner_c_pos >= 0 && (size_t)inner_c_pos < tab_c.size()
                                && inner_c_key >= 0
                                && tab_c[inner_c_pos].count(inner_c_key)) {
                            const BeamEntry& pc = tab_c[inner_c_pos].at(inner_c_key);
                            RecheckStep rs;
                            rs.step_idx = steps;
                            rs.manner_val = ct;
                            rs.pos = pos;
                            rs.e_mfe = var_mfe;
                            rs.pc_mfe = pc.mfe;
                            rs.inner_a = pc.a; rs.inner_b = pc.b;
                            rs.inner_i = pc.i; rs.inner_j = pc.j;
                            rs.inner_x = pc.x;  rs.inner_y = pc.y;
                            rs.outer_a = a; rs.outer_b = b;
                            rs.outer_i = i; rs.outer_j = j;
                            rs.outer_x = x; rs.outer_y = y;
                            rs.seg_len_left = slen_left;
                            rs.seg_len_right = slen_right;
                            recheck.push_back(rs);
                        }
                    }
                }
                else if (is_stk) {
                    // Stacking: child.mfe = parent.mfe + stackE. Need parent C.
                    // bt_info = {c_key, inner_x, inner_y}; parent C at tab_c[pos-1].
                    int c_close = pos - 1;
                    if (c_close >= 0 && (size_t)c_close < tab_c.size()
                            && (int)bt.size() >= 1
                            && tab_c[c_close].count(bt[0])) {
                        const BeamEntry& pc = tab_c[c_close].at(bt[0]);
                        int want_cx = (int)bt.size() >= 2 ? bt[1] : -1;
                        int want_cy = (int)bt.size() >= 3 ? bt[2] : -1;
                        double pc_mfe = pc.mfe;
                        for (const auto& v : pc.variants) {
                            if ((want_cx < 0 || v.x == want_cx) &&
                                (want_cy < 0 || v.y == want_cy)) { pc_mfe = v.mfe; break; }
                        }
                        int cvx = (want_cx >= 0 ? want_cx : pc.x);
                        int cvy = (want_cy >= 0 ? want_cy : pc.y);
                        RecheckStep rs;
                        rs.step_idx = steps;
                        rs.manner_val = ct;
                        rs.pos = pos;
                        rs.e_mfe = var_mfe;
                        rs.pc_mfe = pc_mfe;
                        rs.inner_a = pc.a; rs.inner_b = pc.b;
                        rs.inner_i = pc.i; rs.inner_j = pc.j;
                        rs.inner_x = cvx;  rs.inner_y = cvy;
                        rs.outer_a = a; rs.outer_b = b;
                        rs.outer_i = i; rs.outer_j = j;
                        rs.outer_x = x; rs.outer_y = y;
                        rs.seg_len_left = 0;
                        rs.seg_len_right = 0;
                        recheck.push_back(rs);
                    }
                }
                else if (is_hp) {
                    // Hairpin close: the child's mfe should equal Zuker::hairpin_loop(...)
                    // (or the special tri/tetra/hexaloop energy if applicable).
                    // Save for pass-2 recomputation with full nucle_seq.
                    RecheckStep rs;
                    rs.step_idx = steps;
                    rs.manner_val = ct;
                    rs.pos = pos;
                    rs.e_mfe = var_mfe;
                    rs.pc_mfe = 0;  // no parent C for hairpin
                    rs.inner_a = -1; rs.inner_b = -1;
                    rs.inner_i = -1; rs.inner_j = -1;
                    rs.inner_x = -1; rs.inner_y = -1;
                    rs.outer_a = a; rs.outer_b = b;
                    rs.outer_i = i; rs.outer_j = j;
                    rs.outer_x = x; rs.outer_y = y;
                    rs.seg_len_left = 0; rs.seg_len_right = 0;
                    recheck.push_back(rs);
                }
                else if (is_rb || is_lb || is_il) {
                    int c_slot;
                    int seg_len_left, seg_len_right;
                    // The inner C's position in tab_c differs per manner:
                    //   C_StoC (right bulge)    → inner right = cs_inner_left
                    //   S_C_StoC (internal loop) → inner right = cs_inner_left
                    //   S_CtoC (left bulge)      → inner right = pos - 1 (adjacent to outer right)
                    int closed_right;
                    if (is_rb) { c_slot = 0; seg_len_left = 0;                   seg_len_right = var_cs_right_len; closed_right = var_cs_inner_left; }
                    else if (is_lb) { c_slot = 1; seg_len_left = var_cs_right_len; seg_len_right = 0;                   closed_right = pos - 1; }
                    else             { c_slot = 1; seg_len_left = var_cs_inner_right; seg_len_right = var_cs_right_len; closed_right = var_cs_inner_left; }
                    if (closed_right >= 0 && (size_t)closed_right < tab_c.size()
                            && (int)bt.size() > c_slot
                            && tab_c[closed_right].count(bt[c_slot])) {
                        const BeamEntry& pc = tab_c[closed_right].at(bt[c_slot]);
                        int want_cx = -1, want_cy = -1;
                        if (is_rb && (int)bt.size() >= 4) { want_cx = bt[2]; want_cy = bt[3]; }
                        if (is_lb && (int)bt.size() >= 4) { want_cx = bt[2]; want_cy = bt[3]; }
                        if (is_il && (int)bt.size() >= 7) { want_cx = bt[5]; want_cy = bt[6]; }
                        double pc_mfe = pc.mfe;
                        for (const auto& v : pc.variants) {
                            if ((want_cx < 0 || v.x == want_cx) &&
                                (want_cy < 0 || v.y == want_cy)) { pc_mfe = v.mfe; break; }
                        }
                        int cvx = (want_cx >= 0 ? want_cx : pc.x);
                        int cvy = (want_cy >= 0 ? want_cy : pc.y);
                        RecheckStep rs;
                        rs.step_idx = steps;
                        rs.manner_val = ct;
                        rs.pos = pos;
                        rs.e_mfe = var_mfe;
                        rs.pc_mfe = pc_mfe;
                        rs.inner_a = pc.a; rs.inner_b = pc.b;
                        rs.inner_i = pc.i; rs.inner_j = pc.j;
                        rs.inner_x = cvx;  rs.inner_y = cvy;
                        rs.outer_a = a; rs.outer_b = b;
                        rs.outer_i = i; rs.outer_j = j;
                        rs.outer_x = x; rs.outer_y = y;
                        rs.seg_len_left = seg_len_left;
                        rs.seg_len_right = seg_len_right;
                        recheck.push_back(rs);
                    }
                }
            }
        }

        int li = sigma(a, i), rj = sigma(b, j);
        // Codon-consistency check: detect overwrites of codon_selection[aa].
        {
            static const bool tb_trace_enabled = (std::getenv("DERNA_TB_TRACE") != nullptr);
            if (tb_trace_enabled) {
                if (codon_selection[a] != -1 && codon_selection[a] != x) {
                    std::cerr << "[codon_conflict] step=" << steps
                              << " manner=" << manner_name(e.backtrace_type)
                              << " aa_a=" << a << " prior=" << codon_selection[a] << " new=" << x << "\n";
                }
                if (a != b && codon_selection[b] != -1 && codon_selection[b] != y) {
                    std::cerr << "[codon_conflict] step=" << steps
                              << " manner=" << manner_name(e.backtrace_type)
                              << " aa_b=" << b << " prior=" << codon_selection[b] << " new=" << y << "\n";
                }
            }
        }
        codon_selection[a] = x;
        codon_selection[b] = y;
        if (li >= 0 && li < nuc_len && protein[a] >= 0 && protein[a] < 20 && x >= 0 && x < n_codon[protein[a]])
            nucle_seq[li] = nucleotides[protein[a]][x][i];
        if (rj >= 0 && rj < nuc_len && protein[b] >= 0 && protein[b] < 20 && y >= 0 && y < n_codon[protein[b]])
            nucle_seq[rj] = nucleotides[protein[b]][y][j];

        // CtoF is NOT in is_closed: the base pair is already recorded when the underlying
        // CtoC/NtoC entry is traced; including CtoF here causes double-recording of bp_bond.
        bool is_closed = (manner == (int)DernaManner::MANNER_NtoC || manner == (int)DernaManner::MANNER_CStoC ||
                          manner == (int)DernaManner::MANNER_CtoC || manner == (int)DernaManner::MANNER_S_CtoC ||
                          manner == (int)DernaManner::MANNER_S_CStoC || manner == (int)DernaManner::MANNER_MultitoC ||
                          manner == (int)DernaManner::MANNER_C_StoC || manner == (int)DernaManner::MANNER_S_C_StoC ||
                          manner == (int)DernaManner::MANNER_SpecialHP);
        if (is_closed && li >= 0 && li < nuc_len && rj >= 0 && rj < nuc_len && li < rj) {
            if (!paired[li] && !paired[rj]) {
                bp_bond.push_back({li, rj});
                paired[li] = true;
                paired[rj] = true;
            } else {
                // Pair-drop: DP accumulated energy for this pair but the position
                // is already occupied by an earlier-traced pair. This is a prime
                // suspect for the DP-vs-evaluate phantom-energy discrepancy.
                static const bool tb_trace_enabled = (std::getenv("DERNA_TB_TRACE") != nullptr);
                if (tb_trace_enabled) {
                    std::cerr << "[pair_dropped] step=" << steps
                              << " manner=" << manner_name(manner)
                              << " created_m=" << manner_name(e.backtrace_type)
                              << " li=" << li << " rj=" << rj
                              << " paired[li]=" << (paired[li] ? "Y" : "N")
                              << " paired[rj]=" << (paired[rj] ? "Y" : "N")
                              << " e.mfe=" << e.mfe
                              << "\n";
                }
            }
        }

        // Helper: push an item; for same-table extension transitions propagate (x,y) as hint
        // so the variant-resolution logic at the source chooses the right codon pair.
        auto push_ext = [&](int p, int k, int m, int sl) {
            stack.push_back({p, k, m, sl, false, 0LL, -1LL, x, y});
        };

        if (mt == (int)DernaManner::MANNER_F_EtoF && bt.size() >= 1 && pos > 0) {
            int inner_x = (bt.size() >= 2) ? bt[1] : x;
            int inner_y = (bt.size() >= 3) ? bt[2] : y;
            stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_F_EtoF, -1, false, 0LL, -1LL, inner_x, inner_y});
        }

        else if (mt == (int)DernaManner::MANNER_F_CtoF && bt.size() >= 2) {
            // bt_info = {f_key, c_key, seam_pos, f_var.y, cvar.x, cvar.y}
            int seam = (bt.size() >= 3) ? bt[2] : (sigma(e.a, e.i) - 1);
            int inner_f_y = (bt.size() >= 4) ? bt[3] : -1;
            int cx = (bt.size() >= 6) ? bt[4] : -1;
            int cy = (bt.size() >= 6) ? bt[5] : -1;
            if (seam >= 0 && (size_t)seam < tab_f.size() && tab_f[seam].count(bt[0]))
                stack.push_back({seam, bt[0], (int)DernaManner::MANNER_F_EtoF, -1, false, 0LL, -1LL, x, inner_f_y});
            int c_close = nuc_key_c_to_close(bt[1], nuc_len);
            if (c_close >= 0 && (size_t)c_close < tab_c.size())
                stack.push_back({c_close, bt[1], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
        }
        else if (mt == (int)DernaManner::MANNER_CtoF) {
            // The F-table key is index_struct; recompute the C-table key (nuc_key_c_nuc) from entry fields.
            // nucleotides[] is a global array; protein[] is the const ref parameter.
            int nuc_lo_tb = nucleotides[protein[a]][x][i];
            int nuc_ro_tb = nucleotides[protein[b]][y][j];
            int c_key = nuc_key_c_nuc(sigma(a, i), sigma(b, j), nuc_lo_tb, nuc_ro_tb, nuc_len);
            int c_close = nuc_key_c_to_close(c_key, nuc_len);
            if (c_close >= 0 && (size_t)c_close < tab_c.size())
                push_ext(c_close, c_key, (int)DernaManner::MANNER_CtoC, -1);
        }

        else if (mt == (int)DernaManner::MANNER_N_EtoN && bt.size() >= 1 && pos > 0) {
            int inner_x = (bt.size() >= 2) ? bt[1] : x;
            int inner_y = (bt.size() >= 3) ? bt[2] : y;
            stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_N_EtoN, -1, false, 0LL, -1LL, inner_x, inner_y});
        }
        else if (mt == (int)DernaManner::MANNER_NONEtoN || mt == (int)DernaManner::MANNER_NONEtoS) { /* terminal */ }

        else if (mt == (int)DernaManner::MANNER_NtoC && bt.size() >= 1 && pos > 0) {
            int inner_x = (bt.size() >= 2) ? bt[1] : -1;
            int inner_y = (bt.size() >= 3) ? bt[2] : -1;
            stack.push_back({pos - 1, bt[0], (int)DernaManner::MANNER_N_EtoN, -1, false, 0LL, -1LL, inner_x, inner_y});
        }

        else if (mt == (int)DernaManner::MANNER_SpecialHP) {
            // Special-hairpin terminal: bt_info stores codon choices for aa's in [a..b].
            // Apply them directly to codon_selection (main codon_selection of boundary aa's was already
            // set above; interior aa codons must be recorded here). No further stack entries needed.
            for (int kk = 0; kk < (int)bt.size(); ++kk) {
                int aa = e.a + kk;
                if (aa < 0 || aa >= n) continue;
                int pa2 = protein[aa];
                int cc = bt[kk];
                if (cc < 0 || cc >= n_codon[pa2]) continue;
                codon_selection[aa] = cc;
                for (int s = 0; s < 3; ++s) {
                    int nt = 3 * aa + s;
                    if (nt >= 0 && nt < nuc_len) nucle_seq[nt] = nucleotides[pa2][cc][s];
                }
            }
        }

        else if (mt == (int)DernaManner::MANNER_S_EtoS && bt.size() >= 1 && pos > 0 && seg_len > 1) {
            if (seg_len - 1 < (int)tab_s[pos - 1].size() && tab_s[pos - 1][seg_len - 1].count(bt[0])) {
                int inner_x = (bt.size() >= 2) ? bt[1] : x;
                int inner_y = (bt.size() >= 3) ? bt[2] : y;
                stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_S_EtoS, seg_len-1, false, 0LL, -1LL, inner_x, inner_y});
            }
        }
        else if (mt == (int)DernaManner::MANNER_S_EtoS && bt.size() >= 1 && pos > 0 && seg_len == 1) { /* S len 1 terminal */ }

        else if (mt == (int)DernaManner::MANNER_CStoC && pos > 0)
            stack.push_back({pos - 1, 0, (int)DernaManner::MANNER_C_StoCS, -1, true, var_cs_pack_outer, 0LL, -1, -1});

        else if (mt == (int)DernaManner::MANNER_C_StoCS && bt.size() >= 2) {
            // bt_info = {c_key, s_key, cvar.x, cvar.y, svar.x, svar.y}
            int closed_right = sigma(e.a, e.i) - 1;
            int seg = pos - closed_right;
            int cx = (bt.size() >= 4) ? bt[2] : -1;
            int cy = (bt.size() >= 4) ? bt[3] : -1;
            int sx = (bt.size() >= 6) ? bt[4] : -1;
            int sy = (bt.size() >= 6) ? bt[5] : -1;
            int c_close = nuc_key_c_to_close(bt[0], nuc_len);  // bt[0] is a nuc_key_c_nuc C-table key
            if (c_close >= 0 && (size_t)c_close < tab_c.size())
                stack.push_back({c_close, bt[0], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
            if (seg > 0 && seg < (int)tab_s[pos].size() && tab_s[pos][seg].count(bt[1]))
                stack.push_back({pos, bt[1], (int)DernaManner::MANNER_S_EtoS, seg, false, 0LL, -1LL, sx, sy});
        }
        else if (mt == (int)DernaManner::MANNER_CtoC && bt.size() >= 1 && pos > 0) {
            // Use stored inner (var.x, var.y) from bt_info to correctly resolve the predecessor C variant.
            int inner_x = (bt.size() >= 2) ? bt[1] : x;
            int inner_y = (bt.size() >= 3) ? bt[2] : y;
            stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, inner_x, inner_y});
        }

        else if (mt == (int)DernaManner::MANNER_SCLefttoC) {
            // SCLefttoC: outer pair closes at (outer_left, pos). With dual-keyed SCLeft,
            // the C entry's bt_info stores the SCLeft variant identifier {sleft_x_v, inner_y_v}
            // (so traceback can find the correct merged-table variant).
            int scl_pos = pos - 1;
            int scl_key = (int)var_cs_pack_outer;
            if (scl_pos >= 0 && (size_t)scl_pos < tab_sc_left.size() &&
                tab_sc_left[scl_pos].count(scl_key)) {
                const BeamEntry& scl = tab_sc_left[scl_pos].at(scl_key);
                int scl_var_x = (bt.size() >= 1) ? bt[0] : -1;
                int scl_var_y = (bt.size() >= 2) ? bt[1] : -1;
                // Resolve specific variant inside the merged SCLeft entry.
                const XYVariant* sv = nullptr;
                if (scl_var_x >= 0 && scl_var_x < 6 && scl_var_y >= 0 && scl_var_y < 6) {
                    int vi = scl.variant_idx[scl_var_x][scl_var_y];
                    if (vi >= 0 && vi < (int)scl.variants.size()) sv = &scl.variants[vi];
                }
                if (!sv && !scl.variants.empty()) sv = &scl.variants.front();  // fallback (shouldn't be needed)
                // SCLeft variant bt_info = {c_key, s_key, cvar.x, svar.y}
                const BtInfo& var_bt = sv ? sv->bt_info : scl.bt_info;
                int c_key_inner = (var_bt.size() >= 1) ? var_bt[0] : scl.cs_inner_key;
                int s_key       = (var_bt.size() >= 2) ? var_bt[1] : scl.cs_right_s_key;
                int cx          = (var_bt.size() >= 3) ? var_bt[2] : -1;       // cvar.x (internal)
                int sy          = (var_bt.size() >= 4) ? var_bt[3] : -1;       // svar.y (internal)
                int sx          = sv ? sv->x : -1;                              // svar.x = variant's x
                int cy          = sv ? sv->y : -1;                              // cvar.y = variant's y
                // Inner C at scl_pos (same position as SCLeft)
                if (c_key_inner >= 0 && (size_t)scl_pos < tab_c.size() && tab_c[scl_pos].count(c_key_inner))
                    stack.push_back({scl_pos, c_key_inner, (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
                // S_left at seam_pos = cs_inner_left - 1
                int seam_pos = scl.cs_inner_left - 1;
                int seg = scl.cs_right_len;
                if (seam_pos >= 0 && seg > 0 && (size_t)seam_pos < tab_s.size() &&
                    seg < (int)tab_s[seam_pos].size() && tab_s[seam_pos][seg].count(s_key))
                    stack.push_back({seam_pos, s_key, (int)DernaManner::MANNER_S_EtoS, seg, false, 0LL, -1LL, sx, sy});
            }
        }

        else if (mt == (int)DernaManner::MANNER_S_CtoC && bt.size() >= 2) {
            // bt_info = {s_key, c_key, cvar.x, cvar.y, s_x, s_y}; cs_inner_left = inner_left; cs_right_len = seg_len (left_unpaired)
            // Left bulge (n2=0): inner_right = pos-1 (directly adjacent to outer_right=pos)
            int cx = (bt.size() >= 4) ? bt[2] : -1;
            int cy = (bt.size() >= 4) ? bt[3] : -1;
            int sx = (bt.size() >= 6) ? bt[4] : -1;
            int sy = (bt.size() >= 6) ? bt[5] : -1;
            int c_close = pos - 1;
            if (c_close >= 0 && (size_t)c_close < tab_c.size() && tab_c[c_close].count(bt[1]))
                stack.push_back({c_close, bt[1], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
            int s_inner_left = var_cs_inner_left;
            int s_seg_len    = var_cs_right_len;
            int s_seam_pos   = s_inner_left - 1;
            if (s_seam_pos >= 0 && s_seg_len > 0 && s_seam_pos < (int)tab_s.size() &&
                s_seg_len < (int)tab_s[s_seam_pos].size() && tab_s[s_seam_pos][s_seg_len].count(bt[0]))
                stack.push_back({s_seam_pos, bt[0], (int)DernaManner::MANNER_S_EtoS, s_seg_len, false, 0LL, -1LL, sx, sy});
        }
        else if (mt == (int)DernaManner::MANNER_S_CStoC && bt.size() >= 1) {
            // bt_info = {s_left_key, 0, slv.x, slv.y, best_csvar_y, seg_len_left}
            int cs_pos = pos - 1;
            int cs_want_y = ((int)bt.size() >= 5) ? bt[4] : -1;
            if (cs_pos >= 0 && (size_t)cs_pos < tab_cs.size())
                stack.push_back({cs_pos, 0, (int)DernaManner::MANNER_C_StoCS, -1, true, var_cs_pack_outer, 0LL, -1, cs_want_y});

            int sl_x = (bt.size() >= 4) ? bt[2] : -1;
            int sl_y = (bt.size() >= 4) ? bt[3] : -1;
            int seg_len_left = (bt.size() >= 6) ? bt[5] : -1;
            int seam_pos_left = (bt.size() >= 7) ? bt[6] : -1;
            if (seg_len_left > 0 && seam_pos_left >= 0 && (size_t)seam_pos_left < tab_s.size() &&
                seg_len_left < (int)tab_s[seam_pos_left].size() &&
                tab_s[seam_pos_left][seg_len_left].count(bt[0])) {
                stack.push_back({seam_pos_left, bt[0], (int)DernaManner::MANNER_S_EtoS, seg_len_left, false, 0LL, -1LL, sl_x, sl_y});
            } else {
                // Legacy fallback: prior (buggy) scan. Kept for safety so states with
                // the old 4-slot bt_info don't crash. Should not fire once all states
                // are filled by the new code.
                int s_close = codon_beam_key_left_pos(bt[0]);
                if (s_close >= 0 && (size_t)s_close < tab_s.size()) {
                    for (int seg = 1; seg < (int)tab_s[s_close].size(); ++seg) {
                        if (tab_s[s_close][seg].count(bt[0])) {
                            stack.push_back({s_close, bt[0], (int)DernaManner::MANNER_S_EtoS, seg, false, 0LL, -1LL, sl_x, sl_y});
                            break;
                        }
                    }
                }
            }
        }

        else if (mt == (int)DernaManner::MANNER_C_StoC && bt.size() >= 2) {
            // bt_info = {c_key, s_right_key, c_var.x, c_var.y, sr_var.x, sr_var.y}
            int closed_right = var_cs_inner_left;
            int seg_len_right = var_cs_right_len;
            int s_seam_pos = pos - 1;
            int cx = (bt.size() >= 4) ? bt[2] : -1;
            int cy = (bt.size() >= 4) ? bt[3] : -1;
            int srx = (bt.size() >= 6) ? bt[4] : -1;
            int sry = (bt.size() >= 6) ? bt[5] : -1;
            if (closed_right >= 0 && (size_t)closed_right < tab_c.size() && tab_c[closed_right].count(bt[0]))
                stack.push_back({closed_right, bt[0], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
            if (s_seam_pos >= 0 && seg_len_right > 0 && (size_t)s_seam_pos < tab_s.size() &&
                seg_len_right < (int)tab_s[s_seam_pos].size() && tab_s[s_seam_pos][seg_len_right].count(bt[1]))
                stack.push_back({s_seam_pos, bt[1], (int)DernaManner::MANNER_S_EtoS, seg_len_right, false, 0LL, -1LL, srx, sry});
        }
        else if (mt == (int)DernaManner::MANNER_S_C_StoC && bt.size() >= 3) {
            // bt_info = {s_left_key, c_key, s_right_key, sl_var.x, sl_var.y, c_var.x, c_var.y, sr_var.x, sr_var.y}
            int closed_right = var_cs_inner_left;
            int seg_len_right = var_cs_right_len;
            int seam_pos_left = var_cs_single_start;
            int seg_len_left = var_cs_inner_right;
            int s_seam_pos = pos - 1;
            int slx = (bt.size() >= 9) ? bt[3] : -1;
            int sly = (bt.size() >= 9) ? bt[4] : -1;
            int cx  = (bt.size() >= 9) ? bt[5] : -1;
            int cy  = (bt.size() >= 9) ? bt[6] : -1;
            int srx = (bt.size() >= 9) ? bt[7] : -1;
            int sry = (bt.size() >= 9) ? bt[8] : -1;
            if (closed_right >= 0 && (size_t)closed_right < tab_c.size() && tab_c[closed_right].count(bt[1]))
                stack.push_back({closed_right, bt[1], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, cx, cy});
            if (s_seam_pos >= 0 && seg_len_right > 0 && (size_t)s_seam_pos < tab_s.size() &&
                seg_len_right < (int)tab_s[s_seam_pos].size() && tab_s[s_seam_pos][seg_len_right].count(bt[2]))
                stack.push_back({s_seam_pos, bt[2], (int)DernaManner::MANNER_S_EtoS, seg_len_right, false, 0LL, -1LL, srx, sry});
            if (seam_pos_left >= 0 && seg_len_left > 0 && (size_t)seam_pos_left < tab_s.size() &&
                seg_len_left < (int)tab_s[seam_pos_left].size() && tab_s[seam_pos_left][seg_len_left].count(bt[0]))
                stack.push_back({seam_pos_left, bt[0], (int)DernaManner::MANNER_S_EtoS, seg_len_left, false, 0LL, -1LL, slx, sly});
        }

        else if (mt == (int)DernaManner::MANNER_Multi_EtoMulti && bt.size() >= 1 && pos > 0) {
            int inner_x = (bt.size() >= 2) ? bt[1] : x;
            int inner_y = (bt.size() >= 3) ? bt[2] : y;
            stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_Multi_EtoMulti, -1, false, 0LL, -1LL, inner_x, inner_y});
        }
        else if (mt == (int)DernaManner::MANNER_MultitoC && bt.size() >= 1 && pos > 0) {
            if ((size_t)(pos - 1) < tab_multi.size() && tab_multi[pos - 1].count(bt[0])) {
                int inner_x = (bt.size() >= 2) ? bt[1] : x;
                int inner_y = (bt.size() >= 3) ? bt[2] : y;
                stack.push_back({pos - 1, bt[0], (int)DernaManner::MANNER_Multi_EtoMulti, -1, false, 0LL, -1LL, inner_x, inner_y});
            }
        }

        else if (mt == (int)DernaManner::MANNER_CtoM1 && bt.size() >= 1) {
            // CtoM1 bt_info enriched to {c_key, c_var.x, c_var.y} so we recover the
            // exact C variant even when M1's top-level (x,y) was overwritten by a later
            // M1_EtoM1 update at the same key.
            int c_want_x = (bt.size() >= 2) ? bt[1] : x;
            int c_want_y = (bt.size() >= 3) ? bt[2] : y;
            int c_close = nuc_key_c_to_close(bt[0], nuc_len);
            if (c_close >= 0 && (size_t)c_close < tab_c.size())
                stack.push_back({c_close, bt[0], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, c_want_x, c_want_y});
        }
        else if (mt == (int)DernaManner::MANNER_M1_CtoM2 && bt.size() >= 2) {
            int c_left = nuc_key_c_to_left(bt[1], nuc_len);
            int seam_pos_m1 = c_left - 1;
            int m1_want_x = (bt.size() >= 3) ? bt[2] : -1;
            int m1_want_y = (bt.size() >= 4) ? bt[3] : -1;
            if (seam_pos_m1 >= 0 && (size_t)seam_pos_m1 < tab_m1.size() && tab_m1[seam_pos_m1].count(bt[0]))
                stack.push_back({seam_pos_m1, bt[0], (int)DernaManner::MANNER_M1_EtoM1, -1, false, 0LL, -1LL, m1_want_x, m1_want_y});
            int c_close = nuc_key_c_to_close(bt[1], nuc_len);
            int c_want_x = (bt.size() >= 5) ? bt[4] : -1;
            int c_want_y = (bt.size() >= 6) ? bt[5] : -1;
            if (c_close >= 0 && (size_t)c_close < tab_c.size())
                stack.push_back({c_close, bt[1], (int)DernaManner::MANNER_CtoC, -1, false, 0LL, -1LL, c_want_x, c_want_y});
        }
        else if (mt == (int)DernaManner::MANNER_M1_EtoM1 && bt.size() >= 1 && pos > 0) {
            int inner_x = (bt.size() >= 2) ? bt[1] : x;
            int inner_y = (bt.size() >= 3) ? bt[2] : y;
            stack.push_back({pos-1, bt[0], (int)DernaManner::MANNER_M1_EtoM1, -1, false, 0LL, -1LL, inner_x, inner_y});
        }

        else if (mt == (int)DernaManner::MANNER_M2toM1 && bt.size() >= 1) {
            // bt[0] is the M2 key at the SAME position; M1 inherits M2's (x,y) directly.
            push_ext(pos, bt[0], (int)DernaManner::MANNER_M1_CtoM2, -1);
        }
        else if (mt == (int)DernaManner::MANNER_M2toMulti && bt.size() >= 1) {
            // bt[0] is the M2 key at the SAME position; Multi inherits M2's (x,y) directly.
            push_ext(pos, bt[0], (int)DernaManner::MANNER_M1_CtoM2, -1);
        }
        else if (mt == (int)DernaManner::MANNER_S_M2toMulti && bt.size() >= 2) {
            // bt[0] = S key, bt[1] = M2 key at same pos; bt[2] = s_var.y (= m2var.x via seam constraint)
            // bt[3] = seg_len, bt[4] = seam_pos (direct index into tab_s, avoids ambiguous seg scan)
            // Multi.y = m2var.y, so pass (m2var.x=bt[2], m2var.y=y) to find the exact M2 variant.
            int m2_want_x = (bt.size() >= 3) ? bt[2] : -1;
            stack.push_back({pos, bt[1], (int)DernaManner::MANNER_M1_CtoM2, -1, false, 0LL, -1LL, m2_want_x, y});

            int s_inner_y = (bt.size() >= 3) ? bt[2] : -1;  // stored s_var.y for exact inner-S variant lookup
            int seg_len_s = (bt.size() >= 4) ? bt[3] : -1;
            int seam_pos_s = (bt.size() >= 5) ? bt[4] : -1;
            if (seg_len_s > 0 && seam_pos_s >= 0 && (size_t)seam_pos_s < tab_s.size() &&
                seg_len_s < (int)tab_s[seam_pos_s].size() &&
                tab_s[seam_pos_s][seg_len_s].count(bt[0])) {
                stack.push_back({seam_pos_s, bt[0], (int)DernaManner::MANNER_S_EtoS, seg_len_s, false, 0LL, -1LL, x, s_inner_y});
            } else {
                // Legacy fallback — only used for chart items produced before this fix, if any.
                int s_close = codon_beam_key_left_pos(bt[0]);
                if (s_close >= 0 && (size_t)s_close < tab_s.size()) {
                    for (int seg = 1; seg < (int)tab_s[s_close].size(); ++seg) {
                        if (tab_s[s_close][seg].count(bt[0])) {
                            stack.push_back({s_close, bt[0], (int)DernaManner::MANNER_S_EtoS, seg, false, 0LL, -1LL, x, s_inner_y});
                            break;
                        }
                    }
                }
            }
        }
    }

    std::cerr << "PositionBeamDP traceback: skip_key_not_found=" << skip_key_not_found << ", steps=" << steps << std::endl;

    for (int p = 0; p < nuc_len; ++p) {
        if (nucle_seq[p] >= 0) continue;
        int aa = p / 3, slot = p % 3;
        if (aa >= 0 && aa < n && codon_selection[aa] >= 0) {
            int pa = protein[aa];
            if (pa >= 0 && pa < 20 && codon_selection[aa] < n_codon[pa])
                nucle_seq[p] = nucleotides[pa][codon_selection[aa]][slot];
        }
    }

    // === DERNA_TB_TRACE 2-pass bug_check (pass 2: recompute with real adjacent nucleotides) ===
    // nucle_seq is now populated from codon_selection. Re-run score_single_loop_with_mismatch
    // on each saved bulge/IL close using the *actual* adjacent nucleotides the fill would
    // have seen, and report mismatches.
    {
        static const bool tb_trace_enabled = (std::getenv("DERNA_TB_TRACE") != nullptr);
        if (tb_trace_enabled && !recheck.empty()) {
            std::cerr << "\n=== DERNA_TB_TRACE pass-2 bug_check (n=" << recheck.size() << ") ===\n";
            auto nuc_at = [&](int p) -> int {
                if (p < 0 || p >= nuc_len) return -1;
                return nucle_seq[p];
            };
            double total_delta = 0.0;
            int mismatches = 0;
            for (const auto& rs : recheck) {
                int outer_left  = sigma(rs.outer_a, rs.outer_i);
                int outer_right = sigma(rs.outer_b, rs.outer_j);
                int nuc_outer_L = nuc_at(outer_left);
                int nuc_outer_R = nuc_at(outer_right);
                double loop_e = 0.0, expected = 0.0;
                int inner_left = -1, inner_right = -1;
                int nuc_inner_L = -1, nuc_inner_R = -1;
                int nuc_i1 = -1, nuc_jm1 = -1, nuc_pm1 = -1, nuc_qp1 = -1;
                if (rs.manner_val == (int)DernaManner::MANNER_NtoC) {
                    // Hairpin: loop_e = Zuker::hairpin_loop(nuc_outer_L, nuc_outer_R, nuc_i1, nuc_jm1, loop_len)
                    int loop_len = outer_right - outer_left - 1;
                    nuc_i1  = nuc_at(outer_left + 1);
                    nuc_jm1 = nuc_at(outer_right - 1);
                    // Check for special hairpins (triloop/tetraloop/hexaloop)
                    loop_e = (double)Zuker::hairpin_loop(nuc_outer_L, nuc_outer_R, nuc_i1, nuc_jm1, loop_len);
                    // Compare to what evaluate_structure_energy would compute:
                    // it uses hairpinE lookup for len 3/4/6 with actual loop string.
                    expected = loop_e;
                } else {
                    inner_left  = sigma(rs.inner_a, rs.inner_i);
                    inner_right = sigma(rs.inner_b, rs.inner_j);
                    nuc_inner_L = nuc_at(inner_left);
                    nuc_inner_R = nuc_at(inner_right);
                    nuc_i1  = nuc_at(outer_left + 1);
                    nuc_jm1 = nuc_at(outer_right - 1);
                    nuc_pm1 = nuc_at(inner_left - 1);
                    nuc_qp1 = nuc_at(inner_right + 1);
                    loop_e = score_single_loop_with_mismatch(
                        1.0, nuc_outer_L, nuc_outer_R,
                        nuc_inner_L, nuc_inner_R,
                        rs.seg_len_left, rs.seg_len_right,
                        nuc_i1, nuc_jm1, nuc_pm1, nuc_qp1);
                    expected = rs.pc_mfe + loop_e;
                }
                double delta = rs.e_mfe - expected;
                total_delta += delta;
                bool mism = (std::abs(delta) > 0.5);
                if (mism) mismatches++;
                std::cerr << "[bug_check2] step=" << rs.step_idx
                          << " manner=" << manner_name(rs.manner_val)
                          << " pos=" << rs.pos
                          << " outer=(" << outer_left << "," << outer_right << ")"
                          << " inner=(" << inner_left << "," << inner_right << ")"
                          << " seg_L=" << rs.seg_len_left
                          << " seg_R=" << rs.seg_len_right
                          << " nucs_outer=(" << nuc_outer_L << "," << nuc_outer_R << ")"
                          << " nucs_inner=(" << nuc_inner_L << "," << nuc_inner_R << ")"
                          << " adj=(" << nuc_i1 << "," << nuc_jm1 << "," << nuc_pm1 << "," << nuc_qp1 << ")"
                          << " parent_C.mfe=" << rs.pc_mfe
                          << " loop_e=" << loop_e
                          << " expected=" << expected
                          << " e.mfe=" << rs.e_mfe
                          << " delta=" << delta
                          << (mism ? "  <-- MISMATCH" : "")
                          << "\n";
            }
            std::cerr << "=== bug_check2 summary: total_delta=" << total_delta
                      << " mismatches=" << mismatches << "/" << recheck.size() << " ===\n";
        }
    }

    double tb_ms = duration<double, milli>(high_resolution_clock::now() - t_tb_start).count();
    cerr << "[PositionBeamDP profile] traceback_ms=" << tb_ms << endl;
}
