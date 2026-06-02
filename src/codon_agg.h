//
// codon_agg.h — codon-pair pre-aggregation infrastructure for PositionBeamDP.
//
// Step 1-2 of the planned closure refactor: builds a 36-cell "best (mfe, cai)
// per (cvar.x, cvar.y) tuple" lookup over each S/C entry's variants[] at a
// given (i, j). The closure rewrite (Step 3) is NOT included here — these
// tables are populated but not yet consumed by Block 1.
//
// Key insight (from MEMORY.md / refactor plan): each XYVariant already carries
// cumulative (mfe, cai) per codon-pair tuple. Multiple state-instances at the
// same position-key may have variants with the same (cvar.x, cvar.y) but
// different boundary nucleotides — boundary nucs are functionally determined
// by (codon_assignment, position_slot), so collapsing by codon-pair is
// information-preserving for downstream closure energy.
//
// Diagnostic instrumentation is gated on the DERNA_AGG_DIAG env var. Diagnostic
// data is written to output/agg_diagnostics.log when enabled.
//

#ifndef DERNA_CODON_AGG_H
#define DERNA_CODON_AGG_H

#include "default.h"
#include "PositionBeamDP.h"
#include <array>
#include <cstdint>
#include <vector>
#include <unordered_map>

// NCOD: codon count constant. Verified from default.h: BeamEntry::variant_idx[6][6]
// has 36 slots, so cvar.x ∈ [0,6), cvar.y ∈ [0,6). 6 is the max codons per AA
// (e.g. Leu, Arg, Ser have 6 synonyms).
//
// Note: some AAs have fewer codons (e.g. Met has 1, Trp has 1). For such AAs,
// only a subset of the 36 cells will ever be populated — we do NOT consider
// that "low fill rate" as a failure of the aggregation: it's a structural
// reflection of the codon DFA's branching factor.
static constexpr int NCOD = 6;
static constexpr int CODON_PAIR_GRID = NCOD * NCOD;  // = 36

// One cell of the codon-pair-aggregated lookup table.
//
// score: same beam-ranking metric used elsewhere (combined_score = lambda*mfe
// + (lambda-1)*cai, stored on XYVariant.score).
//
// src_state / src_var_idx: back-pointers for traceback. valid only while the
// underlying tab_s[i][j] / tab_c[i][j] entry is alive. Aggregation tables
// must be invalidated/rebuilt if the source map is mutated.
struct CodonAggEntry {
    float mfe;
    float cai;
    float score;
    BeamEntry* src_state;  // back-pointer to source state-instance
    int src_var_idx;       // index into src_state->variants[]
    bool valid;

    static CodonAggEntry invalid() {
        return CodonAggEntry{0.0f, 0.0f, 0.0f, nullptr, -1, false};
    }
};

// 36-cell array indexed as cvar.x * NCOD + cvar.y.
//
// `populated_indices` lists the indices into `cells` that have been populated
// (cell.valid == true) — a sparse iteration view. The closure refactor (Step 3)
// iterates only populated cells (~5/36 typical), so the dense 36-cell scan is
// avoided. Indices are in insertion order (first-write-wins).
struct CodonAggGrid {
    std::array<CodonAggEntry, CODON_PAIR_GRID> cells;
    std::vector<uint8_t> populated_indices;
};

// Per-table containers. Keyed by structural position only (left_pos at this
// fill position) — NOT by the full state-key. This is the load-bearing design
// choice: the closure refactor's collapse hypothesis is that all state-instances
// at the same (left_pos, right_pos) carrying the same (cvar.x, cvar.y) tuple
// are interchangeable for downstream closure energy (because boundary nucs are
// functionally determined by codon assignment + slot). Multiple state-instances
// in tab_c[pos] with DIFFERENT keys but SAME (left_pos, right_pos) collapse
// into one 36-cell grid here.
//
// For C: tab_c[pos] entries differ in (left_pos, nuc_lo, nuc_ro). At fixed
// pos (= right_pos), grouping by left_pos yields up to 16 distinct keys per
// (left_pos), each contributing variants — these all compete to populate the
// same 36-cell grid.
//
// For S: tab_s[pos][seg_len] entries differ in (left_pos, nuc_left, nuc_right).
// At fixed (pos, seg_len), grouping by left_pos collapses the 16 nuc-key
// variants into one 36-cell grid per left_pos.
//
// Aggregation keying:
//   s_codon_agg[pos][seg_len][left_pos] = CodonAggGrid
//   c_codon_agg[pos][left_pos]          = CodonAggGrid
//
// This is the same coordinate Block 1 iterates over (left_pos = outer-left
// nucleotide position, pos = right-closing position).
using CodonAggMap = std::unordered_map<int, CodonAggGrid>;

struct CodonAggTables {
    // Mirrors tab_s shape: [pos][seg_len] -> map<key, grid>
    std::vector<std::vector<CodonAggMap>> s_agg;
    // Mirrors tab_c shape: [pos] -> map<key, grid>
    std::vector<CodonAggMap> c_agg;

    // Diagnostics (only populated when DERNA_AGG_DIAG=1).
    //
    // s_total_states: sum over seg_len of |tab_s[pos][seg_len]| — how many
    //                 state-instances Block 1 currently iterates.
    // s_total_variants: sum over states of |state.variants| — total variant slots.
    // s_cells_populated: # of (left_pos, cvar.x, cvar.y) cells that received a
    //                    variant. Equivalent to "# of distinct codon-pair tuples
    //                    surfaced in the agg table at this pos".
    // s_cells_possible: sum over distinct left_pos at this pos of 36
    //                   = (# distinct grids) * 36.
    //                   Fill rate = s_cells_populated / s_cells_possible.
    // s_score_collisions: # of variants whose target cell was already populated
    //                     (= "redundant" variants from the closure's perspective —
    //                     they offered no new codon-pair tuple). HIGH collision
    //                     rate => the codon-pair collapse is information-preserving
    //                     for many redundant emits.
    // s_score_overwrites: subset of collisions where the new variant beat the
    //                     existing cell's score (= a Pareto distinction would have
    //                     been lost if the agg consumer rejected the second hit
    //                     instead of keeping the better one).
    struct PerPosDiag {
        // S-table diagnostics aggregated across all seg_lens at this pos.
        uint64_t s_total_states = 0;
        uint64_t s_total_variants = 0;
        uint64_t s_cells_populated = 0;
        uint64_t s_cells_possible = 0;
        uint64_t s_score_collisions = 0;
        uint64_t s_score_overwrites = 0;
        // C-table diagnostics at this pos.
        uint64_t c_total_states = 0;
        uint64_t c_total_variants = 0;
        uint64_t c_cells_populated = 0;
        uint64_t c_cells_possible = 0;
        uint64_t c_score_collisions = 0;
        uint64_t c_score_overwrites = 0;
    };
    std::vector<PerPosDiag> diag;

    void resize(int nuc_len, int s_seg_len_cap) {
        s_agg.assign(nuc_len, std::vector<CodonAggMap>(s_seg_len_cap));
        c_agg.assign(nuc_len, CodonAggMap{});
        diag.assign(nuc_len, PerPosDiag{});
    }
};

// True iff DERNA_AGG_DIAG=1 in env. Cached on first call.
bool agg_diag_enabled();

// Build the codon-pair aggregation grid for tab_s[pos][seg_len], keyed by
// left_pos (= sigma(entry.a, entry.i) = 3*entry.a + entry.i). Must be called
// AFTER tab_s[pos][seg_len] has been finally pruned for this position.
// Skips entries with empty variants[]. Skips variants with sentinel x/y
// outside [0,NCOD).
//
// `tab_s_cell` is the underlying DernaBeamMap (i.e., tab_s[pos][seg_len]).
// `agg_cell` is the destination CodonAggMap; cleared on entry.
// `diag` is the per-pos diagnostic accumulator (updated in place if non-null).
void build_s_codon_agg(DernaBeamMap& tab_s_cell,
                       CodonAggMap& agg_cell,
                       CodonAggTables::PerPosDiag* diag);

// Build the codon-pair aggregation grid for tab_c[pos], keyed by left_pos.
// Must be called AFTER tab_c[pos] has been finally pruned for this position.
void build_c_codon_agg(DernaBeamMap& tab_c_cell,
                       CodonAggMap& agg_cell,
                       CodonAggTables::PerPosDiag* diag);

// Write per-pos diagnostic summary to output/agg_diagnostics.log.
// Header is written on first call; subsequent calls append rows.
// No-op if DERNA_AGG_DIAG is not set.
void agg_diag_write_pos(int pos, const CodonAggTables::PerPosDiag& diag);

// Close the diagnostics log file (if open). Call at end of fill_position_beam_tables.
void agg_diag_close();

#endif  // DERNA_CODON_AGG_H
