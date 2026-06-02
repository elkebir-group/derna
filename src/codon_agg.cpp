//
// codon_agg.cpp — implementation of codon-pair pre-aggregation infrastructure.
//
// See codon_agg.h for design rationale. These functions build per-(i,j)
// aggregation grids from the variants[] inside each BeamEntry of tab_s/tab_c.
// They are populated AFTER final pruning at each position. They are NOT yet
// consumed by the closure logic — Block 1 still iterates state-instances.
//

#include "codon_agg.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>

namespace {

// ---- Diagnostic infrastructure ----

bool s_diag_initialized = false;
bool s_diag_enabled = false;
std::ofstream s_diag_out;
std::mutex s_diag_mu;

}  // namespace

bool agg_diag_enabled() {
    if (!s_diag_initialized) {
        const char* e = std::getenv("DERNA_AGG_DIAG");
        s_diag_enabled = (e && e[0] == '1');
        s_diag_initialized = true;
    }
    return s_diag_enabled;
}

static void open_diag_log_if_needed() {
    if (!agg_diag_enabled()) return;
    if (s_diag_out.is_open()) return;
    // Try output/agg_diagnostics.log first (relative to invocation cwd).
    s_diag_out.open("output/agg_diagnostics.log", std::ios::out);
    if (!s_diag_out.is_open()) {
        // Fallback to current working directory.
        s_diag_out.open("agg_diagnostics.log", std::ios::out);
    }
    if (s_diag_out.is_open()) {
        s_diag_out << "# DERNA codon-agg diagnostics\n";
        s_diag_out << "# Columns: pos table total_states total_variants cells_populated"
                      " cells_possible fill_rate score_collisions score_overwrites"
                      " avg_variants_per_state avg_filled_per_state\n";
    }
}

// Update one cell with a candidate (variant). Returns:
//   0 = inserted (was invalid)
//   1 = collision but new score did NOT beat existing
//   2 = collision AND new score beat existing (overwrite)
static inline int update_agg_cell(CodonAggGrid& grid, int idx,
                                  float mfe, float cai, float score,
                                  BeamEntry* src_state, int src_var_idx) {
    CodonAggEntry& cell = grid.cells[idx];
    if (!cell.valid) {
        cell.mfe = mfe;
        cell.cai = cai;
        cell.score = score;
        cell.src_state = src_state;
        cell.src_var_idx = src_var_idx;
        cell.valid = true;
        grid.populated_indices.push_back((uint8_t)idx);
        return 0;  // first population, no collision
    }
    // Collision: another variant has already populated this codon-pair cell.
    if (score < cell.score) {
        cell.mfe = mfe;
        cell.cai = cai;
        cell.score = score;
        cell.src_state = src_state;
        cell.src_var_idx = src_var_idx;
        return 2;  // overwrite
    }
    return 1;  // no overwrite
}

// Initialize all 36 cells of a fresh grid to invalid.
static inline void init_grid(CodonAggGrid& grid) {
    for (auto& c : grid.cells) c = CodonAggEntry::invalid();
    grid.populated_indices.clear();
}

// sigma(a, i) = 3*a + i = nucleotide-position index in the [0, nuc_len) range.
// Defined inline here to keep codon_agg decoupled from PositionBeamDP internals.
static inline int agg_sigma(int a, int i) { return 3 * a + i; }

void build_s_codon_agg(DernaBeamMap& tab_s_cell,
                       CodonAggMap& agg_cell,
                       CodonAggTables::PerPosDiag* diag) {
    // Clear any prior aggregation data for this cell. Must be a fresh build:
    // beam pruning may have replaced state-instances, so a stale grid would
    // hold dangling back-pointers.
    agg_cell.clear();

    if (tab_s_cell.empty()) return;

    for (auto& kv : tab_s_cell) {
        BeamEntry& entry = kv.second;
        if (entry.variants.empty()) continue;

        // Aggregation key: left_pos = sigma(entry.a, entry.i).
        // All states with the same left_pos but different boundary nucleotides
        // collapse into the same 36-cell grid here.
        const int left_pos = agg_sigma(entry.a, entry.i);

        auto [it, inserted] = agg_cell.try_emplace(left_pos);
        CodonAggGrid& grid = it->second;
        if (inserted) init_grid(grid);

        if (diag) {
            diag->s_total_states++;
            diag->s_total_variants += (uint64_t)entry.variants.size();
        }

        for (size_t vi = 0; vi < entry.variants.size(); ++vi) {
            const XYVariant& v = entry.variants[vi];
            // Skip sentinel/uninitialized variants defensively.
            if (v.x < 0 || v.x >= NCOD || v.y < 0 || v.y >= NCOD) continue;
            // Skip variants with non-finite score (treated as no candidate).
            if (!std::isfinite(v.score)) continue;

            int idx = v.x * NCOD + v.y;
            int rc = update_agg_cell(grid, idx,
                                     (float)v.mfe, (float)v.cai, (float)v.score,
                                     &entry, (int)vi);
            if (diag) {
                if (rc == 0) {
                    diag->s_cells_populated++;
                } else {
                    // Collision: cell was already populated.
                    diag->s_score_collisions++;
                    if (rc == 2) diag->s_score_overwrites++;
                }
            }
        }
    }

    if (diag) {
        // cells_possible reports the total grid capacity exposed by the agg
        // table (one 36-cell grid per distinct left_pos).
        diag->s_cells_possible += (uint64_t)agg_cell.size() * (uint64_t)CODON_PAIR_GRID;
    }
}

void build_c_codon_agg(DernaBeamMap& tab_c_cell,
                       CodonAggMap& agg_cell,
                       CodonAggTables::PerPosDiag* diag) {
    agg_cell.clear();

    if (tab_c_cell.empty()) return;

    for (auto& kv : tab_c_cell) {
        BeamEntry& entry = kv.second;
        if (entry.variants.empty()) continue;

        const int left_pos = agg_sigma(entry.a, entry.i);

        auto [it, inserted] = agg_cell.try_emplace(left_pos);
        CodonAggGrid& grid = it->second;
        if (inserted) init_grid(grid);

        if (diag) {
            diag->c_total_states++;
            diag->c_total_variants += (uint64_t)entry.variants.size();
        }

        for (size_t vi = 0; vi < entry.variants.size(); ++vi) {
            const XYVariant& v = entry.variants[vi];
            if (v.x < 0 || v.x >= NCOD || v.y < 0 || v.y >= NCOD) continue;
            if (!std::isfinite(v.score)) continue;

            int idx = v.x * NCOD + v.y;
            int rc = update_agg_cell(grid, idx,
                                     (float)v.mfe, (float)v.cai, (float)v.score,
                                     &entry, (int)vi);
            if (diag) {
                if (rc == 0) {
                    diag->c_cells_populated++;
                } else {
                    diag->c_score_collisions++;
                    if (rc == 2) diag->c_score_overwrites++;
                }
            }
        }
    }

    if (diag) {
        diag->c_cells_possible += (uint64_t)agg_cell.size() * (uint64_t)CODON_PAIR_GRID;
    }
}

void agg_diag_write_pos(int pos, const CodonAggTables::PerPosDiag& diag) {
    if (!agg_diag_enabled()) return;
    std::lock_guard<std::mutex> lk(s_diag_mu);
    open_diag_log_if_needed();
    if (!s_diag_out.is_open()) return;

    // S row.
    {
        double fill_rate = (diag.s_cells_possible > 0)
            ? (double)diag.s_cells_populated / (double)diag.s_cells_possible : 0.0;
        double avg_var = (diag.s_total_states > 0)
            ? (double)diag.s_total_variants / (double)diag.s_total_states : 0.0;
        // avg_filled_per_state is over distinct populated cells per state; with
        // 36 cells and one state per (key) we can derive it as cells_populated /
        // total_states (across this pos and all seg_lens).
        double avg_filled = (diag.s_total_states > 0)
            ? (double)diag.s_cells_populated / (double)diag.s_total_states : 0.0;
        s_diag_out << pos << "\tS\t"
                   << diag.s_total_states << "\t"
                   << diag.s_total_variants << "\t"
                   << diag.s_cells_populated << "\t"
                   << diag.s_cells_possible << "\t"
                   << fill_rate << "\t"
                   << diag.s_score_collisions << "\t"
                   << diag.s_score_overwrites << "\t"
                   << avg_var << "\t"
                   << avg_filled << "\n";
    }
    // C row.
    {
        double fill_rate = (diag.c_cells_possible > 0)
            ? (double)diag.c_cells_populated / (double)diag.c_cells_possible : 0.0;
        double avg_var = (diag.c_total_states > 0)
            ? (double)diag.c_total_variants / (double)diag.c_total_states : 0.0;
        double avg_filled = (diag.c_total_states > 0)
            ? (double)diag.c_cells_populated / (double)diag.c_total_states : 0.0;
        s_diag_out << pos << "\tC\t"
                   << diag.c_total_states << "\t"
                   << diag.c_total_variants << "\t"
                   << diag.c_cells_populated << "\t"
                   << diag.c_cells_possible << "\t"
                   << fill_rate << "\t"
                   << diag.c_score_collisions << "\t"
                   << diag.c_score_overwrites << "\t"
                   << avg_var << "\t"
                   << avg_filled << "\n";
    }
    s_diag_out.flush();
}

void agg_diag_close() {
    std::lock_guard<std::mutex> lk(s_diag_mu);
    if (s_diag_out.is_open()) s_diag_out.close();
}
