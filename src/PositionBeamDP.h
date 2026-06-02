//
// Position-based beam DP filling AllTablesDerna using BeamEntry (from default.h) and (a,b,i,j,x,y) indexing; energy from Zuker.h.
//

#ifndef DERNA_POSITIONBEAMDP_H
#define DERNA_POSITIONBEAMDP_H

#include "default.h"
#include <cstddef>
#include <functional>                                                                                                                                                           
#include <vector>                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
#include <unordered_map>
#include "third_party/unordered_dense.h"

using DernaBeamMap = ankerl::unordered_dense::map<int, BeamEntry>;

#ifndef HAIRPIN_GAP
#define HAIRPIN_GAP 3
#endif
// Cap on single-stranded segment length inside interior loops, hairpin loops,
// and multi-loops. The energy model permits up to MAXLOOP=30 (Turner standard,
// used by LinearFold and CONTRAfold). Lower values prune more aggressively at
// the cost of dropping optimal structures whose segments exceed the cap.
#ifndef SINGLE_MAX_LEN
#define SINGLE_MAX_LEN 30
#endif

// DP tables: indexed by nucleotide position; each holds states keyed by index(a,b,i,j,x,y).
// bestCS uses an int key from cs_get_index: encodes (cs_inner_left, nuc_inner_L, nuc_outer_R, nuc_single_start, single_len).
struct AllTablesDerna {
    std::vector<DernaBeamMap> bestN;
    std::vector<std::vector<DernaBeamMap>> bestS;
    std::vector<DernaBeamMap> bestF;
    std::vector<DernaBeamMap> bestC;
    std::vector<DernaBeamMap> bestCS;
    std::vector<DernaBeamMap> bestSCLeft; // SC_left = S_left + C composites (for left-bulge closing)
    std::vector<DernaBeamMap> bestM1;
    std::vector<DernaBeamMap> bestM2;
    std::vector<DernaBeamMap> bestMulti;
};

// Backtrace manner constants (stored in BeamEntry.backtrace_type).
enum DernaManner {
    MANNER_NONEtoN = 1,
    MANNER_NONEtoF = 2,
    MANNER_NONEtoS = 3,
    MANNER_N_EtoN = 4,
    MANNER_NtoC = 5,
    MANNER_S_EtoS = 6,
    MANNER_CStoC = 7,
    MANNER_S_CStoC = 8,
    MANNER_C_StoCS = 9,
    MANNER_CtoC = 10,
    MANNER_S_CtoC = 11,
    MANNER_Multi_EtoMulti = 12,
    MANNER_MultitoC = 13,
    MANNER_CtoM1 = 14,
    MANNER_M1_CtoM2 = 15,
    MANNER_M1_EtoM1 = 16,
    MANNER_M2toM1 = 17,
    MANNER_M2toMulti = 18,
    MANNER_S_M2toMulti = 19,
    MANNER_F_EtoF = 20,
    MANNER_CtoF = 21,
    MANNER_F_CtoF = 22,
    MANNER_C_StoSCLeft = 23, // C + S_left → SCLeft composite (for left-bulge)
    MANNER_SCLefttoC = 24,   // SCLeft → C (close left-bulge with outer pair)
    MANNER_C_StoC = 25,      // Direct right-bulge: C[q] + S_right[j-1] → C[j] (Block 1, no CS composite)
    MANNER_S_C_StoC = 26,    // Direct internal loop: S_left + C[q] + S_right[j-1] → C[j] (Block 1)
    MANNER_SpecialHP = 27    // Seeded special-hairpin closure (LCDSfold initialize_Special_HP_LD)
};

/**
 * Fill AllTablesDerna by position-based beam search (Zuker stacking/bulge energy).
 * Maximizes score = (1-lambda)*CAI - lambda*MFE (stored as mfe + cai).
 *
 * @param n number of amino acids
 * @param protein amino acid sequence (indices 0..19)
 * @param lambda weight for structure (1-lambda for CAI)
 * @param beamsize beam width for pruning
 * @param tables tables to fill; resized internally
 */
void fill_position_beam_tables(int n, std::vector<int>& protein, double lambda,
                               int beamsize, AllTablesDerna& tables);

/**
 * Traceback from filled AllTablesDerna to reconstruct RNA sequence and structure.
 * Fills nucle_seq (nucleotide codes 0..3), codon_selection (codon index per aa), and bp_bond (base pairs).
 *
 * @param tables tables filled by fill_position_beam_tables
 * @param n number of amino acids
 * @param protein amino acid sequence (indices 0..19)
 * @param nucle_seq output: length 3*n, nucleotide at each position (0=A,1=C,2=G,3=U)
 * @param codon_selection output: length n, codon index chosen per amino acid
 * @param bp_bond output: list of (i,j) base pairs
 */
void traceback_position_beam_tables(const AllTablesDerna& tables, int n, const std::vector<int>& protein,
                                    std::vector<int>& nucle_seq, std::vector<int>& codon_selection,
                                    std::vector<bond>& bp_bond);

#endif // DERNA_POSITIONBEAMDP_H
