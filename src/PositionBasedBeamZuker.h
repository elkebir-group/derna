//
// Position-Based Beam Search using BeamZuker's (a,b,i,j,x,y) structure
// Implements LinearCDSfold's position-based DP strategy with BeamZuker's indexing
//

#ifndef DERNA_POSITIONBASEDBEAMZUKER_H
#define DERNA_POSITIONBASEDBEAMZUKER_H

#include "default.h"
#include "utils.h"
#include <tuple>

// Position-based state organized by ending nucleotide position
// Uses BeamZuker's (a,b,i,j,x,y) indexing but organized positionally
class PositionBasedBeamZuker {
    // Position-based tables: O[pos], E[pos], E_single[pos][len], M1[pos], M2[pos], M[pos]
    // Each table at position pos contains states ending at nucleotide position pos
    // States are indexed by (a,b,i,j,x,y) but organized by ending position sigma(b,j)
    
    // N[pos] = open/unpaired stretch ending at position pos (used for hairpin N->C)
    vector<unordered_map<int, BeamEntry>> N;
    
    // O[pos] = best free structure ending at position pos (matches BeamZuker's O)
    vector<unordered_map<int, BeamEntry>> O;
    
    // E[pos] = best closed structure (base pair) ending at position pos (matches BeamZuker's E)
    vector<unordered_map<int, BeamEntry>> E;
    
    // E_single[pos][len] = best single-stranded base pair structure (hairpin/internal loop) 
    // of length len ending at position pos (matches BeamZuker's E for hairpins/internal loops)
    vector<vector<unordered_map<int, BeamEntry>>> E_single;

    // CS[pos]: auxiliary table for right-bulge / right-unpaired extensions after a closed pair.
    // Keyed by the S-segment end-state index (same index space as E_single entries).
    vector<unordered_map<int, BeamEntry>> CS;
    
    // M1[pos] = best multiloop with 1 unpaired region ending at position pos (matches BeamZuker's M)
    vector<unordered_map<int, BeamEntry>> M1;
    
    // M2[pos] = best multiloop with 2 unpaired regions ending at position pos
    vector<unordered_map<int, BeamEntry>> M2;
    
    // M[pos] = best multiloop with 3+ unpaired regions ending at position pos (matches BeamZuker's M)
    vector<unordered_map<int, BeamEntry>> M;
    
    // Backtrace stored in BeamEntry.bt_info for O, E, CS, M1, M2, M (no separate _bt tables)

    // Common data structures
    vector<int> protein, nucle_seq, basepair;
    vector<int> codon_selection;
    vector<bond> bp_bond;
    vector<stack_> sector;
    int n, k;
    int nuc_len;  // Total nucleotide length = 3*n
    double last_lambda_ = 1.0;  // Last lambda used in calculate_position_based (for get_final_mfe/cai display scaling)
    
    // CRITICAL FIX: Removed enable_stage_a_compression member.
    // Lossy Stage A compression has been removed - we now only use score-based pruning
    // to match LCDSfold and pseudocode algorithm semantics.
    
    // Helper: Convert (a,b,i,j) to nucleotide position
    inline int sigma(int a, int i) const {
        return 3*a + i;
    }
    
    // Helper: Convert nucleotide position to (a, i)
    inline pair<int, int> pos_to_ai(int pos) const {
        return {pos / 3, pos % 3};
    }
    
    // Helper: Get index for (a,b,i,j,x,y). Must be collision-free for a,b in [0,n-1], i,j in [0,2], x,y in [0,5].
    // Encoding: 324 = 9*36 (9 for (i,j), 36 for (x,y) with x,y up to 5). Ensures uniqueness for any n.
    inline int index(int a, int b, int i, int j, int x, int y) const {
        const int ij = 3 * i + j;  // 0..8
        const int xy = 6 * x + y;  // 0..35 for x,y in [0,5]
        return a * (n * 324) + b * 324 + ij * 36 + xy;
    }
    
    // Helper: Extract (a,b,i,j,x,y) from index (inverse of index())
    inline tuple<int,int,int,int,int,int> index_to_tuple(int idx) const {
        int rest = idx % 324;
        int y = rest % 6; rest /= 6;
        int x = rest % 6; rest /= 6;
        int ij = rest;  // 0..8
        int j = ij % 3;
        int i = ij / 3;
        idx /= 324;
        int b = idx % n;
        int a = idx / n;
        return {a, b, i, j, x, y};
    }
    
public:
    // Public accessor for index_to_tuple (for comparison)
    tuple<int,int,int,int,int,int> decode_index(int idx) const {
        return index_to_tuple(idx);
    }
    
    // Beam pruning: Keep top-k states at position pos
    void prune_beam(unordered_map<int, BeamEntry>& states, int pos, double lambda);
    
    // QuickSelect for efficient pruning (O(n) average case)
    static unsigned long QuickselectPartition(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper);
    static double QuickSelect(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k);
    
    // Two-pass approach: Compute span by length (Pass 1)
    void compute_span_by_length(int a, int b, int i, int j, int x, int y,
                               int start_pos, int end_pos, int len, double lambda);
    
    // Extend state by one nucleotide (codon extension)
    void extend_by_nucleotide(int j, double lambda);
    
    // Form closed structure from single strand
    void form_closed_structure(int j, double lambda);
    
    // Form single strand (hairpin, internal loop, etc.)
    void form_single_strand(int j, double lambda);
    
    // Form multiloop structures
    void form_multiloop(int j, double lambda);
    
    // Compute hairpin energy. Optional loop-flank (a_loop_start,...,y_loop_end) from N entry for correct loop sequence.
    std::tuple<double,double,double> compute_hairpin_energy(int a, int b, int i, int j, int x, int y, double lambda,
        int a_loop_start = -1, int i_loop_start = 0, int x_loop_start = -1,
        int b_loop_end = -1, int j_loop_end = 0, int y_loop_end = -1);
    
    // Compute internal loop energy
    double compute_internal_loop_energy(int a, int b, int i, int j, int x, int y, 
                                         int c, int d, int i1, int j1, int xh, int xk, double lambda);
    
    // Compute multiloop closing energy
    double compute_multiloop_closing(int a, int b, int i, int j, int x, int y, double lambda);
    
    // Helper functions
    bool rightCodon(int l1, int l2, int x, int y) const;
    // Availability mask for flanking bases (match Zuker ava_nucle_p/m). Returns bitmask of allowed bases.
    int ava_nucleotides_mask(int a, int x, int i, int dir) const;
    double add_hairpin_CAI_2(int a, int b, int x, int y, int a1, int x1, int b1, int y1) const;
    double add_hairpin_CAI_2(int a, int b, int x, int y, int a1, int x1) const;
    double add_hairpin_CAI_2(int a, int b, int x, int y) const;
    double add_CAI(int a, int x) const;
    double add_interior_CAI_2(int a, int c, int x, int na, int x1, int pc, int h1) const;
    
public:
    // Constructor
    // CRITICAL FIX: Removed enable_stage_a_compression parameter.
    // Lossy Stage A compression has been removed entirely - we now only use score-based pruning
    // to match LCDSfold and pseudocode algorithm semantics.
    PositionBasedBeamZuker(int n, vector<int>& protein_seq, int k = 10);
    ~PositionBasedBeamZuker();
    void init_values();
    
    // Cleanup function to fix corrupted E_single before destruction
    // Call this BEFORE the object goes out of scope to prevent crashes
    void cleanup_before_destruction();
    
    /**
     * Main computation: Position-based beam search
     * Computes best structure ending at each position j from left to right
     * 
     * @param lambda weight between MFE and CAI
     * @return best score at final position
     */
    double calculate_position_based(double lambda);
    
    /**
     * Backtrace: Reconstruct sequence and structure
     */
    void traceback(double lambda);
    
    // Getters
    vector<int> get_nucle_seq() const { return nucle_seq; }
    vector<bond> get_bp_bond() const { return bp_bond; }
    vector<int> get_codon_selection() const { return codon_selection; }
    
    // Get RNA sequence and base pairs (same interface as Zuker/BeamZuker)
    void get_bp(string & bp);
    void get_rna(string & rna);
    void get_rna_cai(string & rna);
    void get_rna_X(string & rna);
    
    // Get final score (for comparison)
    double get_final_score() const;
    double get_final_mfe() const;
    double get_final_cai() const;
    
    // Comparison helpers (index_to_tuple is already declared above)
    const vector<unordered_map<int, BeamEntry>>& get_O_table() const { return O; }
    const vector<unordered_map<int, BeamEntry>>& get_E_table() const { return E; }
    const vector<unordered_map<int, BeamEntry>>& get_N_table()     const { return N; }
    const vector<unordered_map<int, BeamEntry>>& get_CS_table()    const { return CS; }
    const vector<unordered_map<int, BeamEntry>>& get_M1_table()    const { return M1; }
    const vector<unordered_map<int, BeamEntry>>& get_M2_table()    const { return M2; }
    const vector<unordered_map<int, BeamEntry>>& get_M_table()     const { return M; }
    const vector<vector<unordered_map<int, BeamEntry>>>& get_E_single_table() const { return E_single; }
};

#endif // DERNA_POSITIONBASEDBEAMZUKER_H
