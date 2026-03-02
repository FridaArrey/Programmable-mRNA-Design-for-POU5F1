import re

# --- BIOLOGICAL GATE LOGIC ---

def check_GC_windows(seq, window_size=50):
    """Hard Gate: Ensures GC content stays between 30-70% in every 50bp window."""
    # This ensures the mRNA doesn't have extreme local stability issues
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        gc = (window.count('G') + window.count('C')) / window_size
        if gc < 0.3 or gc > 0.7: return False
    return True

def check_homopolymers(seq, limit=8):
    """Hard Gate: Prevents synthesis 'stuttering' by banning runs > 8nt."""
    # Long runs of a single nucleotide (e.g., AAAAAAAAA) are difficult to synthesize
    for base in "ACGT":
        if base * (limit + 1) in seq: return False
    return True

# --- COMPONENT SCORING (PROXIES) ---

def compute_DES(construct):
    """Selectivity Proxy: Counts miR-122 sites for liver detargeting."""
    # miR-122 binding sites act as 'molecular scissors' in the liver
    mir122_seed = "CAAACACCA"
    sites = construct.count(mir122_seed)
    return min(sites / 3.0, 1.0) # Normalized: 1.0 score at 3+ sites

def compute_tau(construct):
    """Selectivity Proxy: Rewards Uridine-depletion for immune stealth."""
    # Lower U-content helps the mRNA bypass TLR7/8 immune sensors in monocytes
    u_count = construct.count('T') / len(construct)
    return max(0, 1 - (u_count / 0.3)) # Ideal is < 30% U

def helix_mRNA_score(construct):
    """Potency Proxy: Measures stability via GC-content normalization."""
    # Target GC content of ~55% is ideal for mRNA half-life
    gc = (construct.count('G') + construct.count('C')) / len(construct)
    return 1.0 - abs(0.55 - gc)

def riboNN_score(construct):
    """Potency Proxy: Estimates Ribosome Loading (MRL)."""
    # High GC at the start (5' end) correlates with better translation initiation
    head_gc = (construct[:50].count('G') + construct[:50].count('C')) / 50
    return head_gc

# --- MAIN FITNESS ENGINE ---

def fitness(construct_seq, weights={'selectivity': 0.6, 'stability': 0.4}):
    """
    The Master Evaluator:
    1. Manufacturability Hard Gate (0 or 1)
    2. Selectivity (60% Weighting)
    3. Stability/Potency (40% Weighting)
    """
    # 1. THE BOUNCER
    if not (check_GC_windows(construct_seq) and check_homopolymers(construct_seq)):
        print(">>> FAILED: Manufacturability Hard Gate (GC Windows or Homopolymers)")
        return 0.0

    # 2. THE SAFETY LOCK (60% Weight)
    # Balanced between miR-122 sites (60% of component) and U-depletion (40% of component)
    selectivity_score = (compute_DES(construct_seq) * 0.6) + (compute_tau(construct_seq) * 0.4)

    # 3. THE POTENCY LOCK (40% Weight)
    # Balanced between half-life and ribosome loading
    stability_score = (helix_mRNA_score(construct_seq) * 0.4) + (riboNN_score(construct_seq) * 0.4) + 0.2

    final_score = (weights['selectivity'] * selectivity_score) + (weights['stability'] * stability_score)
    
    return round(final_score, 4)

# --- STANDALONE TESTER ---
if __name__ == "__main__":
    # Example: A test sequence with a miR-122 site
    test_seq = "GGGAAATAAG" + "ATGGAGACTGCAACCGAGACTGCGCAGCAGCAG" + "CAAACACCA"
    print(f"Testing Sequence Analysis...")
    print(f"Final Fitness Score: {fitness(test_seq)}")