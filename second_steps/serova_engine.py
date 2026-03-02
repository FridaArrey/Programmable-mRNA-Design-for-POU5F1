import re

def check_GC_windows(seq, window_size=50):
    """Hard Gate: Ensures GC content stays between 30-70% in every 50bp window."""
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        gc = (window.count('G') + window.count('C')) / window_size
        if gc < 0.3 or gc > 0.7: return False
    return True

def check_homopolymers(seq, limit=8):
    """Hard Gate: Prevents synthesis 'stuttering' by banning runs > 8nt."""
    for base in "ACGT":
        if base * (limit + 1) in seq: return False
    return True

def compute_DES(construct):
    """Proxy for miR-122 silencing. Normalizes site count 0-1."""
    sites = construct.count("CAAACACCA")
    return min(sites / 3.0, 1.0) # Max score at 3 sites

def compute_tau(construct):
    """Proxy for PARADE: Rewards U-depletion (Stealth/Specificity)."""
    u_count = construct.count('T') / len(construct)
    return max(0, 1 - (u_count / 0.3)) # Penalty if U > 30%

def helix_mRNA_score(construct):
    """Proxy for Stability: Rewards GC content closer to 55%."""
    gc = (construct.count('G') + construct.count('C')) / len(construct)
    return 1.0 - abs(0.55 - gc)

def riboNN_score(construct):
    """Proxy for MRL: Rewards Codon optimality."""
    # Simplified: High GC in the first 50bp usually improves ribosome loading
    head_gc = (construct[:50].count('G') + construct[:50].count('C')) / 50
    return head_gc

def fitness(construct_seq, weights={'selectivity': 0.6, 'stability': 0.4}):
    """The Final Triple-Lock Fitness Evaluator."""
    # 1. THE BOUNCER
    if not (check_GC_windows(construct_seq) and check_homopolymers(construct_seq)):
        return 0.0

    # 2. THE SAFETY LOCK (60%)
    selectivity_score = (compute_DES(construct_seq) * 0.6) + (compute_tau(construct_seq) * 0.4)

    # 3. THE POTENCY LOCK (40%)
    stability_score = (helix_mRNA_score(construct_seq) * 0.4) + (riboNN_score(construct_seq) * 0.4) + 0.2

    return (weights['selectivity'] * selectivity_score) + (weights['stability'] * stability_score)