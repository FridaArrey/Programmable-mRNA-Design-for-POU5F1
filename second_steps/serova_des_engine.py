import math

class mRNAConstruct:
    def __init__(self, five_prime, cds, three_prime):
        self.five_prime = five_prime
        self.cds = cds
        self.three_prime = three_prime

# --- PROXY MODELS FOR SEROVA BIOTOOLS ---

def riboNN_score(sequence):
    """Proxy for Ribosome Load: Rewards high GC at the start of translation."""
    # Higher scores represent open initiation sites in Dendritic Cells.
    return (sequence[:50].count('G') + sequence[:50].count('C')) / 50

def predict_MFE(sequence, context="liver"):
    """Proxy for DIGs: Returns structural stability (Negative kcal/mol)."""
    # More negative = more structured = lower expression in that tissue.
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    return -25.0 * gc_content if context == "liver" else -15.0 * gc_content

def compute_miRNA_repression(three_prime, miRNA="hsa-miR-122-5p"):
    """Proxy for REPRESS: Returns fold-reduction in hepatocytes."""
    # Each miR-122 site provides exponential degradation in the liver.
    mre_count = three_prime.count("CAAACACCA")
    repression = 1.0 / (25.0 * mre_count) if mre_count > 0 else 1.0
    return repression # e.g., 0.02 = 98% degraded

def compute_CIRBP_boost(three_prime):
    """Proxy for DIGs: Returns expression durability in Dendritic Cells."""
    # CIRBP motifs sustain expression specifically in immune cells.
    return 3.5 if "TGTGT" in three_prime else 1.0

# --- THE DES ENGINE ---

def compute_DES(construct):
    """
    Calculates the Differential Expression Score (DES).
    A high DES ensures POU5F1 is ON in DC/Monocytes but OFF in Hepatocytes.
    """
    # 1. ARM 1: 5' UTR — RiboNN MRL ratio (uORF Gating)
    # Target: Dendritic Cell (DC) / Detarget: Hepatocyte
    mrl_target = riboNN_score(construct.five_prime + construct.cds[:100])
    mrl_baseline = 0.4 # Reference wild-type MRL
    utr_component = mrl_target / mrl_baseline

    # 2. ARM 3: CDS — Structural Differential (DIGs)
    # Liver MFE more negative = more structured = lower expression
    mfe_dc = predict_MFE(construct.cds, context="immune")
    mfe_hepatocyte = predict_MFE(construct.cds, context="liver")
    structure_component = mfe_hepatocyte / mfe_dc

    # 3. ARM 2: 3' UTR — miR-122 Silencing + CIRBP Stability
    miRNA_repression_hep = compute_miRNA_repression(construct.three_prime)
    cirbp_stability_dc = compute_CIRBP_boost(construct.three_prime)

    # FINAL DES CALCULATION
    # Numerator (Target DC) / Denominator (Detarget Liver)
    numerator = utr_component * structure_component * cirbp_stability_dc
    denominator = 1.0 * 1.0 * miRNA_repression_hep 

    des_score = numerator / denominator
    
    # 4. τ index (Specificity Proxy)
    tau = 1.0 - (denominator / numerator) if numerator > denominator else 0.0

    return {
        "DES_Score": round(des_score, 2),
        "Tau_Index": round(tau, 4),
        "Hepatocyte_Suppression": f"{round((1-miRNA_repression_hep)*100, 1)}%",
        "DC_Boost": f"{round(cirbp_stability_dc, 1)}x"
    }

if __name__ == "__main__":
    # Test Construct: Optimized POU5F1 with miR-122 and CIRBP motifs
    test_construct = mRNAConstruct(
        five_prime="GGGAAATAAG", 
        cds="ATGGAGACTGCAACCGAGACTGCGCAGCAGCAG", 
        three_prime="CAAACACCATGTGTCAAACACCA"
    )
    
    report = compute_DES(test_construct)
    print("\n--- SEROVA PRECISION EVALUATION ---")
    print(f"DES Score: {report['DES_Score']} (Target/Detarget Ratio)")
    print(f"Tau Specificity: {report['Tau_Index']}")
    print(f"Liver Safety: {report['Hepatocyte_Suppression']} degradation")
    print(f"Immune Potency: {report['DC_Boost']} sustained expression")