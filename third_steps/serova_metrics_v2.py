"""
serova_metrics_v2.py — Fixed Biological Metrics
================================================
Addresses:
  ⚠️  Stealth (U-content) showing NO CHANGE between wild-type and optimised
  ⚠️  Safety (miR-122 sites) showing NO CHANGE between wild-type and optimised

Root cause of original bug:
  The original optimizer was performing codon substitution but selecting
  synonymous codons at random (uniform), not preferentially selecting low-U
  synonyms. miR-122 sites were never actually inserted into the 3'UTR.

This module provides:
  1. Verified metric computation with before/after diff checking
  2. Active U-minimising codon selection
  3. Active miR-122 site insertion into 3'UTR
  4. Diagnostic reporting that flags NO CHANGE as an error
"""

import json
from typing import Tuple, Dict, List

# ── Codon tables sorted by U-content (ascending) ─────────────────────────────
# For each amino acid, synonyms are sorted: lowest U-count first
LOW_U_PREFERRED = {
    "Ala": ["GCC", "GCA", "GCG", "GCU"],   # GCU has 1 U; GCC has 0
    "Arg": ["CGC", "CGG", "CGA", "CGU", "AGG", "AGA"],
    "Asn": ["AAC", "AAU"],
    "Asp": ["GAC", "GAU"],
    "Cys": ["UGC", "UGU"],                  # both have U; UGC preferred (less U weight)
    "Gln": ["CAG", "CAA"],
    "Glu": ["GAG", "GAA"],
    "Gly": ["GGC", "GGG", "GGA", "GGU"],
    "His": ["CAC", "CAU"],
    "Ile": ["AUC", "AUA", "AUU"],           # AUC < AUA < AUU in U-density contribution
    "Leu": ["CUC", "CUG", "CUA", "CUU", "UUG", "UUA"],
    "Lys": ["AAG", "AAA"],
    "Met": ["AUG"],
    "Phe": ["UUC", "UUU"],
    "Pro": ["CCC", "CCG", "CCA", "CCU"],
    "Ser": ["AGC", "UCC", "UCG", "AGU", "UCA", "UCU"],
    "Thr": ["ACC", "ACG", "ACA", "ACU"],
    "Trp": ["UGG"],
    "Tyr": ["UAC", "UAU"],
    "Val": ["GUC", "GUG", "GUA", "GUU"],
    "Stop": ["UAG", "UAA", "UGA"],
}

# Reverse map: codon → amino acid
CODON_TO_AA = {}
for aa, codons in LOW_U_PREFERRED.items():
    for codon in codons:
        CODON_TO_AA[codon] = aa

MIR122_SITE = "CAAACACCAUUGUCACACUCCA"   # full reverse complement of miR-122
MIR122_SEED = "CAAACACC"                  # 8-mer seed (sufficient for Ago2 loading)


# ── Core metric functions ─────────────────────────────────────────────────────
def compute_u_content(seq: str) -> float:
    rna = seq.upper().replace("T", "U")
    return rna.count("U") / len(rna) if rna else 0.0

def compute_gc_content(seq: str) -> float:
    rna = seq.upper().replace("T", "U")
    return (rna.count("G") + rna.count("C")) / len(rna) if rna else 0.0

def count_mir122_sites(seq: str, full_site: bool = False) -> int:
    """Count miR-122 binding sites (seed match by default)."""
    rna = seq.upper().replace("T", "U")
    query = MIR122_SITE if full_site else MIR122_SEED
    return rna.count(query)

def compute_cai(seq: str) -> float:
    """Codon Adaptation Index (human codon frequencies, Kazusa DB)."""
    FREQ = {
        "GCU": 0.27, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
        "CGU": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.20,
        "AAU": 0.47, "AAC": 0.53, "GAU": 0.46, "GAC": 0.54,
        "UGU": 0.46, "UGC": 0.54, "CAA": 0.27, "CAG": 0.73,
        "GAA": 0.42, "GAG": 0.58, "GGU": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
        "CAU": 0.42, "CAC": 0.58, "AUU": 0.36, "AUC": 0.48, "AUA": 0.16,
        "UUA": 0.08, "UUG": 0.13, "CUU": 0.13, "CUC": 0.20, "CUA": 0.07, "CUG": 0.41,
        "AAA": 0.43, "AAG": 0.57, "AUG": 1.00,
        "UUU": 0.46, "UUC": 0.54, "CCU": 0.29, "CCC": 0.32, "CCA": 0.28, "CCG": 0.11,
        "UCU": 0.15, "UCC": 0.22, "UCA": 0.15, "UCG": 0.06, "AGU": 0.15, "AGC": 0.24,
        "ACU": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11,
        "UGG": 1.00, "UAU": 0.44, "UAC": 0.56,
        "GUU": 0.18, "GUC": 0.24, "GUA": 0.12, "GUG": 0.46,
    }
    rna = seq.upper().replace("T", "U")
    codons = [rna[i:i+3] for i in range(0, len(rna) - 2, 3)]
    vals = [FREQ.get(c, 0.10) for c in codons[1:-1] if len(c) == 3]
    return sum(vals) / len(vals) if vals else 0.0

def full_metrics(seq: str) -> Dict:
    return {
        "u_content": round(compute_u_content(seq), 4),
        "gc_content": round(compute_gc_content(seq), 4),
        "mir122_seed_sites": count_mir122_sites(seq, full_site=False),
        "mir122_full_sites": count_mir122_sites(seq, full_site=True),
        "cai": round(compute_cai(seq), 4),
        "length": len(seq),
    }


# ── Active optimisers ─────────────────────────────────────────────────────────
def minimise_u_content(cds: str) -> str:
    """
    Replace each codon with the synonymous alternative that has the
    lowest U content. This WILL change U-content — verifiably.
    """
    rna = cds.upper().replace("T", "U")
    codons = [rna[i:i+3] for i in range(0, len(rna) - 2, 3) if len(rna[i:i+3]) == 3]
    result = []
    for i, codon in enumerate(codons):
        aa = CODON_TO_AA.get(codon)
        if aa is None or i == 0:  # preserve start codon
            result.append(codon)
            continue
        # Pick synonymous codon with fewest U nucleotides
        synonyms = LOW_U_PREFERRED.get(aa, [codon])
        best = min(synonyms, key=lambda c: (c.count("U"), -_human_freq(c)))
        result.append(best)
    return "".join(result)

def _human_freq(codon: str) -> float:
    FREQ = {"GCC": 0.40, "GAC": 0.54, "CAG": 0.73, "GAG": 0.58, "GGC": 0.34,
            "CAC": 0.58, "AUC": 0.48, "CUG": 0.41, "AAG": 0.57, "AUG": 1.00,
            "UUC": 0.54, "CCC": 0.32, "AGC": 0.24, "ACC": 0.36, "UGG": 1.00,
            "UAC": 0.56, "GUG": 0.46}
    return FREQ.get(codon, 0.15)

def insert_mir122_sites(utr3: str, n_sites: int = 3) -> str:
    """
    Append miR-122 binding site(s) to the 3'UTR.
    Sites are spaced 10nt apart to avoid steric clash (Jain et al. 2018).
    """
    spacer = "GCAUAUGCAU"  # neutral spacer sequence
    sites = []
    for _ in range(n_sites):
        sites.append(MIR122_SITE + spacer)
    insertion = "".join(sites)
    return utr3 + insertion

def build_optimised_mrna(cds: str, utr5: str = "", utr3: str = "", n_mir122: int = 3) -> Dict:
    """
    Full pipeline: minimise U in CDS + insert miR-122 sites in 3'UTR.
    Returns dict with sequences and metric comparisons.
    """
    # Build wild-type full mRNA
    wt_full = utr5 + cds + utr3

    # Step 1: Minimise U in CDS
    opt_cds = minimise_u_content(cds)

    # Step 2: Insert miR-122 sites into 3'UTR
    opt_utr3 = insert_mir122_sites(utr3, n_sites=n_mir122)

    # Build optimised full mRNA
    opt_full = utr5 + opt_cds + opt_utr3

    return {
        "wild_type_cds": cds,
        "optimised_cds": opt_cds,
        "wild_type_utr3": utr3,
        "optimised_utr3": opt_utr3,
        "wild_type_full": wt_full,
        "optimised_full": opt_full,
    }


# ── Validated comparison reporter ─────────────────────────────────────────────
def run_validated_comparison(
    wild_type: str,
    optimised: str,
    label: str = "SEROVA",
) -> None:
    """
    Compare wild-type vs optimised and VERIFY that each optimised metric
    actually changed. Flag any ⚠️ NO CHANGE as an error, not a status.
    """
    wt = full_metrics(wild_type)
    opt = full_metrics(optimised)

    TARGETS = {
        "u_content":          ("↓ lower is better",  lambda a, b: b < a),
        "gc_content":         ("↑ 0.45–0.60 range",  lambda a, b: abs(b - 0.52) < abs(a - 0.52)),
        "mir122_seed_sites":  ("↑ more sites = safer", lambda a, b: b > a),
        "cai":                ("↑ higher is better",  lambda a, b: b > a),
    }

    print(f"\n{'='*72}")
    print(f"  VALIDATED BIOLOGICAL METRICS: {label}")
    print(f"{'='*72}")
    print(f"  {'Metric':<28} | {'Wild-Type':^12} | {'Optimised':^12} | {'Delta':^10} | Status")
    print(f"  {'-'*68}")

    all_pass = True
    results = {}

    for key, (direction, is_better) in TARGETS.items():
        wt_val = wt[key]
        opt_val = opt[key]
        delta = opt_val - wt_val

        if abs(delta) < 1e-6:
            status = "🔴 ERROR: NO CHANGE — optimizer is not modifying this feature"
            all_pass = False
        elif is_better(wt_val, opt_val):
            status = "✅ IMPROVED"
        else:
            status = "⚠️  REGRESSED"

        print(f"  {key:<28} | {wt_val:^12} | {opt_val:^12} | {delta:^+10.4f} | {status}")
        results[key] = {"wt": wt_val, "opt": opt_val, "delta": delta, "status": status}

    print(f"  {'-'*68}")
    overall = "✅ ALL METRICS CHANGED AS EXPECTED" if all_pass else "🔴 SOME METRICS DID NOT CHANGE — CHECK OPTIMIZER"
    print(f"  Overall: {overall}")
    print(f"{'='*72}\n")

    return results


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Demo: Wild-type CDS vs properly optimised version
    WT_CDS = "ATGGAGACUGCAACCGAGACUGCGGAUCGCUUGCAGAACGAAUG"
    WT_UTR3 = "GCAUAUGCAUAUGCAUAUGCAUAUGCAUAUGCAUAUGCAUAUGCA"

    sequences = build_optimised_mrna(WT_CDS, utr3=WT_UTR3, n_mir122=3)

    print(f"\n--- CDS U-content reduction ---")
    wt_u = compute_u_content(sequences["wild_type_cds"])
    opt_u = compute_u_content(sequences["optimised_cds"])
    print(f"  Wild-type CDS U-content:  {wt_u:.4f}")
    print(f"  Optimised CDS U-content:  {opt_u:.4f}")
    assert opt_u < wt_u, "🔴 BUG: U-content did not decrease!"
    print(f"  ✅ U-content decreased by {wt_u - opt_u:.4f}")

    print(f"\n--- miR-122 site insertion ---")
    wt_mir = count_mir122_sites(sequences["wild_type_utr3"])
    opt_mir = count_mir122_sites(sequences["optimised_utr3"])
    print(f"  Wild-type 3'UTR miR-122 seed sites: {wt_mir}")
    print(f"  Optimised 3'UTR miR-122 seed sites: {opt_mir}")
    assert opt_mir > wt_mir, "🔴 BUG: miR-122 sites were not inserted!"
    print(f"  ✅ Added {opt_mir - wt_mir} miR-122 seed site(s)")

    # Full validated comparison
    run_validated_comparison(
        wild_type=sequences["wild_type_full"],
        optimised=sequences["optimised_full"],
        label="SEROVA POU5F1 v3",
    )

    # Save results
    with open("metrics_v2_results.json", "w") as f:
        json.dump({
            "wild_type_metrics": full_metrics(sequences["wild_type_full"]),
            "optimised_metrics": full_metrics(sequences["optimised_full"]),
            "sequences": sequences,
        }, f, indent=2)
    print("Saved: metrics_v2_results.json")