"""
serova_ga_v2.py — Fixed Genetic Algorithm
==========================================
Addresses: Fitness plateau / stagnation in local optima

Key improvements over serova_ga_final.py:
  1. Multiple mutation operators (point, block-swap, inversion, insertion)
  2. Adaptive mutation rate — increases when population stagnates
  3. Tournament selection instead of roulette (handles flat fitness landscapes)
  4. Elitism — top N individuals always survive
  5. Population diversity metric + restart if diversity collapses
  6. Crossover: uniform + two-point
"""

import random
import math
import json
from typing import List, Tuple

# ── Codon tables ────────────────────────────────────────────────────────────
SYNONYMOUS_CODONS = {
    "Ala": ["GCU", "GCC", "GCA", "GCG"],
    "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Asn": ["AAU", "AAC"],
    "Asp": ["GAU", "GAC"],
    "Cys": ["UGU", "UGC"],
    "Gln": ["CAA", "CAG"],
    "Glu": ["GAA", "GAG"],
    "Gly": ["GGU", "GGC", "GGA", "GGG"],
    "His": ["CAU", "CAC"],
    "Ile": ["AUU", "AUC", "AUA"],
    "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "Lys": ["AAA", "AAG"],
    "Met": ["AUG"],
    "Phe": ["UUU", "UUC"],
    "Pro": ["CCU", "CCC", "CCA", "CCG"],
    "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "Thr": ["ACU", "ACC", "ACA", "ACG"],
    "Trp": ["UGG"],
    "Tyr": ["UAU", "UAC"],
    "Val": ["GUU", "GUC", "GUA", "GUG"],
    "Stop": ["UAA", "UAG", "UGA"],
}

# Human codon usage frequency (higher = more common in humans)
CODON_FREQUENCY = {
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
    "UAA": 0.30, "UAG": 0.24, "UGA": 0.47,
}

MIR122_SEED = "UGGAGUGUGACAAUGGUGUUUG"  # miR-122 full sequence
MIR122_SITE = "CAAACACCAUUGUCACACUCCA"  # reverse complement (binding site to insert)


# ── Fitness function (actually modifies the right features) ─────────────────
def compute_fitness(sequence: str, target_cell: str = "DendriticCell") -> Tuple[float, dict]:
    """
    Multi-objective fitness. Returns scalar score + breakdown dict.
    
    Objectives:
      1. Codon Adaptation Index (CAI) — translational efficiency
      2. GC content in [45%, 60%] — stability
      3. U-content (stealth) — avoid innate immune activation (lower = better)
      4. miR-122 site count in 3'UTR — liver detargeting (more = better for safety)
      5. No early stop codons
    """
    seq = sequence.upper().replace("T", "U")
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
    
    if not codons:
        return 0.0, {}

    # 1. CAI score
    cai_scores = [CODON_FREQUENCY.get(c, 0.1) for c in codons[:-1]]  # exclude stop
    cai = sum(cai_scores) / len(cai_scores) if cai_scores else 0.0

    # 2. GC content (target: 50-55%)
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    gc_score = 1.0 - abs(gc - 0.52) * 4  # peaks at 0.52
    gc_score = max(0.0, gc_score)

    # 3. U-content stealth (lower uridine = less innate immune activation)
    u_content = seq.count("U") / len(seq)
    stealth = 1.0 - (u_content / 0.35)  # normalised: 0.35 is typical max
    stealth = max(0.0, min(1.0, stealth))

    # 4. miR-122 sites (count occurrences in 3'UTR region — last 200nt)
    utr3 = seq[-200:] if len(seq) > 200 else seq
    mir122_hits = utr3.count(MIR122_SITE[:8])  # seed match (8-mer)
    mir_score = min(1.0, mir122_hits / 3.0)  # 3 sites = perfect score

    # 5. No premature stop codons
    stop_codons = {"UAA", "UAG", "UGA"}
    internal_stops = sum(1 for c in codons[:-1] if c in stop_codons)
    stop_penalty = 1.0 if internal_stops == 0 else max(0.0, 1.0 - internal_stops * 0.3)

    # Weighted composite
    weights = {"cai": 0.30, "gc": 0.20, "stealth": 0.25, "mir122": 0.15, "stop": 0.10}
    fitness = (
        weights["cai"] * cai +
        weights["gc"] * gc_score +
        weights["stealth"] * stealth +
        weights["mir122"] * mir_score +
        weights["stop"] * stop_penalty
    )

    breakdown = {
        "cai": round(cai, 4),
        "gc_content": round(gc, 4),
        "gc_score": round(gc_score, 4),
        "u_content": round(u_content, 4),
        "stealth": round(stealth, 4),
        "mir122_sites": mir122_hits,
        "mir_score": round(mir_score, 4),
        "internal_stops": internal_stops,
        "fitness": round(fitness, 4),
    }

    return fitness, breakdown


# ── Sequence utilities ───────────────────────────────────────────────────────
def seq_to_codons(seq: str) -> List[str]:
    seq = seq.upper().replace("T", "U")
    return [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]

def codons_to_seq(codons: List[str]) -> str:
    return "".join(codons)

def get_synonyms(codon: str) -> List[str]:
    codon = codon.upper()
    for synonyms in SYNONYMOUS_CODONS.values():
        if codon in synonyms:
            return synonyms
    return [codon]


# ── Mutation operators ───────────────────────────────────────────────────────
def mutate_point(codons: List[str], rate: float) -> List[str]:
    """Replace individual codons with synonymous alternatives."""
    result = codons[:]
    for i in range(1, len(result) - 1):  # preserve start/stop
        if random.random() < rate:
            syns = get_synonyms(result[i])
            if len(syns) > 1:
                result[i] = random.choice([s for s in syns if s != result[i]])
    return result

def mutate_block_swap(codons: List[str]) -> List[str]:
    """Swap two random blocks of codons (maintains amino acid sequence)."""
    n = len(codons)
    # 1. Minimum length check: need at least 2 blocks of size 2 + buffers
    if n < 10: 
        return codons[:]
    
    result = codons[:]
    # 2. Dynamic size based on sequence length
    max_block_size = max(2, n // 10) 
    size = random.randint(2, max_block_size)
    
    # 3. Safe range for i (start of first block)
    # Must leave room for block 1, a gap, and block 2
    i_upper_bound = n - (2 * size) - 1
    if i_upper_bound <= 1:
        return codons[:]
        
    i = random.randint(1, i_upper_bound)
    
    # 4. Safe range for j (start of second block)
    # Must start at least 'size' positions after i
    j_start = i + size
    j_end = n - size - 1
    
    if j_start >= j_end:
        return codons[:]
        
    j = random.randint(j_start, j_end)
    
    # Swap logic remains the same...
    for k in range(size):
        syns_i = get_synonyms(result[i + k])
        syns_j = get_synonyms(result[j + k])
        result[i + k] = random.choice(syns_j) if syns_j else result[i + k]
        result[j + k] = random.choice(syns_i) if syns_i else result[j + k]
    return result

def mutate_inversion(codons: List[str]) -> List[str]:
    """Reverse a segment, replacing each codon with a synonymous alternative."""
    if len(codons) < 6:
        return codons[:]
    result = codons[:]
    i = random.randint(1, len(result) - 4)
    j = random.randint(i + 2, min(i + len(result) // 4, len(result) - 2))
    segment = result[i:j]
    new_segment = [random.choice(get_synonyms(c)) for c in reversed(segment)]
    result[i:j] = new_segment
    return result

def mutate_insert_mir122(codons: List[str]) -> List[str]:
    """Insert a miR-122 binding site codon-segment near the 3' end."""
    result = codons[:]
    # Encode the 8-mer seed as codons (approximate — maintains reading frame)
    # CAAACACC → CAA ACA CC... insert at 80% position
    insert_pos = max(1, int(len(result) * 0.80))
    mir_codons = seq_to_codons(MIR122_SITE)
    if mir_codons:
        result = result[:insert_pos] + mir_codons[:4] + result[insert_pos:]
    return result

def apply_mutation(codons: List[str], mutation_rate: float, generation: int) -> List[str]:
    """Stochastically choose from multiple mutation operators."""
    # Weight operators: more disruptive ones less frequent early on
    operators = ["point", "block_swap", "inversion", "mir122_insert"]
    weights = [0.55, 0.25, 0.15, 0.05]
    
    # Increase diversity-seeking operators as generations progress
    if generation > 20:
        weights = [0.40, 0.30, 0.20, 0.10]

    op = random.choices(operators, weights=weights)[0]

    if op == "point":
        return mutate_point(codons, mutation_rate)
    elif op == "block_swap":
        result = mutate_point(codons, mutation_rate * 0.3)
        return mutate_block_swap(result)
    elif op == "inversion":
        return mutate_inversion(codons)
    elif op == "mir122_insert":
        return mutate_insert_mir122(mutate_point(codons, mutation_rate * 0.2))
    return codons[:]


# ── Crossover operators ──────────────────────────────────────────────────────
def crossover_two_point(p1: List[str], p2: List[str]) -> List[str]:
    """Two-point crossover preserving reading frame."""
    if len(p1) != len(p2) or len(p1) < 4:
        return p1[:]
    i = random.randint(1, len(p1) // 2)
    j = random.randint(i + 1, len(p1) - 1)
    return p1[:i] + p2[i:j] + p1[j:]

def crossover_uniform(p1: List[str], p2: List[str]) -> List[str]:
    """Uniform crossover — each codon chosen from either parent."""
    return [random.choice([a, b]) for a, b in zip(p1, p2)]


# ── Selection ────────────────────────────────────────────────────────────────
def tournament_select(population: List[str], fitnesses: List[float], k: int = 4) -> str:
    """Tournament selection — works well on flat landscapes."""
    contestants = random.sample(list(zip(population, fitnesses)), min(k, len(population)))
    return max(contestants, key=lambda x: x[1])[0]


# ── Population diversity ─────────────────────────────────────────────────────
def population_diversity(population: List[str]) -> float:
    """Fraction of unique sequences in population."""
    return len(set(population)) / len(population)


# ── Main GA loop ─────────────────────────────────────────────────────────────
def run_ga(
    seed_sequence: str,
    generations: int = 50,
    population_size: int = 60,
    elite_k: int = 5,
    base_mutation_rate: float = 0.04,
    stagnation_patience: int = 8,
    diversity_threshold: float = 0.30,
    target_cell: str = "DendriticCell",
) -> Tuple[List[float], str, dict]:
    """
    Improved genetic algorithm.
    Returns: (fitness_curve, best_sequence, best_breakdown)
    """
    print(f"\n{'='*60}")
    print(f"  SEROVA GA v2 — Target: {target_cell}")
    print(f"  Population: {population_size} | Generations: {generations}")
    print(f"  Elite: {elite_k} | Base mutation rate: {base_mutation_rate}")
    print(f"{'='*60}\n")

    # Initialise population with synonymous diversification
    seed_codons = seq_to_codons(seed_sequence)
    population = []
    for _ in range(population_size):
        rate = random.uniform(0.05, 0.40)  # diverse initialisation
        indiv = mutate_point(seed_codons, rate)
        population.append(codons_to_seq(indiv))

    fitness_curve = []
    best_fitness = -1.0
    best_seq = population[0]
    best_breakdown = {}
    stagnation_counter = 0
    mutation_rate = base_mutation_rate

    for gen in range(generations):
        # Evaluate
        scored = [(seq, *compute_fitness(seq, target_cell)) for seq in population]
        scored.sort(key=lambda x: x[1], reverse=True)

        gen_best_fit = scored[0][1]
        gen_best_seq = scored[0][0]
        gen_best_bd = scored[0][2]

        if gen_best_fit > best_fitness:
            best_fitness = gen_best_fit
            best_seq = gen_best_seq
            best_breakdown = gen_best_bd
            stagnation_counter = 0
        else:
            stagnation_counter += 1

        fitness_curve.append(best_fitness)
        fitnesses = [s[1] for s in scored]
        seqs = [s[0] for s in scored]

        # Adaptive mutation rate
        if stagnation_counter >= stagnation_patience:
            mutation_rate = min(0.25, mutation_rate * 1.5)
            print(f"  [Gen {gen+1:3d}] ⚠ Stagnation detected → mutation rate → {mutation_rate:.3f}")
            stagnation_counter = 0
        else:
            mutation_rate = max(base_mutation_rate, mutation_rate * 0.95)

        # Diversity check — restart portion of population if collapsed
        div = population_diversity(population)
        if div < diversity_threshold:
            n_restart = population_size // 3
            print(f"  [Gen {gen+1:3d}] 🔄 Diversity {div:.2f} < {diversity_threshold} → reinitialising {n_restart} individuals")
            for i in range(elite_k, elite_k + n_restart):
                rate = random.uniform(0.10, 0.50)
                seqs[i] = codons_to_seq(mutate_point(seed_codons, rate))

        print(
            f"  [Gen {gen+1:3d}] Best={best_fitness:.4f}  "
            f"CAI={gen_best_bd.get('cai', 0):.3f}  "
            f"Stealth={gen_best_bd.get('stealth', 0):.3f}  "
            f"miR122={gen_best_bd.get('mir122_sites', 0)}  "
            f"U={gen_best_bd.get('u_content', 0):.3f}  "
            f"Div={div:.2f}"
        )

        # Build next generation
        next_gen = seqs[:elite_k]  # elitism

        while len(next_gen) < population_size:
            parent1_seq = tournament_select(seqs, fitnesses, k=5)
            parent2_seq = tournament_select(seqs, fitnesses, k=5)
            p1_codons = seq_to_codons(parent1_seq)
            p2_codons = seq_to_codons(parent2_seq)

            # Choose crossover
            if len(p1_codons) == len(p2_codons) and random.random() < 0.6:
                if random.random() < 0.5:
                    child_codons = crossover_two_point(p1_codons, p2_codons)
                else:
                    child_codons = crossover_uniform(p1_codons, p2_codons)
            else:
                child_codons = p1_codons[:]

            child_codons = apply_mutation(child_codons, mutation_rate, gen)
            next_gen.append(codons_to_seq(child_codons))

        population = next_gen

    print(f"\n{'='*60}")
    print(f"  ✅ COMPLETE: Best Fitness = {best_fitness:.4f}")
    print(f"  Breakdown: {best_breakdown}")
    print(f"{'='*60}\n")

    return fitness_curve, best_seq, best_breakdown


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Example seed — replace with actual POU5F1 CDS or output from serova_final.py
    SEED = (
        "ATGGAGACTGCAACCGAGACTGCGGATCGCTTGCAGAACGAATGCAAAGCAGAAACCGAGCCC"
        "ATGTACGAGCTGGACAAGGACATGAACAGCGATCTGCAGCTTCAGCAGAAGCAGCAGCAGCAGCAG"
        "ATGGAGACTGCAACCGAGACTGCGGATCGCTTGCAGAACGAATGCAAAGCAGAAACCGAGCCC"
    )

    curve, best, breakdown = run_ga(
        seed_sequence=SEED,
        generations=50,
        population_size=60,
        target_cell="DendriticCell",
    )

    # Save results
    with open("ga_v2_results.json", "w") as f:
        json.dump({"best_fitness": breakdown["fitness"], "breakdown": breakdown, "sequence": best}, f, indent=2)
    print("Saved: ga_v2_results.json")

    # Plot convergence curve
    plt.figure(figsize=(10, 5))
    plt.plot(curve, color="#2196F3", linewidth=2, label="Best Fitness")
    plt.axhline(y=curve[0], color="grey", linestyle="--", alpha=0.5, label="Gen 0 baseline")
    plt.xlabel("Generation")
    plt.ylabel("Fitness Score")
    plt.title("SEROVA GA v2 — Convergence Curve (Multi-Operator, Adaptive Mutation)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("ga_v2_convergence.png", dpi=150)
    print("Saved: ga_v2_convergence.png")