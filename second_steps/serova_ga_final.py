import random
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from serova_engine import fitness

# Standard Codon Table for Synonymous Mutations
CODON_TABLE = {
    'ATT': ['ATT', 'ATC', 'ATA'], 'CTT': ['CTT', 'CTC', 'CTA', 'CTG'],
    'GTT': ['GTT', 'GTC', 'GTA', 'GTG'], 'TTT': ['TTT', 'TTC'],
    'ACT': ['ACT', 'ACC', 'ACA', 'ACG'], 'GCT': ['GCT', 'GCC', 'GCA', 'GCG'],
    # ... (Logic can be expanded with full table)
}

def mutate(sequence):
    """Performs a single synonymous codon swap to maintain protein identity."""
    sequence = str(sequence)
    pos = random.randint(0, (len(sequence)//3) - 1) * 3
    codon = sequence[pos:pos+3]
    synonyms = CODON_TABLE.get(codon, [codon])
    new_codon = random.choice(synonyms)
    return sequence[:pos] + new_codon + sequence[pos+3:]

def score_construct(cds, utr_5, utr_3):
    # SELECTIVITY (miR-122 Silencing)
    mir122_seed = "CAAACACCA"
    mre_count = utr_3.count(mir122_seed)
    selectivity = mre_count / 5.0  

    # STABILITY (GC Content Gate)
    gc = (cds.count('G') + cds.count('C')) / len(cds)
    stability = 1.0 if 0.30 <= gc <= 0.70 else 0.0

    # MANUFACTURABILITY (Homopolymer Check)
    has_homopolymer = any(c*9 in cds for c in 'ACGT')
    manufacturability = 0.0 if has_homopolymer else 1.0

    # FINAL FITNESS
    fitness = manufacturability * (0.6 * selectivity + 0.4 * stability)
    return fitness, mre_count, gc

def run_evolutionary_pipeline(n_generations=5):
    # 1. Fetch from NCBI
    Entrez.email = "frida.arreytakubetang@gmail.com" # Required by NCBI
    handle = Entrez.efetch(db="nucleotide", id="NM_002701.6", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    wt_cds = ""
    for f in record.features:
        if f.type == "CDS": wt_cds = str(f.location.extract(record).seq)
    
    # SEEDING: Manually add a miR-122 site to the 3' UTR placeholder
    # This simulates a "best-practice" starting construct
    optimized_utr3 = "CAAACACCA" + "A" * 20
    
    # 2. Initial Population (10 variants of WT)
    population = [mutate(wt_cds) for _ in range(10)]
    fitness_curve = []

    print(f"--- Starting Evolution for 5 Generations ---")

    for gen in range(n_generations):
        # Score and Sort
        # (Assuming your score_construct function is imported/defined)
        scored = []
        for cand in population:
            # We use your Triple-Lock logic: Selectivity, Stability, Manufacturability
            fit_val, mres, gc = score_construct(cand, "GGGAAATAAG", optimized_utr3)
            scored.append((fit_val, cand))
        
        scored.sort(reverse=True, key=lambda x: x[0])
        best_score = scored[0][0]
        fitness_curve.append(best_score)
        
        print(f"Generation {gen+1}: Best Fitness = {best_score:.4f}")

        # Selection: Top 3 survive
        survivors = [c for _, c in scored[:3]]
        
        # Repopulate: Mutate survivors to get back to 10
        new_population = survivors[:]
        while len(new_population) < 10:
            parent = random.choice(survivors)
            new_population.append(mutate(parent))
        population = new_population

    return fitness_curve, scored[0][1]

# Visualization call
curve, best_seq = run_evolutionary_pipeline()