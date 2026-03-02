import random
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 1. FETCH ACTUAL POU5F1 FROM NCBI
def fetch_POU5F1():
    Entrez.email = "frida.arreytakubetang@gmail.com" # Required by NCBI
    # NM_002701.6 is the official RefSeq for Human POU5F1
    handle = Entrez.efetch(db="nucleotide", id="NM_002701.6", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    # Extract just the CDS (Coding Sequence)
    for feature in record.features:
        if feature.type == "CDS":
            return feature.location.extract(record).seq
    return record.seq

# 2. DESIGN MOCK UTRs (Best Practice)
def design_5utr_simple():
    return "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAG" # Standard high-flow UTR

def design_3utr_simple():
    # We bake in the "Safety Lock" (miR-122 sites) here
    return "CAAACACCAtagCAAACACCAtagCAAACACCAtag" # 3x miR-122 MREs

# 3. SYNONYMOUS SHUFFLE (The Engine)
def synonymous_shuffle(seq):
    # This is a simplified version for the hackathon
    # In a real run, you'd use a full codon frequency table
    seq_list = list(str(seq))
    # Randomly flip a few codons to synonymous versions
    # Placeholder: real synonymous shuffling is more complex
    return "".join(seq_list) 

# 4. SCORING (Your Logic)
def score_construct(cds, utr_5, utr_3):
    mir122_seed = "CAAACACCA"
    mre_count = utr_3.count(mir122_seed)
    selectivity = mre_count / 5.0  

    gc = (cds.count('G') + cds.count('C')) / len(cds)
    stability = 1.0 if 0.30 <= gc <= 0.70 else 0.0

    has_homopolymer = any(c*9 in cds for c in 'ACGT')
    manufacturability = 0.0 if has_homopolymer else 1.0

    fitness = manufacturability * (0.6 * selectivity + 0.4 * stability)
    return fitness, mre_count, gc

# 5. EXECUTION
def run_pipeline():
    print("--- STEP 1: Fetching POU5F1 from NCBI ---")
    raw_cds = fetch_POU5F1()
    
    results = []
    print(f"--- STEP 2: Scoring Variants ---")
    for i in range(10):
        variant_cds = synonymous_shuffle(raw_cds)
        utr_5 = design_5utr_simple()
        utr_3 = design_3utr_simple()
        score, mres, gc = score_construct(variant_cds, utr_5, utr_3)
        results.append((score, variant_cds, utr_5, utr_3, mres, gc))
    
    results.sort(reverse=True, key=lambda x: x[0])
    best = results[0]
    
    # Save as FASTA
    final_seq = best[2] + best[1] + best[3] # 5'UTR + CDS + 3'UTR
    record = SeqRecord(Seq(final_seq), id="SEROVA-001", description="Optimized POU5F1")
    SeqIO.write(record, "POU5F1_optimized.fasta", "fasta")
    
    print(f"--- SUCCESS: Best Fitness {best[0]:.2f} saved to FASTA ---")
    return best

if __name__ == "__main__":
    run_pipeline()