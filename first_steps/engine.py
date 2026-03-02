import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def save_output(winning_sequence, target_name):
    # No Alphabet needed in modern Biopython
    record = SeqRecord(
        Seq(winning_sequence),
        id="POU5F1-OPT",
        name="Oct4_Optimized",
        description=f"POU5F1 optimized for {target_name} vs Liver",
        annotations={"molecule_type": "mRNA"} # This replaces Alphabet
    )
    filename = f"POU5F1_{target_name}_winning.gb"
    with open(filename, "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")
    print(f"\n--- SUCCESS: {filename} generated ---")

def get_synonymous_codons():
    # Simplified map for the hackathon
    return {
        'GCA': ['GCT', 'GCC', 'GCG'], 'AGA': ['AGG', 'CGT', 'CGC'],
        # Add more from a standard codon table...
    }

def mutate_cds(sequence):
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    pos = random.randint(0, len(codons) - 1)
    table = get_synonymous_codons()
    
    if codons[pos] in table:
        codons[pos] = random.choice(table[codons[pos]])
    return "".join(codons)

def score_sequence(cds, target="DC"):
    # LIVER DETARGETING (The Negative Arm)
    mir122_site = "CAAACACCA"
    liver_penalty = cds.count(mir122_site) * 0.5
    
    # TARGET SPECIFICITY (The Positive Arm)
    if target == "HSC":
        # HSC: Focus on GC balance for stability in quiescence
        gc_content = (cds.count('G') + cds.count('C')) / len(cds)
        target_score = 1.0 - abs(0.5 - gc_content) 
    elif target == "Monocyte":
        # Monocyte: Focus on Uridine depletion (Stealth)
        u_content = cds.count('T') / len(cds)
        target_score = 1.0 - u_content
    else: # Default: DC
        # DC: High translation (CAI proxy)
        target_score = 0.8 # Placeholder for RiboNN
        
    return target_score - liver_penalty