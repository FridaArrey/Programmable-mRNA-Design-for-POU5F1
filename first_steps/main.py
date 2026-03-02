from engine import mutate_cds, score_sequence, save_output

# Full Human POU5F1 CDS snippet (or use the full 1080bp sequence)
full_pou5f1 = "ATGGAGACTGCAACCGAGACTGCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" # ... etc
target_cell = "Monocyte"

print(f"--- Launching Evolution for {target_cell} ---")

current_best = full_pou5f1
for generation in range(1, 51): # 50 generations looks better in the terminal
    candidate = mutate_cds(current_best)
    if score_sequence(candidate, target_cell) > score_sequence(current_best, target_cell):
        current_best = candidate
        print(f"Gen {generation}: Improvement found! Score: {score_sequence(current_best, target_cell):.4f}")

# THE DELIVERABLE
save_output(current_best, target_cell)