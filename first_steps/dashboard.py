from engine import mutate_cds, score_sequence, save_output
from metrics import evaluate_biologic_performance

# 1. SETUP
target = "Monocyte"
wt_sequence = "ATGGAGACTGCAACCGAGACTGCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" # Placeholder
current_best = wt_sequence

print(f"--- RUNNING FINAL VALIDATION FOR {target} ---")

# 2. RUN MINI-GA (50 Generations)
for gen in range(1, 51):
    candidate = mutate_cds(current_best)
    if score_sequence(candidate, target) > score_sequence(current_best, target):
        current_best = candidate

# 3. OUTPUT THE "BIOLOGIC REPORT"
evaluate_biologic_performance(wt_sequence, current_best, target)

# 4. SAVE THE FINAL FILE
save_output(current_best, target)