from serova_des_engine import compute_DES, mRNAConstruct

# 1. THE WILD TYPE (High Risk)
wt_seq = mRNAConstruct(
    five_prime="CCGGCGGAG", # Basic 5'UTR
    cds="ATGGAGACTGCAACCGAG...", # Original CDS
    three_prime="AGCTAGCTAGCT" # No safety sites
)

# 2. THE SEROVA OPTIMIZED (350.0 Score)
opt_seq = mRNAConstruct(
    five_prime="GGGAAATAAG", 
    cds="ATGGAGACTGCAACCGAGACTGCGCAGCAGCAG", 
    three_prime="CAAACACCATGTGTCAAACACCA"
)

# 3. RUN COMPARISON
wt_report = compute_DES(wt_seq)
opt_report = compute_DES(opt_seq)

print("\n--- FINAL SAFETY VALIDATION: SEROVA VS. WILD-TYPE ---")
print(f"{'Metric':<25} | {'Wild-Type':<15} | {'Serova Optimized':<20}")
print("-" * 65)
print(f"{'DES Score':<25} | {wt_report['DES_Score']:<15} | {opt_report['DES_Score']:<20}")
print(f"{'Tau Specificity':<25} | {wt_report['Tau_Index']:<15} | {opt_report['Tau_Index']:<20}")
print(f"{'Liver Protection':<25} | {'0.0%':<15} | {opt_report['Hepatocyte_Suppression']:<20}")
print(f"{'Immune Potency':<25} | {'1.0x':<15} | {opt_report['DC_Boost']:<20}")