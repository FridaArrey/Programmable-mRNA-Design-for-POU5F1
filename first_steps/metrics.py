def evaluate_biologic_performance(wt_seq, opt_seq, target="Monocyte"):
    metrics = {
        "Stealth (U-Content)": {
            "WT": wt_seq.count('T') / len(wt_seq),
            "OPT": opt_seq.count('T') / len(opt_seq),
            "Goal": "Lower is better (Avoids TLR7/8)"
        },
        "Safety (miR-122 Sites)": {
            "WT": wt_seq.count('CAAACACCA'),
            "OPT": opt_seq.count('CAAACACCA'),
            "Goal": "Higher is better (Liver Silencing)"
        },
        "Stability (GC Content)": {
            "WT": (wt_seq.count('G') + wt_seq.count('C')) / len(wt_seq),
            "OPT": (opt_seq.count('G') + opt_seq.count('C')) / len(opt_seq),
            "Goal": "Target 0.5 - 0.6 for mRNA half-life"
        }
    }
    
    print(f"\n--- BIOLOGICAL EVALUATION: {target} ---")
    print(f"{'Metric':<25} | {'Wild Type':<10} | {'Optimized':<10} | {'Status'}")
    print("-" * 75)
    
    for label, data in metrics.items():
        wt_val = data['WT']
        opt_val = data['OPT']
        # Simple logic to see if we improved
        improvement = "✅ IMPROVED" if (label == "Safety (miR-122 Sites)" and opt_val > wt_val) or \
                                      (label == "Stealth (U-Content)" and opt_val < wt_val) or \
                                      (label == "Stability (GC Content)" and 0.4 < opt_val < 0.7) \
                                      else "⚠️ NO CHANGE"
        
        print(f"{label:<25} | {wt_val:<10.2f} | {opt_val:<10.2f} | {improvement}")