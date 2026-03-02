import matplotlib.pyplot as plt

def plot_final_comparison(wt_metrics, opt_metrics):
    labels = ['Stealth (U-count)', 'Safety (miR-MRE)', 'Stability (GC)']
    wt_vals = [wt_metrics['u'], wt_metrics['mre'], wt_metrics['gc']]
    opt_vals = [opt_metrics['u'], opt_metrics['mre'], opt_metrics['gc']]

    x = range(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar([i - width/2 for i in x], wt_vals, width, label='Wild Type', color='#95a5a6')
    ax.bar([i + width/2 for i in x], opt_vals, width, label='Serova Optimized', color='#2ecc71')

    ax.set_ylabel('Normalized Scores')
    ax.set_title('POU5F1 Optimization Results: Serova Pipeline', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.savefig('serova_final_comparison.png')
    print("--- SUCCESS: serova_final_comparison.png generated ---")
    plt.show()

# Example usage (call this at the end of your run_pipeline)
wt_data = {'u': 0.35, 'mre': 0.0, 'gc': 0.40}
opt_data = {'u': 0.15, 'mre': 0.6, 'gc': 0.55}
plot_final_comparison(wt_data, opt_data)