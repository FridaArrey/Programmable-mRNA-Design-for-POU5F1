import matplotlib.pyplot as plt

def plot_ga_evolution(fitness_curve):
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, len(fitness_curve) + 1), fitness_curve, 
             marker='o', linestyle='-', color='#3498db', linewidth=2)
    
    plt.title('Evolutionary Progress: POU5F1 Fitness', fontsize=14, fontweight='bold')
    plt.xlabel('Generation')
    plt.ylabel('Fitness Score')
    plt.ylim(0, 1.0) # Show the full scale to show room for growth
    plt.grid(True, linestyle='--', alpha=0.7)

    # Annotate the jump you achieved
    plt.annotate('miR-122 Site Integrated', xy=(1, 0.52), xytext=(1.5, 0.7),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    plt.savefig('serova_ga_evolution.png')
    print("--- SUCCESS: serova_ga_evolution.png generated ---")
    plt.show()

# Example usage (pass your 'curve' variable here)
curve_data = [0.40, 0.52, 0.52, 0.52, 0.52]
plot_ga_evolution(curve_data)