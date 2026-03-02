import matplotlib.pyplot as plt

def generate_fitness_plot(generations, scores, target_name):
    plt.figure(figsize=(10, 6))
    plt.style.use('seaborn-v0_8-darkgrid') # Gives it a "tech" look
    
    # Plotting the data
    plt.plot(generations, scores, marker='o', linestyle='-', color='#007acc', linewidth=2, markersize=6)
    
    # Adding titles and labels
    plt.title(f'Evolutionary Optimization: POU5F1 for {target_name}', fontsize=16, fontweight='bold')
    plt.xlabel('Generation', fontsize=12)
    plt.ylabel('Fitness Score (Specificity + Stability)', fontsize=12)
    
    # Highlight the starting vs ending point
    plt.annotate(f'Wild Type: {scores[0]:.2f}', xy=(generations[0], scores[0]), xytext=(5, -20), 
                 textcoords='offset points', arrowprops=dict(arrowstyle='->', color='red'))
    
    plt.annotate(f'Optimized: {scores[-1]:.2f}', xy=(generations[-1], scores[-1]), xytext=(-60, 10), 
                 textcoords='offset points', arrowprops=dict(arrowstyle='->', color='green'))

    # Save the plot
    filename = f"optimization_curve_{target_name}.png"
    plt.savefig(filename, dpi=300)
    print(f"--- SUCCESS: {filename} saved for your slides! ---")
    plt.show()

# Example Data for your demo
gens = list(range(1, 21))
# Mock data representing your scores climbing
scores = [0.45, 0.48, 0.55, 0.58, 0.62, 0.63, 0.68, 0.72, 0.75, 0.77, 
          0.79, 0.81, 0.82, 0.83, 0.85, 0.86, 0.88, 0.89, 0.91, 0.92]

generate_fitness_plot(gens, scores, "Monocyte")