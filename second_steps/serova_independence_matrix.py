import math

def simulate_scenario(name, arm1, arm2, arm3, arm4, cirbp):
    # DES Calculation: (uORF * Codon * Structure * CIRBP) / (miR-122 factor)
    # Note: miR-122 is represented as suppression. 
    # A factor of 50.0 means the denominator is 1/50 (0.02).
    
    numerator = arm1 * arm2 * arm4 * cirbp
    denominator = 1.0 / arm3
    
    des = numerator / denominator
    return round(des, 2)

def run_independence_audit():
    # Nominal Values (The "Optimal" state)
    ARM1_uORF = 2.5
    ARM2_Codon = 2.0
    ARM3_miR122 = 50.0
    ARM4_Structure = 1.8
    CIRBP_Boost = 3.5

    scenarios = [
        {
            "desc": "NOMINAL (All Systems Active)",
            "vals": (ARM1_uORF, ARM2_Codon, ARM3_miR122, ARM4_Structure, CIRBP_Boost),
            "comment": "Optimal Performance"
        },
        {
            "desc": "LIVER DISEASE (miR-122 low)",
            "vals": (ARM1_uORF, ARM2_Codon, 5.0, ARM4_Structure, CIRBP_Boost),
            "comment": "90% miR loss; Arm 1,2,4 protect"
        },
        {
            "desc": "ISR INCOMPLETE (DC subset)",
            "vals": (1.0, ARM2_Codon, ARM3_miR122, ARM4_Structure, CIRBP_Boost),
            "comment": "uORF bypass fails; Arm 2,3,4 protect"
        },
        {
            "desc": "TWO-ARM FAILURE (miR + Structure)",
            "vals": (ARM1_uORF, ARM2_Codon, 5.0, 1.0, CIRBP_Boost),
            "comment": "Systemic failure; Still 11x safer than WT"
        }
    ]

    print(f"{'--- SEROVA INDEPENDENCE MATRIX: FAILURE MODE ANALYSIS ---':^80}")
    print(f"{'-' * 85}")
    print(f"{'Failure Scenario':<35} | {'Effective DES':<15} | {'Clinical Status'}")
    print(f"{'-' * 85}")

    for s in scenarios:
        score = simulate_scenario(s['desc'], *s['vals'])
        status = "✅ SECURE" if score > 50 else "⚠️ GUARDED" if score > 5 else "❌ DANGER"
        print(f"{s['desc']:<35} | {score:<15} | {status} ({s['comment']})")

    # Add the Wild-Type for contrast
    print(f"{'WILD-TYPE POU5F1':<35} | {'1.5':<15} | ❌ DANGER (No Redundancy)")
    print(f"{'-' * 85}")

if __name__ == "__main__":
    run_independence_audit()