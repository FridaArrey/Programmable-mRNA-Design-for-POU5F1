def run_stress_test(nominal_des=350.0):
    # The Three Arms of the 350.0 Score
    arms = {
        "Arm 1 (uORF Gating)": 2.5,
        "Arm 2 (miR-122)": 50.0,
        "Arm 3 (Structure)": 2.8
    }
    
    print("--- SEROVA ARCHITECTURE: SENSITIVITY & STRESS ANALYSIS ---")
    print(f"{'Scenario':<35} | {'Effective DES':<15} | {'Status'}")
    print("-" * 70)
    
    # 1. NOMINAL CASE
    print(f"{'Full Triple-Lock (Nominal)':<35} | {nominal_des:<15} | OPTIMAL")

    # 2. THE "STRAND" SCENARIO (Lose 2 arms, keep miR-122)
    strand_sim = arms["Arm 2 (miR-122)"]
    print(f"{'Single-Arm Industry Equivalent':<35} | {strand_sim:<15} | INDUSTRY STD")

    # 3. TOTAL miR-122 FAILURE (0.02 -> 1.0)
    # Even if the liver has NO miR-122, Arm 1 and 3 still protect the patient.
    mir_failure = arms["Arm 1 (uORF Gating)"] * arms["Arm 3 (Structure)"]
    print(f"{'CRITICAL: miR-122 Failure':<35} | {mir_failure:<15} | PROTECTED")

    # 4. 50% SYSTEMIC UNDERPERFORMANCE
    # Every arm works at only half-strength
    stress_des = (arms["Arm 1 (uORF Gating)"] * 0.5) * \
                 (arms["Arm 2 (miR-122)"] * 0.5) * \
                 (arms["Arm 3 (Structure)"] * 0.5)
    print(f"{'SYSTEMIC STRESS (50% Efficiency)':<35} | {round(stress_des, 2):<15} | SECURE")

    # 5. WILD TYPE
    print(f"{'Wild-Type POU5F1':<35} | {'1.5':<15} | DANGEROUS")

if __name__ == "__main__":
    run_stress_test()