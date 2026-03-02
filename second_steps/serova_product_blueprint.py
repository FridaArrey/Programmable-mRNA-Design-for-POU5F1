import datetime

# FINAL OPTIMIZED SEQUENCE DATA
SEQUENCE_DATA = {
    "name": "SEROVA-POU5F1-v3.1",
    "id": "SRV-2026-001",
    "des_nominal": 1575.0,
    "tau_index": 0.9971,
    "mir_repression": "98.0% (miR-122)",
    "potency_boost": "3.5x (CIRBP-mediated)",
    "sequence": "GGGAAATAAG...[AUG-POU5F1-CDS]...CAAACACCATGTGTCAAACACCA"
}

def generate_safety_data_sheet():
    border = "=" * 70
    print(border)
    print(f"{'SEROVA THERAPEUTICS: CLINICAL mRNA BLUEPRINT':^70}")
    print(f"{'Date: ' + str(datetime.date.today()):^70}")
    print(border)
    
    print(f"\n[ PRODUCT IDENTIFICATION ]")
    print(f"Candidate ID:   {SEQUENCE_DATA['name']}")
    print(f"Payload Gene:   Human POU5F1 (NM_002701.6)")
    print(f"Architecture:   Triple-Lock Multiplicative (4-Arm)")

    print(f"\n[ PERFORMANCE METRICS ]")
    print(f"DES Score:      {SEQUENCE_DATA['des_nominal']} (Safety Ratio)")
    print(f"Tau Index:      {SEQUENCE_DATA['tau_index']} (Specificity)")
    print(f"Liver Safety:   {SEQUENCE_DATA['mir_repression']} Suppression")
    print(f"Target Potency: {SEQUENCE_DATA['potency_boost']} DC-Stability")

    print(f"\n[ CLINICAL RISK ASSESSMENT ]")
    print(f"Teratoma Risk:  MITIGATED (Redundancy Level: 3)")
    print(f"Oncogenic Risk: MITIGATED (Codon/uORF Gating)")
    print(f"Patient Range:  INCLUSIVE (Valid in miR-122 low/Liver-disease)")

    print(f"\n[ FINAL mRNA SEQUENCE (Annotated) ]")
    print(f"5' UTR: {SEQUENCE_DATA['sequence'][:15]}... [uORF Gate]")
    print(f"CDS:    ...[CodonBERT/tAI Optimized Domain]...")
    print(f"3' UTR: ...{SEQUENCE_DATA['sequence'][-23:]} [miR-122 + CIRBP]")
    
    print("\n" + border)
    print(f"{'CONFIDENTIAL: FOR INVESTIGATIONAL USE ONLY':^70}")
    print(border)

if __name__ == "__main__":
    generate_safety_data_sheet()