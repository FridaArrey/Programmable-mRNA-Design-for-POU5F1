Project CUSTO for SEROVA : Multi-Objective mRNA Design Platform
Target Cargo: POU5F1 (Oct4) for Precision Reprogramming | Status: ⚠️ Computationally Validated — Wet-Lab Synthesis Pending
Team: Dr. Andrew Newman (Neuronal Genomics) · Nico Trummer (Bioinformatics / Oncology) · Filippo Conforto (Ph.D. Candidate, Biophysics) · Selin Abdullazade (Software) · Dr. Frida Arrey (Immunologist)
Chain of Custody is a computational design suite for engineering tissue-specific mRNA therapeutics. By integrating evolutionary algorithms, a 196-cell-type miRNA expression atlas, and multi-objective scoring, the platform designs POU5F1 constructs that are biophysically stable, immunologically stealthy, and hepatically detargeted. All sequence-level claims are in silico predictions pending wet-lab validation.

The "Triple-Lock" Architecture
The core problem is that systemic LNPs sequester in the liver. POU5F1 expression there risks teratoma formation. The platform addresses this through three independently generated layers concatenated into the final construct: 5'UTR + CDS + 3'UTR.
5' UTR — ISR-Responsive uORF Gating. The 5' UTR is evolved via NSGA-II to maximise translation efficiency while embedding weak-Kozak uORFs that trap ribosomes in normal liver conditions. In target cells undergoing an Integrated Stress Response, ribosomes bypass these traps to translate POU5F1. Note that ISR is not cell-type exclusive — selectivity depends on the combined effect of all three locks, not this arm alone.
CDS — Multi-Objective Codon Optimization. The canonical CDS is fetched from Ensembl and optimised simultaneously across stealth (U-content < 5%), efficiency (CAI), stability (GC windows 30–70%), and differential expression proxy. NSGA-II produces a Pareto-optimal library rather than a single weighted-sum solution, allowing researchers to select a construct suited to their specific LNP vehicle.
3' UTR — Database-Driven miRNA Site Selection. The 3' UTR is built by querying a curated miRNA–cell type expression atlas (Patil et al. 2022, 196 cell types; Kozomara et al. 2019) to identify miRNAs highly expressed in off-target tissues and absent in the target. Binding site efficacy follows McGeary et al. 2019. The dashboard reports coverage explicitly — e.g., 93.3% of off-targets covered, 13 uncovered — rather than a single suppression figure. Coverage is a database prediction; actual silencing efficiency depends on Ago2 loading kinetics and LNP uptake not captured by RPM values alone.

What Is Innovative
Generalizable detargeting. Rather than hardcoding miR-122 as the liver off-switch, the platform selects the minimal miRNA set covering the maximum off-target tissue panel for any user-specified target/off-target cell pair.
Transparent coverage gaps. Explicitly reporting uncovered off-targets is a stronger epistemic position than a single suppression percentage, and directly useful for a researcher deciding whether a construct is ready for synthesis.
Pareto library with configurable weights. The five evaluation metrics (codon quality, cell-type selectivity, miRNA detargeting, structure, manufacturability) are normalized to 0–1 and combined with configurable weights — surfacing trade-off assumptions to the user rather than embedding them invisibly.
Deep learning ready. The scoring interface is designed to plug directly into foundation models (Evo2, helix-mRNA, RiboNN), replacing proxy scores with experimentally-trained predictions when available.


Performance Benchmarking
All benchmarks validate the scoring models against published experimental datasets — not the final optimized sequences directly.
Benchmark Dataset
Metric
Result
Status
Sample et al. 2019
Ribosome Load (MRL)
Pearson r = 0.965
✅ Model Validated
Presnyak et al. 2015
mRNA Half-Life
Pearson r = 0.966
✅ Model Validated
Jain et al. 2018
miR-122 Silencing
Pearson r = 0.904
✅ Model Validated

These correlations confirm that the underlying scoring models track known experimental outcomes well. They do not constitute direct experimental validation of Chain of Custody output sequences, which would require wet-lab synthesis and reporter assays.

Roadmap to Scale
The post-hackathon roadmap runs: Phase 1 (W1–2) code refinement → Phase 2 (W2–4) synthesis of top Pareto candidates → Phase 3 (W5–6) MPRA assays with empirical feedback into kinetic constants → Phase 4 (W7+) multi-gene, multi-tissue expansion. Phase 3 is the inflection point at which computational predictions become experimentally grounded claims.
The commercial target is the longevity biotech sector — Calico, Altos Labs, Retro Biosciences, NewLimit, Shift Bioscience — all pursuing OSKM partial reprogramming programmes where tissue-specific dosage control is the central unsolved safety problem (Browder et al. 2022, Nature Aging; Singh & Newman 2019, Epigenetics & Chromatin).


Why the Team Chose the Final chainofcustody Implementation: https://github.com/retr0ever/chainofcustody

After five iterative development phases, the team chose a CLI-first pipeline over the more elaborated research prototype. This was the right call for five reasons.
Earlier phases accumulated overclaimed metrics — DES Score 350.0, Tau Index 0.997, >99% hepatocyte suppression — that were artifacts of the scoring formula's construction, not experimental measurements. The final implementation replaces these with normalized 0–1 scores per metric. Two specific mechanisms were also quietly retired: the CIRBP 3.5x stability multiplier (temperature- and context-dependent, not reliably engineerable) and the uORF ISR-bypass as cell-type-specific (ISR is not exclusive to the target cell). The final tool scores structural accessibility and manufacturability as metrics, leaving biological interpretation to the researcher.
The benchmark Pearson correlations validate the scoring models against literature datasets, not the output sequences. The final implementation makes this distinction explicit — the evaluation module ranks candidates without asserting any sequence is synthesis-ready or clinically validated.
Usability is itself a form of scientific integrity. The composable fetch | evaluate CLI, batch ranking, multiple output formats, and interactive dashboard are immediately usable by a wet-lab collaborator. And configurable fitness weights mean that a researcher with a different LNP vehicle or cell target can adapt the scoring to their context without modifying the codebase — surfacing assumptions rather than hiding them.

References
Kozomara A et al. (2019). Nucleic Acids Research.
Patil AH et al. A curated human cellular microRNAome based on 196 primary cell types. (2022). GigaScience.
McGeary SE et al. The biochemical basis of microRNA targeting efficacy. (2019). Science.
Singh & Newman. (2019). Epigenetics & Chromatin.
Browder et al. (2022). Nature Aging.

Repository Roadmap
Step 1: The Engine — Initial GA logic and codon swapping.
Step 2: The Safety Lock — Introduction of hard-gates and detargeting arms.
Step 3: Deep Learning & Benchmarks — Scoring model validation against literature datasets.
Step 4: Industrial Pipeline — End-to-end synthesis-ready sequence generation.
Step 5: Multi-Objective Optimization — Pareto-front engineering across four simultaneous objectives.
