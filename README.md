# Project CUSTO for SEROVA: Multi-Objective mRNA Design Platform

**Target Cargo:** *POU5F1* (Oct4) for Precision Reprogramming  
**Status:** ⚠️ Computationally Validated — Wet-Lab Synthesis Pending  
**Team:** Dr. Andrew Newman · Nico Trummer · Filippo Conforto · Selin Abdullazade · Dr. Frida Arrey

<img width="1408" height="768" alt="CUSTO" src="https://github.com/user-attachments/assets/368d1b9e-182d-4f93-b21b-31d76f0ce920" />

---

## Overview
**Chain of Custody** is a computational design suite for engineering tissue-specific mRNA therapeutics. By integrating evolutionary algorithms, a 196-cell-type miRNA expression atlas, and multi-objective scoring, the platform designs *POU5F1* constructs that are biophysically stable, immunologically stealthy, and hepatically detargeted.

> [!IMPORTANT]
> All sequence-level claims are *in silico* predictions pending wet-lab validation.

---

##  The "Triple-Lock" Architecture
To solve "The Liver Problem" (LNP sequestration), the platform addresses safety through three independently generated layers concatenated into the final construct:

### 1. 5' UTR — ISR-Responsive uORF Gating
Evolved via **NSGA-II** to maximize translation efficiency while embedding weak-Kozak uORFs. 
* **Mechanism:** Traps ribosomes in normal liver conditions. 
* **Bypass:** In target cells undergoing an Integrated Stress Response (ISR), ribosomes bypass these traps to translate the cargo.

### 2. CDS — Multi-Objective Codon Optimization
The canonical CDS is optimized simultaneously across:
* **Stealth:** U-content < 5% (to evade immune detection).
* **Efficiency:** Codon Adaptation Index (CAI).
* **Stability:** GC windows kept between 30–70%.
* **Output:** Produces a **Pareto-optimal library** rather than a single solution.

### 3. 3' UTR — Database-Driven miRNA Site Selection
Built by querying a curated miRNA–cell type expression atlas (Patil et al. 2022) to identify miRNAs highly expressed in off-target tissues.
* **Binding Efficacy:** Follows McGeary et al. 2019.
* **Transparency:** Reports explicit coverage (e.g., 93.3% of off-targets covered) rather than a flat suppression figure.

---

## Innovation Highlights
* **Generalizable Detargeting:** Selects minimal miRNA sets covering the maximum off-target panel for any user-specified cell pair.
* **Transparent Gaps:** Reporting uncovered off-targets provides a stronger epistemic position for researchers.
* **Pareto Library:** Evaluates five normalized metrics (0–1), surfacing trade-off assumptions to the user.
* **DL Ready:** Designed to plug directly into foundation models like **Evo2**, **helix-mRNA**, and **RiboNN**.

---

##  Performance Benchmarking
*Scoring models validated against published experimental datasets (not output sequences).*

| Benchmark Dataset | Metric | Result | Status |
| :--- | :--- | :--- | :--- |
| **Sample et al. 2019** | Ribosome Load (MRL) | Pearson $r = 0.965$ | ✅ Validated |
| **Presnyak et al. 2015** | mRNA Half-Life | Pearson $r = 0.966$ | ✅ Validated |
| **Jain et al. 2018** | miR-122 Silencing | Pearson $r = 0.904$ | ✅ Validated |

---

## Roadmap to Scale
* **Phase 1 (W1–2):** Code refinement and pipeline stabilization.
* **Phase 2 (W2–4):** Synthesis of top Pareto candidates.
* **Phase 3 (W5–6):** **MPRA assays** with empirical feedback into kinetic constants.
* **Phase 4 (W7+):** Multi-gene, multi-tissue expansion for the longevity sector.

---

##  Repository Roadmap
- [x] **Step 1: The Engine** — Initial GA logic and codon swapping.
- [x] **Step 2: The Safety Lock** — Introduction of hard-gates and detargeting arms.
- [x] **Step 3: Deep Learning & Benchmarks** — Scoring model validation against literature.
- [x] **Step 4: Industrial Pipeline** — End-to-end synthesis-ready sequence generation.
- [x] **Step 5: Multi-Objective Optimization** — Pareto-front engineering across four simultaneous objectives.

## 🛠 The "Triple-Lock" Architecture
The platform addresses "The Liver Problem" through three orthogonal engineering layers:


graph TD
    A[Gene Input: POU5F1] --> B{Triple-Lock Engine}
    subgraph "The Three Locks"
    B --> C["5' UTR: ISR-Responsive uORF Gating"]
    B --> D["CDS: Multi-Objective Codon Opt (U < 5%)"]
    B --> E["3' UTR: Cooperative miR-122 Sites"]
    end
    C & D & E --> F[Pareto-Optimal Library]
    F --> G[Synthesis-Ready FASTA/GenBank]
---
## Project Structure
* `serova_pipeline_final.py`: Master entry point for sequence generation.
* `serova_ago2_model.py`: Cooperative silencing biophysics.
* `serova_uorf_designer.py`: 5'UTR gating logic and sequence generation.
* `SEROVA_sequence_final.fasta`: Final optimized nucleotide sequence.


##  References
* **Kozomara A et al. (2019).** *Nucleic Acids Research.*
* **Patil AH et al. (2022).** *GigaScience.*
* **McGeary SE et al. (2019).** *Science.*
* **Browder et al. (2022).** *Nature Aging.*

---
**Official Implementation:** [github.com/retr0ever/chainofcustody](https://github.com/retr0ever/chainofcustody)


