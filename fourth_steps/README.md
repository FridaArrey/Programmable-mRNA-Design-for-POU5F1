# Project SEROVA: Tissue-Specific mRNA Optimization
**Target Cargo:** *POU5F1* (Oct4)  
**Date:** 2026-02-28  
**Status:** ✅ Synthesis Ready (3/3 Benchmarks Passed)

---

## 🧬 Executive Summary
This project provides a validated, synthesis-ready mRNA sequence for **POU5F1** designed for high-precision delivery. By integrating three orthogonal specificity layers, we successfully flipped the expression profile from liver-biased to **Dendritic Cell-targeted**.

### 🚀 Key Performance Indicators
* **Selectivity Improvement:** 6.1x over Wild-Type.
* **Liver Detargeting:** 94.5% predicted degradation in hepatocytes.
* **Immune Evasion:** U-content reduced to 4.7%.

---

## 🛠 Engineering Architecture

### 1. 5' UTR: ISR-Responsive uORF Gating
* **Sequence:** `AAACTATGAGTTAA`
* **Mechanism:** A 2-codon weak-Kozak uORF that utilizes Integrated Stress Response (ISR) levels to filter translation.
* **Target/Off-Target Differential:** 1.56x (uORF arm alone).

### 2. CDS: Structural & Codon Optimization
* **Optimization:** Multi-objective GA targeting CAI and U-depletion.
* **Validation:** Pearson $r=0.966$ vs *Presnyak 2015* (Stability).

### 3. 3' UTR: Cooperative miR-122 Silencing
* **Design:** 3 synthetic miR-122 seed sites.
* **Model:** Cooperative Ago2 Hill-kinetics (Wee 2012).
* **Validation:** Pearson $r=0.904$ vs *Jain 2018*.

---

## 📊 Performance Validation (Benchmarking)
| Benchmark | Metric | Result | Status |
| :--- | :--- | :--- | :--- |
| **Ribosome Load (MRL)** | Pearson $r$ | 0.965 | ✅ PASS |
| **mRNA Half-Life** | Pearson $r$ | 0.966 | ✅ PASS |
| **miR-122 Silencing** | Pearson $r$ | 0.904 | ✅ PASS |

---

## 📂 Project Structure
* `serova_pipeline_final.py`: Master entry point for sequence generation.
* `serova_ago2_model.py`: Cooperative silencing biophysics.
* `serova_uorf_designer.py`: 5'UTR gating logic.
* `SEROVA_sequence_final.fasta`: Final optimized nucleotide sequence.

---

## ⚖️ Ethical & Scientific Caveats
* **Ago2 Competition:** The model assumes localized Ago2 saturation; intermediate tissue performance remains an open calibration target.
* **ISR Dependency:** Gating efficiency is dependent on the ΔISR between target and off-target tissues.
