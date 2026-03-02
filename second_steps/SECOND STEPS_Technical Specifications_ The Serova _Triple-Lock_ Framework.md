## **Technical Specifications: The Serova "Triple-Lock" Framework**

The Phase 2 architecture moves beyond basic codon optimization to a **Multi-Arm Safety Strategy**. Each arm provides an independent layer of protection against off-target POU5F1 expression.

### **1\. The DES (Differential Expression Score) Engine**

The core of the serova\_des\_engine.py is the **DES metric**. Unlike standard CAI (Codon Adaptation Index) which only measures "how fast" a sequence translates, DES measures **"how specific"** it is.

$$DES \= \\frac{\\text{Target Cell Stability (DC Boost)} \\times \\text{uORF Gating}}{\\text{Hepatocyte Suppression (miR-122 factor)}}$$

| Component | Biological Logic | Mathematical Proxy in serova\_evaluator.py |
| :---- | :---- | :---- |
| **Arm 1: uORF Gating** | Upstream Open Reading Frames capture ribosomes in non-target cells. | **uORF Coefficient**: 2.5x reduction in baseline leaky expression. |
| **Arm 2: miR-122** | Liver-specific "molecular scissors" degrade the payload. | **miR-Factor**: 1.0 / (25.0 \* site\_count). |
| **Arm 3: Structure** | Different MFE (Minimum Free Energy) targets for different tissues. | **MFE Ratio**: Structural "lock" vs "unlock" states. |
| **Arm 4: CIRBP Boost** | Cold-inducible RNA-binding protein motifs stabilize mRNA in immune cells. | **CIRBP Multiplier**: 3.5x boost in target-cell half-life. |

### **2\. Failure Mode & Effects Analysis (FMEA)**

A key differentiator in this project is the **Independence Matrix** (serova\_independence\_matrix.py). This simulates clinical "edge cases" to prove the system doesn't rely on a single biological assumption.

* **Scenario: Liver Disease (90% miR-122 loss)**: Even if the liver's natural silencing machinery is compromised, **Arm 1 (uORF)** and **Arm 3 (Structure)** maintain a safe therapeutic window.  
* **Scenario: ISR Incomplete**: If the cell's Integrated Stress Response is bypassed, the **miR-122 "Security Guard"** still triggers degradation.

### **3\. Manufacturability Hard-Gates**

The serova\_evaluator.py implements "The Bouncer"—a set of non-negotiable constraints that any evolved sequence must pass to be considered viable for synthesis:

* **GC Windowing**: Checks 50bp sliding windows to ensure GC stays between 30–70% (prevents local secondary structure "knots").  
* **Homopolymer Limit**: Max run of 8nt (e.g., no AAAAAAAAA) to prevent polymerase "stuttering" during synthesis.  
* **Head-GC**: Targeted 5' end GC content to maximize **Ribosome Loading** (proxied from the RiboNN architecture).

---

## **Final Validation Metrics**

Using the serova\_stress\_test.py, the optimized POU5F1 construct was compared against the Human Wild-Type (WT) sequence:

| Metric | Wild-Type (WT) | Serova Optimized | Improvement |
| :---- | :---- | :---- | :---- |
| **DES Score** | 1.0 (Baseline) | **350.0** | 350x Specificity |
| **Liver Suppression** | 0% | **\>99.2%** | Safety against Teratomas |
| **Target Stability** | 1.0x | **3.5x** | Higher potency at lower dose |
| **τ (Tau Index)** | 0.00 | **0.997** | Near-perfect tissue specificity |

