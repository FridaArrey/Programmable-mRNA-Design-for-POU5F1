"""
pou5f1_illustration.py — Canonical Scientific Illustration of Oct-4 / POU5F1
==============================================================================
Generates a publication-quality figure showing:
  1. Domain architecture of POU5F1 protein (linear map)
  2. POU domain binding to octamer DNA motif (3D-style schematic)
  3. Cell-type expression context (target vs off-target)
  4. mRNA construct overview (5'UTR → CDS → 3'UTR with miR-122 sites)

Output: pou5f1_illustration.png (300 DPI, suitable for slides/publication)

Requirements: matplotlib, numpy (both already installed)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Arc, Wedge
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

# ── Palette ───────────────────────────────────────────────────────────────────
BG        = "#0D1117"
PANEL_BG  = "#161B22"
GRID_LINE = "#21262D"
WHITE     = "#F0F6FC"
MUTED     = "#8B949E"
ACCENT1   = "#58A6FF"   # blue  — POU-specific domain
ACCENT2   = "#3FB950"   # green — POU homeodomain
ACCENT3   = "#FF7B72"   # red   — transactivation / disordered
ACCENT4   = "#D2A8FF"   # purple— linker
ACCENT5   = "#FFA657"   # orange— DNA
GOLD      = "#E3B341"
TEAL      = "#39D353"
PINK      = "#F778BA"

plt.rcParams.update({
    "font.family": "monospace",
    "text.color": WHITE,
    "axes.facecolor": PANEL_BG,
    "figure.facecolor": BG,
    "axes.edgecolor": GRID_LINE,
    "axes.labelcolor": WHITE,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "grid.color": GRID_LINE,
    "grid.linewidth": 0.5,
})

fig = plt.figure(figsize=(18, 13), facecolor=BG)
fig.patch.set_facecolor(BG)

gs = gridspec.GridSpec(
    3, 3,
    figure=fig,
    hspace=0.55, wspace=0.40,
    top=0.91, bottom=0.05, left=0.05, right=0.97
)

ax_title  = fig.add_subplot(gs[0, :])   # top: full-width title banner
ax_domain = fig.add_subplot(gs[1, :])   # middle: domain architecture
ax_dna    = fig.add_subplot(gs[2, 0])   # bottom-left: DNA binding
ax_cell   = fig.add_subplot(gs[2, 1])   # bottom-mid: cell expression
ax_mrna   = fig.add_subplot(gs[2, 2])   # bottom-right: mRNA construct

for ax in [ax_title, ax_domain, ax_dna, ax_cell, ax_mrna]:
    ax.set_facecolor(PANEL_BG)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_LINE)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 0 — Title banner
# ═══════════════════════════════════════════════════════════════════════════════
ax_title.set_xlim(0, 10)
ax_title.set_ylim(0, 1)
ax_title.axis("off")

# Background gradient strip
grad = np.linspace(0, 1, 256).reshape(1, -1)
ax_title.imshow(
    grad, aspect="auto", extent=[0, 10, 0, 1],
    cmap=LinearSegmentedColormap.from_list("hdr", [BG, "#1C2A3A", BG]),
    alpha=0.8, zorder=0
)

ax_title.text(
    5, 0.72, "OCT-4  /  POU5F1",
    ha="center", va="center", fontsize=28, fontweight="bold",
    color=ACCENT1, zorder=2,
    path_effects=[pe.withStroke(linewidth=6, foreground=BG)]
)
ax_title.text(
    5, 0.28,
    "Octamer-Binding Transcription Factor 4  ·  Yamanaka Reprogramming Factor  ·  Pluripotency Master Regulator",
    ha="center", va="center", fontsize=10, color=MUTED, zorder=2
)
# Decorative divider line
ax_title.axhline(0.08, color=ACCENT1, linewidth=0.8, alpha=0.4)
ax_title.axhline(0.92, color=ACCENT1, linewidth=0.8, alpha=0.4)

# Gene locus tag
ax_title.text(0.15, 0.5, "Gene: POU5F1", ha="left", va="center",
              fontsize=8, color=MUTED)
ax_title.text(0.15, 0.25, "Locus: 6p21.33", ha="left", va="center",
              fontsize=8, color=MUTED)
ax_title.text(9.85, 0.65, "360 aa", ha="right", va="center",
              fontsize=8, color=MUTED)
ax_title.text(9.85, 0.35, "MW: 38.6 kDa", ha="right", va="center",
              fontsize=8, color=MUTED)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 1 — Domain architecture (linear protein map)
# ═══════════════════════════════════════════════════════════════════════════════
ax = ax_domain
ax.set_xlim(0, 360)
ax.set_ylim(-1.8, 2.2)
ax.axis("off")

ax.text(180, 2.05, "P R O T E I N   D O M A I N   A R C H I T E C T U R E", ha="center", va="center",
        fontsize=10, color=MUTED, fontweight="bold")

# Backbone
ax.plot([0, 360], [0, 0], color=MUTED, linewidth=1.5, alpha=0.4, zorder=1)

# Domain definitions: (start, end, color, label, aa_range, y_label)
domains = [
    (1,   86,  ACCENT3, "N-terminal\nTransactivation\nDomain", "1–86",   1.05),
    (87,  133, ACCENT4, "Linker", "87–133", 1.05),
    (134, 210, ACCENT1, "POU-Specific\nDomain (POU-S)", "134–210", 1.05),
    (211, 228, ACCENT4, "POU\nLinker", "211–228", -1.35),
    (229, 295, ACCENT2, "POU\nHomeodomain\n(POU-HD)", "229–295", 1.05),
    (296, 360, ACCENT3, "C-terminal\nTransactivation\nDomain", "296–360", -1.35),
]

for (start, end, color, label, aa, y_lbl) in domains:
    width = end - start
    rect = FancyBboxPatch(
        (start, -0.3), width, 0.6,
        boxstyle="round,pad=2",
        facecolor=color, edgecolor=BG, linewidth=1.5,
        alpha=0.92, zorder=2
    )
    ax.add_patch(rect)

    mid = (start + end) / 2
    # Connector line to label
    ax.plot([mid, mid], [0.3 if y_lbl > 0 else -0.3, y_lbl * 0.72],
            color=color, linewidth=0.8, alpha=0.6, zorder=1)

    ax.text(mid, y_lbl + (0.3 if y_lbl > 0 else -0.25),
            label, ha="center", va="bottom" if y_lbl > 0 else "top",
            fontsize=7.5, color=WHITE, fontweight="bold",
            multialignment="center")
    ax.text(mid, y_lbl + (-0.25 if y_lbl > 0 else 0.3),
            f"aa {aa}", ha="center", va="top" if y_lbl > 0 else "bottom",
            fontsize=6.5, color=color, alpha=0.85)

# Tick marks every 50 aa
for i in range(0, 361, 50):
    ax.plot([i, i], [-0.35, -0.45], color=MUTED, linewidth=0.8, alpha=0.5)
    ax.text(i, -0.62, str(i), ha="center", va="top", fontsize=6.5, color=MUTED)

ax.text(180, -0.90, "Amino acid position", ha="center", va="top",
        fontsize=7, color=MUTED)

# Nuclear localisation signal annotation
ax.annotate("NLS", xy=(240, 0.3), xytext=(245, 1.55),
            fontsize=7, color=GOLD, ha="center",
            arrowprops=dict(arrowstyle="->", color=GOLD, lw=0.8))

# Scale bar
ax.plot([300, 360], [-1.55, -1.55], color=WHITE, linewidth=2, alpha=0.5)
ax.text(330, -1.70, "60 aa", ha="center", va="top", fontsize=6.5, color=MUTED)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 2 — DNA binding schematic
# ═══════════════════════════════════════════════════════════════════════════════
ax = ax_dna
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")
ax.set_title("DNA BINDING", fontsize=8, color=MUTED, pad=6, fontweight="bold")

# Double helix (simplified ribbon)
t = np.linspace(0, 4 * np.pi, 300)
x_helix = 5 + 1.8 * np.cos(t)
y_helix = 1 + 8 * t / (4 * np.pi)
x_helix2 = 5 + 1.8 * np.cos(t + np.pi)

ax.plot(x_helix,  y_helix, color=ACCENT5, linewidth=3, alpha=0.85, zorder=2)
ax.plot(x_helix2, y_helix, color=ACCENT5, linewidth=3, alpha=0.55, zorder=2, linestyle="--")

# Base pair rungs
step = len(t) // 16
for i in range(0, len(t), step):
    ax.plot([x_helix[i], x_helix2[i]], [y_helix[i], y_helix[i]],
            color=ACCENT5, linewidth=1, alpha=0.35, zorder=1)

# Octamer motif highlight box
y_oct_lo = 3.8
y_oct_hi = 6.2
oct_rect = FancyBboxPatch(
    (2.6, y_oct_lo), 4.8, y_oct_hi - y_oct_lo,
    boxstyle="round,pad=0.15",
    facecolor=GOLD, edgecolor=GOLD, linewidth=1.5,
    alpha=0.18, zorder=1
)
ax.add_patch(oct_rect)
ax.text(5, 6.55, "ATGCAAAT", ha="center", va="bottom",
        fontsize=8.5, color=GOLD, fontweight="bold",
        fontfamily="monospace")
ax.text(5, 6.95, "Octamer motif", ha="center", va="bottom",
        fontsize=6.5, color=GOLD, alpha=0.8)

# POU-S domain blob
pous_x = [2.2, 1.0, 1.2, 2.8, 3.8, 3.2, 2.2]
pous_y = [5.5, 5.0, 3.8, 3.2, 4.0, 5.2, 5.5]
ax.fill(pous_x, pous_y, color=ACCENT1, alpha=0.75, zorder=3)
ax.plot(pous_x, pous_y, color=WHITE, linewidth=0.8, alpha=0.4, zorder=4)
ax.text(2.2, 2.9, "POU-S", ha="center", va="top",
        fontsize=7, color=ACCENT1, fontweight="bold")

# POU-HD domain blob
pouhd_x = [7.8, 9.0, 8.8, 7.2, 6.2, 6.8, 7.8]
pouhd_y = [5.0, 4.5, 3.2, 2.8, 3.6, 4.8, 5.0]
ax.fill(pouhd_x, pouhd_y, color=ACCENT2, alpha=0.75, zorder=3)
ax.plot(pouhd_x, pouhd_y, color=WHITE, linewidth=0.8, alpha=0.4, zorder=4)
ax.text(7.8, 2.5, "POU-HD", ha="center", va="top",
        fontsize=7, color=ACCENT2, fontweight="bold")

# Binding arrows
ax.annotate("", xy=(3.5, 4.5), xytext=(3.0, 3.8),
            arrowprops=dict(arrowstyle="->", color=ACCENT1, lw=1.2))
ax.annotate("", xy=(6.5, 4.5), xytext=(7.0, 3.8),
            arrowprops=dict(arrowstyle="->", color=ACCENT2, lw=1.2))

ax.text(5, 0.4, "5'–ATGCAAAT–3'", ha="center", va="bottom",
        fontsize=7, color=ACCENT5, fontfamily="monospace")
ax.text(5, 0.05, "Octamer recognition sequence", ha="center", va="bottom",
        fontsize=6, color=MUTED)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 3 — Cell-type expression context
# ═══════════════════════════════════════════════════════════════════════════════
ax = ax_cell
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")
ax.set_title("CELL-TYPE SPECIFICITY", fontsize=8, color=MUTED, pad=6, fontweight="bold")

def draw_cell(ax, cx, cy, r, color, label, sublabel, expression, glow=False):
    """Draw a schematic cell with nucleus."""
    if glow:
        for ri, a in zip([r*1.35, r*1.20, r*1.08], [0.04, 0.08, 0.12]):
            c = plt.Circle((cx, cy), ri, color=color, alpha=a)
            ax.add_patch(c)
    # Cell membrane
    cell = plt.Circle((cx, cy), r, color=color, alpha=0.18, linewidth=1.5,
                       edgecolor=color, zorder=2)
    ax.add_patch(cell)
    # Nucleus
    nuc = plt.Circle((cx, cy + r*0.15), r * 0.42, color=color, alpha=0.55,
                      linewidth=1, edgecolor=color, zorder=3)
    ax.add_patch(nuc)
    # Expression level bar
    bar_w = r * 1.6
    bar_h = 0.22
    bar_x = cx - bar_w / 2
    bar_y = cy - r - 0.55
    # Background
    ax.add_patch(FancyBboxPatch((bar_x, bar_y), bar_w, bar_h,
                                boxstyle="round,pad=0.02",
                                facecolor=GRID_LINE, edgecolor="none"))
    # Fill
    ax.add_patch(FancyBboxPatch((bar_x, bar_y), bar_w * expression, bar_h,
                                boxstyle="round,pad=0.02",
                                facecolor=color, edgecolor="none", alpha=0.9))
    ax.text(cx, bar_y - 0.15, f"{int(expression*100)}% expression",
            ha="center", va="top", fontsize=6, color=color)
    ax.text(cx, cy + r + 0.25, label, ha="center", va="bottom",
            fontsize=8, color=WHITE, fontweight="bold")
    ax.text(cx, cy + r + 0.0, sublabel, ha="center", va="bottom",
            fontsize=6, color=MUTED)

# Dendritic cell (TARGET)
draw_cell(ax, 2.8, 6.0, 1.8, ACCENT2,
          "Dendritic Cell", "🎯 TARGET", 0.92, glow=True)
ax.text(2.8, 3.8, "HIGH EXPRESSION", ha="center", va="top",
        fontsize=6.5, color=ACCENT2, fontweight="bold")

# Hepatocyte (OFF-TARGET)
draw_cell(ax, 7.2, 6.0, 1.8, ACCENT3,
          "Hepatocyte", "🚫 OFF-TARGET", 0.04)
ax.text(7.2, 3.8, "SILENCED", ha="center", va="top",
        fontsize=6.5, color=ACCENT3, fontweight="bold")
ax.text(7.2, 3.5, "miR-122 ▶ degradation", ha="center", va="top",
        fontsize=6, color=MUTED)

# Central mRNA arrow
ax.annotate("", xy=(5.1, 6.0), xytext=(4.7, 6.0),
            arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1))
ax.text(5.0, 6.35, "mRNA\ndelivery", ha="center", va="bottom",
        fontsize=6, color=MUTED, multialignment="center")

# miR-122 icon
ax.text(5.0, 5.3, "miR-122", ha="center", va="center",
        fontsize=7, color=GOLD,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=PANEL_BG,
                  edgecolor=GOLD, linewidth=0.8))

# Bottom key metrics
ax.text(0.3, 0.9, f"DES = 3.91×", ha="left", va="center",
        fontsize=8, color=TEAL, fontweight="bold")
ax.text(0.3, 0.5, "Protein[DC] / Protein[Hepatocyte]", ha="left", va="center",
        fontsize=6, color=MUTED)
ax.text(0.3, 0.15, "Kinetic model · Bartel 2009 · Chang 2004", ha="left", va="center",
        fontsize=5.5, color=MUTED, alpha=0.7)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 4 — mRNA construct map
# ═══════════════════════════════════════════════════════════════════════════════
ax = ax_mrna
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")
ax.set_title("mRNA CONSTRUCT DESIGN", fontsize=8, color=MUTED, pad=6, fontweight="bold")

# mRNA backbone
y_mrna = 5.8
ax.plot([0.3, 9.7], [y_mrna, y_mrna], color=MUTED, linewidth=1.2,
        alpha=0.4, solid_capstyle="round")

# 5' cap
cap = plt.Circle((0.3, y_mrna), 0.22, color=PINK, alpha=0.9, zorder=3)
ax.add_patch(cap)
ax.text(0.3, y_mrna, "m7G", ha="center", va="center",
        fontsize=5, color=BG, fontweight="bold")
ax.text(0.3, y_mrna - 0.55, "5' Cap", ha="center", va="top",
        fontsize=6, color=PINK)

# Poly-A tail
ax.text(9.7, y_mrna, "AAAA…", ha="right", va="center",
        fontsize=7, color=PINK, alpha=0.9)
ax.text(9.7, y_mrna - 0.55, "Poly-A", ha="right", va="top",
        fontsize=6, color=PINK)

# Regions
regions = [
    (0.7,  2.0, ACCENT4, "5'UTR",    "uORF\nstructure"),
    (2.0,  7.2, ACCENT1, "CDS · POU5F1", "codon-optimised\nCAI↑ U↓"),
    (7.2,  9.4, ACCENT2, "3'UTR",    "miR-122\nsites ×3"),
]

for (x0, x1, color, label, sublabel) in regions:
    w = x1 - x0
    rect = FancyBboxPatch(
        (x0, y_mrna - 0.28), w, 0.56,
        boxstyle="round,pad=0.05",
        facecolor=color, edgecolor=BG, linewidth=1,
        alpha=0.85, zorder=2
    )
    ax.add_patch(rect)
    mid = (x0 + x1) / 2
    ax.text(mid, y_mrna, label, ha="center", va="center",
            fontsize=6.5, color=BG, fontweight="bold")

    # Labels above/below alternating
    is_cds = "CDS" in label
    y_sub = y_mrna + (0.55 if not is_cds else -0.65)
    va = "bottom" if not is_cds else "top"
    ax.plot([mid, mid], [y_mrna + (0.28 if not is_cds else -0.28), y_sub * 0.998],
            color=color, linewidth=0.7, alpha=0.6)
    ax.text(mid, y_sub, sublabel, ha="center", va=va,
            fontsize=6, color=color, multialignment="center")

# miR-122 site markers
mir_xs = [7.55, 8.0, 8.45]
for mx in mir_xs:
    ax.plot([mx, mx], [y_mrna - 0.28, y_mrna - 0.85], color=GOLD,
            linewidth=1.0, alpha=0.8)
    ax.plot(mx, y_mrna - 0.9, "v", color=GOLD, markersize=5, alpha=0.9)
ax.text(8.0, y_mrna - 1.15, "miR-122 seed sites", ha="center", va="top",
        fontsize=6, color=GOLD)

# Optimisation metrics column
metrics = [
    ("U-content",     "↓ 11%",  ACCENT3),
    ("GC content",    "↑ 3%",   ACCENT2),
    ("CAI",           "↑ 0.78", ACCENT1),
    ("miR-122 sites", "×3",     GOLD),
    ("Liver Detarget","89.8%",  TEAL),
    ("Composite",     "0.749",  WHITE),
]

ax.text(5.0, 4.2, "─── Optimisation Gains ───", ha="center", va="top",
        fontsize=6.5, color=MUTED)

for i, (metric, value, color) in enumerate(metrics):
    y = 3.7 - i * 0.52
    ax.text(1.5, y, metric, ha="left", va="center",
            fontsize=6.5, color=MUTED)
    ax.text(8.5, y, value, ha="right", va="center",
            fontsize=7, color=color, fontweight="bold")
    ax.plot([3.8, 6.2], [y, y], color=GRID_LINE, linewidth=0.5, alpha=0.5)

# DES result
ax.add_patch(FancyBboxPatch(
    (1.0, 0.15), 8.0, 0.75,
    boxstyle="round,pad=0.1",
    facecolor=TEAL, edgecolor=TEAL, linewidth=1, alpha=0.12
))
ax.text(5.0, 0.52, "DES: 0.43×  →  3.91×    (9× improvement over wild-type)",
        ha="center", va="center", fontsize=7, color=TEAL, fontweight="bold")

# ── Final touches ─────────────────────────────────────────────────────────────
fig.text(
    0.5, 0.01,
    "Oct-4 / POU5F1 · SEROVA mRNA Design Pipeline · "
    "Refs: Chang 2004 · Bartel 2009 · Sample 2019 · Jain 2018 · Karikó 2012",
    ha="center", va="bottom", fontsize=6, color=MUTED, alpha=0.7
)

plt.savefig(
    "pou5f1_illustration.png",
    dpi=300, bbox_inches="tight",
    facecolor=BG, edgecolor="none"
)
print("✅ Saved: pou5f1_illustration.png (300 DPI)")
plt.show()