#!/usr/bin/env python3
"""
Generate Polished Overview, Methodology & Dashboard User Guide Figures:
- Master Figure 0: Unified 4-Panel Master Poster (summary_figure0_overview_and_methodology.png / .pdf)
- Standalone Figure 0A: The Human Globin Paradigm & Evolutionary Inquiry (summary_figure0a_human_paradigm.png / .pdf)
- Standalone Figure 0B: Synteny-Anchored Locus Retrieval Pipeline & 4 Archetypes (summary_figure0b_synteny_pipeline.png / .pdf)
- Standalone Figure 0C: Anatomy of Genomic Metrics, Biophysics & Synonymous GC3 Codon Shift (summary_figure0c_metric_anatomy_and_gc3.png / .pdf)
- Standalone Figure 0D: Interactive Dashboard Architecture & Analytical User Guide (summary_figure0d_dashboard_guide.png / .pdf)
"""

import os
import sys
import shutil
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.lines import Line2D

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
SYNC_DIR = os.path.join(BASE_DIR, "1_data_from_mike", "results")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(SYNC_DIR, exist_ok=True)


def draw_card(ax, x, y, w, h, title="", fc="#ffffff", ec="#cbd5e1", lw=1.2, radius=0.012, zorder=1):
    """Draw a clean card container with an optional header."""
    box = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.002,rounding_size={radius}",
                         facecolor=fc, edgecolor=ec, lw=lw, zorder=zorder)
    ax.add_patch(box)
    if title:
        header_h = 0.040
        hdr = FancyBboxPatch((x, y + h - header_h), w, header_h,
                             boxstyle=f"round,pad=0.002,rounding_size={radius}",
                             facecolor="#f1f5f9", edgecolor=ec, lw=lw, zorder=zorder + 1)
        ax.add_patch(hdr)
        ax.text(x + w / 2.0, y + h - header_h / 2.0, title,
                ha='center', va='center', fontsize=9.8, fontweight='bold', color='#0f172a', zorder=zorder + 2)


def draw_gene_arrow(ax, x, y, w, h, name, fc="#0284c7", ec="#0369a1", text_col="white", zorder=4, fontsize=7.6):
    """Draw a directed gene arrow box."""
    head_w = min(0.015, w * 0.25)
    body_w = w - head_w
    verts = [
        (x, y - h/2.0),
        (x + body_w, y - h/2.0),
        (x + w, y),
        (x + body_w, y + h/2.0),
        (x, y + h/2.0),
        (x, y - h/2.0)
    ]
    poly = Polygon(verts, closed=True, facecolor=fc, edgecolor=ec, lw=1.1, zorder=zorder)
    ax.add_patch(poly)
    ax.text(x + body_w/2.0, y, name, ha='center', va='center', fontsize=fontsize, fontweight='bold', color=text_col, zorder=zorder+1)


# =========================================================================
# PANEL A: THE HUMAN GLOBIN PARADIGM & EVOLUTIONARY INQUIRY
# =========================================================================
def draw_panel_a(ax, standalone=False):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    t_size = 14.5 if standalone else 12.5
    s_size = 10.0 if standalone else 8.5

    ax.text(0.02, 0.965, "A. The Human Globin Paradigm & The Core Evolutionary Inquiry",
            fontsize=t_size, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Canonical dual-cluster arrangement on chromosomes 16 and 11, and the hypothesis of universal isochore polarity",
            fontsize=s_size, fontstyle='italic', color='#475569', va='top')

    # Card 1: Chr 16 - Alpha Globin Cluster (Subtelomeric Heavy Isochore)
    draw_card(ax, 0.02, 0.525, 0.96, 0.39, title="Human Alpha-Globin Cluster — Subtelomeric Heavy Isochore (16p13.3)",
              fc="#fdfefe", ec="#93c5fd", lw=1.4)
    
    iso_a = FancyBboxPatch((0.035, 0.535), 0.93, 0.33, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f0f9ff", edgecolor="none", zorder=1)
    ax.add_patch(iso_a)
    ax.text(0.95, 0.84, "Regional Heavy Isochore H3 (GC > 54%)", ha='right', fontsize=8.0, fontweight='bold', color='#0284c7', zorder=2)

    # Chr 16 Ideogram
    ax.plot([0.06, 0.92], [0.805, 0.805], color="#cbd5e1", lw=11, solid_capstyle='round', zorder=2)
    ax.plot([0.42, 0.42], [0.785, 0.825], color="#475569", lw=3.0, zorder=3)
    ax.text(0.42, 0.832, "Cen (35.3 Mb)", fontsize=7.5, ha='center', color='#64748b', fontweight='bold')
    ax.plot([0.06, 0.06], [0.785, 0.825], color="#dc2626", lw=3.5, zorder=4)
    ax.text(0.06, 0.832, "16p Telomere", fontsize=8.0, ha='left', color='#dc2626', fontweight='bold')
    
    a_x = 0.075
    ax.plot([a_x, a_x], [0.805, 0.745], color="#0284c7", lw=1.8, zorder=5)
    ax.scatter([a_x], [0.805], color="#0284c7", s=36, zorder=6)
    ax.text(a_x, 0.732, "Alpha Locus: 0.15–0.22 Mb (0.16% from telomere — Subtelomeric)", fontsize=8.0, ha='left', va='top',
            color='#0369a1', fontweight='bold')

    ax.plot([0.06, 0.92], [0.635, 0.635], color="#cbd5e1", lw=1.5, zorder=2)
    
    anchors_alpha = [
        ("NPRL3", 0.07, 0.075, "#7c3aed"),
        ("MPG", 0.16, 0.060, "#7c3aed"),
        ("RHBDF1", 0.235, 0.075, "#7c3aed"),
        ("LUC7L", 0.81, 0.075, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_alpha:
        draw_gene_arrow(ax, gx, 0.635, gw, 0.048, gname, fc=gcol, ec="#5b21b6", text_col="white")
        ax.text(gx + gw/2.0, 0.595, "Anchor", ha='center', va='top', fontsize=7.0, color='#6b21a8')

    globins_alpha = [
        ("HBZ (ζ)", 0.33, 0.075, "#0284c7"),
        ("HBM (μ)", 0.42, 0.070, "#0284c7"),
        ("HBA2 (α2)", 0.505, 0.075, "#1d4ed8"),
        ("HBA1 (α1)", 0.595, 0.075, "#1d4ed8"),
        ("HBQ1 (θ)", 0.685, 0.075, "#0284c7"),
    ]
    for gname, gx, gw, gcol in globins_alpha:
        draw_gene_arrow(ax, gx, 0.635, gw, 0.054, gname, fc=gcol, ec="#1e3a8a", text_col="white")

    m_box_a = FancyBboxPatch((0.06, 0.540), 0.86, 0.042, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#eff6ff", edgecolor="#bfdbfe", lw=1, zorder=3)
    ax.add_patch(m_box_a)
    ax.text(0.49, 0.561,
            "α-Metrics:  CDS GC = 66.2%  |  GC3 = 91.1%  |  Gene GC = 68.2%  |  Mean Gene Length = 810 bp  |  Isochore: Heavy H3",
            ha='center', va='center', fontsize=8.2, fontweight='bold', color='#1e40af')

    # Card 2: Chr 11 - Beta Globin Cluster (Interstitial Light Isochore)
    draw_card(ax, 0.02, 0.115, 0.96, 0.39, title="Human Beta-Globin Cluster — Interstitial Light Isochore (11p15.4)",
              fc="#fdfefe", ec="#86efac", lw=1.4)
    
    iso_b = FancyBboxPatch((0.035, 0.125), 0.93, 0.33, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f0fdf4", edgecolor="none", zorder=1)
    ax.add_patch(iso_b)
    ax.text(0.95, 0.43, "Regional Light Isochore L1/L2 (GC ~40%)", ha='right', fontsize=8.0, fontweight='bold', color='#16a34a', zorder=2)

    # Chr 11 Ideogram
    ax.plot([0.06, 0.92], [0.395, 0.395], color="#cbd5e1", lw=11, solid_capstyle='round', zorder=2)
    ax.plot([0.45, 0.45], [0.375, 0.415], color="#475569", lw=3.0, zorder=3)
    ax.text(0.45, 0.422, "Cen (51.6 Mb)", fontsize=7.5, ha='center', color='#64748b', fontweight='bold')
    ax.plot([0.06, 0.06], [0.375, 0.415], color="#dc2626", lw=3.5, zorder=4)
    ax.text(0.06, 0.422, "11p Telomere", fontsize=8.0, ha='left', color='#dc2626', fontweight='bold')
    
    b_x = 0.20
    ax.plot([b_x, b_x], [0.395, 0.335], color="#16a34a", lw=1.8, zorder=5)
    ax.scatter([b_x], [0.395], color="#16a34a", s=36, zorder=6)
    ax.text(b_x, 0.322, "Beta Locus: 5.22–5.27 Mb (3.85% from telomere — Interstitial Desert)", fontsize=8.0, ha='left', va='top',
            color='#15803d', fontweight='bold')

    ax.plot([0.06, 0.92], [0.225, 0.225], color="#cbd5e1", lw=1.5, zorder=2)

    anchors_beta = [
        ("OR51V1", 0.065, 0.070, "#7c3aed"),
        ("OR51B4", 0.145, 0.070, "#7c3aed"),
        ("TRIM21", 0.770, 0.070, "#7c3aed"),
        ("RRM1", 0.850, 0.070, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_beta:
        draw_gene_arrow(ax, gx, 0.225, gw, 0.048, gname, fc=gcol, ec="#5b21b6", text_col="white")
        ax.text(gx + gw/2.0, 0.185, "Anchor", ha='center', va='top', fontsize=7.0, color='#6b21a8')

    globins_beta = [
        ("HBE1 (ε)", 0.245, 0.075, "#15803d"),
        ("HBG2 (Gγ)", 0.335, 0.075, "#15803d"),
        ("HBG1 (Aγ)", 0.425, 0.075, "#15803d"),
        ("HBD (δ)", 0.535, 0.070, "#16a34a"),
        ("HBB (β)", 0.620, 0.075, "#16a34a"),
    ]
    for gname, gx, gw, gcol in globins_beta:
        draw_gene_arrow(ax, gx, 0.225, gw, 0.054, gname, fc=gcol, ec="#14532d", text_col="white")

    m_box_b = FancyBboxPatch((0.06, 0.130), 0.86, 0.042, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#f0fdf4", edgecolor="#bbf7d0", lw=1, zorder=3)
    ax.add_patch(m_box_b)
    ax.text(0.49, 0.151,
            "β-Metrics:  CDS GC = 54.4%  |  GC3 = 65.8%  |  Gene GC = 43.4%  |  Mean Gene Length = 1,608 bp  |  Isochore: Light L1",
            ha='center', va='center', fontsize=8.2, fontweight='bold', color='#166534')

    # Core Evolutionary Inquiry Callout Box
    q_box = FancyBboxPatch((0.02, 0.015), 0.96, 0.085, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f8fafc", edgecolor="#94a3b8", lw=1.2, linestyle='--', zorder=3)
    ax.add_patch(q_box)
    ax.text(0.04, 0.072, "The Core Macroevolutionary Inquiry:", fontsize=8.8, fontweight='bold', color='#0f172a', zorder=4)
    ax.text(0.04, 0.048, "1. Is this extreme GC disparity (subtelomeric high-GC α vs. interstitial low-GC β) an invariant vertebrate feature?",
            fontsize=7.8, color='#334155', zorder=4)
    ax.text(0.04, 0.026, "2. Did physical distance to the telomere dictate locus GC content via GC-biased gene conversion (gBGC) across 500 Myr?",
            fontsize=7.8, color='#334155', zorder=4)


# =========================================================================
# PANEL B: SYNTENIC RETRIEVAL PIPELINE & 4 ARCHETYPES
# =========================================================================
def draw_panel_b(ax, standalone=False):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    t_size = 14.5 if standalone else 12.5
    s_size = 10.0 if standalone else 8.5

    ax.text(0.02, 0.965, "B. Synteny-Anchored Locus Retrieval Pipeline Across 68 Genomes",
            fontsize=t_size, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Identification of homologous intervals via conserved syntenic anchors and discovery of 4 macroevolutionary archetypes",
            fontsize=s_size, fontstyle='italic', color='#475569', va='top')

    draw_card(ax, 0.02, 0.770, 0.96, 0.150, title="Synteny Discovery Pipeline & Chromosome Contiguity Filtering",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)
    
    steps = [
        ("1. Syntenic Search", "Query 68 genomes for flanking markers:\nα: NPRL3, MPG, RHBDF1, LUC7L\nβ: OR51 cluster, CDCA8, TRIM21, RRM1", "#3b82f6"),
        ("2. Interval Bounds", "Define syntenic window containing all\nflankers and globin duplications.\nExtract verified coding models.", "#8b5cf6"),
        ("3. Scaffold vs Chr", "Discriminate assembled chromosomes\nfrom unplaced scaffolds to avoid\nartificial telomeric misassignment.", "#f59e0b"),
        ("4. Baseline Calibration", "Compute whole-genome CDS & gene\nbaselines (mean, SD, percentiles)\nto calculate locus z-score deviations.", "#10b981")
    ]
    for idx, (stitle, sdesc, scol) in enumerate(steps):
        sx = 0.04 + idx * 0.235
        sbox = FancyBboxPatch((sx, 0.782), 0.220, 0.095, boxstyle="round,pad=0.002,rounding_size=0.006",
                              facecolor="#f8fafc", edgecolor=scol, lw=1.3, zorder=3)
        ax.add_patch(sbox)
        ax.text(sx + 0.012, 0.862, stitle, fontsize=8.0, fontweight='bold', color=scol, zorder=4)
        ax.text(sx + 0.012, 0.840, sdesc, fontsize=6.8, color='#334155', va='top', zorder=4)

    draw_card(ax, 0.02, 0.015, 0.96, 0.740, title="The 4 Core Macroevolutionary Archetypes Across Vertebrate Evolution",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    archetypes = [
        ("1. Linked Primordial Archetype (Basal Actinopterygii & Amphibians — Spotted Gar, Xenopus)",
         "Ancestral tandem arrangement: α and β globin genes reside on the same chromosome in a single cluster.",
         "Both genes share identical genomic environment and regional recombination rate. ΔGC3 (α - β) ≈ 0%.",
         "#1565c0", "#f0fdf4", 0.575, "Single Chromosome (Gar Chr 19 / Xenopus Chr 1)", "LINKED"),

        ("2. Canonical Decoupled Subtelomeric-α (Eutherian & Metatherian Mammals — Human, Mouse)",
         "Amniote translocation decouples α and β. α captured by subtelomere (<1 Mb to telomere); β remains interstitial.",
         "Runaway subtelomeric gBGC supercharges α-GC3 (>85–95%) while β-GC3 remains moderate (~60–65%). ΔGC3 >> 0%.",
         "#e65100", "#fff7ed", 0.390, "Two Chromosomes (α Subtelomeric / β Interstitial)", "MAMMALIAN"),

        ("3. Dual Microchromosomal Partitioning (Archosauria & Aves — Chicken, Eagle, Hummingbird)",
         "Both α and β clusters translocate onto distinct, highly recombinogenic microchromosomes (Chr 14 and Chr 25).",
         "High crossover density drives synchronous hyper-elevation of both clusters: α-GC3 > 85%, β-GC3 > 85%. ΔGC3 ≈ 0%.",
         "#00796b", "#f0fdfa", 0.205, "Two Microchromosomes (Micro-14 and Micro-25)", "AVIAN"),

        ("4. The Great Inversion (Squamate Reptiles — Lizards & Snakes — Wall Lizard, Garter Snake)",
         "Karyotypic inversion translocates β-locus to microchromosome/subtelomere, leaving α in an interstitial desert.",
         "Inverted gBGC polarity: β-GC3 skyrockets (>90–96%) while α-GC3 languishes (~48–55%). ΔGC3 << 0% (β >> α)!",
         "#d32f2f", "#fef2f2", 0.020, "Inverted Polarity (β Subtelomeric / α Interstitial)", "SQUAMATE"),
    ]

    for atitle, adesc1, adesc2, acol, afc, ay, achr, atype in archetypes:
        abox = FancyBboxPatch((0.04, ay), 0.92, 0.165, boxstyle="round,pad=0.002,rounding_size=0.008",
                              facecolor=afc, edgecolor=acol, lw=1.2, zorder=3)
        ax.add_patch(abox)
        ax.text(0.06, ay + 0.145, atitle, fontsize=8.6, fontweight='bold', color=acol, va='top', zorder=4)
        ax.text(0.06, ay + 0.122, adesc1, fontsize=7.6, color='#334155', va='top', zorder=4)
        ax.text(0.06, ay + 0.100, adesc2, fontsize=7.6, fontweight='bold', color='#1e293b', va='top', zorder=4)
        ax.text(0.06, ay + 0.078, f"Cytogenetic State: {achr}", fontsize=7.3, color='#64748b', va='top', zorder=4)

        mini_y = ay + 0.035
        ax.plot([0.06, 0.88], [mini_y, mini_y], color="#cbd5e1", lw=8, solid_capstyle='round', zorder=4)
        
        if atype == "LINKED":
            ax.plot([0.42, 0.42], [mini_y - 0.012, mini_y + 0.012], color="#475569", lw=2, zorder=5)
            ax.add_patch(FancyBboxPatch((0.70, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.add_patch(FancyBboxPatch((0.745, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.722, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.767, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.81, mini_y, "◀ Linked α-β (ΔGC3 ≈ 0)", fontsize=7.2, fontweight='bold', color='#1565c0', va='center', zorder=7)
            
        elif atype == "MAMMALIAN":
            ax.plot([0.06, 0.06], [mini_y - 0.012, mini_y + 0.012], color="#dc2626", lw=2.5, zorder=5)
            ax.add_patch(FancyBboxPatch((0.08, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.102, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.14, mini_y, "Subtelomeric α (GC3 > 90%)", fontsize=7.2, fontweight='bold', color='#0369a1', va='center', zorder=7)
            
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.622, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial β (GC3 ~ 65%)", fontsize=7.2, fontweight='bold', color='#15803d', va='center', zorder=7)

        elif atype == "AVIAN":
            ax.add_patch(FancyBboxPatch((0.15, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.172, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.add_patch(FancyBboxPatch((0.55, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.572, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.21, mini_y, "Micro-14 (GC3 ~ 88%)", fontsize=7.2, fontweight='bold', color='#00796b', va='center', zorder=7)
            ax.text(0.61, mini_y, "Micro-25 (GC3 ~ 92%)", fontsize=7.2, fontweight='bold', color='#00796b', va='center', zorder=7)

        elif atype == "SQUAMATE":
            ax.plot([0.06, 0.06], [mini_y - 0.012, mini_y + 0.012], color="#dc2626", lw=2.5, zorder=5)
            ax.add_patch(FancyBboxPatch((0.08, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.102, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.14, mini_y, "Subtelomeric β (GC3 > 94%)", fontsize=7.2, fontweight='bold', color='#15803d', va='center', zorder=7)
            
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.622, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial α (GC3 ~ 50%) — INVERTED POLARITY!", fontsize=7.2, fontweight='bold', color='#b91c1c', va='center', zorder=7)


# =========================================================================
# PANEL C: METRIC ANATOMY, BIOPHYSICS & SYNONYMOUS GC3 CODON ZOOM
# =========================================================================
def draw_panel_c(ax, standalone=False):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    t_size = 14.5 if standalone else 12.5
    s_size = 10.0 if standalone else 8.5

    ax.text(0.02, 0.965, "C. Anatomy of Calculated Metrics: Gene Architecture & Codon Physics",
            fontsize=t_size, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Biophysical basis of CDS GC%, GC3%, Whole-Gene GC%, Intron Compaction, and Direct Synonymous Shifts",
            fontsize=s_size, fontstyle='italic', color='#475569', va='top')

    # Card 1: Gene Architecture & Codon Decomposition
    draw_card(ax, 0.02, 0.625, 0.96, 0.300, title="1. Gene Architecture: Full Gene GC% vs. CDS GC% & Triplet Decomposition",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    # Gene schematic diagram
    ax.plot([0.06, 0.92], [0.850, 0.850], color="#94a3b8", lw=2.2, zorder=2)
    # UTRs
    ax.add_patch(FancyBboxPatch((0.06, 0.835), 0.055, 0.030, boxstyle="square,pad=0", fc="#e2e8f0", ec="#64748b", lw=1, zorder=3))
    ax.text(0.087, 0.850, "5' UTR", ha='center', va='center', fontsize=6.8, color='#475569')
    ax.add_patch(FancyBboxPatch((0.880, 0.835), 0.055, 0.030, boxstyle="square,pad=0", fc="#e2e8f0", ec="#64748b", lw=1, zorder=3))
    ax.text(0.907, 0.850, "3' UTR", ha='center', va='center', fontsize=6.8, color='#475569')

    # Exons
    exons = [
        ("Exon 1 (~95 bp)", 0.125, 0.130, "#1d4ed8"),
        ("Exon 2 (~220 bp)", 0.365, 0.150, "#1d4ed8"),
        ("Exon 3 (~130 bp)", 0.740, 0.130, "#1d4ed8"),
    ]
    for ename, ex, ew, ecol in exons:
        ax.add_patch(FancyBboxPatch((ex, 0.825), ew, 0.050, boxstyle="round,pad=0.002,rounding_size=0.006",
                                   fc=ecol, ec="#1e3a8a", lw=1.1, zorder=4))
        ax.text(ex + ew/2.0, 0.850, ename, ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=5)

    ax.text(0.310, 0.885, "Intron 1 (~120 bp)", ha='center', fontsize=6.8, color='#475569', fontweight='bold')
    ax.text(0.628, 0.885, "Intron 2: Compact in α (140 bp) vs Expanded in β (850 bp)", ha='center', fontsize=7.0, color='#b91c1c', fontweight='bold')

    # Definitions
    ax.text(0.05, 0.795, "Whole-Gene GC%: Exons + Introns + UTRs (Reflects regional chromatin & isochore baseline GC)",
            fontsize=7.8, fontweight='bold', color='#1e293b')
    ax.text(0.05, 0.775, "CDS GC%: Coding sequence only (Concatenated Exons 1, 2, 3; reflects amino acid & mutational bias)",
            fontsize=7.8, fontweight='bold', color='#1d4ed8')

    # Codon Triplet & Selection Corridor (Split on 2 clean lines)
    ax.text(0.05, 0.742, "Codon Triplet Breakdown:", fontsize=8.0, fontweight='bold', color='#0f172a')
    
    # Line 1
    ax.scatter([0.055], [0.720], color="#0284c7", s=18, zorder=4)
    ax.text(0.065, 0.720, "GC1: Moderate purifying constraint (~50–65% GC)", fontsize=6.8, color='#334155', va='center')
    ax.scatter([0.485], [0.720], color="#d97706", s=18, zorder=4)
    ax.text(0.495, 0.720, "GC2: Locked in Purifying Corridor (45–55%)", fontsize=6.8, color='#334155', va='center')

    # Line 2
    ax.scatter([0.055], [0.698], color="#dc2626", s=18, zorder=4)
    ax.text(0.065, 0.698, "GC3: Synonymous site! Barometer of gBGC (35–98% GC)", fontsize=6.8, color='#334155', va='center')
    ax.scatter([0.485], [0.698], color="#7c3aed", s=18, zorder=4)
    ax.text(0.495, 0.698, "GC4d: 4-Fold degenerate (100% synonymous; pure gBGC)", fontsize=6.8, color='#334155', va='center')

    # Mini Purifying Selection Corridor Axis
    corr_y = 0.655
    ax.plot([0.55, 0.92], [corr_y, corr_y], color="#94a3b8", lw=1.2, zorder=4)
    ax.add_patch(patches.Rectangle((0.63, corr_y - 0.007), 0.09, 0.014, facecolor="#bbf7d0", edgecolor="#16a34a", lw=0.8, zorder=5))
    ax.text(0.675, corr_y + 0.012, "Purifying Selection Corridor (45–55%)", ha='center', fontsize=6.5, color='#15803d', fontweight='bold', zorder=6)
    ax.scatter([0.65], [corr_y], color="#d97706", s=20, zorder=7)
    ax.scatter([0.70], [corr_y], color="#0284c7", s=20, zorder=7)
    ax.scatter([0.88], [corr_y], color="#dc2626", s=24, zorder=7)
    ax.text(0.65, corr_y - 0.012, "GC2", ha='center', fontsize=6.2, color='#d97706', fontweight='bold')
    ax.text(0.70, corr_y - 0.012, "GC1", ha='center', fontsize=6.2, color='#0284c7', fontweight='bold')
    ax.text(0.88, corr_y - 0.012, "GC3 ➔ (+gBGC)", ha='center', fontsize=6.2, color='#dc2626', fontweight='bold')

    # Card 2: Conserved Globin Heme Pocket: Proximal Histidine (His-F8) Motif & Synonymous Codon Shift
    draw_card(ax, 0.02, 0.190, 0.96, 0.420,
              title="2. Conserved Globin Heme Pocket: Proximal Histidine (His-F8) Motif & Synonymous GC3 Divergence",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    # Sub-header context
    ax.text(0.04, 0.575, "Comparison across 350+ Myr: 100% Identical Amino Acids, Yet 100% of Divergent 3rd Positions Shifted to G/C in Human!",
            fontsize=7.8, fontweight='bold', color='#0f172a')

    # Data for the 10 conserved amino acids around Proximal His-F8
    codons_human = [
        ('Leu', 'CTG', 'L'),
        ('Ser', 'AGC', 'S'),
        ('Asp', 'GAC', 'D'),
        ('Leu', 'CTG', 'L'),
        ('His', 'CAC', 'H'), # Proximal His-F8
        ('Ala', 'GCG', 'A'),
        ('Arg', 'CGG', 'R'),
        ('Val', 'GTG', 'V'),
        ('Asp', 'GAC', 'D'),
        ('Pro', 'CCG', 'P'),
    ]

    codons_xenopus = [
        ('Leu', 'CTG', 'L'),
        ('Ser', 'AGT', 'S'),
        ('Asp', 'GAC', 'D'),
        ('Leu', 'CTC', 'L'),
        ('His', 'CAT', 'H'), # Proximal His-F8
        ('Ala', 'GCC', 'A'),
        ('Arg', 'AGA', 'R'),
        ('Val', 'GTG', 'V'),
        ('Asp', 'GAT', 'D'),
        ('Pro', 'CCA', 'P'),
    ]

    x_start = 0.190
    x_step = 0.076
    y_human = 0.510
    y_xen = 0.285

    # Human Track Info Card on Left
    h_info = FancyBboxPatch((0.035, y_human - 0.025), 0.135, 0.065, boxstyle="round,pad=0.001,rounding_size=0.004",
                            fc="#eff6ff", ec="#bfdbfe", lw=1.1, zorder=3)
    ax.add_patch(h_info)
    ax.text(0.1025, y_human + 0.022, "Human HBA1 (α)", ha='center', fontsize=7.6, fontweight='bold', color='#1d4ed8', zorder=4)
    ax.text(0.1025, y_human + 0.004, "16p13.3 (Subtelo)", ha='center', fontsize=6.8, color='#2563eb', zorder=4)
    ax.text(0.1025, y_human - 0.014, "Window GC3: 100%", ha='center', fontsize=7.2, fontweight='bold', color='#1e40af', zorder=4)

    # Xenopus Track Info Card on Left
    x_info = FancyBboxPatch((0.035, y_xen - 0.025), 0.135, 0.065, boxstyle="round,pad=0.001,rounding_size=0.004",
                            fc="#f0fdf4", ec="#bbf7d0", lw=1.1, zorder=3)
    ax.add_patch(x_info)
    ax.text(0.1025, y_xen + 0.022, "Xenopus hba1 (α)", ha='center', fontsize=7.6, fontweight='bold', color='#15803d', zorder=4)
    ax.text(0.1025, y_xen + 0.004, "Chr 1 (Interstitial)", ha='center', fontsize=6.8, color='#16a34a', zorder=4)
    ax.text(0.1025, y_xen - 0.014, "Window GC3: 50%", ha='center', fontsize=7.2, fontweight='bold', color='#14532d', zorder=4)

    for i in range(len(codons_human)):
        xh = x_start + i * x_step
        ah, ch, sh = codons_human[i]
        ax_aa, cx, sx = codons_xenopus[i]
        
        is_f8 = (ah == 'His')
        is_shifted = (cx[2] in 'AT' and ch[2] in 'GC')
        
        # Human AA Badge
        aa_bg = '#fef3c7' if is_f8 else '#f1f5f9'
        aa_ec = '#d97706' if is_f8 else '#cbd5e1'
        ax.add_patch(FancyBboxPatch((xh - 0.030, y_human + 0.015), 0.060, 0.022, boxstyle="round,pad=0.001,rounding_size=0.004",
                                   fc=aa_bg, ec=aa_ec, lw=1.1, zorder=3))
        label = f"{ah}-{84+i}" if not is_f8 else "His-F8 (87)"
        ax.text(xh, y_human + 0.026, label, ha='center', va='center', fontsize=6.8, fontweight='bold', color='#0f172a', zorder=4)
        
        # Human Codon Triplets (Pos 1, 2, 3)
        for b_idx in range(2):
            bx = xh - 0.027 + b_idx * 0.018
            ax.add_patch(FancyBboxPatch((bx, y_human - 0.018), 0.016, 0.024, boxstyle="round,pad=0.001,rounding_size=0.002",
                                       fc="#f1f5f9", ec="#94a3b8", lw=0.7, zorder=3))
            ax.text(bx + 0.008, y_human - 0.006, ch[b_idx], ha='center', va='center', fontsize=7.8, fontweight='bold', color='#334155', zorder=4)
        
        # Human Pos 3 (GC3) - Highlighted
        bx3 = xh - 0.027 + 2 * 0.018
        gc3_col = '#1d4ed8' if ch[2] in 'GC' else '#dc2626'
        gc3_bg = '#dbeafe' if ch[2] in 'GC' else '#fee2e2'
        ax.add_patch(FancyBboxPatch((bx3, y_human - 0.018), 0.020, 0.024, boxstyle="round,pad=0.001,rounding_size=0.003",
                                   fc=gc3_bg, ec=gc3_col, lw=1.3, zorder=5))
        ax.text(bx3 + 0.010, y_human - 0.006, ch[2], ha='center', va='center', fontsize=8.8, fontweight='bold', color=gc3_col, zorder=6)
        
        # Xenopus AA Badge
        ax.add_patch(FancyBboxPatch((xh - 0.030, y_xen - 0.030), 0.060, 0.022, boxstyle="round,pad=0.001,rounding_size=0.004",
                                   fc=aa_bg, ec=aa_ec, lw=1.1, zorder=3))
        ax.text(xh, y_xen - 0.019, label, ha='center', va='center', fontsize=6.8, fontweight='bold', color='#0f172a', zorder=4)
        
        # Xenopus Codon Triplets (Pos 1, 2, 3)
        for b_idx in range(2):
            bx = xh - 0.027 + b_idx * 0.018
            ax.add_patch(FancyBboxPatch((bx, y_xen + 0.002), 0.016, 0.024, boxstyle="round,pad=0.001,rounding_size=0.002",
                                       fc="#f1f5f9", ec="#94a3b8", lw=0.7, zorder=3))
            ax.text(bx + 0.008, y_xen + 0.014, cx[b_idx], ha='center', va='center', fontsize=7.8, fontweight='bold', color='#334155', zorder=4)
        
        # Xenopus Pos 3 (GC3)
        bx3_x = xh - 0.027 + 2 * 0.018
        gc3_col_x = '#1d4ed8' if cx[2] in 'GC' else '#dc2626'
        gc3_bg_x = '#dbeafe' if cx[2] in 'GC' else '#fee2e2'
        ax.add_patch(FancyBboxPatch((bx3_x, y_xen + 0.002), 0.020, 0.024, boxstyle="round,pad=0.001,rounding_size=0.003",
                                   fc=gc3_bg_x, ec=gc3_col_x, lw=1.3, zorder=5))
        ax.text(bx3_x + 0.010, y_xen + 0.014, cx[2], ha='center', va='center', fontsize=8.8, fontweight='bold', color=gc3_col_x, zorder=6)

        # Connector between tracks
        mid_y = (y_human - 0.018 + y_xen + 0.026) / 2.0
        if is_shifted:
            ax.plot([bx3 + 0.010, bx3 + 0.010], [y_human - 0.022, y_xen + 0.028], color='#dc2626', lw=1.4, linestyle='--', zorder=4)
            ax.add_patch(FancyBboxPatch((xh - 0.025, mid_y - 0.014), 0.050, 0.028, boxstyle="round,pad=0.001,rounding_size=0.003",
                                       fc="#fef2f2", ec="#ef4444", lw=0.9, zorder=6))
            ax.text(xh, mid_y, f"{cx[2]} ➔ {ch[2]}\n(+gBGC)", ha='center', va='center', fontsize=6.2, fontweight='bold', color='#b91c1c', zorder=7)
        else:
            ax.plot([xh, xh], [y_human - 0.022, y_xen + 0.028], color='#cbd5e1', lw=1.0, zorder=2)
            ax.text(xh, mid_y, "==", ha='center', va='center', fontsize=6.8, color='#94a3b8', zorder=4)

    # Card 3: Bernardi Compaction Law & Distance-to-Telomere Driver
    draw_card(ax, 0.02, 0.015, 0.96, 0.165, title="3. Bernardi's Compaction Law & Telomeric Distance Driver",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    ax.text(0.04, 0.138, "Bernardi's Compaction Law: High-GC isochores suppress intron expansion; low-GC accumulate massive introns.",
            fontsize=7.3, color='#1e293b')
    ax.text(0.04, 0.120, "Human HBA1 (α): 843 bp (Intron 2 = 140 bp) [Compact High-GC]  vs.  Human HBB (β): 1,608 bp (Intron 2 = 850 bp) [Expanded Low-GC]",
            fontsize=7.3, color='#1e40af', fontweight='bold')

    ax.text(0.04, 0.088, "Distance-to-Telomere (dtelo / Lchr) Drives gBGC:", fontsize=7.8, fontweight='bold', color='#0f172a')
    ax.text(0.04, 0.070, "• Meiotic crossovers are physically biased toward chromosome ends, accelerating GC-biased gene conversion (gBGC).",
            fontsize=7.0, color='#334155')
    ax.text(0.04, 0.052, "• Subtelomeric loci (dtelo < 1 Mb) undergo runaway GC3 elevation; internal/centromeric loci remain buffered in AT-rich baselines.",
            fontsize=7.0, color='#b91c1c', fontweight='bold')
    ax.text(0.04, 0.030, "Relative Telomeric Distance: dtelo_rel = min(Pos, Chr_Len - Pos) / Chr_Len    (0.000 Subtelomeric  ➔  0.500 Metacentric)",
            fontsize=6.8, family='monospace', color='#047857', fontweight='bold')


# =========================================================================
# PANEL D: DASHBOARD USER GUIDE & WORKFLOW MANUAL
# =========================================================================
def draw_panel_d(ax, standalone=False):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    t_size = 14.5 if standalone else 12.5
    s_size = 10.0 if standalone else 8.5

    ax.text(0.02, 0.965, "D. Interactive Dashboard Architecture & Analytical User Guide",
            fontsize=t_size, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "How to navigate and interpret the 5 interconnected exploratory modules in dashboard.html",
            fontsize=s_size, fontstyle='italic', color='#475569', va='top')

    modules = [
        ("Tab 1: Global Macroevolutionary Figures",
         "The Big-Picture Evolutionary Synthesis",
         ["• Fig 1A: Curated benchmark summary across 33 representative species.",
          "• Fig 1B: Exhaustive distribution across all 68 vertebrate genomes.",
          "• Fig 2: Mechanistic drivers (Bivariate GC regressions, gene length compaction, telomere distance).",
          "• Fig 3: Time-calibrated phylogenetic chronogram, chromosome synteny maps, and codon ΔGC3 divergence."],
         "#1e40af", "#eff6ff", 0.690, "[ Fig 1A ]  [ Fig 1B ]  [ Fig 2 ]  [ Fig 3 ]"),

        ("Tab 2: 1-on-1 Head-to-Head Comparison",
         "Direct Pairwise Syntenic & Structural Contrast",
         ["• Select Species A & B from searchable comboboxes (e.g., Human vs. Spotted Gar).",
          "• Synchronized views: Linear Chromosome tracks, Bivariate with Mini Inset, or 6-Panel Diagnostic suites.",
          "• Instant differential readout: ΔCDS GC%, ΔWhole-Gene GC%, and chromosomal context."],
         "#065f46", "#ecfdf5", 0.465, "[ Linear Tracks ]  [ Inset Map ]  [ Diagnostic Suite ]"),

        ("Tab 3: Single Species Deep Dive",
         "Locus Anatomy & Architectural Metadata",
         ["• Comprehensive profile of any individual species in the 68-genome database.",
          "• Exact genomic interval bounds, flanking syntenic markers, and gene count.",
          "• High-resolution chromosome ideograms displaying exact physical coordinates."],
         "#7c2d12", "#fff7ed", 0.240, "[ Interval Bounds ]  [ Flanking Markers ]  [ Ideograms ]"),

        ("Tab 4 & 5: Interval Database & Statistical Contrasts",
         "Curated Gene Models & Hypothesis Testing",
         ["• Tab 4 (Database): Filterable, searchable table of all 1,530 curated gene models with z-scores.",
          "• Tab 5 (Statistics): Two-sided Wilcoxon rank-sum tests, Hodges-Lehmann effect sizes, and FDR q-values."],
         "#581c87", "#faf5ff", 0.015, "[ 1,530 Gene Models ]  [ Wilcoxon Tests ]  [ Effect Sizes ]"),
    ]

    if standalone:
        # Standalone layout with taller cards and dedicated workflow footer
        card_h = 0.180
        spacing = 0.020
        start_y = 0.720
        for idx, (mtitle, msub, mpts, mcol, mfc, _, mbtns) in enumerate(modules):
            cy = start_y - idx * (card_h + spacing)
            card = FancyBboxPatch((0.02, cy), 0.96, card_h, boxstyle="round,pad=0.002,rounding_size=0.008",
                                  facecolor=mfc, edgecolor=mcol, lw=1.3, zorder=2)
            ax.add_patch(card)
            ax.text(0.04, cy + card_h - 0.025, mtitle, fontsize=10.5, fontweight='bold', color=mcol, zorder=3)
            ax.text(0.04, cy + card_h - 0.048, msub, fontsize=9.0, fontstyle='italic', color='#475569', zorder=3)
            for p_idx, pt in enumerate(mpts):
                ax.text(0.04, cy + card_h - 0.075 - p_idx * 0.026, pt, fontsize=8.2, color='#1e293b', zorder=3)
            ax.text(0.96, cy + card_h - 0.025, mbtns, ha='right', fontsize=8.2, family='monospace', fontweight='bold', color=mcol, zorder=3)

        # Dedicated Workflow & Shortcuts Footer Card
        foot_box = FancyBboxPatch((0.02, 0.020), 0.96, 0.095, boxstyle="round,pad=0.002,rounding_size=0.008",
                                  facecolor="#f8fafc", edgecolor="#94a3b8", lw=1.2, zorder=2)
        ax.add_patch(foot_box)
        ax.text(0.04, 0.088, "Recommended Exploratory Workflow & Keyboard Shortcuts:", fontsize=9.2, fontweight='bold', color='#0f172a', zorder=3)
        ax.text(0.04, 0.065, "• Step 1 (Tab 1): Inspect macroevolutionary distributions ➔ Step 2 (Tab 2): 1-on-1 syntenic contrasts ➔ Step 3 (Tab 3): Deep dive interval architecture", fontsize=7.8, color='#334155', zorder=3)
        ax.text(0.04, 0.042, "• Interactive Lightbox: Click any image to open fullscreen • Shortcuts: [Esc] Close • [0] Fit to screen • [+] / [-] Zoom • [1]–[6] Diagnostic panels", fontsize=7.8, color='#047857', fontweight='bold', zorder=3)
    else:
        for mtitle, msub, mpts, mcol, mfc, my, mbtns in modules:
            card = FancyBboxPatch((0.02, my), 0.96, 0.210, boxstyle="round,pad=0.002,rounding_size=0.008",
                                  facecolor=mfc, edgecolor=mcol, lw=1.2, zorder=2)
            ax.add_patch(card)
            ax.text(0.04, my + 0.185, mtitle, fontsize=8.8, fontweight='bold', color=mcol, zorder=3)
            ax.text(0.04, my + 0.165, msub, fontsize=7.6, fontstyle='italic', color='#475569', zorder=3)
            for p_idx, pt in enumerate(mpts):
                ax.text(0.04, my + 0.138 - p_idx * 0.024, pt, fontsize=7.2, color='#1e293b', zorder=3)
            ax.text(0.96, my + 0.185, mbtns, ha='right', fontsize=7.0, family='monospace', fontweight='bold', color=mcol, zorder=3)


# =========================================================================
# GENERATORS FOR MASTER AND STANDALONE FIGURES
# =========================================================================
def generate_master_figure():
    fig = plt.figure(figsize=(28, 20), dpi=300)
    fig.patch.set_facecolor('#ffffff')

    fig.text(0.5, 0.985,
             "Figure 0: Conceptual Paradigm, Syntenic Discovery Pipeline, Metric Anatomy & Interactive Dashboard User Guide",
             ha='center', va='top', fontsize=18, fontweight='bold', color='#0f172a')
    fig.text(0.5, 0.968,
             "A Macroevolutionary Atlas of Vertebrate Globin Clusters across 68 Genomes • Resolving Isochore Polarity, Syntenic Translocations, and Codon Biophysics",
             ha='center', va='top', fontsize=12, color='#475569')

    ax_a = fig.add_axes([0.025, 0.505, 0.465, 0.445])
    ax_b = fig.add_axes([0.510, 0.505, 0.465, 0.445])
    ax_c = fig.add_axes([0.025, 0.025, 0.465, 0.455])
    ax_d = fig.add_axes([0.510, 0.025, 0.465, 0.455])

    draw_panel_a(ax_a, standalone=False)
    draw_panel_b(ax_b, standalone=False)
    draw_panel_c(ax_c, standalone=False)
    draw_panel_d(ax_d, standalone=False)

    out_png = os.path.join(RESULTS_DIR, "summary_figure0_overview_and_methodology.png")
    out_pdf = os.path.join(RESULTS_DIR, "summary_figure0_overview_and_methodology.pdf")
    fig.savefig(out_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig.savefig(out_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Master Figure 0 generated: {out_png}")


def generate_standalone_subpanels():
    # Figure 0A
    fig_a, ax_a = plt.subplots(figsize=(18, 12), dpi=300)
    fig_a.patch.set_facecolor('#ffffff')
    draw_panel_a(ax_a, standalone=True)
    out_a_png = os.path.join(RESULTS_DIR, "summary_figure0a_human_paradigm.png")
    out_a_pdf = os.path.join(RESULTS_DIR, "summary_figure0a_human_paradigm.pdf")
    fig_a.savefig(out_a_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_a.savefig(out_a_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_a)
    print(f"✓ Standalone Figure 0A generated: {out_a_png}")

    # Figure 0B
    fig_b, ax_b = plt.subplots(figsize=(18, 12), dpi=300)
    fig_b.patch.set_facecolor('#ffffff')
    draw_panel_b(ax_b, standalone=True)
    out_b_png = os.path.join(RESULTS_DIR, "summary_figure0b_synteny_pipeline.png")
    out_b_pdf = os.path.join(RESULTS_DIR, "summary_figure0b_synteny_pipeline.pdf")
    fig_b.savefig(out_b_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_b.savefig(out_b_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_b)
    print(f"✓ Standalone Figure 0B generated: {out_b_png}")

    # Figure 0C
    fig_c, ax_c = plt.subplots(figsize=(18, 13), dpi=300)
    fig_c.patch.set_facecolor('#ffffff')
    draw_panel_c(ax_c, standalone=True)
    out_c_png = os.path.join(RESULTS_DIR, "summary_figure0c_metric_anatomy_and_gc3.png")
    out_c_pdf = os.path.join(RESULTS_DIR, "summary_figure0c_metric_anatomy_and_gc3.pdf")
    fig_c.savefig(out_c_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_c.savefig(out_c_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_c)
    print(f"✓ Standalone Figure 0C generated: {out_c_png}")

    # Figure 0D
    fig_d, ax_d = plt.subplots(figsize=(18, 12), dpi=300)
    fig_d.patch.set_facecolor('#ffffff')
    draw_panel_d(ax_d, standalone=True)
    out_d_png = os.path.join(RESULTS_DIR, "summary_figure0d_dashboard_guide.png")
    out_d_pdf = os.path.join(RESULTS_DIR, "summary_figure0d_dashboard_guide.pdf")
    fig_d.savefig(out_d_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_d.savefig(out_d_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_d)
    print(f"✓ Standalone Figure 0D generated: {out_d_png}")


def sync_files():
    files_to_sync = [
        "summary_figure0_overview_and_methodology.png",
        "summary_figure0_overview_and_methodology.pdf",
        "summary_figure0_overview_and_methodology.html",
        "overview_guide.html",
        "summary_figure0a_human_paradigm.png",
        "summary_figure0a_human_paradigm.pdf",
        "summary_figure0b_synteny_pipeline.png",
        "summary_figure0b_synteny_pipeline.pdf",
        "summary_figure0c_metric_anatomy_and_gc3.png",
        "summary_figure0c_metric_anatomy_and_gc3.pdf",
        "summary_figure0d_dashboard_guide.png",
        "summary_figure0d_dashboard_guide.pdf",
    ]
    for fn in files_to_sync:
        src = os.path.join(RESULTS_DIR, fn)
        dst = os.path.join(SYNC_DIR, fn)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"  Synchronized: {fn}")


def main():
    print("=== Generating Polished Master Figure 0 and Standalone Subpanel Figures ===")
    generate_master_figure()
    generate_standalone_subpanels()
    
    # Generate interactive overview and methodology guide
    try:
        from generate_overview_guide_html import main as generate_guide
        generate_guide()
    except ImportError:
        import subprocess
        script_path = os.path.join(os.path.dirname(__file__), "generate_overview_guide_html.py")
        if os.path.exists(script_path):
            subprocess.run([sys.executable, script_path], check=True)

    sync_files()
    print("\nAll overview figures and interactive guides successfully created and synchronized!")


if __name__ == "__main__":
    main()
