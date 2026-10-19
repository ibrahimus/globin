#!/usr/bin/env python3
"""
Generate Polished Overview, Methodology & Dashboard User Guide Figures:
- Figure 0: Interactive Dashboard Architecture & Analytical User Guide (summary_figure0_dashboard_guide.png / .pdf)
- Figure 0A: The Human Globin Paradigm & Evolutionary Inquiry (summary_figure0a_human_paradigm.png / .pdf)
- Figure 0B: Synteny-Anchored Locus Retrieval Pipeline & 4 Archetypes (summary_figure0b_synteny_pipeline.png / .pdf)
- Figure 0C: Anatomy of Genomic Metrics, Biophysics & Synonymous GC3 Codon Dynamics (summary_figure0c_metric_anatomy_and_gc3.png / .pdf)
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


def draw_card(ax, x, y, w, h, title="", subtitle="", fc="#ffffff", ec="#cbd5e1", lw=1.2, radius=0.010, zorder=1):
    """Draw a clean, modern card container with an integrated header (no awkward inner box)."""
    box = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.002,rounding_size={radius}",
                         facecolor=fc, edgecolor=ec, lw=lw, zorder=zorder)
    ax.add_patch(box)
    if title:
        ax.text(x + 0.020, y + h - 0.030, title,
                ha='left', va='center', fontsize=11.5, fontweight='bold', color='#0f172a', zorder=zorder + 2)
        if subtitle:
            ax.text(x + 0.020, y + h - 0.052, subtitle,
                    ha='left', va='center', fontsize=8.8, color='#64748b', zorder=zorder + 2)
            div_y = y + h - 0.066
        else:
            div_y = y + h - 0.046
        ax.plot([x + 0.015, x + w - 0.015], [div_y, div_y], color='#e2e8f0', lw=1.0, zorder=zorder + 1)


def draw_gene_arrow(ax, x, y, w, h, name, fc="#0284c7", ec="#0369a1", text_col="white", zorder=4, fontsize=8.5):
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
# FIGURE 0A: THE HUMAN GLOBIN PARADIGM & EVOLUTIONARY INQUIRY
# =========================================================================
def draw_panel_a(ax, standalone=True):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    # Title without subtitle as requested
    ax.text(0.02, 0.970, "Figure 0A: The Human Globin Paradigm & Evolutionary Inquiry",
            fontsize=15.5, fontweight='bold', color='#0f172a', va='top')

    # Card 1: Chr 16 - Alpha Globin Cluster (Subtelomeric Heavy Isochore)
    draw_card(ax, 0.02, 0.520, 0.96, 0.420,
              title="Human Alpha-Globin Cluster — Subtelomeric Heavy Isochore (16p13.3)",
              fc="#fdfefe", ec="#93c5fd", lw=1.4)
    
    # Regional Isochore Badge
    ax.text(0.95, 0.895, "Regional Heavy Isochore H3 (GC > 54%)", ha='right', fontsize=9.2, fontweight='bold', color='#0284c7', zorder=2)

    # Chr 16 Ideogram (NO telomere markings as requested)
    ax.plot([0.06, 0.92], [0.840, 0.840], color="#cbd5e1", lw=12, solid_capstyle='round', zorder=2)
    ax.plot([0.42, 0.42], [0.815, 0.865], color="#475569", lw=3.2, zorder=3)
    ax.text(0.42, 0.875, "Cen (35.3 Mb)", fontsize=8.5, ha='center', color='#64748b', fontweight='bold')
    
    a_x = 0.075
    ax.plot([a_x, a_x], [0.840, 0.770], color="#0284c7", lw=2.2, zorder=5)
    ax.scatter([a_x], [0.840], color="#0284c7", s=45, zorder=6)
    
    # Increased font size for Alpha locus as requested
    ax.text(a_x, 0.755, "Alpha Locus: 0.15–0.22 Mb (0.16% from telomere — Subtelomeric)",
            fontsize=10.5, ha='left', va='top', color='#0369a1', fontweight='bold')

    ax.plot([0.06, 0.92], [0.655, 0.655], color="#cbd5e1", lw=1.8, zorder=2)
    
    anchors_alpha = [
        ("NPRL3", 0.07, 0.075, "#7c3aed"),
        ("MPG", 0.16, 0.060, "#7c3aed"),
        ("RHBDF1", 0.235, 0.075, "#7c3aed"),
        ("LUC7L", 0.81, 0.075, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_alpha:
        draw_gene_arrow(ax, gx, 0.655, gw, 0.052, gname, fc=gcol, ec="#5b21b6", text_col="white", fontsize=8.5)
        ax.text(gx + gw/2.0, 0.612, "Anchor", ha='center', va='top', fontsize=7.8, color='#6b21a8')

    globins_alpha = [
        ("HBZ (ζ)", 0.33, 0.075, "#0284c7"),
        ("HBM (μ)", 0.42, 0.070, "#0284c7"),
        ("HBA2 (α2)", 0.505, 0.075, "#1d4ed8"),
        ("HBA1 (α1)", 0.595, 0.075, "#1d4ed8"),
        ("HBQ1 (θ)", 0.685, 0.075, "#0284c7"),
    ]
    for gname, gx, gw, gcol in globins_alpha:
        draw_gene_arrow(ax, gx, 0.655, gw, 0.058, gname, fc=gcol, ec="#1e3a8a", text_col="white", fontsize=8.5)

    m_box_a = FancyBboxPatch((0.06, 0.540), 0.86, 0.046, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#eff6ff", edgecolor="#bfdbfe", lw=1.1, zorder=3)
    ax.add_patch(m_box_a)
    ax.text(0.49, 0.563,
            "α-Metrics:  CDS GC = 66.2%  |  GC3 = 91.1%  |  Gene GC = 68.2%  |  Mean Gene Length = 810 bp  |  Isochore: Heavy H3",
            ha='center', va='center', fontsize=9.2, fontweight='bold', color='#1e40af')

    # Card 2: Chr 11 - Beta Globin Cluster (Interstitial Light Isochore)
    draw_card(ax, 0.02, 0.100, 0.96, 0.400,
              title="Human Beta-Globin Cluster — Interstitial Light Isochore (11p15.4)",
              fc="#fdfefe", ec="#86efac", lw=1.4)
    
    # Regional Isochore Badge
    ax.text(0.95, 0.455, "Regional Light Isochore L1/L2 (GC ~40%)", ha='right', fontsize=9.2, fontweight='bold', color='#16a34a', zorder=2)

    # Chr 11 Ideogram (NO telomere markings as requested)
    ax.plot([0.06, 0.92], [0.405, 0.405], color="#cbd5e1", lw=12, solid_capstyle='round', zorder=2)
    ax.plot([0.45, 0.45], [0.380, 0.430], color="#475569", lw=3.2, zorder=3)
    ax.text(0.45, 0.440, "Cen (51.6 Mb)", fontsize=8.5, ha='center', color='#64748b', fontweight='bold')
    
    b_x = 0.20
    ax.plot([b_x, b_x], [0.405, 0.335], color="#16a34a", lw=2.2, zorder=5)
    ax.scatter([b_x], [0.405], color="#16a34a", s=45, zorder=6)
    
    # Increased font size for Beta locus as requested
    ax.text(b_x, 0.320, "Beta Locus: 5.22–5.27 Mb (3.85% from telomere — Interstitial)",
            fontsize=10.5, ha='left', va='top', color='#15803d', fontweight='bold')

    ax.plot([0.06, 0.92], [0.220, 0.220], color="#cbd5e1", lw=1.8, zorder=2)

    anchors_beta = [
        ("OR51V1", 0.065, 0.070, "#7c3aed"),
        ("OR51B4", 0.145, 0.070, "#7c3aed"),
        ("TRIM21", 0.770, 0.070, "#7c3aed"),
        ("RRM1", 0.850, 0.070, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_beta:
        draw_gene_arrow(ax, gx, 0.220, gw, 0.052, gname, fc=gcol, ec="#5b21b6", text_col="white", fontsize=8.5)
        ax.text(gx + gw/2.0, 0.177, "Anchor", ha='center', va='top', fontsize=7.8, color='#6b21a8')

    globins_beta = [
        ("HBE1 (ε)", 0.245, 0.075, "#15803d"),
        ("HBG2 (Gγ)", 0.335, 0.075, "#15803d"),
        ("HBG1 (Aγ)", 0.425, 0.075, "#15803d"),
        ("HBD (δ)", 0.535, 0.070, "#16a34a"),
        ("HBB (β)", 0.620, 0.075, "#16a34a"),
    ]
    for gname, gx, gw, gcol in globins_beta:
        draw_gene_arrow(ax, gx, 0.220, gw, 0.058, gname, fc=gcol, ec="#14532d", text_col="white", fontsize=8.5)

    m_box_b = FancyBboxPatch((0.06, 0.120), 0.86, 0.046, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#f0fdf4", edgecolor="#bbf7d0", lw=1.1, zorder=3)
    ax.add_patch(m_box_b)
    ax.text(0.49, 0.143,
            "β-Metrics:  CDS GC = 54.4%  |  GC3 = 65.8%  |  Gene GC = 43.4%  |  Mean Gene Length = 1,608 bp  |  Isochore: Light L1",
            ha='center', va='center', fontsize=9.2, fontweight='bold', color='#166534')

    # Core Evolutionary Inquiry Callout Box
    q_box = FancyBboxPatch((0.02, 0.010), 0.96, 0.075, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f8fafc", edgecolor="#94a3b8", lw=1.2, zorder=3)
    ax.add_patch(q_box)
    ax.text(0.04, 0.062, "Core Macroevolutionary Inquiries:", fontsize=9.8, fontweight='bold', color='#0f172a', zorder=4)
    ax.text(0.04, 0.040, "1. Is this extreme GC disparity (subtelomeric high-GC α vs. interstitial low-GC β) an invariant vertebrate feature, or a mammalian derived state?",
            fontsize=8.8, color='#334155', zorder=4)
    ax.text(0.04, 0.020, "2. Did physical distance to the telomere dictate locus GC content via GC-biased gene conversion (gBGC) across 500 Myr of vertebrate evolution?",
            fontsize=8.8, color='#334155', zorder=4)


# =========================================================================
# FIGURE 0B: SYNTENIC RETRIEVAL PIPELINE & 4 ARCHETYPES
# =========================================================================
def draw_panel_b(ax, standalone=True):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title without subtitle as requested
    ax.text(0.02, 0.970, "Figure 0B: Synteny-Anchored Locus Retrieval Pipeline & 4 Archetypes",
            fontsize=15.5, fontweight='bold', color='#0f172a', va='top')

    # Pipeline Box: 5 Comprehensive Stages including HMMER proteome scanning and de novo modeling
    draw_card(ax, 0.02, 0.740, 0.96, 0.200, title="Synteny Discovery, HMMER Scanning & Chromosome Validation Pipeline",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)
    
    # 5 Detailed Pipeline Steps with increased font size
    steps = [
        ("1. Syntenic Flank Search",
         "Query 68 genomes for anchor markers:\nα: NPRL3, MPG, RHBDF1, LUC7L\nβ: OR51 cluster, TRIM21, RRM1",
         "#2563eb"),
        ("2. HMMER / TBLASTN Scan",
         "When annotations missing/divergent:\nProfile HMM scan (Pfam PF00042)\nCross-species seeded TBLASTN",
         "#7c3aed"),
        ("3. De Novo Gene Models",
         "Extract canonical 3-exon / 2-intron\nglobin coding sequences (GT-AG);\nvalidate invariant heme pocket",
         "#059669"),
        ("4. Scaffold vs. Chromosome",
         "Discriminate fully assembled arms\nfrom unplaced contigs to guarantee\ntrue telomere distance measurement",
         "#d97706"),
        ("5. Baseline Calibration",
         "Compute whole-genome CDS baselines\n(N = 15k–35k genes) to quantify\nisochore z-score deviations (ΔGC3)",
         "#0d9488")
    ]
    for idx, (stitle, sdesc, scol) in enumerate(steps):
        sx = 0.035 + idx * 0.188
        sbox = FancyBboxPatch((sx, 0.755), 0.178, 0.130, boxstyle="round,pad=0.002,rounding_size=0.006",
                              facecolor="#f8fafc", edgecolor=scol, lw=1.4, zorder=3)
        ax.add_patch(sbox)
        # Increased font size as requested
        ax.text(sx + 0.010, 0.865, stitle, fontsize=9.6, fontweight='bold', color=scol, zorder=4)
        ax.text(sx + 0.010, 0.840, sdesc, fontsize=8.2, color='#334155', va='top', zorder=4, linespacing=1.25)

    # Archetypes Card (No exclamation marks, calm scientific tone)
    draw_card(ax, 0.02, 0.015, 0.96, 0.705, title="The 4 Core Macroevolutionary Archetypes Across Vertebrate Evolution",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    archetypes = [
        ("1. Linked Primordial Archetype (Basal Actinopterygii & Amphibians — Spotted Gar, Xenopus)",
         "Ancestral tandem arrangement: α and β globin genes reside on the same chromosome in a single cluster.",
         "Both genes share identical genomic environment and regional recombination rate. ΔGC3 (α - β) ≈ 0%.",
         "#1565c0", "#f0fdf4", 0.505, "Single Chromosome (Gar Chr 19 / Xenopus Chr 1)", "LINKED"),

        ("2. Canonical Decoupled Subtelomeric-α (Eutherian & Metatherian Mammals — Human, Mouse)",
         "Amniote translocation decouples α and β. α captured by subtelomere (<1 Mb to telomere); β remains interstitial.",
         "Runaway subtelomeric gBGC supercharges α-GC3 (>85–95%) while β-GC3 remains moderate (~60–65%). ΔGC3 > 0%.",
         "#e65100", "#fff7ed", 0.345, "Two Chromosomes (α Subtelomeric / β Interstitial)", "MAMMALIAN"),

        ("3. Dual Microchromosomal Partitioning (Archosauria & Aves — Chicken, Eagle, Hummingbird)",
         "Both α and β clusters translocate onto distinct, highly recombinogenic microchromosomes (Chr 14 and Chr 25).",
         "High crossover density drives synchronous hyper-elevation of both clusters: α-GC3 > 85%, β-GC3 > 85%. ΔGC3 ≈ 0%.",
         "#00796b", "#f0fdfa", 0.185, "Two Microchromosomes (Micro-14 and Micro-25)", "AVIAN"),

        ("4. Squamate Inverted Polarity (Lizards & Snakes — Wall Lizard, Garter Snake)",
         "Squamates exhibit inverted chromosomal polarity: the β locus is translocated to a microchromosome or subtelomere,",
         "driving elevated β-GC3 (85–95%), while the α locus resides in an interstitial genomic desert (GC3 48–55%). ΔGC3 < 0%.",
         "#b91c1c", "#fef2f2", 0.025, "Inverted Polarity (β Subtelomeric / α Interstitial)", "SQUAMATE"),
    ]

    for atitle, adesc1, adesc2, acol, afc, ay, achr, atype in archetypes:
        abox = FancyBboxPatch((0.035, ay), 0.93, 0.138, boxstyle="round,pad=0.002,rounding_size=0.008",
                              facecolor=afc, edgecolor=acol, lw=1.2, zorder=3)
        ax.add_patch(abox)
        ax.text(0.055, ay + 0.123, atitle, fontsize=9.4, fontweight='bold', color=acol, va='top', zorder=4)
        ax.text(0.055, ay + 0.103, adesc1, fontsize=8.3, color='#334155', va='top', zorder=4)
        ax.text(0.055, ay + 0.083, adesc2, fontsize=8.3, fontweight='bold', color='#1e293b', va='top', zorder=4)
        ax.text(0.055, ay + 0.063, f"Cytogenetic State: {achr}", fontsize=7.8, color='#64748b', va='top', zorder=4)

        mini_y = ay + 0.026
        ax.plot([0.055, 0.91], [mini_y, mini_y], color="#cbd5e1", lw=7, solid_capstyle='round', zorder=4)
        
        if atype == "LINKED":
            ax.plot([0.42, 0.42], [mini_y - 0.010, mini_y + 0.010], color="#475569", lw=2, zorder=5)
            ax.add_patch(FancyBboxPatch((0.70, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.add_patch(FancyBboxPatch((0.745, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.722, mini_y, "α", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.767, mini_y, "β", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.81, mini_y, "◀ Linked α-β (ΔGC3 ≈ 0)", fontsize=8.0, fontweight='bold', color='#1565c0', va='center', zorder=7)
            
        elif atype == "MAMMALIAN":
            ax.plot([0.055, 0.055], [mini_y - 0.010, mini_y + 0.010], color="#dc2626", lw=2.5, zorder=5)
            ax.add_patch(FancyBboxPatch((0.075, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.097, mini_y, "α", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.135, mini_y, "Subtelomeric α (GC3 > 90%)", fontsize=8.0, fontweight='bold', color='#0369a1', va='center', zorder=7)
            
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.622, mini_y, "β", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial β (GC3 ~ 65%)", fontsize=8.0, fontweight='bold', color='#15803d', va='center', zorder=7)

        elif atype == "AVIAN":
            ax.add_patch(FancyBboxPatch((0.15, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.172, mini_y, "α", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.add_patch(FancyBboxPatch((0.55, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.572, mini_y, "β", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.21, mini_y, "Micro-14 (GC3 ~ 88%)", fontsize=8.0, fontweight='bold', color='#00796b', va='center', zorder=7)
            ax.text(0.61, mini_y, "Micro-25 (GC3 ~ 92%)", fontsize=8.0, fontweight='bold', color='#00796b', va='center', zorder=7)

        elif atype == "SQUAMATE":
            ax.plot([0.055, 0.055], [mini_y - 0.010, mini_y + 0.010], color="#dc2626", lw=2.5, zorder=5)
            ax.add_patch(FancyBboxPatch((0.075, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.097, mini_y, "β", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.135, mini_y, "Subtelomeric β (GC3 > 94%)", fontsize=8.0, fontweight='bold', color='#15803d', va='center', zorder=7)
            
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.014), 0.045, 0.028, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.622, mini_y, "α", ha='center', va='center', fontsize=7.5, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial α (GC3 ~ 50%) — Inverted Polarity", fontsize=8.0, fontweight='bold', color='#b91c1c', va='center', zorder=7)
            


# =========================================================================
# FIGURE 0C: METRIC ANATOMY, BIOPHYSICS & SYNONYMOUS GC3 CODON DYNAMICS
# =========================================================================
def draw_panel_c(ax, standalone=True):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title with "Codon Dynamics" as requested, ample space above Card 1
    ax.text(0.02, 0.978, "Figure 0C: Anatomy of Calculated Metrics: Gene Architecture & Codon Dynamics",
            fontsize=15.5, fontweight='bold', color='#0f172a', va='top')

    # Card 1: Gene Architecture & Splicing Demarcation (starts at y=0.555, height 0.370 -> top at 0.925)
    draw_card(ax, 0.02, 0.555, 0.96, 0.370,
              title="1. Gene Architecture: Whole-Gene GC% vs. CDS GC% Splicing Demarcation",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    # Demarcation Line 1: Whole-Gene GC% Range (Exons 1–3 + Introns 1–2, strictly EXCLUDING UTRs)
    y_g_dem = 0.852
    ax.plot([0.220, 0.860], [y_g_dem, y_g_dem], color='#047857', lw=1.8, marker='|', markersize=8, zorder=6)
    ax.text(0.540, y_g_dem + 0.012, "◀── Whole-Gene GC% Range: Exons 1–3 + Introns 1–2 (UTRs Excluded to Avoid Aberrant GC Skew: ~0.84 kb in α vs ~1.61 kb in β) ──▶",
            ha='center', fontsize=8.2, fontweight='bold', color='#047857')

    # Genomic Track (Upper)
    y_gen = 0.812
    ax.text(0.035, y_gen, "Genomic DNA Locus:", fontsize=9.0, fontweight='bold', color='#0f172a', va='center')
    ax.plot([0.140, 0.940], [y_gen, y_gen], color="#94a3b8", lw=2.5, zorder=2)
    
    # UTRs (Explicitly Excluded from Whole-Gene GC% to avoid strange GC/CpG skew)
    ax.add_patch(FancyBboxPatch((0.145, y_gen - 0.016), 0.060, 0.032, boxstyle="square,pad=0",
                               fc="#f1f5f9", ec="#94a3b8", lw=1.2, hatch="//", zorder=3))
    ax.text(0.175, y_gen, "5' UTR", ha='center', va='center', fontsize=7.5, color='#475569', fontweight='bold', zorder=4)
    ax.text(0.175, y_gen + 0.024, "Excluded", ha='center', fontsize=7.0, color='#dc2626', fontweight='bold', zorder=5)

    ax.add_patch(FancyBboxPatch((0.875, y_gen - 0.016), 0.060, 0.032, boxstyle="square,pad=0",
                               fc="#f1f5f9", ec="#94a3b8", lw=1.2, hatch="//", zorder=3))
    ax.text(0.905, y_gen, "3' UTR", ha='center', va='center', fontsize=7.5, color='#475569', fontweight='bold', zorder=4)
    ax.text(0.905, y_gen + 0.024, "Excluded", ha='center', fontsize=7.0, color='#dc2626', fontweight='bold', zorder=5)

    # Exons on Genomic Track
    exons_gen = [
        ("Exon 1 (~95 bp)", 0.220, 0.110, "#1d4ed8"),
        ("Exon 2 (~220 bp)", 0.430, 0.140, "#1d4ed8"),
        ("Exon 3 (~130 bp)", 0.740, 0.120, "#1d4ed8"),
    ]
    for ename, ex, ew, ecol in exons_gen:
        ax.add_patch(FancyBboxPatch((ex, y_gen - 0.020), ew, 0.040, boxstyle="round,pad=0.002,rounding_size=0.005",
                                   fc=ecol, ec="#1e3a8a", lw=1.2, zorder=4))
        ax.text(ex + ew/2.0, y_gen, ename, ha='center', va='center', fontsize=8.0, fontweight='bold', color='white', zorder=5)

    ax.text(0.380, y_gen + 0.026, "Intron 1 (~120 bp)", ha='center', fontsize=7.4, color='#475569', fontweight='bold')
    ax.text(0.655, y_gen + 0.026, "Intron 2 (140 bp in α vs 850 bp in β)", ha='center', fontsize=7.6, color='#b91c1c', fontweight='bold')

    # Splicing Connectors: Pulling exons down to form merged CDS
    y_cds = 0.722
    ax.text(0.035, y_cds, "Spliced mRNA (CDS):", fontsize=9.0, fontweight='bold', color='#1d4ed8', va='center')
    
    # Spliced Exons concatenated together
    exons_spliced = [
        ("Exon 1", 0.260, 0.110, "#1d4ed8"),
        ("Exon 2", 0.370, 0.180, "#1d4ed8"),
        ("Exon 3", 0.550, 0.140, "#1d4ed8"),
    ]
    for ename, ex, ew, ecol in exons_spliced:
        ax.add_patch(FancyBboxPatch((ex, y_cds - 0.018), ew, 0.036, boxstyle="round,pad=0.001,rounding_size=0.004",
                                   fc=ecol, ec="#1e3a8a", lw=1.2, zorder=4))
        ax.text(ex + ew/2.0, y_cds, ename, ha='center', va='center', fontsize=8.2, fontweight='bold', color='white', zorder=5)

    # Demarcation Line 2: CDS GC% Range (Merged coding sequence only)
    y_c_dem = y_cds - 0.028
    ax.plot([0.260, 0.690], [y_c_dem, y_c_dem], color='#1d4ed8', lw=1.8, marker='|', markersize=8, zorder=6)
    ax.text(0.475, y_c_dem - 0.013, "◀── CDS GC% Calculation Range (Concatenated Spliced Coding Sequence Only: 426–444 bp) ──▶",
            ha='center', fontsize=8.4, fontweight='bold', color='#1d4ed8')

    # Splicing dashed guideline paths
    ax.plot([0.220 + 0.055, 0.260 + 0.055], [y_gen - 0.022, y_cds + 0.020], color='#3b82f6', lw=1.1, linestyle='--', zorder=2)
    ax.plot([0.430 + 0.070, 0.370 + 0.090], [y_gen - 0.022, y_cds + 0.020], color='#3b82f6', lw=1.1, linestyle='--', zorder=2)
    ax.plot([0.740 + 0.060, 0.550 + 0.070], [y_gen - 0.022, y_cds + 0.020], color='#3b82f6', lw=1.1, linestyle='--', zorder=2)
    ax.text(0.760, y_cds, "Splicing removes introns;\nCDS feeds codon decomposition", fontsize=7.6, color='#475569', va='center')

    # Codon Position Breakdown & Purifying Corridor
    div_c1 = 0.655
    ax.plot([0.035, 0.945], [div_c1, div_c1], color='#e2e8f0', lw=1.0, zorder=2)
    
    ax.text(0.035, 0.636, "Codon Position Decomposition & Selection Constraints:", fontsize=9.4, fontweight='bold', color='#0f172a')
    
    ax.scatter([0.045], [0.612], color="#0284c7", s=32, zorder=10)
    ax.text(0.058, 0.612, "GC1: 1st position (~50–65% GC) — Constrained by amino acid polarity", fontsize=8.5, color='#1e293b', va='center', fontweight='bold')
    
    ax.scatter([0.045], [0.590], color="#d97706", s=32, zorder=10)
    ax.text(0.058, 0.590, "GC2: 2nd position (45–55% GC) — Strictly locked in Purifying Selection Corridor across all 68 vertebrates", fontsize=8.5, color='#1e293b', va='center', fontweight='bold')

    ax.scatter([0.045], [0.568], color="#dc2626", s=32, zorder=10)
    ax.text(0.058, 0.568, "GC3: 3rd position (35–98% GC) — Synonymous site! Neutral barometer of GC-biased gene conversion (gBGC)", fontsize=8.5, color='#b91c1c', va='center', fontweight='bold')

    # Purifying Selection Corridor Axis (Right side)
    corr_y = 0.602
    ax.plot([0.62, 0.94], [corr_y, corr_y], color="#94a3b8", lw=1.8, zorder=4)
    ax.add_patch(patches.Rectangle((0.68, corr_y - 0.015), 0.10, 0.030, facecolor="#bbf7d0", edgecolor="#16a34a", lw=1.2, zorder=5))
    ax.text(0.730, corr_y + 0.022, "Purifying Selection Corridor (45–55%)", ha='center', fontsize=8.2, color='#15803d', fontweight='bold', zorder=6)
    
    ax.scatter([0.70], [corr_y], color="#d97706", s=45, zorder=10)
    ax.scatter([0.75], [corr_y], color="#0284c7", s=45, zorder=10)
    ax.scatter([0.91], [corr_y], color="#dc2626", s=50, zorder=10)
    ax.text(0.70, corr_y - 0.018, "GC2", ha='center', fontsize=8.2, color='#d97706', fontweight='bold', zorder=11)
    ax.text(0.75, corr_y - 0.018, "GC1", ha='center', fontsize=8.2, color='#0284c7', fontweight='bold', zorder=11)
    ax.text(0.91, corr_y - 0.018, "GC3 (+gBGC)", ha='center', fontsize=8.2, color='#dc2626', fontweight='bold', zorder=11)

    # Card 2: Conserved Globin Heme Pocket: Proximal Histidine (His-F8) Motif
    draw_card(ax, 0.02, 0.175, 0.96, 0.365,
              title="2. Conserved Globin Heme Pocket: Proximal Histidine (His-F8) Motif & Synonymous GC3 Divergence",
              subtitle="Comparison across 350+ Myr: 100% Identical Amino Acids, Yet 100% of Divergent 3rd Positions Shifted to G/C in Human",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

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
    y_human = 0.435
    y_xen = 0.245

    # Human Track Info Card on Left
    h_info = FancyBboxPatch((0.035, y_human - 0.025), 0.135, 0.065, boxstyle="round,pad=0.001,rounding_size=0.004",
                            fc="#eff6ff", ec="#bfdbfe", lw=1.2, zorder=3)
    ax.add_patch(h_info)
    ax.text(0.1025, y_human + 0.022, "Human HBA1 (α)", ha='center', fontsize=8.6, fontweight='bold', color='#1d4ed8', zorder=4)
    ax.text(0.1025, y_human + 0.004, "16p13.3 (Subtelo)", ha='center', fontsize=7.6, color='#2563eb', zorder=4)
    ax.text(0.1025, y_human - 0.014, "Window GC3: 100%", ha='center', fontsize=8.2, fontweight='bold', color='#1e40af', zorder=4)

    # Xenopus Track Info Card on Left
    x_info = FancyBboxPatch((0.035, y_xen - 0.025), 0.135, 0.065, boxstyle="round,pad=0.001,rounding_size=0.004",
                            fc="#f0fdf4", ec="#bbf7d0", lw=1.2, zorder=3)
    ax.add_patch(x_info)
    ax.text(0.1025, y_xen + 0.022, "Xenopus hba1 (α)", ha='center', fontsize=8.6, fontweight='bold', color='#15803d', zorder=4)
    ax.text(0.1025, y_xen + 0.004, "Chr 1 (Interstitial)", ha='center', fontsize=7.6, color='#16a34a', zorder=4)
    ax.text(0.1025, y_xen - 0.014, "Window GC3: 50%", ha='center', fontsize=8.2, fontweight='bold', color='#14532d', zorder=4)

    for i in range(len(codons_human)):
        xh = x_start + i * x_step
        ah, ch, sh = codons_human[i]
        ax_aa, cx, sx = codons_xenopus[i]
        
        is_f8 = (ah == 'His')
        is_shifted = (cx[2] in 'AT' and ch[2] in 'GC')
        
        # Human AA Badge
        aa_bg = '#fef3c7' if is_f8 else '#f1f5f9'
        aa_ec = '#d97706' if is_f8 else '#cbd5e1'
        ax.add_patch(FancyBboxPatch((xh - 0.030, y_human + 0.015), 0.060, 0.024, boxstyle="round,pad=0.001,rounding_size=0.004",
                                   fc=aa_bg, ec=aa_ec, lw=1.2, zorder=3))
        label = f"{ah}-{84+i}" if not is_f8 else "His-F8 (87)"
        ax.text(xh, y_human + 0.027, label, ha='center', va='center', fontsize=7.6, fontweight='bold', color='#0f172a', zorder=4)
        
        # Human Codon Triplets (Pos 1, 2, 3)
        for b_idx in range(2):
            bx = xh - 0.027 + b_idx * 0.018
            ax.add_patch(FancyBboxPatch((bx, y_human - 0.018), 0.016, 0.026, boxstyle="round,pad=0.001,rounding_size=0.002",
                                       fc="#f1f5f9", ec="#94a3b8", lw=0.8, zorder=3))
            ax.text(bx + 0.008, y_human - 0.005, ch[b_idx], ha='center', va='center', fontsize=8.5, fontweight='bold', color='#334155', zorder=4)
        
        # Human Pos 3 (GC3) - Highlighted
        bx3 = xh - 0.027 + 2 * 0.018
        gc3_col = '#1d4ed8' if ch[2] in 'GC' else '#dc2626'
        gc3_bg = '#dbeafe' if ch[2] in 'GC' else '#fee2e2'
        ax.add_patch(FancyBboxPatch((bx3, y_human - 0.018), 0.020, 0.026, boxstyle="round,pad=0.001,rounding_size=0.003",
                                   fc=gc3_bg, ec=gc3_col, lw=1.4, zorder=5))
        ax.text(bx3 + 0.010, y_human - 0.005, ch[2], ha='center', va='center', fontsize=9.2, fontweight='bold', color=gc3_col, zorder=6)
        
        # Xenopus AA Badge
        ax.add_patch(FancyBboxPatch((xh - 0.030, y_xen - 0.032), 0.060, 0.024, boxstyle="round,pad=0.001,rounding_size=0.004",
                                   fc=aa_bg, ec=aa_ec, lw=1.2, zorder=3))
        ax.text(xh, y_xen - 0.020, label, ha='center', va='center', fontsize=7.6, fontweight='bold', color='#0f172a', zorder=4)
        
        # Xenopus Codon Triplets (Pos 1, 2, 3)
        for b_idx in range(2):
            bx = xh - 0.027 + b_idx * 0.018
            ax.add_patch(FancyBboxPatch((bx, y_xen + 0.002), 0.016, 0.026, boxstyle="round,pad=0.001,rounding_size=0.002",
                                       fc="#f1f5f9", ec="#94a3b8", lw=0.8, zorder=3))
            ax.text(bx + 0.008, y_xen + 0.015, cx[b_idx], ha='center', va='center', fontsize=8.5, fontweight='bold', color='#334155', zorder=4)
        
        # Xenopus Pos 3 (GC3)
        bx3_x = xh - 0.027 + 2 * 0.018
        gc3_col_x = '#1d4ed8' if cx[2] in 'GC' else '#dc2626'
        gc3_bg_x = '#dbeafe' if cx[2] in 'GC' else '#fee2e2'
        ax.add_patch(FancyBboxPatch((bx3_x, y_xen + 0.002), 0.020, 0.026, boxstyle="round,pad=0.001,rounding_size=0.003",
                                   fc=gc3_bg_x, ec=gc3_col_x, lw=1.4, zorder=5))
        ax.text(bx3_x + 0.010, y_xen + 0.015, cx[2], ha='center', va='center', fontsize=9.2, fontweight='bold', color=gc3_col_x, zorder=6)

        # Connector between tracks
        mid_y = (y_human - 0.018 + y_xen + 0.028) / 2.0
        if is_shifted:
            ax.plot([bx3 + 0.010, bx3 + 0.010], [y_human - 0.022, y_xen + 0.030], color='#dc2626', lw=1.6, linestyle='--', zorder=4)
            ax.add_patch(FancyBboxPatch((xh - 0.025, mid_y - 0.015), 0.050, 0.030, boxstyle="round,pad=0.001,rounding_size=0.003",
                                       fc="#fef2f2", ec="#ef4444", lw=1.0, zorder=6))
            ax.text(xh, mid_y, f"{cx[2]} ➔ {ch[2]}\n(+gBGC)", ha='center', va='center', fontsize=7.2, fontweight='bold', color='#b91c1c', zorder=7)
        else:
            ax.plot([xh, xh], [y_human - 0.022, y_xen + 0.030], color='#cbd5e1', lw=1.2, zorder=2)
            ax.text(xh, mid_y, "==", ha='center', va='center', fontsize=7.5, color='#94a3b8', zorder=4)

    # Card 3: Bernardi Compaction Law & Distance-to-Telomere Driver
    draw_card(ax, 0.02, 0.008, 0.96, 0.156,
              title="3. Bernardi's Compaction Law & Telomeric Distance Driver",
              subtitle="High-GC isochores suppress intron expansion; low-GC regions accumulate expanded introns.",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    ax.text(0.04, 0.076, "Human HBA1 (α): 843 bp (Intron 2 = 140 bp) [Compact High-GC]  vs.  Human HBB (β): 1,608 bp (Intron 2 = 850 bp) [Expanded Low-GC]",
            fontsize=8.4, color='#1e40af', fontweight='bold')

    ax.text(0.04, 0.056, "Distance-to-Telomere (dtelo / Lchr) Drives gBGC: Subtelomeric loci (<1 Mb) undergo runaway GC3 elevation; internal loci remain AT-buffered.",
            fontsize=8.2, color='#0f172a', fontweight='bold')
    ax.text(0.04, 0.036, "• Meiotic crossovers are physically biased toward chromosome ends, accelerating GC-biased gene conversion (gBGC).",
            fontsize=8.0, color='#334155')
    ax.text(0.04, 0.017, "Relative Telomeric Distance: dtelo_rel = min(Pos, Chr_Len - Pos) / Chr_Len    (0.000 Subtelomeric  ➔  0.500 Metacentric)",
            fontsize=7.8, family='monospace', color='#047857', fontweight='bold')


# =========================================================================
# FIGURE 0: DASHBOARD USER GUIDE & WORKFLOW MANUAL (FIRST FIGURE)
# =========================================================================
def draw_panel_d(ax, standalone=True):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    ax.text(0.02, 0.970, "Figure 0: Interactive Dashboard Architecture & Analytical User Guide",
            fontsize=15.5, fontweight='bold', color='#0f172a', va='top')

    modules = [
        ("Tab 1: Global Macroevolutionary Figures",
         "The Big-Picture Evolutionary Synthesis across 68 Vertebrates",
         ["• Fig 0, 0A–0C: Conceptual framework, syntenic retrieval pipeline, metric anatomy, and codon dynamics.",
          "• Fig 1A & 1B: Curated summary across 33 species (1A) and exhaustive 68-vertebrate macroevolutionary atlas (1B).",
          "• Fig 2 & 3: Mechanistic drivers of GC polarity (Fig 2) and time-calibrated phylogenetic synteny chronogram (Fig 3)."],
         "#1e40af", "#eff6ff", "[ Fig 0 ]  [ Fig 1A ]  [ Fig 1B ]  [ Fig 2 ]  [ Fig 3 ]"),

        ("Tab 2: 1-on-1 Head-to-Head Comparison",
         "Direct Pairwise Syntenic & Structural Contrast Between Any Two Species",
         ["• Select Species A & B from searchable autocomplete comboboxes (e.g., Human vs. Spotted Gar).",
          "• Synchronized views: Linear Chromosome tracks, Bivariate with Mini Inset, or 6-Panel Diagnostic suites.",
          "• Instant differential readout: ΔCDS GC%, ΔWhole-Gene GC%, chromosome coordinates, and syntenic context."],
         "#065f46", "#ecfdf5", "[ Linear Tracks ]  [ Inset Map ]  [ Diagnostic Suite ]"),

        ("Tab 3: Single Species Deep Dive",
         "Locus Anatomy & Architectural Metadata for Any Individual Genome",
         ["• Comprehensive profile of any individual species in the 68-genome database.",
          "• Exact genomic interval bounds, flanking syntenic markers, and curated gene count.",
          "• High-resolution chromosome ideograms displaying exact physical coordinates and centromere locations."],
         "#7c2d12", "#fff7ed", "[ Interval Bounds ]  [ Flanking Markers ]  [ Ideograms ]"),

        ("Tabs 4 & 5: Interval Database & Statistical Contrasts",
         "Curated Gene Models & Hypothesis Testing Suite",
         ["• Tab 4 (Database): Filterable, searchable table of all 1,180 verified globin genes with z-scores and coordinates.",
          "• Tab 5 (Statistics): Two-sided Wilcoxon signed-rank and Mann-Whitney U tests with Hodges-Lehmann effect sizes."],
         "#581c87", "#faf5ff", "[ 1,180 Gene Models ]  [ Wilcoxon Tests ]  [ Effect Sizes ]"),
    ]

    card_h = 0.170
    spacing = 0.020
    start_y = 0.735
    for idx, (mtitle, msub, mpts, mcol, mfc, mbtns) in enumerate(modules):
        cy = start_y - idx * (card_h + spacing)
        card = FancyBboxPatch((0.02, cy), 0.96, card_h, boxstyle="round,pad=0.002,rounding_size=0.008",
                              facecolor=mfc, edgecolor=mcol, lw=1.3, zorder=2)
        ax.add_patch(card)
        ax.text(0.04, cy + card_h - 0.028, mtitle, fontsize=11.2, fontweight='bold', color=mcol, zorder=3)
        ax.text(0.04, cy + card_h - 0.050, msub, fontsize=9.4, fontstyle='italic', color='#475569', zorder=3)
        for p_idx, pt in enumerate(mpts):
            ax.text(0.04, cy + card_h - 0.078 - p_idx * 0.026, pt, fontsize=8.6, color='#1e293b', zorder=3)
        ax.text(0.96, cy + card_h - 0.028, mbtns, ha='right', fontsize=8.6, family='monospace', fontweight='bold', color=mcol, zorder=3)

    # Dedicated Workflow & Shortcuts Footer Card
    foot_box = FancyBboxPatch((0.02, 0.020), 0.96, 0.095, boxstyle="round,pad=0.002,rounding_size=0.008",
                              facecolor="#f8fafc", edgecolor="#94a3b8", lw=1.2, zorder=2)
    ax.add_patch(foot_box)
    ax.text(0.04, 0.088, "Recommended Exploratory Workflow & Keyboard Shortcuts:", fontsize=9.8, fontweight='bold', color='#0f172a', zorder=3)
    ax.text(0.04, 0.065, "• Step 1 (Tab 1): Inspect macroevolutionary distributions ➔ Step 2 (Tab 2): 1-on-1 syntenic contrasts ➔ Step 3 (Tab 3): Deep dive interval architecture", fontsize=8.4, color='#334155', zorder=3)
    ax.text(0.04, 0.042, "• Interactive Lightbox: Click any image to open fullscreen • Shortcuts: [Esc] Close • [0] Fit to screen • [+] / [-] Zoom • [1]–[6] Diagnostic panels", fontsize=8.4, color='#047857', fontweight='bold', zorder=3)


# =========================================================================
# GENERATORS FOR STANDALONE FIGURES
# =========================================================================
def generate_all_figures():
    # Figure 0: Dashboard Guide (First Figure)
    fig_0, ax_0 = plt.subplots(figsize=(18, 12), dpi=300)
    fig_0.patch.set_facecolor('#ffffff')
    draw_panel_d(ax_0, standalone=True)
    out_0_png = os.path.join(RESULTS_DIR, "summary_figure0_dashboard_guide.png")
    out_0_pdf = os.path.join(RESULTS_DIR, "summary_figure0_dashboard_guide.pdf")
    fig_0.savefig(out_0_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_0.savefig(out_0_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_0)
    print(f"✓ Figure 0 (Dashboard Guide) generated: {out_0_png}")

    # Also save as summary_figure0_overview_and_methodology and summary_figure0d for full compatibility
    shutil.copy2(out_0_png, os.path.join(RESULTS_DIR, "summary_figure0_overview_and_methodology.png"))
    shutil.copy2(out_0_pdf, os.path.join(RESULTS_DIR, "summary_figure0_overview_and_methodology.pdf"))
    shutil.copy2(out_0_png, os.path.join(RESULTS_DIR, "summary_figure0d_dashboard_guide.png"))
    shutil.copy2(out_0_pdf, os.path.join(RESULTS_DIR, "summary_figure0d_dashboard_guide.pdf"))

    # Figure 0A: Human Paradigm
    fig_a, ax_a = plt.subplots(figsize=(18, 12), dpi=300)
    fig_a.patch.set_facecolor('#ffffff')
    draw_panel_a(ax_a, standalone=True)
    out_a_png = os.path.join(RESULTS_DIR, "summary_figure0a_human_paradigm.png")
    out_a_pdf = os.path.join(RESULTS_DIR, "summary_figure0a_human_paradigm.pdf")
    fig_a.savefig(out_a_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_a.savefig(out_a_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_a)
    print(f"✓ Standalone Figure 0A generated: {out_a_png}")

    # Figure 0B: Synteny Pipeline
    fig_b, ax_b = plt.subplots(figsize=(18, 12), dpi=300)
    fig_b.patch.set_facecolor('#ffffff')
    draw_panel_b(ax_b, standalone=True)
    out_b_png = os.path.join(RESULTS_DIR, "summary_figure0b_synteny_pipeline.png")
    out_b_pdf = os.path.join(RESULTS_DIR, "summary_figure0b_synteny_pipeline.pdf")
    fig_b.savefig(out_b_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_b.savefig(out_b_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_b)
    print(f"✓ Standalone Figure 0B generated: {out_b_png}")

    # Figure 0C: Metrics & Codon Dynamics
    fig_c, ax_c = plt.subplots(figsize=(18, 13), dpi=300)
    fig_c.patch.set_facecolor('#ffffff')
    draw_panel_c(ax_c, standalone=True)
    out_c_png = os.path.join(RESULTS_DIR, "summary_figure0c_metric_anatomy_and_gc3.png")
    out_c_pdf = os.path.join(RESULTS_DIR, "summary_figure0c_metric_anatomy_and_gc3.pdf")
    fig_c.savefig(out_c_png, dpi=300, facecolor='#ffffff', bbox_inches='tight')
    fig_c.savefig(out_c_pdf, facecolor='#ffffff', bbox_inches='tight')
    plt.close(fig_c)
    print(f"✓ Standalone Figure 0C generated: {out_c_png}")


def sync_files():
    files_to_sync = [
        "summary_figure0_dashboard_guide.png",
        "summary_figure0_dashboard_guide.pdf",
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
    print("=== Generating Polished Overview, Methodology & Guide Figures ===")
    generate_all_figures()
    
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
