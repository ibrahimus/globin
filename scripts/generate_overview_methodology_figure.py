#!/usr/bin/env python3
"""
Generate Polished Overview, Methodology & Dashboard User Guide Figure:
"Figure 0: Conceptual Framework, Syntenic Retrieval Pipeline, Metric Anatomy & Interactive Dashboard Guide"

Sub-panels:
A. The Human Globin Paradigm & The Macroevolutionary Question
B. Cross-Vertebrate Synteny Retrieval & Karyotypic Tracking Pipeline (68 Genomes)
C. Anatomy of Genomic Metrics: Gene Architecture, Codon Triplet Physics & Compaction
D. The Analytical Framework & Interactive Dashboard Navigation Guide
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch, Polygon
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DATA_DIR = os.path.join(BASE_DIR, "data")
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
                ha='center', va='center', fontsize=10.5, fontweight='bold', color='#0f172a', zorder=zorder + 2)


def draw_gene_arrow(ax, x, y, w, h, name, fc="#0284c7", ec="#0369a1", text_col="white", zorder=4):
    """Draw a directed gene arrow box."""
    # Arrow head occupies 20% of width
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
    ax.text(x + body_w/2.0, y, name, ha='center', va='center', fontsize=7.8, fontweight='bold', color=text_col, zorder=zorder+1)


def draw_panel_a(ax):
    """PANEL A: The Human Globin Paradigm & The Macroevolutionary Question."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    # Panel Title
    ax.text(0.02, 0.965, "A. The Human Globin Paradigm & The Core Evolutionary Inquiry",
            fontsize=13.5, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Canonical dual-cluster arrangement on chromosomes 16 and 11, and the hypothesis of universal isochore polarity",
            fontsize=9.0, fontstyle='italic', color='#475569', va='top')

    # Card 1: Chr 16 - Alpha Globin Cluster (Subtelomeric Heavy Isochore)
    draw_card(ax, 0.02, 0.525, 0.96, 0.39, title="Human Alpha-Globin Cluster — Subtelomeric Heavy Isochore (16p13.3)",
              fc="#fdfefe", ec="#93c5fd", lw=1.4)
    
    # Isochore background tint
    iso_a = FancyBboxPatch((0.035, 0.535), 0.93, 0.33, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f0f9ff", edgecolor="none", zorder=1)
    ax.add_patch(iso_a)
    ax.text(0.95, 0.84, "Regional Heavy Isochore H3 (GC > 54%)", ha='right', fontsize=7.8, fontweight='bold', color='#0284c7', zorder=2)

    # Chr 16 Ideogram
    ax.plot([0.06, 0.92], [0.805, 0.805], color="#cbd5e1", lw=11, solid_capstyle='round', zorder=2)
    # Centromere
    ax.plot([0.42, 0.42], [0.785, 0.825], color="#475569", lw=3.0, zorder=3)
    ax.text(0.42, 0.832, "Cen (35.3 Mb)", fontsize=7.5, ha='center', color='#64748b', fontweight='bold')
    # Telomere cap (16p)
    ax.plot([0.06, 0.06], [0.785, 0.825], color="#dc2626", lw=3.5, zorder=4)
    ax.text(0.06, 0.832, "16p Telomere", fontsize=8.0, ha='left', color='#dc2626', fontweight='bold')
    
    # Locus Pin on Chr 16
    a_x = 0.075
    ax.plot([a_x, a_x], [0.805, 0.745], color="#0284c7", lw=1.8, zorder=5)
    ax.scatter([a_x], [0.805], color="#0284c7", s=36, zorder=6)
    ax.text(a_x, 0.732, "Alpha Locus: 0.15–0.22 Mb (0.16% from telomere — Subtelomeric)", fontsize=8.0, ha='left', va='top',
            color='#0369a1', fontweight='bold')

    # Detailed Alpha Gene Cluster Track
    ax.plot([0.06, 0.92], [0.635, 0.635], color="#cbd5e1", lw=1.5, zorder=2)
    
    # Syntenic Anchors (Alpha)
    anchors_alpha = [
        ("NPRL3", 0.07, 0.075, "#7c3aed"),
        ("MPG", 0.16, 0.060, "#7c3aed"),
        ("RHBDF1", 0.235, 0.075, "#7c3aed"),
        ("LUC7L", 0.81, 0.075, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_alpha:
        draw_gene_arrow(ax, gx, 0.635, gw, 0.048, gname, fc=gcol, ec="#5b21b6", text_col="white")
        ax.text(gx + gw/2.0, 0.595, "Anchor", ha='center', va='top', fontsize=7.0, color='#6b21a8')

    # Alpha Globin Genes
    globins_alpha = [
        ("HBZ (ζ)", 0.33, 0.075, "#0284c7"),
        ("HBM (μ)", 0.42, 0.070, "#0284c7"),
        ("HBA2 (α2)", 0.505, 0.075, "#1d4ed8"),
        ("HBA1 (α1)", 0.595, 0.075, "#1d4ed8"),
        ("HBQ1 (θ)", 0.685, 0.075, "#0284c7"),
    ]
    for gname, gx, gw, gcol in globins_alpha:
        draw_gene_arrow(ax, gx, 0.635, gw, 0.054, gname, fc=gcol, ec="#1e3a8a", text_col="white")

    # Metrics Badge for Alpha
    m_box_a = FancyBboxPatch((0.06, 0.540), 0.86, 0.042, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#eff6ff", edgecolor="#bfdbfe", lw=1, zorder=3)
    ax.add_patch(m_box_a)
    ax.text(0.49, 0.561,
            "α-Metrics:  CDS GC = 66.2%  |  GC3 = 91.1%  |  Gene GC = 68.2%  |  Mean Gene Length = 810 bp  |  Isochore: Heavy H3",
            ha='center', va='center', fontsize=8.2, fontweight='bold', color='#1e40af')

    # Card 2: Chr 11 - Beta Globin Cluster (Interstitial Light Isochore)
    draw_card(ax, 0.02, 0.115, 0.96, 0.39, title="Human Beta-Globin Cluster — Interstitial Light Isochore (11p15.4)",
              fc="#fdfefe", ec="#86efac", lw=1.4)
    
    # Isochore background tint
    iso_b = FancyBboxPatch((0.035, 0.125), 0.93, 0.33, boxstyle="round,pad=0.002,rounding_size=0.008",
                           facecolor="#f0fdf4", edgecolor="none", zorder=1)
    ax.add_patch(iso_b)
    ax.text(0.95, 0.43, "Regional Light Isochore L1/L2 (GC ~40%)", ha='right', fontsize=7.8, fontweight='bold', color='#16a34a', zorder=2)

    # Chr 11 Ideogram
    ax.plot([0.06, 0.92], [0.395, 0.395], color="#cbd5e1", lw=11, solid_capstyle='round', zorder=2)
    # Centromere
    ax.plot([0.45, 0.45], [0.375, 0.415], color="#475569", lw=3.0, zorder=3)
    ax.text(0.45, 0.422, "Cen (51.6 Mb)", fontsize=7.5, ha='center', color='#64748b', fontweight='bold')
    # Telomere cap (11p)
    ax.plot([0.06, 0.06], [0.375, 0.415], color="#dc2626", lw=3.5, zorder=4)
    ax.text(0.06, 0.422, "11p Telomere", fontsize=8.0, ha='left', color='#dc2626', fontweight='bold')
    
    # Locus Pin on Chr 11 (Interstitial, 5.2 Mb internal)
    b_x = 0.20
    ax.plot([b_x, b_x], [0.395, 0.335], color="#16a34a", lw=1.8, zorder=5)
    ax.scatter([b_x], [0.395], color="#16a34a", s=36, zorder=6)
    ax.text(b_x, 0.322, "Beta Locus: 5.22–5.27 Mb (3.85% from telomere — Interstitial Desert)", fontsize=8.0, ha='left', va='top',
            color='#15803d', fontweight='bold')

    # Detailed Beta Gene Cluster Track
    ax.plot([0.06, 0.92], [0.225, 0.225], color="#cbd5e1", lw=1.5, zorder=2)
    
    # Syntenic Anchors (Beta)
    anchors_beta = [
        ("FOLR1", 0.07, 0.075, "#7c3aed"),
        ("TRIM21", 0.16, 0.075, "#7c3aed"),
        ("OR51", 0.25, 0.060, "#7c3aed"),
        ("CDCA8", 0.81, 0.075, "#7c3aed"),
    ]
    for gname, gx, gw, gcol in anchors_beta:
        draw_gene_arrow(ax, gx, 0.225, gw, 0.048, gname, fc=gcol, ec="#5b21b6", text_col="white")
        ax.text(gx + gw/2.0, 0.185, "Anchor", ha='center', va='top', fontsize=7.0, color='#6b21a8')

    # Beta Globin Genes
    globins_beta = [
        ("HBE1 (ε)", 0.33, 0.075, "#16a34a"),
        ("HBG2 (Gγ)", 0.42, 0.075, "#16a34a"),
        ("HBG1 (Aγ)", 0.51, 0.075, "#16a34a"),
        ("HBD (δ)", 0.60, 0.070, "#15803d"),
        ("HBB (β)", 0.685, 0.070, "#15803d"),
    ]
    for gname, gx, gw, gcol in globins_beta:
        draw_gene_arrow(ax, gx, 0.225, gw, 0.054, gname, fc=gcol, ec="#14532d", text_col="white")

    # Metrics Badge for Beta
    m_box_b = FancyBboxPatch((0.06, 0.130), 0.86, 0.042, boxstyle="round,pad=0.002,rounding_size=0.006",
                             facecolor="#f0fdf4", edgecolor="#bbf7d0", lw=1, zorder=3)
    ax.add_patch(m_box_b)
    ax.text(0.49, 0.151,
            "β-Metrics:  CDS GC = 54.4%  |  GC3 = 65.8%  |  Gene GC = 43.4%  |  Mean Gene Length = 1,608 bp  |  Isochore: Light L1",
            ha='center', va='center', fontsize=8.2, fontweight='bold', color='#166534')

    # Core Evolutionary Question Banner
    q_box = FancyBboxPatch((0.02, 0.015), 0.96, 0.085, boxstyle="round,pad=0.002,rounding_size=0.010",
                           facecolor="#fffbeb", edgecolor="#f59e0b", lw=1.3, zorder=3)
    ax.add_patch(q_box)
    ax.text(0.04, 0.078, "The Core Evolutionary Inquiry:", fontsize=9.2, fontweight='bold', color='#b45309')
    ax.text(0.04, 0.048, "Is this dual-cluster, isochore-polarized architecture an ancient vertebrate feature, or a derived mammalian specialization?",
            fontsize=8.5, fontweight='bold', color='#1e293b')
    ax.text(0.04, 0.026, "How did 500 million years of chromosomal rearrangements, whole-genome duplications (2R/3R), and meiotic recombination shape globin synteny?",
            fontsize=8.0, color='#475569')


def draw_panel_b(ax):
    """PANEL B: Cross-Vertebrate Synteny Retrieval & Karyotypic Tracking Pipeline."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.02, 0.965, "B. Synteny-Anchored Discovery & Karyotypic Tracking (68 Vertebrates)",
            fontsize=13.5, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Systematic discovery pipeline across 10 vertebrate clades spanning 500 million years of divergence",
            fontsize=9.0, fontstyle='italic', color='#475569', va='top')

    # Step 1: Synteny Screening Funnel
    draw_card(ax, 0.02, 0.775, 0.96, 0.145, title="Step 1: Automated Syntenic Anchor Screening & Orthology Mapping",
              fc="#ffffff", ec="#e2e8f0", lw=1.2)
    
    ax.text(0.04, 0.852, "68 Genomes Screened across Tree of Life:", fontsize=8.8, fontweight='bold', color='#1e293b')
    ax.text(0.04, 0.828, "Agnatha • Chondrichthyes • Actinopterygii • Sarcopterygii • Amphibia • Testudines • Crocodylia • Aves • Squamata • Mammalia",
            fontsize=7.8, color='#0284c7')
    ax.text(0.04, 0.796, "Target Syntenic Markers:", fontsize=8.5, fontweight='bold', color='#334155')
    ax.text(0.23, 0.796, "α: NPRL3, MPG, RHBDF1, LUC7L, POLR3K  |  β: OR51, OR52, CDCA8, TRIM21, FOLR1, RRM1, STIM1",
            fontsize=7.8, color='#475569', family='monospace')

    # Step 2: Four Evolutionary Structural Archetypes
    draw_card(ax, 0.02, 0.015, 0.96, 0.745, title="Step 2: Four Core Macroevolutionary Chromosomal Archetypes",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    archetypes = [
        ("1. Primordial Linked Single Locus (Gar, Teleosts, Amphibians)",
         "α and β globin genes physically paired in head-to-head tandem arrays on a single chromosome arm.",
         "Identical chromatin context ensures parallel compositional evolution: ΔGC3 ≈ 0%  (α = β).",
         "#1565c0", "#eff6ff", 0.585, "Single Chromosome (e.g. Gar LG13 / Frog Chr 3)", "LINKED"),

        ("2. Canonical Decoupled Subtelomeric-α (Eutherian & Metatherian Mammals)",
         "Amniote translocation decouples α and β. α captured by subtelomere (<1 Mb to telomere); β remains interstitial.",
         "Runaway subtelomeric gBGC supercharges α-GC3 (>85–95%) while β-GC3 remains moderate (~60–65%). ΔGC3 >> 0%.",
         "#e65100", "#fff7ed", 0.400, "Two Chromosomes (α Subtelomeric / β Interstitial)", "MAMMALIAN"),

        ("3. Dual Microchromosomal Partitioning (Archosauria & Aves — Chicken, Eagle)",
         "Both α and β clusters translocate onto distinct, highly recombinogenic microchromosomes (Chr 14 and Chr 25).",
         "High crossover density drives synchronous hyper-elevation of both clusters: α-GC3 > 85%, β-GC3 > 85%. ΔGC3 ≈ 0%.",
         "#00796b", "#f0fdfa", 0.215, "Two Microchromosomes (Micro-14 and Micro-25)", "AVIAN"),

        ("4. The Great Inversion (Squamate Reptiles — Lizards & Snakes)",
         "Karyotypic inversion translocates β-locus to microchromosome/subtelomere, leaving α in an interstitial desert.",
         "Inverted gBGC polarity: β-GC3 skyrockets (>90–96%) while α-GC3 languishes (~48–55%). ΔGC3 << 0% (β >> α)!",
         "#d32f2f", "#fef2f2", 0.030, "Inverted Polarity (β Subtelomeric / α Interstitial)", "SQUAMATE"),
    ]

    for atitle, adesc1, adesc2, acol, afc, ay, achr, atype in archetypes:
        abox = FancyBboxPatch((0.04, ay), 0.92, 0.170, boxstyle="round,pad=0.002,rounding_size=0.008",
                              facecolor=afc, edgecolor=acol, lw=1.2, zorder=3)
        ax.add_patch(abox)
        ax.text(0.06, ay + 0.148, atitle, fontsize=8.8, fontweight='bold', color=acol, va='top', zorder=4)
        ax.text(0.06, ay + 0.124, adesc1, fontsize=7.8, color='#334155', va='top', zorder=4)
        ax.text(0.06, ay + 0.100, adesc2, fontsize=7.8, fontweight='bold', color='#1e293b', va='top', zorder=4)
        ax.text(0.06, ay + 0.076, f"Cytogenetic State: {achr}", fontsize=7.5, color='#64748b', va='top', zorder=4)

        # Draw Mini Chromosome Diagram for each archetype
        mini_y = ay + 0.035
        ax.plot([0.06, 0.88], [mini_y, mini_y], color="#cbd5e1", lw=8, solid_capstyle='round', zorder=4)
        
        if atype == "LINKED":
            # Single linked alpha-beta block on arm
            ax.plot([0.42, 0.42], [mini_y - 0.012, mini_y + 0.012], color="#475569", lw=2, zorder=5) # Cen
            ax.add_patch(FancyBboxPatch((0.70, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.add_patch(FancyBboxPatch((0.745, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.722, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.767, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.81, mini_y, "◀ Linked α-β (ΔGC3 ≈ 0)", fontsize=7.2, fontweight='bold', color='#1565c0', va='center', zorder=7)
            
        elif atype == "MAMMALIAN":
            # Top Chr: alpha subtelo; Bottom Chr: beta interstitial
            ax.plot([0.06, 0.06], [mini_y - 0.012, mini_y + 0.012], color="#dc2626", lw=2.5, zorder=5) # Telomere
            ax.add_patch(FancyBboxPatch((0.08, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.102, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.14, mini_y, "Subtelomeric α (GC3 > 90%)", fontsize=7.2, fontweight='bold', color='#0369a1', va='center', zorder=7)
            
            # Beta interstitial pin
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.622, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial β (GC3 ~ 65%)", fontsize=7.2, fontweight='bold', color='#15803d', va='center', zorder=7)

        elif atype == "AVIAN":
            # Two microchromosomes
            ax.add_patch(FancyBboxPatch((0.15, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.172, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.add_patch(FancyBboxPatch((0.55, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.572, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.21, mini_y, "Micro-14 (GC3 ~ 88%)", fontsize=7.2, fontweight='bold', color='#00796b', va='center', zorder=7)
            ax.text(0.61, mini_y, "Micro-25 (GC3 ~ 92%)", fontsize=7.2, fontweight='bold', color='#00796b', va='center', zorder=7)

        elif atype == "SQUAMATE":
            # Inverted polarity: Beta subtelomeric, Alpha interstitial
            ax.plot([0.06, 0.06], [mini_y - 0.012, mini_y + 0.012], color="#dc2626", lw=2.5, zorder=5) # Telomere
            ax.add_patch(FancyBboxPatch((0.08, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#43a047", ec="#1b5e20", lw=1, zorder=6))
            ax.text(0.102, mini_y, "β", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.14, mini_y, "Subtelomeric β (GC3 > 94%)", fontsize=7.2, fontweight='bold', color='#15803d', va='center', zorder=7)
            
            # Alpha interstitial pin
            ax.add_patch(FancyBboxPatch((0.60, mini_y - 0.018), 0.045, 0.036, boxstyle="square,pad=0", fc="#1e88e5", ec="#0d47a1", lw=1, zorder=6))
            ax.text(0.622, mini_y, "α", ha='center', va='center', fontsize=7, fontweight='bold', color='white', zorder=7)
            ax.text(0.66, mini_y, "Interstitial α (GC3 ~ 50%) — INVERTED POLARITY!", fontsize=7.2, fontweight='bold', color='#b91c1c', va='center', zorder=7)


def draw_panel_c(ax):
    """PANEL C: Anatomy of Genomic Metrics: Gene Architecture, Codon Triplet Physics & Compaction."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.02, 0.965, "C. Anatomy of Calculated Metrics: Gene Architecture & Codon Physics",
            fontsize=13.5, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "Biophysical basis of CDS GC%, GC3%, Whole-Gene GC%, Intron Expansion, and Telomeric Gradients",
            fontsize=9.0, fontstyle='italic', color='#475569', va='top')

    # Section 1: Gene Structure Decomposition
    draw_card(ax, 0.02, 0.525, 0.96, 0.395, title="1. Gene Architecture: Full Gene GC% vs. CDS GC% vs. Codon Positions",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    # Gene schematic diagram
    ax.plot([0.06, 0.92], [0.810, 0.810], color="#94a3b8", lw=2.5, zorder=2)
    # UTRs
    ax.add_patch(FancyBboxPatch((0.06, 0.790), 0.055, 0.040, boxstyle="square,pad=0", fc="#e2e8f0", ec="#64748b", lw=1, zorder=3))
    ax.text(0.087, 0.810, "5' UTR", ha='center', va='center', fontsize=7.0, color='#475569')
    ax.add_patch(FancyBboxPatch((0.880, 0.790), 0.055, 0.040, boxstyle="square,pad=0", fc="#e2e8f0", ec="#64748b", lw=1, zorder=3))
    ax.text(0.907, 0.810, "3' UTR", ha='center', va='center', fontsize=7.0, color='#475569')

    # Exons (CDS) - carefully proportioned with wide intron gaps
    exons = [
        ("Exon 1 (~95 bp)", 0.125, 0.130, "#1d4ed8"),
        ("Exon 2 (~220 bp)", 0.365, 0.150, "#1d4ed8"),
        ("Exon 3 (~130 bp)", 0.740, 0.130, "#1d4ed8"),
    ]
    for ename, ex, ew, ecol in exons:
        ax.add_patch(FancyBboxPatch((ex, 0.775), ew, 0.070, boxstyle="round,pad=0.002,rounding_size=0.006",
                                   fc=ecol, ec="#1e3a8a", lw=1.1, zorder=4))
        ax.text(ex + ew/2.0, 0.810, ename, ha='center', va='center', fontsize=7.8, fontweight='bold', color='white', zorder=5)

    # Introns - perfectly centered in wide gaps with zero overlap
    ax.text(0.310, 0.852, "Intron 1", ha='center', fontsize=7.2, color='#475569', fontweight='bold')
    ax.text(0.310, 0.835, "(~120 bp)", ha='center', fontsize=6.8, color='#64748b')

    ax.text(0.628, 0.852, "Intron 2 (Compacted in α vs. Expanded in β)", ha='center', fontsize=7.2, color='#b91c1c', fontweight='bold')
    ax.text(0.628, 0.835, "(140 bp in HBA1  vs.  850 bp in HBB)", ha='center', fontsize=6.8, color='#991b1b')

    # Metric Definitions
    ax.text(0.06, 0.735, "Whole-Gene GC%: Exons + Introns + UTRs (Reflects regional chromatin & isochore baseline GC)",
            fontsize=8.2, fontweight='bold', color='#1e293b')
    ax.text(0.06, 0.710, "CDS GC%: Coding sequence only (Concatenated Exons 1, 2, 3; reflects amino acid & mutational bias)",
            fontsize=8.2, fontweight='bold', color='#1d4ed8')

    # Codon Triplet Zoom In Card with Purifying Corridor Visual
    trip_box = FancyBboxPatch((0.05, 0.535), 0.90, 0.160, boxstyle="round,pad=0.002,rounding_size=0.006",
                              fc="#f8fafc", ec="#cbd5e1", lw=1, zorder=3)
    ax.add_patch(trip_box)
    
    # 3 Codon Positions
    ax.text(0.07, 0.672, "Codon Triplet Decomposition & Biophysical Roles:", fontsize=8.8, fontweight='bold', color='#0f172a')
    
    positions = [
        ("Position 1 (GC1)", "Moderate purifying constraint; conservative amino acid substitutions (~50–65% GC)", "#0284c7"),
        ("Position 2 (GC2)", "Extreme purifying constraint! Strictly locked in Purifying Selection Corridor (45–55%)", "#d97706"),
        ("Position 3 (GC3)", "70% synonymous/neutral! Pure biophysical barometer of recombination & gBGC (35–98% GC)", "#dc2626"),
        ("4-Fold Degenerate (GC4d)", "100% synonymous (Gly, Ala, Val, Pro, Thr, Arg, Ser, Leu). Pure mutational/gBGC ground truth", "#7c3aed"),
    ]
    for idx, (pname, pdesc, pcol) in enumerate(positions):
        py = 0.648 - idx * 0.024
        ax.scatter([0.075], [py + 0.003], color=pcol, s=26, zorder=4)
        ax.text(0.095, py, f"{pname}: {pdesc}", fontsize=7.5, color='#334155', va='center')

    # Mini Purifying Selection Corridor Axis Graphic
    corr_y = 0.548
    ax.plot([0.62, 0.92], [corr_y, corr_y], color="#94a3b8", lw=1.5, zorder=4)
    # Shaded Corridor
    ax.add_patch(patches.Rectangle((0.68, corr_y - 0.008), 0.08, 0.016, facecolor="#bbf7d0", edgecolor="#16a34a", lw=0.8, zorder=5))
    ax.text(0.72, corr_y + 0.012, "Purifying Corridor (45–55%)", ha='center', fontsize=6.8, color='#15803d', fontweight='bold', zorder=6)
    # Points
    ax.scatter([0.70], [corr_y], color="#d97706", s=22, zorder=7) # GC2
    ax.scatter([0.74], [corr_y], color="#0284c7", s=22, zorder=7) # GC1
    ax.scatter([0.90], [corr_y], color="#dc2626", s=28, zorder=7) # GC3
    ax.text(0.70, corr_y - 0.012, "GC2", ha='center', fontsize=6.5, color='#d97706', fontweight='bold')
    ax.text(0.74, corr_y - 0.012, "GC1", ha='center', fontsize=6.5, color='#0284c7', fontweight='bold')
    ax.text(0.90, corr_y - 0.012, "GC3 ➔", ha='center', fontsize=6.5, color='#dc2626', fontweight='bold')

    # Section 2: Bernardi's Compaction & Telomeric Gradient
    draw_card(ax, 0.02, 0.015, 0.96, 0.485, title="2. Isochore Compaction Theory & Distance-to-Telomere Drivers",
              fc="#ffffff", ec="#cbd5e1", lw=1.2)

    # Gene Compaction Graphic
    ax.text(0.05, 0.450, "Bernardi's Genome Compaction Law:", fontsize=9.0, fontweight='bold', color='#1e293b')
    ax.text(0.05, 0.425, "High-GC isochores suppress intron expansion; low-GC isochores accumulate massive introns.", fontsize=7.8, color='#475569')

    # Human HBA1 vs HBB length bars
    ax.text(0.05, 0.380, "Human HBA1 (α-globin):", fontsize=8.0, fontweight='bold', color='#1d4ed8')
    ax.add_patch(FancyBboxPatch((0.05, 0.345), 0.22, 0.022, boxstyle="round,pad=0.001,rounding_size=0.004", fc="#1d4ed8", ec="none"))
    ax.text(0.285, 0.356, "843 bp total (Intron 2 = 140 bp)  [COMPACT HIGH-GC]", fontsize=7.5, color='#1e40af', fontweight='bold', va='center')

    ax.text(0.05, 0.310, "Human HBB (β-globin):", fontsize=8.0, fontweight='bold', color='#15803d')
    ax.add_patch(FancyBboxPatch((0.05, 0.275), 0.42, 0.022, boxstyle="round,pad=0.001,rounding_size=0.004", fc="#15803d", ec="none"))
    ax.text(0.485, 0.286, "1,608 bp total (Intron 2 = 850 bp)  [EXPANDED LOW-GC]", fontsize=7.5, color='#166534', fontweight='bold', va='center')

    # Telomeric Distance Driver
    ax.text(0.05, 0.215, "Distance-to-Telomere (dtelo / Lchr) as the Primary Engine of GC Polarity:",
            fontsize=9.0, fontweight='bold', color='#0f172a')
    ax.text(0.05, 0.180, "• Meiotic crossovers are physically biased toward chromosome ends (especially in male meiosis).",
            fontsize=7.8, color='#334155')
    ax.text(0.05, 0.155, "• Recombination-associated double-strand break repair preferentially resolves with GC alleles (gBGC).",
            fontsize=7.8, color='#334155')
    ax.text(0.05, 0.130, "• Subtelomeric loci experience runaway GC3 elevation; interstitial loci remain buffered in AT-rich deserts.",
            fontsize=7.8, color='#b91c1c', fontweight='bold')
    
    # Distance Formula box
    f_box = FancyBboxPatch((0.05, 0.045), 0.90, 0.065, boxstyle="round,pad=0.002,rounding_size=0.006",
                           facecolor="#f1f5f9", edgecolor="#cbd5e1", lw=1)
    ax.add_patch(f_box)
    ax.text(0.50, 0.088, "Relative Telomeric Distance Formula:", ha='center', fontsize=8.0, fontweight='bold', color='#334155')
    ax.text(0.50, 0.062, "dtelo_rel = min(Locus_Pos, Chr_Len - Locus_Pos) / Chr_Len    (Range: 0.000 Subtelomeric  to  0.500 Metacentric/Internal)",
            ha='center', fontsize=7.5, family='monospace', color='#0f172a', fontweight='bold')


def draw_panel_d(ax):
    """PANEL D: The Analytical Framework & Interactive Dashboard Navigation Guide."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title - lowered to avoid overlapping card
    ax.text(0.02, 0.965, "D. Interactive Dashboard Architecture & Analytical User Guide",
            fontsize=13.5, fontweight='bold', color='#0f172a', va='top')
    ax.text(0.02, 0.935, "How to navigate and interpret the 5 interconnected exploratory modules in dashboard.html",
            fontsize=9.0, fontstyle='italic', color='#475569', va='top')

    # Dashboard Modules Overview
    modules = [
        ("Tab 1: Global Macroevolutionary Figures",
         "The Big-Picture Evolutionary Synthesis",
         ["• Fig 1A: Curated benchmark summary across 33 representative species.",
          "• Fig 1B: Exhaustive distribution across all 68 vertebrate genomes.",
          "• Fig 2: Mechanistic drivers (Bivariate GC regressions, gene length compaction, telomere distance).",
          "• Fig 3: Time-calibrated phylogenetic chronogram, chromosome synteny maps, and codon ΔGC3 divergence."],
         "#1e40af", "#eff6ff", 0.690, "[ Fig 1A ]  [ Fig 1B ]  [ Fig 2 ]  [ Fig 3 ]"),

        ("Tab 2: 1-on-1 Head-to-Head Comparison",
         "Pairwise Species Contrast & Polarity Scorecard",
         ["• Select Species A and Species B from full 68-genome searchable dropdowns.",
          "• Live Comparative Scorecard: Instant ΔCDS GC, ΔGC3, and Isochore Polarity classification.",
          "• Synchronized views: Linear Chromosome Tracks, Bivariate with Mini Inset, 6-Panel Diagnostic Distributions.",
          "• Click any image or press Space/Enter to launch full-screen interactive pan & zoom modal."],
         "#065f46", "#ecfdf5", 0.465, "[ Species A: Human ▾ ]  vs  [ Species B: Corn Snake ▾ ]   [ Inversion: ΔGC3 = -38% ]"),

        ("Tab 3: Single Species Deep Dive",
         "High-Resolution Multi-Track Locus Inspection",
         ["• Inspect individual chromosome coordinates, centromere locations, and telomeric distances.",
          "• View 6-panel diagnostic distributions contrasting globin loci against the whole-genome baseline.",
          "• Focus on individual subpanels (Panels A–F) using keyboard shortcuts 1–6 in lightbox mode."],
         "#854d0e", "#fefce8", 0.240, "[ 6-Panel Diagnostic Suite ]  [ Interactive Lightbox Pan & Zoom ]"),

        ("Tabs 4 & 5: Interval Database & Statistical Rigor",
         "Raw Gene Data & Hypothesis Testing",
         ["• Tab 4 (Interval Genes): Searchable table of all 1,180 globin and flanking syntenic genes with GC metrics.",
          "• Tab 5 (Statistical Contrasts): Non-parametric Wilcoxon tests, Z-score distributions, and empirical percentiles.",
          "• Downloadable CSV exports for independent bioinformatics replication."],
         "#581c87", "#faf5ff", 0.015, "[ 1,180 Genes DB ]  [ Non-Parametric Wilcoxon Tests ]  [ Z-Scores ]"),
    ]

    for mtitle, msub, mbullets, mcol, mfc, my, mui in modules:
        mbox = FancyBboxPatch((0.02, my), 0.96, 0.205, boxstyle="round,pad=0.002,rounding_size=0.010",
                             facecolor=mfc, edgecolor=mcol, lw=1.2, zorder=2)
        ax.add_patch(mbox)
        ax.text(0.04, my + 0.182, mtitle, fontsize=9.8, fontweight='bold', color=mcol, va='top', zorder=3)
        ax.text(0.04, my + 0.156, msub, fontsize=8.0, fontstyle='italic', color='#475569', va='top', zorder=3)
        for bidx, bullet in enumerate(mbullets):
            by = my + 0.128 - bidx * 0.023
            ax.text(0.04, by, bullet, fontsize=7.6, color='#1e293b', va='top', zorder=3)
            
        # UI mockup button bar at bottom of card
        ax.text(0.04, my + 0.022, f"UI Elements: {mui}", fontsize=7.2, family='monospace', fontweight='bold', color=mcol, va='bottom', zorder=3)


def plot_overview_publication(out_dir=RESULTS_DIR):
    """Generate high-DPI publication vector PDF and raster PNG for Figure 0."""
    png_path = os.path.join(out_dir, "summary_figure0_overview_and_methodology.png")
    pdf_path = os.path.join(out_dir, "summary_figure0_overview_and_methodology.pdf")
    
    png_alt = os.path.join(out_dir, "overview_methodology_figure.png")
    pdf_alt = os.path.join(out_dir, "overview_methodology_figure.pdf")

    fig = plt.figure(figsize=(28, 20), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.14, wspace=0.12,
                           left=0.03, right=0.97, top=0.95, bottom=0.03)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)
    draw_panel_d(ax_d)

    # Master Figure Supertitle
    fig.suptitle("Figure 0: Conceptual Paradigm, Syntenic Discovery Pipeline, Metric Anatomy & Interactive Dashboard User Guide",
                 fontsize=17, fontweight='bold', y=0.985, color='#0f172a')

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.savefig(png_alt, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_alt, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Generated Overview Figure (PNG/PDF): {png_path} and {pdf_path}")
    return png_path, pdf_path


def plot_overview_interactive(out_dir=RESULTS_DIR):
    """Generate interactive HTML preview for Figure 0."""
    html_path = os.path.join(out_dir, "summary_figure0_overview_and_methodology.html")
    html_alt = os.path.join(out_dir, "overview_methodology_figure.html")

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Overview &amp; Methodology Guide — Vertebrate Globin Genomics</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 20px; line-height: 1.5; }
        .container { max-width: 1400px; margin: 0 auto; background: #1e293b; border-radius: 12px; padding: 25px; border: 1px solid #334155; box-shadow: 0 10px 25px rgba(0,0,0,0.5); }
        h1 { margin-top: 0; color: #38bdf8; font-size: 1.6rem; border-bottom: 1px solid #334155; padding-bottom: 12px; }
        .subtitle { color: #94a3b8; font-size: 0.95rem; margin-top: -6px; margin-bottom: 20px; }
        .img-box { background: white; border-radius: 8px; overflow: hidden; border: 1px solid #475569; text-align: center; }
        .img-box img { width: 100%; height: auto; display: block; }
        .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 25px; }
        .card { background: #0f172a; border-radius: 8px; padding: 18px; border: 1px solid #334155; }
        .card h3 { margin-top: 0; font-size: 1.1rem; color: #38bdf8; }
        ul { padding-left: 20px; color: #cbd5e1; font-size: 0.9rem; }
        li { margin-bottom: 8px; }
        code { background: #334155; padding: 2px 6px; border-radius: 4px; font-size: 0.85rem; color: #f1f5f9; }
        .back-bar { margin-top: 25px; padding-top: 15px; border-top: 1px solid #334155; display: flex; justify-content: space-between; align-items: center; }
        .btn { background: #0284c7; color: white; border: none; padding: 8px 16px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 0.9rem; cursor: pointer; }
        .btn:hover { background: #0369a1; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Overview &amp; Methodology Guide</h1>
        <div class="subtitle">Conceptual Framework, Syntenic Retrieval Pipeline, Metric Anatomy &amp; Interactive Dashboard User Guide</div>
        
        <div class="img-box">
            <img src="summary_figure0_overview_and_methodology.png" alt="Overview and Methodology Figure">
        </div>

        <div class="grid-2">
            <div class="card">
                <h3>Key Scientific Insights</h3>
                <ul>
                    <li><b>The Human Paradigm:</b> Human α (16p) is subtelomeric in an extreme heavy isochore (GC3 = 91%), while β (11p) is interstitial in a light isochore (GC3 = 66%).</li>
                    <li><b>Not Universal:</b> Teleosts, Gar, and Amphibians maintain a single linked locus with identical GC context (ΔGC3 ≈ 0%).</li>
                    <li><b>The Great Inversion:</b> Squamate reptiles invert this polarity completely: β is subtelomeric/microchromosomal (GC3 &gt; 90%), while α is interstitial!</li>
                    <li><b>Purifying Corridor:</b> GC1 and GC2 are constrained to 45–55% across all 68 species, while GC3 ranges unconstrained from 35% to 98%.</li>
                </ul>
            </div>
            <div class="card">
                <h3>How to Navigate the Dashboard</h3>
                <ul>
                    <li><b>Tab 1 (Macro Figures):</b> Explore clade-level divergence across 500 Ma (Fig 1A, 1B, 2, 3).</li>
                    <li><b>Tab 2 (1-on-1 Comparison):</b> Select any two species to generate side-by-side chromosome tracks, bivariate insets, and diagnostic distributions.</li>
                    <li><b>Tab 3 (Deep Dive):</b> Inspect individual locus coordinates, centromeres, and telomeric distances.</li>
                    <li><b>Tabs 4 &amp; 5 (Database &amp; Stats):</b> Filter 1,180 genes and review non-parametric statistical tests.</li>
                </ul>
            </div>
        </div>

        <div class="back-bar">
            <span style="color: #64748b; font-size: 0.85rem;">Vertebrate Globin Genomics Explorer • High-Resolution Diagnostic Suite</span>
            <a href="dashboard.html" class="btn">Open Main Dashboard →</a>
        </div>
    </div>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    with open(html_alt, "w", encoding="utf-8") as f:
        f.write(html_content)
        
    print(f"Generated Overview Interactive HTML: {html_path}")
    return html_path


def main():
    print("=== Generating Polished Figure 0: Overview, Methodology & Dashboard Guide ===")
    png_path, pdf_path = plot_overview_publication(out_dir=RESULTS_DIR)
    html_path = plot_overview_interactive(out_dir=RESULTS_DIR)

    # Synchronize to 1_data_from_mike/results/
    for fname in ["summary_figure0_overview_and_methodology.png",
                  "summary_figure0_overview_and_methodology.pdf",
                  "summary_figure0_overview_and_methodology.html",
                  "overview_methodology_figure.png",
                  "overview_methodology_figure.pdf",
                  "overview_methodology_figure.html"]:
        src = os.path.join(RESULTS_DIR, fname)
        dst = os.path.join(SYNC_DIR, fname)
        if os.path.exists(src):
            with open(src, "rb") as f_in, open(dst, "wb") as f_out:
                f_out.write(f_in.read())
            print(f"  Synchronized: {fname}")

    print("\nFigure 0 generation and synchronization complete!")


if __name__ == "__main__":
    main()
