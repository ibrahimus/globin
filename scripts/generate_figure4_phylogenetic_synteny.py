"""Generate Figure 4: Macroevolutionary Synteny and Codon Divergence (GC1/GC2/GC3) Map
across the Vertebrate Tree (Polished Publication Version).

Combines:
1. Time-calibrated phylogenetic tree (520 Ma -> 0 Ma) with orthogonal branch layout.
2. Syntenic chromosome schematics:
   - Bichromosomal species display TWO distinct chromosome ribbons (one for α-chromosome, one for β-chromosome)
     with realistic chromosomal coordinates and physical megabase spans.
   - Non-amniote tetrapods and fish display ONE single chromosome ribbon representing the ancestral linked α-β locus.
3. Paired codon position decomposition (GC1, GC2, GC3):
   - Purifying selection corridor (GC1 and GC2: 45–55%).
   - Subunit Divergence Bar connecting α-GC3 and β-GC3, colored by a continuous diverging colormap (coolwarm)
     depicting ΔGC3 (α vs β polarity).
   - Dedicated legend positioned cleanly below data points without any overlap.
4. Landmark macroevolutionary transition callouts (N1 through N8).

Outputs:
- results/summary_figure4_phylogenetic_synteny.png (300 dpi)
- results/summary_figure4_phylogenetic_synteny.pdf (vector)
- results/summary_figure4_phylogenetic_synteny.html (interactive Plotly)
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# -------------------------------------------------------------
# 22 Representative Taxa with Realistic Chromosomal Architecture
# -------------------------------------------------------------
REPRESENTATIVE_TAXA = [
    # 0. Agnatha
    {
        'species': 'petromyzon_marinus',
        'display': 'Sea Lamprey (P. marinus)',
        'clade': 'Agnatha',
        'divergence_mya': 500,
        'state': 'CYCLOSTOME_PROTO',
        'color': '#7b1fa2',
        'is_bichromosomal': False,
        'chr_label': 'GL478636',
        'pos': 0.03,
        'mb_range': '[0.04–0.09 Mb]',
        'anchors': 'Ancestral Proto-Globin',
        'notes': 'Single ancestral proto-globin locus; monomeric/homomeric',
    },
    # 1. Chondrichthyes - Grey Reef Shark
    {
        'species': 'carcharhinus_amblyrhynchos',
        'display': 'Grey Reef Shark (C. amblyrhynchos)',
        'clade': 'Chondrichthyes',
        'divergence_mya': 435,
        'state': 'CHONDRICHTHYAN_GRADIENT',
        'color': '#00838f',
        'is_bichromosomal': False,
        'chr_label': 'Chr 18',
        'pos': 0.04,
        'mb_range': '[3.8–4.0 Mb]',
        'anchors': '◀ NPRL3-MPG (linked)',
        'notes': 'Single linked locus on Chr 18; subtelomeric distance governs GC',
    },
    # 2. Chondrichthyes - Catshark
    {
        'species': 'scyliorhinus_canicula',
        'display': 'Catshark (S. canicula)',
        'clade': 'Chondrichthyes',
        'divergence_mya': 430,
        'state': 'CHONDRICHTHYAN_GRADIENT',
        'color': '#00838f',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 15',
        'alpha_pos': 0.02,
        'alpha_mb': '[3.8 Mb]',
        'alpha_anchors': 'NPRL3-RHBDF1',
        'beta_chr': 'Chr 21',
        'beta_pos': 0.45,
        'beta_mb': '[73.6 Mb]',
        'beta_anchors': 'Deep Interstitial',
        'notes': 'Bichromosomal; α subtelomeric, β deep interstitial',
    },
    # 3. Basal Actinopterygii - Spotted Gar
    {
        'species': 'lepisosteus_oculatus',
        'display': 'Spotted Gar (L. oculatus)',
        'clade': 'Holostei',
        'divergence_mya': 350,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1565c0',
        'is_bichromosomal': False,
        'chr_label': 'LG13',
        'pos': 0.43,
        'mb_range': '[14.7–14.9 Mb]',
        'anchors': '◀ NPRL3-MPG-RHBDF1 (linked)',
        'notes': 'Pre-teleost 3R state; single linked α-β cluster on LG13',
    },
    # 4. Teleostei - Zebrafish
    {
        'species': 'danio_rerio',
        'display': 'Zebrafish (D. rerio)',
        'clade': 'Teleostei',
        'divergence_mya': 310,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1976d2',
        'is_bichromosomal': False,
        'chr_label': 'Chr 3 & 12',
        'pos': 0.03,
        'mb_range': '[20.2–20.5 Mb]',
        'anchors': '◀ NPRL3-MPG (linked)',
        'notes': '3R duplicate loci (LA Chr 3, LB Chr 12); tandem linked α-β pairs',
    },
    # 5. Teleostei - Fugu
    {
        'species': 'takifugu_rubripes',
        'display': 'Fugu (T. rubripes)',
        'clade': 'Teleostei',
        'divergence_mya': 310,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1976d2',
        'is_bichromosomal': False,
        'chr_label': 'Chr 5 & 1',
        'pos': 0.43,
        'mb_range': '[10.3–10.5 Mb]',
        'anchors': '◀ NPRL3-MPG (linked)',
        'notes': 'Ultra-compact genome; linked α-β pairs retain parallel GC',
    },
    # 6. Sarcopterygii - Coelacanth
    {
        'species': 'latimeria_chalumnae',
        'display': 'Coelacanth (L. chalumnae)',
        'clade': 'Sarcopterygii',
        'divergence_mya': 415,
        'state': 'LINKED_ANCESTRAL',
        'color': '#2e7d32',
        'is_bichromosomal': False,
        'chr_label': 'Chr 5',
        'pos': 0.03,
        'mb_range': '[166.2–166.7 Mb]',
        'anchors': '◀ NPRL3-LUC7L (linked)',
        'notes': 'Basal sarcopterygian; retains ancestral linked α-β locus on Chr 5',
    },
    # 7. Amphibia - Clawed Frog
    {
        'species': 'xenopus_tropicalis',
        'display': 'Clawed Frog (X. tropicalis)',
        'clade': 'Amphibia',
        'divergence_mya': 270,
        'state': 'LINKED_ANCESTRAL',
        'color': '#388e3c',
        'is_bichromosomal': False,
        'chr_label': 'Chr 9',
        'pos': 0.35,
        'mb_range': '[40.3–41.3 Mb]',
        'anchors': '◀ NPRL3-RHBDF1-LUC7L (linked)',
        'notes': 'Non-amniote tetrapod; retains single linked α-β locus on Chr 9 (hba + hbg/hbe)',
    },
    # 8. Amphibia - Axolotl
    {
        'species': 'ambystoma_mexicanum',
        'display': 'Axolotl (A. mexicanum)',
        'clade': 'Amphibia',
        'divergence_mya': 270,
        'state': 'LINKED_ANCESTRAL',
        'color': '#388e3c',
        'is_bichromosomal': False,
        'chr_label': 'Chr 2',
        'pos': 0.45,
        'mb_range': '[681–689 Mb]',
        'anchors': '◀ NPRL3-RHBDF1 (linked)',
        'notes': 'Giant salamander genome (32 Gb); retains ancestral linked α-β locus on Chr 2',
    },
    # 9. Testudines - Green Sea Turtle
    {
        'species': 'chelonia_mydas',
        'display': 'Green Sea Turtle (C. mydas)',
        'clade': 'Testudines',
        'divergence_mya': 245,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr NC_051250',
        'alpha_pos': 0.37,
        'alpha_mb': '[31.2 Mb]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'Chr NC_057849',
        'beta_pos': 0.01,
        'beta_mb': '[3.5 Mb subtelo]',
        'beta_anchors': 'OR51-RRM1',
        'notes': 'Basal archelosaur; microchromosomal partitioning of α and β',
    },
    # 10. Crocodylia - American Alligator
    {
        'species': 'alligator_mississippiensis',
        'display': 'American Alligator (A. mississippiensis)',
        'clade': 'Crocodylia',
        'divergence_mya': 240,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr NC_081836',
        'alpha_pos': 0.37,
        'alpha_mb': '[21.2 Mb]',
        'alpha_anchors': 'NPRL3-LUC7L',
        'beta_chr': 'Chr NC_081824',
        'beta_pos': 0.01,
        'beta_mb': '[4.4 Mb subtelo]',
        'beta_anchors': 'RRM1-STIM1',
        'notes': 'Archosaur karyotype; slow morphological evolution, stable GC',
    },
    # 11. Aves - Chicken
    {
        'species': 'Gallus_gallus',
        'display': 'Chicken (G. gallus)',
        'clade': 'Aves',
        'divergence_mya': 100,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'is_bichromosomal': True,
        'alpha_chr': 'micro-Chr 14',
        'alpha_pos': 0.22,
        'alpha_mb': '[12.0 Mb]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'macro-Chr 1',
        'beta_pos': 0.003,
        'beta_mb': '[196.0 Mb subtelo]',
        'beta_anchors': 'FOLR1-OR51',
        'notes': 'α on micro-Chr 14; β on macro-Chr 1 subtelomere',
    },
    # 12. Aves - Hummingbird
    {
        'species': 'calypte_anna',
        'display': "Anna's Hummingbird (C. anna)",
        'clade': 'Aves',
        'divergence_mya': 100,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'is_bichromosomal': True,
        'alpha_chr': 'micro-Chr 14',
        'alpha_pos': 0.07,
        'alpha_mb': '[14.2 Mb]',
        'alpha_anchors': 'NPRL3-LUC7L',
        'beta_chr': 'micro-Chr 25',
        'beta_pos': 0.05,
        'beta_mb': '[Subtelomeric]',
        'beta_anchors': 'OR51-RRM1',
        'notes': 'Ultra-high metabolic rate; both loci on microchromosomes',
    },
    # 13. Squamata - Wall Lizard
    {
        'species': 'podarcis_muralis',
        'display': 'Wall Lizard (P. muralis)',
        'clade': 'Squamata',
        'divergence_mya': 165,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#d32f2f',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr NC_135668',
        'alpha_pos': 0.02,
        'alpha_mb': '[1.2 Mb subtelo]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'Chr NC_135658',
        'beta_pos': 0.01,
        'beta_mb': '[1.6 Mb subtelo]',
        'beta_anchors': 'CDCA8 (subtelo)',
        'notes': 'Both α and β in terminal subtelomeres; parallel hyper-GC surges',
    },
    # 14. Squamata - Corn Snake
    {
        'species': 'pantherophis_guttatus',
        'display': 'Corn Snake (P. guttatus)',
        'clade': 'Squamata',
        'divergence_mya': 110,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#c2185b',
        'is_bichromosomal': True,
        'alpha_chr': 'macro-Chr 4',
        'alpha_pos': 0.50,
        'alpha_mb': '[22.9 Mb internal]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'micro-Chr',
        'beta_pos': 0.006,
        'beta_mb': '[0.7 Mb subtelo]',
        'beta_anchors': '◀ CDCA8-OR51-RRM1',
        'notes': 'The Squamate Inversion: β translocated to microchromosome subtelomere',
    },
    # 15. Squamata - Common Adder
    {
        'species': 'vipera_berus',
        'display': 'Common Adder (V. berus)',
        'clade': 'Squamata',
        'divergence_mya': 110,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#c2185b',
        'is_bichromosomal': True,
        'alpha_chr': 'macro-Chr 4',
        'alpha_pos': 0.50,
        'alpha_mb': '[143 kb internal]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'micro-Chr',
        'beta_pos': 0.005,
        'beta_mb': '[563 kb subtelo]',
        'beta_anchors': '◀ RRM1-STIM1',
        'notes': 'Extreme squamate polarity: β-GC3 reaches 93.3% while α remains low',
    },
    # 16. Mammalia - Platypus
    {
        'species': 'ornithorhynchus_anatinus',
        'display': 'Platypus (O. anatinus)',
        'clade': 'Monotremata',
        'divergence_mya': 165,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#f57c00',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 21',
        'alpha_pos': 0.36,
        'alpha_mb': '[15.7 Mb]',
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_chr': 'Chr 2',
        'beta_pos': 0.19,
        'beta_mb': '[139.4 Mb internal]',
        'beta_anchors': 'OR51',
        'notes': 'Basal mammal; α and β decoupled, transitional isochore structure',
    },
    # 17. Mammalia - Gray Opossum
    {
        'species': 'Monodelphis_domestica',
        'display': 'Gray Opossum (M. domestica)',
        'clade': 'Marsupialia',
        'divergence_mya': 140,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#f57c00',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 6',
        'alpha_pos': 0.49,
        'alpha_mb': '[149.1 Mb]',
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_chr': 'Chr 4',
        'beta_pos': 0.19,
        'beta_mb': '[352.2 Mb internal]',
        'beta_anchors': 'FOLR1-OR51',
        'notes': 'Metatherian model; intermediate GC levels, stable autosomes',
    },
    # 18. Mammalia - Human
    {
        'species': 'Homo_sapiens',
        'display': 'Human (H. sapiens)',
        'clade': 'Primates',
        'divergence_mya': 85,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 16',
        'alpha_pos': 0.002,
        'alpha_mb': '[0.15 Mb subtelo]',
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_chr': 'Chr 11',
        'beta_pos': 0.039,
        'beta_mb': '[5.2 Mb interstitial]',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Canonical mammalian paradigm: α at extreme 16p telomere (GC3=91.1%), β at 11p (65.8%)',
    },
    # 19. Mammalia - Mouse
    {
        'species': 'Mus_musculus',
        'display': 'House Mouse (M. musculus)',
        'clade': 'Rodentia',
        'divergence_mya': 85,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 11',
        'alpha_pos': 0.26,
        'alpha_mb': '[32.2 Mb]',
        'alpha_anchors': 'NPRL3-MPG',
        'beta_chr': 'Chr 7',
        'beta_pos': 0.29,
        'beta_mb': '[103.5 Mb internal]',
        'beta_anchors': 'OR51-OR52',
        'notes': 'Telocentric chromosomes; karyotypic shuffling attenuates polarity',
    },
    # 20. Mammalia - Cattle
    {
        'species': 'Bos_taurus',
        'display': 'Cattle (B. taurus)',
        'clade': 'Artiodactyla',
        'divergence_mya': 60,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 25',
        'alpha_pos': 0.004,
        'alpha_mb': '[0.17 Mb subtelo]',
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_chr': 'Chr 15',
        'beta_pos': 0.43,
        'beta_mb': '[48.3 Mb internal]',
        'beta_anchors': 'OR51-OR52',
        'notes': 'Intense GC polarity: α subtelomeric (89.5%), β deep interstitial (65.5%)',
    },
    # 21. Mammalia - Sperm Whale
    {
        'species': 'physeter_catodon',
        'display': 'Sperm Whale (P. catodon)',
        'clade': 'Cetacea',
        'divergence_mya': 60,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'is_bichromosomal': True,
        'alpha_chr': 'Chr 15 scaf',
        'alpha_pos': 0.01,
        'alpha_mb': '[0.02 Mb subtelo]',
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_chr': 'Chr 15',
        'beta_pos': 0.42,
        'beta_mb': '[60.3 Mb internal]',
        'beta_anchors': 'TRIM21-FOLR1-OR51',
        'notes': 'Marine diving adaptation; high myoglobin, high α-globin GC',
    },
]


def load_data():
    """Load locus summaries, species histories, and gene decompositions."""
    loc_file = os.path.join(RESULTS_DIR, "locus_codon_and_synteny_summary.csv")
    sp_file = os.path.join(RESULTS_DIR, "species_locus_evolutionary_history.csv")
    node_file = os.path.join(RESULTS_DIR, "ancestral_node_reconstruction.csv")
    codon_file = os.path.join(RESULTS_DIR, "codon_decomposition_per_gene.csv")
    
    df_loc = pd.read_csv(loc_file)
    df_sp = pd.read_csv(sp_file)
    df_nodes = pd.read_csv(node_file)
    df_codons = pd.read_csv(codon_file)
    return df_loc, df_sp, df_nodes, df_codons


def plot_figure4_publication(df_loc, df_sp, df_nodes, df_codons, out_dir=RESULTS_DIR):
    """Generate high-DPI vector & raster publication Figure 4."""
    png_path = os.path.join(out_dir, "summary_figure4_phylogenetic_synteny.png")
    pdf_path = os.path.join(out_dir, "summary_figure4_phylogenetic_synteny.pdf")
    
    n_taxa = len(REPRESENTATIVE_TAXA)
    fig_h = max(18.0, n_taxa * 0.95)
    
    fig, (ax_tree, ax_syn, ax_codon) = plt.subplots(
        1, 3, figsize=(32, fig_h),
        gridspec_kw={'width_ratios': [4.7, 4.6, 4.2]},
        dpi=300
    )
    plt.subplots_adjust(wspace=0.10, left=0.03, right=0.97, top=0.94, bottom=0.06)
    
    y_positions = np.arange(n_taxa)[::-1]  # Top to bottom (index 0 at top, y=21)
    
    # -------------------------------------------------------------
    # PANEL 1: Strictly Orthogonal Chronogram Tree & Timeline
    # -------------------------------------------------------------
    ax_tree.set_title("A. Calibrated Vertebrate Tree & Evolutionary Milestones", fontsize=15, fontweight='bold', pad=18)
    ax_tree.set_xlim(670, -215)
    ax_tree.set_ylim(-2.2, n_taxa - 0.2)
    ax_tree.set_xlabel("Divergence Time (Million Years Ago - Ma)", fontsize=13, fontweight='bold', labelpad=10)
    
    # Time Gridlines (500 to 0 Ma)
    for t_mya in [500, 400, 300, 200, 100, 0]:
        ax_tree.axvline(t_mya, color="#ebebeb", linestyle="--", linewidth=1.0, zorder=1)
        ax_tree.text(t_mya, -0.6, f"{t_mya} Ma", fontsize=9.5, ha='center', va='top', color='#757575')

    def draw_orthogonal(t_parent, t_child, y_p, y_c, col, lw=2.4):
        """Draw standard orthogonal chronogram branch: vertical at t_parent, horizontal to t_child."""
        ax_tree.plot([t_parent, t_parent], [y_p, y_c], color=col, linewidth=lw, solid_capstyle='round', zorder=2)
        ax_tree.plot([t_parent, t_child], [y_c, y_c], color=col, linewidth=lw, solid_capstyle='round', zorder=2)

    # Calculate exact internal node heights (strictly bottom-up midpoint geometry)
    y_teleo = (y_positions[4] + y_positions[5]) / 2.0         # 16.5
    y_actino = (y_positions[3] + y_teleo) / 2.0               # 17.25
    y_chondro = (y_positions[1] + y_positions[2]) / 2.0       # 19.5
    y_aves = (y_positions[11] + y_positions[12]) / 2.0        # 9.5
    y_archo = (y_positions[10] + y_aves) / 2.0                # 10.25
    y_archosauro = (y_positions[9] + y_archo) / 2.0           # 11.125
    y_serpentes = (y_positions[14] + y_positions[15]) / 2.0   # 6.5
    y_squamata = (y_positions[13] + y_serpentes) / 2.0        # 7.25
    y_sauropsida = (y_archosauro + y_squamata) / 2.0          # 9.1875
    y_euarchonto = (y_positions[18] + y_positions[19]) / 2.0   # 2.5
    y_laurasia = (y_positions[20] + y_positions[21]) / 2.0     # 0.5
    y_placentalia = (y_euarchonto + y_laurasia) / 2.0         # 1.5
    y_theria = (y_positions[17] + y_placentalia) / 2.0         # 2.75
    y_mammalia = (y_positions[16] + y_theria) / 2.0           # 3.875
    y_amniota = (y_sauropsida + y_mammalia) / 2.0             # 6.531
    y_amphibia = (y_positions[7] + y_positions[8]) / 2.0       # 13.5
    y_tetrapoda = (y_amphibia + y_amniota) / 2.0              # 10.015
    y_sarco = (y_positions[6] + y_tetrapoda) / 2.0            # 12.508
    y_osteo = (y_actino + y_sarco) / 2.0                      # 14.879
    y_gnatho = (y_chondro + y_osteo) / 2.0                    # 17.189
    y_root = (y_positions[0] + y_gnatho) / 2.0                # 19.095

    # 1. Root / Vertebrata (500 Ma)
    ax_tree.plot([520, 500], [y_root, y_root], color='#424242', linewidth=2.6, zorder=2)
    draw_orthogonal(500, 0, y_root, y_positions[0], '#7b1fa2', lw=2.6) # Lamprey
    draw_orthogonal(500, 450, y_root, y_gnatho, '#424242', lw=2.6) # Gnathostome trunk

    # 2. Gnathostomata (450 Ma)
    draw_orthogonal(450, 435, y_gnatho, y_chondro, '#00838f', lw=2.6) # Chondrichthyes stem
    draw_orthogonal(450, 430, y_gnatho, y_osteo, '#1565c0', lw=2.6) # Osteichthyes stem

    # 3. Chondrichthyes (435 Ma)
    draw_orthogonal(435, 0, y_chondro, y_positions[1], '#00838f', lw=2.4) # Grey Reef Shark
    draw_orthogonal(435, 0, y_chondro, y_positions[2], '#00838f', lw=2.4) # Catshark

    # 4. Osteichthyes (430 Ma)
    draw_orthogonal(430, 350, y_osteo, y_actino, '#1565c0', lw=2.6) # Actinopterygii stem
    draw_orthogonal(430, 415, y_osteo, y_sarco, '#2e7d32', lw=2.6) # Sarcopterygii stem

    # 5. Actinopterygii (350 Ma)
    draw_orthogonal(350, 0, y_actino, y_positions[3], '#1565c0', lw=2.4) # Spotted Gar
    draw_orthogonal(350, 310, y_actino, y_teleo, '#1976d2', lw=2.6) # Teleost trunk
    draw_orthogonal(310, 0, y_teleo, y_positions[4], '#1976d2', lw=2.4) # Zebrafish
    draw_orthogonal(310, 0, y_teleo, y_positions[5], '#1976d2', lw=2.4) # Fugu

    # 6. Sarcopterygii (415 Ma)
    draw_orthogonal(415, 0, y_sarco, y_positions[6], '#2e7d32', lw=2.4) # Coelacanth
    draw_orthogonal(415, 355, y_sarco, y_tetrapoda, '#388e3c', lw=2.6) # Tetrapoda trunk

    # 7. Tetrapoda (355 Ma)
    draw_orthogonal(355, 270, y_tetrapoda, y_amphibia, '#388e3c', lw=2.6) # Amphibia trunk
    draw_orthogonal(355, 318, y_tetrapoda, y_amniota, '#616161', lw=2.6) # Amniota trunk
    draw_orthogonal(270, 0, y_amphibia, y_positions[7], '#388e3c', lw=2.4) # Xenopus
    draw_orthogonal(270, 0, y_amphibia, y_positions[8], '#388e3c', lw=2.4) # Axolotl

    # 8. Amniota (318 Ma) - THE AUTOSOMAL SPLIT
    draw_orthogonal(318, 280, y_amniota, y_sauropsida, '#00796b', lw=2.6) # Sauropsida stem
    draw_orthogonal(318, 165, y_amniota, y_mammalia, '#e65100', lw=2.6) # Mammalia stem

    # 9. Sauropsida (280 Ma)
    draw_orthogonal(280, 245, y_sauropsida, y_archosauro, '#00796b', lw=2.6) # Archosauromorpha stem
    draw_orthogonal(280, 165, y_sauropsida, y_squamata, '#d32f2f', lw=2.6) # Squamata stem

    # 10. Archosauromorpha (245 Ma)
    draw_orthogonal(245, 0, y_archosauro, y_positions[9], '#00796b', lw=2.4) # Green Sea Turtle
    draw_orthogonal(245, 240, y_archosauro, y_archo, '#00796b', lw=2.6) # Archosauria trunk
    draw_orthogonal(240, 0, y_archo, y_positions[10], '#00796b', lw=2.4) # American Alligator
    draw_orthogonal(240, 100, y_archo, y_aves, '#00796b', lw=2.6) # Aves trunk
    draw_orthogonal(100, 0, y_aves, y_positions[11], '#00796b', lw=2.4) # Chicken
    draw_orthogonal(100, 0, y_aves, y_positions[12], '#00796b', lw=2.4) # Hummingbird

    # 11. Squamata (165 Ma)
    draw_orthogonal(165, 0, y_squamata, y_positions[13], '#d32f2f', lw=2.4) # Wall Lizard
    draw_orthogonal(165, 110, y_squamata, y_serpentes, '#c2185b', lw=2.6) # Serpentes trunk
    draw_orthogonal(110, 0, y_serpentes, y_positions[14], '#c2185b', lw=2.4) # Corn Snake
    draw_orthogonal(110, 0, y_serpentes, y_positions[15], '#c2185b', lw=2.4) # Common Adder

    # 12. Mammalia (165 Ma crown)
    draw_orthogonal(165, 0, y_mammalia, y_positions[16], '#f57c00', lw=2.4) # Platypus
    draw_orthogonal(165, 140, y_mammalia, y_theria, '#e65100', lw=2.6) # Theria trunk

    # 13. Theria (140 Ma)
    draw_orthogonal(140, 0, y_theria, y_positions[17], '#f57c00', lw=2.4) # Gray Opossum
    draw_orthogonal(140, 95, y_theria, y_placentalia, '#e65100', lw=2.6) # Placentalia trunk

    # 14. Placentalia (95 Ma)
    draw_orthogonal(95, 85, y_placentalia, y_euarchonto, '#e65100', lw=2.6) # Euarchontoglires
    draw_orthogonal(95, 60, y_placentalia, y_laurasia, '#e65100', lw=2.6) # Laurasiatheria

    # Euarchontoglires (85 Ma)
    draw_orthogonal(85, 0, y_euarchonto, y_positions[18], '#e65100', lw=2.4) # Human
    draw_orthogonal(85, 0, y_euarchonto, y_positions[19], '#e65100', lw=2.4) # Mouse

    # Laurasiatheria (60 Ma)
    draw_orthogonal(60, 0, y_laurasia, y_positions[20], '#e65100', lw=2.4) # Cattle
    draw_orthogonal(60, 0, y_laurasia, y_positions[21], '#e65100', lw=2.4) # Sperm Whale

    # Tips and Taxon Display Names
    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        col = t['color']
        ax_tree.scatter(0, y, color=col, s=48, edgecolors='black', linewidths=0.9, zorder=5)
        ax_tree.text(-12, y, t['display'], fontsize=9.8, fontweight='bold', va='center', ha='left', color='#212121', zorder=5)

    # Milestone Badges on Tree Nodes
    milestones = [
        (500, y_positions[0], "N1", "500 Ma: Proto-Globin Single Locus (Agnatha)"),
        (450, y_gnatho, "N2", "450 Ma: Gnathostome α-β Tandem Duplication"),
        (435, y_chondro, "N3", "435 Ma: Chondrichthyan Telomeric Segregation"),
        (310, y_teleo, "N4", "310 Ma: Teleost 3R Duplication & Tandem Retention"),
        (318, y_amniota, "N5", "318 Ma: Amniote Autosomal Split (α: NPRL3, β translocates to CDCA8-OR51)"),
        (245, y_archosauro, "N6", "245 Ma: Archosaur Microchromosomal Partitioning"),
        (220, y_mammalia, "N7", "220 Ma: Mammalian Subtelomeric-α Capture"),
        (165, y_squamata, "N8", "165 Ma: Squamate Subtelomeric-β Inversion"),
    ]
    
    for mx, my, tag, desc in milestones:
        bbox_props = dict(boxstyle="circle,pad=0.26", fc="#fff9c4", ec="#f57f17", lw=1.3)
        ax_tree.text(mx, my, tag, fontsize=9.5, fontweight='bold', ha='center', va='center', bbox=bbox_props, zorder=6)

    # Clade Label Banners along the left margin
    clade_spans = [
        ("Agnatha", y_positions[0], y_positions[0], "#7b1fa2"),
        ("Chondrichthyes", y_positions[1], y_positions[2], "#00838f"),
        ("Actinopterygii", y_positions[3], y_positions[5], "#1565c0"),
        ("Sarcopterygii", y_positions[6], y_positions[6], "#2e7d32"),
        ("Amphibia (Linked)", y_positions[7], y_positions[8], "#388e3c"),
        ("Archosauromorpha", y_positions[9], y_positions[12], "#00796b"),
        ("Squamata (Inverted)", y_positions[13], y_positions[15], "#d32f2f"),
        ("Mammalia (Canonical)", y_positions[16], y_positions[21], "#e65100"),
    ]
    
    for cname, y1, y2, col in clade_spans:
        ymid = (y1 + y2) / 2.0
        ax_tree.plot([650, 650], [y1 - 0.28, y2 + 0.28], color=col, linewidth=4.5, solid_capstyle='round')
        ax_tree.text(640, ymid, cname, fontsize=10.5, fontweight='bold', ha='left', va='center', color=col, rotation=0)

    ax_tree.set_yticks([])
    ax_tree.grid(False)
    ax_tree.spines['top'].set_visible(False)
    ax_tree.spines['right'].set_visible(False)
    ax_tree.spines['left'].set_visible(False)

    # -------------------------------------------------------------
    # PANEL 2: Conserved Syntenic Blocks & Chromosomal Architecture
    # -------------------------------------------------------------
    ax_syn.set_title("B. Conserved Syntenic Blocks & Chromosomal Architecture", fontsize=15, fontweight='bold', pad=18)
    ax_syn.set_xlim(-0.04, 1.04)
    ax_syn.set_ylim(-2.2, n_taxa - 0.2)
    ax_syn.set_xlabel("Normalized Chromosome Coordinate (0.0 = p-Telomere, 0.5 = Centromere, 1.0 = q-Telomere)", fontsize=13, fontweight='bold', labelpad=10)
    
    # Telomere & Centromere guide markers
    ax_syn.axvline(0.12, color="#bdbdbd", linestyle=":", linewidth=1.0)
    ax_syn.axvline(0.55, color="#e0e0e0", linestyle=":", linewidth=0.8)
    ax_syn.axvline(0.98, color="#bdbdbd", linestyle=":", linewidth=1.0)
    ax_syn.text(0.12, -0.6, "p-Telomere (0.0)", fontsize=8.5, ha='center', va='top', color='#757575', fontweight='bold')
    ax_syn.text(0.55, -0.6, "Centromere (~0.5)", fontsize=8.5, ha='center', va='top', color='#757575', fontweight='bold')
    ax_syn.text(0.98, -0.6, "q-Telomere (1.0)", fontsize=8.5, ha='center', va='top', color='#757575', fontweight='bold')

    # Function to draw a clean chromosome ribbon with centromeric constriction
    def draw_chromosome_ribbon(ax, x0, x1, y_mid, height, chr_name, fc="#f8f9fa", ec="#cccccc", cen_rel=0.48):
        width = x1 - x0
        cen_x = x0 + cen_rel * width
        
        # Left arm
        left_arm = FancyBboxPatch((x0, y_mid - height/2.0), cen_x - x0 - 0.015, height,
                                  boxstyle="round,pad=0.005,rounding_size=0.03",
                                  fc=fc, ec=ec, lw=1.1, zorder=2)
        ax.add_patch(left_arm)
        
        # Right arm
        right_arm = FancyBboxPatch((cen_x + 0.015, y_mid - height/2.0), x1 - cen_x - 0.015, height,
                                   boxstyle="round,pad=0.005,rounding_size=0.03",
                                   fc=fc, ec=ec, lw=1.1, zorder=2)
        ax.add_patch(right_arm)
        
        # Centromere constriction
        ax.plot([cen_x - 0.01, cen_x + 0.01], [y_mid, y_mid], color="#9e9e9e", lw=2.0, zorder=3)
        ax.scatter(cen_x, y_mid, s=20, marker='o', fc='#bdbdbd', ec='#757575', lw=0.6, zorder=4)
        
        # Chromosome label on the far left
        ax.text(x0 - 0.012, y_mid, chr_name, fontsize=8.0, fontweight='bold', ha='right', va='center',
                color='#37474f', bbox=dict(boxstyle="round,pad=0.2", fc="#eceff1", ec="#b0bec5", lw=0.8), zorder=5)

    X_START = 0.12
    X_END = 0.98
    X_SPAN = X_END - X_START

    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        is_bi = t['is_bichromosomal']
        
        if is_bi:
            # -------------------------------------------------------------
            # BICHROMOSOMAL: Two distinct ribbons (Top = α-Chr, Bottom = β-Chr)
            # -------------------------------------------------------------
            y_alpha = y + 0.17
            y_beta = y - 0.17
            
            # 1. Top Ribbon: Alpha Chromosome
            draw_chromosome_ribbon(ax_syn, X_START, X_END, y_alpha, height=0.15, chr_name=t['alpha_chr'])
            a_rel = np.clip(t['alpha_pos'], 0.01, 0.94)
            a_x = X_START + a_rel * X_SPAN
            ax_syn.add_patch(FancyBboxPatch((a_x, y_alpha - 0.065), 0.045, 0.13,
                                            boxstyle="round,pad=0.005,rounding_size=0.02",
                                            fc="#1e88e5", ec="#0d47a1", lw=1.2, zorder=6))
            ax_syn.text(a_x + 0.0225, y_alpha, "α", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=7)
            # Label placement: left if in interior, right if near p-telomere
            if a_rel >= 0.28:
                ax_syn.text(a_x - 0.012, y_alpha, f"{t['alpha_anchors']} {t.get('alpha_mb', '')}",
                            fontsize=7.2, va='center', ha='right', color="#0d47a1", zorder=7)
            else:
                ax_syn.text(a_x + 0.052, y_alpha, f"{t['alpha_anchors']} {t.get('alpha_mb', '')}",
                            fontsize=7.2, va='center', ha='left', color="#0d47a1", zorder=7)
            
            # 2. Bottom Ribbon: Beta Chromosome
            draw_chromosome_ribbon(ax_syn, X_START, X_END, y_beta, height=0.15, chr_name=t['beta_chr'])
            b_rel = np.clip(t['beta_pos'], 0.01, 0.94)
            b_x = X_START + b_rel * X_SPAN
            b_fc = "#d32f2f" if t['state'] == 'SQUAMATE_BETA_INVERTED' else "#43a047"
            b_ec = "#b71c1c" if t['state'] == 'SQUAMATE_BETA_INVERTED' else "#1b5e20"
            b_txt_col = "#b71c1c" if t['state'] == 'SQUAMATE_BETA_INVERTED' else "#1b5e20"
            
            ax_syn.add_patch(FancyBboxPatch((b_x, y_beta - 0.065), 0.045, 0.13,
                                            boxstyle="round,pad=0.005,rounding_size=0.02",
                                            fc=b_fc, ec=b_ec, lw=1.2, zorder=6))
            ax_syn.text(b_x + 0.0225, y_beta, "β", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=7)
            # Label placement: left if in interior, right if near p-telomere
            if b_rel >= 0.28:
                ax_syn.text(b_x - 0.012, y_beta, f"{t['beta_anchors']} {t.get('beta_mb', '')}",
                            fontsize=7.2, va='center', ha='right', color=b_txt_col, zorder=7)
            else:
                ax_syn.text(b_x + 0.052, y_beta, f"{t['beta_anchors']} {t.get('beta_mb', '')}",
                            fontsize=7.2, va='center', ha='left', color=b_txt_col, zorder=7)
            
        else:
            # -------------------------------------------------------------
            # SINGLE CHROMOSOME: Tandem linked or ancestral proto locus
            # -------------------------------------------------------------
            draw_chromosome_ribbon(ax_syn, X_START, X_END, y, height=0.24, chr_name=t['chr_label'])
            pos_rel = np.clip(t['pos'], 0.01, 0.90)
            x_pos = X_START + pos_rel * X_SPAN
            
            if t['state'] == 'CYCLOSTOME_PROTO':
                ax_syn.add_patch(FancyBboxPatch((x_pos, y - 0.09), 0.065, 0.18,
                                                boxstyle="round,pad=0.005,rounding_size=0.02",
                                                fc="#ab47bc", ec="#4a148c", lw=1.2, zorder=6))
                ax_syn.text(x_pos + 0.0325, y, "Proto-Hb", fontsize=8.0, fontweight='bold', ha='center', va='center', color='white', zorder=7)
                ax_syn.text(x_pos + 0.075, y, f"{t['anchors']} {t['mb_range']}", fontsize=7.5, va='center', ha='left', color="#4a148c", zorder=7)
            else:
                # Tandem linked α-β block
                ax_syn.add_patch(FancyBboxPatch((x_pos, y - 0.09), 0.045, 0.18,
                                                boxstyle="round,pad=0.005,rounding_size=0.02",
                                                fc="#1e88e5", ec="#0d47a1", lw=1.2, zorder=6))
                ax_syn.add_patch(FancyBboxPatch((x_pos + 0.045, y - 0.09), 0.045, 0.18,
                                                boxstyle="round,pad=0.005,rounding_size=0.02",
                                                fc="#43a047", ec="#1b5e20", lw=1.2, zorder=6))
                ax_syn.text(x_pos + 0.0225, y, "α", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=7)
                ax_syn.text(x_pos + 0.0675, y, "β", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=7)
                if pos_rel >= 0.32:
                    ax_syn.text(x_pos - 0.012, y, f"{t['anchors']} {t['mb_range']}", fontsize=7.5, va='center', ha='right', color="#1565c0", zorder=7)
                else:
                    ax_syn.text(x_pos + 0.10, y, f"{t['anchors']} {t['mb_range']}", fontsize=7.5, va='center', ha='left', color="#1565c0", zorder=7)

    ax_syn.set_yticks([])
    ax_syn.spines['top'].set_visible(False)
    ax_syn.spines['right'].set_visible(False)
    ax_syn.spines['left'].set_visible(False)

    # -------------------------------------------------------------
    # PANEL 3: Codon Position Decomposition & Diverging GC3 Contrast
    # -------------------------------------------------------------
    ax_codon.set_title("C. Codon Decomposition (GC1, GC2, GC3) & Subunit Divergence (ΔGC3)", fontsize=15, fontweight='bold', pad=18)
    ax_codon.set_xlim(30, 118)
    ax_codon.set_ylim(-2.2, n_taxa - 0.2)
    ax_codon.set_xlabel("GC Content (%) [Corridor: Conserved GC1 & GC2 | Bar: ΔGC3 Subunit Contrast]", fontsize=13, fontweight='bold', labelpad=10)
    ax_codon.set_xticks([40, 50, 60, 70, 80, 90, 100])
    
    # Purifying Selection Corridor (45% to 55%) across all species
    ax_codon.axvspan(45, 55, color='#eceff1', alpha=0.75, zorder=1)
    ax_codon.axvline(45, color='#cfd8dc', linestyle=':', linewidth=1.1, zorder=2)
    ax_codon.axvline(55, color='#cfd8dc', linestyle=':', linewidth=1.1, zorder=2)
    ax_codon.text(50, -0.6, "Purifying Selection (GC1/2: 45–55%)", fontsize=8.0, ha='center', va='top', color='#546e7a', fontweight='bold')

    # Background guides
    for gc_mark in [40, 60, 70, 80, 90, 100]:
        ax_codon.axvline(gc_mark, color="#f5f5f5", linestyle="--", linewidth=0.9, zorder=1)

    # Continuous Diverging Colormap (coolwarm: Blue for β>α, Grey for Balanced, Red for α>β)
    div_norm = mcolors.Normalize(vmin=-25.0, vmax=25.0)
    coolwarm_cmap = plt.cm.coolwarm

    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        sp = t['species']
        sp_rows = df_sp[df_sp['species'] == sp]
        loc_rows = df_loc[df_loc['species'] == sp]
        
        # Retrieve codon metrics
        if t['state'] == 'LINKED_ANCESTRAL':
            sp_g = df_codons[(df_codons['species'] == sp) & (df_codons['is_globin'])]
            a_genes = sp_g[sp_g['gene_name'].str.contains(r'hba|hbz|alpha|hbaa|hbae|LOC102354|LOC138569', case=False, regex=True)]
            b_genes = sp_g[sp_g['gene_name'].str.contains(r'hbb|hbe|hbg|beta|hbba|hbbe|LOC102356|LOC102355|LOC138579|LOC138580', case=False, regex=True)]
            a_gc1 = a_genes['gc1'].mean() if not a_genes.empty and 'gc1' in a_genes and pd.notna(a_genes['gc1'].mean()) else 52.0
            a_gc2 = a_genes['gc2'].mean() if not a_genes.empty and 'gc2' in a_genes and pd.notna(a_genes['gc2'].mean()) else 42.0
            a_gc3 = a_genes['gc3'].mean() if not a_genes.empty and 'gc3' in a_genes and pd.notna(a_genes['gc3'].mean()) else 60.0
            b_gc1 = b_genes['gc1'].mean() if not b_genes.empty and 'gc1' in b_genes and pd.notna(b_genes['gc1'].mean()) else 50.0
            b_gc2 = b_genes['gc2'].mean() if not b_genes.empty and 'gc2' in b_genes and pd.notna(b_genes['gc2'].mean()) else 38.0
            b_gc3 = b_genes['gc3'].mean() if not b_genes.empty and 'gc3' in b_genes and pd.notna(b_genes['gc3'].mean()) else 60.0
            has_a = not a_genes.empty
            has_b = not b_genes.empty
        else:
            a_row = loc_rows[loc_rows['locus_class'].isin(['alpha_globin', 'linked_alpha_beta'])]
            b_row = loc_rows[loc_rows['locus_class'].isin(['beta_globin', 'linked_alpha_beta'])]
            has_a = not a_row.empty
            has_b = not b_row.empty
            a_gc1 = a_row['mean_gc1'].iloc[0] if has_a and pd.notna(a_row['mean_gc1'].iloc[0]) else 52.0
            a_gc2 = a_row['mean_gc2'].iloc[0] if has_a and pd.notna(a_row['mean_gc2'].iloc[0]) else 42.0
            a_gc3 = a_row['mean_gc3'].iloc[0] if has_a and pd.notna(a_row['mean_gc3'].iloc[0]) else 65.0
            b_gc1 = b_row['mean_gc1'].iloc[0] if has_b and pd.notna(b_row['mean_gc1'].iloc[0]) else 50.0
            b_gc2 = b_row['mean_gc2'].iloc[0] if has_b and pd.notna(b_row['mean_gc2'].iloc[0]) else 38.0
            b_gc3 = b_row['mean_gc3'].iloc[0] if has_b and pd.notna(b_row['mean_gc3'].iloc[0]) else 58.0
        
        # Calculate Delta GC3
        d_gc3 = (a_gc3 - b_gc3) if (has_a and has_b and t['state'] != 'CYCLOSTOME_PROTO') else 0.0
        div_color = coolwarm_cmap(div_norm(d_gc3))

        # 1. Purifying selection markers (GC1 and GC2)
        if has_a:
            ax_codon.scatter(a_gc1, y + 0.14, marker='o', s=30, color='#1565c0', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(a_gc2, y + 0.14, marker='^', s=30, color='#42a5f5', edgecolors='black', linewidths=0.6, zorder=4)

        if has_b and t['state'] != 'CYCLOSTOME_PROTO':
            ax_codon.scatter(b_gc1, y - 0.14, marker='s', s=30, color='#2e7d32', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(b_gc2, y - 0.14, marker='v', s=30, color='#81c784', edgecolors='black', linewidths=0.6, zorder=4)

        # 2. Subunit Divergence Bar connecting alpha-GC3 and beta-GC3
        if has_a and has_b and t['state'] != 'CYCLOSTOME_PROTO':
            min_gc3 = min(a_gc3, b_gc3)
            max_gc3 = max(a_gc3, b_gc3)
            # Chunky rounded contrast bar
            ax_codon.plot([min_gc3, max_gc3], [y, y],
                          color=div_color, linewidth=5.5, zorder=3, solid_capstyle='round')
            
            # Alpha GC3 marker (Blue Diamond)
            ax_codon.scatter(a_gc3, y, marker='D', s=70, color='#1565c0', edgecolors='black', linewidths=0.9, zorder=5)
            # Beta GC3 marker (Green or Crimson Square)
            b_mcol = '#d32f2f' if t['state'] == 'SQUAMATE_BETA_INVERTED' else '#2e7d32'
            ax_codon.scatter(b_gc3, y, marker='s', s=70, color=b_mcol, edgecolors='black', linewidths=0.9, zorder=5)
            
        elif t['state'] == 'CYCLOSTOME_PROTO':
            # Single proto-globin
            ax_codon.scatter(a_gc3, y, marker='D', s=75, color='#ab47bc', edgecolors='black', linewidths=1.0, zorder=5)

        # 3. Delta GC3 Badge on the right
        if t['state'] == 'CYCLOSTOME_PROTO':
            badge_txt = "Proto-Hb"
            badge_bg = "#f3e5f5"
            badge_ec = "#ab47bc"
            badge_col = "#6a1b9a"
        elif d_gc3 > 5.0:
            badge_txt = f"ΔGC3: +{d_gc3:.1f}% (α > β)"
            badge_bg = "#fff3e0" if d_gc3 < 15 else "#ffebee"
            badge_ec = div_color
            badge_col = "#b71c1c"
        elif d_gc3 < -5.0:
            badge_txt = f"ΔGC3: {d_gc3:.1f}% (β > α)"
            badge_bg = "#e8eaf6"
            badge_ec = div_color
            badge_col = "#1a237e"
        else:
            badge_txt = f"ΔGC3: {d_gc3:+.1f}% (Bal)"
            badge_bg = "#f5f5f5"
            badge_ec = "#9e9e9e"
            badge_col = "#424242"
            
        ax_codon.text(103.5, y, badge_txt, fontsize=8.2, fontweight='bold', va='center', ha='left',
                      color=badge_col, bbox=dict(boxstyle="round,pad=0.28", fc=badge_bg, ec=badge_ec, lw=1.1))

    # Dedicated Legend below the plot (absolutely zero overlap with data points)
    leg_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1565c0', markeredgecolor='black', markersize=6.5, label='α-Globin GC1'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='#42a5f5', markeredgecolor='black', markersize=6.5, label='α-Globin GC2'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#1565c0', markeredgecolor='black', markersize=7.5, label='α-Globin GC3'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#2e7d32', markeredgecolor='black', markersize=6.5, label='β-Globin GC1'),
        Line2D([0], [0], marker='v', color='w', markerfacecolor='#81c784', markeredgecolor='black', markersize=6.5, label='β-Globin GC2'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#2e7d32', markeredgecolor='black', markersize=7.5, label='β-Globin GC3 (Canonical/Linked)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#d32f2f', markeredgecolor='black', markersize=7.5, label='β-Globin GC3 (Squamate Inverted)'),
        Line2D([0], [0], color='#b40426', linewidth=4.0, label='ΔGC3 > 0 Divergence Bar (α-dominant)'),
        Line2D([0], [0], color='#dddddd', linewidth=4.0, label='ΔGC3 ≈ 0 Divergence Bar (Balanced)'),
        Line2D([0], [0], color='#3b4cc0', linewidth=4.0, label='ΔGC3 < 0 Divergence Bar (β-dominant)'),
    ]
    ax_codon.legend(handles=leg_elements, loc='upper center', bbox_to_anchor=(0.5, -0.06),
                    fontsize=8.5, frameon=True, facecolor='#ffffff', edgecolor='#cccccc', ncol=4)

    ax_codon.set_yticks([])
    ax_codon.spines['top'].set_visible(False)
    ax_codon.spines['right'].set_visible(False)
    ax_codon.spines['left'].set_visible(False)

    plt.suptitle("Figure 4: Macroevolutionary Synteny Architecture and Codon Divergence Across Vertebrates",
                 fontsize=18, fontweight='bold', y=0.975)

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 4 (PNG/PDF): {png_path} and {pdf_path}")
    return png_path, pdf_path


def plot_figure4_interactive(df_loc, df_sp, df_nodes, df_codons, out_dir=RESULTS_DIR):
    """Generate high-fidelity interactive Plotly HTML for Figure 4."""
    html_path = os.path.join(out_dir, "summary_figure4_phylogenetic_synteny.html")
    
    n_taxa = len(REPRESENTATIVE_TAXA)
    y_positions = list(range(n_taxa))[::-1]
    
    fig = make_subplots(
        rows=1, cols=3,
        shared_yaxes=True,
        column_widths=[0.32, 0.38, 0.30],
        subplot_titles=[
            "<b>A. Calibrated Vertebrate Chronogram & Transitions</b>",
            "<b>B. Conserved Syntenic Architecture (Separate Ribbons for Separate Autosomes)</b>",
            "<b>C. Codon Decomposition & Diverging ΔGC3 Contrast</b>"
        ],
        horizontal_spacing=0.06
    )
    
    # 1. Panel A: Tree Trunks and Species Tips
    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        col = t['color']
        div_t = t['divergence_mya']
        
        # Horizontal branch
        fig.add_trace(go.Scatter(
            x=[div_t, 0], y=[y, y],
            mode='lines',
            line=dict(color=col, width=3.5),
            hoverinfo='skip',
            showlegend=False
        ), row=1, col=1)
        
        # Tip marker with tooltip
        tip_text = (
            f"<b>{t['display']}</b><br>"
            f"Clade: {t['clade']}<br>"
            f"Divergence: {div_t} Ma<br>"
            f"State: {t['state']}<br>"
            f"<i>{t['notes']}</i>"
        )
        fig.add_trace(go.Scatter(
            x=[0], y=[y],
            mode='markers+text',
            marker=dict(size=11, color=col, line=dict(color='black', width=1)),
            text=[f"  {t['display']}"],
            textposition="middle right",
            hovertext=tip_text,
            hoverinfo='text',
            showlegend=False
        ), row=1, col=1)

    # Milestones in Panel A
    milestones = [
        (500, 21.0, "N1", "500 Ma: Proto-Globin Single Locus (Agnatha)"),
        (450, 17.2, "N2", "450 Ma: Gnathostome α-β Tandem Duplication"),
        (435, 19.5, "N3", "435 Ma: Chondrichthyan Telomeric Segregation"),
        (310, 16.5, "N4", "310 Ma: Teleost 3R Duplication & Tandem Retention"),
        (318, 6.5, "N5", "318 Ma: Amniote Autosomal Split (α: NPRL3, β translocates to CDCA8-OR51)"),
        (245, 11.1, "N6", "245 Ma: Archosaur Microchromosomal Partitioning"),
        (220, 3.9, "N7", "220 Ma: Mammalian Subtelomeric-α Capture"),
        (165, 7.3, "N8", "165 Ma: Squamate Subtelomeric-β Inversion"),
    ]
    for mx, my, tag, desc in milestones:
        fig.add_trace(go.Scatter(
            x=[mx], y=[my],
            mode='markers+text',
            marker=dict(size=14, color='#fff9c4', line=dict(color='#f57f17', width=2)),
            text=[tag],
            textposition="middle center",
            hovertext=f"<b>Milestone {tag}</b><br>{desc}",
            hoverinfo='text',
            showlegend=False
        ), row=1, col=1)

    # 2. Panel B: Syntenic Chromosome Schematics
    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        is_bi = t['is_bichromosomal']
        
        if is_bi:
            # Top Ribbon: Alpha Chromosome
            y_a = y + 0.17
            fig.add_trace(go.Scatter(
                x=[0.12, 0.98], y=[y_a, y_a],
                mode='lines',
                line=dict(color='#e0e0e0', width=10),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=2)
            
            a_pos = 0.12 + np.clip(t['alpha_pos'], 0.01, 0.94) * 0.86
            fig.add_trace(go.Scatter(
                x=[a_pos], y=[y_a],
                mode='markers+text',
                marker=dict(size=12, color='#1e88e5', symbol='square', line=dict(color='#0d47a1', width=1)),
                text=[f"  α: {t['alpha_chr']} {t.get('alpha_mb', '')}"],
                textposition="middle right",
                hovertext=f"<b>α-Globin</b><br>Chr: {t['alpha_chr']}<br>Pos: {t.get('alpha_mb', '')}<br>Anchors: {t['alpha_anchors']}",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
            
            # Bottom Ribbon: Beta Chromosome
            y_b = y - 0.17
            fig.add_trace(go.Scatter(
                x=[0.12, 0.98], y=[y_b, y_b],
                mode='lines',
                line=dict(color='#e0e0e0', width=10),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=2)
            
            b_pos = 0.12 + np.clip(t['beta_pos'], 0.01, 0.94) * 0.86
            b_col = '#d32f2f' if t['state'] == 'SQUAMATE_BETA_INVERTED' else '#43a047'
            fig.add_trace(go.Scatter(
                x=[b_pos], y=[y_b],
                mode='markers+text',
                marker=dict(size=12, color=b_col, symbol='square', line=dict(color='black', width=1)),
                text=[f"  β: {t['beta_chr']} {t.get('beta_mb', '')}"],
                textposition="middle right",
                hovertext=f"<b>β-Globin</b><br>Chr: {t['beta_chr']}<br>Pos: {t.get('beta_mb', '')}<br>Anchors: {t['beta_anchors']}",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
            
        else:
            # Single Ribbon: Linked or Proto
            fig.add_trace(go.Scatter(
                x=[0.12, 0.98], y=[y, y],
                mode='lines',
                line=dict(color='#e0e0e0', width=14),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=2)
            
            p_pos = 0.12 + np.clip(t['pos'], 0.01, 0.90) * 0.86
            if t['state'] == 'CYCLOSTOME_PROTO':
                fig.add_trace(go.Scatter(
                    x=[p_pos], y=[y],
                    mode='markers+text',
                    marker=dict(size=14, color='#ab47bc', symbol='square'),
                    text=[f"  Proto-Hb ({t['chr_label']})"],
                    textposition="middle right",
                    hovertext=f"Ancestral Proto-Globin: {t['chr_label']} {t['mb_range']}",
                    hoverinfo='text',
                    showlegend=False
                ), row=1, col=2)
            else:
                fig.add_trace(go.Scatter(
                    x=[p_pos], y=[y],
                    mode='markers+text',
                    marker=dict(size=14, color='#1e88e5', symbol='square'),
                    text=[f"  α-β linked ({t['chr_label']}) {t['mb_range']}"],
                    textposition="middle right",
                    hovertext=f"<b>Linked α-β Locus</b><br>Chr: {t['chr_label']}<br>Pos: {t['mb_range']}<br>Anchors: {t['anchors']}",
                    hoverinfo='text',
                    showlegend=False
                ), row=1, col=2)

    # 3. Panel C: Codon Decomposition & Diverging Contrast
    div_norm = mcolors.Normalize(vmin=-25.0, vmax=25.0)
    coolwarm_cmap = plt.cm.coolwarm

    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        sp = t['species']
        loc_rows = df_loc[df_loc['species'] == sp]
        
        if t['state'] == 'LINKED_ANCESTRAL':
            sp_g = df_codons[(df_codons['species'] == sp) & (df_codons['is_globin'])]
            a_genes = sp_g[sp_g['gene_name'].str.contains(r'hba|hbz|alpha|hbaa|hbae|LOC102354|LOC138569', case=False, regex=True)]
            b_genes = sp_g[sp_g['gene_name'].str.contains(r'hbb|hbe|hbg|beta|hbba|hbbe|LOC102356|LOC102355|LOC138579|LOC138580', case=False, regex=True)]
            a_gc3 = a_genes['gc3'].mean() if not a_genes.empty and 'gc3' in a_genes and pd.notna(a_genes['gc3'].mean()) else 60.0
            b_gc3 = b_genes['gc3'].mean() if not b_genes.empty and 'gc3' in b_genes and pd.notna(b_genes['gc3'].mean()) else 60.0
            has_a = not a_genes.empty
            has_b = not b_genes.empty
        else:
            a_row = loc_rows[loc_rows['locus_class'].isin(['alpha_globin', 'linked_alpha_beta'])]
            b_row = loc_rows[loc_rows['locus_class'].isin(['beta_globin', 'linked_alpha_beta'])]
            has_a = not a_row.empty
            has_b = not b_row.empty
            a_gc3 = a_row['mean_gc3'].iloc[0] if has_a and pd.notna(a_row['mean_gc3'].iloc[0]) else 65.0
            b_gc3 = b_row['mean_gc3'].iloc[0] if has_b and pd.notna(b_row['mean_gc3'].iloc[0]) else 58.0
        
        d_gc3 = (a_gc3 - b_gc3) if (has_a and has_b and t['state'] != 'CYCLOSTOME_PROTO') else 0.0
        rgba = coolwarm_cmap(div_norm(d_gc3))
        hex_col = mcolors.to_hex(rgba)
        
        if has_a and has_b and t['state'] != 'CYCLOSTOME_PROTO':
            # Subunit Divergence Bar
            fig.add_trace(go.Scatter(
                x=[min(a_gc3, b_gc3), max(a_gc3, b_gc3)], y=[y, y],
                mode='lines',
                line=dict(color=hex_col, width=8),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=3)
            
            # Alpha GC3 Diamond
            fig.add_trace(go.Scatter(
                x=[a_gc3], y=[y],
                mode='markers',
                marker=dict(size=12, color='#1565c0', symbol='diamond', line=dict(color='black', width=1)),
                hovertext=f"<b>{t['display']} (α-Globin)</b><br>GC3: {a_gc3:.1f}%<br>ΔGC3 (α - β): {d_gc3:+.1f}%",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=3)
            
            # Beta GC3 Square
            b_c = '#d32f2f' if t['state'] == 'SQUAMATE_BETA_INVERTED' else '#2e7d32'
            fig.add_trace(go.Scatter(
                x=[b_gc3], y=[y],
                mode='markers',
                marker=dict(size=12, color=b_c, symbol='square', line=dict(color='black', width=1)),
                hovertext=f"<b>{t['display']} (β-Globin)</b><br>GC3: {b_gc3:.1f}%<br>ΔGC3 (α - β): {d_gc3:+.1f}%",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=3)
            
        elif t['state'] == 'CYCLOSTOME_PROTO':
            fig.add_trace(go.Scatter(
                x=[a_gc3], y=[y],
                mode='markers',
                marker=dict(size=13, color='#ab47bc', symbol='diamond', line=dict(color='black', width=1)),
                hovertext=f"<b>{t['display']} (Proto-Hb)</b><br>GC3: {a_gc3:.1f}%",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=3)

    # Layout styling
    fig.update_xaxes(title_text="<b>Divergence Time (Ma)</b>", range=[560, -180], autorange="reversed", row=1, col=1)
    fig.update_xaxes(title_text="<b>Normalized Position (0 = p-Telomere, 1 = q-Telomere)</b>", range=[-0.05, 1.05], row=1, col=2)
    fig.update_xaxes(title_text="<b>GC Content (%) [Diamond = α-GC3, Square = β-GC3]</b>", range=[32, 115], row=1, col=3)
    
    fig.update_yaxes(showticklabels=False, range=[-1.5, n_taxa - 0.2])
    
    fig.update_layout(
        title="<b>Figure 4: Macroevolutionary Synteny Architecture and Codon Divergence Across Vertebrates</b>",
        height=max(900, n_taxa * 44),
        width=1800,
        template="plotly_white",
        margin=dict(l=40, r=40, t=80, b=50)
    )
    
    fig.write_html(html_path)
    print(f"Generated Figure 4 Interactive HTML: {html_path}")
    return html_path


def main():
    print("=== Generating Polished Figure 4: Separate Chromosomes & Diverging Codon Map ===")
    df_loc, df_sp, df_nodes, df_codons = load_data()
    
    png_path, pdf_path = plot_figure4_publication(df_loc, df_sp, df_nodes, df_codons)
    html_path = plot_figure4_interactive(df_loc, df_sp, df_nodes, df_codons)
    
    # Synchronize to 1_data_from_mike/results/
    target_dir = os.path.join(BASE_DIR, "1_data_from_mike/results")
    if os.path.exists(target_dir):
        for f in [png_path, pdf_path, html_path]:
            dest = os.path.join(target_dir, os.path.basename(f))
            shutil.copy2(f, dest)
            print(f"  Synchronized: {os.path.basename(f)}")
            
    print("\nFigure 4 generation and synchronization complete!")

if __name__ == "__main__":
    main()
