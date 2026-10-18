"""Generate Figure 4: Macroevolutionary Synteny and Codon Divergence (GC1/GC2/GC3) Map
across the Vertebrate Tree (Polished Publication Version).

Combines:
1. Time-calibrated phylogenetic tree (520 Ma -> 0 Ma) with orthogonal branch layout.
2. Syntenic chromosome schematics with cytogenetic centromeric constrictions,
   telomeres, and cleanly positioned syntenic anchor labels.
3. Paired codon position decomposition (GC1, GC2, GC3) illustrating the synonymous gBGC surge
   and the polarity reversal between Mammalia (alpha-dominant) and Squamata (beta-dominant).
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
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")

REPRESENTATIVE_TAXA = [
    # 0. Agnatha
    {
        'species': 'petromyzon_marinus',
        'display': 'Sea Lamprey (P. marinus)',
        'clade': 'Agnatha',
        'divergence_mya': 500,
        'state': 'CYCLOSTOME_PROTO',
        'color': '#7b1fa2',
        'alpha_chr': 'GL478636',
        'alpha_pos': 0.08,
        'beta_chr': None,
        'beta_pos': None,
        'alpha_anchors': 'Proto-Hb',
        'beta_anchors': '',
        'notes': 'Single ancestral globin locus; monomeric/homomeric',
    },
    # 1-2. Chondrichthyes
    {
        'species': 'carcharhinus_amblyrhynchos',
        'display': 'Grey Reef Shark (C. amblyrhynchos)',
        'clade': 'Chondrichthyes',
        'divergence_mya': 435,
        'state': 'CHONDRICHTHYAN_GRADIENT',
        'color': '#00838f',
        'alpha_chr': 'Chr 18',
        'alpha_pos': 0.06,
        'beta_chr': 'Chr 18',
        'beta_pos': 0.16,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': 'Flank',
        'notes': 'Separate loci; GC strictly governed by telomeric distance gradient',
    },
    {
        'species': 'scyliorhinus_canicula',
        'display': 'Catshark (S. canicula)',
        'clade': 'Chondrichthyes',
        'divergence_mya': 430,
        'state': 'CHONDRICHTHYAN_GRADIENT',
        'color': '#00838f',
        'alpha_chr': 'Chr 15',
        'alpha_pos': 0.04,
        'beta_chr': 'Chr 21',
        'beta_pos': 0.68,
        'alpha_anchors': 'NPRL3-RHBDF1',
        'beta_anchors': 'Flank',
        'notes': 'Bichromosomal; alpha subtelomeric, beta deep interstitial',
    },
    # 3. Basal Actinopterygii
    {
        'species': 'lepisosteus_oculatus',
        'display': 'Spotted Gar (L. oculatus)',
        'clade': 'Holostei',
        'divergence_mya': 350,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1565c0',
        'alpha_chr': 'LG13',
        'alpha_pos': 0.12,
        'beta_chr': 'LG13',
        'beta_pos': 0.12,
        'alpha_anchors': 'NPRL3-MPG-RHBDF1',
        'beta_anchors': '',
        'notes': 'Pre-teleost 3R state; single linked alpha-beta cluster on LG13',
    },
    # 4-5. Teleostei
    {
        'species': 'danio_rerio',
        'display': 'Zebrafish (D. rerio)',
        'clade': 'Teleostei',
        'divergence_mya': 310,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1976d2',
        'alpha_chr': 'Chr 3 & 12',
        'alpha_pos': 0.10,
        'beta_chr': 'Chr 3 & 12',
        'beta_pos': 0.10,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': '',
        'notes': '3R duplicate loci (LA Chr 3, LB Chr 12); tandem head-to-head pairs',
    },
    {
        'species': 'takifugu_rubripes',
        'display': 'Fugu (T. rubripes)',
        'clade': 'Teleostei',
        'divergence_mya': 310,
        'state': 'LINKED_ANCESTRAL',
        'color': '#1976d2',
        'alpha_chr': 'Chr 2',
        'alpha_pos': 0.18,
        'beta_chr': 'Chr 2',
        'beta_pos': 0.18,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': '',
        'notes': 'Ultra-compact genome; linked alpha-beta pairs retain parallel GC',
    },
    # 6. Sarcopterygii
    {
        'species': 'latimeria_chalumnae',
        'display': 'Coelacanth (L. chalumnae)',
        'clade': 'Sarcopterygii',
        'divergence_mya': 415,
        'state': 'LINKED_ANCESTRAL',
        'color': '#2e7d32',
        'alpha_chr': 'Chr 5',
        'alpha_pos': 0.14,
        'beta_chr': 'Chr 5',
        'beta_pos': 0.14,
        'alpha_anchors': 'NPRL3-LUC7L',
        'beta_anchors': '',
        'notes': 'Basal lobe-finned fish; retains primordial linked alpha-beta locus flanked by NPRL3-LUC7L',
    },
    # 7-8. Amphibia
    {
        'species': 'xenopus_tropicalis',
        'display': 'Clawed Frog (X. tropicalis)',
        'clade': 'Amphibia',
        'divergence_mya': 270,
        'state': 'LINKED_ANCESTRAL',
        'color': '#388e3c',
        'alpha_chr': 'Chr 9',
        'alpha_pos': 0.10,
        'beta_chr': 'Chr 9',
        'beta_pos': 0.10,
        'alpha_anchors': 'NPRL3-RHBDF1-LUC7L',
        'beta_anchors': '',
        'notes': 'Non-amniote tetrapod; retains single linked alpha-beta locus on Chr 9 (hba + hbg/hbe)',
    },
    {
        'species': 'ambystoma_mexicanum',
        'display': 'Axolotl (A. mexicanum)',
        'clade': 'Amphibia',
        'divergence_mya': 270,
        'state': 'LINKED_ANCESTRAL',
        'color': '#388e3c',
        'alpha_chr': 'Chr 2',
        'alpha_pos': 0.12,
        'beta_chr': 'Chr 2',
        'beta_pos': 0.12,
        'alpha_anchors': 'NPRL3-RHBDF1',
        'beta_anchors': '',
        'notes': 'Giant salamander genome (32 Gb); retains primordial linked alpha-beta locus on Chr 2',
    },
    # 9-12. Archosauromorpha (Testudines, Crocodylia, Aves)
    {
        'species': 'chelonia_mydas',
        'display': 'Green Sea Turtle (C. mydas)',
        'clade': 'Testudines',
        'divergence_mya': 245,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'alpha_chr': 'Chr 14',
        'alpha_pos': 0.04,
        'beta_chr': 'Chr 21',
        'beta_pos': 0.03,
        'alpha_anchors': 'NPRL3-LUC7L',
        'beta_anchors': 'OR51-RRM1',
        'notes': 'Archosauromorph sister; balanced microchromosomal placement',
    },
    {
        'species': 'alligator_mississippiensis',
        'display': 'American Alligator (A. mississippiensis)',
        'clade': 'Crocodylia',
        'divergence_mya': 240,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'alpha_chr': 'Chr 15',
        'alpha_pos': 0.02,
        'beta_chr': 'Chr 18',
        'beta_pos': 0.02,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'OR51-NUMA1',
        'notes': 'Both clusters situated near microchromosome ends',
    },
    {
        'species': 'Gallus_gallus',
        'display': 'Chicken (G. gallus)',
        'clade': 'Aves',
        'divergence_mya': 100,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'alpha_chr': 'Chr 14',
        'alpha_pos': 0.02,
        'beta_chr': 'Chr 25',
        'beta_pos': 0.01,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'OR51-NUMA1',
        'notes': 'Classic avian microchromosomes (Chr 14 & 25); parallel high GC3 (>75%)',
    },
    {
        'species': 'calypte_anna',
        'display': "Anna's Hummingbird (C. anna)",
        'clade': 'Aves',
        'divergence_mya': 100,
        'state': 'ARCHOSAUR_BALANCED',
        'color': '#00796b',
        'alpha_chr': 'Chr 14',
        'alpha_pos': 0.02,
        'beta_chr': 'Chr 25',
        'beta_pos': 0.01,
        'alpha_anchors': 'NPRL3-LUC7L',
        'beta_anchors': 'OR51-NUMA1',
        'notes': 'High metabolic rate; both microchromosomal clusters maintain elevated GC',
    },
    # 13-15. Squamata - The Inversion
    {
        'species': 'podarcis_muralis',
        'display': 'Wall Lizard (P. muralis)',
        'clade': 'Lacertoidea',
        'divergence_mya': 165,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#d32f2f',
        'alpha_chr': 'Chr 5',
        'alpha_pos': 0.023,
        'beta_chr': 'Chr 18',
        'beta_pos': 0.014,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': 'CDCA8-OR51-RRM1',
        'notes': 'Squamate microchromosomal β; both loci capture subtelomeric domains',
    },
    {
        'species': 'pantherophis_guttatus',
        'display': 'Corn Snake (P. guttatus)',
        'clade': 'Serpentes',
        'divergence_mya': 110,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#c2185b',
        'alpha_chr': 'Chr 4',
        'alpha_pos': 0.75,
        'beta_chr': 'Chr 16',
        'beta_pos': 0.01,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': 'CDCA8-OR51-RRM1-STIM1',
        'notes': 'Runaway β-GC3 (92.2%) driven by subtelomeric microchromosome capture',
    },
    {
        'species': 'vipera_berus',
        'display': 'Common Adder (V. berus)',
        'clade': 'Serpentes',
        'divergence_mya': 110,
        'state': 'SQUAMATE_BETA_INVERTED',
        'color': '#c2185b',
        'alpha_chr': 'Chr 4',
        'alpha_pos': 0.74,
        'beta_chr': 'Chr 16',
        'beta_pos': 0.01,
        'alpha_anchors': 'NPRL3-MPG',
        'beta_anchors': 'CDCA8-OR51-RRM1',
        'notes': 'Viperid extreme polarity reversal: β-GC3 93.3%, α-GC3 70.6%',
    },
    # 16-17. Mammalia - Monotremata & Marsupialia
    {
        'species': 'ornithorhynchus_anatinus',
        'display': 'Platypus (O. anatinus)',
        'clade': 'Monotremata',
        'divergence_mya': 165,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#f57c00',
        'alpha_chr': 'Chr 21',
        'alpha_pos': 0.03,
        'beta_chr': 'Chr 18',
        'beta_pos': 0.65,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'OR51-NUMA1',
        'notes': 'Basal mammal; α already captures subtelomeric position on small chr',
    },
    {
        'species': 'Monodelphis_domestica',
        'display': 'Gray Opossum (M. domestica)',
        'clade': 'Marsupialia',
        'divergence_mya': 140,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#f57c00',
        'alpha_chr': 'Chr 1',
        'alpha_pos': 0.04,
        'beta_chr': 'Chr 3',
        'beta_pos': 0.68,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Marsupial canonical arrangement: α subtelomeric, β interstitial',
    },
    # 18-21. Mammalia - Eutheria (Placentalia)
    {
        'species': 'Homo_sapiens',
        'display': 'Human (H. sapiens)',
        'clade': 'Primates',
        'divergence_mya': 85,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'alpha_chr': 'Chr 16p',
        'alpha_pos': 0.01,
        'beta_chr': 'Chr 11p',
        'beta_pos': 0.65,
        'alpha_anchors': 'NPRL3-MPG-RHBDF1-LUC7L',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Canonical mammalian paradigm: α at extreme 16p telomere (GC3=91.1%), β at 11p (65.8%)',
    },
    {
        'species': 'Mus_musculus',
        'display': 'House Mouse (M. musculus)',
        'clade': 'Rodentia',
        'divergence_mya': 85,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'alpha_chr': 'Chr 11',
        'alpha_pos': 0.02,
        'beta_chr': 'Chr 7',
        'beta_pos': 0.70,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Murine chromosome reshuffling; α remains subtelomeric, β deep internal',
    },
    {
        'species': 'Bos_taurus',
        'display': 'Cattle (B. taurus)',
        'clade': 'Artiodactyla',
        'divergence_mya': 60,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'alpha_chr': 'Chr 25',
        'alpha_pos': 0.03,
        'beta_chr': 'Chr 15',
        'beta_pos': 0.68,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Ruminant tandem β expansion; α-GC3 89.5%, β-GC3 65.5%',
    },
    {
        'species': 'physeter_catodon',
        'display': 'Sperm Whale (P. catodon)',
        'clade': 'Cetacea',
        'divergence_mya': 60,
        'state': 'MAMMALIAN_ALPHA_SUBTELO',
        'color': '#e65100',
        'alpha_chr': 'Scaffold NW_021146805.1',
        'alpha_pos': 0.01,
        'beta_chr': 'Chr 14',
        'beta_pos': 0.72,
        'alpha_anchors': 'NPRL3-MPG-LUC7L',
        'beta_anchors': 'FOLR1-OR51-NUMA1',
        'notes': 'Abyssal diver; intact 5-gene α-cluster on GC-rich scaffold, β deep interstitial',
    },
]

def load_data():
    df_loc = pd.read_csv(os.path.join(RESULTS_DIR, "locus_codon_and_synteny_summary.csv"))
    df_sp = pd.read_csv(os.path.join(RESULTS_DIR, "species_locus_evolutionary_history.csv"))
    df_nodes = pd.read_csv(os.path.join(RESULTS_DIR, "ancestral_node_reconstruction.csv"))
    df_codons = pd.read_csv(os.path.join(RESULTS_DIR, "codon_decomposition_per_gene.csv"))
    return df_loc, df_sp, df_nodes, df_codons

def plot_figure4_publication(df_loc, df_sp, df_nodes, df_codons, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, "summary_figure4_phylogenetic_synteny.png")
    pdf_path = os.path.join(out_dir, "summary_figure4_phylogenetic_synteny.pdf")
    
    n_taxa = len(REPRESENTATIVE_TAXA)
    fig_h = max(20, n_taxa * 0.95)
    
    # 3 major panels:
    # 1. Phylogenetic Tree & Time Scale (width 4.8)
    # 2. Syntenic Chromosome Schematics (width 4.3)
    # 3. Codon Decomposition & Synonymous Divergence (width 4.3)
    fig, (ax_tree, ax_syn, ax_codon) = plt.subplots(
        1, 3, figsize=(31, fig_h),
        gridspec_kw={'width_ratios': [4.8, 4.3, 4.3]},
        dpi=300
    )
    plt.subplots_adjust(wspace=0.10, left=0.03, right=0.97, top=0.94, bottom=0.04)
    
    y_positions = np.arange(n_taxa)[::-1]  # Top to bottom (index 0 at top, y=21)
    
    # -------------------------------------------------------------
    # PANEL 1: Strictly Orthogonal Chronogram Tree & Timeline
    # -------------------------------------------------------------
    ax_tree.set_title("A. Calibrated Vertebrate Tree & Evolutionary Milestones", fontsize=15, fontweight='bold', pad=18)
    ax_tree.set_xlim(670, -215)
    ax_tree.set_ylim(-0.8, n_taxa - 0.2)
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
    ax_tree.plot([520, 500], [y_root, y_root], color='#424242', linewidth=2.6, zorder=2) # Stem
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

    # 8. Amniota (318 Ma)
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
        # Tip marker
        ax_tree.scatter(0, y, color=col, s=48, edgecolors='black', linewidths=0.9, zorder=5)
        # Taxon name text placed cleanly inside ax_tree (never clipped)
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

    # Clade Label Banners along the left margin (cleanly inside xlim, no overlap)
    clade_spans = [
        ("Agnatha", y_positions[0], y_positions[0], "#7b1fa2"),
        ("Chondrichthyes", y_positions[1], y_positions[2], "#00838f"),
        ("Actinopterygii", y_positions[3], y_positions[5], "#1565c0"),
        ("Sarcopterygii", y_positions[6], y_positions[6], "#2e7d32"),
        ("Amphibia", y_positions[7], y_positions[8], "#388e3c"),
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
    ax_syn.set_xlim(-0.02, 1.02)
    ax_syn.set_ylim(-0.8, n_taxa - 0.2)
    ax_syn.set_xlabel("Relative Chromosome Position (0.0 = 5' Telomere, 1.0 = 3' Telomere)", fontsize=13, fontweight='bold', labelpad=10)
    
    # 5' and 3' Telomere markers
    ax_syn.axvline(0.0, color="#bdbdbd", linestyle=":", linewidth=1.0)
    ax_syn.axvline(1.0, color="#bdbdbd", linestyle=":", linewidth=1.0)
    ax_syn.text(0.0, -0.6, "5' Telomere", fontsize=8.5, ha='center', va='top', color='#757575', fontweight='bold')
    ax_syn.text(1.0, -0.6, "3' Telomere", fontsize=8.5, ha='center', va='top', color='#757575', fontweight='bold')

    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        state = t['state']
        
        # Centromere position
        cen_x = 0.48 if state != 'MAMMALIAN_ALPHA_SUBTELO' else 0.42
        
        # Left arm ribbon
        left_arm = FancyBboxPatch((0.01, y - 0.14), cen_x - 0.03, 0.28,
                                  boxstyle="round,pad=0.01,rounding_size=0.04",
                                  fc="#f8f9fa", ec="#cccccc", lw=1.1, zorder=2)
        ax_syn.add_patch(left_arm)
        
        # Right arm ribbon
        right_arm = FancyBboxPatch((cen_x + 0.02, y - 0.14), 0.98 - cen_x - 0.02, 0.28,
                                   boxstyle="round,pad=0.01,rounding_size=0.04",
                                   fc="#f8f9fa", ec="#cccccc", lw=1.1, zorder=2)
        ax_syn.add_patch(right_arm)
        
        # Centromere constriction circle
        ax_syn.plot([cen_x - 0.015, cen_x + 0.015], [y, y], color="#9e9e9e", lw=3.0, zorder=3)
        ax_syn.scatter(cen_x, y, s=35, marker='o', fc='#757575', ec='black', lw=0.6, zorder=4)

        # Draw Globin Loci & Syntenic Anchors
        if state == 'CYCLOSTOME_PROTO':
            ax_syn.add_patch(Rectangle((0.04, y - 0.12), 0.07, 0.24, fc="#ab47bc", ec="black", lw=1.0, zorder=5))
            ax_syn.text(0.075, y, "Proto-Hb", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(0.13, y, "Ancestral Agnathan Single Locus", fontsize=8, va='center', ha='left', color="#424242", zorder=6)
            
        elif state == 'LINKED_ANCESTRAL':
            pos = t['alpha_pos']
            ax_syn.add_patch(Rectangle((pos, y - 0.12), 0.055, 0.24, fc="#1e88e5", ec="black", lw=1.0, zorder=5))
            ax_syn.add_patch(Rectangle((pos + 0.055, y - 0.12), 0.055, 0.24, fc="#43a047", ec="black", lw=1.0, zorder=5))
            ax_syn.text(pos + 0.027, y, "α", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(pos + 0.082, y, "β", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(pos + 0.12, y, f"◀ {t['alpha_anchors']} (linked)", fontsize=7.5, va='center', ha='left', color="#1565c0", zorder=6)
            
        elif state == 'CHONDRICHTHYAN_GRADIENT':
            a_pos = t['alpha_pos']
            b_pos = t['beta_pos']
            ax_syn.add_patch(Rectangle((a_pos, y - 0.12), 0.05, 0.24, fc="#1e88e5", ec="black", lw=1.0, zorder=5))
            ax_syn.text(a_pos + 0.025, y, "α", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.add_patch(Rectangle((b_pos, y - 0.12), 0.05, 0.24, fc="#43a047", ec="black", lw=1.0, zorder=5))
            ax_syn.text(b_pos + 0.025, y, "β", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(a_pos + 0.06, y, t['alpha_anchors'], fontsize=7.5, va='center', ha='left', color="#00838f")
            
        elif state == 'MAMMALIAN_ALPHA_SUBTELO':
            a_pos = 0.03
            b_pos = 0.68
            ax_syn.add_patch(Rectangle((a_pos, y - 0.12), 0.055, 0.24, fc="#1e88e5", ec="#0d47a1", lw=1.2, zorder=5))
            ax_syn.text(a_pos + 0.027, y, "α", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(a_pos + 0.065, y, "NPRL3-MPG-LUC7L", fontsize=7.5, va='center', ha='left', color="#0d47a1", zorder=6)
            
            ax_syn.add_patch(Rectangle((b_pos, y - 0.12), 0.055, 0.24, fc="#43a047", ec="#1b5e20", lw=1.2, zorder=5))
            ax_syn.text(b_pos + 0.027, y, "β", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(b_pos + 0.065, y, "FOLR1-OR51-NUMA1", fontsize=7.5, va='center', ha='left', color="#1b5e20", zorder=6)
            
        elif state == 'SQUAMATE_BETA_INVERTED':
            b_pos = 0.03
            a_pos = 0.72
            # Crimson Beta
            ax_syn.add_patch(Rectangle((b_pos, y - 0.12), 0.055, 0.24, fc="#d32f2f", ec="#b71c1c", lw=1.3, zorder=5))
            ax_syn.text(b_pos + 0.027, y, "β", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(b_pos + 0.065, y, "CDCA8-OR51-RRM1 ◀ INVERSION", fontsize=7.5, fontweight='bold', va='center', ha='left', color="#b71c1c", zorder=6)
            
            # Alpha internal
            ax_syn.add_patch(Rectangle((a_pos, y - 0.12), 0.055, 0.24, fc="#1e88e5", ec="#0d47a1", lw=1.0, zorder=5))
            ax_syn.text(a_pos + 0.027, y, "α", fontsize=9, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(a_pos + 0.065, y, "NPRL3", fontsize=7.5, va='center', ha='left', color="#0d47a1", zorder=6)

        elif state == 'ARCHOSAUR_BALANCED':
            a_pos = 0.03
            b_pos = 0.12
            ax_syn.add_patch(Rectangle((a_pos, y - 0.12), 0.05, 0.24, fc="#1e88e5", ec="black", lw=1.0, zorder=5))
            ax_syn.text(a_pos + 0.025, y, "α", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.add_patch(Rectangle((b_pos, y - 0.12), 0.05, 0.24, fc="#43a047", ec="black", lw=1.0, zorder=5))
            ax_syn.text(b_pos + 0.025, y, "β", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.text(b_pos + 0.065, y, "Microchromosomes (Chr 14 & 25)", fontsize=7.5, va='center', ha='left', color="#004d40", zorder=6)

        else: # General fallback
            a_pos = t.get('alpha_pos', 0.10)
            b_pos = t.get('beta_pos', 0.65)
            ax_syn.add_patch(Rectangle((a_pos, y - 0.12), 0.05, 0.24, fc="#1e88e5", ec="black", lw=1.0, zorder=5))
            ax_syn.text(a_pos + 0.025, y, "α", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)
            ax_syn.add_patch(Rectangle((b_pos, y - 0.12), 0.05, 0.24, fc="#43a047", ec="black", lw=1.0, zorder=5))
            ax_syn.text(b_pos + 0.025, y, "β", fontsize=8.5, fontweight='bold', ha='center', va='center', color='white', zorder=6)

    ax_syn.set_yticks([])
    ax_syn.spines['top'].set_visible(False)
    ax_syn.spines['right'].set_visible(False)
    ax_syn.spines['left'].set_visible(False)

    # -------------------------------------------------------------
    # PANEL 3: Codon Position Decomposition & Synonymous Divergence
    # -------------------------------------------------------------
    ax_codon.set_title("C. Codon Decomposition (GC1, GC2, GC3) & Divergence", fontsize=15, fontweight='bold', pad=18)
    ax_codon.set_xlim(30, 118)
    ax_codon.set_ylim(-0.8, n_taxa - 0.2)
    ax_codon.set_xlabel("GC Content (%) [Diamond = GC3; Shaded Bar = GC12 to GC3 Surge]", fontsize=13, fontweight='bold', labelpad=10)
    ax_codon.set_xticks([40, 50, 60, 70, 80, 90, 100])
    
    # Background guides
    for gc_mark in [40, 50, 60, 70, 80, 90, 100]:
        ax_codon.axvline(gc_mark, color="#f0f0f0", linestyle="--", linewidth=0.9, zorder=1)

    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        sp = t['species']
        sp_rows = df_sp[df_sp['species'] == sp]
        loc_rows = df_loc[df_loc['species'] == sp]
        
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
        
        # 1. Alpha Codons (Circle/Triangle/Diamond)
        if has_a:
            # Shaded bar from GC12 to GC3 highlighting synonymous surge
            a_gc12 = (a_gc1 + a_gc2) / 2.0
            ax_codon.plot([a_gc12, a_gc3], [y + 0.12, y + 0.12],
                          color="#ffcc80" if a_gc3 > 80 else "#bbdefb", linewidth=3.4, zorder=3, solid_capstyle='round')
            ax_codon.scatter(a_gc1, y + 0.12, marker='o', s=34, color='#1565c0', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(a_gc2, y + 0.12, marker='^', s=34, color='#42a5f5', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(a_gc3, y + 0.12, marker='D', s=65, color='#e65100' if a_gc3 > 80 else '#1e88e5', edgecolors='black', linewidths=0.8, zorder=5)

        # 2. Beta Codons (Square/Inverted-Triangle/Diamond)
        if has_b and t['state'] != 'CYCLOSTOME_PROTO':
            b_gc12 = (b_gc1 + b_gc2) / 2.0
            ax_codon.plot([b_gc12, b_gc3], [y - 0.12, y - 0.12],
                          color="#ef9a9a" if b_gc3 > 80 else "#c8e6c9", linewidth=3.4, zorder=3, solid_capstyle='round')
            ax_codon.scatter(b_gc1, y - 0.12, marker='s', s=34, color='#2e7d32', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(b_gc2, y - 0.12, marker='v', s=34, color='#81c784', edgecolors='black', linewidths=0.6, zorder=4)
            ax_codon.scatter(b_gc3, y - 0.12, marker='D', s=65, color='#d32f2f' if b_gc3 > 80 else '#43a047', edgecolors='black', linewidths=0.8, zorder=5)

        # Delta GC3 Badge on the right (positioned at 103, plenty of room before 118)
        if not sp_rows.empty and pd.notna(sp_rows['delta_gc3_a_minus_b'].iloc[0]):
            d_gc3 = sp_rows['delta_gc3_a_minus_b'].iloc[0]
            if d_gc3 > 5.0:
                txt = f"ΔGC3: +{d_gc3:.1f}% (α > β)"
                col = "#e65100"
                bg = "#fff3e0"
            elif d_gc3 < -5.0:
                txt = f"ΔGC3: {d_gc3:.1f}% (β > α)"
                col = "#b71c1c"
                bg = "#ffebee"
            else:
                txt = f"ΔGC3: {d_gc3:+.1f}% (Bal)"
                col = "#424242"
                bg = "#f5f5f5"
            ax_codon.text(103, y, txt, fontsize=8.5, fontweight='bold', va='center', ha='left', color=col,
                          bbox=dict(boxstyle="round,pad=0.25", fc=bg, ec=col, lw=0.9))

    # Legend for Panel 3
    leg_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1565c0', markeredgecolor='black', markersize=6.5, label='α-Globin GC1'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='#42a5f5', markeredgecolor='black', markersize=6.5, label='α-Globin GC2'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#e65100', markeredgecolor='black', markersize=7.5, label='α-Globin GC3 (Synonymous Surge)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#2e7d32', markeredgecolor='black', markersize=6.5, label='β-Globin GC1'),
        Line2D([0], [0], marker='v', color='w', markerfacecolor='#81c784', markeredgecolor='black', markersize=6.5, label='β-Globin GC2'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#d32f2f', markeredgecolor='black', markersize=7.5, label='β-Globin GC3 (Synonymous Surge)'),
    ]
    ax_codon.legend(handles=leg_elements, loc='lower left', fontsize=8.5, frameon=True, facecolor='#fafafa', edgecolor='#cccccc', ncol=2)

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
        column_widths=[0.32, 0.36, 0.32],
        subplot_titles=[
            "<b>A. Calibrated Vertebrate Chronogram & Transitions</b>",
            "<b>B. Conserved Syntenic Architecture (0 = Telomere)</b>",
            "<b>C. Codon Decomposition (GC1, GC2, GC3)</b>"
        ],
        horizontal_spacing=0.06
    )
    
    # 1. Panel A: Tree Trunks and Species Tips
    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        sp = t['species']
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
        state = t['state']
        cen_x = 0.48 if state != 'MAMMALIAN_ALPHA_SUBTELO' else 0.42
        
        # Chromosome ribbon
        fig.add_trace(go.Scatter(
            x=[0.01, 0.99], y=[y, y],
            mode='lines',
            line=dict(color='#e0e0e0', width=16),
            hoverinfo='skip',
            showlegend=False
        ), row=1, col=2)
        
        # Centromere
        fig.add_trace(go.Scatter(
            x=[cen_x], y=[y],
            mode='markers',
            marker=dict(size=10, color='#757575', symbol='circle', line=dict(color='black', width=1)),
            hovertext=f"Centromere (cen ~ {cen_x:.2f})",
            hoverinfo='text',
            showlegend=False
        ), row=1, col=2)
        
        # Globin Loci markers
        if state == 'CYCLOSTOME_PROTO':
            fig.add_trace(go.Scatter(
                x=[0.08], y=[y],
                mode='markers+text',
                marker=dict(size=15, color='#ab47bc', symbol='square'),
                text=["Proto-Hb"],
                textposition="middle right",
                hovertext="Proto-Globin Single Locus (GL478636)",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
        elif state == 'LINKED_ANCESTRAL':
            pos = t['alpha_pos']
            fig.add_trace(go.Scatter(
                x=[pos], y=[y],
                mode='markers',
                marker=dict(size=14, color='#1e88e5', symbol='square', line=dict(color='black', width=1)),
                hovertext=f"Linked α-Globin ({t['alpha_chr']}: {t['alpha_anchors']})",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
            fig.add_trace(go.Scatter(
                x=[pos + 0.035], y=[y],
                mode='markers+text',
                marker=dict(size=14, color='#43a047', symbol='square', line=dict(color='black', width=1)),
                text=[f"  α-β linked ({t['alpha_anchors']})"],
                textposition="middle right",
                hovertext=f"Linked β-Globin ({t['alpha_chr']}: {t['alpha_anchors']})",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
        elif state == 'SQUAMATE_BETA_INVERTED':
            fig.add_trace(go.Scatter(
                x=[0.03], y=[y],
                mode='markers+text',
                marker=dict(size=15, color='#d32f2f', symbol='square'),
                text=["β (subtelo)"],
                textposition="middle right",
                hovertext="β-Globin Inverted to Subtelomere (CDCA8-OR51)",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
            fig.add_trace(go.Scatter(
                x=[0.72], y=[y],
                mode='markers+text',
                marker=dict(size=15, color='#1e88e5', symbol='square'),
                text=["α (internal)"],
                textposition="middle right",
                hovertext="α-Globin Retained Internally (NPRL3-MPG)",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
        else: # Mammalia, Amphibia, Archosauria, Chondrichthyes
            a_p = t['alpha_pos'] if t['alpha_pos'] is not None else 0.03
            b_p = t['beta_pos'] if t['beta_pos'] is not None else 0.68
            fig.add_trace(go.Scatter(
                x=[a_p], y=[y],
                mode='markers',
                marker=dict(size=14, color='#1e88e5', symbol='square'),
                hovertext=f"α-Globin: {t['alpha_chr']} ({t['alpha_anchors']})",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)
            fig.add_trace(go.Scatter(
                x=[b_p], y=[y],
                mode='markers',
                marker=dict(size=14, color='#43a047', symbol='square'),
                hovertext=f"β-Globin: {t['beta_chr']} ({t['beta_anchors']})",
                hoverinfo='text',
                showlegend=False
            ), row=1, col=2)

    # 3. Panel C: Codon Decomposition
    for i, t in enumerate(REPRESENTATIVE_TAXA):
        y = y_positions[i]
        sp = t['species']
        sp_rows = df_sp[df_sp['species'] == sp]
        loc_rows = df_loc[df_loc['species'] == sp]
        
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
        
        if has_a:
            a_gc12 = (a_gc1 + a_gc2) / 2.0
            
            # Connector bar
            fig.add_trace(go.Scatter(
                x=[a_gc12, a_gc3], y=[y + 0.12, y + 0.12],
                mode='lines',
                line=dict(color='#ffcc80' if a_gc3 > 80 else '#bbdefb', width=6),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=3)
            
            # GC3 marker
            hover_a = (
                f"<b>{t['display']} (α-Globin)</b><br>"
                f"GC1: {a_gc1:.1f}% | GC2: {a_gc2:.1f}%<br>"
                f"GC12: {a_gc12:.1f}%<br>"
                f"<b>GC3: {a_gc3:.1f}%</b> (Surge: {a_gc3 - a_gc12:+.1f}%)"
            )
            fig.add_trace(go.Scatter(
                x=[a_gc3], y=[y + 0.12],
                mode='markers',
                marker=dict(size=12, color='#e65100' if a_gc3 > 80 else '#1e88e5', symbol='diamond', line=dict(color='black', width=1)),
                hovertext=hover_a,
                hoverinfo='text',
                showlegend=False
            ), row=1, col=3)

        if has_b and t['state'] != 'CYCLOSTOME_PROTO':
            b_gc12 = (b_gc1 + b_gc2) / 2.0
            
            fig.add_trace(go.Scatter(
                x=[b_gc12, b_gc3], y=[y - 0.12, y - 0.12],
                mode='lines',
                line=dict(color='#ef9a9a' if b_gc3 > 80 else '#c8e6c9', width=6),
                hoverinfo='skip',
                showlegend=False
            ), row=1, col=3)
            
            hover_b = (
                f"<b>{t['display']} (β-Globin)</b><br>"
                f"GC1: {b_gc1:.1f}% | GC2: {b_gc2:.1f}%<br>"
                f"GC12: {b_gc12:.1f}%<br>"
                f"<b>GC3: {b_gc3:.1f}%</b> (Surge: {b_gc3 - b_gc12:+.1f}%)"
            )
            fig.add_trace(go.Scatter(
                x=[b_gc3], y=[y - 0.12],
                mode='markers',
                marker=dict(size=12, color='#d32f2f' if b_gc3 > 80 else '#43a047', symbol='diamond', line=dict(color='black', width=1)),
                hovertext=hover_b,
                hoverinfo='text',
                showlegend=False
            ), row=1, col=3)

    # Layout styling
    fig.update_xaxes(title_text="<b>Divergence Time (Ma)</b>", range=[560, -180], autorange="reversed", row=1, col=1)
    fig.update_xaxes(title_text="<b>Relative Position (0 = 5' Telomere, 1 = 3' Telomere)</b>", range=[-0.05, 1.05], row=1, col=2)
    fig.update_xaxes(title_text="<b>GC Content (%) [Diamond = GC3]</b>", range=[32, 105], row=1, col=3)
    
    fig.update_yaxes(showticklabels=False, range=[-0.8, n_taxa - 0.2])
    
    fig.update_layout(
        title="<b>Figure 4: Macroevolutionary Synteny Architecture and Codon Divergence Across Vertebrates</b>",
        height=max(900, n_taxa * 42),
        width=1750,
        template="plotly_white",
        margin=dict(l=40, r=40, t=80, b=50)
    )
    
    fig.write_html(html_path)
    print(f"Generated Figure 4 Interactive HTML: {html_path}")
    return html_path


def main():
    print("=== Generating Polished Figure 4: Phylogenetic Synteny & Codon Divergence Map ===")
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
            
    print("\nPolished Figure 4 generation and synchronization complete!")

if __name__ == "__main__":
    main()
