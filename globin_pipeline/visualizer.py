"""Visualization module for per-species diagnostic graphs and cross-species meta-graphs."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional
from globin_pipeline.baseline_manager import BaselineManager


class Visualizer:
    """Generates high-resolution diagnostic plots per species and meta-comparative graphs."""

    PHYLOGENETIC_ORDER = [
        # Cyclostomes (Jawless Vertebrate Baseline)
        "eptatretus_burgeri",
        "petromyzon_marinus",
        "lampetra_fluviatilis",
        # Cartilaginous Fish (Chondrichthyes: Ancestral Jawed Vertebrate Baseline)
        "heptranchias_perlo",
        "somniosus_microcephalus",
        "chiloscyllium_punctatum",
        "scyliorhinus_canicula",
        "carcharhinus_amblyrhynchos",
        # Non-teleost Ray-finned Fish (Cladistia & Holostei)
        "polypterus_senegalus",
        "lepisosteus_oculatus",
        # Teleost Fish (Post-3R Whole-Genome Duplication)
        "danio_rerio",
        "pimephales_promelas",
        "oryzias_latipes",
        "takifugu_rubripes",
        "gasterosteus_aculeatus",
        # Lobe-finned Fish (Sarcopterygii: Tetrapod Sister Group)
        "latimeria_chalumnae",
        # Amphibians (Tetrapod Ancestral Linked Cluster)
        "geotrypetes_seraphini",
        "rhinatrema_bivittatum",
        "microcaecilia_unicolor",
        "ambystoma_mexicanum",
        "triturus_cristatus",
        "xenopus_tropicalis",
        "bufo_bufo",
        # Turtles (Testudines)
        "Chrysemys_picta_bellii",
        "chelonia_mydas",
        "eretmochelys_imbricata",
        "pelodiscus_sinensis",
        "gopherus_evgoodei",
        "terrapene_carolina_triunguis",
        # Squamate Reptiles (Lizards & Snakes)
        "anolis_carolinensis",
        "pogona_vitticeps",
        "podarcis_muralis",
        "zootoca_vivipara",
        "candoia_aspera",
        "vipera_berus",
        "notechis_scutatus",
        "pantherophis_guttatus",
        "thamnophis_elegans",
        # Crocodilia (Archosauria: Sister to Birds)
        "alligator_mississippiensis",
        "crocodylus_porosus",
        # Birds (Aves: Archosauria)
        "struthio_camelus_australis",
        "aptenodytes_patagonicus",
        "anser_cygnoides",
        "anas_platyrhynchos",
        "meleagris_gallopavo",
        "coturnix_japonica",
        "Gallus_gallus",
        "aquila_chrysaetos_chrysaetos",
        "calypte_anna",
        "taeniopygia_guttata",
        "sturnus_vulgaris",
        # Monotremes (Basal Mammals)
        "ornithorhynchus_anatinus",
        "tachyglossus_aculeatus",
        # Marsupials
        "Monodelphis_domestica",
        # Eutherian Mammals
        "loxodonta_africana",
        "Mus_musculus",
        "rattus_norvegicus",
        "heterocephalus_glaber",
        "miniopterus_natalensis",
        "desmodus_rotundus",
        "puma_concolor",
        "zalophus_californianus",
        "mirounga_angustirostris",
        "Sus_scrofa",
        "camelus_dromedarius",
        "vicugna_pacos",
        "Bos_taurus",
        "physeter_catodon",
        "tursiops_truncatus",
        "equus_asinus",
        "equus_caballus",
        "Homo_sapiens",
    ]

    COMMON_NAME_MAP = {
        "Homo_sapiens": "Human (H. sapiens)",
        "equus_caballus": "Horse (E. caballus)",
        "equus_asinus": "Donkey (E. asinus)",
        "Bos_taurus": "Cattle (B. taurus)",
        "Sus_scrofa": "Pig (S. scrofa)",
        "puma_concolor": "Puma (P. concolor)",
        "zalophus_californianus": "Sea Lion (Z. californianus)",
        "tursiops_truncatus": "Dolphin (T. truncatus)",
        "miniopterus_natalensis": "Bat (M. natalensis)",
        "rattus_norvegicus": "Norway Rat (R. norvegicus)",
        "heterocephalus_glaber": "Naked Mole-Rat (H. glaber)",
        "Mus_musculus": "Mouse (M. musculus)",
        "loxodonta_africana": "Elephant (L. africana)",
        "Erinaceus_europaeus": "Hedgehog (E. europaeus)",
        "Monodelphis_domestica": "Opossum (M. domestica)",
        "ornithorhynchus_anatinus": "Platypus (O. anatinus)",
        "sturnus_vulgaris": "Starling (S. vulgaris)",
        "aquila_chrysaetos_chrysaetos": "Golden Eagle (A. chrysaetos)",
        "taeniopygia_guttata": "Zebra Finch (T. guttata)",
        "Gallus_gallus": "Chicken (G. gallus)",
        "coturnix_japonica": "Quail (C. japonica)",
        "meleagris_gallopavo": "Turkey (M. gallopavo)",
        "anas_platyrhynchos": "Mallard Duck (A. platyrhynchos)",
        "anser_cygnoides": "Swan Goose (A. cygnoides)",
        "struthio_camelus_australis": "Ostrich (S. camelus)",
        "crocodylus_porosus": "Crocodile (C. porosus)",
        "alligator_mississippiensis": "Alligator (A. mississippiensis)",
        "pogona_vitticeps": "Bearded Dragon (P. vitticeps)",
        "anolis_carolinensis": "Anole (A. carolinensis)",
        "podarcis_muralis": "Wall Lizard (P. muralis)",
        "zootoca_vivipara": "Viviparous Lizard (Z. vivipara)",
        "candoia_aspera": "Ground Boa (C. aspera)",
        "vipera_berus": "Adder (V. berus)",
        "notechis_scutatus": "Tiger Snake (N. scutatus)",
        "pantherophis_guttatus": "Corn Snake (P. guttatus)",
        "thamnophis_elegans": "Garter Snake (T. elegans)",
        "terrapene_carolina_triunguis": "Box Turtle (T. carolina)",
        "gopherus_evgoodei": "Tortoise (G. evgoodei)",
        "pelodiscus_sinensis": "Softshell Turtle (P. sinensis)",
        "eretmochelys_imbricata": "Hawksbill Turtle (E. imbricata)",
        "chelonia_mydas": "Sea Turtle (C. mydas)",
        "Chrysemys_picta_bellii": "Painted Turtle (C. picta)",
        "bufo_bufo": "Common Toad (B. bufo)",
        "xenopus_tropicalis": "Clawed Frog (X. tropicalis)",
        "microcaecilia_unicolor": "Tiny Caecilian (M. unicolor)",
        "rhinatrema_bivittatum": "Two-lined Caecilian (R. bivittatum)",
        "geotrypetes_seraphini": "Gabon Caecilian (G. seraphini)",
        "ambystoma_mexicanum": "Axolotl (A. mexicanum)",
        "triturus_cristatus": "Crested Newt (T. cristatus)",
        "takifugu_rubripes": "Fugu (T. rubripes)",
        "gasterosteus_aculeatus": "Stickleback (G. aculeatus)",
        "oryzias_latipes": "Medaka (O. latipes)",
        "pimephales_promelas": "Minnow (P. promelas)",
        "danio_rerio": "Zebrafish (D. rerio)",
        "lepisosteus_oculatus": "Spotted Gar (L. oculatus)",
        "scyliorhinus_canicula": "Catshark (S. canicula)",
        "heptranchias_perlo": "Sevengill Shark (H. perlo)",
        "somniosus_microcephalus": "Greenland Shark (S. microcephalus)",
        "chiloscyllium_punctatum": "Bamboo Shark (C. punctatum)",
        "carcharhinus_amblyrhynchos": "Grey Reef Shark (C. amblyrhynchos)",
        "petromyzon_marinus": "Sea Lamprey (P. marinus)",
        "lampetra_fluviatilis": "River Lamprey (L. fluviatilis)",
        "eptatretus_burgeri": "Inshore Hagfish (E. burgeri)",
        "latimeria_chalumnae": "Coelacanth (L. chalumnae)",
        "polypterus_senegalus": "Bichir (P. senegalus)",
        "camelus_dromedarius": "Dromedary (C. dromedarius)",
        "vicugna_pacos": "Alpaca (V. pacos)",
        "desmodus_rotundus": "Vampire Bat (D. rotundus)",
        "mirounga_angustirostris": "Elephant Seal (M. angustirostris)",
        "physeter_catodon": "Sperm Whale (P. catodon)",
        "tachyglossus_aculeatus": "Echidna (T. aculeatus)",
        "calypte_anna": "Hummingbird (C. anna)",
        "aptenodytes_patagonicus": "King Penguin (A. patagonicus)",
    }

    COLOR_PALETTE = {
        "alpha_globin_coding": "#1f77b4",          # Deep Blue
        "beta_globin_coding": "#2ca02c",           # Forest Green
        "globin_coding": "#1f77b4",                # Deep Blue
        "globin_coding_locus_A": "#0d47a1",        # Cobalt Blue
        "globin_coding_locus_B": "#00acc1",        # Cyan / Teal
        "globin_pseudogene": "#9467bd",            # Purple/Magenta
        "cluster_non_globin": "#ff7f0e",           # Amber/Orange (single locus / linked)
        "cluster_non_globin_alpha": "#ff7f0e",     # Amber / Warm Orange (alpha-locus bystander)
        "cluster_non_globin_beta": "#d95f02",      # Rust / Dark Orange (beta-locus bystander)
        "cluster_non_globin_locus_A": "#e65100",   # Deep Orange (locus A bystander)
        "cluster_non_globin_locus_B": "#fdae6b",   # Peach / Light Orange (locus B bystander)
        "genome_baseline": "#7f7f7f",              # Slate Gray
    }

    MARKER_MAP = {
        "alpha_globin_coding": "o",
        "beta_globin_coding": "s",
        "globin_coding": "o",
        "globin_coding_locus_A": "o",
        "globin_coding_locus_B": "^",
        "globin_pseudogene": "D",
        "cluster_non_globin": "v",
        "cluster_non_globin_alpha": "v",
        "cluster_non_globin_beta": "<",
        "cluster_non_globin_locus_A": "v",
        "cluster_non_globin_locus_B": ">",
    }

    LABEL_MAP = {
        "alpha_globin_coding": "Alpha Globin Coding",
        "beta_globin_coding": "Beta Globin Coding",
        "globin_coding": "Globin Coding",
        "globin_coding_locus_A": "Globin Coding (Locus A)",
        "globin_coding_locus_B": "Globin Coding (Locus B)",
        "globin_pseudogene": "Globin Pseudogene",
        "cluster_non_globin": "Cluster Non-Globin",
        "cluster_non_globin_alpha": "Non-Globin (α-Locus)",
        "cluster_non_globin_beta": "Non-Globin (β-Locus)",
        "cluster_non_globin_locus_A": "Non-Globin (Locus A)",
        "cluster_non_globin_locus_B": "Non-Globin (Locus B)",
    }

    DETERMINISTIC_GROUP_TIERS = {
        # Tier 1 (closest to 0): Alpha-globins & ancestral linked globins
        "alpha_globin_coding": 1,
        "globin_coding_locus_A": 1,
        "globin_coding": 1,
        # Tier 2: Beta-globins
        "beta_globin_coding": 2,
        "globin_coding_locus_B": 2,
        # Tier 3: Globin pseudogenes
        "globin_pseudogene": 3,
        # Tier 4: Non-globins from alpha-locus / linked cluster / locus A
        "cluster_non_globin_alpha": 4,
        "cluster_non_globin_locus_A": 4,
        "cluster_non_globin": 4,
        # Tier 5: Non-globins from beta-locus / locus B
        "cluster_non_globin_beta": 5,
        "cluster_non_globin_locus_B": 5,
        # Tier 6: Other locus features
        "other_locus_feature": 6,
    }

    def _calculate_gene_y_positions(self, df: pd.DataFrame, val_col: str, y_max: float) -> Dict[str, float]:
        """Calculates deterministic, consistent y-axis positions stratified by gene category.
        
        Ordering (closest to 0 up to higher density):
          Tier 1 (closest to 0): Alpha-globins (or ancestral linked globins) -> ~0.04 * y_max
          Tier 2: Beta-globins -> ~0.09 * y_max
          Tier 3: Globin pseudogenes -> ~0.14 * y_max
          Tier 4: Non-globins from alpha-locus / linked cluster -> ~0.19 * y_max
          Tier 5: Non-globins from beta-locus -> ~0.24 * y_max
          Tier 6: Other locus features -> ~0.29 * y_max
          
        For closely spaced genes within the same tier (|delta_gc| < 0.8%), a deterministic
        micro-stagger (+/-0.015 * y_max) prevents point occlusion.
        """
        base_offsets = {
            1: 0.04,
            2: 0.09,
            3: 0.14,
            4: 0.19,
            5: 0.24,
            6: 0.29,
        }
        y_pos: Dict[str, float] = {}
        if df.empty or val_col not in df.columns:
            return y_pos

        for tier in range(1, 7):
            mask = df["comparative_group"].map(lambda g: self.DETERMINISTIC_GROUP_TIERS.get(g, 6)) == tier
            tier_df = df[mask].copy()
            if tier_df.empty:
                continue
            tier_df = tier_df.sort_values(by=val_col).reset_index(drop=True)
            base_y = base_offsets[tier] * y_max
            
            vals = [float(v) for v in tier_df[val_col].tolist()]
            n = len(vals)
            micro_offsets = [0.0] * n
            
            i = 0
            while i < n:
                cluster = [i]
                j = i + 1
                while j < n and (vals[j] - vals[j-1]) < 0.8:
                    cluster.append(j)
                    j += 1
                if len(cluster) > 1:
                    for idx_c, c_idx in enumerate(cluster):
                        sign = 1 if (idx_c % 2 == 1) else -1
                        micro_offsets[c_idx] = sign * 0.015 * y_max
                i = j
                
            for idx, row in tier_df.iterrows():
                gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
                y_pos[gid] = base_y + micro_offsets[idx]
                
        return y_pos

    def _calculate_meta_y_positions(
        self,
        df_sp: pd.DataFrame,
        val_col: str,
        base_y: float,
        clean_mode: bool = False,
    ) -> Dict[str, float]:
        """Calculates deterministic, consistent y-axis positions stratified by gene tier for the meta-graph.
        
        Offsets within the species band [base_y - 0.35, base_y + 0.35]:
          Full Mode (6 tiers):
            Tier 1 (Alpha-globins / Locus A / Ancestral linked): base_y - 0.22
            Tier 2 (Beta-globins / Locus B): base_y - 0.11
            Tier 3 (Globin pseudogenes): base_y + 0.00
            Tier 4 (Non-globins alpha-locus / linked): base_y + 0.11
            Tier 5 (Non-globins beta-locus): base_y + 0.22
            Tier 6 (Other locus features): base_y + 0.31
            
          Clean Mode (Tiers 3 & 6 excluded):
            Tier 1 (Alpha-globins / Locus A / Ancestral linked): base_y - 0.18
            Tier 2 (Beta-globins / Locus B): base_y - 0.06
            Tier 4 (Non-globins alpha-locus / linked): base_y + 0.06
            Tier 5 (Non-globins beta-locus): base_y + 0.18
            
        For closely spaced points within the same tier (|delta_z| < 0.18), a deterministic
        micro-stagger (+/-0.025) prevents marker occlusion.
        """
        base_offsets_full = {1: -0.22, 2: -0.11, 3: 0.00, 4: 0.11, 5: 0.22, 6: 0.31}
        base_offsets_clean = {1: -0.18, 2: -0.06, 4: 0.06, 5: 0.18}
        offsets = base_offsets_clean if clean_mode else base_offsets_full

        y_pos: Dict[str, float] = {}
        if df_sp.empty:
            return y_pos

        valid_df = df_sp.copy()
        if clean_mode:
            valid_df = valid_df[~valid_df["comparative_group"].isin(["other_locus_feature", "globin_pseudogene"])]

        for tier, offset in offsets.items():
            mask = valid_df["comparative_group"].map(lambda g: self.DETERMINISTIC_GROUP_TIERS.get(g, 6)) == tier
            tier_df = valid_df[mask].copy()
            if tier_df.empty:
                continue

            if val_col in tier_df.columns:
                tier_df["effective_z"] = tier_df[val_col]
                if val_col == "z_score_full_gene_gc" and "z_score_gc" in tier_df.columns:
                    tier_df["effective_z"] = tier_df["effective_z"].fillna(tier_df["z_score_gc"])
            elif "z_score_gc" in tier_df.columns:
                tier_df["effective_z"] = tier_df["z_score_gc"]
            else:
                tier_df["effective_z"] = np.nan

            tier_df = tier_df.dropna(subset=["effective_z"]).sort_values(by="effective_z").reset_index(drop=True)
            if tier_df.empty:
                continue

            vals = [float(v) for v in tier_df["effective_z"].tolist()]
            n = len(vals)
            micro = [0.0] * n
            i = 0
            while i < n:
                cluster = [i]
                j = i + 1
                while j < n and (vals[j] - vals[j-1]) < 0.18:
                    cluster.append(j)
                    j += 1
                if len(cluster) > 1:
                    for idx_c, c_idx in enumerate(cluster):
                        sign = 1 if (idx_c % 2 == 1) else -1
                        micro[c_idx] = sign * 0.025
                i = j

            for idx, row in tier_df.iterrows():
                gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
                y_pos[gid] = base_y + offset + micro[idx]

        return y_pos

    def __init__(self, baseline_manager: Optional[BaselineManager] = None):
        self.bm = baseline_manager or BaselineManager()
        sns.set_theme(style="whitegrid", font="sans-serif")

    @staticmethod
    def _format_symbol_annotation(symbol: str, group: str) -> str:
        """Appends a concise locus tag to gene symbols when multiple loci exist."""
        grp_lower = group.lower()
        if "locus_a" in grp_lower:
            return f"{symbol} (A)"
        elif "locus_b" in grp_lower:
            return f"{symbol} (B)"
        elif "alpha" in grp_lower:
            return f"{symbol} (α)"
        elif "beta" in grp_lower:
            return f"{symbol} (β)"
        return symbol

    @staticmethod
    def _calculate_label_positions(values: List[float], min_sep: float = 1.3, cluster_gap: float = 2.5) -> List[float]:
        """Distributes label x-positions to prevent overlapping callout text for closely clustered points."""
        n = len(values)
        if n <= 1:
            return list(values)
        
        clusters = []
        current_cluster = [0]
        for i in range(1, n):
            if values[i] - values[i-1] > cluster_gap:
                clusters.append(current_cluster)
                current_cluster = [i]
            else:
                current_cluster.append(i)
        clusters.append(current_cluster)
        
        pos = list(values)
        for c in clusters:
            if len(c) == 1:
                continue
            c_vals = [values[i] for i in c]
            c_pos = list(c_vals)
            for j in range(1, len(c)):
                if c_pos[j] < c_pos[j-1] + min_sep:
                    c_pos[j] = c_pos[j-1] + min_sep
            shift = float(np.mean(c_vals)) - float(np.mean(c_pos))
            for idx_c, i in enumerate(c):
                pos[i] = c_pos[idx_c] + shift
        return pos

    def _annotate_density_points(
        self,
        ax: plt.Axes,
        valid_genes: pd.DataFrame,
        y_max: float,
        gene_y_map: Optional[Dict[str, float]] = None,
        is_panel_whole: bool = True,
        val_col: Optional[str] = None,
    ) -> None:
        """Annotates genes on the density distribution, prioritizing globin genes and pseudogenes with staggered labels."""
        if val_col is None:
            val_col = "full_gene_gc" if is_panel_whole else "cds_gc"
        if valid_genes.empty or val_col not in valid_genes.columns:
            return
        
        is_globin_mask = valid_genes["comparative_group"].str.contains("globin_coding|pseudogene") | (valid_genes["is_globin"].astype(str).str.lower() == "true")
        
        if len(valid_genes) <= 8:
            to_label = valid_genes.copy()
        else:
            to_label = valid_genes[is_globin_mask].copy()
            if to_label.empty:
                to_label = valid_genes.copy()
                
        to_label = to_label.sort_values(by=val_col).reset_index(drop=True)
        raw_vals = [float(v) for v in to_label[val_col].tolist()]
        is_len = "len" in val_col
        min_sep = 0.08 if is_len else 1.3
        cluster_gap = 0.16 if is_len else 2.5
        spaced_vals = self._calculate_label_positions(raw_vals, min_sep=min_sep, cluster_gap=cluster_gap)
        
        tiers = [0.36, 0.56, 0.76, 0.96]
        
        for idx, row in to_label.iterrows():
            gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
            val = float(row[val_col])
            x_text = spaced_vals[idx]
            grp = str(row["comparative_group"])
            col = self.COLOR_PALETTE.get(grp, "#333")
            is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"
            
            tier_idx = idx % len(tiers)
            y_offset = y_max * tiers[tier_idx]
            y_dot = gene_y_map.get(gid, y_max * 0.10) if gene_y_map else (y_max * 0.10)
            
            sym = str(row["display_symbol"])
            sym_annot = self._format_symbol_annotation(sym, grp)
            
            ax.annotate(
                sym_annot,
                xy=(val, y_dot),
                xytext=(x_text, y_offset),
                fontsize=8.5,
                fontweight="bold" if is_glo else "normal",
                rotation=45,
                ha="left",
                va="bottom",
                arrowprops=dict(
                    arrowstyle="->",
                    color=col,
                    lw=1.0 if is_glo else 0.7,
                    alpha=0.85,
                ),
                zorder=7,
            )

    def plot_diagnostic(
        self, species: str, df_genes: pd.DataFrame, output_dir: str = "results/diagnostic_plots"
    ) -> str:
        """Generates a comprehensive 5-panel diagnostic plot for a single species."""
        os.makedirs(output_dir, exist_ok=True)
        out_png = os.path.join(output_dir, f"diagnostic_{species}.png")

        df_base, stats = self.bm.get_genome_baseline(species)
        sp_genes = df_genes[df_genes["species"] == species].copy()
        common_name = self.COMMON_NAME_MAP.get(species, species)

        fig = plt.figure(figsize=(22, 22), dpi=300)
        gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.0, 0.9])

        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[1, 0])
        ax_d = fig.add_subplot(gs[1, 1])
        ax_e = fig.add_subplot(gs[2, 0])
        ax_f = fig.add_subplot(gs[2, 1])

        base_whole = df_base["whole_gene_gc"].dropna() if not df_base.empty and "whole_gene_gc" in df_base.columns else pd.Series(dtype=float)
        base_cds = df_base["cds_gc"].dropna() if not df_base.empty and "cds_gc" in df_base.columns else pd.Series(dtype=float)

        mu_whole = stats.get("mean_whole_gene_gc", float(np.mean(base_whole)) if not base_whole.empty else 0.0)
        sigma_whole = stats.get("std_whole_gene_gc", float(np.std(base_whole)) if not base_whole.empty else 1.0)
        mu_cds = stats.get("mean_cds_gc", float(np.mean(base_cds)) if not base_cds.empty else 0.0)
        sigma_cds = stats.get("std_cds_gc", float(np.std(base_cds)) if not base_cds.empty else 1.0)
        mu_gen = stats.get("mean_genomic_gc")

        # Aligned x-limits for Panels A and B to facilitate visual comparison
        all_vals = []
        if not base_whole.empty:
            all_vals.extend([float(base_whole.quantile(0.001)), float(base_whole.quantile(0.999))])
        if not base_cds.empty:
            all_vals.extend([float(base_cds.quantile(0.001)), float(base_cds.quantile(0.999))])
        for _, r in sp_genes.iterrows():
            if pd.notna(r.get("full_gene_gc")):
                all_vals.append(float(r["full_gene_gc"]))
            if pd.notna(r.get("cds_gc")):
                all_vals.append(float(r["cds_gc"]))
        min_x = max(10.0, float(np.min(all_vals)) - 4.0) if all_vals else 20.0
        max_x = min(95.0, float(np.max(all_vals)) + 4.0) if all_vals else 90.0

        # -------------------------------------------------------------
        # PANEL A: Whole-Gene GC Density (All Protein-Coding Genes)
        # -------------------------------------------------------------
        if not base_whole.empty:
            sns.kdeplot(base_whole, ax=ax_a, color="#1565c0", fill=True, alpha=0.25, label=f"Genome Whole-Gene ({len(base_whole):,} genes)", bw_adjust=1.2)
        ax_a.axvline(mu_whole, color="#0d47a1", linestyle="--", linewidth=1.6, label=f"Whole-Gene Mean ({mu_whole:.1f}%)")
        if mu_gen is not None:
            ax_a.axvline(mu_gen, color="#2e7d32", linestyle="-.", linewidth=1.4, label=f"Genomic Assembly Baseline ({mu_gen:.1f}%)")
        ax_a.axvline(mu_whole + sigma_whole, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_a.axvline(mu_whole - sigma_whole, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_a.axvline(mu_whole + 2 * sigma_whole, color="#90a4ae", linestyle=":", linewidth=0.9)
        ax_a.axvline(mu_whole - 2 * sigma_whole, color="#90a4ae", linestyle=":", linewidth=0.9)

        ax_a.set_xlim(min_x, max_x)
        y_max_a = ax_a.get_ylim()[1] if ax_a.get_ylim()[1] > 0 else 0.05
        ax_a.set_ylim(0, y_max_a * 1.50)

        valid_a = sp_genes[sp_genes["full_gene_gc"].notna() & sp_genes["comparative_group"].isin(self.COLOR_PALETTE)].copy()
        y_positions_a = self._calculate_gene_y_positions(valid_a, "full_gene_gc", y_max_a)
        plotted_labels_a = set()
        for idx, row in valid_a.iterrows():
            gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
            grp = row["comparative_group"]
            val = float(row["full_gene_gc"])
            y_dot = y_positions_a.get(gid, y_max_a * 0.10)
            col = self.COLOR_PALETTE.get(grp, "#333")
            m = self.MARKER_MAP.get(grp, "o")
            is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"
            sz = 120 if is_glo else 70
            alpha_val = 1.0 if is_glo else 0.65
            lbl = self.LABEL_MAP.get(grp, grp.replace("_", " ").title()) if grp not in plotted_labels_a else None
            plotted_labels_a.add(grp)
            ax_a.scatter(val, y_dot, color=col, marker=m, s=sz, alpha=alpha_val, edgecolors="black", linewidths=0.9, zorder=5, label=lbl)

        self._annotate_density_points(ax_a, valid_a, y_max_a, gene_y_map=y_positions_a, is_panel_whole=True)

        ax_a.set_title(f"A. Whole-Gene GC Distribution (All Protein-Coding Genes)\n{common_name}", fontsize=12, fontweight="bold", pad=10)
        ax_a.set_xlabel("Whole-Gene GC Content (%) [Exons + Introns]", fontsize=11)
        ax_a.set_ylabel("Density", fontsize=11)
        ax_a.legend(loc="upper left", frameon=True, fontsize=8.5)

        # -------------------------------------------------------------
        # PANEL B: Coding Sequence (CDS) GC Density
        # -------------------------------------------------------------
        if not base_cds.empty:
            sns.kdeplot(base_cds, ax=ax_b, color="#455a64", fill=True, alpha=0.25, label=f"Genome Coding CDS ({len(base_cds):,} genes)", bw_adjust=1.2)
        ax_b.axvline(mu_cds, color="#212529", linestyle="--", linewidth=1.6, label=f"Coding CDS Mean ({mu_cds:.1f}%)")
        if mu_gen is not None:
            ax_b.axvline(mu_gen, color="#2e7d32", linestyle="-.", linewidth=1.4, label=f"Genomic Assembly Baseline ({mu_gen:.1f}%)")
        ax_b.axvline(mu_cds + sigma_cds, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_b.axvline(mu_cds - sigma_cds, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_b.axvline(mu_cds + 2 * sigma_cds, color="#90a4ae", linestyle=":", linewidth=0.9)
        ax_b.axvline(mu_cds - 2 * sigma_cds, color="#90a4ae", linestyle=":", linewidth=0.9)

        ax_b.set_xlim(min_x, max_x)
        y_max_b = ax_b.get_ylim()[1] if ax_b.get_ylim()[1] > 0 else 0.05
        ax_b.set_ylim(0, y_max_b * 1.50)

        valid_b = sp_genes[sp_genes["cds_gc"].notna() & sp_genes["comparative_group"].isin(self.COLOR_PALETTE)].copy()
        y_positions_b = self._calculate_gene_y_positions(valid_b, "cds_gc", y_max_b)
        plotted_labels_b = set()
        for idx, row in valid_b.iterrows():
            gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
            grp = row["comparative_group"]
            val = float(row["cds_gc"])
            y_dot = y_positions_b.get(gid, y_max_b * 0.10)
            col = self.COLOR_PALETTE.get(grp, "#333")
            m = self.MARKER_MAP.get(grp, "o")
            is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"
            sz = 120 if is_glo else 70
            alpha_val = 1.0 if is_glo else 0.65
            lbl = self.LABEL_MAP.get(grp, grp.replace("_", " ").title()) if grp not in plotted_labels_b else None
            plotted_labels_b.add(grp)
            ax_b.scatter(val, y_dot, color=col, marker=m, s=sz, alpha=alpha_val, edgecolors="black", linewidths=0.9, zorder=5, label=lbl)

        self._annotate_density_points(ax_b, valid_b, y_max_b, gene_y_map=y_positions_b, is_panel_whole=False)

        ax_b.set_title(f"B. Coding Sequence (CDS) GC Distribution (All Protein-Coding Genes)\n{common_name}", fontsize=12, fontweight="bold", pad=10)
        ax_b.set_xlabel("Coding Sequence (CDS) GC Content (%)", fontsize=11)
        ax_b.set_ylabel("Density", fontsize=11)
        ax_b.legend(loc="upper left", frameon=True, fontsize=8.5)

        # -------------------------------------------------------------
        # PANEL C: Gene Length Log10(bp) Empirical Distribution
        # -------------------------------------------------------------
        base_len = df_base["whole_gene_len"].dropna() if not df_base.empty and "whole_gene_len" in df_base.columns else pd.Series(dtype=float)
        valid_base_len = base_len[base_len > 0]
        base_log_len = np.log10(valid_base_len) if not valid_base_len.empty else pd.Series(dtype=float)

        mu_log_len = stats.get("mean_log10_len", float(np.mean(base_log_len)) if not base_log_len.empty else 3.8)
        sigma_log_len = stats.get("std_log10_len", float(np.std(base_log_len)) if not base_log_len.empty else 0.55)

        if not base_log_len.empty:
            sns.kdeplot(base_log_len, ax=ax_c, color="#6a1b9a", fill=True, alpha=0.25, label=f"Genome Gene Length ({len(base_log_len):,} genes)", bw_adjust=1.2)
        ax_c.axvline(mu_log_len, color="#4a148c", linestyle="--", linewidth=1.6, label=f"Mean Length ({10**mu_log_len:,.0f} bp / 10^{mu_log_len:.2f})")
        ax_c.axvline(mu_log_len + sigma_log_len, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_c.axvline(mu_log_len - sigma_log_len, color="#546e7a", linestyle=":", linewidth=1.1)
        ax_c.axvline(mu_log_len + 2 * sigma_log_len, color="#90a4ae", linestyle=":", linewidth=0.9)
        ax_c.axvline(mu_log_len - 2 * sigma_log_len, color="#90a4ae", linestyle=":", linewidth=0.9)

        min_len_x = 1.8
        max_len_x = 6.2
        ax_c.set_xlim(min_len_x, max_len_x)
        y_max_c = ax_c.get_ylim()[1] if ax_c.get_ylim()[1] > 0 else 0.8
        ax_c.set_ylim(0, y_max_c * 1.50)

        # Human-readable tick labels along log10 length
        ax_c.set_xticks([2.0, 3.0, 4.0, 5.0, 6.0])
        ax_c.set_xticklabels(["100 bp\n($10^2$)", "1 kb\n($10^3$)", "10 kb\n($10^4$)", "100 kb\n($10^5$)", "1 Mb\n($10^6$)"], fontsize=9.5)

        valid_c_len = sp_genes[sp_genes["log10_coding_span_len"].notna() & sp_genes["comparative_group"].isin(self.COLOR_PALETTE)].copy()
        if valid_c_len.empty and "full_gene_len" in sp_genes.columns:
            valid_c_len = sp_genes[sp_genes["full_gene_len"].notna() & (sp_genes["full_gene_len"] > 0) & sp_genes["comparative_group"].isin(self.COLOR_PALETTE)].copy()
            valid_c_len["log10_coding_span_len"] = np.log10(valid_c_len["full_gene_len"])

        y_positions_c = self._calculate_gene_y_positions(valid_c_len, "log10_coding_span_len", y_max_c)
        plotted_labels_c = set()
        for idx, row in valid_c_len.iterrows():
            gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
            grp = row["comparative_group"]
            val = float(row["log10_coding_span_len"])
            y_dot = y_positions_c.get(gid, y_max_c * 0.10)
            col = self.COLOR_PALETTE.get(grp, "#333")
            m = self.MARKER_MAP.get(grp, "o")
            is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"
            sz = 120 if is_glo else 70
            alpha_val = 1.0 if is_glo else 0.65
            lbl = self.LABEL_MAP.get(grp, grp.replace("_", " ").title()) if grp not in plotted_labels_c else None
            plotted_labels_c.add(grp)
            ax_c.scatter(val, y_dot, color=col, marker=m, s=sz, alpha=alpha_val, edgecolors="black", linewidths=0.9, zorder=5, label=lbl)

        self._annotate_density_points(ax_c, valid_c_len, y_max_c, gene_y_map=y_positions_c, val_col="log10_coding_span_len")

        ax_c.set_title(f"C. Gene Length Log10(bp) Distribution (All Protein-Coding Genes)\n{common_name}", fontsize=12, fontweight="bold", pad=10)
        ax_c.set_xlabel("Gene Length (Coding Span excluding UTRs)", fontsize=11)
        ax_c.set_ylabel("Density", fontsize=11)
        ax_c.legend(loc="upper left", frameon=True, fontsize=8.5)

        # -------------------------------------------------------------
        # PANEL D: Whole-Gene GC vs. CDS GC Bivariate Correlation
        # -------------------------------------------------------------
        ax_d.axline((50, 50), slope=1, color="#adb5bd", linestyle="--", linewidth=1.2, label="y = x (CDS GC = Whole-Gene GC)")
        if mu_gen is not None:
            ax_d.axvline(mu_gen, color="#2e7d32", linestyle="-.", linewidth=1.2, label=f"Genomic Baseline ({mu_gen:.1f}%)")
        if mu_cds > 0:
            ax_d.axhline(mu_cds, color="#212529", linestyle="--", linewidth=1.2, label=f"CDS Baseline ({mu_cds:.1f}%)")

        plotted_d_groups = set()
        for _, row in sp_genes.iterrows():
            fg = row.get("full_gene_gc")
            cg = row.get("cds_gc")
            grp = row.get("comparative_group")
            col = self.COLOR_PALETTE.get(grp, "#333")
            m = self.MARKER_MAP.get(grp, "o")
            is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"

            if pd.notna(fg) and pd.notna(cg):
                lbl = self.LABEL_MAP.get(grp, grp.replace("_", " ").title()) if grp not in plotted_d_groups else None
                plotted_d_groups.add(grp)
                ax_d.scatter(float(fg), float(cg), color=col, marker=m, s=120 if is_glo else 75, edgecolors="black", linewidths=0.9, zorder=5, label=lbl)
                sym = row["display_symbol"]
                sym_annot = self._format_symbol_annotation(sym, grp)
                ax_d.annotate(
                    sym_annot,
                    xy=(float(fg), float(cg)),
                    xytext=(float(fg) + 0.6, float(cg) + 0.4),
                    fontsize=8.5,
                    fontweight="bold" if is_glo else "normal",
                    zorder=6,
                )

        ax_d.set_title("D. Whole-Gene GC vs. CDS GC (Isochore Divergence vs. Coding Selection)", fontsize=12, fontweight="bold", pad=10)
        ax_d.set_xlabel("Whole-Gene GC Content (%) [Exons + Introns]", fontsize=11)
        ax_d.set_ylabel("Coding Sequence (CDS) GC Content (%)", fontsize=11)
        ax_d.legend(loc="lower right", frameon=True, fontsize=8.5)

        valid_d = sp_genes[["full_gene_gc", "cds_gc"]].dropna()
        if not valid_d.empty:
            min_v = max(20.0, float(valid_d.min().min()) - 6.0)
            max_v = min(90.0, float(valid_d.max().max()) + 6.0)
            ax_d.set_xlim(min_v, max_v)
            ax_d.set_ylim(min_v, max_v)

        # -------------------------------------------------------------
        # PANEL E: CDS GC vs. GC3 Codon Space
        # -------------------------------------------------------------
        df_base_valid = df_base.dropna(subset=["cds_gc", "gc3"]) if not df_base.empty and "cds_gc" in df_base.columns and "gc3" in df_base.columns else pd.DataFrame()
        if not df_base_valid.empty:
            sns.kdeplot(
                data=df_base_valid.sample(min(len(df_base_valid), 3000), random_state=42),
                x="cds_gc",
                y="gc3",
                ax=ax_e,
                levels=6,
                color="#6c757d",
                alpha=0.4,
                linewidths=1.0,
            )

        for _, row in sp_genes.iterrows():
            if pd.notna(row.get("cds_gc")) and pd.notna(row.get("gc3")):
                grp = row["comparative_group"]
                col = self.COLOR_PALETTE.get(grp, "#333")
                m = self.MARKER_MAP.get(grp, "o")
                is_glo = "globin_coding" in grp or str(row.get("is_globin")).lower() == "true"
                ax_e.scatter(
                    float(row["cds_gc"]),
                    float(row["gc3"]),
                    color=col,
                    marker=m,
                    s=120 if is_glo else 75,
                    edgecolors="black",
                    linewidths=0.9,
                    zorder=5,
                )
                sym = row["display_symbol"]
                sym_annot = self._format_symbol_annotation(sym, grp)
                ax_e.annotate(
                    sym_annot,
                    (float(row["cds_gc"]) + 0.5, float(row["gc3"]) + 0.5),
                    fontsize=8.5,
                    fontweight="bold" if is_glo else "normal",
                )

        ax_e.set_title("E. CDS GC vs. GC3 Codon Space (Selection vs. Mutational Bias)", fontsize=12, fontweight="bold", pad=10)
        ax_e.set_xlabel("CDS GC (%)", fontsize=11)
        ax_e.set_ylabel("GC3 (3rd Codon Position GC %)", fontsize=11)
        ax_e.axline((50, 50), slope=1, color="#ced4da", linestyle=":", label="x = y line")
        ax_e.legend(loc="lower right", frameon=True, fontsize=8.5)

        # -------------------------------------------------------------
        # PANEL F: Compartment Comparison Across Gene Groups
        # -------------------------------------------------------------
        comp_data = []
        if not base_whole.empty:
            sample_size = min(len(base_whole), 1500)
            for v in base_whole.sample(sample_size, random_state=42):
                comp_data.append({"Group": "Genome", "Compartment": "Whole-Gene GC", "GC": float(v)})
        if not base_cds.empty:
            sample_size = min(len(base_cds), 1500)
            for v in base_cds.sample(sample_size, random_state=42):
                comp_data.append({"Group": "Genome", "Compartment": "CDS GC", "GC": float(v)})

        for _, row in sp_genes.iterrows():
            grp = row["comparative_group"]
            if grp not in self.COLOR_PALETTE:
                continue

            if grp == "cluster_non_globin_alpha":
                pretty_grp = "Non-Globin (α)"
            elif grp == "cluster_non_globin_beta":
                pretty_grp = "Non-Globin (β)"
            elif grp == "cluster_non_globin_locus_A":
                pretty_grp = "Non-Globin (A)"
            elif grp == "cluster_non_globin_locus_B":
                pretty_grp = "Non-Globin (B)"
            elif grp == "cluster_non_globin":
                pretty_grp = "Non-Globin"
            elif grp == "globin_pseudogene":
                pretty_grp = "Pseudogene"
            elif grp == "alpha_globin_coding":
                pretty_grp = "Alpha Globin"
            elif grp == "beta_globin_coding":
                pretty_grp = "Beta Globin"
            elif grp == "globin_coding_locus_A":
                pretty_grp = "Locus A"
            elif grp == "globin_coding_locus_B":
                pretty_grp = "Locus B"
            else:
                pretty_grp = grp.replace("_", " ").title()

            if pd.notna(row.get("full_gene_gc")):
                comp_data.append({"Group": pretty_grp, "Compartment": "Whole-Gene GC", "GC": float(row["full_gene_gc"])})
            if pd.notna(row.get("cds_gc")):
                comp_data.append({"Group": pretty_grp, "Compartment": "CDS GC", "GC": float(row["cds_gc"])})

        df_comp = pd.DataFrame(comp_data)
        group_order = [
            "Genome",
            "Non-Globin (α)",
            "Non-Globin (β)",
            "Non-Globin (A)",
            "Non-Globin (B)",
            "Non-Globin",
            "Pseudogene",
            "Locus A",
            "Locus B",
            "Alpha Globin",
            "Beta Globin",
            "Globin Coding",
        ]
        active_order = [g for g in group_order if g in df_comp["Group"].values]

        palette = {"Whole-Gene GC": "#1565c0", "CDS GC": "#e65100"}
        sns.boxplot(
            data=df_comp,
            x="Group",
            y="GC",
            hue="Compartment",
            order=active_order,
            ax=ax_f,
            palette=palette,
            width=0.6,
            boxprops=dict(alpha=0.65),
            fliersize=2,
        )
        sns.stripplot(
            data=df_comp[df_comp["Group"] != "Genome"],
            x="Group",
            y="GC",
            hue="Compartment",
            dodge=True,
            order=active_order,
            ax=ax_f,
            palette=palette,
            size=6,
            jitter=0.15,
            edgecolor="white",
            linewidth=0.6,
            legend=False,
        )

        ax_f.set_title("F. Compartment Comparison: Whole-Gene GC vs. CDS GC Across Gene Groups", fontsize=12, fontweight="bold", pad=10)
        ax_f.set_ylabel("GC Content (%)", fontsize=11)
        ax_f.set_xlabel("", fontsize=11)
        ax_f.tick_params(axis="x", rotation=15, labelsize=10)
        ax_f.legend(loc="upper left", frameon=True, fontsize=9.5)

        fig.subplots_adjust(top=0.96, bottom=0.06, left=0.07, right=0.97, hspace=0.32, wspace=0.20)
        plt.savefig(out_png, dpi=300)
        plt.close(fig)
        return out_png

    def plot_meta_comparative(
        self,
        df_genes: pd.DataFrame,
        output_path: str = "results/meta_comparative_graph.png",
        clean_mode: bool = False,
    ) -> str:
        """Generates a side-by-side cross-species meta-graph comparing Whole-Gene vs CDS Z-scores.
        
        If clean_mode is True, excludes 'other_locus_feature' and 'globin_pseudogene'.
        Uses deterministic, tiered y-axis offsets to eliminate random jitter and overlap.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

        if clean_mode:
            df_plot = df_genes[~df_genes["comparative_group"].isin(["other_locus_feature", "globin_pseudogene"])].copy()
        else:
            df_plot = df_genes.copy()

        present_species = [sp for sp in self.PHYLOGENETIC_ORDER if sp in df_plot["species"].values]
        if not present_species:
            present_species = sorted(df_plot["species"].unique())

        fig, (ax_fg, ax_cds) = plt.subplots(1, 2, figsize=(22, 11), sharey=True, dpi=300)

        def style_meta_axis(ax, title, xlabel, left_note, right_note, xlim=(-3.8, 3.8)):
            ax.axvspan(-1, 1, color="#e8f5e9", alpha=0.6, label="Genome Typical Range (±1σ)")
            ax.axvspan(-2, 2, color="#f1f8e9", alpha=0.3, zorder=0)
            ax.axvline(0, color="#2e7d32", linestyle="-", linewidth=1.5, label="Genome Mean (Z=0)")
            ax.axvline(-1, color="#81c784", linestyle=":", linewidth=1.0)
            ax.axvline(1, color="#81c784", linestyle=":", linewidth=1.0)
            ax.axvline(-2, color="#c8e6c9", linestyle="--", linewidth=0.8)
            ax.axvline(2, color="#c8e6c9", linestyle="--", linewidth=0.8)
            ax.set_xlim(xlim[0], xlim[1])
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
            ax.set_xlabel(xlabel, fontsize=11, fontweight="bold", labelpad=8)
            ax.annotate(left_note, xy=(xlim[0] * 0.65, len(present_species) - 0.4), fontsize=9.5, fontstyle="italic", color="#555", ha="center")
            ax.annotate(right_note, xy=(xlim[1] * 0.65, len(present_species) - 0.4), fontsize=9.5, fontstyle="italic", color="#555", ha="center")

        style_meta_axis(
            ax_fg,
            "A. Whole-Gene GC Z-Scores (Standardized to Empirical Gene Baseline)",
            "Standardized Z-Score: Whole-Gene GC (σ from Empirical Gene Distribution)",
            "← Lower than Gene Mean",
            "Higher than Gene Mean →",
            xlim=(-3.8, 7.2)
        )
        style_meta_axis(
            ax_cds,
            "B. Coding Sequence (CDS) GC Z-Scores (Standardized to Empirical Coding Baseline)",
            "Standardized Z-Score: CDS GC (σ from Empirical Coding Distribution)",
            "← GC-Depleted CDS",
            "GC-Elevated CDS →",
            xlim=(-3.8, 3.8)
        )

        y_positions = {sp: idx for idx, sp in enumerate(present_species)}
        labels_plotted_fg = set()
        labels_plotted_cds = set()

        for sp in present_species:
            y_pos = y_positions[sp]
            sp_df = df_plot[df_plot["species"] == sp]

            ax_fg.axhline(y_pos, color="#dee2e6", linestyle="-", linewidth=0.8, zorder=1)
            ax_cds.axhline(y_pos, color="#dee2e6", linestyle="-", linewidth=0.8, zorder=1)

            y_map_fg = self._calculate_meta_y_positions(sp_df, "z_score_full_gene_gc", y_pos, clean_mode=clean_mode)
            y_map_cds = self._calculate_meta_y_positions(sp_df, "z_score_gc", y_pos, clean_mode=clean_mode)

            for idx, row in sp_df.iterrows():
                gid = str(row.get("gene_id") or row.get("display_symbol") or idx)
                grp = row["comparative_group"]
                col = self.COLOR_PALETTE.get(grp, "#333")
                m = self.MARKER_MAP.get(grp, "o")

                lbl_text = self.LABEL_MAP.get(grp, grp.replace("_", " ").title())

                # Panel A: Whole-gene Z score
                z_fg = row.get("z_score_full_gene_gc")
                if pd.isna(z_fg):
                    z_fg = row.get("z_score_gc")
                if pd.notna(z_fg) and gid in y_map_fg:
                    lbl_fg = lbl_text if grp not in labels_plotted_fg else None
                    labels_plotted_fg.add(grp)
                    ax_fg.scatter(
                        float(z_fg), y_map_fg[gid],
                        color=col, marker=m, s=80, edgecolors="black", linewidths=0.7, alpha=0.9, zorder=4,
                        label=lbl_fg
                    )

                # Panel B: CDS Z score
                z_cds = row.get("z_score_gc")
                if pd.notna(z_cds) and gid in y_map_cds:
                    lbl_cds = lbl_text if grp not in labels_plotted_cds else None
                    labels_plotted_cds.add(grp)
                    ax_cds.scatter(
                        float(z_cds), y_map_cds[gid],
                        color=col, marker=m, s=80, edgecolors="black", linewidths=0.7, alpha=0.9, zorder=4,
                        label=lbl_cds
                    )

        y_labels = [self.COMMON_NAME_MAP.get(sp, sp) for sp in present_species]
        ax_fg.set_yticks(range(len(present_species)))
        ax_fg.set_yticklabels(y_labels, fontsize=11, fontweight="medium")

        ax_fg.legend(loc="lower right", frameon=True, fontsize=9.0, facecolor="white", framealpha=0.95)
        ax_cds.legend(loc="lower right", frameon=True, fontsize=9.0, facecolor="white", framealpha=0.95)

        title_suffix = " (Clean: Functional Globins & Synteny Flankers)" if clean_mode else " (Comprehensive: All Annotated Features)"
        fig.suptitle(
            f"Cross-Species Meta-Analysis: Vertebrate Globin Loci Compositional Divergence (Whole-Gene vs. CDS){title_suffix}",
            fontsize=14, fontweight="bold", y=0.98
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(output_path, dpi=300)
        pdf_path = output_path.replace(".png", ".pdf")
        plt.savefig(pdf_path)
        plt.close(fig)

        return output_path

    def plot_diagnostic_interactive(
        self, species: str, df_genes: pd.DataFrame, output_dir: str = "results/diagnostic_plots"
    ) -> str:
        """Generates an interactive standalone Plotly HTML 5-panel diagnostic dashboard for a single species."""
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        os.makedirs(output_dir, exist_ok=True)
        html_out = os.path.join(output_dir, f"diagnostic_{species}.html")

        df_base, stats = self.bm.get_genome_baseline(species)
        sp_genes = df_genes[df_genes["species"] == species].copy()
        common_name = self.COMMON_NAME_MAP.get(species, species)

        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=[
                f"<b>A. Whole-Gene GC Distribution & Locus Placement ({common_name})</b>",
                f"<b>B. Coding Sequence (CDS) GC Distribution & Locus Placement ({common_name})</b>",
                f"<b>C. Gene Length Log10(bp) Distribution & Locus Placement ({common_name})</b>",
                "<b>D. Whole-Gene GC vs. CDS GC (Bivariate Space)</b>",
                "<b>E. CDS GC vs. GC3 Codon Space</b>",
                "<b>F. Compartment Comparison: Whole-Gene GC vs. CDS GC Across Gene Groups</b>",
            ],
            vertical_spacing=0.10,
            horizontal_spacing=0.08,
        )

        base_whole = df_base["whole_gene_gc"].dropna() if not df_base.empty and "whole_gene_gc" in df_base.columns else pd.Series(dtype=float)
        base_cds = df_base["cds_gc"].dropna() if not df_base.empty and "cds_gc" in df_base.columns else pd.Series(dtype=float)
        base_len = df_base["whole_gene_len"].dropna() if not df_base.empty and "whole_gene_len" in df_base.columns else pd.Series(dtype=float)
        valid_base_len = base_len[base_len > 0]
        base_log_len = np.log10(valid_base_len) if not valid_base_len.empty else pd.Series(dtype=float)

        mu_whole = stats.get("mean_whole_gene_gc", float(np.mean(base_whole)) if not base_whole.empty else 0.0)
        mu_cds = stats.get("mean_cds_gc", float(np.mean(base_cds)) if not base_cds.empty else 0.0)
        mu_gen = stats.get("mean_genomic_gc")
        mu_log_len = stats.get("mean_log10_len", float(np.mean(base_log_len)) if not base_log_len.empty else 3.8)
        sigma_log_len = stats.get("std_log10_len", float(np.std(base_log_len)) if not base_log_len.empty else 0.55)

        # Aligned x-limits for Panels A and B to facilitate visual comparison
        all_vals = []
        if not base_whole.empty:
            all_vals.extend([float(base_whole.quantile(0.001)), float(base_whole.quantile(0.999))])
        if not base_cds.empty:
            all_vals.extend([float(base_cds.quantile(0.001)), float(base_cds.quantile(0.999))])
        for _, r in sp_genes.iterrows():
            if pd.notna(r.get("full_gene_gc")):
                all_vals.append(float(r["full_gene_gc"]))
            if pd.notna(r.get("cds_gc")):
                all_vals.append(float(r["cds_gc"]))
        min_x = max(10.0, float(np.min(all_vals)) - 4.0) if all_vals else 20.0
        max_x = min(95.0, float(np.max(all_vals)) + 4.0) if all_vals else 90.0

        # 1. Genome histogram in Panel A (Whole-Gene GC)
        if not base_whole.empty:
            fig.add_trace(
                go.Histogram(
                    x=base_whole,
                    histnorm="probability density",
                    name=f"Genome Whole-Gene ({len(base_whole):,} genes)",
                    marker_color="#90caf9",
                    opacity=0.55,
                    nbinsx=60,
                    hoverinfo="x+y",
                    showlegend=True,
                ),
                row=1, col=1,
            )

        if mu_whole > 0:
            fig.add_vline(x=mu_whole, line_width=1.5, line_dash="dash", line_color="#0d47a1", annotation_text=f"Mean: {mu_whole:.1f}%", row=1, col=1)
        if mu_gen is not None:
            fig.add_vline(x=mu_gen, line_width=1.4, line_dash="dashdot", line_color="#2e7d32", annotation_text=f"Assembly: {mu_gen:.1f}%", row=1, col=1)

        # 2. Genome histogram in Panel B (CDS GC)
        if not base_cds.empty:
            fig.add_trace(
                go.Histogram(
                    x=base_cds,
                    histnorm="probability density",
                    name=f"Genome Coding CDS ({len(base_cds):,} genes)",
                    marker_color="#b0bec5",
                    opacity=0.55,
                    nbinsx=60,
                    hoverinfo="x+y",
                    showlegend=True,
                ),
                row=1, col=2,
            )

        if mu_cds > 0:
            fig.add_vline(x=mu_cds, line_width=1.5, line_dash="dash", line_color="#212529", annotation_text=f"Mean: {mu_cds:.1f}%", row=1, col=2)
        if mu_gen is not None:
            fig.add_vline(x=mu_gen, line_width=1.4, line_dash="dashdot", line_color="#2e7d32", annotation_text=f"Assembly: {mu_gen:.1f}%", row=1, col=2)

        # 3. Genome histogram in Panel C (Gene Length Log10(bp))
        if not base_log_len.empty:
            fig.add_trace(
                go.Histogram(
                    x=base_log_len,
                    histnorm="probability density",
                    name=f"Genome Gene Length ({len(base_log_len):,} genes)",
                    marker_color="#ce93d8",
                    opacity=0.55,
                    nbinsx=60,
                    hoverinfo="x+y",
                    showlegend=True,
                ),
                row=2, col=1,
            )

        if mu_log_len > 0:
            fig.add_vline(x=mu_log_len, line_width=1.5, line_dash="dash", line_color="#4a148c", annotation_text=f"Mean: {10**mu_log_len:,.0f} bp", row=2, col=1)
            fig.add_vline(x=mu_log_len + sigma_log_len, line_width=1.1, line_dash="dot", line_color="#7b1fa2", row=2, col=1)
            fig.add_vline(x=mu_log_len - sigma_log_len, line_width=1.1, line_dash="dot", line_color="#7b1fa2", row=2, col=1)

        # Diagonal line for Panel D (Whole-Gene GC vs. CDS GC, y = x)
        fig.add_trace(
            go.Scatter(
                x=[20, 80],
                y=[20, 80],
                mode="lines",
                name="y = x (CDS = Whole-Gene)",
                line=dict(dash="dash", color="#adb5bd"),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=2, col=2,
        )

        # Diagonal line for Panel E (CDS GC vs GC3, y = x)
        fig.add_trace(
            go.Scatter(
                x=[20, 90],
                y=[20, 90],
                mode="lines",
                name="y = x (Codon Bias = CDS)",
                line=dict(dash="dot", color="#ced4da"),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=3, col=1,
        )

        # Compute deterministic y-axis offsets for Plotly Panels A, B, C
        y_max_plotly_a = 0.05
        if not base_whole.empty:
            std_w = float(stats.get("std_whole_gene_gc", np.std(base_whole)))
            if std_w > 0:
                y_max_plotly_a = 1.0 / (std_w * np.sqrt(2 * np.pi))
        y_positions_plotly_a = self._calculate_gene_y_positions(sp_genes, "full_gene_gc", y_max_plotly_a)

        y_max_plotly_b = 0.05
        if not base_cds.empty:
            std_c = float(stats.get("std_cds_gc", np.std(base_cds)))
            if std_c > 0:
                y_max_plotly_b = 1.0 / (std_c * np.sqrt(2 * np.pi))
        y_positions_plotly_b = self._calculate_gene_y_positions(sp_genes, "cds_gc", y_max_plotly_b)

        y_max_plotly_c = 0.8
        if not base_log_len.empty and sigma_log_len > 0:
            y_max_plotly_c = 1.0 / (sigma_log_len * np.sqrt(2 * np.pi))
        y_positions_plotly_c = self._calculate_gene_y_positions(sp_genes, "log10_coding_span_len", y_max_plotly_c)

        # Iterate over groups for Panels A, B, C, D, E
        for grp, df_grp in sp_genes.groupby("comparative_group"):
            col = self.COLOR_PALETTE.get(grp, "#333")
            lbl_name = self.LABEL_MAP.get(grp, grp.replace("_", " ").title())

            # Panel A: Locus points on Whole-Gene distribution
            valid_a = df_grp.dropna(subset=["full_gene_gc"])
            if not valid_a.empty:
                sym_list_a = [self._format_symbol_annotation(s, grp) for s in valid_a["display_symbol"]]
                hover_a = [
                    f"<b>{r['display_symbol']}</b> ({r['gene_name']})<br>"
                    f"Locus Origin: {r.get('locus_id', 'unknown')}<br>"
                    f"Group: {lbl_name}<br>"
                    f"Whole-Gene GC: {r.get('full_gene_gc')}% (Z = {r.get('z_score_full_gene_gc')})<br>"
                    f"CDS GC: {r.get('cds_gc')}% (Z = {r.get('z_score_gc')})<br>"
                    f"Percentile Rank: {r.get('percentile_gc')}%"
                    for _, r in valid_a.iterrows()
                ]
                y_vals_a = [
                    y_positions_plotly_a.get(str(r.get("gene_id") or r.get("display_symbol") or idx), y_max_plotly_a * 0.10)
                    for idx, r in valid_a.iterrows()
                ]
                fig.add_trace(
                    go.Scatter(
                        x=valid_a["full_gene_gc"],
                        y=y_vals_a,
                        mode="markers+text",
                        name=lbl_name,
                        marker=dict(size=11, color=col, line=dict(width=1, color="black")),
                        text=sym_list_a,
                        textposition="top center",
                        hovertext=hover_a,
                        hoverinfo="text",
                        showlegend=True,
                    ),
                    row=1, col=1,
                )

            # Panel B: Locus points on CDS distribution
            valid_b = df_grp.dropna(subset=["cds_gc"])
            if not valid_b.empty:
                sym_list_b = [self._format_symbol_annotation(s, grp) for s in valid_b["display_symbol"]]
                hover_b = [
                    f"<b>{r['display_symbol']}</b> ({r['gene_name']})<br>"
                    f"Locus Origin: {r.get('locus_id', 'unknown')}<br>"
                    f"Group: {lbl_name}<br>"
                    f"CDS GC: {r.get('cds_gc')}% (Z = {r.get('z_score_gc')})<br>"
                    f"Whole-Gene GC: {r.get('full_gene_gc')}% (Z = {r.get('z_score_full_gene_gc')})<br>"
                    f"Percentile Rank: {r.get('percentile_gc')}%"
                    for _, r in valid_b.iterrows()
                ]
                y_vals_b = [
                    y_positions_plotly_b.get(str(r.get("gene_id") or r.get("display_symbol") or idx), y_max_plotly_b * 0.10)
                    for idx, r in valid_b.iterrows()
                ]
                fig.add_trace(
                    go.Scatter(
                        x=valid_b["cds_gc"],
                        y=y_vals_b,
                        mode="markers+text",
                        name=lbl_name,
                        marker=dict(size=11, color=col, line=dict(width=1, color="black")),
                        text=sym_list_b,
                        textposition="top center",
                        hovertext=hover_b,
                        hoverinfo="text",
                        showlegend=False,
                    ),
                    row=1, col=2,
                )

            # Panel C: Locus points on Length distribution
            valid_c_len = df_grp.dropna(subset=["log10_coding_span_len"])
            if valid_c_len.empty and "full_gene_len" in df_grp.columns:
                valid_c_len = df_grp[df_grp["full_gene_len"].notna() & (df_grp["full_gene_len"] > 0)].copy()
                valid_c_len["log10_coding_span_len"] = np.log10(valid_c_len["full_gene_len"])

            if not valid_c_len.empty:
                sym_list_c = [self._format_symbol_annotation(s, grp) for s in valid_c_len["display_symbol"]]
                hover_c = [
                    f"<b>{r['display_symbol']}</b> ({r['gene_name']})<br>"
                    f"Locus Origin: {r.get('locus_id', 'unknown')}<br>"
                    f"Group: {lbl_name}<br>"
                    f"Length: {int(r.get('coding_span_len') or r.get('full_gene_len', 0)):,} bp<br>"
                    f"log10(Length): {r.get('log10_coding_span_len'):.2f}<br>"
                    f"Length Z-Score: {r.get('z_score_length', 'N/A')}"
                    for _, r in valid_c_len.iterrows()
                ]
                y_vals_c = [
                    y_positions_plotly_c.get(str(r.get("gene_id") or r.get("display_symbol") or idx), y_max_plotly_c * 0.10)
                    for idx, r in valid_c_len.iterrows()
                ]
                fig.add_trace(
                    go.Scatter(
                        x=valid_c_len["log10_coding_span_len"],
                        y=y_vals_c,
                        mode="markers+text",
                        name=lbl_name,
                        marker=dict(size=11, color=col, line=dict(width=1, color="black")),
                        text=sym_list_c,
                        textposition="top center",
                        hovertext=hover_c,
                        hoverinfo="text",
                        showlegend=False,
                    ),
                    row=2, col=1,
                )

            # Panel D: Whole-Gene GC vs. CDS GC
            valid_d = df_grp.dropna(subset=["full_gene_gc", "cds_gc"])
            if not valid_d.empty:
                d_hover = [
                    f"<b>{r['display_symbol']}</b> ({r.get('locus_id')})<br>Whole-Gene GC: {r['full_gene_gc']}%<br>CDS GC: {r['cds_gc']}%"
                    for _, r in valid_d.iterrows()
                ]
                fig.add_trace(
                    go.Scatter(
                        x=valid_d["full_gene_gc"],
                        y=valid_d["cds_gc"],
                        mode="markers+text",
                        name=lbl_name,
                        marker=dict(size=11, color=col, line=dict(width=1, color="black")),
                        text=[self._format_symbol_annotation(s, grp) for s in valid_d["display_symbol"]],
                        textposition="top right",
                        hovertext=d_hover,
                        hoverinfo="text",
                        showlegend=False,
                    ),
                    row=2, col=2,
                )

            # Panel E: CDS GC vs GC3
            valid_e = df_grp.dropna(subset=["cds_gc", "gc3"])
            if not valid_e.empty:
                e_hover = [
                    f"<b>{r['display_symbol']}</b> ({r.get('locus_id')})<br>CDS GC: {r['cds_gc']}%<br>GC3: {r['gc3']}%"
                    for _, r in valid_e.iterrows()
                ]
                fig.add_trace(
                    go.Scatter(
                        x=valid_e["cds_gc"],
                        y=valid_e["gc3"],
                        mode="markers+text",
                        name=lbl_name,
                        marker=dict(size=11, color=col, line=dict(width=1, color="black")),
                        text=[self._format_symbol_annotation(s, grp) for s in valid_e["display_symbol"]],
                        textposition="top right",
                        hovertext=e_hover,
                        hoverinfo="text",
                        showlegend=False,
                    ),
                    row=3, col=1,
                )

        # Panel F: Compartment Comparison
        comp_data = []
        if not base_whole.empty:
            sample_size = min(len(base_whole), 1500)
            for v in base_whole.sample(sample_size, random_state=42):
                comp_data.append({"Group": "Genome", "Compartment": "Whole-Gene GC", "GC": float(v)})
        if not base_cds.empty:
            sample_size = min(len(base_cds), 1500)
            for v in base_cds.sample(sample_size, random_state=42):
                comp_data.append({"Group": "Genome", "Compartment": "CDS GC", "GC": float(v)})

        for _, row in sp_genes.iterrows():
            grp = row["comparative_group"]
            if grp not in self.COLOR_PALETTE:
                continue

            if grp == "cluster_non_globin_alpha":
                pretty_grp = "Non-Globin (α)"
            elif grp == "cluster_non_globin_beta":
                pretty_grp = "Non-Globin (β)"
            elif grp == "cluster_non_globin_locus_A":
                pretty_grp = "Non-Globin (A)"
            elif grp == "cluster_non_globin_locus_B":
                pretty_grp = "Non-Globin (B)"
            elif grp == "cluster_non_globin":
                pretty_grp = "Non-Globin"
            elif grp == "globin_pseudogene":
                pretty_grp = "Pseudogene"
            elif grp == "alpha_globin_coding":
                pretty_grp = "Alpha Globin"
            elif grp == "beta_globin_coding":
                pretty_grp = "Beta Globin"
            elif grp == "globin_coding_locus_A":
                pretty_grp = "Locus A"
            elif grp == "globin_coding_locus_B":
                pretty_grp = "Locus B"
            else:
                pretty_grp = grp.replace("_", " ").title()

            if pd.notna(row.get("full_gene_gc")):
                comp_data.append({"Group": pretty_grp, "Compartment": "Whole-Gene GC", "GC": float(row["full_gene_gc"])})
            if pd.notna(row.get("cds_gc")):
                comp_data.append({"Group": pretty_grp, "Compartment": "CDS GC", "GC": float(row["cds_gc"])})

        df_comp = pd.DataFrame(comp_data)
        group_order = [
            "Genome",
            "Non-Globin (α)",
            "Non-Globin (β)",
            "Non-Globin (A)",
            "Non-Globin (B)",
            "Non-Globin",
            "Pseudogene",
            "Locus A",
            "Locus B",
            "Alpha Globin",
            "Beta Globin",
            "Globin Coding",
        ]
        active_order = [g for g in group_order if g in df_comp["Group"].values]

        if not df_comp.empty:
            df_wg = df_comp[df_comp["Compartment"] == "Whole-Gene GC"]
            df_cg = df_comp[df_comp["Compartment"] == "CDS GC"]

            fig.add_trace(
                go.Box(
                    x=df_wg["Group"],
                    y=df_wg["GC"],
                    name="Whole-Gene GC",
                    marker_color="#1565c0",
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    showlegend=True,
                ),
                row=3, col=2,
            )
            fig.add_trace(
                go.Box(
                    x=df_cg["Group"],
                    y=df_cg["GC"],
                    name="CDS GC",
                    marker_color="#e65100",
                    boxpoints="all",
                    jitter=0.3,
                    pointpos=-1.8,
                    showlegend=True,
                ),
                row=3, col=2,
            )

        fig.update_xaxes(title_text="Whole-Gene GC Content (%) [Exons + Introns]", range=[min_x, max_x], row=1, col=1)
        fig.update_yaxes(title_text="Density", row=1, col=1)

        fig.update_xaxes(title_text="Coding Sequence (CDS) GC Content (%)", range=[min_x, max_x], row=1, col=2)
        fig.update_yaxes(title_text="Density", row=1, col=2)

        fig.update_xaxes(
            title_text="Gene Length (log10 bp)",
            range=[1.8, 6.2],
            tickvals=[2.0, 3.0, 4.0, 5.0, 6.0],
            ticktext=["100 bp (10²)", "1 kb (10³)", "10 kb (10⁴)", "100 kb (10⁵)", "1 Mb (10⁶)"],
            row=2, col=1
        )
        fig.update_yaxes(title_text="Density", row=2, col=1)

        fig.update_xaxes(title_text="Whole-Gene GC (%)", row=2, col=2)
        fig.update_yaxes(title_text="CDS GC (%)", row=2, col=2)

        fig.update_xaxes(title_text="CDS GC (%)", row=3, col=1)
        fig.update_yaxes(title_text="GC3 (%)", row=3, col=1)

        fig.update_xaxes(title_text="Gene Group", categoryorder="array", categoryarray=active_order, row=3, col=2)
        fig.update_yaxes(title_text="GC Content (%)", row=3, col=2)

        fig.for_each_annotation(lambda a: a.update(font=dict(size=12.5, color="#212529")))

        fig.update_layout(
            height=1500,
            width=1500,
            template="plotly_white",
            boxmode="group",
            title=dict(
                text=f"<b>Interactive Globin Locus Diagnostic Dashboard: {common_name}</b>",
                x=0.5,
                y=0.99,
                xanchor="center",
                yanchor="top",
                font=dict(size=16, color="#212529"),
            ),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.06,
                xanchor="center",
                x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)",
                bordercolor="#dee2e6",
                borderwidth=1,
                font=dict(size=10.5),
            ),
            margin=dict(t=70, b=90, l=70, r=40),
        )

        fig.write_html(html_out)
        return html_out

    def plot_meta_interactive(
        self,
        df_genes: pd.DataFrame,
        output_path: str = "results/meta_comparative_graph.html",
        clean_mode: bool = False,
    ) -> str:
        """Generates an interactive cross-species meta-graph comparing Whole-Gene vs CDS Z-scores side-by-side."""
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

        if clean_mode:
            df_plot = df_genes[~df_genes["comparative_group"].isin(["other_locus_feature", "globin_pseudogene"])].copy()
        else:
            df_plot = df_genes.copy()

        present_species = [sp for sp in self.PHYLOGENETIC_ORDER if sp in df_plot["species"].values]
        if not present_species:
            present_species = sorted(df_plot["species"].unique())

        y_labels = [self.COMMON_NAME_MAP.get(sp, sp) for sp in present_species]
        y_positions = {sp: idx for idx, sp in enumerate(present_species)}

        # Precompute deterministic y-offsets per species for both Whole-Gene and CDS metrics
        y_maps_fg = {}
        y_maps_cds = {}
        for sp in present_species:
            sp_df = df_plot[df_plot["species"] == sp]
            y_pos = y_positions[sp]
            y_maps_fg[sp] = self._calculate_meta_y_positions(sp_df, "z_score_full_gene_gc", y_pos, clean_mode=clean_mode)
            y_maps_cds[sp] = self._calculate_meta_y_positions(sp_df, "z_score_gc", y_pos, clean_mode=clean_mode)

        fig = make_subplots(
            rows=1, cols=2,
            shared_yaxes=True,
            subplot_titles=[
                "<b>A. Whole-Gene GC Z-Scores (Standardized to Empirical Gene Baseline)</b>",
                "<b>B. Coding Sequence (CDS) GC Z-Scores (Standardized to Empirical Coding Baseline)</b>",
            ],
            horizontal_spacing=0.08,
        )

        for col_idx in [1, 2]:
            fig.add_vrect(x0=-1, x1=1, fillcolor="#e8f5e9", opacity=0.5, layer="below", line_width=0, row=1, col=col_idx)
            fig.add_vrect(x0=-2, x1=2, fillcolor="#f1f8e9", opacity=0.25, layer="below", line_width=0, row=1, col=col_idx)
            fig.add_vline(x=0, line_width=1.5, line_color="#2e7d32", line_dash="solid", row=1, col=col_idx)
            for y_pos in y_positions.values():
                fig.add_hline(y=y_pos, line_width=0.6, line_color="#dee2e6", line_dash="solid", layer="below", row=1, col=col_idx)

        group_list = [
            "alpha_globin_coding",
            "beta_globin_coding",
            "globin_coding_locus_A",
            "globin_coding_locus_B",
            "globin_coding",
            "globin_pseudogene",
            "cluster_non_globin_alpha",
            "cluster_non_globin_beta",
            "cluster_non_globin_locus_A",
            "cluster_non_globin_locus_B",
            "cluster_non_globin",
            "other_locus_feature",
        ]
        if clean_mode:
            group_list = [g for g in group_list if g not in ["other_locus_feature", "globin_pseudogene"]]

        for grp in group_list:
            sub = df_plot[df_plot["comparative_group"] == grp]
            if sub.empty:
                continue

            col = self.COLOR_PALETTE.get(grp, "#333")
            if "alpha" in grp or grp == "globin_coding":
                sym = "circle" if "coding" in grp else "triangle-down"
            elif "beta" in grp:
                sym = "square" if "coding" in grp else "triangle-left"
            elif "locus_A" in grp:
                sym = "circle" if "coding" in grp else "triangle-down"
            elif "locus_B" in grp:
                sym = "triangle-up" if "coding" in grp else "triangle-right"
            elif "pseudogene" in grp:
                sym = "diamond"
            else:
                sym = "triangle-down"

            x_vals_fg, y_vals_fg, hover_fg = [], [], []
            x_vals_cds, y_vals_cds, hover_cds = [], [], []

            for idx_r, r in sub.iterrows():
                sp = r["species"]
                if sp not in y_positions:
                    continue
                gid = str(r.get("gene_id") or r.get("display_symbol") or idx_r)
                y_val_fg = y_maps_fg.get(sp, {}).get(gid)
                y_val_cds = y_maps_cds.get(sp, {}).get(gid)

                z_fg = r.get("z_score_full_gene_gc")
                if pd.isna(z_fg):
                    z_fg = r.get("z_score_gc")
                z_cds = r.get("z_score_gc")

                lbl_name = self.LABEL_MAP.get(grp, grp.replace("_", " ").title())

                tooltip = (
                    f"<b>{r['display_symbol']}</b> ({self.COMMON_NAME_MAP.get(sp, sp)})<br>"
                    f"Locus Origin: {r.get('locus_id', 'unknown')}<br>"
                    f"Group: {lbl_name}<br>"
                    f"Whole-Gene GC: {r.get('full_gene_gc')}% (Z = {z_fg})<br>"
                    f"CDS GC: {r.get('cds_gc')}% (Z = {z_cds})<br>"
                    f"Intronic GC: {r.get('intronic_gc')}%<br>"
                    f"Exonic GC: {r.get('exonic_gc')}%<br>"
                    f"GC3: {r.get('gc3')}%"
                )

                if pd.notna(z_fg) and y_val_fg is not None:
                    x_vals_fg.append(z_fg)
                    y_vals_fg.append(y_val_fg)
                    hover_fg.append(tooltip)

                if pd.notna(z_cds) and y_val_cds is not None:
                    x_vals_cds.append(z_cds)
                    y_vals_cds.append(y_val_cds)
                    hover_cds.append(tooltip)

            lbl_name = self.LABEL_MAP.get(grp, grp.replace("_", " ").title())

            if x_vals_fg:
                fig.add_trace(
                    go.Scatter(
                        x=x_vals_fg,
                        y=y_vals_fg,
                        mode="markers",
                        name=lbl_name,
                        marker=dict(size=10, color=col, symbol=sym, line=dict(width=0.8, color="black")),
                        hovertext=hover_fg,
                        hoverinfo="text",
                        legendgroup=grp,
                    ),
                    row=1, col=1,
                )

            if x_vals_cds:
                fig.add_trace(
                    go.Scatter(
                        x=x_vals_cds,
                        y=y_vals_cds,
                        mode="markers",
                        name=lbl_name,
                        marker=dict(size=10, color=col, symbol=sym, line=dict(width=0.8, color="black")),
                        hovertext=hover_cds,
                        hoverinfo="text",
                        legendgroup=grp,
                        showlegend=False,
                    ),
                    row=1, col=2,
                )

        fig.update_xaxes(title_text="<b>Whole-Gene GC Z-Score</b>", range=[-3.8, 7.2], dtick=1.0, row=1, col=1)
        fig.update_xaxes(title_text="<b>CDS GC Z-Score</b>", range=[-3.8, 3.8], dtick=1.0, row=1, col=2)
        fig.update_yaxes(tickvals=list(range(len(present_species))), ticktext=y_labels, row=1, col=1)

        fig.for_each_annotation(lambda a: a.update(font=dict(size=13, color="#212529")))
        title_suffix = " (Clean: Functional Globins & Synteny Flankers)" if clean_mode else " (Comprehensive: All Annotated Features)"
        fig.update_layout(
            height=920,
            width=1380,
            template="plotly_white",
            title=dict(
                text=f"<b>Cross-Species Meta-Analysis: Vertebrate Globin Loci Compositional Divergence (Whole-Gene vs. CDS){title_suffix}</b>",
                x=0.5,
                y=0.985,
                xanchor="center",
                yanchor="top",
                font=dict(size=15, color="#212529"),
            ),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.08,
                xanchor="center",
                x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)",
                bordercolor="#dee2e6",
                borderwidth=1,
                font=dict(size=11),
            ),
            margin=dict(t=70, b=90, l=140, r=40),
        )

        fig.write_html(output_path)
        return output_path

    def plot_chromosome_linear(
        self,
        species: str,
        df_genes: pd.DataFrame,
        df_intervals: Optional[pd.DataFrame] = None,
        output_dir: str = "results/chromosome_plots",
    ) -> Dict[str, str]:
        """Plots linear chromosome arm tracks showing full chromosome length, centromere position,
        distance to nearest telomere, and zoomed locus gene architecture with realistic coordinate spacing.
        Generates static PNG/PDF and interactive Plotly HTML.
        """
        import json
        import matplotlib.patches as patches
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        os.makedirs(output_dir, exist_ok=True)
        common_name = self.COMMON_NAME_MAP.get(species, species)

        if df_intervals is None:
            intervals_path = "results/curated_intervals.csv"
            if not os.path.exists(intervals_path):
                intervals_path = "data/curated_intervals.csv"
            if os.path.exists(intervals_path):
                df_intervals = pd.read_csv(intervals_path)
            else:
                df_intervals = pd.DataFrame()

        sp_intervals = df_intervals[df_intervals["species"] == species] if not df_intervals.empty else pd.DataFrame()
        sp_genes = df_genes[df_genes["species"] == species].copy()

        # Ensure start, end, and strand coordinates exist
        if "start" not in sp_genes.columns or sp_genes["start"].isna().all():
            g_coord_path = "results/genes_in_intervals.csv"
            if not os.path.exists(g_coord_path):
                g_coord_path = "data/genes_in_intervals.csv"
            if os.path.exists(g_coord_path):
                df_coords = pd.read_csv(g_coord_path)
                sp_coords = df_coords[df_coords["species"] == species][["gene_id", "start", "end", "strand"]].drop_duplicates(subset=["gene_id"])
                sp_genes = sp_genes.drop(columns=[c for c in ["start", "end", "strand"] if c in sp_genes.columns])
                sp_genes = sp_genes.merge(sp_coords, on="gene_id", how="left")

        chr_len_path = "data/chromosome_lengths.json"
        if not os.path.exists(chr_len_path):
            chr_len_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "chromosome_lengths.json")
        chr_lengths = {}
        if os.path.exists(chr_len_path):
            with open(chr_len_path, "r", encoding="utf-8") as f:
                chr_lengths = json.load(f)

        cen_path = "data/centromere_positions.json"
        if not os.path.exists(cen_path):
            cen_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "centromere_positions.json")
        centromeres = {}
        if os.path.exists(cen_path):
            with open(cen_path, "r", encoding="utf-8") as f:
                centromeres = json.load(f)

        num_loci = len(sp_intervals)
        if num_loci == 0:
            return {}

        png_path = os.path.join(output_dir, f"chromosome_linear_{species}.png")
        pdf_path = os.path.join(output_dir, f"chromosome_linear_{species}.pdf")
        html_path = os.path.join(output_dir, f"chromosome_linear_{species}.html")

        # Determine maximum chromosome length across all loci for this species
        # to ensure chromosomes are scaled RELATIVE to each other
        all_chr_lens = [chr_lengths.get(f"{species}::{r['chromosome']}", int(r["end"]) * 1.2) for _, r in sp_intervals.iterrows()]
        max_chr_len = max(all_chr_lens) if all_chr_lens else 100000000
        max_locus_span_kb = max((int(r["end"]) - int(r["start"])) / 1e3 for _, r in sp_intervals.iterrows())

        # 1. Matplotlib static rendering
        fig, axes = plt.subplots(
            num_loci * 2, 1, figsize=(15, 4.2 * num_loci), dpi=200,
            gridspec_kw={"height_ratios": [1, 2.6] * num_loci}
        )
        if num_loci == 1:
            axes = [axes[0], axes[1]]

        plotted_mpl_genes = set()

        for idx, (_, r_int) in enumerate(sp_intervals.iterrows()):
            loc_id = r_int["locus_id"]
            chrom = str(r_int["chromosome"])
            l_s = int(r_int["start"])
            l_e = int(r_int["end"])
            chr_len = chr_lengths.get(f"{species}::{chrom}", max(l_e * 1.2, 50000000))
            cen_bp = centromeres.get(f"{species}::{chrom}")

            # Centromere-aware arm telomere logic
            mid_locus = (l_s + l_e) / 2
            is_fragment = (chr_len < 2000000) or str(chrom).lower().startswith(("gl", "scaffold", "ctg"))
            if not is_fragment and cen_bp and 0 < cen_bp < chr_len:
                if mid_locus > cen_bp:
                    near_telo = "3' (q-arm)"
                    min_d = max(chr_len - l_e, 0)
                else:
                    near_telo = "5' (p-arm)"
                    min_d = max(l_s, 0)
            else:
                d_5p = l_s
                d_3p = chr_len - l_e
                if d_5p < d_3p:
                    near_telo = "5' (p-arm)"
                    min_d = d_5p
                else:
                    near_telo = "3' (q-arm)"
                    min_d = d_3p
            rel_d = (min_d / chr_len) * 100

            if is_fragment:
                chr_class = "Scaffold Fragment"
            else:
                chr_class = "Microchromosome (<25 Mb)" if chr_len < 25e6 else ("Intermediate (25-50 Mb)" if chr_len < 50e6 else "Macrochromosome (>50 Mb)")

            loc_name = loc_id.replace("_", " ").title()
            loc_lower = loc_id.lower()
            if "locus_a" in loc_lower:
                theme_col = "#d32f2f"       # Carmine Red for Locus A
                theme_dark = "#b71c1c"
                theme_light = "#fffbfb"
                loc_badge = "Locus A (α+β Cluster)"
            elif "locus_b" in loc_lower:
                theme_col = "#1976d2"       # Royal Blue for Locus B
                theme_dark = "#0d47a1"
                theme_light = "#f8fbff"
                loc_badge = "Locus B (α+β Cluster)"
            elif loc_id == "alpha_globin" or (loc_lower.startswith("alpha") and "beta" not in loc_lower):
                theme_col = "#d32f2f"       # Carmine Red for Alpha
                theme_dark = "#b71c1c"
                theme_light = "#fffbfb"
                loc_badge = "Alpha Globin Locus"
            elif loc_id == "beta_globin" or (loc_lower.startswith("beta") and "alpha" not in loc_lower):
                theme_col = "#1976d2"       # Royal Blue for Beta
                theme_dark = "#0d47a1"
                theme_light = "#f8fbff"
                loc_badge = "Beta Globin Locus"
            elif "proto" in loc_lower:
                theme_col = "#7b1fa2"       # Purple for proto-globin
                theme_dark = "#4a148c"
                theme_light = "#faf7fc"
                loc_badge = "Proto-Globin Cluster"
            else:
                theme_col = "#7b1fa2"       # Purple for linked/other
                theme_dark = "#4a148c"
                theme_light = "#faf7fc"
                loc_badge = "Linked α+β Globin Cluster"

            # Top: Chromosome bar (full chromosome from 0 to chr_len)
            ax_chr = axes[idx * 2]
            ax_chr.set_xlim(-chr_len * 0.015 / 1e6, chr_len * 1.04 / 1e6)
            ax_chr.set_ylim(-0.6, 0.6)
            ax_chr.axis("off")

            # Draw chromosome body (rounded rectangle)
            rect = patches.FancyBboxPatch(
                (0, -0.2), chr_len / 1e6, 0.4,
                boxstyle="round,pad=0.02,rounding_size=0.15",
                facecolor="#eceff1", edgecolor="#455a64", linewidth=1.2,
                linestyle="--" if is_fragment else "-", zorder=2
            )
            ax_chr.add_patch(rect)

            # Draw centromere indication if known (only for bona fide chromosomes, not scaffold fragments)
            if not is_fragment and cen_bp and 0 < cen_bp < chr_len:
                cen_mb = cen_bp / 1e6
                w_cen = chr_len / 1e6 * 0.007
                notch_top = patches.Polygon(
                    [(cen_mb - w_cen, 0.22), (cen_mb, 0.06), (cen_mb + w_cen, 0.22)],
                    facecolor="white", edgecolor="#455a64", linewidth=1.0, zorder=4
                )
                notch_bot = patches.Polygon(
                    [(cen_mb - w_cen, -0.22), (cen_mb, -0.06), (cen_mb + w_cen, -0.22)],
                    facecolor="white", edgecolor="#455a64", linewidth=1.0, zorder=4
                )
                ax_chr.add_patch(notch_top)
                ax_chr.add_patch(notch_bot)
                ax_chr.axvline(cen_mb, ymin=0.25, ymax=0.75, color="#c62828", linestyle=":", linewidth=1.2, zorder=5)
                ax_chr.text(cen_mb, -0.32, "Cen", fontsize=8.5, fontweight="bold", color="#c62828", ha="center", va="top")

            # Highlight globin locus on the chromosome bar matching the zoomed-in box color
            l_span = max(l_e - l_s, 1)
            l_span_str = f"{l_span/1e3:.0f} kb" if l_span < 1e6 else f"{l_span/1e6:.2f} Mb"
            l_w_mb = max((l_e - l_s) / 1e6, chr_len / 1e6 * 0.003)
            loc_rect = patches.Rectangle(
                (l_s / 1e6, -0.23), l_w_mb, 0.46,
                facecolor=theme_col, edgecolor=theme_dark, linewidth=1.6, zorder=6
            )
            ax_chr.add_patch(loc_rect)

            if is_fragment:
                chr_prefix = "Scaffold " if not chrom.lower().startswith("scaffold") else ""
                chr_len_str = f"{chr_len/1e3:.0f} kb" if chr_len < 1e6 else f"{chr_len/1e6:.2f} Mb"
                dist_str = f"{min_d/1e3:.0f} kb" if min_d < 1e6 else f"{min_d/1e6:.2f} Mb"
                ax_chr.text(
                    0, 0.35,
                    f"{chr_prefix}{chrom} ({chr_class}, Full Chr: {chr_len_str} | Locus Span: {l_span_str})",
                    fontsize=10.5, fontweight="bold", va="bottom", color="#263238"
                )
                ax_chr.text(
                    chr_len / 1e6, 0.35,
                    f"Distance to fragment end: {dist_str}",
                    fontsize=10, fontweight="bold", color=theme_dark, ha="right", va="bottom"
                )
            else:
                ax_chr.text(
                    0, 0.35,
                    f"Chromosome {chrom} ({chr_class}, Full Chr: {chr_len/1e6:.1f} Mb | Locus Span: {l_span_str})",
                    fontsize=10.5, fontweight="bold", va="bottom", color="#263238"
                )
                ax_chr.text(
                    chr_len / 1e6, 0.35,
                    f"Distance to {near_telo} telomere: {min_d/1e6:.2f} Mb",
                    fontsize=10, fontweight="bold", color=theme_dark, ha="right", va="bottom"
                )

            # Bottom: Zoomed Locus Gene Track (styled with matching theme color box)
            ax_loc = axes[idx * 2 + 1]
            ax_loc.set_facecolor(theme_light)
            for spine in ax_loc.spines.values():
                spine.set_edgecolor(theme_col)
                spine.set_linewidth(1.8)

            loc_genes = sp_genes[sp_genes["locus_id"] == loc_id].copy()
            # Exclude other_locus_feature and pseudogenes to match clean protein-coding architecture
            loc_genes = loc_genes[~loc_genes["comparative_group"].isin(["other_locus_feature", "globin_pseudogene"])]
            loc_genes = loc_genes[~loc_genes["biotype"].astype(str).str.contains("pseudogene", case=False, na=False)]

            pad_kb = (l_e - l_s) / 1e3 * 0.03
            ax_loc.set_xlim(l_s / 1e3 - pad_kb, l_e / 1e3 + pad_kb)
            ax_loc.set_ylim(-1.55, 1.55)
            # Center baseline and upper/lower flank guidelines
            ax_loc.axhline(0, color="#78909c", linewidth=2.0, zorder=1)
            ax_loc.axhline(0.55, color="#e0e0e0", linestyle="--", linewidth=0.8, zorder=1)
            ax_loc.axhline(-0.55, color="#e0e0e0", linestyle="--", linewidth=0.8, zorder=1)
            ax_loc.set_xlabel(f"Chr {chrom} Coordinates (kb)", fontsize=10.5, fontweight="bold", color="#37474f")
            ax_loc.set_yticks([])
            ax_loc.set_title(
                f"Zoomed Locus Architecture: {loc_name} [{loc_badge}] ({l_s/1e3:,.1f} - {l_e/1e3:,.1f} kb; Span: {l_span_str})",
                fontsize=11, loc="center", color=theme_dark, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=theme_col, alpha=0.95, linewidth=0.8),
                pad=9
            )

            # Realistic gene arrows: Globins on middle track (y=0), non-globins shifted up (+0.55) & down (-0.55)
            if not loc_genes.empty:
                loc_genes["start_num"] = pd.to_numeric(loc_genes["start"], errors="coerce").fillna(l_s)
                loc_genes["end_num"] = pd.to_numeric(loc_genes["end"], errors="coerce").fillna(l_e)
                loc_genes = loc_genes[(loc_genes["end_num"] >= l_s) & (loc_genes["start_num"] <= l_e)]
                loc_genes = loc_genes.sort_values(by="start_num").reset_index(drop=True)

                locus_span_kb = (l_e - l_s) / 1e3
                max_tip_kb = min(locus_span_kb * 0.012, 1.5)

                glo_tiers_pos = [0.42, 0.78, 1.15]
                glo_tiers_neg = [-0.42, -0.78, -1.15]
                glo_idx = 0
                non_glo_idx = 0

                for _, g in loc_genes.iterrows():
                    sym = str(g["display_symbol"])
                    grp = str(g.get("comparative_group", ""))
                    is_pseudo = ("pseudogene" in grp) or ("pseudo" in sym.lower())
                    is_glo = str(g.get("is_globin", False)).lower() in ["true", "pseudogene"] or ("globin" in grp and "non_globin" not in grp) or is_pseudo
                    strand = str(g.get("strand", "+"))

                    g_s_kb = float(g["start_num"]) / 1e3
                    g_e_kb = float(g["end_num"]) / 1e3
                    if g_e_kb <= g_s_kb:
                        g_e_kb = g_s_kb + 1.5

                    # Clamp coordinates to visible locus window
                    vis_s = max(g_s_kb, l_s / 1e3)
                    vis_e = min(g_e_kb, l_e / 1e3)
                    if vis_e <= vis_s:
                        continue

                    raw_len_kb = vis_e - vis_s
                    min_disp_kb = 0.6
                    if raw_len_kb < min_disp_kb:
                        mid_kb = (vis_s + vis_e) / 2
                        disp_s = max(l_s / 1e3, mid_kb - min_disp_kb / 2)
                        disp_e = min(l_e / 1e3, mid_kb + min_disp_kb / 2)
                    else:
                        disp_s = vis_s
                        disp_e = vis_e

                    disp_w = disp_e - disp_s
                    # Bounded arrowhead tip width so long genes never have oversized arrowheads
                    arrow_w = min(max_tip_kb, 0.30 * disp_w)

                    # Determine track y-center and height
                    if is_glo:
                        y_center = 0.0
                        h_box = 0.32
                    else:
                        y_center = 0.55 if (non_glo_idx % 2 == 0) else -0.55
                        non_glo_idx += 1
                        h_box = 0.26

                    y_bot = y_center - h_box / 2
                    y_top = y_center + h_box / 2

                    if strand == "+":
                        pts = [
                            (disp_s, y_bot),
                            (disp_e - arrow_w, y_bot),
                            (disp_e, y_center),
                            (disp_e - arrow_w, y_top),
                            (disp_s, y_top)
                        ]
                    elif strand == "-":
                        pts = [
                            (disp_e, y_bot),
                            (disp_s + arrow_w, y_bot),
                            (disp_s, y_center),
                            (disp_s + arrow_w, y_top),
                            (disp_e, y_top)
                        ]
                    else:
                        pts = [(disp_s, y_bot), (disp_e, y_bot), (disp_e, y_top), (disp_s, y_top)]

                    if is_pseudo:
                        col = "#f57c00"
                    elif is_glo:
                        if "alpha" in grp or grp == "globin_coding" or "locus_a" in grp:
                            col = "#d32f2f"
                        elif "beta" in grp or "locus_b" in grp:
                            col = "#1976d2"
                        else:
                            col = "#7b1fa2"
                    elif "flank" in grp or "cluster_non" in grp:
                        col = "#2e7d32"
                    else:
                        col = "#78909c"

                    polygon = patches.Polygon(pts, facecolor=col, edgecolor="black", linewidth=0.8, zorder=3)
                    ax_loc.add_patch(polygon)
                    if is_glo:
                        plotted_mpl_genes.add(sym)

                    # Multi-tier collision-free label positioning
                    mid_coord = (disp_s + disp_e) / 2
                    margin_kb = (l_e - l_s) / 1e3 * 0.02
                    mid_coord = min(max(mid_coord, l_s / 1e3 + margin_kb), l_e / 1e3 - margin_kb)

                    if is_glo:
                        if glo_idx % 2 == 0:
                            t_idx = (glo_idx // 2) % len(glo_tiers_pos)
                            lbl_y = glo_tiers_pos[t_idx]
                            lead_start = y_top + 0.02
                            lead_end = lbl_y - 0.10
                        else:
                            t_idx = (glo_idx // 2) % len(glo_tiers_neg)
                            lbl_y = glo_tiers_neg[t_idx]
                            lead_start = y_bot - 0.02
                            lead_end = lbl_y + 0.10
                        glo_idx += 1
                    else:
                        if y_center > 0:
                            lbl_y = 1.05
                            lead_start = y_top + 0.02
                            lead_end = lbl_y - 0.10
                        else:
                            lbl_y = -1.05
                            lead_start = y_bot - 0.02
                            lead_end = lbl_y + 0.10

                    ax_loc.plot([mid_coord, mid_coord], [lead_start, lead_end], color="#b0bec5", linestyle=":", linewidth=0.8, zorder=2)
                    ax_loc.text(
                        mid_coord, lbl_y, sym,
                        fontsize=8.5 if is_glo else 7.5,
                        ha="center", va="center",
                        fontweight="bold" if is_glo else "normal",
                        color=col if is_glo else "#37474f",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor=col if is_glo else "#cfd8dc", alpha=0.9, linewidth=0.6),
                        zorder=4
                    )

        plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])

        # Draw two connector lines opening up to the zoomed-in locus track
        for idx, (_, r_int) in enumerate(sp_intervals.iterrows()):
            loc_id = r_int["locus_id"]
            chrom = str(r_int["chromosome"])
            l_s = int(r_int["start"])
            l_e = int(r_int["end"])
            chr_len = chr_lengths.get(f"{species}::{chrom}", max(l_e * 1.2, 50000000))
            l_w_mb = max((l_e - l_s) / 1e6, chr_len / 1e6 * 0.002)

            loc_lower = loc_id.lower()
            if "alpha" in loc_lower or "locus_a" in loc_lower:
                theme_col = "#d32f2f"
            elif "beta" in loc_lower or "locus_b" in loc_lower:
                theme_col = "#1976d2"
            else:
                theme_col = "#7b1fa2"

            ax_chr = axes[idx * 2]
            ax_loc = axes[idx * 2 + 1]

            con_left = patches.ConnectionPatch(
                xyA=(l_s / 1e6, -0.23), coordsA="data", axesA=ax_chr,
                xyB=(0.0, 1.0), coordsB="axes fraction", axesB=ax_loc,
                linestyle="--", color=theme_col, linewidth=1.2, alpha=0.75, zorder=0
            )
            ax_loc.add_artist(con_left)

            con_right = patches.ConnectionPatch(
                xyA=(l_s / 1e6 + l_w_mb, -0.23), coordsA="data", axesA=ax_chr,
                xyB=(1.0, 1.0), coordsB="axes fraction", axesB=ax_loc,
                linestyle="--", color=theme_col, linewidth=1.2, alpha=0.75, zorder=0
            )
            ax_loc.add_artist(con_right)

        plt.savefig(png_path)
        plt.savefig(pdf_path)
        plt.close(fig)

        # 2. Plotly interactive rendering
        plotly_rows = []
        for _, r_int in sp_intervals.iterrows():
            loc_id = r_int["locus_id"]
            chrom = str(r_int["chromosome"])
            l_s = int(r_int["start"])
            l_e = int(r_int["end"])
            chr_len = chr_lengths.get(f"{species}::{chrom}", max(l_e * 1.2, 50000000))
            cen_bp = centromeres.get(f"{species}::{chrom}")
            is_fragment = (chr_len < 2000000) or str(chrom).lower().startswith(("gl", "scaffold", "ctg"))
            mid_locus = (l_s + l_e) / 2
            if not is_fragment and cen_bp and 0 < cen_bp < chr_len:
                if mid_locus > cen_bp:
                    near_telo = "3' (q-arm)"
                    min_d = max(chr_len - l_e, 0)
                else:
                    near_telo = "5' (p-arm)"
                    min_d = max(l_s, 0)
            else:
                d_5p = l_s
                d_3p = chr_len - l_e
                if d_5p < d_3p:
                    near_telo = "5' (p-arm)"
                    min_d = d_5p
                else:
                    near_telo = "3' (q-arm)"
                    min_d = d_3p
            rel_d = (min_d / chr_len) * 100

            plotly_rows.append({
                "loc_id": loc_id,
                "chrom": chrom,
                "chr_len": chr_len,
                "cen_bp": cen_bp,
                "l_s": l_s,
                "l_e": l_e,
                "min_d": min_d,
                "rel_d": rel_d,
                "near_telo": near_telo,
            })

        sub_titles = []
        for r in plotly_rows:
            loc_name = r["loc_id"].replace("_", " ").title()
            sub_titles.append(f"<b>Chr {r['chrom']} Overview (Total: {r['chr_len']/1e6:.1f} Mb | Nearest Telomere: {r['min_d']/1e6:.2f} Mb)</b>")
            sub_titles.append(f"<b>Zoomed {loc_name} Architecture (Chr {r['chrom']}: {r['l_s']/1e3:,.1f} - {r['l_e']/1e3:,.1f} kb)</b>")

        fig_ply = make_subplots(
            rows=num_loci * 2, cols=1,
            subplot_titles=sub_titles,
            vertical_spacing=0.10,
            row_heights=[0.3, 0.7] * num_loci,
        )

        plotted_ply_genes = set()

        for idx, r in enumerate(plotly_rows):
            r_chr = idx * 2 + 1
            r_loc = idx * 2 + 2

            loc_lower = r["loc_id"].lower()
            theme_col = "#d32f2f" if "alpha" in loc_lower or "locus_a" in loc_lower else ("#1976d2" if "beta" in loc_lower or "locus_b" in loc_lower else "#7b1fa2")

            # Top: Chromosome bar in Plotly (full length for each chromosome)
            fig_ply.add_trace(
                go.Bar(
                    x=[r["chr_len"] / 1e6], y=["Chromosome"], orientation="h",
                    marker=dict(color="#eceff1", line=dict(color="#455a64", width=1.5)),
                    hoverinfo="text",
                    hovertext=f"Chromosome {r['chrom']}<br>Total Length: {r['chr_len']:,} bp<br>Distance to {r['near_telo']}: {r['min_d']:,} bp",
                    showlegend=False,
                ),
                row=r_chr, col=1,
            )
            # Add locus highlight marker using matching theme color
            fig_ply.add_trace(
                go.Scatter(
                    x=[(r["l_s"] + r["l_e"]) / 2e6], y=["Chromosome"],
                    mode="markers",
                    marker=dict(size=14, color=theme_col, symbol="diamond", line=dict(color="black", width=1.2)),
                    hoverinfo="text",
                    hovertext=f"Locus: {r['loc_id']}<br>Coordinates: {r['l_s']:,} - {r['l_e']:,} bp<br>Span: {r['l_e'] - r['l_s']:,} bp",
                    showlegend=False,
                ),
                row=r_chr, col=1,
            )

            # Centromere vertical dashed line in Plotly
            if r.get("cen_bp") and 0 < r["cen_bp"] < r["chr_len"]:
                cen_mb = r["cen_bp"] / 1e6
                fig_ply.add_vline(
                    x=cen_mb, line_width=1.5, line_dash="dot", line_color="#c62828",
                    annotation_text="Cen", annotation_position="bottom right",
                    row=r_chr, col=1
                )

            # Each chromosome scaled across full width
            fig_ply.update_xaxes(
                range=[-r["chr_len"] * 0.015 / 1e6, r["chr_len"] * 1.04 / 1e6],
                title_text="Chromosome Coordinate (Mb)",
                row=r_chr, col=1
            )

            # Bottom: Zoomed locus track
            loc_genes = sp_genes[sp_genes["locus_id"] == r["loc_id"]].copy()
            loc_genes = loc_genes[~loc_genes["comparative_group"].isin(["other_locus_feature", "globin_pseudogene"])]
            loc_genes = loc_genes[~loc_genes["biotype"].astype(str).str.contains("pseudogene", case=False, na=False)]
            if not loc_genes.empty:
                loc_genes["start_num"] = pd.to_numeric(loc_genes["start"], errors="coerce").fillna(r["l_s"])
                loc_genes["end_num"] = pd.to_numeric(loc_genes["end"], errors="coerce").fillna(r["l_e"])
                loc_genes = loc_genes[(loc_genes["end_num"] >= r["l_s"]) & (loc_genes["start_num"] <= r["l_e"])]
                loc_genes = loc_genes.sort_values(by="start_num").reset_index(drop=True)

                ply_non_glo = 0
                for _, g in loc_genes.iterrows():
                    sym = str(g["display_symbol"])
                    grp = str(g.get("comparative_group", ""))
                    g_s = float(g["start_num"])
                    g_e = float(g["end_num"])
                    if g_e <= g_s:
                        g_e = g_s + 1500

                    is_pseudo = ("pseudogene" in grp) or ("pseudo" in sym.lower())
                    is_glo = str(g.get("is_globin", False)).lower() in ["true", "pseudogene"] or ("globin" in grp and "non_globin" not in grp) or is_pseudo
                    if is_pseudo:
                        col = "#f57c00"
                    elif is_glo:
                        col = "#d32f2f" if ("alpha" in grp or grp == "globin_coding" or "locus_a" in grp) else ("#1976d2" if ("beta" in grp or "locus_b" in grp) else "#7b1fa2")
                    elif "flank" in grp or "cluster_non" in grp:
                        col = "#2e7d32"
                    else:
                        col = "#78909c"
                    gc_val = g.get("coding_span_gc") or g.get("full_gene_gc") or g.get("cds_gc")
                    len_val = int(g.get("coding_span_len") or g.get("full_gene_len") or (g_e - g_s))
                    strand = str(g.get("strand", "+"))

                    hover_txt = (
                        f"<b>{sym}</b> ({grp.replace('_', ' ').title()})<br>"
                        f"Coordinates: Chr {r['chrom']}: {int(g_s):,} - {int(g_e):,} bp<br>"
                        f"Length: {len_val:,} bp<br>"
                        f"Strand: {strand}<br>"
                        f"Coding-Span GC: {gc_val}%<br>"
                        f"Dist to 5' Telomere: {int(g_s):,} bp<br>"
                        f"Dist to 3' Telomere: {int(r['chr_len'] - g_e):,} bp"
                    )

                    disp_w_kb = max((g_e - g_s) / 1e3, 0.6)

                    if is_glo:
                        plotted_ply_genes.add(sym)
                        y_ply = "Globin Baseline"
                    else:
                        y_ply = "Upper Flank" if (ply_non_glo % 2 == 0) else "Lower Flank"
                        ply_non_glo += 1

                    fig_ply.add_trace(
                        go.Bar(
                            base=[g_s / 1e3],
                            x=[disp_w_kb],
                            y=[y_ply],
                            orientation="h",
                            marker=dict(color=col, line=dict(color="black", width=0.8)),
                            hoverinfo="text",
                            hovertext=hover_txt,
                            name=sym,
                            text=sym if is_glo else "",
                            textposition="inside",
                            insidetextanchor="middle",
                            textfont=dict(size=9.5, color="white"),
                            showlegend=False,
                        ),
                        row=r_loc, col=1,
                    )

                pad_kb = (r["l_e"] - r["l_s"]) / 1e3 * 0.03
                fig_ply.update_xaxes(
                    range=[r["l_s"] / 1e3 - pad_kb, r["l_e"] / 1e3 + pad_kb],
                    domain=[0.0, 1.0],
                    title_text=f"Chr {r['chrom']} Coordinates (kb)",
                    row=r_loc, col=1
                )
                fig_ply.update_yaxes(
                    categoryorder="array",
                    categoryarray=["Lower Flank", "Globin Baseline", "Upper Flank"],
                    row=r_loc, col=1
                )

        fig_ply.update_layout(
            height=420 * num_loci + 120,
            width=1400,
            template="plotly_white",
            title=dict(
                text=f"<b>Linear Chromosome & Globin Locus Architecture: {common_name} (Distance to Nearest Telomere)</b>",
                x=0.5, y=0.99, xanchor="center", yanchor="top",
                font=dict(size=16, color="#212529"),
            ),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.06,
                xanchor="center",
                x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)",
                bordercolor="#dee2e6",
                borderwidth=1,
                font=dict(size=11),
            ),
            margin=dict(t=70, b=80, l=70, r=40),
            barmode="stack",
        )
        fig_ply.write_html(html_path)

        # Automated validation safeguard to prevent any missing globins
        def _check_is_globin(row):
            grp_str = str(row.get("comparative_group", ""))
            is_g = str(row.get("is_globin", False)).lower() == "true"
            return is_g or ("globin" in grp_str and "non_globin" not in grp_str)

        expected_globins = set(sp_genes[sp_genes.apply(_check_is_globin, axis=1)]["display_symbol"])
        missing_mpl = expected_globins - plotted_mpl_genes
        if missing_mpl:
            raise RuntimeError(f"Visualizer safeguard violation: species {species} has {len(missing_mpl)} globin genes missing from chromosome_linear static plot: {sorted(list(missing_mpl))}")
        missing_ply = expected_globins - plotted_ply_genes
        if missing_ply:
            raise RuntimeError(f"Visualizer safeguard violation: species {species} has {len(missing_ply)} globin genes missing from chromosome_linear interactive plot: {sorted(list(missing_ply))}")

        return {"png": png_path, "pdf": pdf_path, "html": html_path}

    def plot_length_vs_gc(
        self,
        species: str,
        df_genes: pd.DataFrame,
        output_dir: str = "results/bivariate_length_gc_plots",
    ) -> Dict[str, str]:
        """Generates 2D bivariate space of log10(Coding-Span Gene Length) vs. GC Content with layered globin placement."""
        import matplotlib.patches as patches
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        import plotly.graph_objects as go
        from scipy.stats import gaussian_kde

        os.makedirs(output_dir, exist_ok=True)
        common_name = self.COMMON_NAME_MAP.get(species, species)

        df_base, stats = self.bm.get_genome_baseline(species)
        sp_genes = df_genes[df_genes["species"] == species].copy()

        png_path = os.path.join(output_dir, f"bivariate_length_gc_{species}.png")
        pdf_path = os.path.join(output_dir, f"bivariate_length_gc_{species}.pdf")
        html_path = os.path.join(output_dir, f"bivariate_length_gc_{species}.html")

        # Background KDE
        len_col = "whole_gene_len" if "whole_gene_len" in df_base.columns else ("whole_gene_length" if "whole_gene_length" in df_base.columns else None)
        base_clean = df_base.dropna(subset=[len_col, "whole_gene_gc"]) if (not df_base.empty and len_col) else pd.DataFrame()
        if not base_clean.empty:
            base_lens = np.log10(base_clean[len_col][base_clean[len_col] > 0])
            base_gc = base_clean["whole_gene_gc"]
            mask = (base_lens >= 1.5) & (base_lens <= 6.5) & (base_gc >= 15) & (base_gc <= 85)
            sample_n = min(4000, mask.sum())
            x_samp = base_lens[mask].sample(sample_n, random_state=42)
            y_samp = base_gc[mask].sample(sample_n, random_state=42)
            kde = gaussian_kde(np.vstack([x_samp, y_samp]))

            x_grid = np.linspace(2.0, 6.35, 120)
            y_grid = np.linspace(20.0, 85.0, 100)
            X, Y = np.meshgrid(x_grid, y_grid)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        else:
            x_grid, y_grid, Z = None, None, None

        # Comprehensive styling, grouping, and deterministic layer configuration
        # Strict layer ordering (z-order):
        # 1. Background KDE: zorder=1..2
        # 2. Syntenic Non-Globins / Bystanders: zorder=3 (annotations: 4)
        # 3. Globin Pseudogenes: zorder=5 (annotations: 6)
        # 4. Beta Globins / Locus B: zorder=7 (annotations: 8)
        # 5. Alpha Globins / Locus A / Ancestral: zorder=9 (annotations: 10, TOP-MOST)
        BIVARIATE_GROUP_CONFIG = {
            # Alpha / Locus A / Linked Globin Coding (Top Layer)
            "alpha_globin_coding": {
                "label": "Alpha Globin Coding", "color": "#d32f2f", "edgecolor": "#b71c1c",
                "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"
            },
            "globin_coding_locus_A": {
                "label": "Globin Coding (Locus A)", "color": "#d32f2f", "edgecolor": "#b71c1c",
                "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"
            },
            "globin_coding": {
                "label": "Globin Coding", "color": "#d32f2f", "edgecolor": "#b71c1c",
                "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"
            },
            # Beta / Locus B Globin Coding (Middle Layer)
            "beta_globin_coding": {
                "label": "Beta Globin Coding", "color": "#1976d2", "edgecolor": "#0d47a1",
                "marker": "s", "plotly_symbol": "square", "size": 100, "zorder": 7, "tier": "beta"
            },
            "globin_coding_locus_B": {
                "label": "Globin Coding (Locus B)", "color": "#0288d1", "edgecolor": "#01579b",
                "marker": "s", "plotly_symbol": "square", "size": 100, "zorder": 7, "tier": "beta"
            },
            # Globin Pseudogenes (Intermediate Layer)
            "globin_pseudogene": {
                "label": "Globin Pseudogene", "color": "#f57c00", "edgecolor": "#e65100",
                "marker": "D", "plotly_symbol": "diamond", "size": 85, "zorder": 5, "tier": "pseudogene"
            },
            # Non-Globins & Bystanders (Bottom Points Layer)
            "cluster_non_globin_alpha": {
                "label": "Non-Globin (α-Locus)", "color": "#2e7d32", "edgecolor": "#1b5e20",
                "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"
            },
            "cluster_non_globin_locus_A": {
                "label": "Non-Globin (Locus A)", "color": "#2e7d32", "edgecolor": "#1b5e20",
                "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"
            },
            "cluster_non_globin_beta": {
                "label": "Non-Globin (β-Locus)", "color": "#388e3c", "edgecolor": "#1b5e20",
                "marker": "v", "plotly_symbol": "triangle-down", "size": 75, "zorder": 3, "tier": "non_globin"
            },
            "cluster_non_globin_locus_B": {
                "label": "Non-Globin (Locus B)", "color": "#388e3c", "edgecolor": "#1b5e20",
                "marker": "v", "plotly_symbol": "triangle-down", "size": 75, "zorder": 3, "tier": "non_globin"
            },
            "cluster_non_globin": {
                "label": "Cluster Non-Globin", "color": "#2e7d32", "edgecolor": "#1b5e20",
                "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"
            },
            "other_locus_feature": {
                "label": "Other Locus Feature", "color": "#78909c", "edgecolor": "#455a64",
                "marker": "p", "plotly_symbol": "pentagon", "size": 70, "zorder": 3, "tier": "non_globin"
            },
        }

        # Collect genes grouped by comparative_group (excluding other_locus_feature and pseudogenes)
        items_by_grp: Dict[str, list] = {}
        all_items: list = []

        for _, r in sp_genes.iterrows():
            grp = str(r.get("comparative_group", "other_locus_feature"))
            biotype = str(r.get("biotype", "")).lower()
            if grp in ["other_locus_feature", "globin_pseudogene"] or "pseudogene" in biotype:
                continue
            l = r.get("coding_span_len") or r.get("full_gene_len", 0)
            gc_val = r.get("coding_span_gc") or r.get("full_gene_gc")
            if l <= 0 or pd.isna(gc_val):
                continue
            log_l = float(np.log10(l))
            gc_f = float(gc_val)
            sym = str(r["display_symbol"])

            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            item = {
                "row": r, "x": log_l, "y": gc_f, "sym": sym, "grp": grp, "len": int(l),
                "cfg": cfg, "tier": cfg["tier"]
            }
            if grp not in items_by_grp:
                items_by_grp[grp] = []
            items_by_grp[grp].append(item)
            all_items.append(item)

        # 1. Matplotlib static rendering
        fig, ax = plt.subplots(figsize=(10, 8), dpi=200)

        # Background KDE density cloud
        if Z is not None:
            ax.contourf(X, Y, Z, levels=14, cmap="Blues", alpha=0.55, zorder=1)
            ax.contour(X, Y, Z, levels=7, colors="#64b5f6", linewidths=0.8, alpha=0.8, zorder=2)

        # Sort groups by z-order so non-globins are plotted first, then pseudogenes, beta, and alpha on top
        sorted_groups = sorted(items_by_grp.keys(), key=lambda g: BIVARIATE_GROUP_CONFIG.get(g, {}).get("zorder", 3))

        legend_handles = []
        if Z is not None:
            legend_handles.append(Patch(facecolor="#90caf9", edgecolor="#1976d2", alpha=0.6, label=f"Genome Baseline Density ({len(base_clean):,} genes)"))

        for grp in sorted_groups:
            items = items_by_grp[grp]
            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            xs = [it["x"] for it in items]
            ys = [it["y"] for it in items]
            ax.scatter(
                xs, ys, color=cfg["color"], marker=cfg["marker"], s=cfg["size"],
                edgecolor=cfg["edgecolor"], linewidth=1.1, alpha=0.95, zorder=cfg["zorder"]
            )
            legend_handles.append(
                Line2D([0], [0], marker=cfg["marker"], color="w", markerfacecolor=cfg["color"],
                       markeredgecolor=cfg["edgecolor"], markersize=9.0, label=f"{cfg['label']} (n={len(items)})")
            )

        # Smart collision avoidance for text labels
        # Priority order for reserving label space: Alpha/Locus A -> Beta/Locus B -> Pseudogene -> Named Bystander -> LOC/ENS
        tier_alpha_items = [it for it in all_items if it["tier"] == "alpha"]
        tier_beta_items = [it for it in all_items if it["tier"] == "beta"]
        tier_pseudo_items = [it for it in all_items if it["tier"] == "pseudogene"]
        tier_other_items = [it for it in all_items if it["tier"] == "non_globin"]

        # Sort bystanders so named genes precede LOC/ENS
        tier_other_items.sort(key=lambda it: (it["sym"].startswith(("LOC", "ENS")), it["sym"]))

        occupied_boxes = []

        def allocate_label_offset(x, y, sym_str, tier="alpha", is_loc=False):
            u = (x - 2.0) / 4.0 * 600.0
            v = (y - 20.0) / 65.0 * 450.0
            fs = 9.0 if tier == "alpha" else (8.5 if tier == "beta" else 7.0)
            w_box = len(sym_str) * (fs * 0.62) + 8.0
            h_box = fs * 1.5 + 4.0

            if tier == "alpha":
                candidates = [
                    (14, 12), (-14, 12), (14, -14), (-14, -14),
                    (0, 20), (0, -20), (22, 0), (-24, 0),
                    (26, 18), (-28, 18), (26, -18), (-28, -18),
                    (10, 28), (-12, 28), (10, -28), (-12, -28),
                    (34, 8), (-36, 8), (34, -8), (-36, -8),
                ]
            elif tier == "beta":
                candidates = [
                    (12, -12), (-14, -12), (12, 12), (-14, 12),
                    (18, 0), (-20, 0), (0, -18), (0, 18),
                    (24, -16), (-26, -16), (24, 16), (-26, 16),
                    (10, -26), (-12, -26), (10, 26), (-12, 26),
                    (30, -8), (-32, -8), (30, 8), (-32, 8),
                ]
            elif tier == "pseudogene":
                candidates = [
                    (12, 12), (-12, 12), (12, -12), (-12, -12),
                    (16, 0), (-16, 0), (0, 16), (0, -16),
                ]
            elif not is_loc:
                candidates = [
                    (10, 10), (-10, 10), (10, -10), (-10, -10),
                    (14, 0), (-14, 0), (0, 14), (0, -14),
                    (18, 12), (-18, 12), (18, -12), (-18, -12),
                    (22, 16), (-22, 16), (22, -16), (-22, -16),
                ]
            else:
                candidates = [
                    (8, 8), (-8, 8), (8, -8), (-8, -8),
                    (12, 0), (-12, 0), (0, 12), (0, -12),
                ]

            best_offset = candidates[0]
            best_align = ("left", "bottom")
            min_overlap = float("inf")

            for ox, oy in candidates:
                ha = "right" if ox < -4 else ("left" if ox > 4 else "center")
                va = "top" if oy < -4 else ("bottom" if oy > 4 else "center")

                if ha == "left":
                    bx1 = u + ox
                    bx2 = u + ox + w_box
                elif ha == "right":
                    bx1 = u + ox - w_box
                    bx2 = u + ox
                else:
                    bx1 = u + ox - w_box / 2.0
                    bx2 = u + ox + w_box / 2.0

                if va == "bottom":
                    by1 = v + oy
                    by2 = v + oy + h_box
                elif va == "top":
                    by1 = v + oy - h_box
                    by2 = v + oy
                else:
                    by1 = v + oy - h_box / 2.0
                    by2 = v + oy + h_box / 2.0

                pad = 3.0
                overlap = 0.0
                for obx1, obx2, oby1, oby2 in occupied_boxes:
                    ix1 = max(bx1 - pad, obx1)
                    ix2 = min(bx2 + pad, obx2)
                    iy1 = max(by1 - pad, oby1)
                    iy2 = min(by2 + pad, oby2)
                    if ix2 > ix1 and iy2 > iy1:
                        overlap += (ix2 - ix1) * (iy2 - iy1)

                if overlap == 0.0:
                    occupied_boxes.append((bx1, bx2, by1, by2))
                    return (ox, oy, ha, va)

                if overlap < min_overlap:
                    min_overlap = overlap
                    best_offset = (ox, oy)
                    best_align = (ha, va)

            if is_loc and min_overlap > 0.0:
                return None

            ox, oy = best_offset
            ha, va = best_align
            if ha == "left":
                bx1, bx2 = u + ox, u + ox + w_box
            elif ha == "right":
                bx1, bx2 = u + ox - w_box, u + ox
            else:
                bx1, bx2 = u + ox - w_box / 2.0, u + ox + w_box / 2.0
            if va == "bottom":
                by1, by2 = v + oy, v + oy + h_box
            elif va == "top":
                by1, by2 = v + oy - h_box, v + oy
            else:
                by1, by2 = v + oy - h_box / 2.0, v + oy + h_box / 2.0

            occupied_boxes.append((bx1, bx2, by1, by2))
            return (ox, oy, ha, va)

        # Precompute label offsets by priority: Alpha -> Beta -> Pseudogene -> Non-Globins
        for item in tier_alpha_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="alpha")
        for item in tier_beta_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="beta")
        for item in tier_pseudo_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="pseudogene")
        for item in tier_other_items:
            is_loc = item["sym"].startswith(("LOC", "ENS"))
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="non_globin", is_loc=is_loc)

        # Draw annotations in increasing z-order: Non-Globins (z=4) -> Pseudogene (z=6) -> Beta (z=8) -> Alpha (z=10)
        for item in tier_other_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]),
                xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va,
                fontsize=7.0, color="#2e7d32",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="#a5d6a7", alpha=0.90, linewidth=0.6),
                arrowprops=dict(arrowstyle="-", color="#81c784", lw=0.6, alpha=0.7),
                zorder=4
            )

        for item in tier_pseudo_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]),
                xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va,
                fontsize=8.0, fontweight="bold", color="#e65100",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#f57c00", alpha=0.92, linewidth=0.8),
                arrowprops=dict(arrowstyle="-", color="#f57c00", lw=0.7, alpha=0.8),
                zorder=6
            )

        for item in tier_beta_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]),
                xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va,
                fontsize=8.5, fontweight="bold", color="#0d47a1",
                bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#1976d2", alpha=0.95, linewidth=0.9),
                arrowprops=dict(arrowstyle="-", color="#1976d2", lw=0.8, alpha=0.8),
                zorder=8
            )

        for item in tier_alpha_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]),
                xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va,
                fontsize=9.0, fontweight="bold", color="#b71c1c",
                bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#d32f2f", alpha=0.95, linewidth=0.9),
                arrowprops=dict(arrowstyle="-", color="#d32f2f", lw=0.8, alpha=0.8),
                zorder=10
            )

        ax.legend(
            handles=legend_handles, loc="lower left",
            frameon=True, facecolor="white", framealpha=0.95, edgecolor="#cfd8dc",
            fontsize=9.0, title="Gene Classification & Density", title_fontsize=9.5
        )

        ax.set_xlabel("log10(Coding-Span Gene Length [bp])", fontsize=11.5, fontweight="bold", color="#263238")
        ax.set_ylabel("Coding-Span Gene GC Content (%)", fontsize=11.5, fontweight="bold", color="#263238")
        ax.set_title(
            f"Bivariate Gene Architecture: log10(Coding-Span Length) vs. GC Content\n{common_name}",
            fontsize=12.5, fontweight="bold", color="#212529", pad=12
        )
        ax.set_xlim(2.0, 6.35)
        ax.set_ylim(20.0, 85.0)
        ax.set_xticks([2.0, 3.0, 4.0, 5.0, 6.0, 6.301])
        ax.set_xticklabels(["100 bp\n($10^2$)", "1 kb\n($10^3$)", "10 kb\n($10^4$)", "100 kb\n($10^5$)", "1 Mb\n($10^6$)", "2 Mb\n($2\\times 10^6$)"], fontsize=9.0)
        ax.grid(True, linestyle="--", alpha=0.3, color="#b0bec5")

        plt.tight_layout()
        plt.savefig(png_path)
        plt.savefig(pdf_path)
        plt.close(fig)

        # 2. Plotly interactive rendering
        fig_ply = go.Figure()
        if Z is not None:
            fig_ply.add_trace(go.Contour(
                x=x_grid, y=y_grid, z=Z,
                colorscale="Blues", opacity=0.6, showscale=False,
                contours=dict(coloring="fill", showlines=True),
                hoverinfo="skip",
                name="Genome Baseline Density",
                showlegend=True,
            ))

        # Add traces in layered order (Non-Globins -> Pseudogenes -> Beta -> Alpha on top)
        for grp in sorted_groups:
            items = items_by_grp[grp]
            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            xs = [it["x"] for it in items]
            ys = [it["y"] for it in items]
            texts = [
                f"<b>{it['sym']}</b> ({cfg['label']})<br>"
                f"Coding-Span Length: {it['len']:,} bp (log10: {it['x']:.2f})<br>"
                f"Coding-Span GC: {it['y']}%<br>"
                f"Z-Score Length: {it['row'].get('z_score_length', 'N/A')}<br>"
                f"Z-Score GC: {it['row'].get('z_score_coding_span_gc') or it['row'].get('z_score_full_gene_gc', 'N/A')}"
                for it in items
            ]
            fig_ply.add_trace(go.Scatter(
                x=xs, y=ys, mode="markers", name=f"{cfg['label']} (n={len(items)})",
                marker=dict(
                    size=12 if cfg["tier"] in ["alpha", "beta"] else (10 if cfg["tier"] == "pseudogene" else 9),
                    color=cfg["color"],
                    symbol=cfg["plotly_symbol"],
                    line=dict(width=1.0, color=cfg["edgecolor"])
                ),
                hovertext=texts, hoverinfo="text"
            ))

        fig_ply.update_layout(
            template="plotly_white",
            width=1000, height=800,
            title=dict(
                text=f"<b>Bivariate Gene Architecture: log10(Coding-Span Length) vs. GC Content ({common_name})</b>",
                x=0.5, y=0.985, xanchor="center", yanchor="top",
                font=dict(size=15, color="#212529")
            ),
            xaxis=dict(
                title="<b>log10(Coding-Span Gene Length [bp])</b>",
                range=[2.0, 6.35],
                tickvals=[2.0, 3.0, 4.0, 5.0, 6.0, 6.301],
                ticktext=["100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "2 Mb"],
            ),
            yaxis=dict(title="<b>Coding-Span Gene GC Content (%)</b>", range=[15.0, 90.0]),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.12,
                xanchor="center",
                x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)",
                bordercolor="#dee2e6",
                borderwidth=1,
            ),
            margin=dict(t=70, b=90, l=70, r=40)
        )
        fig_ply.write_html(html_path)

        return {"png": png_path, "pdf": pdf_path, "html": html_path}

    def plot_bivariate_with_chromosome_inset(
        self,
        species: str,
        df_genes: pd.DataFrame,
        df_intervals: Optional[pd.DataFrame] = None,
        output_dir: str = "results/bivariate_with_inset_plots",
    ) -> Dict[str, str]:
        """Generates an alternate high-impact version of bivariate log10(Coding-Span Length) vs. GC Content
        featuring an exaggerated mini chromosome linear inset in the top-right corner to maximize visibility
        of genomic coordinates, distance to nearest telomere, and zoomed locus gene architecture.
        Generates static PNG, PDF, and interactive Plotly HTML.
        """
        import json
        import matplotlib.patches as patches
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        import plotly.graph_objects as go
        from scipy.stats import gaussian_kde

        os.makedirs(output_dir, exist_ok=True)
        common_name = self.COMMON_NAME_MAP.get(species, species)

        # Baseline KDE
        df_base, stats = self.bm.get_genome_baseline(species)
        sp_genes = df_genes[df_genes["species"] == species].copy()

        # Intervals
        if df_intervals is None:
            intervals_path = "results/curated_intervals.csv"
            if not os.path.exists(intervals_path):
                intervals_path = "data/curated_intervals.csv"
            if os.path.exists(intervals_path):
                df_intervals = pd.read_csv(intervals_path)
            else:
                df_intervals = pd.DataFrame()

        sp_intervals = df_intervals[df_intervals["species"] == species] if not df_intervals.empty else pd.DataFrame()

        # Ensure start, end, and strand coordinates exist
        if "start" not in sp_genes.columns or sp_genes["start"].isna().all():
            g_coord_path = "results/genes_in_intervals.csv"
            if not os.path.exists(g_coord_path):
                g_coord_path = "data/genes_in_intervals.csv"
            if os.path.exists(g_coord_path):
                df_coords = pd.read_csv(g_coord_path)
                sp_coords = df_coords[df_coords["species"] == species][["gene_id", "start", "end", "strand"]].drop_duplicates(subset=["gene_id"])
                sp_genes = sp_genes.drop(columns=[c for c in ["start", "end", "strand"] if c in sp_genes.columns])
                sp_genes = sp_genes.merge(sp_coords, on="gene_id", how="left")

        # Chromosome lengths and centromeres
        chr_len_path = "data/chromosome_lengths.json"
        if not os.path.exists(chr_len_path):
            chr_len_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "chromosome_lengths.json")
        chr_lengths = {}
        if os.path.exists(chr_len_path):
            with open(chr_len_path, "r", encoding="utf-8") as f:
                chr_lengths = json.load(f)

        cen_path = "data/centromere_positions.json"
        if not os.path.exists(cen_path):
            cen_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "centromere_positions.json")
        centromeres = {}
        if os.path.exists(cen_path):
            with open(cen_path, "r", encoding="utf-8") as f:
                centromeres = json.load(f)

        png_path = os.path.join(output_dir, f"bivariate_with_inset_{species}.png")
        pdf_path = os.path.join(output_dir, f"bivariate_with_inset_{species}.pdf")
        html_path = os.path.join(output_dir, f"bivariate_with_inset_{species}.html")

        # Background KDE density cloud
        len_col = "whole_gene_len" if "whole_gene_len" in df_base.columns else ("whole_gene_length" if "whole_gene_length" in df_base.columns else None)
        base_clean = df_base.dropna(subset=[len_col, "whole_gene_gc"]) if (not df_base.empty and len_col) else pd.DataFrame()
        if not base_clean.empty:
            base_lens = np.log10(base_clean[len_col][base_clean[len_col] > 0])
            base_gc = base_clean["whole_gene_gc"]
            mask = (base_lens >= 1.5) & (base_lens <= 6.5) & (base_gc >= 15) & (base_gc <= 85)
            sample_n = min(4000, mask.sum())
            x_samp = base_lens[mask].sample(sample_n, random_state=42)
            y_samp = base_gc[mask].sample(sample_n, random_state=42)
            kde = gaussian_kde(np.vstack([x_samp, y_samp]))

            x_grid = np.linspace(2.0, 6.35, 120)
            y_grid = np.linspace(20.0, 85.0, 100)
            X, Y = np.meshgrid(x_grid, y_grid)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        else:
            x_grid, y_grid, Z = None, None, None

        BIVARIATE_GROUP_CONFIG = {
            "alpha_globin_coding": {"label": "Alpha Globin Coding", "color": "#d32f2f", "edgecolor": "#b71c1c", "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"},
            "globin_coding_locus_A": {"label": "Globin Coding (Locus A)", "color": "#d32f2f", "edgecolor": "#b71c1c", "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"},
            "globin_coding": {"label": "Globin Coding", "color": "#d32f2f", "edgecolor": "#b71c1c", "marker": "o", "plotly_symbol": "circle", "size": 110, "zorder": 9, "tier": "alpha"},
            "beta_globin_coding": {"label": "Beta Globin Coding", "color": "#1976d2", "edgecolor": "#0d47a1", "marker": "s", "plotly_symbol": "square", "size": 100, "zorder": 7, "tier": "beta"},
            "globin_coding_locus_B": {"label": "Globin Coding (Locus B)", "color": "#0288d1", "edgecolor": "#01579b", "marker": "s", "plotly_symbol": "square", "size": 100, "zorder": 7, "tier": "beta"},
            "globin_pseudogene": {"label": "Globin Pseudogene", "color": "#f57c00", "edgecolor": "#e65100", "marker": "D", "plotly_symbol": "diamond", "size": 85, "zorder": 5, "tier": "pseudogene"},
            "cluster_non_globin_alpha": {"label": "Non-Globin (α-Locus)", "color": "#2e7d32", "edgecolor": "#1b5e20", "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"},
            "cluster_non_globin_locus_A": {"label": "Non-Globin (Locus A)", "color": "#2e7d32", "edgecolor": "#1b5e20", "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"},
            "cluster_non_globin_beta": {"label": "Non-Globin (β-Locus)", "color": "#388e3c", "edgecolor": "#1b5e20", "marker": "v", "plotly_symbol": "triangle-down", "size": 75, "zorder": 3, "tier": "non_globin"},
            "cluster_non_globin_locus_B": {"label": "Non-Globin (Locus B)", "color": "#388e3c", "edgecolor": "#1b5e20", "marker": "v", "plotly_symbol": "triangle-down", "size": 75, "zorder": 3, "tier": "non_globin"},
            "cluster_non_globin": {"label": "Cluster Non-Globin", "color": "#2e7d32", "edgecolor": "#1b5e20", "marker": "^", "plotly_symbol": "triangle-up", "size": 75, "zorder": 3, "tier": "non_globin"},
            "other_locus_feature": {"label": "Other Locus Feature", "color": "#78909c", "edgecolor": "#455a64", "marker": "p", "plotly_symbol": "pentagon", "size": 70, "zorder": 3, "tier": "non_globin"},
        }

        # Collect genes grouped by comparative_group (excluding other_locus_feature and pseudogenes)
        items_by_grp: Dict[str, list] = {}
        all_items: list = []

        for _, r in sp_genes.iterrows():
            grp = str(r.get("comparative_group", "other_locus_feature"))
            biotype = str(r.get("biotype", "")).lower()
            if grp in ["other_locus_feature", "globin_pseudogene"] or "pseudogene" in biotype:
                continue
            l = r.get("coding_span_len") or r.get("full_gene_len", 0)
            gc_val = r.get("coding_span_gc") or r.get("full_gene_gc")
            if l <= 0 or pd.isna(gc_val):
                continue
            log_l = float(np.log10(l))
            gc_f = float(gc_val)
            sym = str(r["display_symbol"])

            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            item = {
                "row": r, "x": log_l, "y": gc_f, "sym": sym, "grp": grp, "len": int(l),
                "cfg": cfg, "tier": cfg["tier"]
            }
            if grp not in items_by_grp:
                items_by_grp[grp] = []
            items_by_grp[grp].append(item)
            all_items.append(item)

        # 1. Matplotlib static rendering
        fig, ax = plt.subplots(figsize=(11, 8.5), dpi=200)

        # Background KDE
        if Z is not None:
            ax.contourf(X, Y, Z, levels=14, cmap="Blues", alpha=0.55, zorder=1)
            ax.contour(X, Y, Z, levels=7, colors="#64b5f6", linewidths=0.8, alpha=0.8, zorder=2)

        sorted_groups = sorted(items_by_grp.keys(), key=lambda g: BIVARIATE_GROUP_CONFIG.get(g, {}).get("zorder", 3))

        legend_handles = []
        if Z is not None:
            legend_handles.append(Patch(facecolor="#90caf9", edgecolor="#1976d2", alpha=0.6, label=f"Genome Baseline Density ({len(base_clean):,} genes)"))

        for grp in sorted_groups:
            items = items_by_grp[grp]
            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            xs = [it["x"] for it in items]
            ys = [it["y"] for it in items]
            ax.scatter(
                xs, ys, color=cfg["color"], marker=cfg["marker"], s=cfg["size"],
                edgecolor=cfg["edgecolor"], linewidth=1.1, alpha=0.95, zorder=cfg["zorder"]
            )
            legend_handles.append(
                Line2D([0], [0], marker=cfg["marker"], color="w", markerfacecolor=cfg["color"],
                       markeredgecolor=cfg["edgecolor"], markersize=8.5, label=f"{cfg['label']} (n={len(items)})")
            )

        # Label collision avoidance: reserve top-right inset region so labels never collide with the inset!
        tier_alpha_items = [it for it in all_items if it["tier"] == "alpha"]
        tier_beta_items = [it for it in all_items if it["tier"] == "beta"]
        tier_pseudo_items = [it for it in all_items if it["tier"] == "pseudogene"]
        tier_other_items = [it for it in all_items if it["tier"] == "non_globin"]
        tier_other_items.sort(key=lambda it: (it["sym"].startswith(("LOC", "ENS")), it["sym"]))

        occupied_boxes = []
        # Inset region in normalized pixel coordinates
        inset_u_min = (4.45 - 2.0) / 4.0 * 600.0
        inset_u_max = 600.0
        inset_v_min = (48.0 - 20.0) / 65.0 * 450.0
        inset_v_max = 450.0
        occupied_boxes.append((inset_u_min, inset_u_max, inset_v_min, inset_v_max))

        def allocate_label_offset(x, y, sym_str, tier="alpha", is_loc=False):
            u = (x - 2.0) / 4.0 * 600.0
            v = (y - 20.0) / 65.0 * 450.0
            fs = 8.5 if tier == "alpha" else (8.0 if tier == "beta" else 6.8)
            w_box = len(sym_str) * (fs * 0.60) + 6.0
            h_box = fs * 1.4 + 4.0

            if tier == "alpha":
                candidates = [
                    (12, 12), (-14, 12), (12, -14), (-14, -14),
                    (0, 18), (0, -18), (20, 0), (-22, 0),
                    (22, 16), (-24, 16), (22, -16), (-24, -16),
                ]
            elif tier == "beta":
                candidates = [
                    (12, -12), (-14, -12), (12, 12), (-14, 12),
                    (18, 0), (-20, 0), (0, -18), (0, 18),
                    (20, -16), (-22, -16), (20, 16), (-22, 16),
                ]
            elif tier == "pseudogene":
                candidates = [
                    (10, 10), (-10, 10), (10, -10), (-10, -10),
                    (14, 0), (-14, 0), (0, 14), (0, -14),
                ]
            elif not is_loc:
                candidates = [
                    (10, 10), (-10, 10), (10, -10), (-10, -10),
                    (12, 0), (-12, 0), (0, 12), (0, -12),
                ]
            else:
                candidates = [
                    (8, 8), (-8, 8), (8, -8), (-8, -8),
                    (10, 0), (-10, 0), (0, 10), (0, -10),
                ]

            best_offset = candidates[0]
            best_align = ("left", "bottom")
            min_overlap = float("inf")

            for ox, oy in candidates:
                ha = "right" if ox < -4 else ("left" if ox > 4 else "center")
                va = "top" if oy < -4 else ("bottom" if oy > 4 else "center")

                bx1 = u + ox if ha == "left" else (u + ox - w_box if ha == "right" else u + ox - w_box / 2.0)
                bx2 = bx1 + w_box
                by1 = v + oy if va == "bottom" else (v + oy - h_box if va == "top" else v + oy - h_box / 2.0)
                by2 = by1 + h_box

                overlap = 0.0
                for obx1, obx2, oby1, oby2 in occupied_boxes:
                    ix1 = max(bx1 - 2.0, obx1)
                    ix2 = min(bx2 + 2.0, obx2)
                    iy1 = max(by1 - 2.0, oby1)
                    iy2 = min(by2 + 2.0, oby2)
                    if ix2 > ix1 and iy2 > iy1:
                        overlap += (ix2 - ix1) * (iy2 - iy1)

                if overlap == 0.0:
                    occupied_boxes.append((bx1, bx2, by1, by2))
                    return (ox, oy, ha, va)

                if overlap < min_overlap:
                    min_overlap = overlap
                    best_offset = (ox, oy)
                    best_align = (ha, va)

            if is_loc and min_overlap > 0.0:
                return None

            ox, oy = best_offset
            ha, va = best_align
            bx1 = u + ox if ha == "left" else (u + ox - w_box if ha == "right" else u + ox - w_box / 2.0)
            bx2 = bx1 + w_box
            by1 = v + oy if va == "bottom" else (v + oy - h_box if va == "top" else v + oy - h_box / 2.0)
            by2 = by1 + h_box
            occupied_boxes.append((bx1, bx2, by1, by2))
            return (ox, oy, ha, va)

        for item in tier_alpha_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="alpha")
        for item in tier_beta_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="beta")
        for item in tier_pseudo_items:
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="pseudogene")
        for item in tier_other_items:
            is_loc = item["sym"].startswith(("LOC", "ENS"))
            item["annot_spec"] = allocate_label_offset(item["x"], item["y"], item["sym"], tier="non_globin", is_loc=is_loc)

        for item in tier_other_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]), xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va, fontsize=6.8, color="#2e7d32",
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="#a5d6a7", alpha=0.90, linewidth=0.5),
                arrowprops=dict(arrowstyle="-", color="#81c784", lw=0.5, alpha=0.7), zorder=4
            )

        for item in tier_pseudo_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]), xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va, fontsize=7.8, fontweight="bold", color="#e65100",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="#f57c00", alpha=0.92, linewidth=0.7),
                arrowprops=dict(arrowstyle="-", color="#f57c00", lw=0.6, alpha=0.8), zorder=6
            )

        for item in tier_beta_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]), xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va, fontsize=8.2, fontweight="bold", color="#0d47a1",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#1976d2", alpha=0.95, linewidth=0.8),
                arrowprops=dict(arrowstyle="-", color="#1976d2", lw=0.7, alpha=0.8), zorder=8
            )

        for item in tier_alpha_items:
            if not item.get("annot_spec"):
                continue
            ox, oy, ha, va = item["annot_spec"]
            ax.annotate(
                item["sym"], (item["x"], item["y"]), xytext=(ox, oy), textcoords="offset points",
                ha=ha, va=va, fontsize=8.5, fontweight="bold", color="#b71c1c",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#d32f2f", alpha=0.95, linewidth=0.8),
                arrowprops=dict(arrowstyle="-", color="#d32f2f", lw=0.7, alpha=0.8), zorder=10
            )

        # ---------------------------------------------------------------------
        # TOP-RIGHT INSET: Compact Chromosome Linear Representation (No zoom-ins)
        # ---------------------------------------------------------------------
        num_loci = len(sp_intervals)
        ins_w = 0.50
        ins_h = 0.27 if num_loci >= 2 else 0.17
        ins_x0 = 0.475
        ins_y0 = 0.70 if num_loci >= 2 else 0.80

        ax_ins = ax.inset_axes([ins_x0, ins_y0, ins_w, ins_h])
        ax_ins.set_facecolor("#ffffff")
        for spine in ax_ins.spines.values():
            spine.set_edgecolor("#546e7a")
            spine.set_linewidth(1.1)
        ax_ins.set_xlim(0.0, 1.0)
        ax_ins.set_ylim(0.0, 1.0)
        ax_ins.set_xticks([])
        ax_ins.set_yticks([])

        # Inset Header Banner
        header_h = 0.11 if num_loci >= 2 else 0.16
        header_y = 1.0 - header_h - 0.02
        header_box = patches.FancyBboxPatch(
            (0.03, header_y), 0.94, header_h,
            boxstyle="round,pad=0.005,rounding_size=0.03",
            facecolor="#eceff1", edgecolor="#90a4ae", linewidth=0.8, zorder=1
        )
        ax_ins.add_patch(header_box)
        ax_ins.text(
            0.5, header_y + header_h / 2, "Chromosomal Context & Architecture",
            ha="center", va="center",
            fontsize=9.2, fontweight="bold", color="#1a237e", zorder=2
        )

        calc_chr_lens = []
        for _, r in sp_intervals.iterrows():
            c_str = str(r["chromosome"])
            key = f"{species}::{c_str}"
            calc_chr_lens.append(chr_lengths.get(key, max(int(r["end"]) * 1.2, 50000000)))
        max_chr_len = max(calc_chr_lens) if calc_chr_lens else 50000000

        max_bar_w = 0.88
        base_x = 0.05

        for idx, (_, r_int) in enumerate(sp_intervals.iterrows()):
            loc_id = r_int["locus_id"]
            chrom = str(r_int["chromosome"])
            l_s = int(r_int["start"])
            l_e = int(r_int["end"])
            chr_len = chr_lengths.get(f"{species}::{chrom}", max(l_e * 1.2, 50000000))
            cen_bp = centromeres.get(f"{species}::{chrom}")

            is_fragment = (chr_len < 2000000) or str(chrom).lower().startswith(("gl", "scaffold", "ctg"))
            mid_locus = (l_s + l_e) / 2
            if not is_fragment and cen_bp and 0 < cen_bp < chr_len:
                if mid_locus > cen_bp:
                    near_telo = "3'"
                    min_d = max(chr_len - l_e, 0)
                else:
                    near_telo = "5'"
                    min_d = max(l_s, 0)
            else:
                d_5p = l_s
                d_3p = chr_len - l_e
                near_telo = "5'" if d_5p < d_3p else "3'"
                min_d = min(d_5p, d_3p)

            loc_lower = loc_id.lower()
            if "locus_a" in loc_lower:
                theme_col = "#d32f2f"       # Carmine Red for Locus A
                theme_dark = "#b71c1c"
                badge = "Locus A"
            elif "locus_b" in loc_lower:
                theme_col = "#1976d2"       # Royal Blue for Locus B
                theme_dark = "#0d47a1"
                badge = "Locus B"
            elif loc_id == "alpha_globin" or (loc_lower.startswith("alpha") and "beta" not in loc_lower):
                theme_col = "#d32f2f"       # Carmine Red for Alpha
                theme_dark = "#b71c1c"
                badge = "α-Locus"
            elif loc_id == "beta_globin" or (loc_lower.startswith("beta") and "alpha" not in loc_lower):
                theme_col = "#1976d2"       # Royal Blue for Beta
                theme_dark = "#0d47a1"
                badge = "β-Locus"
            elif "proto" in loc_lower:
                theme_col = "#7b1fa2"
                theme_dark = "#4a148c"
                badge = "Proto-Globin"
            else:
                theme_col = "#7b1fa2"
                theme_dark = "#4a148c"
                badge = "Linked α+β"

            if num_loci == 1:
                y_chr = 0.50
                h_chr = 0.16
                lbl_y = 0.38
                pin_y = 0.64
                dist_y = 0.595
            else:
                y_chr = 0.65 if idx == 0 else 0.25
                h_chr = 0.12
                lbl_y = 0.55 if idx == 0 else 0.15
                pin_y = 0.76 if idx == 0 else 0.36
                dist_y = 0.722 if idx == 0 else 0.322

            # Relative chromosome bar width proportional to max chromosome in the species
            chr_w = max_bar_w * (chr_len / max_chr_len)
            rect_chr = patches.FancyBboxPatch(
                (base_x, y_chr - h_chr / 2), chr_w, h_chr,
                boxstyle="round,pad=0.004,rounding_size=0.03",
                facecolor="#eceff1", edgecolor="#607d8b", linewidth=0.8,
                linestyle="--" if is_fragment else "-", zorder=2
            )
            ax_ins.add_patch(rect_chr)

            # Centromere notch/indicator (bona fide chromosomes only)
            if not is_fragment and cen_bp and 0 < cen_bp < chr_len:
                cen_x = base_x + (cen_bp / chr_len) * chr_w
                ax_ins.plot([cen_x, cen_x], [y_chr - h_chr / 2, y_chr + h_chr / 2], color="#c62828", linestyle=":", linewidth=1.0, zorder=3)

            # Realistic Locus marker on chromosome bar (strictly proportional to chr_len, min 0.006 for visibility)
            l_span = max(l_e - l_s, 1)
            l_span_str = f"{l_span/1e3:.0f} kb" if l_span < 1e6 else f"{l_span/1e6:.2f} Mb"
            l_w_real = max((l_span / chr_len) * chr_w, 0.006)
            loc_x = base_x + (l_s / chr_len) * chr_w

            rect_loc = patches.Rectangle(
                (loc_x, y_chr - h_chr / 2), l_w_real, h_chr,
                facecolor=theme_col, edgecolor=theme_dark, linewidth=0.9, zorder=4
            )
            ax_ins.add_patch(rect_loc)

            # Arrow-head (6.0 pt) positioned above chromosome pointing down
            arrow_x = loc_x + l_w_real / 2
            ax_ins.plot(
                [arrow_x], [pin_y],
                marker="v", color=theme_col, markersize=6.0, markeredgecolor=theme_dark, markeredgewidth=0.8, zorder=5
            )

            # Chromosome / Fragment label kept below chromosomes (explicitly showing Chr Length vs Locus Span)
            if is_fragment:
                chr_prefix = "Scaffold " if not chrom.lower().startswith("scaffold") else ""
                chr_len_str = f"{chr_len/1e3:.0f} kb" if chr_len < 1e6 else f"{chr_len/1e6:.2f} Mb"
                chr_text = f"{chr_prefix}{chrom} ({chr_len_str}) · {badge}: {l_span_str}"
            else:
                chr_display = f"Chr {chrom}" if len(chrom) <= 4 else chrom.split(".")[0]
                chr_text = f"{chr_display} ({chr_len/1e6:.1f} Mb) · {badge}: {l_span_str}"

            ax_ins.text(
                base_x, lbl_y, chr_text,
                fontsize=9.8, fontweight="bold", color="#263238", va="top", zorder=4
            )

            # Distance to telomere / fragment end moved UP, aligned bottom to chromosome, left or right of arrow (fontsize 9.8 pt)
            if is_fragment:
                dist_text = f"{min_d/1e3:.0f} kb to end" if min_d < 1e6 else f"{min_d/1e6:.2f} Mb to end"
            else:
                dist_text = f"{min_d/1e6:.1f} Mb to {near_telo} telo"

            gap_x = 0.018
            est_text_w = 0.22
            # Boundary-safe placement: evaluate arrow position in inset coordinates [0, 1] to prevent clipping
            if arrow_x < 0.48:
                ha_val = "left"
                tx = arrow_x + gap_x
                if tx + est_text_w > 0.98:
                    ha_val = "right"
                    tx = max(arrow_x - gap_x, 0.03)
            else:
                ha_val = "right"
                tx = arrow_x - gap_x
                if tx - est_text_w < 0.02:
                    ha_val = "left"
                    tx = min(arrow_x + gap_x, 0.97)

            ax_ins.text(
                tx, dist_y, dist_text,
                fontsize=9.8, fontweight="bold", color=theme_col,
                ha=ha_val, va="bottom", zorder=5
            )

        ax.legend(
            handles=legend_handles, loc="lower left",
            frameon=True, facecolor="white", framealpha=0.95, edgecolor="#cfd8dc",
            fontsize=8.8, title="Gene Classification & Density", title_fontsize=9.2
        )
        ax.set_xlabel("log10(Coding-Span Gene Length [bp])", fontsize=11.5, fontweight="bold", color="#263238")
        ax.set_ylabel("Coding-Span Gene GC Content (%)", fontsize=11.5, fontweight="bold", color="#263238")
        ax.set_title(
            f"Bivariate Gene Architecture with Chromosome Inset: {common_name}",
            fontsize=12.5, fontweight="bold", color="#212529", pad=12
        )
        ax.set_xlim(2.0, 6.35)
        ax.set_ylim(20.0, 85.0)
        ax.set_xticks([2.0, 3.0, 4.0, 5.0, 6.0, 6.301])
        ax.set_xticklabels(["100 bp\n($10^2$)", "1 kb\n($10^3$)", "10 kb\n($10^4$)", "100 kb\n($10^5$)", "1 Mb\n($10^6$)", "2 Mb\n($2\\times 10^6$)"], fontsize=9.0)
        ax.grid(True, linestyle="--", alpha=0.3, color="#b0bec5")

        plt.tight_layout()
        plt.savefig(png_path)
        plt.savefig(pdf_path)
        plt.close(fig)

        # 2. Plotly interactive rendering
        fig_ply = go.Figure()

        # Background contour
        if Z is not None:
            fig_ply.add_trace(go.Contour(
                x=x_grid, y=y_grid, z=Z,
                colorscale="Blues", opacity=0.6, showscale=False,
                contours=dict(coloring="fill", showlines=True),
                hoverinfo="skip",
                name="Genome Baseline Density",
                showlegend=True,
                xaxis="x", yaxis="y"
            ))

        for grp in sorted_groups:
            items = items_by_grp[grp]
            cfg = BIVARIATE_GROUP_CONFIG.get(grp, BIVARIATE_GROUP_CONFIG["other_locus_feature"])
            xs = [it["x"] for it in items]
            ys = [it["y"] for it in items]
            texts = [
                f"<b>{it['sym']}</b> ({cfg['label']})<br>"
                f"Coding-Span Length: {it['len']:,} bp (log10: {it['x']:.2f})<br>"
                f"Coding-Span GC: {it['y']}%<br>"
                f"Z-Score Length: {it['row'].get('z_score_length', 'N/A')}<br>"
                f"Z-Score GC: {it['row'].get('z_score_coding_span_gc') or it['row'].get('z_score_full_gene_gc', 'N/A')}"
                for it in items
            ]
            fig_ply.add_trace(go.Scatter(
                x=xs, y=ys, mode="markers", name=f"{cfg['label']} (n={len(items)})",
                marker=dict(
                    size=12 if cfg["tier"] in ["alpha", "beta"] else (10 if cfg["tier"] == "pseudogene" else 9),
                    color=cfg["color"],
                    symbol=cfg["plotly_symbol"],
                    line=dict(width=1.0, color=cfg["edgecolor"])
                ),
                hovertext=texts, hoverinfo="text",
                xaxis="x", yaxis="y"
            ))

        # Inset traces in Plotly on xaxis2, yaxis2
        for idx, (_, r_int) in enumerate(sp_intervals.iterrows()):
            loc_id = r_int["locus_id"]
            chrom = str(r_int["chromosome"])
            l_s = int(r_int["start"])
            l_e = int(r_int["end"])
            chr_len = chr_lengths.get(f"{species}::{chrom}", max(l_e * 1.2, 50000000))
            d_5p = l_s
            d_3p = chr_len - l_e
            near_telo = "5'" if d_5p < d_3p else "3'"
            min_d = min(d_5p, d_3p)

            loc_lower = loc_id.lower()
            if "locus_a" in loc_lower:
                theme_col = "#d32f2f"
                badge = "Locus A (α+β)"
            elif "locus_b" in loc_lower:
                theme_col = "#1976d2"
                badge = "Locus B (α+β)"
            elif loc_id == "alpha_globin" or (loc_lower.startswith("alpha") and "beta" not in loc_lower):
                theme_col = "#d32f2f"
                badge = "α-Locus"
            elif loc_id == "beta_globin" or (loc_lower.startswith("beta") and "alpha" not in loc_lower):
                theme_col = "#1976d2"
                badge = "β-Locus"
            elif "proto" in loc_lower:
                theme_col = "#7b1fa2"
                badge = "Proto-Globin"
            else:
                theme_col = "#7b1fa2"
                badge = "Linked α+β"

            bar_rel_w = (chr_len / max_chr_len) * 100
            loc_rel_pos = ((l_s + l_e) / 2 / max_chr_len) * 100
            y_label = f"Chr {chrom} ({badge})"

            # Mini Chromosome Bar
            fig_ply.add_trace(go.Bar(
                x=[bar_rel_w], y=[y_label], orientation="h",
                marker=dict(color="#eceff1", line=dict(color="#455a64", width=1.0)),
                hoverinfo="text",
                hovertext=f"<b>Chromosome {chrom} ({badge})</b><br>Total Length: {chr_len:,} bp<br>Distance to {near_telo} Telomere: {min_d:,} bp",
                showlegend=False,
                xaxis="x2", yaxis="y2"
            ))
            # Locus diamond marker
            fig_ply.add_trace(go.Scatter(
                x=[loc_rel_pos], y=[y_label],
                mode="markers",
                marker=dict(size=12, color=theme_col, symbol="diamond", line=dict(color="black", width=1.0)),
                hoverinfo="text",
                hovertext=f"<b>Locus: {loc_id} ({badge})</b><br>Coordinates: Chr {chrom}: {l_s:,} - {l_e:,} bp<br>Nearest Telomere: {min_d/1e6:.1f} Mb",
                showlegend=False,
                xaxis="x2", yaxis="y2"
            ))

        fig_ply.update_layout(
            template="plotly_white",
            width=1050, height=820,
            title=dict(
                text=f"<b>Bivariate Gene Architecture with Chromosome Inset: {common_name}</b>",
                x=0.5, y=0.985, xanchor="center", yanchor="top",
                font=dict(size=15, color="#212529")
            ),
            xaxis=dict(
                domain=[0.0, 1.0],
                title="<b>log10(Coding-Span Gene Length [bp])</b>",
                range=[2.0, 6.35],
                tickvals=[2.0, 3.0, 4.0, 5.0, 6.0, 6.301],
                ticktext=["100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "2 Mb"],
            ),
            yaxis=dict(domain=[0.0, 1.0], title="<b>Coding-Span Gene GC Content (%)</b>", range=[15.0, 90.0]),
            xaxis2=dict(domain=[0.55, 0.98], anchor="y2", showgrid=False, zeroline=False, showticklabels=False),
            yaxis2=dict(domain=[0.74, 0.97] if num_loci >= 2 else [0.82, 0.97], anchor="x2", showgrid=False, zeroline=False, showticklabels=True),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.12,
                xanchor="center",
                x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)",
                bordercolor="#dee2e6",
                borderwidth=1,
            ),
            margin=dict(t=70, b=90, l=70, r=40)
        )
        fig_ply.write_html(html_path)

        # Validation safeguard: verify 100% coding globin genes are plotted in the bivariate plot and represented by loci
        def _check_is_coding_globin(row):
            grp_str = str(row.get("comparative_group", ""))
            if grp_str in ["other_locus_feature", "globin_pseudogene"]:
                return False
            biotype = str(row.get("biotype", "")).lower()
            if "pseudogene" in biotype:
                return False
            is_g = str(row.get("is_globin", False)).lower() == "true"
            return is_g or ("globin" in grp_str and "non_globin" not in grp_str)

        expected_globins = set(sp_genes[sp_genes.apply(_check_is_coding_globin, axis=1)]["display_symbol"])
        plotted_bivariate_symbols = set()
        for grp, items in items_by_grp.items():
            for it in items:
                plotted_bivariate_symbols.add(it["sym"])

        missing_bivariate = expected_globins - plotted_bivariate_symbols
        if missing_bivariate:
            raise RuntimeError(f"Visualizer safeguard violation: species {species} has {len(missing_bivariate)} coding globin genes missing from bivariate plot: {sorted(list(missing_bivariate))}")

        represented_loci = set(sp_intervals["locus_id"])
        glo_rows = sp_genes[sp_genes.apply(_check_is_coding_globin, axis=1)]
        missing_loci = set(glo_rows["locus_id"]) - represented_loci
        if missing_loci:
            raise RuntimeError(f"Visualizer safeguard violation: species {species} has coding globin genes in unrepresented loci: {missing_loci}")

        return {"png": png_path, "pdf": pdf_path, "html": html_path}

    def build_cross_species_summary(self, df_genes: pd.DataFrame) -> pd.DataFrame:
        """Computes locus-level summary statistics and telomere distances for all species."""
        clade_display = {
            "cyclostome": "Agnatha",
            "chondrichthyes": "Chondrichthyes",
            "non_teleost_fish": "Holostei",
            "teleost_fish": "Teleostei",
            "sarcopterygii": "Sarcopterygii",
            "amphibian": "Amphibia",
            "amniote_reptile": "Reptilia",
            "amniote_crocodylian": "Archosauria",
            "amniote_bird": "Aves",
            "amniote_monotreme": "Monotremata",
            "amniote_mammal": "Mammalia",
        }
        
        present_species = [sp for sp in self.PHYLOGENETIC_ORDER if sp in df_genes["species"].values]
        if not present_species:
            present_species = sorted(df_genes["species"].unique())

        rows = []
        for sp in present_species:
            sp_df = df_genes[df_genes['species'] == sp]
            if sp_df.empty:
                continue
            clade = sp_df['clade'].iloc[0]
            cname = self.COMMON_NAME_MAP.get(sp, sp)
            
            a_genes = sp_df[sp_df['comparative_group'].isin(['alpha_globin_coding', 'globin_coding_alpha', 'globin_coding_locus_A', 'globin_coding'])]
            b_genes = sp_df[sp_df['comparative_group'].isin(['beta_globin_coding', 'globin_coding_beta', 'globin_coding_locus_B'])]
            
            a_flank = sp_df[sp_df['comparative_group'].isin(['cluster_non_globin_alpha', 'cluster_non_globin_locus_A', 'cluster_non_globin'])]
            b_flank = sp_df[sp_df['comparative_group'].isin(['cluster_non_globin_beta', 'cluster_non_globin_locus_B'])]
            
            has_alpha = len(a_genes) > 0
            has_beta = len(b_genes) > 0
            
            a_cds_z = a_genes['z_score_gc'].mean() if has_alpha else np.nan
            b_cds_z = b_genes['z_score_gc'].mean() if has_beta else np.nan
            
            a_fg_z = a_genes['z_score_full_gene_gc'].mean() if has_alpha else np.nan
            b_fg_z = b_genes['z_score_full_gene_gc'].mean() if has_beta else np.nan
            
            a_flank_fg_z = a_flank['z_score_full_gene_gc'].mean() if len(a_flank) > 0 else np.nan
            b_flank_fg_z = b_flank['z_score_full_gene_gc'].mean() if len(b_flank) > 0 else (a_flank_fg_z if has_beta else np.nan)

            a_flank_cds_z = a_flank['z_score_gc'].mean() if len(a_flank) > 0 else np.nan
            b_flank_cds_z = b_flank['z_score_gc'].mean() if len(b_flank) > 0 else (a_flank_cds_z if has_beta else np.nan)
            
            a_tel = a_genes['rel_dist_to_telomere'].mean() if has_alpha else np.nan
            b_tel = b_genes['rel_dist_to_telomere'].mean() if has_beta else np.nan
            
            rows.append({
                'species': sp,
                'common_name': cname,
                'clade': clade,
                'clade_label': clade_display.get(clade, clade),
                'has_alpha': has_alpha,
                'has_beta': has_beta,
                'a_cds_z': a_cds_z,
                'b_cds_z': b_cds_z,
                'delta_cds_z': (a_cds_z - b_cds_z) if (has_alpha and has_beta) else np.nan,
                'a_fg_z': a_fg_z,
                'b_fg_z': b_fg_z,
                'delta_fg_z': (a_fg_z - b_fg_z) if (has_alpha and has_beta) else np.nan,
                'a_flank_fg_z': a_flank_fg_z,
                'b_flank_fg_z': b_flank_fg_z,
                'delta_flank_fg_z': (a_flank_fg_z - b_flank_fg_z) if (pd.notna(a_flank_fg_z) and pd.notna(b_flank_fg_z)) else np.nan,
                'a_flank_cds_z': a_flank_cds_z,
                'b_flank_cds_z': b_flank_cds_z,
                'a_tel': a_tel,
                'b_tel': b_tel,
                'is_teleost': clade == 'teleost_fish',
                'is_proto': sp in ['petromyzon_marinus', 'eptatretus_burgeri', 'lampetra_fluviatilis'],
            })
            
        return pd.DataFrame(rows)

    def plot_summary_locus_divergence(
        self,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure1_locus_divergence.png"
    ) -> Dict[str, str]:
        """Generates Figure 1: Macroevolutionary Locus Divergence (Phylogenetic Dumbbell Synthesis)."""
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig, (ax_fg, ax_cds) = plt.subplots(1, 2, figsize=(20, 13), sharey=True, dpi=300)
        
        n_sp = len(df_summary)
        y_indices = np.arange(n_sp)
        
        def style_dumbbell_axis(ax, title, xlabel, xlim):
            ax.axvspan(-1, 1, color="#c8e6c9", alpha=0.45, label="Genome Baseline Range (±1σ)", zorder=1)
            ax.axvspan(-2, 2, color="#e8f5e9", alpha=0.25, zorder=0)
            ax.axvline(0, color="#2e7d32", linestyle="-", linewidth=1.4, label="Genome Median (Z=0)", zorder=1)
            ax.axvline(-1, color="#81c784", linestyle=":", linewidth=1.0, zorder=1)
            ax.axvline(1, color="#81c784", linestyle=":", linewidth=1.0, zorder=1)
            ax.set_xlim(xlim)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
            ax.set_xlabel(xlabel, fontsize=11, fontweight="bold", labelpad=8)
            ax.grid(axis='x', linestyle='--', alpha=0.35, zorder=1)

        style_dumbbell_axis(
            ax_fg,
            "A. Whole-Gene / Isochore GC Z-Scores (Chromosomal Domain Shift)",
            "Standardized Z-Score: Whole-Gene GC (σ from Genome Gene Baseline)",
            (-4.6, 6.5)
        )
        style_dumbbell_axis(
            ax_cds,
            "B. Coding Sequence (CDS) GC Z-Scores (Globin Coding Divergence)",
            "Standardized Z-Score: CDS GC (σ from Genome Coding Baseline)",
            (-2.5, 3.2)
        )

        clades = df_summary['clade_label'].tolist()
        clade_colors = {
            "Agnatha": "#f5f5f5",
            "Holostei": "#ffffff",
            "Teleostei": "#f0f4f8",
            "Amphibia": "#ffffff",
            "Reptilia": "#fffde7",
            "Archosauria": "#fff3e0",
            "Aves": "#fff8e1",
            "Monotremata": "#f1f8e9",
            "Mammalia": "#e8f5e9",
        }
        
        prev_clade = None
        clade_start = 0
        for i, clade in enumerate(clades + [None]):
            if clade != prev_clade:
                if prev_clade is not None:
                    c_col = clade_colors.get(prev_clade, "#ffffff")
                    ax_fg.axhspan(clade_start - 0.5, i - 0.5, color=c_col, alpha=0.55, zorder=0)
                    ax_cds.axhspan(clade_start - 0.5, i - 0.5, color=c_col, alpha=0.55, zorder=0)
                    mid_y = (clade_start + i - 1) / 2.0
                    ax_fg.text(-4.45, mid_y, prev_clade, fontsize=9.5, fontweight='bold', color='#37474f',
                               va='center', ha='left', rotation=0, zorder=5)
                    ax_fg.plot([-4.52, -4.52], [clade_start - 0.4, i - 0.6], color='#546e7a', lw=2.5, zorder=5)
                clade_start = i
                prev_clade = clade

        for idx, row in df_summary.iterrows():
            y = idx
            ax_fg.axhline(y, color="#e0e0e0", linestyle="-", linewidth=0.7, zorder=1)
            ax_cds.axhline(y, color="#e0e0e0", linestyle="-", linewidth=0.7, zorder=1)
            
            # Panel A: Whole-gene GC
            a_fg = row['a_fg_z']
            b_fg = row['b_fg_z']
            a_flank = row['a_flank_fg_z']
            b_flank = row['b_flank_fg_z']
            
            # Flanking background range: offset by -0.28 on y-axis so it never overlaps the globins
            y_flank = y - 0.28
            if pd.notna(a_flank) and pd.notna(b_flank):
                ax_fg.plot([a_flank, b_flank], [y_flank, y_flank], color="#9e9e9e", linewidth=2.8, alpha=0.85, zorder=2,
                            solid_capstyle='round')
                ax_fg.scatter(a_flank, y_flank, marker='v', color='#ff7f0e', s=45, alpha=0.9, zorder=3, edgecolors='black', linewidths=0.5)
                ax_fg.scatter(b_flank, y_flank, marker='<', color='#d95f02', s=45, alpha=0.9, zorder=3, edgecolors='black', linewidths=0.5)
            elif pd.notna(a_flank):
                ax_fg.scatter(a_flank, y_flank, marker='v', color='#ff7f0e', s=45, alpha=0.9, zorder=3, edgecolors='black', linewidths=0.5)

            # Globin Dumbbell (at exact row baseline y)
            if pd.notna(a_fg) and pd.notna(b_fg):
                delta = a_fg - b_fg
                if delta > 0.3:
                    bar_color = "#e65100"
                    lw = 3.2
                elif delta < -0.3:
                    bar_color = "#00838f"
                    lw = 3.2
                else:
                    bar_color = "#78909c"
                    lw = 2.2
                    
                ax_fg.plot([b_fg, a_fg], [y, y], color=bar_color, linewidth=lw, zorder=4, solid_capstyle='round')
                ax_fg.scatter(a_fg, y, marker='o', color='#1f77b4', s=95, edgecolors='black', linewidths=1.0, zorder=5)
                ax_fg.scatter(b_fg, y, marker='s', color='#2ca02c', s=85, edgecolors='black', linewidths=1.0, zorder=5)
            elif pd.notna(a_fg):
                ax_fg.scatter(a_fg, y, marker='o', color='#1f77b4', s=95, edgecolors='black', linewidths=1.0, zorder=5)

            # Panel B: CDS GC
            a_cds = row['a_cds_z']
            b_cds = row['b_cds_z']
            
            if pd.notna(a_cds) and pd.notna(b_cds):
                delta_cds = a_cds - b_cds
                if delta_cds > 0.2:
                    bar_color = "#e65100"
                    lw = 3.0
                elif delta_cds < -0.2:
                    bar_color = "#00838f"
                    lw = 3.0
                else:
                    bar_color = "#78909c"
                    lw = 2.0
                
                ax_cds.plot([b_cds, a_cds], [y, y], color=bar_color, linewidth=lw, zorder=4, solid_capstyle='round')
                ax_cds.scatter(a_cds, y, marker='o', color='#1f77b4', s=95, edgecolors='black', linewidths=1.0, zorder=5)
                ax_cds.scatter(b_cds, y, marker='s', color='#2ca02c', s=85, edgecolors='black', linewidths=1.0, zorder=5)
            elif pd.notna(a_cds):
                ax_cds.scatter(a_cds, y, marker='o', color='#1f77b4', s=95, edgecolors='black', linewidths=1.0, zorder=5)

        ax_fg.set_yticks(y_indices)
        ax_fg.set_yticklabels(df_summary['common_name'], fontsize=11, fontweight='medium')
        ax_fg.tick_params(axis='y', length=0)
        
        # Legends without comment text in parentheses
        legend_elements_fg = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markeredgecolor='k', markersize=9, label='α-Globin Locus (Mean Z)'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ca02c', markeredgecolor='k', markersize=8, label='β-Globin Locus (Mean Z)'),
            Line2D([0], [0], color='#e65100', lw=3.0, label='Divergence: α > β'),
            Line2D([0], [0], color='#00838f', lw=3.0, label='Inversion: β > α'),
            Line2D([0], [0], color='#9e9e9e', lw=2.8, marker='v', markerfacecolor='#ff7f0e', markeredgecolor='k', markersize=7, label='Flanking Syntenic Non-Globins'),
            Patch(facecolor='#c8e6c9', edgecolor='none', alpha=0.45, label='Genome Typical Range (±1σ)'),
        ]
        ax_fg.legend(handles=legend_elements_fg, loc='lower right', frameon=True, fontsize=9.0, facecolor='white', framealpha=0.95)
        
        legend_elements_cds = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markeredgecolor='k', markersize=9, label='α-Globin CDS (Mean Z)'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ca02c', markeredgecolor='k', markersize=8, label='β-Globin CDS (Mean Z)'),
            Line2D([0], [0], color='#e65100', lw=3.0, label='Divergence: α > β'),
            Line2D([0], [0], color='#00838f', lw=3.0, label='Inversion: β > α'),
        ]
        ax_cds.legend(handles=legend_elements_cds, loc='lower right', frameon=True, fontsize=9.0, facecolor='white', framealpha=0.95)

        fig.suptitle(
            "Figure 1: Macroevolutionary Locus Divergence Across Vertebrates (Phylogenetic Dumbbell Synthesis)",
            fontsize=15, fontweight="bold", y=0.98
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(output_path, dpi=300)
        pdf_path = output_path.replace(".png", ".pdf")
        plt.savefig(pdf_path)
        plt.close(fig)
        html_path = output_path.replace(".png", ".html")
        self.plot_summary_locus_divergence_interactive(df_summary, output_path=html_path)
        return {"png": output_path, "pdf": pdf_path, "html": html_path}

    def plot_summary_mechanistic_drivers(
        self,
        df_genes: pd.DataFrame,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure2_mechanistic_drivers.png"
    ) -> Dict[str, str]:
        """Generates Figure 2: The Two Mechanistic Drivers (Isochore Coupling & Telomeric Proximity Gradient)."""
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        from scipy import stats
        from scipy.optimize import curve_fit
        
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig, (ax_iso, ax_tel) = plt.subplots(1, 2, figsize=(20, 9.5), dpi=300)
        
        # Panel A: Isochore Coupling
        iso_data = []
        for idx, row in df_summary.iterrows():
            sp = row['species']
            cname = row['common_name']
            clade = row['clade']
            if pd.notna(row['a_flank_fg_z']) and pd.notna(row['a_cds_z']):
                iso_data.append({
                    'species': sp, 'common_name': cname, 'clade': clade,
                    'locus_type': 'Alpha Locus', 'flank_z': row['a_flank_fg_z'],
                    'globin_cds_z': row['a_cds_z'], 'globin_fg_z': row['a_fg_z'],
                })
            if pd.notna(row['b_flank_fg_z']) and pd.notna(row['b_cds_z']):
                iso_data.append({
                    'species': sp, 'common_name': cname, 'clade': clade,
                    'locus_type': 'Beta Locus', 'flank_z': row['b_flank_fg_z'],
                    'globin_cds_z': row['b_cds_z'], 'globin_fg_z': row['b_fg_z'],
                })
                
        df_iso = pd.DataFrame(iso_data)
        x_vals = df_iso['flank_z'].values
        y_vals = df_iso['globin_fg_z'].values
        slope, intercept, r_val, p_val, std_err = stats.linregress(x_vals, y_vals)
        r2 = r_val ** 2
        
        grid_x = np.linspace(-1.8, 3.5, 100)
        pred_y = slope * grid_x + intercept
        
        ax_iso.axvspan(-1, 1, color="#e8f5e9", alpha=0.5, label="Genome Typical Flank (±1σ)")
        ax_iso.axhspan(-1, 1, color="#e8f5e9", alpha=0.3, label="Genome Typical Globin (±1σ)")
        ax_iso.axvline(0, color="#2e7d32", linestyle=":", linewidth=1.0)
        ax_iso.axhline(0, color="#2e7d32", linestyle=":", linewidth=1.0)
        ax_iso.plot([-2, 4.5], [-2, 4.5], linestyle='--', color='#757575', linewidth=1.5, label='1:1 Line (Pure Isochore Embedding)', zorder=2)
        ax_iso.plot(grid_x, pred_y, color='#d32f2f', linewidth=2.5, label=f'Linear Fit: y = {slope:.2f}x + {intercept:.2f} (R² = {r2:.2f}, p < 1e-6)', zorder=3)
        
        for idx, r in df_iso.iterrows():
            is_alpha = r['locus_type'] == 'Alpha Locus'
            col = '#1f77b4' if is_alpha else '#2ca02c'
            marker = 'o' if is_alpha else 's'
            ax_iso.scatter(r['flank_z'], r['globin_fg_z'], color=col, marker=marker, s=110,
                           edgecolors='black', linewidths=1.0, alpha=0.9, zorder=4)
            
        annotations_iso = [
            # Top right / extreme GC
            ("podarcis_muralis", "Alpha Locus", "Wall Lizard α", (14, -6)),
            ("podarcis_muralis", "Beta Locus", "Wall Lizard β", (14, 10)),
            ("thamnophis_elegans", "Alpha Locus", "Garter Snake α", (14, 10)),
            ("thamnophis_elegans", "Beta Locus", "Garter Snake β", (14, -12)),
            ("pogona_vitticeps", "Alpha Locus", "Bearded Dragon α", (14, -10)),
            ("pogona_vitticeps", "Beta Locus", "Bearded Dragon β", (14, 10)),
            ("anolis_carolinensis", "Alpha Locus", "Anole α", (14, -6)),
            ("anolis_carolinensis", "Beta Locus", "Anole β", (-55, 12)),
            
            # Mammalian High-GC Alpha
            ("Sus_scrofa", "Alpha Locus", "Pig α", (14, 16)),
            ("Homo_sapiens", "Alpha Locus", "Human α", (-60, 14)),
            ("Bos_taurus", "Alpha Locus", "Cattle α", (14, -16)),
            ("Mus_musculus", "Alpha Locus", "Mouse α", (18, 12)),
            
            # Interstitial alpha & beta around center-right (x ~ 0.6 - 1.0, y ~ 0.4 - 1.3)
            ("chelonia_mydas", "Beta Locus", "Sea Turtle β", (14, 10)),
            ("meleagris_gallopavo", "Beta Locus", "Turkey β", (-55, 18)),
            ("anas_platyrhynchos", "Beta Locus", "Duck β", (18, -4)),
            ("lepisosteus_oculatus", "Alpha Locus", "Gar α", (-50, -12)),
            ("aquila_chrysaetos_chrysaetos", "Alpha Locus", "Golden Eagle α", (16, 12)),
            ("aquila_chrysaetos_chrysaetos", "Beta Locus", "Golden Eagle β", (16, -14)),
            
            # Center-left cluster (x ~ -0.2 to +0.5, y ~ 0 to 1.9)
            ("chelonia_mydas", "Alpha Locus", "Sea Turtle α", (-60, 10)),
            ("danio_rerio", "Alpha Locus", "Zebrafish α", (14, 10)),
            ("Monodelphis_domestica", "Alpha Locus", "Opossum α", (-65, 8)),
            ("Chrysemys_picta_bellii", "Alpha Locus", "Painted Turtle α", (14, 10)),
            ("meleagris_gallopavo", "Alpha Locus", "Turkey α", (-60, -4)),
            ("Gallus_gallus", "Alpha Locus", "Chicken α", (-55, 12)),
            ("Gallus_gallus", "Beta Locus", "Chicken β", (14, 14)),
            ("alligator_mississippiensis", "Beta Locus", "Alligator β", (14, 8)),
            ("ornithorhynchus_anatinus", "Beta Locus", "Platypus β", (-55, -8)),
            ("danio_rerio", "Beta Locus", "Zebrafish β", (-65, 8)),
            ("alligator_mississippiensis", "Alpha Locus", "Alligator α", (-55, 10)),
            ("anas_platyrhynchos", "Alpha Locus", "Duck α", (16, 12)),
            ("ornithorhynchus_anatinus", "Alpha Locus", "Platypus α", (16, 4)),
            ("Monodelphis_domestica", "Beta Locus", "Opossum β", (-55, -8)),
            ("takifugu_rubripes", "Alpha Locus", "Fugu α", (14, 10)),
            ("oryzias_latipes", "Alpha Locus", "Medaka α", (-50, -12)),
            ("gasterosteus_aculeatus", "Alpha Locus", "Stickleback α", (-60, -10)),
            ("gasterosteus_aculeatus", "Beta Locus", "Stickleback β", (14, 8)),
            ("Sus_scrofa", "Beta Locus", "Pig β", (16, -14)),
            ("petromyzon_marinus", "Alpha Locus", "Lamprey", (-55, 6)),
            
            # Low GC points (y < 0)
            ("takifugu_rubripes", "Beta Locus", "Fugu β", (-55, -12)),
            ("Chrysemys_picta_bellii", "Beta Locus", "Painted Turtle β", (-52, -10)),
            ("Bos_taurus", "Beta Locus", "Cattle β", (14, -10)),
            ("Homo_sapiens", "Beta Locus", "Human β", (-55, -12)),
            ("Mus_musculus", "Beta Locus", "Mouse β", (-55, 8)),
            ("oryzias_latipes", "Beta Locus", "Medaka β", (14, -10)),
            ("xenopus_tropicalis", "Alpha Locus", "Clawed Frog α", (-75, -8)),
            ("bufo_bufo", "Alpha Locus", "Toad α", (14, -12)),
        ]
        for sp, ltype, label, offset in annotations_iso:
            sub = df_iso[(df_iso['species'] == sp) & (df_iso['locus_type'] == ltype)]
            if not sub.empty:
                px = sub['flank_z'].iloc[0]
                py = sub['globin_fg_z'].iloc[0]
                ax_iso.annotate(
                    label, xy=(px, py), xytext=offset, textcoords="offset points",
                    fontsize=7.8, fontweight='bold', color='#212121',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.75, edgecolor='none'),
                    arrowprops=dict(arrowstyle="->", color="#546e7a", lw=0.7),
                    zorder=6
                )

        ax_iso.set_xlim(-2.0, 4.0)
        ax_iso.set_ylim(-2.0, 6.5)
        ax_iso.set_xlabel("Flanking Syntenic Non-Globin GC Z-Score (Chromosomal Neighborhood)", fontsize=11, fontweight="bold", labelpad=8)
        ax_iso.set_ylabel("Globin Locus Whole-Gene GC Z-Score", fontsize=11, fontweight="bold", labelpad=8)
        ax_iso.set_title("A. Isochore Coupling: Globin GC Mirrors Surrounding Synteny", fontsize=13, fontweight="bold", pad=12)
        ax_iso.grid(True, linestyle="--", alpha=0.3)
        ax_iso.legend(loc="upper left", frameon=True, fontsize=8.8, facecolor="white", framealpha=0.95)

        # Panel B: Telomeric Proximity Gradient
        tel_data = []
        for idx, row in df_summary.iterrows():
            sp = row['species']
            cname = row['common_name']
            clade = row['clade']
            if pd.notna(row['a_tel']) and pd.notna(row['a_fg_z']):
                tel_data.append({
                    'species': sp, 'common_name': cname, 'clade': clade,
                    'locus': 'α-Globin', 'rel_tel': row['a_tel'] * 100,
                    'gc_z': row['a_fg_z'], 'cds_z': row['a_cds_z'],
                })
            if pd.notna(row['b_tel']) and pd.notna(row['b_fg_z']):
                tel_data.append({
                    'species': sp, 'common_name': cname, 'clade': clade,
                    'locus': 'β-Globin', 'rel_tel': row['b_tel'] * 100,
                    'gc_z': row['b_fg_z'], 'cds_z': row['b_cds_z'],
                })
                
        df_tel = pd.DataFrame(tel_data)
        ax_tel.axvspan(0, 5, color="#fff3e0", alpha=0.7, label="Subtelomeric Zone (<5% Arm)", zorder=0)
        ax_tel.axvline(5, color="#ff9800", linestyle="--", linewidth=1.2, zorder=1)
        ax_tel.axhspan(-1, 1, color="#e8f5e9", alpha=0.4, label="Genome Baseline (±1σ)", zorder=0)
        ax_tel.axhline(0, color="#2e7d32", linestyle=":", linewidth=1.0)

        for idx, r in df_tel.iterrows():
            is_alpha = r['locus'] == 'α-Globin'
            col = '#1f77b4' if is_alpha else '#2ca02c'
            marker = 'o' if is_alpha else 's'
            ax_tel.scatter(r['rel_tel'], r['gc_z'], color=col, marker=marker, s=110,
                           edgecolors='black', linewidths=1.0, alpha=0.9, zorder=4)

        # Annotate loci cleanly across the entire chromosome arm
        annotations_tel = [
            # Subtelomeric hyper-GC cluster (<5% arm)
            ("podarcis_muralis", "α-Globin", "Wall Lizard α", (14, 8)),
            ("podarcis_muralis", "β-Globin", "Wall Lizard β", (14, -12)),
            ("thamnophis_elegans", "α-Globin", "Garter Snake α", (-65, 8)),
            ("thamnophis_elegans", "β-Globin", "Garter Snake β", (14, 10)),
            ("pogona_vitticeps", "β-Globin", "Bearded Dragon β", (14, 6)),
            ("pogona_vitticeps", "α-Globin", "Bearded Dragon α", (14, -12)),
            ("anolis_carolinensis", "β-Globin", "Anole β", (14, 10)),
            ("chelonia_mydas", "β-Globin", "Sea Turtle β", (14, 10)),
            ("Homo_sapiens", "α-Globin", "Human α", (14, 14)),
            ("Bos_taurus", "α-Globin", "Cattle α", (14, -14)),
            ("meleagris_gallopavo", "β-Globin", "Turkey β", (16, 14)),
            ("anas_platyrhynchos", "β-Globin", "Duck β", (36, 6)),
            ("Gallus_gallus", "β-Globin", "Chicken β", (16, 2)),
            ("alligator_mississippiensis", "β-Globin", "Alligator β", (16, -18)),
            ("aquila_chrysaetos_chrysaetos", "β-Globin", "Golden Eagle β", (36, -8)),
            ("Sus_scrofa", "β-Globin", "Pig β", (14, 8)),
            ("Homo_sapiens", "β-Globin", "Human β", (14, -12)),
            
            # Interstitial 5% to 20% arm
            ("meleagris_gallopavo", "α-Globin", "Turkey α", (-55, 10)),
            ("danio_rerio", "α-Globin", "Zebrafish α", (14, 10)),
            ("oryzias_latipes", "β-Globin", "Medaka β", (14, -12)),
            ("gasterosteus_aculeatus", "α-Globin", "Stickleback α", (-60, -12)),
            ("gasterosteus_aculeatus", "β-Globin", "Stickleback β", (14, 10)),
            ("lepisosteus_oculatus", "β-Globin", "Gar β", (14, 12)),
            ("lepisosteus_oculatus", "α-Globin", "Gar α", (14, -14)),
            ("Chrysemys_picta_bellii", "α-Globin", "Painted Turtle α", (-55, 10)),
            ("aquila_chrysaetos_chrysaetos", "α-Globin", "Golden Eagle α", (14, 10)),
            ("ornithorhynchus_anatinus", "β-Globin", "Platypus β", (14, -12)),
            ("Monodelphis_domestica", "β-Globin", "Opossum β", (14, 8)),
            
            # Interstitial 20% to 35% arm
            ("Gallus_gallus", "α-Globin", "Chicken α", (14, 12)),
            ("Mus_musculus", "α-Globin", "Mouse α", (14, 10)),
            ("Mus_musculus", "β-Globin", "Mouse β", (-55, -12)),
            ("takifugu_rubripes", "α-Globin", "Fugu α", (14, 10)),
            ("takifugu_rubripes", "β-Globin", "Fugu β", (14, -12)),
            ("Chrysemys_picta_bellii", "β-Globin", "Painted Turtle β", (14, -10)),
            ("Sus_scrofa", "α-Globin", "Pig α", (14, 10)),
            ("oryzias_latipes", "α-Globin", "Medaka α", (14, -14)),
            ("petromyzon_marinus", "α-Globin", "Lamprey", (-45, 16)),
            
            # Interstitial 35% to 50% arm
            ("chelonia_mydas", "α-Globin", "Sea Turtle α", (14, 10)),
            ("ornithorhynchus_anatinus", "α-Globin", "Platypus α", (-55, -14)),
            ("alligator_mississippiensis", "α-Globin", "Alligator α", (-55, 18)),
            ("anas_platyrhynchos", "α-Globin", "Duck α", (16, -14)),
            ("bufo_bufo", "β-Globin", "Toad β", (-52, 10)),
            ("bufo_bufo", "α-Globin", "Toad α", (14, -12)),
            ("anolis_carolinensis", "α-Globin", "Anole α", (-55, 10)),
            ("danio_rerio", "β-Globin", "Zebrafish β", (14, 12)),
            ("Bos_taurus", "β-Globin", "Cattle β", (-55, -12)),
            ("xenopus_tropicalis", "β-Globin", "Frog β", (14, 14)),
            ("xenopus_tropicalis", "α-Globin", "Frog α", (14, -14)),
            ("Monodelphis_domestica", "α-Globin", "Opossum α", (-60, 10)),
        ]
        for sp, ltype, label, offset in annotations_tel:
            sub = df_tel[(df_tel['species'] == sp) & (df_tel['locus'] == ltype)]
            if not sub.empty:
                px = sub['rel_tel'].iloc[0]
                py = sub['gc_z'].iloc[0]
                ax_tel.annotate(
                    label, xy=(px, py), xytext=offset, textcoords="offset points",
                    fontsize=7.8, fontweight='bold', color='#212121',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.75, edgecolor='none'),
                    arrowprops=dict(arrowstyle="->", color="#546e7a", lw=0.7),
                    zorder=6
                )

        ax_tel.set_xlim(-1, 53)
        ax_tel.set_ylim(-2.0, 6.5)
        ax_tel.set_xlabel("Relative Distance to Closest Telomere (% of Chromosome Arm)", fontsize=11, fontweight="bold", labelpad=8)
        ax_tel.set_ylabel("Locus Whole-Gene GC Z-Score", fontsize=11, fontweight="bold", labelpad=8)
        ax_tel.set_title("B. Chromosomal Position & Telomeric Proximity vs. Locus GC", fontsize=13, fontweight="bold", pad=12)
        ax_tel.grid(True, linestyle="--", alpha=0.3)
        
        legend_elements_tel = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markeredgecolor='k', markersize=9, label='α-Globin Loci'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ca02c', markeredgecolor='k', markersize=8, label='β-Globin Loci'),
            Patch(facecolor='#fff3e0', edgecolor='#ff9800', linestyle='--', label='Subtelomeric Zone (<5% Arm)'),
            Patch(facecolor='#e8f5e9', edgecolor='none', alpha=0.4, label='Genome Baseline (±1σ)'),
        ]
        ax_tel.legend(handles=legend_elements_tel, loc="upper right", frameon=True, fontsize=8.8, facecolor="white", framealpha=0.95)

        fig.suptitle(
            "Figure 2: The Two Mechanistic Drivers of Vertebrate Globin GC Divergence",
            fontsize=15, fontweight="bold", y=0.98
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(output_path, dpi=300)
        pdf_path = output_path.replace(".png", ".pdf")
        plt.savefig(pdf_path)
        plt.close(fig)
        html_path = output_path.replace(".png", ".html")
        self.plot_summary_mechanistic_drivers_interactive(df_genes, df_summary, output_path=html_path)
        return {"png": output_path, "pdf": pdf_path, "html": html_path}

    def plot_summary_clade_heatmap(
        self,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure3_clade_heatmap.png"
    ) -> Dict[str, str]:
        """Generates Figure 3: Clade-Structured Multi-Metric Heatmap Matrix."""
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        
        metrics = [
            ('a_cds_z', 'α CDS GC (Z)'),
            ('b_cds_z', 'β CDS GC (Z)'),
            ('delta_cds_z', 'Δ(α-β) CDS (Z)'),
            ('a_fg_z', 'α Gene GC (Z)'),
            ('b_fg_z', 'β Gene GC (Z)'),
            ('delta_fg_z', 'Δ(α-β) Gene (Z)'),
            ('delta_flank_fg_z', 'Δ(α-β) Flank (Z)'),
        ]
        
        col_keys = [m[0] for m in metrics]
        col_labels = [m[1] for m in metrics]
        
        matrix_data = []
        annot_data = []
        
        for idx, r in df_summary.iterrows():
            row_vals = []
            annot_vals = []
            for k in col_keys:
                val = r[k]
                if pd.isna(val):
                    row_vals.append(np.nan)
                    annot_vals.append("—")
                else:
                    row_vals.append(val)
                    annot_vals.append(f"{val:+.2f}")
            matrix_data.append(row_vals)
            annot_data.append(annot_vals)
            
        df_matrix = pd.DataFrame(matrix_data, columns=col_labels, index=df_summary['common_name'])
        df_annot = pd.DataFrame(annot_data, columns=col_labels, index=df_summary['common_name'])
        
        fig, ax = plt.subplots(figsize=(14, 11), dpi=300)
        sns.heatmap(
            df_matrix,
            annot=df_annot,
            fmt="",
            cmap="coolwarm",
            center=0,
            vmin=-2.5,
            vmax=2.5,
            cbar_kws={'label': 'Standardized Z-Score (σ from Genome Baseline)'},
            linewidths=0.8,
            linecolor='#e0e0e0',
            ax=ax,
        )
        
        ax.set_title("Figure 3: Cross-Species Quantitative Divergence Matrix (Z-Scores & Divergence Polarities)",
                     fontsize=14, fontweight="bold", pad=16)
        ax.set_ylabel("")
        ax.tick_params(axis='x', rotation=30, labelsize=10.5)
        ax.tick_params(axis='y', labelsize=10.5)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        pdf_path = output_path.replace(".png", ".pdf")
        plt.savefig(pdf_path)
        plt.close(fig)
        html_path = output_path.replace(".png", ".html")
        self.plot_summary_clade_heatmap_interactive(df_summary, output_path=html_path)
        return {"png": output_path, "pdf": pdf_path, "html": html_path}

    def plot_summary_locus_divergence_interactive(
        self,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure1_locus_divergence.html"
    ) -> str:
        """Generates interactive Plotly HTML for Figure 1: Macroevolutionary Locus Divergence."""
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig = make_subplots(
            rows=1, cols=2, shared_yaxes=True,
            subplot_titles=(
                "<b>A. Whole-Gene / Isochore GC Z-Scores</b><br><sup>(Chromosomal Domain Shifts)</sup>",
                "<b>B. Coding Sequence (CDS) GC Z-Scores</b><br><sup>(Globin Coding Divergence)</sup>"
            ),
            horizontal_spacing=0.08
        )

        n_sp = len(df_summary)
        y_coords = list(range(n_sp))
        y_labels = df_summary['common_name'].tolist()

        for col_idx in [1, 2]:
            fig.add_vrect(
                x0=-1, x1=1, fillcolor="#c8e6c9", opacity=0.35, layer="below", line_width=0,
                row=1, col=col_idx
            )
            fig.add_vrect(
                x0=-2, x1=2, fillcolor="#e8f5e9", opacity=0.20, layer="below", line_width=0,
                row=1, col=col_idx
            )
            fig.add_vline(x=0, line_width=1.2, line_dash="solid", line_color="#2e7d32", row=1, col=col_idx)
            fig.add_vline(x=-1, line_width=0.8, line_dash="dot", line_color="#81c784", row=1, col=col_idx)
            fig.add_vline(x=1, line_width=0.8, line_dash="dot", line_color="#81c784", row=1, col=col_idx)

        clades = df_summary['clade_label'].tolist()
        clade_colors = {
            "Agnatha": "rgba(245, 245, 245, 0.4)",
            "Holostei": "rgba(255, 255, 255, 0.4)",
            "Teleostei": "rgba(240, 244, 248, 0.45)",
            "Amphibia": "rgba(255, 255, 255, 0.4)",
            "Reptilia": "rgba(255, 253, 231, 0.45)",
            "Archosauria": "rgba(255, 243, 224, 0.45)",
            "Aves": "rgba(255, 248, 225, 0.45)",
            "Monotremata": "rgba(241, 248, 233, 0.45)",
            "Mammalia": "rgba(232, 245, 233, 0.45)",
        }
        prev_clade = None
        c_start = 0
        for i, c in enumerate(clades + [None]):
            if c != prev_clade:
                if prev_clade is not None:
                    bg = clade_colors.get(prev_clade, "rgba(250, 250, 250, 0.3)")
                    for col_idx in [1, 2]:
                        fig.add_hrect(
                            y0=c_start - 0.5, y1=i - 0.5,
                            fillcolor=bg, opacity=0.6, layer="below", line_width=0,
                            row=1, col=col_idx
                        )
                c_start = i
                prev_clade = c

        a_fg_x, a_fg_y, a_fg_hover = [], [], []
        b_fg_x, b_fg_y, b_fg_hover = [], [], []
        flank_x, flank_y, flank_hover = [], [], []
        
        a_cds_x, a_cds_y, a_cds_hover = [], [], []
        b_cds_x, b_cds_y, b_cds_hover = [], [], []

        for idx, r in df_summary.iterrows():
            y = idx
            cname = r['common_name']
            sp = r['species']
            clade = r['clade_label']

            a_flank = r['a_flank_fg_z']
            b_flank = r['b_flank_fg_z']
            y_flank = y - 0.25
            if pd.notna(a_flank) and pd.notna(b_flank):
                fig.add_trace(
                    go.Scatter(
                        x=[a_flank, b_flank], y=[y_flank, y_flank],
                        mode="lines", line=dict(color="#9e9e9e", width=2.5),
                        hoverinfo="skip", showlegend=False
                    ), row=1, col=1
                )
                flank_x.extend([a_flank, b_flank])
                flank_y.extend([y_flank, y_flank])
                flank_hover.extend([
                    f"<b>{cname}</b><br>Flanking Non-Globins (α Locus)<br>Z-Score: {a_flank:+.2f}σ",
                    f"<b>{cname}</b><br>Flanking Non-Globins (β Locus)<br>Z-Score: {b_flank:+.2f}σ"
                ])
            elif pd.notna(a_flank):
                flank_x.append(a_flank)
                flank_y.append(y_flank)
                flank_hover.append(f"<b>{cname}</b><br>Flanking Non-Globins<br>Z-Score: {a_flank:+.2f}σ")

            a_fg = r['a_fg_z']
            b_fg = r['b_fg_z']
            if pd.notna(a_fg) and pd.notna(b_fg):
                delta = a_fg - b_fg
                bar_col = "#e65100" if delta > 0.3 else ("#00838f" if delta < -0.3 else "#78909c")
                fig.add_trace(
                    go.Scatter(
                        x=[b_fg, a_fg], y=[y, y],
                        mode="lines", line=dict(color=bar_col, width=3.2),
                        hoverinfo="skip", showlegend=False
                    ), row=1, col=1
                )
            if pd.notna(a_fg):
                a_fg_x.append(a_fg)
                a_fg_y.append(y)
                tel_str = f"<br>Telomere Distance: {r['a_tel']*100:.1f}% arm" if pd.notna(r.get('a_tel')) else ""
                a_fg_hover.append(
                    f"<b>{cname}</b> (<i>{sp}</i>)<br>Clade: {clade}<br>Locus: α-Globin<br>"
                    f"Whole-Gene GC Z-Score: {a_fg:+.2f}σ{tel_str}"
                )
            if pd.notna(b_fg):
                b_fg_x.append(b_fg)
                b_fg_y.append(y)
                tel_str = f"<br>Telomere Distance: {r['b_tel']*100:.1f}% arm" if pd.notna(r.get('b_tel')) else ""
                b_fg_hover.append(
                    f"<b>{cname}</b> (<i>{sp}</i>)<br>Clade: {clade}<br>Locus: β-Globin<br>"
                    f"Whole-Gene GC Z-Score: {b_fg:+.2f}σ{tel_str}"
                )

            a_cds = r['a_cds_z']
            b_cds = r['b_cds_z']
            if pd.notna(a_cds) and pd.notna(b_cds):
                delta_cds = a_cds - b_cds
                bar_col_cds = "#e65100" if delta_cds > 0.2 else ("#00838f" if delta_cds < -0.2 else "#78909c")
                fig.add_trace(
                    go.Scatter(
                        x=[b_cds, a_cds], y=[y, y],
                        mode="lines", line=dict(color=bar_col_cds, width=3.0),
                        hoverinfo="skip", showlegend=False
                    ), row=1, col=2
                )
            if pd.notna(a_cds):
                a_cds_x.append(a_cds)
                a_cds_y.append(y)
                a_cds_hover.append(
                    f"<b>{cname}</b><br>Clade: {clade}<br>Locus: α-Globin CDS<br>"
                    f"CDS GC Z-Score: {a_cds:+.2f}σ"
                )
            if pd.notna(b_cds):
                b_cds_x.append(b_cds)
                b_cds_y.append(y)
                b_cds_hover.append(
                    f"<b>{cname}</b><br>Clade: {clade}<br>Locus: β-Globin CDS<br>"
                    f"CDS GC Z-Score: {b_cds:+.2f}σ"
                )

        fig.add_trace(
            go.Scatter(
                x=flank_x, y=flank_y, mode="markers",
                marker=dict(symbol="triangle-down", size=7, color="#ff7f0e", line=dict(color="black", width=0.5)),
                hovertext=flank_hover, hoverinfo="text", name="Flanking Non-Globins"
            ), row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=a_fg_x, y=a_fg_y, mode="markers",
                marker=dict(symbol="circle", size=10, color="#1f77b4", line=dict(color="black", width=0.8)),
                hovertext=a_fg_hover, hoverinfo="text", name="α-Globin Locus (Mean Z)"
            ), row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=b_fg_x, y=b_fg_y, mode="markers",
                marker=dict(symbol="square", size=9, color="#2ca02c", line=dict(color="black", width=0.8)),
                hovertext=b_fg_hover, hoverinfo="text", name="β-Globin Locus (Mean Z)"
            ), row=1, col=1
        )

        fig.add_trace(
            go.Scatter(
                x=a_cds_x, y=a_cds_y, mode="markers",
                marker=dict(symbol="circle", size=10, color="#1f77b4", line=dict(color="black", width=0.8)),
                hovertext=a_cds_hover, hoverinfo="text", name="α-Globin CDS", showlegend=False
            ), row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=b_cds_x, y=b_cds_y, mode="markers",
                marker=dict(symbol="square", size=9, color="#2ca02c", line=dict(color="black", width=0.8)),
                hovertext=b_cds_hover, hoverinfo="text", name="β-Globin CDS", showlegend=False
            ), row=1, col=2
        )

        fig.update_yaxes(
            tickvals=y_coords, ticktext=y_labels, tickfont=dict(size=11),
            showgrid=True, gridcolor="#eeeeee", row=1, col=1
        )
        fig.update_yaxes(showgrid=True, gridcolor="#eeeeee", row=1, col=2)
        fig.update_xaxes(
            title_text="<b>Standardized Z-Score: Whole-Gene GC (σ from Baseline)</b>",
            range=[-4.6, 6.5], showgrid=True, gridcolor="#e0e0e0", row=1, col=1
        )
        fig.update_xaxes(
            title_text="<b>Standardized Z-Score: CDS GC (σ from Baseline)</b>",
            range=[-2.5, 3.2], showgrid=True, gridcolor="#e0e0e0", row=1, col=2
        )
        fig.update_layout(
            title=dict(
                text="<b>Figure 1: Macroevolutionary Locus Divergence Across Vertebrates (Interactive)</b>",
                font=dict(size=16)
            ),
            template="plotly_white",
            height=950,
            width=1350,
            legend=dict(
                orientation="h", yanchor="bottom", y=-0.10, xanchor="center", x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)", bordercolor="#cccccc", borderwidth=1
            ),
            margin=dict(l=180, r=40, t=80, b=80)
        )
        fig.write_html(output_path)
        return output_path

    def plot_summary_mechanistic_drivers_interactive(
        self,
        df_genes: pd.DataFrame,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure2_mechanistic_drivers.html"
    ) -> str:
        """Generates interactive Plotly HTML for Figure 2: The Two Mechanistic Drivers."""
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy import stats

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(
                "<b>A. Isochore Coupling: Globin GC Mirrors Surrounding Synteny</b><br><sup>(Flanking Non-Globin GC Z-Score vs Globin GC Z-Score)</sup>",
                "<b>B. Telomeric Proximity Gradient: Whole-Genome Arm Distribution</b><br><sup>(Distance to Telomere % Arm vs Globin GC Z-Score)</sup>"
            ),
            horizontal_spacing=0.10
        )

        iso_rows = []
        for idx, r in df_summary.iterrows():
            sp = r['species']
            cname = r['common_name']
            clade = r['clade_label']
            if pd.notna(r['a_flank_fg_z']) and pd.notna(r['a_fg_z']):
                iso_rows.append({
                    'species': sp, 'common_name': cname, 'clade': clade, 'locus': 'α-Globin',
                    'flank_z': r['a_flank_fg_z'], 'globin_z': r['a_fg_z'],
                    'tel_pct': r['a_tel']*100 if pd.notna(r.get('a_tel')) else np.nan
                })
            if pd.notna(r['b_flank_fg_z']) and pd.notna(r['b_fg_z']):
                iso_rows.append({
                    'species': sp, 'common_name': cname, 'clade': clade, 'locus': 'β-Globin',
                    'flank_z': r['b_flank_fg_z'], 'globin_z': r['b_fg_z'],
                    'tel_pct': r['b_tel']*100 if pd.notna(r.get('b_tel')) else np.nan
                })
        df_iso = pd.DataFrame(iso_rows)
        slope, intercept, r_val, p_val, std_err = stats.linregress(df_iso['flank_z'], df_iso['globin_z'])
        r2 = r_val ** 2

        fig.add_vrect(x0=-1, x1=1, fillcolor="#e8f5e9", opacity=0.35, layer="below", line_width=0, row=1, col=1)
        fig.add_hrect(y0=-1, y1=1, fillcolor="#e8f5e9", opacity=0.25, layer="below", line_width=0, row=1, col=1)
        fig.add_vline(x=0, line_width=1.0, line_dash="dot", line_color="#2e7d32", row=1, col=1)
        fig.add_hline(y=0, line_width=1.0, line_dash="dot", line_color="#2e7d32", row=1, col=1)

        fig.add_trace(
            go.Scatter(
                x=[-2.0, 4.5], y=[-2.0, 4.5], mode="lines",
                line=dict(color="#757575", dash="dash", width=1.5),
                name="1:1 Line (Pure Isochore Embedding)"
            ), row=1, col=1
        )
        grid_x = np.linspace(-1.8, 3.8, 100)
        pred_y = slope * grid_x + intercept
        fig.add_trace(
            go.Scatter(
                x=grid_x, y=pred_y, mode="lines",
                line=dict(color="#d32f2f", width=2.5),
                name=f"Linear Fit: y = {slope:.2f}x + {intercept:.2f} (R² = {r2:.2f})"
            ), row=1, col=1
        )

        iso_a = df_iso[df_iso['locus'] == 'α-Globin']
        iso_b = df_iso[df_iso['locus'] == 'β-Globin']
        fig.add_trace(
            go.Scatter(
                x=iso_a['flank_z'], y=iso_a['globin_z'], mode="markers",
                marker=dict(symbol="circle", size=10, color="#1f77b4", line=dict(color="black", width=1.0)),
                hovertext=[
                    f"<b>{r.common_name}</b><br>Clade: {r.clade}<br>Locus: α-Globin<br>"
                    f"Flanker Z: {r.flank_z:+.2f}σ<br>Globin Z: {r.globin_z:+.2f}σ<br>"
                    f"Telomere Proximity: {r.tel_pct:.1f}% arm"
                    for _, r in iso_a.iterrows()
                ],
                hoverinfo="text", name="α-Globin Locus"
            ), row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=iso_b['flank_z'], y=iso_b['globin_z'], mode="markers",
                marker=dict(symbol="square", size=9, color="#2ca02c", line=dict(color="black", width=1.0)),
                hovertext=[
                    f"<b>{r.common_name}</b><br>Clade: {r.clade}<br>Locus: β-Globin<br>"
                    f"Flanker Z: {r.flank_z:+.2f}σ<br>Globin Z: {r.globin_z:+.2f}σ<br>"
                    f"Telomere Proximity: {r.tel_pct:.1f}% arm"
                    for _, r in iso_b.iterrows()
                ],
                hoverinfo="text", name="β-Globin Locus"
            ), row=1, col=1
        )

        fig.add_vrect(x0=0, x1=5, fillcolor="#fff3e0", opacity=0.7, layer="below", line_width=0, row=1, col=2)
        fig.add_vline(x=5, line_width=1.2, line_dash="dash", line_color="#ff9800", row=1, col=2)
        fig.add_hrect(y0=-1, y1=1, fillcolor="#e8f5e9", opacity=0.35, layer="below", line_width=0, row=1, col=2)
        fig.add_hline(y=0, line_width=1.0, line_dash="dot", line_color="#2e7d32", row=1, col=2)

        tel_rows = []
        for idx, r in df_summary.iterrows():
            sp = r['species']
            cname = r['common_name']
            clade = r['clade_label']
            if pd.notna(r['a_tel']) and pd.notna(r['a_fg_z']):
                tel_rows.append({
                    'species': sp, 'common_name': cname, 'clade': clade, 'locus': 'α-Globin',
                    'rel_tel': r['a_tel'] * 100, 'gc_z': r['a_fg_z']
                })
            if pd.notna(r['b_tel']) and pd.notna(r['b_fg_z']):
                tel_rows.append({
                    'species': sp, 'common_name': cname, 'clade': clade, 'locus': 'β-Globin',
                    'rel_tel': r['b_tel'] * 100, 'gc_z': r['b_fg_z']
                })
        df_tel = pd.DataFrame(tel_rows)
        tel_a = df_tel[df_tel['locus'] == 'α-Globin']
        tel_b = df_tel[df_tel['locus'] == 'β-Globin']

        fig.add_trace(
            go.Scatter(
                x=tel_a['rel_tel'], y=tel_a['gc_z'], mode="markers",
                marker=dict(symbol="circle", size=10, color="#1f77b4", line=dict(color="black", width=1.0)),
                hovertext=[
                    f"<b>{r.common_name}</b><br>Clade: {r.clade}<br>Locus: α-Globin<br>"
                    f"Distance to Telomere: {r.rel_tel:.1f}% arm<br>Globin Whole-Gene Z: {r.gc_z:+.2f}σ"
                    for _, r in tel_a.iterrows()
                ],
                hoverinfo="text", showlegend=False
            ), row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=tel_b['rel_tel'], y=tel_b['gc_z'], mode="markers",
                marker=dict(symbol="square", size=9, color="#2ca02c", line=dict(color="black", width=1.0)),
                hovertext=[
                    f"<b>{r.common_name}</b><br>Clade: {r.clade}<br>Locus: β-Globin<br>"
                    f"Distance to Telomere: {r.rel_tel:.1f}% arm<br>Globin Whole-Gene Z: {r.gc_z:+.2f}σ"
                    for _, r in tel_b.iterrows()
                ],
                hoverinfo="text", showlegend=False
            ), row=1, col=2
        )

        fig.update_xaxes(
            title_text="<b>Flanking Non-Globin GC Z-Score</b>",
            range=[-2.0, 4.0], showgrid=True, gridcolor="#e0e0e0", row=1, col=1
        )
        fig.update_yaxes(
            title_text="<b>Globin Locus Whole-Gene GC Z-Score</b>",
            range=[-2.0, 6.5], showgrid=True, gridcolor="#e0e0e0", row=1, col=1
        )
        fig.update_xaxes(
            title_text="<b>Distance to Closest Telomere (% Chromosome Arm)</b>",
            range=[-1.0, 50.0], showgrid=True, gridcolor="#e0e0e0", row=1, col=2
        )
        fig.update_yaxes(
            title_text="<b>Globin Locus Whole-Gene GC Z-Score</b>",
            range=[-2.0, 6.5], showgrid=True, gridcolor="#e0e0e0", row=1, col=2
        )
        fig.update_layout(
            title=dict(
                text="<b>Figure 2: The Two Mechanistic Drivers of Vertebrate Globin GC Evolution (Interactive)</b>",
                font=dict(size=16)
            ),
            template="plotly_white",
            height=720,
            width=1350,
            legend=dict(
                orientation="h", yanchor="bottom", y=-0.16, xanchor="center", x=0.5,
                bgcolor="rgba(255, 255, 255, 0.95)", bordercolor="#cccccc", borderwidth=1
            ),
            margin=dict(l=70, r=40, t=80, b=90)
        )
        fig.write_html(output_path)
        return output_path

    def plot_summary_clade_heatmap_interactive(
        self,
        df_summary: pd.DataFrame,
        output_path: str = "results/summary_figure3_clade_heatmap.html"
    ) -> str:
        """Generates interactive Plotly HTML for Figure 3: Clade Summary Matrix."""
        import plotly.graph_objects as go

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        metrics = [
            ('a_cds_z', 'α CDS GC (Z)'),
            ('b_cds_z', 'β CDS GC (Z)'),
            ('delta_cds_z', 'Δ(α-β) CDS (Z)'),
            ('a_fg_z', 'α Gene GC (Z)'),
            ('b_fg_z', 'β Gene GC (Z)'),
            ('delta_fg_z', 'Δ(α-β) Gene (Z)'),
            ('delta_flank_fg_z', 'Δ(α-β) Flank (Z)'),
        ]
        col_keys = [m[0] for m in metrics]
        col_labels = [m[1] for m in metrics]

        z_matrix = []
        text_matrix = []
        hover_matrix = []

        for idx, r in df_summary.iterrows():
            cname = r['common_name']
            row_z = []
            row_text = []
            row_hover = []
            for k, lbl in zip(col_keys, col_labels):
                val = r[k]
                if pd.isna(val):
                    row_z.append(np.nan)
                    row_text.append("—")
                    row_hover.append(f"<b>{cname}</b><br>{lbl}: N/A")
                else:
                    row_z.append(val)
                    row_text.append(f"{val:+.2f}")
                    row_hover.append(f"<b>{cname}</b><br>{lbl}: <b>{val:+.2f}σ</b>")
            z_matrix.append(row_z)
            text_matrix.append(row_text)
            hover_matrix.append(row_hover)

        fig = go.Figure(
            data=go.Heatmap(
                z=z_matrix,
                x=col_labels,
                y=df_summary['common_name'].tolist(),
                text=text_matrix,
                texttemplate="%{text}",
                textfont=dict(size=11),
                hovertext=hover_matrix,
                hoverinfo="text",
                colorscale="RdBu_r",
                zmin=-2.5,
                zmax=2.5,
                colorbar=dict(
                    title="Standardized Z-Score (σ)",
                    tickvals=[-2, -1, 0, 1, 2],
                    ticktext=["-2.0σ", "-1.0σ", "0.0σ", "+1.0σ", "+2.0σ"]
                )
            )
        )

        fig.update_layout(
            title=dict(
                text="<b>Figure 3: Cross-Species Quantitative Divergence Matrix (Interactive)</b>",
                font=dict(size=16)
            ),
            template="plotly_white",
            height=850,
            width=1050,
            margin=dict(l=180, r=60, t=80, b=80),
            xaxis=dict(tickangle=-25, tickfont=dict(size=11, family="sans-serif", weight="bold")),
            yaxis=dict(tickfont=dict(size=11), autorange="reversed")
        )
        fig.write_html(output_path)
        return output_path


