"""Z-score and relative genomic ranking analyzer module."""

import os
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Any, Optional, Tuple
from globin_pipeline.baseline_manager import BaselineManager


class ZScoreAnalyzer:
    """Computes standardized Z-scores, percentile ranks, and comparative hypothesis tests."""

    def __init__(self, baseline_manager: Optional[BaselineManager] = None):
        self.bm = baseline_manager or BaselineManager()
        self.chr_lengths = {}
        chr_len_path = os.path.join("data", "chromosome_lengths.json")
        if not os.path.exists(chr_len_path):
            chr_len_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "chromosome_lengths.json")
        if os.path.exists(chr_len_path):
            try:
                import json
                with open(chr_len_path, "r", encoding="utf-8") as f:
                    self.chr_lengths = json.load(f)
            except Exception as e:
                print(f"Warning: could not load chromosome lengths from {chr_len_path}: {e}")

    def analyze_genes(
        self, locus_genes: List[Dict[str, Any]], gene_gc_records: List[Dict[str, Any]]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Calculates Z-scores, ranks, and statistical tests for all locus genes."""
        # Index gene_gc_records by (species, gene_id)
        gc_by_gene = {}
        for r in gene_gc_records:
            key = (r["species"], r["gene_id"])
            gc_by_gene[key] = r

        # Cache baseline DFs and stats by species
        baselines_df: Dict[str, pd.DataFrame] = {}
        baselines_stats: Dict[str, Dict[str, Any]] = {}

        results = []

        for g in locus_genes:
            sp = g["species"]
            g_id = g["gene_id"]
            cat = g["category"]
            biotype = str(g.get("biotype", "")).strip()
            locus_id = g.get("locus_id", "alpha_globin")
            is_globin = g.get("is_globin", "False")

            # Load baseline for species
            if sp not in baselines_stats:
                try:
                    df_base, st_base = self.bm.get_genome_baseline(sp)
                    baselines_df[sp] = df_base
                    baselines_stats[sp] = st_base
                except Exception as e:
                    print(f"Warning: could not load baseline for {sp}: {e}")
                    baselines_df[sp] = pd.DataFrame()
                    baselines_stats[sp] = {}

            st = baselines_stats.get(sp, {})
            df_sp = baselines_df.get(sp, pd.DataFrame())

            mu_cds = st.get("mean_cds_gc")
            sigma_cds = st.get("std_cds_gc")
            mu_gc3 = st.get("mean_gc3")
            sigma_gc3 = st.get("std_gc3")
            mu_exon = st.get("mean_exonic_gc", mu_cds)
            sigma_exon = st.get("std_exonic_gc", sigma_cds)
            mu_genomic = st.get("mean_genomic_gc", mu_cds)
            sigma_genomic = st.get("std_genomic_gc", sigma_cds)

            # Get GC metrics for this gene
            gc_rec = gc_by_gene.get((sp, g_id), {})
            cds_gc = gc_rec.get("cds_gc")
            cds_len = gc_rec.get("cds_len", 0)
            gc3 = gc_rec.get("gc3")
            full_gene_gc = gc_rec.get("full_gene_gc") if gc_rec.get("full_gene_gc") is not None else gc_rec.get("genomic_gc")
            full_gene_len = gc_rec.get("full_gene_len", 0)
            exonic_gc = gc_rec.get("exonic_gc")
            exonic_len = gc_rec.get("exonic_len", 0)
            intronic_gc = gc_rec.get("intronic_gc")
            intronic_len = gc_rec.get("intronic_len", 0)

            g_name = str(g.get("gene_name") or "").lower()
            g_desc = str(g.get("description") or "").lower()
            if is_globin == "True" and biotype == "protein_coding":
                if "locus_b" in str(locus_id).lower():
                    comp_group = "globin_coding_locus_B"
                elif "locus_a" in str(locus_id).lower():
                    comp_group = "globin_coding_locus_A"
                elif cat == "alpha_globin":
                    comp_group = "alpha_globin_coding"
                elif cat == "beta_globin":
                    comp_group = "beta_globin_coding"
                elif "alpha" in str(locus_id).lower() and "beta" not in str(locus_id).lower():
                    comp_group = "alpha_globin_coding"
                elif "beta" in str(locus_id).lower() and "alpha" not in str(locus_id).lower():
                    comp_group = "beta_globin_coding"
                else:
                    comp_group = "globin_coding"
            elif cat == "globin_pseudogene" or ("pseudogene" in str(biotype) and ("globin" in g_desc or "hba" in g_name or "hbb" in g_name or "hbz" in g_name or "hbd" in g_name)):
                comp_group = "globin_pseudogene"
            elif cat == "cluster_non_globin" or (biotype == "protein_coding" and is_globin != "True"):
                loc_lower = str(locus_id).lower()
                if "locus_b" in loc_lower:
                    comp_group = "cluster_non_globin_locus_B"
                elif "locus_a" in loc_lower:
                    comp_group = "cluster_non_globin_locus_A"
                elif "alpha_beta" in loc_lower or "cluster" in loc_lower:
                    comp_group = "cluster_non_globin"
                elif "beta" in loc_lower and "alpha" not in loc_lower:
                    comp_group = "cluster_non_globin_beta"
                elif "alpha" in loc_lower and "beta" not in loc_lower:
                    comp_group = "cluster_non_globin_alpha"
                else:
                    comp_group = "cluster_non_globin"
            else:
                comp_group = "other_locus_feature"

            loc_lower = str(locus_id).lower()
            if "locus_b" in loc_lower:
                locus_origin = "locus_B"
            elif "locus_a" in loc_lower:
                locus_origin = "locus_A"
            elif "alpha_beta" in loc_lower or "cluster" in loc_lower:
                locus_origin = "alpha_beta_cluster"
            elif "beta" in loc_lower and "alpha" not in loc_lower:
                locus_origin = "beta_locus"
            elif "alpha" in loc_lower and "beta" not in loc_lower:
                locus_origin = "alpha_locus"
            else:
                locus_origin = "linked_cluster"

            # Primary GC for standardized ranking: CDS GC for coding, exonic/full for pseudogenes
            eval_gc = cds_gc if cds_gc is not None else (exonic_gc if exonic_gc is not None else full_gene_gc)

            # Calculate Z-scores and percentiles
            z_cds_gc = None
            z_gc3 = None
            z_exonic_gc = None
            z_full_gene_gc = None
            delta_full_gene_gc = None
            delta_cds_gc = None
            pct_gc = None
            pct_gc3 = None
            pct_exonic = None

            if eval_gc is not None and mu_cds is not None and sigma_cds and sigma_cds > 0:
                z_cds_gc = round(float((eval_gc - mu_cds) / sigma_cds), 2)
                delta_cds_gc = round(float(eval_gc - mu_cds), 2)
                if not df_sp.empty and "cds_gc" in df_sp.columns:
                    pct_gc = round(float(stats.percentileofscore(df_sp["cds_gc"].dropna(), eval_gc)), 2)

            if gc3 is not None and mu_gc3 is not None and sigma_gc3 and sigma_gc3 > 0:
                z_gc3 = round(float((gc3 - mu_gc3) / sigma_gc3), 2)
                if not df_sp.empty and "gc3" in df_sp.columns:
                    pct_gc3 = round(float(stats.percentileofscore(df_sp["gc3"].dropna(), gc3)), 2)

            if exonic_gc is not None and mu_exon is not None and sigma_exon and sigma_exon > 0:
                z_exonic_gc = round(float((exonic_gc - mu_exon) / sigma_exon), 2)
                if not df_sp.empty and "exonic_gc" in df_sp.columns:
                    pct_exonic = round(float(stats.percentileofscore(df_sp["exonic_gc"].dropna(), exonic_gc)), 2)

            # Empirical Whole-Gene Body baseline (derived from all ~20,000 protein-coding genes)
            mu_whole = st.get("mean_whole_gene_gc") if st.get("mean_whole_gene_gc") is not None else (st.get("mean_exonic_gc") or mu_cds)
            sigma_whole = st.get("std_whole_gene_gc") if st.get("std_whole_gene_gc") is not None else (st.get("std_exonic_gc") or sigma_cds)

            # Whole-gene GC is standardized against the empirical whole-gene baseline
            if full_gene_gc is not None and mu_whole is not None and sigma_whole and sigma_whole > 0:
                z_full_gene_gc = round(float((full_gene_gc - mu_whole) / sigma_whole), 2)
                delta_full_gene_gc = round(float(full_gene_gc - mu_whole), 2)

            pct_full_gene = None
            ref_col = (
                "whole_gene_gc"
                if ("whole_gene_gc" in df_sp.columns and df_sp["whole_gene_gc"].notna().sum() > 0)
                else ("exonic_gc" if ("exonic_gc" in df_sp.columns and df_sp["exonic_gc"].notna().sum() > 0) else "cds_gc")
            )
            if not df_sp.empty and ref_col in df_sp.columns and full_gene_gc is not None:
                pct_full_gene = round(float(stats.percentileofscore(df_sp[ref_col].dropna(), full_gene_gc)), 2)

            # Coding-span GC (excluding UTRs) and gene length
            coding_span_gc = gc_rec.get("coding_span_gc") if gc_rec.get("coding_span_gc") is not None else full_gene_gc
            coding_span_len = gc_rec.get("coding_span_len") if gc_rec.get("coding_span_len") else full_gene_len
            log10_len = round(float(np.log10(coding_span_len)), 4) if coding_span_len > 0 else None

            # Calculate chromosome length and distance to nearest telomere
            chrom = str(g.get("chromosome", ""))
            g_s = int(gc_rec.get("start", g.get("start", 0)))
            g_e = int(gc_rec.get("end", g.get("end", 0)))
            chr_key = f"{sp}::{chrom}"
            chr_len = getattr(self, "chr_lengths", {}).get(chr_key)
            dist_to_telo = None
            rel_dist_to_telo = None
            if chr_len and chr_len > 0:
                dist_to_telo = min(g_s, chr_len - g_e)
                rel_dist_to_telo = round(float(dist_to_telo / chr_len), 4)

            # Gene Length baseline (log10 length across all protein-coding genes)
            z_length = None
            mu_log_len = None
            sigma_log_len = None
            if not df_sp.empty and "whole_gene_len" in df_sp.columns:
                lens = df_sp["whole_gene_len"].dropna()
                valid_lens = lens[lens > 0]
                if not valid_lens.empty:
                    log_lens = np.log10(valid_lens)
                    mu_log_len = round(float(np.mean(log_lens)), 4)
                    sigma_log_len = round(float(np.std(log_lens)), 4)
                    if log10_len is not None and sigma_log_len > 0:
                        z_length = round(float((log10_len - mu_log_len) / sigma_log_len), 2)

            z_coding_span_gc = None
            if coding_span_gc is not None and mu_whole is not None and sigma_whole and sigma_whole > 0:
                z_coding_span_gc = round(float((coding_span_gc - mu_whole) / sigma_whole), 2)

            results.append({
                "species": sp,
                "clade": g.get("clade", "unknown"),
                "locus_id": locus_id,
                "locus_origin": locus_origin,
                "chromosome": chrom,
                "gene_id": g_id,
                "gene_name": g["gene_name"],
                "display_symbol": g.get("display_symbol", g["gene_name"]),
                "biotype": biotype,
                "comparative_group": comp_group,
                "is_globin": is_globin,
                "start": g_s,
                "end": g_e,
                "strand": g.get("strand", "+"),
                "coding_span_gc": coding_span_gc,
                "coding_span_len": coding_span_len,
                "log10_coding_span_len": log10_len,
                "full_gene_gc": full_gene_gc,
                "full_gene_len": full_gene_len,
                "exonic_gc": exonic_gc,
                "exonic_len": exonic_len,
                "intronic_gc": intronic_gc,
                "intronic_len": intronic_len,
                "cds_gc": cds_gc,
                "cds_len": cds_len,
                "gc3": gc3,
                "chrom_length_bp": chr_len,
                "dist_to_telomere_bp": dist_to_telo,
                "rel_dist_to_telomere": rel_dist_to_telo,
                "genome_mean_gene_gc": mu_whole,
                "genome_std_gene_gc": sigma_whole,
                "genome_mean_genomic_gc": mu_genomic,
                "genome_std_genomic_gc": sigma_genomic,
                "genome_mean_cds_gc": mu_cds,
                "genome_std_cds_gc": sigma_cds,
                "genome_mean_gc3": mu_gc3,
                "genome_std_gc3": sigma_gc3,
                "genome_mean_exonic_gc": mu_exon,
                "genome_std_exonic_gc": sigma_exon,
                "genome_mean_log10_len": mu_log_len,
                "genome_std_log10_len": sigma_log_len,
                "z_score_gc": z_cds_gc,
                "z_score_gc3": z_gc3,
                "z_score_exonic_gc": z_exonic_gc,
                "z_score_full_gene_gc": z_full_gene_gc,
                "z_score_coding_span_gc": z_coding_span_gc,
                "z_score_length": z_length,
                "delta_full_gene_gc": delta_full_gene_gc,
                "delta_cds_gc": delta_cds_gc,
                "percentile_gc": pct_gc,
                "percentile_gc3": pct_gc3,
                "percentile_exonic_gc": pct_exonic,
                "percentile_full_gene": pct_full_gene,
            })

        df_genes = pd.DataFrame(results)

        # Statistical hypothesis tests per species
        stat_rows = []
        for sp, grp in df_genes.groupby("species"):
            df_base = baselines_df.get(sp, pd.DataFrame())
            if df_base.empty:
                continue

            base_gc = df_base["cds_gc"].dropna().values
            base_gc3 = df_base["gc3"].dropna().values

            # 1. Globin coding vs Genome
            globins = grp[grp["comparative_group"].str.contains("globin_coding")]
            glob_gc = globins["cds_gc"].dropna().values
            glob_gc3 = globins["gc3"].dropna().values

            p_mw_gc = np.nan
            p_mw_gc3 = np.nan
            if len(glob_gc) > 0 and len(base_gc) > 0:
                _, p_mw_gc = stats.mannwhitneyu(glob_gc, base_gc, alternative="two-sided")
            if len(glob_gc3) > 0 and len(base_gc3) > 0:
                _, p_mw_gc3 = stats.mannwhitneyu(glob_gc3, base_gc3, alternative="two-sided")

            # 2. Alpha Globin vs Beta Globin (in amniotes)
            alpha_glob = grp[grp["comparative_group"] == "alpha_globin_coding"]["cds_gc"].dropna().values
            beta_glob = grp[grp["comparative_group"] == "beta_globin_coding"]["cds_gc"].dropna().values
            p_alpha_beta = np.nan
            diff_alpha_beta = np.nan
            mean_alpha_gc = float(np.mean(alpha_glob)) if len(alpha_glob) > 0 else np.nan
            mean_beta_gc = float(np.mean(beta_glob)) if len(beta_glob) > 0 else np.nan
            if len(alpha_glob) > 0 and len(beta_glob) > 0:
                diff_alpha_beta = round(mean_alpha_gc - mean_beta_gc, 2)
                try:
                    _, p_alpha_beta = stats.mannwhitneyu(alpha_glob, beta_glob, alternative="two-sided")
                except Exception:
                    pass

            # Whole-gene GC comparison for Alpha vs Beta
            alpha_fg = grp[grp["comparative_group"] == "alpha_globin_coding"]["full_gene_gc"].dropna().values
            beta_fg = grp[grp["comparative_group"] == "beta_globin_coding"]["full_gene_gc"].dropna().values
            p_alpha_beta_fg = np.nan
            diff_alpha_beta_fg = np.nan
            mean_alpha_fg = float(np.mean(alpha_fg)) if len(alpha_fg) > 0 else np.nan
            mean_beta_fg = float(np.mean(beta_fg)) if len(beta_fg) > 0 else np.nan
            if len(alpha_fg) > 0 and len(beta_fg) > 0:
                diff_alpha_beta_fg = round(mean_alpha_fg - mean_beta_fg, 2)
                try:
                    _, p_alpha_beta_fg = stats.mannwhitneyu(alpha_fg, beta_fg, alternative="two-sided")
                except Exception:
                    pass

            # 3. Teleost Locus A vs Locus B
            loc_a = grp[grp["comparative_group"] == "globin_coding_locus_A"]["cds_gc"].dropna().values
            loc_b = grp[grp["comparative_group"] == "globin_coding_locus_B"]["cds_gc"].dropna().values
            p_locus_ab = np.nan
            if len(loc_a) > 0 and len(loc_b) > 0:
                try:
                    _, p_locus_ab = stats.mannwhitneyu(loc_a, loc_b, alternative="two-sided")
                except Exception:
                    pass

            # 4. Globin coding vs Globin pseudogenes
            col_pseudo = "full_gene_gc" if "full_gene_gc" in grp.columns else ("genomic_gc" if "genomic_gc" in grp.columns else "cds_gc")
            pseudos = grp[grp["comparative_group"] == "globin_pseudogene"][col_pseudo].dropna().values
            mean_pseudo_gc = float(np.mean(pseudos)) if len(pseudos) > 0 else np.nan
            mean_glob_gc = float(np.mean(glob_gc)) if len(glob_gc) > 0 else np.nan
            diff_glob_pseudo = mean_glob_gc - mean_pseudo_gc if not np.isnan(mean_pseudo_gc) and not np.isnan(mean_glob_gc) else np.nan

            # 5. Cluster Non-globins
            ng_mask = grp["comparative_group"].str.startswith("cluster_non_globin")
            non_globs = grp[ng_mask]["cds_gc"].dropna().values
            mean_nonglob_gc = float(np.mean(non_globs)) if len(non_globs) > 0 else np.nan
            z_nonglob = grp[ng_mask]["z_score_gc"].dropna().values
            mean_z_nonglob = float(np.mean(z_nonglob)) if len(z_nonglob) > 0 else np.nan

            # Locus-specific non-globins (Alpha vs Beta locus)
            ng_alpha_fg = grp[grp["comparative_group"] == "cluster_non_globin_alpha"]["full_gene_gc"].dropna().values
            ng_beta_fg = grp[grp["comparative_group"] == "cluster_non_globin_beta"]["full_gene_gc"].dropna().values
            mean_ng_alpha_fg = float(np.mean(ng_alpha_fg)) if len(ng_alpha_fg) > 0 else np.nan
            mean_ng_beta_fg = float(np.mean(ng_beta_fg)) if len(ng_beta_fg) > 0 else np.nan
            diff_ng_alpha_beta_fg = round(mean_ng_alpha_fg - mean_ng_beta_fg, 2) if not np.isnan(mean_ng_alpha_fg) and not np.isnan(mean_ng_beta_fg) else np.nan

            stat_rows.append({
                "species": sp,
                "clade": grp["clade"].iloc[0],
                "num_globins": len(glob_gc),
                "mean_globin_gc": round(float(np.mean(glob_gc)), 2) if len(glob_gc) > 0 else np.nan,
                "mean_globin_z_gc": round(float(globins["z_score_gc"].mean()), 2) if not globins.empty else np.nan,
                "mean_globin_gc3": round(float(np.mean(glob_gc3)), 2) if len(glob_gc3) > 0 else np.nan,
                "mean_globin_z_gc3": round(float(globins["z_score_gc3"].mean()), 2) if not globins.empty else np.nan,
                "p_globin_vs_genome_gc": p_mw_gc,
                "p_globin_vs_genome_gc3": p_mw_gc3,
                "num_alpha_globins": len(alpha_glob),
                "mean_alpha_globin_gc": round(mean_alpha_gc, 2) if not np.isnan(mean_alpha_gc) else np.nan,
                "num_beta_globins": len(beta_glob),
                "mean_beta_globin_gc": round(mean_beta_gc, 2) if not np.isnan(mean_beta_gc) else np.nan,
                "diff_alpha_minus_beta_gc": diff_alpha_beta,
                "p_alpha_vs_beta_gc": p_alpha_beta,
                "mean_alpha_full_gene_gc": round(mean_alpha_fg, 2) if not np.isnan(mean_alpha_fg) else np.nan,
                "mean_beta_full_gene_gc": round(mean_beta_fg, 2) if not np.isnan(mean_beta_fg) else np.nan,
                "diff_alpha_minus_beta_full_gene_gc": diff_alpha_beta_fg,
                "p_alpha_vs_beta_full_gene_gc": p_alpha_beta_fg,
                "num_pseudogenes": len(pseudos),
                "mean_pseudogene_gc": round(mean_pseudo_gc, 2) if not np.isnan(mean_pseudo_gc) else np.nan,
                "diff_coding_minus_pseudo_gc": round(diff_glob_pseudo, 2) if not np.isnan(diff_glob_pseudo) else np.nan,
                "num_cluster_non_globins": len(non_globs),
                "mean_cluster_non_globin_gc": round(mean_nonglob_gc, 2) if not np.isnan(mean_nonglob_gc) else np.nan,
                "mean_cluster_non_globin_z_gc": round(mean_z_nonglob, 2) if not np.isnan(mean_z_nonglob) else np.nan,
                "mean_alpha_nonglobin_full_gene_gc": round(mean_ng_alpha_fg, 2) if not np.isnan(mean_ng_alpha_fg) else np.nan,
                "mean_beta_nonglobin_full_gene_gc": round(mean_ng_beta_fg, 2) if not np.isnan(mean_ng_beta_fg) else np.nan,
                "diff_alpha_minus_beta_nonglobin_full_gene_gc": diff_ng_alpha_beta_fg,
                "p_locus_a_vs_b_gc": p_locus_ab,
            })

        df_stats = pd.DataFrame(stat_rows)
        return df_genes, df_stats
