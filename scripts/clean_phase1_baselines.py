#!/usr/bin/env python3
"""
Clean Phase 1 Expansion Genome Baselines:
Filter out non-coding RNA features (tRNAs, 5S rRNAs, snoRNAs, etc.) that lack CDS,
retaining strictly protein-coding gene models with valid CDS.
Recalculate summary baseline statistics (mean, std, median) for each species.
"""

import os
import json
import numpy as np
import pandas as pd

BASELINE_DIR = "data/genome_baselines"

PHASE1_SPECIES = [
    "anser_cygnoides",
    "coturnix_japonica",
    "equus_asinus",
    "equus_caballus",
    "eretmochelys_imbricata",
    "geotrypetes_seraphini",
    "microcaecilia_unicolor",
    "miniopterus_natalensis",
    "pimephales_promelas",
    "puma_concolor",
    "rattus_norvegicus",
    "rhinatrema_bivittatum",
    "sturnus_vulgaris",
]


def clean_baselines():
    print(f"Cleaning non-coding features from {len(PHASE1_SPECIES)} species baselines...\n")

    for sp in PHASE1_SPECIES:
        csv_path = os.path.join(BASELINE_DIR, f"{sp}_genes_gc.csv.gz")
        json_path = os.path.join(BASELINE_DIR, f"{sp}_stats.json")

        if not os.path.exists(csv_path):
            print(f"[{sp}] CSV not found: {csv_path}, skipping.")
            continue

        df = pd.read_csv(csv_path)
        n_orig = len(df)

        # Filter strictly for genuine protein-coding genes with valid CDS
        # In RefSeq GFF, non-coding RNAs lack CDS records
        df_clean = df[df["cds_length"].notna() & (df["cds_length"] > 0)].copy()
        df_clean["biotype"] = "protein_coding"
        n_clean = len(df_clean)

        # Write cleaned CSV back to disk
        df_clean.to_csv(csv_path, index=False, compression="gzip")

        # Load existing stats to preserve genomic_gc
        genomic_gc = 42.0
        if os.path.exists(json_path):
            with open(json_path, "r", encoding="utf-8") as f:
                old_stats = json.load(f)
                genomic_gc = old_stats.get("genomic_gc", 42.0)

        # Recalculate baseline statistics from clean protein-coding genes
        w_lens = df_clean["whole_gene_len"].dropna().values
        w_gcs = df_clean["whole_gene_gc"].dropna().values
        c_gcs = df_clean["cds_gc"].dropna().values
        gc3_vals = df_clean["gc3"].dropna().values if "gc3" in df_clean.columns else np.array([])

        stats = {
            "species": sp,
            "total_genes": n_clean,
            "genomic_gc": float(genomic_gc),
            "mean_whole_gene_len": float(np.mean(w_lens)) if len(w_lens) > 0 else 0.0,
            "std_whole_gene_len": float(np.std(w_lens, ddof=1)) if len(w_lens) > 1 else 0.0,
            "median_whole_gene_len": float(np.median(w_lens)) if len(w_lens) > 0 else 0.0,
            "mean_whole_gene_gc": float(np.mean(w_gcs)) if len(w_gcs) > 0 else 0.0,
            "std_whole_gene_gc": float(np.std(w_gcs, ddof=1)) if len(w_gcs) > 1 else 0.0,
            "median_whole_gene_gc": float(np.median(w_gcs)) if len(w_gcs) > 0 else 0.0,
            "mean_cds_gc": float(np.mean(c_gcs)) if len(c_gcs) > 0 else 0.0,
            "std_cds_gc": float(np.std(c_gcs, ddof=1)) if len(c_gcs) > 1 else 0.0,
            "median_cds_gc": float(np.median(c_gcs)) if len(c_gcs) > 0 else 0.0,
            "mean_gc3": float(np.mean(gc3_vals)) if len(gc3_vals) > 0 else None,
            "std_gc3": float(np.std(gc3_vals, ddof=1)) if len(gc3_vals) > 1 else None,
            "median_gc3": float(np.median(gc3_vals)) if len(gc3_vals) > 0 else None,
        }

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)

        print(f"[{sp:25s}] Cleaned: {n_orig:5d} -> {n_clean:5d} genes | Mean Len: {stats['mean_whole_gene_len']:7.0f} bp | Mean Whole GC: {stats['mean_whole_gene_gc']:5.2f}% (std: {stats['std_whole_gene_gc']:4.2f}%) | Mean CDS GC: {stats['mean_cds_gc']:5.2f}%")

    print("\nAll Phase 1 species baselines cleaned and stats recalculated successfully!")


if __name__ == "__main__":
    clean_baselines()
