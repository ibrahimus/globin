"""Unit tests for advanced baseline, zscore, and visualizer modules."""

import os
import pytest
import numpy as np
import pandas as pd
from globin_pipeline.baseline_manager import BaselineManager
from globin_pipeline.zscore_analyzer import ZScoreAnalyzer
from globin_pipeline.visualizer import Visualizer


def test_baseline_manager_calculations():
    bm = BaselineManager(cache_dir="data/test_baselines")

    seq = "ATGGCCGCT"  # Codons: ATG (G), GCC (C), GCT (T) -> len=9
    # GC = (A:1, T:2, G:3, C:3) -> 6/9 = 66.67%
    assert bm.calculate_gc(seq) == 66.67

    # GC3: 3rd bases are G, C, T -> 2/3 = 66.67%
    assert bm.calculate_gc3(seq) == 66.67


def test_zscore_analyzer_logic(tmp_path):
    bm = BaselineManager(cache_dir=str(tmp_path))

    # Mock baseline data
    mock_genes = pd.DataFrame([
        {"species": "Test_sp", "gene_id": f"GENE{i}", "cds_gc": 50.0 + i, "gc3": 55.0 + i}
        for i in range(100)
    ])
    mock_stats = {
        "species": "Test_sp",
        "mean_cds_gc": 50.0,
        "std_cds_gc": 10.0,
        "mean_gc3": 55.0,
        "std_gc3": 15.0,
    }
    # Save mock cache
    t_path, s_path = bm._get_cache_paths("Test_sp")
    mock_genes.to_csv(t_path, index=False, compression="gzip")
    import json
    with open(s_path, "w") as f:
        json.dump(mock_stats, f)

    za = ZScoreAnalyzer(bm)

    locus_genes = [
        {
            "species": "Test_sp",
            "clade": "amniote",
            "locus_id": "alpha_globin",
            "chromosome": "1",
            "gene_id": "GLOB1",
            "gene_name": "HBA1",
            "display_symbol": "HBA1",
            "biotype": "protein_coding",
            "category": "alpha_globin",
            "subtype": "HBA",
            "is_globin": "True",
        },
        {
            "species": "Test_sp",
            "clade": "amniote",
            "locus_id": "alpha_globin",
            "chromosome": "1",
            "gene_id": "PSEUDO1",
            "gene_name": "HBAP1",
            "display_symbol": "HBAP1",
            "biotype": "unprocessed_pseudogene",
            "category": "globin_pseudogene",
            "subtype": "globin_pseudogene",
            "is_globin": "Pseudogene",
        },
    ]

    gene_gc = [
        {"species": "Test_sp", "gene_id": "GLOB1", "cds_gc": 70.0, "gc3": 85.0, "genomic_gc": 68.0},
        {"species": "Test_sp", "gene_id": "PSEUDO1", "cds_gc": None, "gc3": None, "genomic_gc": 45.0},
    ]

    df_res, df_stats = za.analyze_genes(locus_genes, gene_gc)

    # Check Z-score: GLOB1 has CDS GC 70.0 -> (70 - 50) / 10 = +2.0
    glob_row = df_res[df_res["gene_id"] == "GLOB1"].iloc[0]
    assert glob_row["z_score_gc"] == 2.0
    # GC3: 85.0 -> (85 - 55) / 15 = +2.0
    assert glob_row["z_score_gc3"] == 2.0
    assert glob_row["comparative_group"] == "alpha_globin_coding"

    # PSEUDO1 should be in globin_pseudogene
    pseudo_row = df_res[df_res["gene_id"] == "PSEUDO1"].iloc[0]
    assert pseudo_row["comparative_group"] == "globin_pseudogene"


def test_visualizer_plotting(tmp_path):
    viz = Visualizer()
    df_mock = pd.DataFrame([
        {
            "species": "Homo_sapiens",
            "clade": "amniote_mammal",
            "gene_id": "ENSG01",
            "display_symbol": "HBA1",
            "comparative_group": "globin_coding_locus_A",
            "z_score_gc": 1.8,
            "cds_gc": 65.0,
            "gc3": 88.0,
            "genomic_gc": 64.0,
            "full_gene_gc": 64.0,
            "exonic_gc": 65.0,
            "intronic_gc": 62.0,
            "percentile_gc": 95.0,
            "gene_name": "HBA1",
        },
        {
            "species": "danio_rerio",
            "clade": "teleost_fish",
            "gene_id": "ENSD01",
            "display_symbol": "hbaa1",
            "comparative_group": "globin_coding_locus_A",
            "z_score_gc": 0.8,
            "cds_gc": 52.0,
            "gc3": 62.0,
            "genomic_gc": 50.0,
            "full_gene_gc": 50.0,
            "exonic_gc": 52.0,
            "intronic_gc": 48.0,
            "percentile_gc": 60.0,
            "gene_name": "hbaa1",
        },
    ])

    meta_out = str(tmp_path / "test_meta.png")
    out = viz.plot_meta_comparative(df_mock, output_path=meta_out, clean_mode=False)
    assert os.path.exists(out)
    assert os.path.exists(meta_out.replace(".png", ".pdf"))

    meta_clean_out = str(tmp_path / "test_meta_clean.png")
    clean_out = viz.plot_meta_comparative(df_mock, output_path=meta_clean_out, clean_mode=True)
    assert os.path.exists(clean_out)

    meta_html = str(tmp_path / "test_meta.html")
    html_out = viz.plot_meta_interactive(df_mock, output_path=meta_html, clean_mode=False)
    assert os.path.exists(html_out)

    meta_clean_html = str(tmp_path / "test_meta_clean.html")
    clean_html_out = viz.plot_meta_interactive(df_mock, output_path=meta_clean_html, clean_mode=True)
    assert os.path.exists(clean_html_out)

    # Test deterministic y offset calculation
    y_map = viz._calculate_meta_y_positions(df_mock, "z_score_gc", 0.0, clean_mode=True)
    assert "ENSG01" in y_map
    assert "ENSD01" in y_map
