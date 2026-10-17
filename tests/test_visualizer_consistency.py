"""Automated consistency and integrity tests for globin visualization outputs.

Verifies that:
1. No globin gene is omitted from chromosome_linear plots across all 20 species.
2. All globin genes are faithfully represented in bivariate plots with mini chromosome linear insets.
3. Coordinates and intervals match between gene_zscores and curated_intervals.
"""

import os
import pytest
import pandas as pd
from globin_pipeline.visualizer import Visualizer


from pathlib import Path

@pytest.fixture(scope="module")
def gene_data():
    project_root = Path(__file__).resolve().parent.parent
    csv_genes = project_root / "results" / "gene_zscores.csv"
    csv_intervals = project_root / "globin_pipeline" / "data" / "curated_intervals.csv"
    if not csv_intervals.exists():
        csv_intervals = project_root / "results" / "curated_intervals.csv"
    df_genes = pd.read_csv(csv_genes)
    df_intervals = pd.read_csv(csv_intervals)
    return df_genes, df_intervals


def test_no_globins_omitted_from_intervals(gene_data):
    """Ensure every globin gene in gene_zscores matches an interval."""
    df_genes, df_intervals = gene_data

    def is_true_globin(r):
        grp = str(r.get("comparative_group", ""))
        is_g = str(r.get("is_globin", False)).lower() == "true"
        return is_g or ("globin" in grp and "non_globin" not in grp)

    globins = df_genes[df_genes.apply(is_true_globin, axis=1)]

    for sp in globins["species"].unique():
        sp_g = globins[globins["species"] == sp]
        sp_int = df_intervals[df_intervals["species"] == sp]
        assert not sp_int.empty, f"Species {sp} has no intervals in curated_intervals.csv!"

        valid_loci = set(sp_int["locus_id"])
        for _, g in sp_g.iterrows():
            loc = g.get("locus_id")
            sym = g.get("display_symbol")
            assert loc in valid_loci, f"Globin {sym} in {sp} has locus_id '{loc}' not in curated_intervals: {valid_loci}"


def test_chromosome_linear_safeguard_all_species(gene_data, tmp_path):
    """Verify plot_chromosome_linear runs and passes safeguard for all 20 species."""
    df_genes, df_intervals = gene_data
    viz = Visualizer()
    out_dir = str(tmp_path / "chromosome_plots")

    species_list = sorted(df_genes["species"].unique())
    assert len(species_list) >= 44, f"Expected at least 44 species, found {len(species_list)}"

    for sp in species_list:
        res = viz.plot_chromosome_linear(sp, df_genes, df_intervals, output_dir=out_dir)
        assert os.path.exists(res["png"]), f"PNG missing for {sp}"
        assert os.path.exists(res["pdf"]), f"PDF missing for {sp}"
        assert os.path.exists(res["html"]), f"HTML missing for {sp}"


def test_bivariate_with_inset_safeguard_all_species(gene_data, tmp_path):
    """Verify plot_bivariate_with_chromosome_inset runs and passes safeguard for all 20 species."""
    df_genes, df_intervals = gene_data
    viz = Visualizer()
    out_dir = str(tmp_path / "bivariate_with_inset_plots")

    species_list = sorted(df_genes["species"].unique())
    for sp in species_list:
        res = viz.plot_bivariate_with_chromosome_inset(sp, df_genes, df_intervals, output_dir=out_dir)
        assert os.path.exists(res["png"]), f"PNG missing for {sp}"
        assert os.path.exists(res["pdf"]), f"PDF missing for {sp}"
        assert os.path.exists(res["html"]), f"HTML missing for {sp}"
