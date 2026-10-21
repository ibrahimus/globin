#!/usr/bin/env python3
"""
scripts/update_species_artifacts.py

Universal Per-Species Artifact & Dependency Synchronization Pipeline.

This script enforces the complete dependency chain for any modified or updated species:
1. Generates ALL FOUR per-species visualization suites:
   a. Linear Chromosome Track (chromosome_plots/chromosome_linear_{sp}.*) [PNG, PDF, HTML]
   b. Bivariate Length vs GC (bivariate_length_gc_plots/bivariate_length_gc_{sp}.*) [PNG, PDF, HTML]
   c. Bivariate with Chromosome Inset (bivariate_with_inset_plots/bivariate_with_inset_{sp}.*) [PNG, PDF, HTML]
   d. Diagnostic GC Distribution (diagnostic_plots/diagnostic_{sp}.*) [PNG, HTML]
2. Synchronizes all generated artifacts into 1_data_from_mike/results/.
3. Provides smart stale-checking (--check-stale): detects species whose plot files are
   missing or older than results/gene_zscores.csv or results/curated_intervals.csv.
4. Provides full cascade mode (--cascade): updates downstream summary tables,
   Figures 0-4, and results/dashboard.html.

Usage:
  python scripts/update_species_artifacts.py mirounga_angustirostris
  python scripts/update_species_artifacts.py --check-stale
  python scripts/update_species_artifacts.py --all
  python scripts/update_species_artifacts.py mirounga_angustirostris --cascade
"""

import os
import sys
import time
import shutil
import argparse
import subprocess
import pandas as pd
from typing import List, Set
from globin_pipeline.visualizer import Visualizer

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ZSCORES_FILE = os.path.join(RESULTS_DIR, "gene_zscores.csv")
INTERVALS_FILE = os.path.join(RESULTS_DIR, "curated_intervals.csv")


def get_stale_species(df_genes: pd.DataFrame) -> List[str]:
    """Identifies species whose per-species plots are missing or older than the master tables."""
    if not os.path.exists(ZSCORES_FILE):
        return []
    
    table_mtime = max(
        os.path.getmtime(ZSCORES_FILE),
        os.path.getmtime(INTERVALS_FILE) if os.path.exists(INTERVALS_FILE) else 0.0
    )

    stale: Set[str] = set()
    all_species = sorted(df_genes["species"].unique())

    required_artifacts = [
        ("chromosome_plots", "chromosome_linear_{sp}.png"),
        ("bivariate_length_gc_plots", "bivariate_length_gc_{sp}.png"),
        ("bivariate_with_inset_plots", "bivariate_with_inset_{sp}.png"),
        ("diagnostic_plots", "diagnostic_{sp}.png"),
    ]

    for sp in all_species:
        for folder, pattern in required_artifacts:
            path = os.path.join(RESULTS_DIR, folder, pattern.format(sp=sp))
            if not os.path.exists(path) or os.path.getmtime(path) < table_mtime:
                stale.add(sp)
                break

    return sorted(list(stale))


def update_species_plots(species_list: List[str], viz: Visualizer, df_genes: pd.DataFrame, df_intervals: pd.DataFrame):
    """Renders all 4 plot families for each specified species."""
    chr_out_dir = os.path.join(RESULTS_DIR, "chromosome_plots")
    inset_out_dir = os.path.join(RESULTS_DIR, "bivariate_with_inset_plots")
    biv_out_dir = os.path.join(RESULTS_DIR, "bivariate_length_gc_plots")
    diag_out_dir = os.path.join(RESULTS_DIR, "diagnostic_plots")

    os.makedirs(chr_out_dir, exist_ok=True)
    os.makedirs(inset_out_dir, exist_ok=True)
    os.makedirs(biv_out_dir, exist_ok=True)
    os.makedirs(diag_out_dir, exist_ok=True)

    total = len(species_list)
    print(f"\n>> Updating 4-Suite Artifacts for {total} species...")

    for i, sp in enumerate(species_list, 1):
        t0 = time.time()
        print(f"[{i:02d}/{total}] Rendering all suites for: {sp}...", end="", flush=True)

        # 1. Chromosome linear track (PNG, PDF, HTML)
        viz.plot_chromosome_linear(sp, df_genes, df_intervals, output_dir=chr_out_dir)

        # 2. Standard Bivariate Length vs GC (PNG, PDF, HTML)
        viz.plot_length_vs_gc(sp, df_genes, output_dir=biv_out_dir)

        # 3. Bivariate with Chromosome Inset (PNG, PDF, HTML)
        viz.plot_bivariate_with_chromosome_inset(sp, df_genes, df_intervals, output_dir=inset_out_dir)

        # 4. Diagnostic GC Distribution (6-panel PNG, interactive Plotly HTML)
        viz.plot_diagnostic(sp, df_genes, output_dir=diag_out_dir)
        viz.plot_diagnostic_interactive(sp, df_genes, output_dir=diag_out_dir)

        dt = time.time() - t0
        print(f" done in {dt:.1f}s")

    # Synchronize all four suites to 1_data_from_mike/results/
    sync_targets = [
        (chr_out_dir, os.path.join(BASE_DIR, "1_data_from_mike/results/chromosome_plots")),
        (inset_out_dir, os.path.join(BASE_DIR, "1_data_from_mike/results/bivariate_with_inset_plots")),
        (biv_out_dir, os.path.join(BASE_DIR, "1_data_from_mike/results/bivariate_length_gc_plots")),
        (diag_out_dir, os.path.join(BASE_DIR, "1_data_from_mike/results/diagnostic_plots")),
    ]

    print("\n>> Synchronizing artifact suites to 1_data_from_mike/results/...")
    for src, dst in sync_targets:
        os.makedirs(dst, exist_ok=True)
        shutil.copytree(src, dst, dirs_exist_ok=True)
    print(">> Synchronization complete.")


def run_downstream_cascade():
    """Executes the macroevolutionary summary and dashboard cascade."""
    cascade_steps = [
        ("Macroevolutionary Summary Tables", [sys.executable, "scripts/build_evolutionary_history_and_synteny.py"]),
        ("Overview & Guide Figures (Fig 0)", [sys.executable, "scripts/generate_overview_methodology_figure.py"]),
        ("Macroevolutionary Figures (Fig 1, 2)", [sys.executable, "scripts/generate_all_summary_figures.py"]),
        ("Synteny & Codon Heatmaps (Fig 3)", [sys.executable, "scripts/generate_figure3_phylogenetic_synteny.py"]),
        ("Bernardi Compaction Law (Fig 4)", [sys.executable, "scripts/analyze_bernardi_compaction.py"]),
        ("Static HTML Dashboard", [sys.executable, "scripts/generate_static_dashboard.py"]),
    ]

    print("\n" + "=" * 75)
    print(">> EXECUTING DOWNSTREAM DEPENDENCY CASCADE")
    print("=" * 75)

    for step_name, cmd in cascade_steps:
        print(f"\n>> [CASCADE] {step_name}...")
        t0 = time.time()
        res = subprocess.run(cmd, cwd=BASE_DIR)
        if res.returncode != 0:
            print(f"!! [ERROR] Cascade step failed: {step_name}")
            sys.exit(res.returncode)
        print(f">> [SUCCESS] Completed {step_name} in {time.time() - t0:.1f}s")


def main():
    parser = argparse.ArgumentParser(description="Universal Per-Species Artifact & Dependency Synchronization")
    parser.add_argument("species", nargs="*", help="Specific species name(s) to update")
    parser.add_argument("--check-stale", action="store_true", help="Automatically detect and update stale species")
    parser.add_argument("--all", action="store_true", help="Update all species in dataset")
    parser.add_argument("--cascade", action="store_true", help="Cascade downstream to summary figures and dashboard")
    args = parser.parse_args()

    print("=" * 75)
    print("UNIVERSAL SPECIES ARTIFACT & DEPENDENCY MANAGER")
    print("=" * 75)

    df_genes = pd.read_csv(ZSCORES_FILE)
    df_intervals = pd.read_csv(INTERVALS_FILE)
    viz = Visualizer()

    targets: List[str] = []

    if args.all:
        targets = sorted(df_genes["species"].unique())
        print(f"Targeting ALL {len(targets)} species in dataset.")
    elif args.check_stale:
        targets = get_stale_species(df_genes)
        if not targets:
            print("All species artifact suites (chromosome, bivariate, inset, diagnostic) are up to date!")
            if not args.cascade:
                return
        else:
            print(f"Detected {len(targets)} stale species: {', '.join(targets)}")
    elif args.species:
        targets = [sp for sp in args.species if sp in df_genes["species"].unique()]
        unknown = set(args.species) - set(targets)
        if unknown:
            print(f"Warning: Unknown species ignored: {unknown}")
    else:
        # Default behavior if run with no args: check stale
        targets = get_stale_species(df_genes)
        if not targets:
            print("All species artifact suites are currently up to date.")
            print("Run with a species name (e.g. mirounga_angustirostris) or --all to force regeneration.")
            return
        print(f"Auto-detected {len(targets)} stale species: {', '.join(targets)}")

    if targets:
        update_species_plots(targets, viz, df_genes, df_intervals)

    if args.cascade or args.species:
        run_downstream_cascade()

    print("\n" + "=" * 75)
    print("ALL REQUESTED ARTIFACTS AND DEPENDENCIES SUCCESSFULLY SYNCHRONIZED!")
    print("=" * 75)


if __name__ == "__main__":
    main()
