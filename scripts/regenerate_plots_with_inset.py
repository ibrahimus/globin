"""Batch script to regenerate all chromosome_plots and generate bivariate_with_inset_plots for all 20 species,
then synchronize results to 1_data_from_mike/results/.
"""

import os
import sys
import shutil
import pandas as pd
from globin_pipeline.visualizer import Visualizer

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    print("Initializing Visualizer...")
    viz = Visualizer()
    results_dir = os.path.join(BASE_DIR, "results")
    df_genes = pd.read_csv(os.path.join(results_dir, "gene_zscores.csv"))
    df_intervals = pd.read_csv(os.path.join(results_dir, "curated_intervals.csv"))

    species_arg = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if species_arg:
        species_list = [sp for sp in species_arg if sp in df_genes["species"].unique()]
    else:
        species_list = sorted(df_genes["species"].unique())
    print(f"Loaded {len(species_list)} species.")

    chr_out_dir = os.path.join(results_dir, "chromosome_plots")
    inset_out_dir = os.path.join(results_dir, "bivariate_with_inset_plots")
    biv_out_dir = os.path.join(results_dir, "bivariate_length_gc_plots")
    diag_out_dir = os.path.join(results_dir, "diagnostic_plots")
    os.makedirs(chr_out_dir, exist_ok=True)
    os.makedirs(inset_out_dir, exist_ok=True)
    os.makedirs(biv_out_dir, exist_ok=True)
    os.makedirs(diag_out_dir, exist_ok=True)

    print("\n--- Generating chromosome_plots, bivariate_length_gc_plots, bivariate_with_inset_plots, and diagnostic_plots ---")
    for i, sp in enumerate(species_list, 1):
        print(f"[{i:02d}/{len(species_list)}] Processing species: {sp}...")
        
        # 1. Regenerate chromosome_linear (validates 100% globin rendering)
        chr_res = viz.plot_chromosome_linear(sp, df_genes, df_intervals, output_dir=chr_out_dir)
        
        # 2. Generate standard bivariate_length_gc (without other features / pseudogenes)
        biv_res = viz.plot_length_vs_gc(sp, df_genes, output_dir=biv_out_dir)

        # 3. Generate bivariate_with_chromosome_inset (without other features / pseudogenes)
        inset_res = viz.plot_bivariate_with_chromosome_inset(sp, df_genes, df_intervals, output_dir=inset_out_dir)

        # 4. Generate Diagnostic GC distribution (static PNG and interactive Plotly HTML)
        viz.plot_diagnostic(sp, df_genes, output_dir=diag_out_dir)
        viz.plot_diagnostic_interactive(sp, df_genes, output_dir=diag_out_dir)

    print(f"\nAll {len(species_list)} species successfully processed!")

    # Synchronize to 1_data_from_mike/results/
    target_chr = os.path.join(BASE_DIR, "1_data_from_mike/results/chromosome_plots")
    target_inset = os.path.join(BASE_DIR, "1_data_from_mike/results/bivariate_with_inset_plots")
    target_biv = os.path.join(BASE_DIR, "1_data_from_mike/results/bivariate_length_gc_plots")
    target_diag = os.path.join(BASE_DIR, "1_data_from_mike/results/diagnostic_plots")
    
    print(f"Synchronizing {chr_out_dir} -> {target_chr}...")
    shutil.copytree(chr_out_dir, target_chr, dirs_exist_ok=True)
    
    print(f"Synchronizing {inset_out_dir} -> {target_inset}...")
    shutil.copytree(inset_out_dir, target_inset, dirs_exist_ok=True)

    print(f"Synchronizing {biv_out_dir} -> {target_biv}...")
    shutil.copytree(biv_out_dir, target_biv, dirs_exist_ok=True)

    print(f"Synchronizing {diag_out_dir} -> {target_diag}...")
    shutil.copytree(diag_out_dir, target_diag, dirs_exist_ok=True)

    print("Synchronization complete!")

if __name__ == "__main__":
    main()
