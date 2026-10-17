import os
import time
import pandas as pd
from globin_pipeline.visualizer import Visualizer

def main():
    t0 = time.time()
    zscores_path = "results/gene_zscores.csv"
    if not os.path.exists(zscores_path):
        zscores_path = "1_data_from_mike/results/gene_zscores.csv"
    
    df_genes = pd.read_csv(zscores_path)
    viz = Visualizer()
    
    chr_dir = "results/chromosome_plots"
    biv_dir = "results/bivariate_length_gc_plots"
    os.makedirs(chr_dir, exist_ok=True)
    os.makedirs(biv_dir, exist_ok=True)
    
    species_list = sorted(df_genes["species"].unique())
    print(f"Batch generating chromosome linear and bivariate length vs GC plots for {len(species_list)} species...")
    
    for idx, sp in enumerate(species_list, 1):
        sp_t0 = time.time()
        chr_res = viz.plot_chromosome_linear(sp, df_genes, output_dir=chr_dir)
        biv_res = viz.plot_length_vs_gc(sp, df_genes, output_dir=biv_dir)
        elapsed = time.time() - sp_t0
        print(f"[{idx:02d}/{len(species_list):02d}] {sp} ({elapsed:.1f}s)")
        if chr_res:
            print(f"   Chr: {chr_res.get('png')}")
        if biv_res:
            print(f"   Biv: {biv_res.get('png')}")

    print(f"\nAll plots successfully generated in {time.time() - t0:.1f}s.")

if __name__ == "__main__":
    main()
