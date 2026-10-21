#!/usr/bin/env python3
"""
reproduce_all.py

Master One-Command Reproduction Pipeline for the Vertebrate Globin Synteny & Isochore Dynamics Project.

This script executes the entire analytical and figure generation workflow from curated tables:
1. Audits cluster completeness across all 68 vertebrate genomes and 121 loci.
2. Generates Overview & Methodology Figures (Figure 0, 0A, 0B, 0C, 0D).
3. Generates Macroevolutionary Summary Figures (Figure 1A, 1B, Figure 2, Figure 3).
4. Generates Phylogenetic Synteny Heatmap & Architecture (Figure 3/4).
5. Runs the Empirical Proof of Bernardi's Isochore Compaction Law (Figure 4).
6. Builds the Standalone HTML Dashboard (results/dashboard.html).
7. Executes the automated test suite (pytest).
"""

import sys
import os
import time
import subprocess

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def run_step(step_name, command_args):
    print(f"\n{'='*75}")
    print(f">> [STEP] {step_name}")
    print(f"   Command: {' '.join(command_args)}")
    print(f"{'='*75}")
    t0 = time.time()
    res = subprocess.run(command_args, cwd=ROOT_DIR)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"!! [ERROR] Step failed with return code {res.returncode}: {step_name}")
        sys.exit(res.returncode)
    print(f">> [SUCCESS] Completed {step_name} in {dt:.1f}s")

def main():
    print("=" * 75)
    print("VERTEBRATE GLOBIN REPRODUCIBILITY PIPELINE")
    print(f"Working Directory: {ROOT_DIR}")
    print(f"Python: {sys.executable}")
    print("=" * 75)
    
    total_t0 = time.time()

    # Step 1: Audit Cluster Completeness
    run_step("Cluster Completeness Audit", [
        sys.executable, "scripts/verify_cluster_completeness.py"
    ])

    # Step 2: Overview & Methodology Figures (Fig 0)
    run_step("Generate Overview & Methodology Figures (Fig 0A, 0B, 0C, 0D)", [
        sys.executable, "scripts/generate_overview_methodology_figure.py"
    ])
    run_step("Generate Interactive Overview Slide Navigator", [
        sys.executable, "scripts/generate_overview_guide_html.py"
    ])

    # Step 3: Macroevolutionary Divergence & Drivers (Fig 1, 2)
    run_step("Generate Macroevolutionary Figures (Fig 1A, 1B, Fig 2)", [
        sys.executable, "scripts/generate_all_summary_figures.py"
    ])

    # Step 4: Phylogenetic Synteny & Codon Divergence (Fig 3)
    run_step("Generate Phylogenetic Synteny & Codon Heatmaps (Fig 3)", [
        sys.executable, "scripts/generate_figure3_phylogenetic_synteny.py"
    ])

    # Step 5: Empirical Proof of Bernardi's Compaction Law (Fig 4)
    run_step("Analyze Bernardi's Isochore Compaction Law Across Vertebrates (Fig 4)", [
        sys.executable, "scripts/analyze_bernardi_compaction.py"
    ])

    # Step 6: Standalone Static HTML Dashboard
    run_step("Generate Standalone HTML Dashboard", [
        sys.executable, "scripts/generate_static_dashboard.py"
    ])

    # Step 7: Automated Test Suite
    run_step("Run Automated Test Suite (pytest)", [
        sys.executable, "-m", "pytest", "tests"
    ])

    total_dt = time.time() - total_t0
    print("\n" + "=" * 75)
    print(f"ALL STEPS COMPLETED SUCCESSFULLY in {total_dt:.1f}s ({total_dt/60:.2f} min)!")
    print("All summary figures, statistical models, and the dashboard are up to date.")
    print(f"Open the dashboard: file://{os.path.join(ROOT_DIR, 'results', 'dashboard.html')}")
    print("=" * 75)

if __name__ == "__main__":
    main()
