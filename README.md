# globin-pipeline: Comparative Vertebrate Globin Genomics

[![Python](https://img.shields.io/badge/python-%3E%3D3.10-brightgreen.svg)](file:///Users/ibrahim/Documents/globin_project/pyproject.toml)
[![Tests](https://img.shields.io/badge/tests-16%20passed-success.svg)](file:///Users/ibrahim/Documents/globin_project/tests/)
[![Genomes](https://img.shields.io/badge/vertebrate%20genomes-68%20curated-blue.svg)](file:///Users/ibrahim/Documents/globin_project/results/cross_species_summary_table.csv)
[![Loci](https://img.shields.io/badge/curated%20loci-121%20intervals-indigo.svg)](file:///Users/ibrahim/Documents/globin_project/results/curated_intervals.csv)
[![Genes](https://img.shields.io/badge/gene%20models-1%2C544%20verified-purple.svg)](file:///Users/ibrahim/Documents/globin_project/results/gene_zscores.csv)
[![Baselines](https://img.shields.io/badge/genome%20baselines-72%20species%20%28%3E1.4M%20genes%29-teal.svg)](file:///Users/ibrahim/Documents/globin_project/results/bernardi_compaction_genome_baselines.csv)

A standardized, fully reproducible comparative genomics platform and reference atlas for auditing, discovering, curating, and analyzing $\alpha$- and $\beta$-globin loci across vertebrate genomes spanning ~500 million years of evolutionary history.

---

## 🌐 Live Interactive Dashboard (Zero Installation Required)

You can explore the entire project directly in any web browser without installing Python or running any scripts:

👉 **[Launch Live Web Dashboard](https://your-username.github.io/your-repo-name/)** *(served statically via GitHub Pages from [`index.html`](file:///Users/ibrahim/Documents/globin_project/index.html) &rarr; [`results/dashboard.html`](file:///Users/ibrahim/Documents/globin_project/results/dashboard.html))*

### How to Enable GitHub Pages on Your Fork / Repo:
1. Go to your GitHub repository **Settings** &rarr; **Pages** (under "Code and automation").
2. Under **Build and deployment** &rarr; **Source**, select **Deploy from a branch**.
3. Choose Branch: **`main`** and Folder: **`/ (root)`**, then click **Save**.
4. In ~60 seconds, GitHub Pages will deploy your dashboard at `https://<your-username>.github.io/<your-repo-name>/`.

---

## Dataset Scale & Capabilities

- **Taxonomic Breadth**: **68 vertebrate genomes** across 11 evolutionary lineages (Cyclostomata, Chondrichthyes, Basal Actinopterygii, Teleostei, Sarcopterygii, Amphibia, Testudines, Squamata, Crocodylia, Aves, Monotremata, Metatheria, Eutheria).
- **Curated Syntenic Loci**: **121 verified genomic intervals** anchored by ultra-conserved orthologs (*RHBDF1*, *MPG*, *NPRL3*, *LUC7L*, *OR* clusters).
- **Annotated Repertoire**: **1,544 protein-coding genes** classified by structural role, GC compartment, developmental stage, and genomic position.
- **Empirical Baseline Normalization**: Complete distributions across **72 species baselines** totaling **>1,440,000 genes**, providing standardized $Z$-scores ($Z_{\text{gene}}$, $Z_{\text{cds}}$, $Z_{\text{gc3}}$) and genomic percentiles.
- **Empirical Proof of Bernardi's Compaction Law**: Quantitative proof that 98.6% of vertebrate genomes adhere to isochoric compaction, with invariant globin calibrator shrinking up to ~7x.

---

## ⚡ Quickstart & Local Reproduction

Anyone who clones this repository can reproduce the entire analytical workflow, verify cluster completeness, recompute all statistics, and regenerate all publication figures in **under 2 minutes** using the master reproduction script:

```bash
# 1. Clone the repository
git clone https://github.com/your-username/globin_project.git
cd globin_project

# 2. Set up virtual environment and install dependencies
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

# 3. Run the one-command master reproduction pipeline
python scripts/reproduce_all.py
```

### What `scripts/reproduce_all.py` Executes:
1. **Cluster Completeness Audit**: Validates that all 68 species and 121 loci meet strict syntenic and coding standards ([`verify_cluster_completeness.py`](file:///Users/ibrahim/Documents/globin_project/scripts/verify_cluster_completeness.py)).
2. **Overview & Methodology Figures**: Generates Figure 0 Dashboard Guide, 0A Human Paradigm, 0B Synteny Pipeline, 0C Codon Dynamics ([`generate_overview_methodology_figure.py`](file:///Users/ibrahim/Documents/globin_project/scripts/generate_overview_methodology_figure.py)).
3. **Macroevolutionary Figures**: Generates Figure 1A, 1B (Locus Divergence) and Figure 2 (Mechanistic Drivers) ([`generate_all_summary_figures.py`](file:///Users/ibrahim/Documents/globin_project/scripts/generate_all_summary_figures.py)).
4. **Phylogenetic Synteny & Codon Heatmaps**: Generates Figure 3 ([`generate_figure3_phylogenetic_synteny.py`](file:///Users/ibrahim/Documents/globin_project/scripts/generate_figure3_phylogenetic_synteny.py)).
5. **Bernardi's Compaction Law**: Executes global and clade regressions, fold-change ratios, and generates Figure 4 ([`analyze_bernardi_compaction.py`](file:///Users/ibrahim/Documents/globin_project/scripts/analyze_bernardi_compaction.py)).
6. **Standalone HTML Dashboard**: Rebuilds the self-contained static dashboard ([`generate_static_dashboard.py`](file:///Users/ibrahim/Documents/globin_project/scripts/generate_static_dashboard.py)).
7. **Test Suite**: Runs the complete test suite via `pytest`.

---

## Documentation Compendium

Comprehensive scientific background, methodology, and evolutionary syntheses are organized into 8 chapters in [`docs/`](file:///Users/ibrahim/Documents/globin_project/docs/):

| Chapter | Title | Key Scientific Findings |
| :---: | :--- | :--- |
| **[Chapter 01](file:///Users/ibrahim/Documents/globin_project/docs/01_scientific_aims.md)** | **Scientific Aims & Evolutionary Framework** | Macroevolutionary puzzle; selection vs. gBGC vs. isochore inertia; theoretical hypotheses. |
| **[Chapter 02](file:///Users/ibrahim/Documents/globin_project/docs/02_locus_identification_strategy.md)** | **Locus Identification & Curation Strategy** | Syntenic anchor framework (*RHBDF1*–*NPRL3*–*LUC7L*); teleost 3R duplication; Urodele intronic gigantism; CDS filtering. |
| **[Chapter 03](file:///Users/ibrahim/Documents/globin_project/docs/03_species_atlas.md)** | **Vertebrate Globin Species Atlas** | Complete catalog of 68 species and 121 intervals with coordinates, flankers, coding GC, and $Z$-scores. |
| **[Chapter 04](file:///Users/ibrahim/Documents/globin_project/docs/04_analytical_methodology.md)** | **Analytical Methodology & Statistics** | Mathematical formulation of $Z$-scores, intronic dilution, isochore coupling OLS regression, Mann-Whitney subtelomeric tests, and inertia relaxation kinetics ($t_{1/2} \approx 350\text{ My}$). |
| **[Chapter 05](file:///Users/ibrahim/Documents/globin_project/docs/05_user_guide_and_cli.md)** | **User Guide, CLI & Python API** | Installation instructions; `add-genome` and `export-tables` CLI commands; Python API recipes; error diagnostics. |
| **[Chapter 06](file:///Users/ibrahim/Documents/globin_project/docs/06_species_expansion_strategy.md)** | **Species Expansion Strategy** | Systematic scaling to 68 species across all vertebrate evolutionary extrema. |
| **[Chapter 07](file:///Users/ibrahim/Documents/globin_project/docs/07_chondrichthyan_telomeric_gradient_analysis.md)** | **Chondrichthyan Telomeric Gradient Analysis** | Evaluation of cartilaginous fish baselines and lack of isochoric compartmentalization. |
| **[Chapter 08](file:///Users/ibrahim/Documents/globin_project/docs/08_bernardi_compaction_law_and_developmental_dynamics.md)** | **Bernardi's Compaction Law & Developmental Dynamics** | Empirical proof across 72 baselines; paired embryonic vs. adult globin intron contrasts; two-way locus vs. stage decomposition. |

---

## Primary Results & Publication Figures

### Publication Figures
- **Figure 0 (Dashboard Guide & Methodology)**: [`results/summary_figure0_dashboard_guide.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure0_dashboard_guide.png)
  - Subpanel 0A (Human Paradigm): [`results/summary_figure0a_human_paradigm.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure0a_human_paradigm.png)
  - Subpanel 0B (Synteny Pipeline): [`results/summary_figure0b_synteny_pipeline.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure0b_synteny_pipeline.png)
  - Subpanel 0C (Metric Anatomy & GC3): [`results/summary_figure0c_metric_anatomy_and_gc3.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure0c_metric_anatomy_and_gc3.png)
- **Figure 1 (Macroevolutionary Locus Divergence)**: [`results/summary_figure1a_locus_divergence_summary.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure1a_locus_divergence_summary.png) (Curated Summary) & [`results/summary_figure1b_locus_divergence_full.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure1b_locus_divergence_full.png) (All 68 Species).
- **Figure 2 (Mechanistic Drivers of GC Polarity)**: [`results/summary_figure2_mechanistic_drivers.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure2_mechanistic_drivers.png).
- **Figure 3 (Phylogenetic Synteny & Codon Map)**: [`results/summary_figure3_phylogenetic_synteny.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure3_phylogenetic_synteny.png).
- **Figure 4 (Empirical Proof of Bernardi's Compaction Law)**: [`results/summary_figure_bernardi_compaction_law.png`](file:///Users/ibrahim/Documents/globin_project/results/summary_figure_bernardi_compaction_law.png).

### Primary Data Tables
- **Curated Genomic Intervals (121 Loci)**: [`results/curated_intervals.csv`](file:///Users/ibrahim/Documents/globin_project/results/curated_intervals.csv)
- **Verified Gene Models & Z-Scores (1,544 Genes)**: [`results/gene_zscores.csv`](file:///Users/ibrahim/Documents/globin_project/results/gene_zscores.csv)
- **Cross-Species Summary (68 Species)**: [`results/cross_species_summary_table.csv`](file:///Users/ibrahim/Documents/globin_project/results/cross_species_summary_table.csv)
- **Bernardi Compaction Global Stats**: [`results/bernardi_compaction_global_stats.csv`](file:///Users/ibrahim/Documents/globin_project/results/bernardi_compaction_global_stats.csv)
- **Embryonic vs. Adult Globin Compaction**: [`results/bernardi_compaction_embryonic_vs_adult.csv`](file:///Users/ibrahim/Documents/globin_project/results/bernardi_compaction_embryonic_vs_adult.csv)
- **Two-Way Locus vs. Stage Decomposition**: [`results/bernardi_compaction_two_way_locus_stage.csv`](file:///Users/ibrahim/Documents/globin_project/results/bernardi_compaction_two_way_locus_stage.csv)
- **Supplementary Tables S1–S4**: [`results/Supplementary_Table_S1_Species_Baselines.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S1_Species_Baselines.csv) &rarr; [`Supplementary_Table_S4_Statistical_Contrasts.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S4_Statistical_Contrasts.csv).

---

## Generating Single-Species Plot Galleries

The package allows generating 4 per-species publication plot galleries:
- `results/chromosome_plots/`: Linear chromosome tracks with centromeres and telomere distances.
- `results/bivariate_with_inset_plots/`: Bivariate $\log_{10}(\text{Length})$ vs. GC scatter with mini chromosome insets.
- `results/bivariate_length_gc_plots/`: Bivariate length vs. GC density distributions.
- `results/diagnostic_plots/`: Diagnostic chromosome distributions and globin percentiles.

```bash
# Generate all plots for a specific species:
globin-pipeline generate-plots --species Homo_sapiens --output-dir results/

# Generate across all 68 species:
globin-pipeline generate-plots --output-dir results/
```

---

## Automated Testing

```bash
# Run the test suite:
pytest tests/
```
All 16 test suites pass with 100% success.
