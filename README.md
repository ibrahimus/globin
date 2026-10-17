# globin-pipeline: Comparative Vertebrate Globin Genomics

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](file:///Users/ibrahim/Documents/globin_project/pyproject.toml)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-brightgreen.svg)](file:///Users/ibrahim/Documents/globin_project/pyproject.toml)
[![Tests](https://img.shields.io/badge/tests-14%20passed-success.svg)](file:///Users/ibrahim/Documents/globin_project/tests/)

A standardized, reproducible comparative genomics platform and Python library for auditing, discovering, curating, and analyzing $\alpha$- and $\beta$-globin loci across vertebrate genomes spanning ~500 million years of evolutionary history.

---

## Dataset Scale & Capabilities

- **Taxonomic Breadth**: **44 vertebrate genomes** across 9 major evolutionary clades (Cyclostomata, Holostei, Teleostei, Amphibia, Testudines, Lepidosauria, Archosauria, Monotremata, Eutheria).
- **Curated Syntenic Loci**: **80 verified genomic intervals** anchored by ultra-conserved orthologs (*RHBDF1*, *MPG*, *NPRL3*, *LUC7L*, *OR* clusters).
- **Annotated Repertoire**: **1,017 protein-coding genes** classified by structural role, GC compartment, and genomic position.
- **Empirical Baseline Normalization**: Zero-fallback whole-genome baseline distributions providing standardized $Z$-scores ($Z_{\text{gene}}$, $Z_{\text{cds}}$, $Z_{\text{gc3}}$).
- **Autonomous Syntenic Ingestion**: Extensible engine to stream, download, and curate any newly sequenced vertebrate genome from NCBI RefSeq/GenBank.

---

## Documentation Compendium

The project documentation is organized linearly into five comprehensive chapters in [`docs/`](file:///Users/ibrahim/Documents/globin_project/docs/):

| Chapter | Title | Key Contents |
| :---: | :--- | :--- |
| **[Chapter 1](file:///Users/ibrahim/Documents/globin_project/docs/01_scientific_aims.md)** | **Scientific Aims & Evolutionary Framework** | Macroevolutionary puzzle; selection vs. gBGC vs. isochore inertia; theoretical hypotheses. |
| **[Chapter 2](file:///Users/ibrahim/Documents/globin_project/docs/02_locus_identification_strategy.md)** | **Locus Identification & Curation Strategy** | Syntenic anchor framework (*RHBDF1*–*NPRL3*–*LUC7L*); teleost 3R duplication; Urodele intronic gigantism; CDS filtering. |
| **[Chapter 3](file:///Users/ibrahim/Documents/globin_project/docs/03_species_atlas.md)** | **Vertebrate Globin Species Atlas** | Complete catalog of 44 species and 80 intervals with coordinates, flankers, coding GC, and $Z$-scores. |
| **[Chapter 4](file:///Users/ibrahim/Documents/globin_project/docs/04_analytical_methodology.md)** | **Analytical Methodology & Statistics** | Mathematical formulation of $Z$-scores, intronic dilution, isochore coupling OLS regression, Mann-Whitney subtelomeric tests, and inertia relaxation kinetics ($t_{1/2} \approx 350\text{ My}$). |
| **[Chapter 5](file:///Users/ibrahim/Documents/globin_project/docs/05_user_guide_and_cli.md)** | **User Guide, CLI & Python API** | Installation instructions; `add-genome` and `export-tables` CLI commands; Python API recipes; error diagnostics. |

---

## Quickstart & Installation

Install the package in editable mode using `pip` or `uv`:

```bash
# Clone the repository
git clone https://github.com/your-username/globin_project.git
cd globin_project

# Install with dependencies
pip install -e .
# Or using uv:
uv pip install -e .
```

### Command-Line Interface (CLI)

```bash
# 1. Discover and ingest a new vertebrate genome from NCBI:
globin-pipeline add-genome \
  --accession GCF_036496625.1 \
  --species eptatretus_burgeri \
  --clade cyclostome

# 2. Add with manual coordinate override:
globin-pipeline add-genome \
  --accession GCF_036496625.1 \
  --species eptatretus_burgeri \
  --clade cyclostome \
  --chrom NC_140977.1 \
  --start 52050000 \
  --end 52280000

# 3. Export publication-ready supplementary tables:
globin-pipeline export-tables --output-dir results/supplementary_tables/

# 4. Generate all plot galleries (diagnostic, chromosome tracks, bivariate with insets):
globin-pipeline generate-plots --output-dir results/
```

### Python API Usage

```python
from globin_pipeline import GlobinPipeline

# Initialize pipeline with pre-packaged reference data
pipeline = GlobinPipeline()

# Query curated intervals
intervals = pipeline.get_intervals(species="Homo_sapiens")
for locus in intervals:
  print(
      f"{locus.locus_id:<15} | Chr {locus.chromosome}:{locus.start:,}-{locus.end:,} "
      f"({locus.span_bp / 1000:.1f} kb) | {locus.left_gene} ... {locus.right_gene}"
  )

# Access empirical genome baselines
baseline = pipeline.get_baseline(species="Homo_sapiens")
print(f"Human baseline CDS GC: {baseline.mean_cds_gc:.2f}% (std = {baseline.std_cds_gc:.2f}%)")
```

---

## Generating Plot Galleries & Visualizations

The package provides a built-in visualization engine to generate all 4 core publication plot galleries across all 44 species in static **PNG**, publication **PDF**, and interactive **Plotly HTML** formats:

### Command-Line Usage

```bash
# Generate all 4 galleries across all 44 species:
globin-pipeline generate-plots --output-dir results/

# Generate for a single species:
globin-pipeline generate-plots --species Homo_sapiens --output-dir results/

# Generate specific gallery types:
globin-pipeline generate-plots --plot-type chromosome --output-dir results/
globin-pipeline generate-plots --plot-type bivariate-inset --output-dir results/
```

### Gallery Structure & Formats

| Directory | Plot Type | Key Features | Formats |
| :--- | :--- | :--- | :--- |
| `results/diagnostic_plots/` | **Genome Distributions** | Whole-chromosome GC distributions, $\pm 1\sigma, \pm 2\sigma$ standard deviations, and individual globin percentiles. | PNG, HTML |
| `results/chromosome_plots/` | **Linear Chromosome Tracks** | Full chromosome arm length, centromere annotations, telomeric distances, and zoomed locus gene architecture. | PNG, PDF, HTML |
| `results/bivariate_length_gc_plots/` | **Bivariate Length vs. GC** | Bivariate $\log_{10}(\text{Coding-Span Length})$ vs. GC Content scatter with 2D density contours. | PNG, PDF, HTML |
| `results/bivariate_with_inset_plots/` | **Bivariate with Mini Inset** | Exaggerated mini chromosome linear inset in the top-right corner showing chromosome position alongside bivariate distribution. | PNG, PDF, HTML |

### Python API Usage

```python
import pandas as pd
from globin_pipeline.visualizer import Visualizer

viz = Visualizer()
df_genes = pd.read_csv("results/gene_zscores.csv")
df_intervals = pd.read_csv("results/curated_intervals.csv")

# Diagnostic distribution plot
viz.plot_diagnostic("Homo_sapiens", df_genes, output_dir="results/diagnostic_plots")

# Linear chromosome track with centromere & telomere margins
viz.plot_chromosome_linear("Homo_sapiens", df_genes, df_intervals=df_intervals, output_dir="results/chromosome_plots")

# Bivariate log10(Length) vs GC
viz.plot_length_vs_gc("Homo_sapiens", df_genes, output_dir="results/bivariate_length_gc_plots")

# Bivariate with mini chromosome linear inset
viz.plot_bivariate_with_chromosome_inset("Homo_sapiens", df_genes, df_intervals=df_intervals, output_dir="results/bivariate_with_inset_plots")
```

---

## Primary Results & Visualizations

- **Supplementary Tables**:
  - [`Supplementary_Table_S1_Species_Baselines.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S1_Species_Baselines.csv): Whole-genome coding and genomic baselines for all 44 species.
  - [`Supplementary_Table_S2_Globin_Loci_Summary.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S2_Globin_Loci_Summary.csv): Coordinates, gene counts, and flankers for all 80 intervals.
  - [`Supplementary_Table_S3_All_1017_Genes_Annotations_ZScores.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S3_All_1017_Genes_Annotations_ZScores.csv): Annotation database and $Z$-scores for all 1,017 genes.
  - [`Supplementary_Table_S4_Statistical_Contrasts.csv`](file:///Users/ibrahim/Documents/globin_project/results/Supplementary_Table_S4_Statistical_Contrasts.csv): Hypothesis tests, regressions, effect sizes, and biological conclusions.
- **Figures & Visualizations**:
  - Standardized Meta-Comparative Graphs: [`results/meta_comparative_graph.png`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph.png) / [`results/meta_comparative_graph_clean.png`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph_clean.png)
  - Diagnostic Bivariate Distribution Plots: [`results/diagnostic_plots/`](file:///Users/ibrahim/Documents/globin_project/results/diagnostic_plots/) (44 species with centromere annotations and telomeric margins).
