# Implementation Plan: Advanced Genome-Wide Globin GC/GC3 Baseline & Multi-Locus Comparative Analysis

## Goal Description

We want to expand the comparative globin analysis pipeline to investigate fundamental evolutionary and compositional questions:
1. **Genomic Elevation**: Do globin genes possess significantly higher, lower, or normal GC/GC3 levels relative to the overall protein-coding gene repertoire of their host genome?
2. **Multiple Loci Comparison**: In genomes with multiple globin loci (notably teleost fish 3R duplicates: Locus A vs. Locus B), do globin genes in different genomic locations exhibit distinct GC profiles?
3. **Cluster Bystander Effect**: Are non-globin genes located inside or flanking the globin cluster (*RHBDF1, MPG, NPRL3, LUC7L, KANK2*) similarly GC-elevated, indicating a regional isochore / domain effect?
4. **Pseudogene Evolution**: Do globin pseudogenes evolve differently from functional protein-coding globins (e.g., drifting toward lower GC / genomic background due to loss of GC-biased gene conversion and purifying selection)?

To achieve this, we will:
- Download and stream species-specific whole-genome CDS FASTA files from Ensembl to compute whole-genome GC and GC3 distributions across all protein-coding genes (~20,000–30,000 genes per species).
- Persist and cache these genome baselines in `data/genome_baselines/<species>_genes_gc.parquet` (or `.csv.gz`) so they are computed only once per genome.
- Compute standardized **Z-scores** ($Z = \frac{\text{GC} - \mu}{\sigma}$) and percentile ranks for globin coding genes, globin pseudogenes, and intra-cluster non-globin genes against each species' own baseline.
- Generate **high-resolution diagnostic graphs** per species (density distribution, box/strip comparisons, GC vs. GC3 plots).
- Generate a unified **cross-species meta-graph** comparing all 13 vertebrate species along their phylogenetic tree.

---

## User Review Required

> [!IMPORTANT]
> **Whole-Genome Baseline Streaming vs Local Cache**:
> We tested streaming Ensembl CDS FASTA files directly via Python HTTP/gzip streaming. Across our 13 vertebrate species, downloading and parsing all 20,000–30,000 genes takes only **6–20 seconds per genome**!
> Once parsed, we cache each genome's table in compressed parquet/CSV format in `data/genome_baselines/`. Subsequent queries are instantaneous (<0.05 seconds).

> [!NOTE]
> **Pseudogene GC Metric**:
> Functional globin genes have intact coding sequences (CDS), allowing calculation of both CDS GC% and synonymous third-codon position GC% (GC3).
> Globin pseudogenes often have stop codons, frameshifts, or deleted exons; they do not have true functional codons.
> We evaluate globin pseudogenes using **exonic/transcript GC%** and compare them directly against active globin exonic/genomic GC% and the genome-wide baseline.

---

## Open Questions

> [!TIP]
> 1. **Visualization Format**: We propose outputting high-resolution publication-quality PNG and PDF formats (e.g. 300 DPI) using `seaborn` and `matplotlib`. Would you also like interactive HTML plots (via `plotly`)?
> 2. **GC Metric for Primary Z-score**: In isochore genomics, **GC3 (third codon position)** is the most sensitive metric for GC-biased gene conversion (gBGC), while **CDS GC** is the standard whole-coding metric. We propose computing and plotting Z-scores for **both CDS GC and GC3**, featuring both in the diagnostic panels. Does this match your intended focus?

---

## Architecture & Proposed Changes

```
globin_project/
├── data/
│   ├── cache/                           # REST and locus sequence cache
│   └── genome_baselines/                # [NEW] Cached whole-genome gene GC distributions
│       ├── homo_sapiens_genes_gc.parquet
│       ├── danio_rerio_genes_gc.parquet
│       └── ...
├── results/
│   ├── diagnostic_plots/                # [NEW] Diagnostic figures per species
│   │   ├── diagnostic_Homo_sapiens.png
│   │   ├── diagnostic_danio_rerio.png
│   │   └── ...
│   ├── meta_comparative_graph.png       # [NEW] Cross-species meta-graph
│   ├── gene_zscores.csv                 # [NEW] Full table with Z-scores & percentiles
│   └── genome_baseline_stats.csv        # [NEW] Genome means, stds, medians
└── globin_pipeline/
    ├── baseline_manager.py              # [NEW] Computes & caches genome-wide gene GC
    ├── zscore_analyzer.py               # [NEW] Computes Z-scores, ranks, statistical tests
    ├── visualizer.py                    # [NEW] Matplotlib/Seaborn diagnostic & meta plots
    ├── curator.py                       # [UPDATE] Ensure pseudogenes are properly flagged
    ├── annotator.py                     # [UPDATE] Extract globin pseudogenes with exonic coordinates
    └── cli.py                           # [UPDATE] Add baseline, zscore, and plot commands
```

---

### Component Specifications

#### 1. Whole-Genome Baseline Manager ([NEW] `globin_pipeline/baseline_manager.py`)
- Streams `https://ftp.ensembl.org/pub/current_fasta/<species>/cds/<Species>.<Assembly>.cds.all.fa.gz`.
- For each transcript:
  - Parses gene ID, transcript ID, gene symbol, biotype, description from the FASTA header.
  - Keeps the canonical/longest CDS per protein-coding gene.
  - Computes CDS length, CDS GC%, and GC3% (using codon positions 2, 5, 8...).
- Aggregates distribution statistics:
  - Mean CDS GC ($\mu_{\text{GC}}$), standard deviation ($\sigma_{\text{GC}}$), median, IQR.
  - Mean GC3 ($\mu_{\text{GC3}}$), standard deviation ($\sigma_{\text{GC3}}$), median, IQR.
- Caches to `data/genome_baselines/<species>_genes_gc.parquet` and `data/genome_baselines/<species>_stats.json`.
- Implements `get_genome_baseline(species)`: loads cache if present, streams and calculates if not.

#### 2. Z-Score & Statistical Analyzer ([NEW] `globin_pipeline/zscore_analyzer.py`)
- For every gene in the curated loci:
  - Matches gene ID against the species' genome-wide baseline.
  - Computes:
    - $Z_{\text{CDS\_GC}} = \frac{\text{CDS\_GC} - \mu_{\text{GC}}}{\sigma_{\text{GC}}}$
    - $Z_{\text{GC3}} = \frac{\text{GC3} - \mu_{\text{GC3}}}{\sigma_{\text{GC3}}}$
    - Percentile rank within the host genome ($0.0 - 100.0\%$).
- Categorizes each locus gene into 4 comparative classes:
  1. `globin_coding_locus_A`: Protein-coding globins in primary cluster (or single locus).
  2. `globin_coding_locus_B`: Protein-coding globins in secondary cluster (teleost 3R duplicate).
  3. `globin_pseudogene`: Globin pseudogenes (*HBAP1, HBZP1*, etc.).
  4. `cluster_non_globin`: Syntenic non-globin bystander genes (*NPRL3, MPG, RHBDF1, LUC7L, KANK2*).
- Performs hypothesis testing:
  - Mann-Whitney U test & Welch's t-test comparing Globin coding genes vs. Genome background.
  - Comparison between Locus A vs. Locus B in teleosts.
  - Comparison between Globin coding genes vs. Globin pseudogenes.
  - Comparison between Globin coding genes vs. Cluster non-globin genes.

#### 3. Diagnostic & Meta-Plot Visualizer ([NEW] `globin_pipeline/visualizer.py`)
- **Per-Species Diagnostic Graph** (`results/diagnostic_plots/diagnostic_<species>.png`):
  - **Panel A (Whole-Genome Density & Placement)**:
    - Kernel density estimate (KDE) and histogram of genome-wide gene GC% and GC3%.
    - Vertical lines at $\mu$, $\pm 1\sigma$, $\pm 2\sigma$.
    - Labeled markers / rug ticks indicating each globin gene (colored by locus), pseudogenes (diamond), and intra-cluster non-globins (triangle).
  - **Panel B (Group Comparison Boxplot & Jitter)**:
    - Categorical comparison of: Genome Baseline vs. Cluster Non-globins vs. Globin Pseudogenes vs. Globin Coding (Locus A / Locus B).
    - Statistical significance asterisks ($*** p < 0.001$, etc.).
  - **Panel C (GC vs. GC3 Bivariate Space)**:
    - Scatter plot of locus genes overlaid on whole-genome density contours.
- **Cross-Species Meta-Graph** (`results/meta_comparative_graph.png`):
  - Forest / strip plot spanning all 13 vertebrate species arranged phylogenetically:
    - Lamprey $\rightarrow$ Zebrafish $\rightarrow$ Medaka $\rightarrow$ Fugu $\rightarrow$ Xenopus $\rightarrow$ Turtle $\rightarrow$ Chicken $\rightarrow$ Opossum $\rightarrow$ Hedgehog $\rightarrow$ Mouse $\rightarrow$ Pig $\rightarrow$ Cow $\rightarrow$ Human.
  - X-axis: Standardized Z-score ($Z = 0$ is genome baseline, shaded $[-1, +1]$ and $[-2, +2]$ standard deviations).
  - Y-axis: Species.
  - Plotted points:
    - Globin Coding Locus A (Blue circles)
    - Globin Coding Locus B (Teal squares)
    - Globin Pseudogenes (Purple diamonds)
    - Intra-cluster Non-globins (Green triangles)
  - This graph immediately reveals the macro-evolutionary trends across all clades!

#### 4. CLI Extensions ([MODIFY] `globin_pipeline/cli.py`)
- Add CLI commands:
  - `globin-pipeline baseline [--species <name|all>]`
  - `globin-pipeline zscore [--output <path>]`
  - `globin-pipeline plot-diagnostic [--species <name|all>]`
  - `globin-pipeline plot-meta [--output <path>]`
  - `globin-pipeline advanced-analysis [--output-dir <path>]` (runs all steps automatically)

---

## Verification Plan

### Automated Tests
- `tests/test_baseline_manager.py`:
  - Test streaming and parsing mock CDS FASTA.
  - Verify calculation of $\mu$, $\sigma$, median, GC, GC3.
  - Verify parquet/CSV caching and reload speed.
- `tests/test_zscore_analyzer.py`:
  - Test Z-score arithmetic and percentile calculation on mock gene sets.
  - Verify statistical hypothesis tests (Mann-Whitney U).
- `tests/test_visualizer.py`:
  - Test figure generation and save functions without display error (`matplotlib.use('Agg')`).

```bash
uv run --with pytest,matplotlib,seaborn,scipy pytest tests/
```

### Manual Verification
- Review generated per-species diagnostic graphs for Human, Zebrafish, Xenopus, Chicken, and Turtle.
- Inspect the cross-species meta-graph to verify phylogenetic ordering, cluster differentiation (Locus A vs B), pseudogene markers, and bystander gene locations.
- Verify that calculated Z-scores align with established literature (e.g. human globin GC3 $> +1.8\sigma$ above genome mean).
