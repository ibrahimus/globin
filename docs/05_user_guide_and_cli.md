# Chapter 5: User Guide, Command-Line Interface & Python API

## 1. Installation & Environment Setup

`globin-pipeline` is a fully standalone, installable Python package designed for comparative genomics of vertebrate globin repertoires. It bundles curated interval definitions, baseline statistics, and chromosome models for 44 vertebrate species, while providing modular tools to ingest and analyze new genomes directly from NCBI.

### Requirements
- **Python**: Version $\ge 3.10$
- **Supported Platforms**: macOS (x86_64, arm64 / Apple Silicon), Linux (x86_64), Windows (WSL2 recommended)

### Installation
You can install the package in editable mode using `pip` or the high-speed `uv` package manager:

```bash
# Clone the repository
git clone https://github.com/your-username/globin_project.git
cd globin_project

# Option A: Standard pip installation
pip install -e .

# Option B: Fast installation using uv
uv pip install -e .
```

To verify the installation:
```bash
globin-pipeline --version
# globin-pipeline 1.0.0
```

---

## 2. Command-Line Interface (CLI) Reference

The package provides the unified `globin-pipeline` command with multiple subcommands:

```
globin-pipeline [OPTIONS] COMMAND [ARGS]...
```

### 2.1 Discovering and Adding New Genomes (`add-genome`)
The `add-genome` subcommand downloads a vertebrate assembly from NCBI RefSeq/GenBank, parses GFF annotations and CDS sequences, automatically detects candidate globin loci via synteny anchors, and registers the genome into your local manifest.

```bash
# Automated discovery via syntenic anchors:
globin-pipeline add-genome \
  --accession GCF_036496625.1 \
  --species eptatretus_burgeri \
  --clade cyclostome

# Manual coordinate override when locus boundaries are pre-determined:
globin-pipeline add-genome \
  --accession GCF_036496625.1 \
  --species eptatretus_burgeri \
  --clade cyclostome \
  --chrom NC_140977.1 \
  --start 52050000 \
  --end 52280000
```

**Options:**
- `--accession <TEXT>`: NCBI RefSeq or GenBank assembly accession (e.g., `GCF_000001405.40`). *(Required)*
- `--species <TEXT>`: Normalized species name (e.g., `homo_sapiens`). Defaults to accession if omitted.
- `--clade <TEXT>`: Phylogenetic clade (e.g., `amniote_mammal`, `amphibian`, `teleost_fish`).
- `--chrom <TEXT>`: Target chromosome identifier (for manual coordinate override).
- `--start <INT>`: 1-based start coordinate.
- `--end <INT>`: 1-based end coordinate.
- `--output-dir <PATH>`: Working directory for downloaded files (defaults to `data/<species>`).
- `--skip-download`: Skip FTP download if GFF and FASTA files are already present locally.

---

### 2.2 Exporting Standard Supplementary Tables (`export-tables`)
Exports all publication-ready supplementary tables (S1 through S4) as standardized CSV files:

```bash
globin-pipeline export-tables --output-dir results/supplementary_tables/
```

**Exported Artifacts:**
1. `Supplementary_Table_S1_Species_Baselines.csv`: Empirical whole-genome coding and genomic GC baselines for all 44 species.
2. `Supplementary_Table_S2_Globin_Loci_Summary.csv`: Curated interval coordinates, syntenic flanking genes, spans, and gene counts for all 80 loci.
3. `Supplementary_Table_S3_All_1017_Genes_Annotations_ZScores.csv`: Comprehensive gene-by-gene database of all 1,017 genes, including exon count, CDS GC, GC3, whole-gene GC, and standardized $Z$-scores.
4. `Supplementary_Table_S4_Statistical_Contrasts.csv`: Summary of hypothesis tests (paired mammalian $t$-test, avian polarity, squamate hyper-gBGC, isochore coupling OLS regression, Mann-Whitney subtelomeric test).

---

### 2.3 Running Batch Manifest Pipelines (`run-manifest`)
Executes end-to-end data processing across all species defined in a YAML manifest:

```bash
globin-pipeline run-manifest --manifest globin_pipeline/data/species_manifest.yaml
```

---

## 3. Python API Recipes

`globin-pipeline` exposes a rich, typed Python API for programmatic workflows, exploratory data analysis, and Jupyter notebook integration.

### Recipe 1: Initializing the Pipeline and Accessing Bundled Data
```python
from globin_pipeline import GlobinPipeline

# Initialize pipeline (loads pre-packaged 44-species reference data)
pipeline = GlobinPipeline()

# List all available species
species_list = pipeline.list_curated_species()
print(f"Total curated species: {len(species_list)}")

# Inspect curated loci for a specific species
human_intervals = pipeline.get_intervals(species="Homo_sapiens")
for locus in human_intervals:
  print(
      f"Locus: {locus.locus_id:<15} | Chr {locus.chromosome}:{locus.start:,}-{locus.end:,} "
      f"({locus.span_bp / 1000:.1f} kb) | Anchors: {locus.left_gene} ... {locus.right_gene}"
  )
```

---

### Recipe 2: Automated Syntenic Locus Detection in a Custom GFF
```python
from pathlib import Path
from globin_pipeline.synteny import GlobinSyntenyDetector

# Initialize the detector with canonical Gnathostome and Cyclostome anchor dictionaries
detector = GlobinSyntenyDetector()

# Scan an annotation file for candidate globin loci
candidates = detector.detect_candidates(
    gff_path=Path("path/to/genomic.gff"), min_flank_distance_bp=50_000
)

for cand in candidates:
  print(f"Candidate Locus on {cand.chromosome}:")
  print(f"  Coordinates: {cand.start:,} - {cand.end:,}")
  print(f"  Globin genes found: {cand.globin_count}")
  print(f"  Flanking anchors: {cand.flanking_anchors}")
```

---

### Recipe 3: Calculating Standardized $Z$-Scores Against Empirical Baselines
```python
from globin_pipeline import GlobinPipeline

pipeline = GlobinPipeline()

# Retrieve baseline metrics for House Mouse
baseline = pipeline.get_baseline(species="Mus_musculus")
print(
    f"Mouse Baseline: Mean CDS GC = {baseline.mean_cds_gc:.2f}%, Std = {baseline.std_cds_gc:.2f}%"
)

# Compute standardized Z-score for an observed coding GC of 60.5%
observed_gc = 60.5
z_score = (observed_gc - baseline.mean_cds_gc) / baseline.std_cds_gc
print(f"Calculated Z-Score: {z_score:+.2f} sigma")
```

---

## 4. Error Handling & Diagnostics

The package implements a comprehensive hierarchy of domain-specific exceptions in `globin_pipeline.exceptions` to provide clear, actionable diagnostics:

| Exception Class | Trigger Condition | Recommended User Action |
| :--- | :--- | :--- |
| `GenomeNotFoundError` | The requested accession or URL was not found on NCBI FTP servers. | Verify assembly accession format (`GCF_...` or `GCA_...`). |
| `AnnotationMissingError` | The downloaded assembly lacks a valid GFF3 file or contains no protein-coding CDS features. | Ensure the assembly has an NCBI RefSeq annotation release. |
| `ScaffoldOnlyError` | The assembly consists solely of fragmented unplaced scaffolds without chromosome mappings. | Locus detection can still identify gene clusters, but telomeric distance metrics will be skipped. |
| `LocusNotFoundError` | Automated synteny detection did not detect any globin genes or syntenic anchor clusters. | Provide manual locus coordinates via `--chrom`, `--start`, `--end`. |
| `InvalidCoordinateError` | Specified interval coordinates are negative, inverted (`start > end`), or exceed chromosome length. | Check chromosome length using `pipeline.get_chromosome_lengths()`. |
| `BaselineComputationError` | Genome contains fewer than 1,000 verified protein-coding genes. | Verify that GFF parsing included full coding annotations. |

### Example Diagnostic Handling
```python
from globin_pipeline import GlobinPipeline
from globin_pipeline.exceptions import InvalidCoordinateError, LocusNotFoundError

pipeline = GlobinPipeline()

try:
  pipeline.add_genome(
      accession="GCF_999999999.1",
      chrom="chr1",
      start=5000000,
      end=1000000,  # Inverted coordinates
  )
except InvalidCoordinateError as e:
  print(f"Configuration error: {e}")
except LocusNotFoundError as e:
  print(f"Synteny detection warning: {e}")
```
