# Implementation Plan: Vertebrate Globin Loci Comparative Analysis Pipeline

## Goal Description
The objective of this project is to characterize vertebrate alpha (and later beta) globin loci across major evolutionary clades—specifically comparing **fish**, **amphibians**, and **amniotes** (reptiles, birds, mammals). 

The input is provided in periodic collaborative batches from an external bioinformatician ("Mike") in [`1_data_from_mike/<date>/`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930). 
Our goals are to:
1. **Audit and summarize Mike's initial data batch** (`20260930`), identifying syntenic, coordinate, and assembly issues.
2. **Standardize and curate genomic intervals** for the alpha globin cluster across target species and human reference.
3. **Extract all protein-coding genes** within each interval from species-specific GTF/GFF annotations.
4. **Classify genes** into globin family members (alpha-like, beta-like, ancestral/other globins) vs. flanking/intervening non-globin genes.
5. **Extract genomic sequences (FASTA) and calculate GC content** metrics across compartments (whole locus, CDS, GC3 third-codon positions, introns, intergenic) and sliding windows (10 kb, 100 kb isochores), benchmarking against genomic backgrounds.
6. **Build an automated, reproducible Python package/CLI** (`globin-analyzer`) that processes new batches from Mike automatically.

---

## Detailed Audit & Summary of Mike's Files (`20260930`)

The folder [`1_data_from_mike/20260930`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930) contains four primary files:
- [`file_descriptions.txt`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930/file_descriptions.txt): Brief description of attachments.
- [`homologues_config_multi_gene.yaml`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930/homologues_config_multi_gene.yaml): Config defining `query_species: Homo_sapiens`, 9 query genes (*RHBDF1, MPG, NPRL3, HBA2, LUC7L, KANK2, LCMT1, DCHS1, STIM1*), and 12 target species.
- [`Homo_sapiens_selected_flanking_homologues.csv`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930/Homo_sapiens_selected_flanking_homologues.csv): Raw Ensembl Compara orthologue query results.
- [`alpha_globin_domain_intervals.csv`](file:///Users/ibrahim/Documents/globin_project/1_data_from_mike/20260930/alpha_globin_domain_intervals.csv): The derived genomic intervals spanning the alpha globin domain for 11 target species (turtle removed).

### Critical Findings & Problems Identified in Mike's Data

```
Human:      [RHBDF1] --- [MPG] --- [NPRL3] -------- [HBA2 / Globins] -------- [LUC7L]
                                                           |
Xenopus:    [Globins: hba1,3,4, hbz] --- [NPRL3] --- [MPG] --- [RHBDF1] ------------- [LUC7L]
            |----------------- EXCLUDED -----------------| |==== MIKE'S INTERVAL ====|
                                                           (ZERO alpha globins inside!)
```

### 1. Xenopus tropicalis interval is broken (contains 0 alpha globins)
- **What Mike did**: Set left boundary to `rhbdf1` (chr9: 40,712,133) and right boundary to `luc7l` (chr9: 41,211,821), padded with 5 kb -> `chr9: 40,707,133 - 41,216,821`.
- **The Problem**: In *Xenopus tropicalis*, the locus order on chromosome 9 is inverted/rearranged. The alpha globin cluster (*hbz, hba4, hba3, hba1*) is located at **chr9: 40,366,583 – 40,458,276**, followed by *nprl3* (40.54–40.60 Mb) and *mpg* (40.61–40.71 Mb), then *rhbdf1* (40.71–40.82 Mb).
- **Consequence**: Mike's interval **completely excludes every single alpha globin gene** in Xenopus! It only contains *rhbdf1*, an embryonic beta-like globin (*hbe1.2* at 40.87 Mb), and *luc7l*.
- **Fix**: The Xenopus interval must span from *hba1/snrnp25* (~40.31 Mb) through *luc7l* (~41.22 Mb), i.e. `chr9: 40,315,000 – 41,220,000`.

### 2. Teleost Fish Dual Loci (Teleost Genome Duplication / 3R)
- **What Mike did**: In *Danio rerio* (zebrafish), *Oryzias latipes* (medaka), and *Takifugu rubripes* (fugu), Mike used *rhbdf1a* and *kank2* to define intervals on Chr 3, Chr 8, and Chr 5 respectively.
- **The Problem**: Teleosts experienced an ancestral whole-genome duplication (TGD / 3R), resulting in **two distinct paralogous globin loci**:
  - **Locus A (Major Cluster)**: Synteny with *rhbdf1a - mpg - nprl3 - [alpha+beta globins] - kank2* (Zebrafish Chr 3, Medaka Chr 8, Fugu Chr 5).
  - **Locus B (Minor Cluster)**: Synteny with *rhbdf1b - [embryonic alpha+beta globins: hbae4, hbae5, hbbe2, hbbe3] - lcmt1* (Zebrafish Chr 12, Medaka Chr 19, Fugu Chr 1).
- **Consequence**: In `Homo_sapiens_selected_flanking_homologues.csv`, Mike's Ensembl Compara query for HBA2 matched `hbae4` on **Locus B** (Chr 12/19/1), but his intervals table defined **Locus A** (Chr 3/8/5). Furthermore, Locus B is omitted from the intervals table entirely.
- **Fix**: Our pipeline must track and analyze **both** Locus A and Locus B for teleost fish, as both are critical to understanding globin evolution in fish.

### 3. Western Painted Turtle (*Chrysemys picta bellii*) Removed
- **What Mike did**: Removed turtle with note: `(Western painted turtle removed)`.
- **Why it failed**: On turtle scaffold `KK088206.1`, the genes run: Globins (1.42–1.45 Mb) -> *NPRL3* (1.46–1.59 Mb) -> *RHBDF1* (1.70–1.83 Mb). Mike attempted to make RHBDF1 the "left" boundary (1.70 Mb) and tried to set an arbitrary "right" boundary at `NPRL3 + 100kb` (1.57 Mb). Since 1.70 Mb > 1.57 Mb, the coordinates were reversed and the globins (1.42 Mb) were excluded.
- **Fix**: Re-instate turtle with valid boundaries: `KK088206.1: 1,415,000 – 1,865,000` (capturing *ENSCPBG00000026394, ENSCPBG00000026395 [HBA], HBAD, HBZ, NPRL3, MPG, RHBDF1*).

### 4. Hedgehog (*Erinaceus europaeus*) Fragmented Assembly
- **The Problem**: In assembly `HEDGEHOG`, `GeneScaffold_103` only contains *RHBDF1* and *NPRL3*. The globin genes reside on separate unplaced scaffolds (`scaffold_235869` with *GLNC1/HBA* and `scaffold_219931` with *HBM*).
- **Consequence**: Mike's artificial fixed endpoint (+100 kb from NPRL3 on `GeneScaffold_103`) results in an interval containing **zero globin genes**.
- **Fix**: Flag fragmented assemblies where loci span multiple contigs/scaffolds, or update to modern chromosome-level hedgehog assemblies where available.

### 5. Sea Lamprey (*Petromyzon marinus*) Agnathous Globins
- **The Problem**: Lamprey diverged prior to the alpha-beta globin gene duplication. Ensembl Compara returned "Not found" for HBA2.
- **Finding**: On scaffold `GL478636`, directly adjacent to *nprl3* (45–57 kb), reside two unannotated genes: `ENSPMAG00000001587` (70–80 kb) and `ENSPMAG00000001592` (89–93 kb). Protein sequence verification confirms both are true cyclostome hemoglobins! Mike's interval (`20,362 – 150,362`) captures them, but standard annotation tools won't label them as globins without sequence/domain profiling.

### 6. Missing *Homo sapiens* Reference Interval
- Mike's `alpha_globin_domain_intervals.csv` lists 11 target species but omits *Homo sapiens* itself. Human (`GRCh38`, `chr16: 60,000 – 275,000`) must be included as the amniote reference.

### 7. Coordinate Convention Discrepancy
- The columns `Start` and `End` in Mike's CSV are **0-based half-open** (`[Start, End)`), while the `location` text column and Ensembl web displays use **1-based closed** coordinates (`[Start + 1, End]`).
- Our package will enforce explicit coordinate typing to prevent 1-bp off-by-one errors during FASTA slicing or GTF feature overlapping.

---

## User Review Required

> [!IMPORTANT]
> **Action Needed on Xenopus Interval**: Mike's current interval for *Xenopus tropicalis* completely excludes all alpha globins. We propose replacing it with `chr9: 40,315,000 – 41,220,000` to capture all 4 alpha globins (*hbz, hba4, hba3, hba1*), *nprl3*, *mpg*, *rhbdf1*, and *luc7l*.

> [!IMPORTANT]
> **Action Needed on Teleost Fish**: Teleosts have two paralogous globin loci (Locus A on Chr 3/8/5 and Locus B on Chr 12/19/1). We propose analyzing **both** loci for each teleost species rather than choosing only one.

> [!NOTE]
> **Action Needed on Western Painted Turtle**: We can re-include *Chrysemys picta bellii* using coordinates `KK088206.1: 1,415,000 – 1,865,000`.

---

## Open Questions

> [!WARNING]
> 1. **Assembly Versions**: Mike used Ensembl release assemblies (e.g. `ASM229v1` for opossum, `HEDGEHOG` for hedgehog). Some of these older assemblies are fragmented. Should we stick strictly to the exact assembly versions Mike used in `homologues_config_multi_gene.yaml`, or would you like an option in the package to upgrade fragmented assemblies (like hedgehog) to newer chromosome-level assemblies when available?
> 2. **Globin Classification Methodology**: In non-mammalian and marsupial species, many globin genes lack symbols like "HBA" or descriptions mentioning "hemoglobin" (e.g. lamprey `ENSPMAG00000001587`, opossum `ENSMODG00000042487`). To ensure 100% sensitivity, should we classify globins using a hybrid strategy: (a) Ensembl metadata/GO terms + (b) Pfam/InterPro Globin domain (PF00042) match + (c) profile HMM / BLAST search?
> 3. **Isochore & Background GC Window Sizes**: For genomic background comparisons, what window sizes and parameters do you prefer? The standard is 10 kb and 100 kb non-overlapping windows across the whole chromosome and whole genome, plus gene-level GC (overall gene, CDS, GC3, introns, intergenic). Are there other specific windows or metrics you would like included?

---

## Proposed Changes & Package Architecture

We will organize the code into a modular Python package located in the repository:
`/Users/ibrahim/Documents/globin_project/globin_pipeline/`

```
globin_project/
├── 1_data_from_mike/                     # Raw batches from Mike
│   └── 20260930/
├── data/                                 # Downloaded/cached FASTA and GTF/GFF
├── results/                              # Output tables, reports, figures
└── globin_pipeline/                      # Core analysis package
    ├── __init__.py
    ├── parser.py                         # Parses & audits Mike's batches
    ├── curator.py                        # Reconciles intervals & overrides
    ├── data_manager.py                   # Ensembl REST / FTP fetcher for GTF & FASTA
    ├── annotator.py                      # Extracts genes from GTF/GFF & classifies globins
    ├── gc_analyzer.py                    # GC content, GC3, and isochore statistics
    ├── reporter.py                       # Summary tables, comparative metrics & plots
    └── cli.py                            # Unified CLI interface
```

### Component Breakdown

#### 1. Batch Parser & Auditor ([NEW] `globin_pipeline/parser.py`)
- Reads directory batches (`<batch_dir>/` containing YAML, homologues CSV, intervals CSV).
- Validates data integrity: checks for missing species, missing genes, and coordinate math.
- Detects whether flanking genes actually enclose the globin genes by inspecting homology records and synteny order.
- Generates an automated **Audit Report** flagging discrepancies (like the Xenopus inversion or teleost dual loci).

#### 2. Locus Curator ([NEW] `globin_pipeline/curator.py`)
- Maintains curated rules and overrides for edge cases where synteny breaks.
- Re-introduces *Homo sapiens* reference interval (`chr16: 60,000 – 275,000`).
- Re-introduces *Chrysemys picta bellii* (`KK088206.1: 1,415,000 – 1,865,000`).
- Corrects *Xenopus tropicalis* (`chr9: 40,315,000 – 41,220,000`).
- Splits teleost species into `Locus_A` (major cluster) and `Locus_B` (minor cluster).
- Outputs a normalized, unified interval specification table in standard formats (CSV, BED).

#### 3. Data Manager ([NEW] `globin_pipeline/data_manager.py`)
- Connects to Ensembl REST API and Ensembl FTP.
- Automates retrieval of:
  - Genomic sequence of the locus (FASTA)
  - Full chromosome / genome FASTA (for background genomic GC / isochore calculation)
  - Species-specific GTF / GFF3 annotation files
- Caches downloaded files in `data/<species>/<assembly>/` to prevent redundant network traffic.

#### 4. Gene Extraction & Globin Classifier ([NEW] `globin_pipeline/annotator.py`)
- Parses GTF/GFF3 for the target interval using `gffutils` or Ensembl REST.
- Filters for `protein_coding` genes (while recording pseudogenes and ncRNAs).
- Employs a robust 3-tier classification system:
  1. **Annotation Match**: Checks gene symbols and descriptions for globin keywords (`HBA`, `HBB`, `HBZ`, `HBM`, `HBQ`, `globin`, etc.).
  2. **Domain Match**: Checks InterPro / Pfam domain `PF00042` (Globin family).
  3. **HMM Profile / Alignment Match**: Scans translations against vertebrate globin profile HMMs to capture unannotated/novel globins.
- Classifies each gene into:
  - `alpha_globin` (*HBA, HBZ, HBM, HBQ, etc.*)
  - `beta_globin` (*HBB, HBE, etc.* in fish/amphibian shared loci)
  - `other_globin` (*cytoglobin, myoglobin, neuroglobin, unclassified*)
  - `flanking_non_globin` (*RHBDF1, MPG, NPRL3, LUC7L, KANK2, etc.*)
  - `intervening_non_globin` (any non-globin gene located inside the cluster)

#### 5. Sequence & GC Content Analyzer ([NEW] `globin_pipeline/gc_analyzer.py`)
- Computes comprehensive GC metrics:
  - **Locus Overall GC**: %GC across the entire genomic interval.
  - **Genic vs Intergenic GC**: %GC in combined exons vs introns vs intergenic spaces.
  - **Coding Sequence GC (CDS GC) & GC3**: GC content across CDS, and specifically at the third codon position (GC3 / GC at synonymous sites), which is the most sensitive indicator of isochores and biased gene conversion.
  - **Individual Gene GC Profiles**: GC metrics for each individual globin vs non-globin gene.
  - **Isochore Sliding Windows**: 10 kb and 100 kb sliding windows across the locus.
  - **Genomic Baseline Benchmarking**:
    - Whole-chromosome %GC
    - Whole-genome %GC
    - Genome-wide CDS %GC and GC3 distribution
    - Genome-wide 10kb and 100kb isochore distributions (determining whether globin loci sit in GC-rich or GC-poor isochores in each clade).

#### 6. Comparative Clade Reporter ([NEW] `globin_pipeline/reporter.py`)
- Aggregates metrics across clades:
  - **Fish** (*danio_rerio, oryzias_latipes, takifugu_rubripes, petromyzon_marinus*)
  - **Amphibians** (*xenopus_tropicalis*)
  - **Amniotes** (*mammals: Homo_sapiens, Bos_taurus, Sus_scrofa, Mus_musculus, Monodelphis_domestica; birds: Gallus_gallus; reptiles: Chrysemys_picta_bellii*)
- Generates:
  - Publication-ready summary tables (`results/summary_globin_loci.csv`)
  - Comparative plots: GC3 vs Genomic GC across clades, gene organization diagrams, sliding window GC profiles.

#### 7. Unified CLI ([NEW] `globin_pipeline/cli.py`)
- Provides command-line entry points:
  ```bash
  # 1. Audit a new batch from Mike
  globin-pipeline audit --batch 1_data_from_mike/20260930/

  # 2. Curate and standardize intervals
  globin-pipeline curate --batch 1_data_from_mike/20260930/ --output data/curated_intervals.csv

  # 3. Run full gene extraction and GC analysis
  globin-pipeline analyze --intervals data/curated_intervals.csv --outdir results/

  # 4. Generate comparative clade report
  globin-pipeline report --input results/ --outdir results/reports/
  ```

---

## Verification Plan

### Automated Tests
- Unit tests with `pytest`:
  - `test_parser.py`: Verify parsing of Mike's YAML, CSVs, and error detection.
  - `test_coordinates.py`: Verify 0-based vs 1-based coordinate conversions and interval arithmetic.
  - `test_annotator.py`: Verify gene extraction and globin classification on test GTF chunks.
  - `test_gc_analyzer.py`: Verify GC, GC3, intron/exon partitioning, and sliding window calculations on known sequence benchmarks.
- Execution command:
  ```bash
  uv run --with pytest,biopython,pandas pytest tests/
  ```

### Manual Verification
- Inspect the curated intervals against Ensembl genome browser for Human, Xenopus, Chicken, and Zebrafish.
- Confirm all known alpha globin genes in Xenopus (*hbz, hba4, hba3, hba1*) are detected inside the revised interval.
- Confirm both Locus A and Locus B are extracted for Zebrafish and Medaka.
- Validate that calculated genomic GC averages match established literature values (e.g. human genome ~41% GC, warm-blooded amniote GC heterogeneity / isochores vs cold-blooded vertebrates).
