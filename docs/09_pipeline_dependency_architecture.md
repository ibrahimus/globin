# Chapter 09: Pipeline Dependency Architecture & Targeted Species Update Protocol

## 1. Executive Summary

This document establishes the formal dependency hierarchy, automated propagation rules, and targeted update protocol for the Vertebrate Globin Synteny & Isochore Dynamics Project.

To ensure both **computational efficiency** (avoiding costly re-downloads or re-computations of all 68 vertebrate genomes) and **biological/data integrity** (guaranteeing that every upstream update immediately cascades into all downstream tables, publication figures, web dashboards, and test suites), the pipeline is structured into **5 distinct dependency tiers**.

---

## 2. The 5-Tier Dependency Hierarchy

```mermaid
graph TD
    classDef l0 fill:#eceff1,stroke:#607d8b,stroke-width:2px;
    classDef l1 fill:#e3f2fd,stroke:#1976d2,stroke-width:2px;
    classDef l2 fill:#e8f5e9,stroke:#388e3c,stroke-width:2px;
    classDef l3 fill:#fff8e1,stroke:#f57f17,stroke-width:2px;
    classDef l4 fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px;

    subgraph Level0 ["Level 0: Single Source of Truth & Manifests"]
        M[config/species_manifest.yaml]:::l0
        GB[data/genome_baselines/<species>_genes_gc.csv.gz]:::l0
        ST[data/genome_baselines/<species>_stats.json]:::l0
    end

    subgraph Level1 ["Level 1: Targeted Ingestion & Per-Species Caching"]
        FA[data/cache/seq_<species>_<chr>_<start>_<end>.fa]:::l1
        GJ[data/cache/genes_<species>_<chr>_<start>_<end>.json]:::l1
        SJ[data/cache/gene_<gid>.json]:::l1
        P1[results/chromosome_plots/chromosome_linear_<species>.png]:::l1
        P2[results/bivariate_length_gc_plots/bivariate_length_gc_<species>.png]:::l1
        P3[results/bivariate_with_inset_plots/bivariate_with_inset_<species>.png]:::l1
    end

    subgraph Level2 ["Level 2: Global Table Compilation"]
        CI[results/curated_intervals.csv & .bed]:::l2
        GZ[results/gene_zscores.csv]:::l2
        CS[results/cross_species_summary_table.csv]:::l2
        S1[results/Supplementary_Table_S1_Species_Baselines.csv]:::l2
        S2[results/Supplementary_Table_S2_Globin_Loci_Summary.csv]:::l2
        S3[results/Supplementary_Table_S3_All_1497_Genes_Annotations_ZScores.csv]:::l2
        S4[results/Supplementary_Table_S4_Statistical_Contrasts.csv]:::l2
    end

    subgraph Level3 ["Level 3: Macro Figures & Analytical Models"]
        F0[results/summary_figure0_overview_and_methodology.png/pdf]:::l3
        F1[results/summary_figure1a_locus_divergence_summary.png/pdf]:::l3
        F1B[results/summary_figure1b_locus_divergence_full.png/pdf]:::l3
        F2[results/summary_figure2_mechanistic_drivers.png/pdf]:::l3
        F3[results/summary_figure3_phylogenetic_synteny.png/pdf]:::l3
        F4[results/summary_figure_bernardi_compaction_law.png/pdf]:::l3
        CH[results/codon_decomposition_per_gene.csv]:::l3
    end

    subgraph Level4 ["Level 4: Publication Artifacts, Dashboard & Audits"]
        DB[results/dashboard.html]:::l4
        AUD[scripts/verify_cluster_completeness.py]:::l4
        TEST[tests/ (pytest)]:::l4
    end

    M -->|Targeted Fetch| FA
    M -->|Targeted Parse| GJ
    FA & GJ --> SJ
    SJ --> P1 & P2 & P3
    SJ --> GZ
    M --> CI
    GZ & CI --> CS & S1 & S2 & S3 & S4
    GZ & CS --> F1 & F1B & F2
    GZ & CH --> F3 & F4
    GZ & CI & CS & F0 & F1 & F2 & F3 & F4 --> DB
    CI & GZ & CS --> AUD
    DB & GZ & CI --> TEST
```

---

## 3. Targeted Species Update Protocol

When an individual species is upgraded (e.g. replacing a fragmented assembly with PacBio HiFi long-read data or fixing unannotated globin gene models):

### Step 1: Targeted Ingestion (Local Scope)
Only the modified species is queried against NCBI:
1. Update `config/species_manifest.yaml` with the new assembly accession and locus boundaries.
2. Fetch the sequence slice for only the updated locus via NCBI E-utilities (takes ~1-2 seconds) $\rightarrow$ `data/cache/seq_<species>_<chr>_<start>_<end>.fa`.
3. Stream the specific scaffold features from RefSeq GFF3 $\rightarrow$ `data/cache/genes_<species>...json` and `data/cache/gene_<gid>.json`.
4. Render only the 4 single-species plots for that species:
   ```bash
   python scripts/regenerate_plots_with_inset.py <species>
   ```
*Duration: ~10 seconds. Zero recomputation of the other 67 genomes.*

### Step 2: Global Table Compilation
1. Update `results/curated_intervals.csv` from `config/species_manifest.yaml`.
2. Update `results/gene_zscores.csv`: replace only the affected species' rows with the newly computed Z-scores.
3. Recompile the Supplementary Data release package:
   ```bash
   python scripts/build_supplementary_data_package.py
   ```
*Duration: ~3 seconds.*

### Step 3: Cascaded Macro Figures & Statistical Models
1. Re-run `scripts/build_evolutionary_history_and_synteny.py` to refresh codon frequencies and ancestral synteny flags.
2. Re-run `scripts/generate_all_summary_figures.py` (Fig 1A, 1B, Fig 2).
3. Re-run `scripts/generate_figure3_phylogenetic_synteny.py` (Fig 3).
4. Re-run `scripts/analyze_bernardi_compaction.py` (Fig 4 & Bernardi models).
*Duration: ~25 seconds.*

### Step 4: Standalone Dashboard Build & Verification
1. Re-run `scripts/generate_static_dashboard.py` to embed the new tables and link updated plots into `results/dashboard.html`.
2. Run `scripts/verify_cluster_completeness.py` to assert 100% amniote dual-cluster completeness, teleost 3R integrity, and coordinate validity.
3. Run `pytest tests` to verify regression-free state.
*Duration: ~5 seconds.*

**Total turnaround time for a complete, verified, targeted species update: < 45 seconds.**

---

## 4. Master Reproduction Command

To reproduce the entire analytical atlas, all 5 figures, the statistical models, and the web portal in one command:

```bash
python scripts/reproduce_all.py
```

This single command executes the full Level 2 $\rightarrow$ Level 4 cascade with automated pass/fail reporting.
