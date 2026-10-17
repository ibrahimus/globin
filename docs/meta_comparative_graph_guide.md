# Cross-Species Meta-Comparative Visualization Guide

## 1. Overview & Dual-Graph Architecture

To balance comprehensive genomic completeness with visual clarity, the pipeline generates **two synchronized versions** of the cross-species meta-comparative analysis across all 20 gold-standard vertebrate species:

1. **Comprehensive Version (`meta_comparative_graph.png`, `.pdf`, `.html`)**:
   - Contains **all annotated genomic features** across each locus interval:
     - Functional $\alpha$- and $\beta$-globins
     - Tandem duplicated/ancestral globin forms (Loci A & B in teleosts)
     - Globin pseudogenes ($\psi\alpha$, $\psi\beta$)
     - Syntenic non-globin flanking markers (*RHBDF1, MPG, NPRL3, LUC7L, OR51, CORO1A*, etc.)
     - Other locus features (long non-coding RNAs, uncharacterized transcripts, non-globin pseudogenes)
   - Ideal for genomic audits, boundary verification, and locus completeness inspection.

2. **Clean Version (`meta_comparative_graph_clean.png`, `.pdf`, `.html`)**:
   - Strictly focuses on **functional protein-coding genes**:
     - $\alpha$-globin coding genes (circles, blue)
     - $\beta$-globin coding genes (squares, green)
     - Syntenic non-globin flankers from the $\alpha$-locus (down-triangles, orange)
     - Syntenic non-globin flankers from the $\beta$-locus (left-triangles, brown)
   - Filters out clutter from non-coding locus features and pseudogenes to highlight evolutionary trends.

---

## 2. Panels and Compositional Metrics

Each figure presents two side-by-side comparative panels standardized against each species' own whole-genome baseline:

- **Panel A: Whole-Gene Body GC Z-Score ($Z_{\text{Full Gene}}$)**:
  - Computed from the complete transcribed genomic unit (all exons, introns, and UTRs) measured directly against the species' genome-wide whole-gene distribution.
  - Crucially incorporates **intronic base composition**. In lineages with compact GC-rich introns (*Pogona vitticeps*, *Ornithorhynchus anatinus*), whole-gene GC is elevated; in lineages with expanded, AT-rich introns (mammalian $\beta$-globins), whole-gene GC is lower than CDS GC.
- **Panel B: Coding Sequence GC Z-Score ($Z_{\text{CDS}}$)**:
  - Computed strictly from the functional open reading frame (CDS), reflecting both synonymous codon usage (GC3) and amino acid composition under purifying selection.

---

## 3. Deterministic Y-Axis Stratification (Offset Bands)

Both plots use **deterministic, tiered y-offsets** within each species' horizontal band (`y = base_index ± 0.35`). This guarantees identical horizontal tracking across runs, reproducible layouts, and clear alignment between Whole-Gene (Panel A) and CDS (Panel B).

### Tier Mapping

| Tier | Category | Shape / Color | Offset (Clean Mode) | Offset (Full Mode) |
| :--- | :--- | :--- | :---: | :---: |
| **Tier 1** | $\alpha$-Globins, Locus A, Ancestral linked | Blue Circle / Triangle | **`-0.18`** | **`-0.22`** |
| **Tier 2** | $\beta$-Globins, Locus B | Green Square / Triangle | **`-0.06`** | **`-0.11`** |
| **Tier 3** | Globin Pseudogenes ($\psi$) | Purple Diamond | *Filtered out* | **`0.00`** |
| **Tier 4** | Synteny Non-Globins ($\alpha$-locus / linked) | Orange Down-Triangle | **`+0.06`** | **`+0.11`** |
| **Tier 5** | Synteny Non-Globins ($\beta$-locus) | Brown Left-Triangle | **`+0.18`** | **`+0.22`** |
| **Tier 6** | Other Locus Features (lncRNAs, misc) | Gray Circle | *Filtered out* | **`+0.31`** |

### Micro-Staggering for Clustered Genes
When multiple genes within the same tier share nearly identical GC Z-scores ($|\Delta Z| < 0.18$), an alternating micro-stagger ($\pm 0.025$) is applied to prevent marker occlusion while preserving tier separation.

---

## 4. Key Comparative Takeaways

1. **Mammalian $\alpha$- vs. $\beta$-Globin Separation**:
   - In Placental Mammals and Turtles, $\alpha$-globins occupy the high GC right tail ($Z > +1.5\sigma$ to $+3.5\sigma$), while $\beta$-globins sit near or below genome average ($Z \approx -0.5\sigma$ to $+0.5\sigma$).
2. **Avian and Reptilian Subtelomeric & Microchromosomal Drivers**:
   - In Birds (Duck, Turkey, Chicken), Crocodilians (Alligator), and Squamates (Anole, Bearded Dragon), $\beta$-globins match or exceed $\alpha$-globins in GC content. In birds, the distal location of the $\beta$-cluster on Macrochromosome 1 (<150 kb–1.5 Mb from chromosome terminus) subjects it to elevated crossover rates and GC-biased gene conversion (gBGC).
3. **Platypus Whole-Gene Inversion**:
   - In Platypus (*O. anatinus*), $\beta$-globins have higher whole-gene GC than $\alpha$-globins (Panel A), driven by their GC-rich intronic architecture (55.9% intronic GC), reversing the pattern observed in therian mammals.
4. **Squamate Intronic Elevation**:
   - In Bearded Dragon (*P. vitticeps*), both $\alpha$ and $\beta$ globins and their syntenic flankers reside in elevated GC microchromosomes (>66% whole-gene GC, >70% intronic GC, >85% GC3), demonstrating high-GC equilibrium across both coding and non-coding sequence.
5. **Teleost Locus A/B Symmetry**:
   - Teleost fish (Zebrafish, Medaka) retain linked $\alpha$- and $\beta$-globins co-located in shared genomic clusters, showing synchronized compositional evolution.

---

## 5. File References

- **Clean Graph (Publication-Ready)**: [`results/meta_comparative_graph_clean.png`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph_clean.png)
- **Comprehensive Graph**: [`results/meta_comparative_graph.png`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph.png)
- **Vector Graphics**: [`results/meta_comparative_graph_clean.pdf`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph_clean.pdf)
- **Interactive Tooltip HTML**: [`results/meta_comparative_graph.html`](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph.html)
- **Standardized Dataset**: [`results/gene_zscores.csv`](file:///Users/ibrahim/Documents/globin_project/results/gene_zscores.csv)
