# Walkthrough: Empirical Gene GC Standardization & Globin GC Calculation Methodology

This walkthrough details the methodology for calculating globin gene GC content across genomic compartments, documents the transition to **empirical standard deviations and baselines derived directly from all protein-coding genes**, and summarizes the updated cross-species results across **24 vertebrate species (41 curated intervals, 408 genes, 92 functional globins)**.

---

## 1. How Globin Gene GC Content Is Calculated

In [`GCAnalyzer.analyze_locus`](file:///Users/ibrahim/Documents/globin_project/globin_pipeline/gc_analyzer.py#L78-L246), GC content is calculated across **five distinct biological compartments** for every globin gene:

```
┌──────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                GLOBIN GENE COMPARTMENT BREAKDOWN                                 │
├──────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Genomic Locus DNA: 5' ───[Exon 1]───(Intron 1)───[Exon 2]───(Intron 2)───[Exon 3]─── 3'          │
│                                                                                                  │
│ 1. Full Gene Body: ├───────────────────────── Entire Span ────────────────────────┤ (Exons+Introns)│
│ 2. Coding CDS:           ├───CDS 1──┤           ├──CDS 2──┤           ├──CDS 3──┤   (Spliced Coding)│
│ 3. Codon GC3:            [..X ..X]               [..X ..X]             [..X ..X]    (3rd Bases)     │
│ 4. Exonic GC:      ├──Exon 1──┤                 ├──Exon 2─┤           ├──Exon 3─┤   (cDNA Transcr.) │
│ 5. Intronic GC:                 ├────Intron 1───┤         ├────Intron 2───┤         (Introns only)  │
└──────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### A. Coding Sequence GC (`cds_gc`)
- **Sequence Retrieval**: The canonical protein-coding transcript is queried via the Ensembl REST API endpoint:
  ```
  /sequence/id/{transcript_id}?type=cds
  ```
  This returns the spliced coding sequence (only coding exons, from ATG start codon to stop codon; introns and 5'/3' UTRs are completely excluded).
- **Formula**:
  $$\text{CDS\_GC} = \frac{N_G + N_C}{N_A + N_C + N_G + N_T} \times 100$$

### B. Third Codon Position GC (`gc3`)
- **Extraction**: From the canonical CDS sequence, codons are parsed in-frame ($0, 1, 2 \mid 3, 4, 5 \mid \dots$). Bases at 0-indexed positions `2, 5, 8, 11, ...` are extracted.
- **Formula**:
  $$\text{GC3} = \frac{N_{G_3} + N_{C_3}}{N_{\text{valid third codon bases}}} \times 100$$
- Reflects synonymous mutational pressure and GC-biased gene conversion (gBGC), independent of amino acid selection.

### C. Full Gene Body GC (`full_gene_gc`)
- **Sequence Retrieval**: Sliced directly from the unmasked chromosomal locus sequence (`locus_seq`) using the gene's chromosomal coordinates (`start` to `end`):
  $$\text{gene\_genomic\_seq} = \text{locus\_seq}[\text{rel\_start} : \text{rel\_end}]$$
- **Content**: Encompasses the entire genomic transcription unit: **5' UTR + coding exons + all introns + 3' UTR**.
- **Formula**:
  $$\text{Full\_Gene\_GC} = \frac{\text{Count}(G) + \text{Count}(C)}{\text{Count}(A) + \text{Count}(C) + \text{Count}(G) + \text{Count}(T)} \times 100$$

### D. Exonic GC (`exonic_gc`)
- Spliced transcript model derived by concatenating all exon coordinates (`target_transcript["Exon"]`), capturing coding exons and UTRs.

### E. Intronic GC (`intronic_gc`)
- Intervening sequences between consecutive exons within the gene body.

---

## 2. Why Whole-Gene GC Differs from CDS GC

A classic vertebrate globin gene consists of **3 coding exons and 2 introns**. The coding sequence is conserved at ~426 to 444 bp across vertebrates, whereas introns can span from a few hundred base pairs to several kilobases.

In most non-mammalian vertebrates and mammalian $\beta$-globins, **introns are markedly AT-rich**. Because introns comprise >85% of the total gene length, `full_gene_gc` is pulled toward the AT-rich intronic value, while `cds_gc` reflects coding constraints:

| Taxon / Gene | Full Gene GC (`full_gene_gc`) | Coding CDS GC (`cds_gc`) | Intronic GC (`intronic_gc`) | Biological Context |
| :--- | :---: | :---: | :---: | :--- |
| **Human *HBB*** ($\beta$) | **37.64%** | **56.08%** | **31.40%** | AT-rich introns pull full gene body down |
| **Human *HBA1*** ($\alpha$) | **67.86%** | **64.80%** | **69.80%** | Ultradense GC-rich isochore elevates whole locus |
| **Zebrafish *hbae5*** | **31.52%** | **53.01%** | **29.22%** | Coding exons are +0.76$\sigma$, but introns are ~29% GC |
| **Zebrafish *hbbe2*** | **45.23%** | **53.83%** | **30.29%** | Coding exons are +0.93$\sigma$, introns are ~30% GC |
| **Frog *hba3*** (*Xenopus*) | **33.90%** | **50.12%** | **31.80%** | Coding exons are +0.57$\sigma$, introns are ~32% GC |
| **Turtle *ENSCPBG00000021826*** | **42.70%** | **49.12%** | **41.30%** | Intronic GC reduces full gene body GC |
| **Chicken *HBAD*** ($\alpha$) | **61.54%** | **58.65%** | **63.40%** | Avian microchromosome: both exons and introns high GC |

---

## 3. Transition to Empirical Protein-Coding Gene Standard Deviations

### What Was Changed
1. **Elimination of Arbitrary Assembly Baselines**:
   - Previously, Whole-Gene GC was standardized against hardcoded 100 kb genomic assembly window estimates (`GENOMIC_ASSEMBLY_BASELINES`).
   - In [`ZScoreAnalyzer`](file:///Users/ibrahim/Documents/globin_project/globin_pipeline/zscore_analyzer.py#L148-L153), Whole-Gene GC ($Z_{\text{Full\_Gene}}$) is now standardized **strictly against the empirical distribution of all protein-coding genes** ($\mu_{\text{gene\_emp}}, \sigma_{\text{gene\_emp}}$) in that species.
2. **Empirical Standard Deviation Calculation**:
   - For every species ($N \approx 10,415$ to $30,313$ protein-coding genes streamed from Ensembl), the empirical distribution is computed directly:
     $$\mu_{\text{empirical}} = \frac{1}{N} \sum_{i=1}^N \text{Gene\_GC}_i, \quad \sigma_{\text{empirical}} = \sqrt{\frac{1}{N-1} \sum_{i=1}^N (\text{Gene\_GC}_i - \mu_{\text{empirical}})^2}$$
3. **Standardized Z-Score Formulas**:
   $$Z_{\text{CDS}} = \frac{\text{cds\_gc} - \mu_{\text{empirical\_cds}}}{\sigma_{\text{empirical\_cds}}}$$
   $$Z_{\text{Full\_Gene}} = \frac{\text{full\_gene\_gc} - \mu_{\text{empirical\_gene}}}{\sigma_{\text{empirical\_gene}}}$$
   $$Z_{\text{GC3}} = \frac{\text{gc3} - \mu_{\text{empirical\_gc3}}}{\sigma_{\text{empirical\_gc3}}}$$

### Concrete Empirical Baselines Across 24 Species

| Species | Protein-Coding Genes ($N$) | Empirical Coding Mean ($\mu_{\text{cds}}$) | Empirical Coding Std ($\sigma_{\text{cds}}$) | Empirical GC3 Std ($\sigma_{\text{gc3}}$) |
| :--- | :---: | :---: | :---: | :---: |
| **Danio rerio** (Zebrafish) | 30,313 | 49.29% | **4.87%** | 9.40% |
| **Xenopus tropicalis** (Frog) | 22,107 | 46.89% | **5.65%** | 10.82% |
| **Chrysemys picta bellii** (Painted Turtle) | 21,186 | 51.37% | **8.29%** | 16.21% |
| **Pelodiscus sinensis** (Softshell Turtle) | 18,189 | 49.46% | **8.21%** | 16.44% |
| **Crocodylus porosus** (Crocodile) | 15,226 | 49.59% | **7.75%** | 15.36% |
| **Gallus gallus** (Chicken) | 17,007 | 53.12% | **8.59%** | 16.03% |
| **Anas platyrhynchos** (Mallard) | 16,765 | 53.25% | **9.34%** | 17.82% |
| **Pogona vitticeps** (Bearded Dragon) | 19,580 | 49.59% | **6.96%** | 13.82% |
| **Mus musculus** (Mouse) | 22,311 | 51.71% | **6.95%** | 11.71% |
| **Sus scrofa** (Pig) | 22,041 | 53.08% | **9.28%** | 16.47% |
| **Bos taurus** (Cattle) | 20,848 | 53.76% | **8.69%** | 15.88% |
| **Homo sapiens** (Human) | 23,879 | 53.41% | **8.59%** | 15.63% |

---

## 4. Visualizations

### Cross-Species Meta-Comparative Graph
Both panels are now standardized to **empirical protein-coding gene distributions**:

![Cross-Species Meta-Comparative Graph](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/meta_comparative_graph.png)

- **Vector PDF Output**: [meta_comparative_graph.pdf](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph.pdf)
- **Interactive Plotly Dashboard**: [meta_comparative_graph.html](file:///Users/ibrahim/Documents/globin_project/results/meta_comparative_graph.html)
- **Updated Dataset**: [gene_zscores.csv](file:///Users/ibrahim/Documents/globin_project/results/gene_zscores.csv) (408 genes, 92 globins across 24 species).

---

## 5. Redesigned 5-Panel Diagnostic Plots: Empirical Whole-Gene vs. CDS GC Density Distributions

In response to the user's request, the diagnostic graphs have been expanded from a 4-panel 2x2 layout to a publication-grade **5-panel GridSpec layout (3x2)** in both static (`matplotlib`) and interactive (`plotly`) engines:

```
┌───────────────────────────────────────────────┬───────────────────────────────────────────────┐
│ Panel A: Whole-Gene GC Density (Exons+Introns)│ Panel B: Coding Sequence (CDS) GC Density     │
│  - All ~20,000 protein-coding genes           │  - All ~20,000 protein-coding genes           │
│  - Mean (μ_whole), Assembly (μ_gen), ±1σ, ±2σ │  - Mean (μ_cds), Assembly (μ_gen), ±1σ, ±2σ   │
│  - Globin & locus genes marked & labeled      │  - Globin & locus genes marked & labeled      │
│  [Shared x-axis scale min/max for comparison] │  [Shared x-axis scale min/max for comparison] │
├───────────────────────────────────────────────┼───────────────────────────────────────────────┤
│ Panel C: Whole-Gene GC vs. CDS GC             │ Panel D: CDS GC vs. GC3 Codon Space           │
│  - Isochore divergence vs. coding selection   │  - Selection vs. mutational bias/gBGC         │
│  - y = x diagonal line, genomic/CDS baselines │  - Empirical bivariate contour landscape      │
├───────────────────────────────────────────────┴───────────────────────────────────────────────┤
│ Panel E: Compartment Comparison (Spanning Bottom)                                             │
│  - Paired boxplots & stripplots for Whole-Gene GC (blue) vs. CDS GC (orange)                 │
│  - Genome vs Non-Globin (α/β) vs Pseudogenes vs Alpha/Beta Globin Coding                      │
└───────────────────────────────────────────────────────────────────────────────────────────────┘
```

### Key Design Enhancements:
1. **Side-by-Side Whole-Gene vs. CDS Density Curves (Panels A & B)**:
   - **Panel A** plots the empirical distribution of `full_gene_gc` across all ~20,000 protein-coding genes (exons + introns) in the species.
   - **Panel B** plots the empirical distribution of `cds_gc` across all ~20,000 protein-coding genes.
   - **Globin Placement**: Globin genes (and locus non-globin bystander genes) are explicitly marked with large symbols on both density plots.
   - **Aligned X-Axes**: The horizontal axes of Panels A and B share identical minimum and maximum limits, visually demonstrating the genome-wide rightward shift from whole-gene to coding CDS, and highlighting how individual globins shift relative to their genomic background.
   - **Anti-Collision Label Staggering**: Uses 1D cluster detection and 4-tier vertical offsets (`0.32×`, `0.54×`, `0.76×`, `0.96× y_max`) so tandemly duplicated globins with nearly identical GC values (e.g., *HBA1* vs *HBA2*, *HBG1* vs *HBG2*) remain legible without overlapping arrows.
2. **Bottom Spanning Compartment Comparison (Panel E)**:
   - Spanning the full 2-column width (`gs[2, :]`), allowing ample horizontal space for multi-category boxplots and jittered strip markers without label clipping.

### Representative Diagnostic Figures

#### Human (*Homo sapiens*)
- Whole-gene baseline: $\mu = 46.7\%$, assembly baseline: $41.0\%$, CDS baseline: $\mu = 53.1\%$.
- Shows the striking bifurcation: $\beta$-globins sit in an AT-rich whole-gene environment (~38–46% GC) with CDS around 50–56%, while $\alpha$-globins reside in an extreme GC-rich isochore (>60–72% whole-gene GC, >60–70% CDS GC).

![Human Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_Homo_sapiens.png)

#### Zebrafish (*Danio rerio*)
- Whole-gene baseline: $\mu = 37.3\%$, assembly baseline: $36.5\%$, CDS baseline: $\mu = 49.4\%$.
- Illustrates how zebrafish globins on Locus A and Locus B have whole-gene GC values clustered between 32% and 48%, but coding sequence GC values tightly centered around 50–56%.

![Zebrafish Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_danio_rerio.png)

#### Painted Turtle (*Chrysemys picta bellii*)
- Whole-gene baseline: $\mu = 48.0\%$, assembly baseline: $44.0\%$, CDS baseline: $\mu = 51.4\%$.
- Shows *HBA* and *HBAD* elevated into the high-GC tail (>65% whole-gene GC, >66% CDS GC), while *HBB_like* genes remain close to the genome mean.

![Turtle Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_Chrysemys_picta_bellii.png)

---

## 6. Resolution of Missing Whole-Gene GC Distributions & Assembly Coordinates

### Why Whole-Gene GC Was Missing in Four Species
1. **European Hedgehog (*Erinaceus europaeus*)**:
   - The assembly baseline cache (`erinaceus_europaeus_genes_gc.csv.gz`) originated from an early CDS-only extraction script that queried `/sequence/id/` and lacked genomic span coordinates (`whole_gene_gc`).
   - In addition, because the draft assembly (`HEDGEHOG` / `EriEur2.0`) consists of highly fragmented unplaced scaffolds with interrupted globin synteny, Hedgehog is excluded from active curated intervals in `cli.py`.
2. **Spotted Gar (*Lepisosteus oculatus*)**:
   - The species manifest previously specified coordinates on `LG19: 22,100,000–22,450,000`, which was out of bounds (Gar LG19 is only 9.7 Mb in `LepOcu1`).
   - The true intact tandem globin cluster resides on **`LG13: 14,650,000–14,900,000`** (250 kb span flanked by `rhbdf1a` and `luc7l`), harboring linked $\alpha$ (`HBA2`, `hbae5`, `hbae4`) and $\beta$ (`HBE1`) globins.
3. **Platypus (*Ornithorhynchus anatinus*)**:
   - Manifest coordinates referenced legacy unplaced contigs (`Ultra23` and `Ultra5`).
   - In the modern reference assembly `mOrnAna1.p.v1`, the $\alpha$-locus is located on **`Chr 21: 15,600,000–15,860,000`** (`LUC7L` to `RHBDF1`, containing `HBA` and `HBM`), and the $\beta$-locus on **`Chr 2: 139,350,000–139,550,000`** (flanked by olfactory receptors `OR52D1` and `OR51E1`, harboring `HBE`, `HBB`, and `HBE_like`).
4. **Sea Lamprey (*Petromyzon marinus*)**:
   - A typographical error in the scaffold identifier (`GL476331` vs true scaffold **`GL478636: 20,000–150,000`**) previously caused locus sequence retrieval to fail. The cluster contains cyclostome proto-globins `pm_globin1` and `pm_globin2` flanked by `UBE2I` and `nprl3`.

### Extracted Values for Remediated Species

| Species | Gene | Compartment | GC Content | Empirical Gene Mean | Empirical Gene $\sigma$ | Z-Score |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| **Spotted Gar** | *HBA2* | Whole-Gene / CDS | 50.39% / 63.85% | 43.73% / 53.50% | 7.22% / 7.29% | $+0.92$ / $+1.42$ |
| | *hbae5* | Whole-Gene / CDS | 55.08% / 60.92% | 43.73% / 53.50% | 7.22% / 7.29% | $+1.57$ / $+1.02$ |
| | *hbae4* | Whole-Gene / CDS | 47.78% / 59.59% | 43.73% / 53.50% | 7.22% / 7.29% | $+0.56$ / $+0.84$ |
| | *HBE1* | Whole-Gene / CDS | 51.59% / 62.39% | 43.73% / 53.50% | 7.22% / 7.29% | $+1.09$ / $+1.22$ |
| **Platypus** | *HBA* | Whole-Gene / CDS | 54.09% / 63.96% | 52.71% / 57.47% | 8.71% / 9.11% | $+0.16$ / $+0.71$ |
| | *HBM* | Whole-Gene / CDS | 50.90% / 62.44% | 52.71% / 57.47% | 8.71% / 9.11% | $-0.21$ / $+0.55$ |
| | *HBE* | Whole-Gene / CDS | 58.14% / 60.81% | 52.71% / 57.47% | 8.71% / 9.11% | $+0.62$ / $+0.37$ |
| | *HBE_like* | Whole-Gene / CDS | 57.43% / 57.89% | 52.71% / 57.47% | 8.71% / 9.11% | $+0.54$ / $+0.05$ |
| | *HBB* | Whole-Gene | 56.37% | 52.71% | 8.71% | $+0.42$ |
| **Sea Lamprey**| *pm_globin1* | Whole-Gene / CDS | 47.68% / 62.69% | 48.99% / 59.66% | 6.14% / 6.10% | $-0.21$ / $+0.50$ |
| | *pm_globin2* | Whole-Gene / CDS | 49.40% / 61.84% | 48.99% / 59.66% | 6.14% / 6.10% | $+0.07$ / $+0.36$ |

---

## 7. Deterministic Y-Axis Stratification in Diagnostic Density Plots

To eliminate overlapping marker dots while ensuring strict cross-species comparability, a **deterministic vertical stratification scheme** was implemented in [`Visualizer._calculate_gene_y_positions`](file:///Users/ibrahim/Documents/globin_project/globin_pipeline/visualizer.py#L150-L215):

```
Density Max (y_max)
  ▲
  │   [Tier 6] ~0.29 × y_max: Other locus features (e.g. flanking markers)
  │   [Tier 5] ~0.24 × y_max: Non-globin bystander genes from β-locus (or locus B)
  │   [Tier 4] ~0.19 × y_max: Non-globin bystander genes from α-locus (or locus A / linked)
  │   [Tier 3] ~0.14 × y_max: Globin pseudogenes
  │   [Tier 2] ~0.09 × y_max: Beta-globin coding genes (or locus B globins)
  │   [Tier 1] ~0.04 × y_max: Alpha-globin coding genes (or locus A / linked globins) [CLOSEST TO 0]
──┴─────────────────────────────────────────────────────────────────────────────────────────────► GC (%)
```

### Key Technical Properties:
1. **Consistent Between All Graphs**: The vertical baseline for each gene category is strictly fixed relative to $y_{\max}$, making cross-species comparison instantaneous.
2. **Micro-Staggering for Clustered Duplications**: For tandemly duplicated paralogs within the same tier whose GC values differ by $<0.8\%$ (e.g., *HBA1* vs *HBA2* or *HBG1* vs *HBG2*), a deterministic alternating vertical micro-offset ($\pm 0.015 \times y_{\max}$) is applied so dots never occlude one another.
3. **Exact Arrow Pointer Alignment**: Callout annotations (`_annotate_density_points`) point arrow tips directly to $(x, y_{\text{dot}})$ on their respective stratified tiers, with callout text positioned in upper display tiers ($0.36 \times$ to $0.96 \times y_{\max}$).
4. **Parity Across Engines**: Both static publication figures (`matplotlib`) and interactive Plotly HTML dashboards share identical stratification logic.

### Remediated Diagnostic Plots:

#### Spotted Gar (*Lepisosteus oculatus*)
![Gar Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_lepisosteus_oculatus.png)

#### Platypus (*Ornithorhynchus anatinus*)
![Platypus Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_ornithorhynchus_anatinus.png)

#### Sea Lamprey (*Petromyzon marinus*)
![Sea Lamprey Diagnostic Plot](/Users/ibrahim/.gemini/antigravity-cli/brain/9a1c7782-8799-45ce-a82d-191b5d492881/diagnostic_petromyzon_marinus.png)


