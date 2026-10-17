# Feasibility Assessment: 8 Candidate Vertebrate Genomes

This report evaluates the technical feasibility, assembly quality, globin syntenic architecture, and macroevolutionary value of the **8 candidate genomes** queried from NCBI RefSeq.

---

## 1. Executive Summary & Tiered Recommendations

| Accession | Species & Common Name | Clade | Assembly Level | Genome Size | Contiguity (Scaf N50) | Globin Synteny Status | Feasibility Tier | Recommendation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GCF_965207765.1** | *Zootoca vivipara* (Viviparous Lizard) | Squamata (Lacertidae) | **Chromosome** | 1.42 Gb | 95.2 Mb | $\alpha$ on NC_142596.1, $\beta$ on NC_142587.1 | **Tier 1 (High Priority)** | **Integrate**: Directly tests whether *Podarcis* hyper-gBGC ($Z > +5.6\sigma$) is family-wide. |
| **GCF_056783735.1** | *Eptatretus burgeri* (Inshore Hagfish) | Cyclostomata (Myxinidae) | **Chromosome** | 2.00 Gb | 108.8 Mb | Proto-cluster on NC_140977.1 (12 globins + NPRL3) | **Tier 1 (High Priority)** | **Integrate**: Definitive jawless craniate outgroup; resolves pre-1R/2R globin origin. |
| **GCF_964198595.1** | *Lampetra fluviatilis* (European River Lamprey) | Cyclostomata (Petromyzontidae) | **Chromosome** | 1.04 Gb | 13.2 Mb | Proto-cluster on NC_136297.1 (17 globins + NPRL3) | **Tier 1 (High Priority)** | **Integrate**: Upgrades fragmented *Petromyzon* to a modern Darwin Tree of Life chromosome build. |
| **GCF_049634645.1** | *Pelodiscus sinensis* (Chinese Softshell Turtle) | Testudines (Trionychidae) | **Chromosome** | 2.24 Gb | 129.6 Mb | $\alpha$ on NC_134726.1, $\beta$ on NC_134711.1 | **Tier 1 (High Priority)** | **Integrate**: Deepest cryptodiran turtle branch (~150 Mya divergence from Emydidae/Cheloniidae). |
| **GCF_964187855.1** | *Myxine glutinosa* (Atlantic Hagfish) | Cyclostomata (Myxinidae) | **Chromosome** | 3.06 Gb | 192.2 Mb | Proto-cluster on NC_136265.1 (9 globins + NPRL3) | **Tier 1 (Alternative)** | Excellent alternative/companion to *E. burgeri* for hagfish cyclostome baseline. |
| **GCF_059951505.1** | *Python bivittatus* (Burmese Python) | Squamata (Pythonidae) | Contig / Scaffold | 1.43 Gb | 57.1 Mb | $\alpha$ on NW_028277083.1, $\beta$ on NW_028277056.1 | **Tier 2 (Moderate)** | Basal snake anchor (Henophidia), but lacks full chromosome-arm numbering. |
| **GCF_031769625.1** | *Taudactylus pleione* (Kroombit Tinker Frog) | Amphibia (Myobatrachidae) | Contig / Scaffold | 5.51 Gb | 8.85 Mb | Linked $\alpha+\beta$ on NW_028284468.1 (21 globins) | **Tier 2 (Heavy / Scaffold)** | Intact linked cluster, but 5.5 Gb genome and scaffold-level assembly. |
| **GCF_965140235.1** | *Caretta caretta* (Loggerhead Sea Turtle) | Testudines (Cheloniidae) | **Chromosome** | 2.22 Gb | 137.5 Mb | $\alpha$ on NC_134215.1, $\beta$ on NC_134206.1 | **Tier 3 (Redundant)** | Highly contiguous, but functionally identical to *C. mydas* and *E. imbricata* (3rd sea turtle). |

---

## 2. Detailed Evaluation of High-Priority Candidates (Tier 1)

### 1. *Zootoca vivipara* (Viviparous / Common Lizard, GCF_965207765.1)
- **Assembly Quality**: Chromosome-level (`rZooViv2.hap1.1`), 1.42 Gb, Scaffold N50 = 95.2 Mb. RefSeq Annotation Release: September 2026 (20,782 protein-coding genes).
- **Locus Synteny**:
  - **$\alpha$-Globin Locus (`NC_142596.1`: 53,440,000–53,620,000, ~180 kb)**:
    - Anchors: `RHBDF1` (53.56 Mb) — `MPG` (53.59 Mb) — `NPRL3` (53.60 Mb).
    - Globin Genes: `LOC118093266` ($\alpha$), `LOC118093265` ($\alpha-2$), `LOC118093251` ($\alpha-D$).
  - **$\beta$-Globin Locus (`NC_142587.1`: 101,700,000–102,850,000, ~1.15 Mb)**:
    - Globin Genes: 6 tandem adult $\beta$-globins (`LOC132592081`, `LOC118085601`, `LOC118085062`, `LOC118079427`, `LOC118078484`, `LOC118078305`).
- **Macroevolutionary Value**: **Exceptional**. In our current 40-species dataset, *Podarcis muralis* (common wall lizard, also Lacertidae) is the single most extreme GC outlier in the vertebrate kingdom ($\alpha \text{ locus } Z = +5.69\sigma$, $\beta \text{ locus } Z = +5.01\sigma$). Adding *Zootoca vivipara* provides the definitive phylogenetic test of whether extreme subtelomeric gBGC is a private idiosyncrasy of *Podarcis* or a shared synapomorphy across the family Lacertidae.

### 2. *Eptatretus burgeri* (Inshore Hagfish, GCF_056783735.1)
- **Assembly Quality**: Chromosome-level (`Ebu_testis_v1`), 2.00 Gb, Scaffold N50 = 108.8 Mb. RefSeq Annotation Release: July 2026 (19,022 protein-coding genes).
- **Locus Synteny**:
  - **Proto-Globin Locus (`NC_140977.1`: 52,080,000–52,280,000, ~200 kb)**:
    - Contains `NPRL3` (52.08–52.15 Mb) directly adjacent to a tandem array of 12 globin genes (`LOC147010100`, `LOC147010099`, `LOC147010098`, etc.).
    - Flanking syntenic genes `RHBDF1` and `MPG` reside on chromosome `NC_140966.1`.
- **Macroevolutionary Value**: **Monumental**. Hagfishes (Myxiniformes) and lampreys (Petromyzontiformes) together comprise Cyclostomata, the living sister group to all jawed vertebrates (Gnathostomata). Hagfish diverged prior to or at the earliest stages of vertebrate whole-genome duplications (1R/2R). Having a modern, chromosome-level hagfish establishes the absolute ancestral baseline for the vertebrate globin superfamily.

### 3. *Lampetra fluviatilis* (European River Lamprey, GCF_964198595.1)
- **Assembly Quality**: Chromosome-level (`kcLamFluv1.1`, Darwin Tree of Life), 1.04 Gb, Scaffold N50 = 13.2 Mb. RefSeq Annotation Release: January 2026 (17,054 protein-coding genes).
- **Locus Synteny**:
  - **Proto-Globin Locus (`NC_136297.1`: 11,400,000–12,300,000, ~900 kb)**:
    - Contains `NPRL3` (11.41 Mb) adjacent to a tandem array of 17 globin genes (`LOC144940938`, `LOC144940937`, `LOC144940936`, etc.).
- **Macroevolutionary Value**: **High**. Our current cyclostome representative, *Petromyzon marinus*, is an older, fragmented assembly where globin loci reside on an unplaced scaffold (`GL478636`). *Lampetra fluviatilis* provides a complete, modern chromosome-anchored lamprey assembly that enables precise telomeric distance and isochore embedding calculations.

### 4. *Pelodiscus sinensis* (Chinese Softshell Turtle, GCF_049634645.1)
- **Assembly Quality**: Chromosome-level (`ASM4963464v1`), 2.24 Gb, Scaffold N50 = 129.6 Mb. RefSeq Annotation Release: August 2025 (21,463 protein-coding genes).
- **Locus Synteny**:
  - **$\alpha$-Globin Locus (`NC_134726.1`: 4,640,000–4,940,000, ~300 kb)**:
    - Canonical tetrapod cluster: `RHBDF1` (4.83 Mb) — `MPG` (4.79 Mb) — `NPRL3` (4.68 Mb) — `pi` — `alpha_D` — `alpha_A`.
  - **$\beta$-Globin Locus (`NC_134711.1`: 331,060,000–331,100,000, ~40 kb)**:
    - Contains 3 tandem $\beta$-globins (`LOC142826610`, `LOC142818157`, `LOC102453687`).
- **Macroevolutionary Value**: **High**. Softshell turtles (Trionychidae) represent the earliest diverging extant branch of Cryptodira (~150 Mya divergence). While *Chelonia*, *Eretmochelys*, and *Chrysemys* represent hard-shelled turtles, *Pelodiscus* adds morphological and phylogenetic divergence to Testudines.

---

## 3. Secondary & Redundant Candidates

- **`GCF_059951505.1` (*Python bivittatus*, Burmese Python)**:
  - Intact $\alpha$ and $\beta$ clusters on scaffolds `NW_028277083.1` (240 kb) and `NW_028277056.1` (710 kb). However, because these scaffolds are not anchored to numbered chromosome arms, physical distance to the telomere cannot be directly measured.
- **`GCF_031769625.1` (*Taudactylus pleione*, Tinker Frog)**:
  - 5.51 Gb genome size with scaffold-level assembly (`aTauple2`). While it has 21 linked globins on `NW_028284468.1`, its large size and contig status provide lower technical contiguity than our existing amphibians (*Xenopus*, *Bufo*, caecilians, and urodeles).
- **`GCF_965140235.1` (*Caretta caretta*, Loggerhead Sea Turtle)**:
  - Outstanding technical quality, but functionally identical in locus synteny, telomeric distance, and sequence to *Chelonia mydas* and *Eretmochelys imbricata* (3rd sea turtle in the same family, Cheloniidae).

---

## 4. Renaming of *Chrysemys picta bellii*

Per your instruction, *Chrysemys picta bellii* has been renamed across the codebase:
- **Previous Label**: `"Turtle (C. picta)"` / `"Turtle α"` / `"Turtle β"`
- **New Descriptive Label**: **`"Painted Turtle (C. picta)"`** / **`"Painted Turtle α"`** / **`"Painted Turtle β"`**
- **Clarity Gain**: Avoids ambiguity against **"Green Sea Turtle (*C. mydas*)"**, **"Hawksbill Turtle (*E. imbricata*)"**, and **"Softshell Turtle (*P. sinensis*)"**.
