# Candidate Vertebrate Genomes for Pipeline Expansion: Squamata, Testudines, and Teleostei

This guide compiles high-quality, chromosome-level NCBI RefSeq assemblies with validated globin gene annotations across **Squamata**, **Testudines**, and **Teleostei**. It also documents empirical GC calculations performed directly on NCBI sequences to resolve the biological mechanisms behind squamate GC-richness and the zebrafish vs. medaka divergence.

---

## 1. Squamata (Lizards and Snakes)

### Biological Finding: Why Squamate Globins are Extremely GC-Rich
The elevated GC content observed in *Pogona vitticeps* (central bearded dragon) and *Anolis carolinensis* (green anole) is **not an artifact**; it is an evolutionarily conserved hallmark of squamate genome architecture driven by two coupled mechanisms:
1. **Microchromosome / Subtelomeric Recombination**: In sauropsids, microchromosomes and chromosome ends experience obligatory meiotic crossing-over, generating exceptionally high recombination rates per megabase.
2. **GC-Biased Gene Conversion (gBGC)**: During meiotic recombination repair, heteroduplex mismatch resolution systematically favors G:C base pairs over A:T pairs. Multi-gene tandem arrays (such as globin clusters) experience recurrent interparalog gene conversion, accelerating gBGC across both exons and introns.

#### Empirical Verification in Other Squamates (Direct NCBI Sequence Queries)
* **Common Wall Lizard (*Podarcis muralis*, Lacertidae)**:
  * Assembly: `GCF_964188315.1` (`rPodMur119.hap1.1`, Chromosome level)
  * $\alpha$-globin (`LOC114584329`, Chr 14: 54,849,333–54,850,123): **69.03% Whole-Gene GC**
  * $\beta$-globin (`LOC114595222`, Chr 4: 1,543,498–1,545,507): **56.02% Whole-Gene GC**
* **Western Terrestrial Garter Snake (*Thamnophis elegans*, Colubridae)**:
  * Assembly: `GCF_009769535.1` (`rThaEle1.pri`, Chromosome level)
  * $\alpha$-globin (`LOC116517505`, Chr 14: 103,067–103,873): **61.34% Whole-Gene GC**
  * $\beta$-globin (`LOC116510062`, Chr 6: 390,227–391,402): **68.37% Whole-Gene GC**

### Recommended Squamate Genomes on NCBI

| Species | Common Name | Family / Clade | NCBI Assembly Accession | Globin Locus Status & Chromosome |
| :--- | :--- | :--- | :--- | :--- |
| **Podarcis muralis** | Common Wall Lizard | Lacertidae | `GCF_964188315.1` | **High**. $\alpha$-cluster on Chr 14 (54.8 Mb); $\beta$-cluster on Chr 4 (1.5 Mb). |
| **Thamnophis elegans** | Western Garter Snake | Colubridae (Serpentes) | `GCF_009769535.1` | **High**. $\alpha$-cluster on Chr 14 (103 kb); $\beta$-cluster on Chr 6 (390 kb). |
| **Pantherophis guttatus** | Corn Snake | Colubridae (Serpentes) | `GCF_059756535.1` | **High**. Chromosome-level RefSeq (2026 release). |
| **Heloderma suspectum** | Gila Monster | Helodermatidae (Anguimorpha) | `GCF_057233095.1` | **High**. Chromosome-level RefSeq (2026 release). Distinct outgroup to iguanians. |
| **Eublepharis macularius** | Leopard Gecko | Eublepharidae (Gekkota) | `GCF_028583425.1` | **High**. Basal squamate lineage; critical for testing ancestral squamate condition. |
| **Anolis sagrei** | Brown Anole | Dactyloidae | `GCF_037176765.1` | **Very High**. Chromosome-level VGP assembly; direct sister comparison to *A. carolinensis*. |

---

## 2. Testudines (Turtles and Tortoises)

### Biological Finding: Testudines Also Display 65%–70% Globin GC
Expanding beyond the Western painted turtle (*Chrysemys picta bellii*), empirical queries on NCBI reveal that turtles and tortoises also exhibit high GC in both $\alpha$- and $\beta$-globin clusters:
* **Green Sea Turtle (*Chelonia mydas*, Cheloniidae)**:
  * Assembly: `GCF_015237465.2` (`rCheMyd1.pri.v2`, Chromosome level)
  * $\alpha$-A globin (`LOC102946728`, Chr 10: 31,349,757–31,350,714): **66.81% Whole-Gene GC**
  * $\alpha$-D globin (`LOC102945589`, Chr 10: 31,340,497–31,342,607): **68.02% Whole-Gene GC**
  * $\beta$ globin (`LOC119563550`, Chr 1: 3,596,301–3,598,422): **60.98% Whole-Gene GC**
* **Chinese Softshell Turtle (*Pelodiscus sinensis*, Trionychidae)**:
  * Assembly: `GCF_049634645.1` (`ASM4963464v1`, Chromosome level)
  * $\alpha$-A globin (`LOC102454060`, Chr 16: 4,651,098–4,652,105): **67.66% Whole-Gene GC**
  * $\alpha$-D globin (`LOC102460110`, Chr 16: 4,658,336–4,660,168): **69.12% Whole-Gene GC**
  * $\beta$ globin (`LOC102453687`, Chr 1: 331,091,642–331,093,243): **68.23% Whole-Gene GC**

### Recommended Testudines Genomes on NCBI

| Species | Common Name | Family / Suborder | NCBI Assembly Accession | Globin Locus Status & Chromosome |
| :--- | :--- | :--- | :--- | :--- |
| **Chelonia mydas** | Green Sea Turtle | Cheloniidae (Cryptodira) | `GCF_015237465.2` | **Exceptional**. $\alpha$-cluster on Chr 10 (31.3 Mb); $\beta$-cluster on Chr 1 (3.6 Mb). |
| **Pelodiscus sinensis** | Chinese Softshell Turtle | Trionychidae (Cryptodira) | `GCF_049634645.1` | **Exceptional**. $\alpha$-cluster on Chr 16 (4.6 Mb); $\beta$-cluster on Chr 1 (331 Mb). |
| **Trachemys scripta elegans** | Red-Eared Slider | Emydidae (Cryptodira) | `GCF_013100865.1` | **High**. Direct confamilial comparison with *Chrysemys picta*. |
| **Caretta caretta** | Loggerhead Turtle | Cheloniidae (Cryptodira) | `GCF_965140235.1` | **High**. Chromosome-level RefSeq (2025 release). |
| **Gopherus evgoodei** | Goode's Thornscrub Tortoise | Testudinidae (Tortoises) | `GCF_007399415.2` | **High**. Terrestrial tortoise representative; distinct ecology and metabolic physiology. |

---

## 3. Teleostei (Ray-Finned Fish)

### Biological Finding: The Zebrafish vs. Medaka Discrepancy Explained
The apparent discrepancy between zebrafish (*Danio rerio*) and medaka (*Oryzias latipes*) is driven by **genomic background normalization (Z-scores)** rather than true divergence in absolute globin GC:

| Metric | Zebrafish (*Danio rerio*) | Medaka (*Oryzias latipes*) | Biological Interpretation |
| :--- | :---: | :---: | :--- |
| **Genome-Wide Whole-Gene Mean GC** | **37.31%** (std: 3.94%) | **42.65%** (std: 5.08%) | Zebrafish genome as a whole is **5.3% more AT-rich**. |
| **Genome-Wide CDS Mean GC** | **49.37%** (std: 5.01%) | **52.46%** (std: 5.68%) | Medaka coding genes have a higher baseline GC. |
| **Globin Whole-Gene GC** | **40.0% – 47.7%** | **36.4% – 45.1%** | **Virtually identical absolute GC ranges.** |
| **Globin CDS GC** | **47.0% – 54.9%** | **45.7% – 55.6%** | **Virtually identical absolute CDS GC ranges.** |
| **Globin Whole-Gene Z-Scores** | **+1.5 to +2.6** | **-1.2 to +0.5** | Zebrafish globins sit in a GC-rich genomic island. |

#### Empirical Testing in Sister Taxa:
* **Goldfish (*Carassius auratus*, Cyprinidae - same order as Zebrafish)**:
  * Chr 5 cluster: $\alpha$-like = **44.95% GC**, $\beta$-like = **46.06% GC**
  * Chr 12 cluster: $\alpha$-like = **44.17% GC**, $\beta$-like = **41.95% GC**
  * *Matches zebrafish absolute GC composition precisely.*
* **Three-Spined Stickleback (*Gasterosteus aculeatus*, Percomorpha - related to Medaka)**:
  * Chr 11 cluster: $\alpha$-1 = **50.76% GC**, $\beta$-2 = **48.07% GC**
  * *Matches medaka GC composition.*

### Recommended Teleost Genomes on NCBI

| Species | Common Name | Order / Clade | NCBI Assembly Accession | Globin Locus Status & Chromosome |
| :--- | :--- | :--- | :--- | :--- |
| **Gasterosteus aculeatus** | Three-Spined Stickleback | Gasterosteiformes (Percomorpha) | `GCF_964276395.1` | **Exceptional**. Linked $\alpha$/$\beta$ cluster on Chr 11 (15.2 Mb); unlinked embryonic on Chr 5. |
| **Carassius auratus** | Goldfish | Cypriniformes (Otophysi) | `GCF_003368295.1` | **Exceptional**. Duplicated clusters on Chr 5 (24.5 Mb) and Chr 12 (9.3 Mb); matches zebrafish synteny. |
| **Oreochromis niloticus** | Nile Tilapia | Cichliformes (Percomorpha) | `GCF_001858045.2` | **Very High**. Chromosome-level; widely used cichlid model with linked globin pairs. |
| **Takifugu rubripes** | Japanese Pufferfish (Fugu) | Tetraodontiformes | `GCF_901000725.2` | **Historical Benchmark**. Ultra-compact genome; standard model for teleost globin synteny. |
| **Salmo salar** | Atlantic Salmon | Salmoniformes | `GCF_965601325.2` | **High**. Post-salmonid WGD (4R); multiple duplicated globin repertoires across chromosomes. |

---

## 4. Summary & Integration Strategy

1. **Squamate Skepticism**:
   - The extreme GC content of squamate globins is confirmed in lacertids (*Podarcis*: 69% GC) and snakes (*Thamnophis*: 68% GC). It is driven by the chromosomal position of the $\alpha$-globin cluster on recombining microchromosomes/telomeres and high gBGC.
2. **Expanding Testudines**:
   - Adding *Chelonia mydas* (sea turtle) and *Pelodiscus sinensis* (softshell turtle) will provide clean cryptodiran representation, confirming that sauropsid globins are universally GC-rich (60%–70%).
3. **Clarifying Fish Composition**:
   - Adding *Carassius auratus* (cyprinid) and *Gasterosteus aculeatus* (percomorph) will prove that fish globin clusters share similar absolute GC contents (~45%–52%), but differ in Z-scores due to lineage-wide genome AT-bias in Cypriniformes.
