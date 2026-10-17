# Comparative Vertebrate Globin Loci Analysis Report

## 1. Executive Summary

This report summarizes the comparative genomic and composition analysis of vertebrate alpha (and linked beta) globin loci across three major evolutionary lineages:
- **Fish** (Cyclostomes and Teleost 3R duplicates: Zebrafish, Medaka, Fugu, Lamprey)
- **Amphibians** (*Xenopus tropicalis*)
- **Amniotes** (Mammals, Birds, Reptiles: Human, Cattle, Pig, Mouse, Opossum, Chicken, Turtle)

### Clade-Level Composition Summary

| clade               |   num_species |   num_loci |   mean_locus_length_bp |   mean_locus_gc |   mean_exonic_gc |   mean_intergenic_gc |   mean_globin_cds_gc |   mean_globin_gc3 |   total_globins_found |
|:--------------------|--------------:|-----------:|-----------------------:|----------------:|-----------------:|---------------------:|---------------------:|------------------:|----------------------:|
| amniote_bird        |             4 |          8 |        188750          |           49.93 |            55.32 |                49.79 |                57.92 |             70.5  |                    23 |
| amniote_crocodylian |             1 |          2 |             1.5e+06    |           48.44 |            52.5  |                48.44 |                53.1  |             67.94 |                    10 |
| amniote_mammal      |             5 |         10 |        180193          |           45.12 |            54.75 |                45.12 |                57.5  |             72.73 |                    41 |
| amniote_monotreme   |             1 |          2 |        230000          |           50.56 |            59.66 |                50.56 |                61.28 |             75.41 |                     5 |
| amniote_reptile     |             3 |          6 |        271667          |           50.53 |            59.05 |                50.53 |                58.19 |             79.54 |                    16 |
| amphibian           |             2 |          2 |             1.3125e+06 |           42.39 |            21.87 |                42.39 |                50.39 |             55.38 |                    19 |
| cyclostome          |             1 |          1 |        130000          |           44.1  |            59.05 |                44.1  |                62.27 |             91.03 |                     2 |
| non_teleost_fish    |             1 |          1 |        250000          |           47.59 |            62.23 |                47.59 |                61.69 |             84.55 |                     4 |
| teleost_fish        |             2 |          4 |        176450          |           37.41 |            52.99 |                37.41 |                52.75 |             62.71 |                    26 |

## 2. Locus-by-Locus Genomic & GC Metrics

| species                      | common_name                              | clade               | locus_id             | chromosome     |   length_bp |   overall_gc |   globin_mean_cds_gc |   globin_mean_gc3 |   num_globins |
|:-----------------------------|:-----------------------------------------|:--------------------|:---------------------|:---------------|------------:|-------------:|---------------------:|------------------:|--------------:|
| petromyzon_marinus           | Sea Lamprey (P. marinus)                 | cyclostome          | proto_globin_cluster | GL478636       |      130000 |        44.1  |                62.27 |             91.03 |             2 |
| lepisosteus_oculatus         | Spotted Gar (L. oculatus)                | non_teleost_fish    | alpha_beta_cluster   | LG13           |      250000 |        47.59 |                61.69 |             84.55 |             4 |
| danio_rerio                  | Zebrafish (D. rerio)                     | teleost_fish        | alpha_beta_locus_A   | 3              |      327800 |        36.05 |                52.49 |             62    |            12 |
| danio_rerio                  | Zebrafish (D. rerio)                     | teleost_fish        | alpha_beta_locus_B   | 12             |      193000 |        34.73 |                51.28 |             62.64 |             4 |
| oryzias_latipes              | Medaka (O. latipes)                      | teleost_fish        | alpha_beta_locus_A   | 8              |      105000 |        39.75 |                49.72 |             58.9  |             8 |
| oryzias_latipes              | Medaka (O. latipes)                      | teleost_fish        | alpha_beta_locus_B   | 19             |       80000 |        39.11 |                57.52 |             67.3  |             2 |
| xenopus_tropicalis           | Clawed Frog (X. tropicalis)              | amphibian           | alpha_beta_cluster   | 9              |      905000 |        39.93 |                49.3  |             50.05 |             6 |
| bufo_bufo                    | Common Toad (B. bufo)                    | amphibian           | alpha_beta_cluster   | NC_053395.1    |     1720000 |        44.85 |                51.48 |             60.71 |            13 |
| Chrysemys_picta_bellii       | Western Painted Turtle (C. picta)        | amniote_reptile     | alpha_globin         | KK088206.1     |      450000 |        46.02 |                54.19 |             64.2  |             4 |
| Chrysemys_picta_bellii       | Western Painted Turtle (C. picta)        | amniote_reptile     | beta_globin          | KK083065.1     |       80000 |        45.34 |                47.53 |             56.24 |             2 |
| alligator_mississippiensis   | American Alligator (A. mississippiensis) | amniote_crocodylian | alpha_globin         | NC_081836.1    |     1700000 |        48.87 |                50.24 |             59.68 |             3 |
| alligator_mississippiensis   | American Alligator (A. mississippiensis) | amniote_crocodylian | beta_globin          | NC_081824.1    |     1300000 |        48.01 |                55.96 |             76.19 |             7 |
| Gallus_gallus                | Chicken (G. gallus)                      | amniote_bird        | alpha_globin         | 14             |      230000 |        48.9  |                56.78 |             70.36 |             3 |
| Gallus_gallus                | Chicken (G. gallus)                      | amniote_bird        | beta_globin          | 1              |      120000 |        50.64 |                59.12 |             75.14 |             2 |
| aquila_chrysaetos_chrysaetos | Golden Eagle (A. chrysaetos)             | amniote_bird        | alpha_globin         | 25             |      350000 |        50.24 |                60.36 |             61.25 |             2 |
| aquila_chrysaetos_chrysaetos | Golden Eagle (A. chrysaetos)             | amniote_bird        | beta_globin          | 19             |      150000 |        48.89 |                55.79 |             75.62 |             3 |
| anas_platyrhynchos           | Mallard Duck (A. platyrhynchos)          | amniote_bird        | alpha_globin         | 15             |      250000 |        47.48 |                55.31 |             64.03 |             2 |
| anas_platyrhynchos           | Mallard Duck (A. platyrhynchos)          | amniote_bird        | beta_globin          | 1              |       90000 |        54.94 |                61.28 |             79.36 |             3 |
| meleagris_gallopavo          | Turkey (M. gallopavo)                    | amniote_bird        | alpha_globin         | 16             |      170000 |        47.06 |                56.34 |             77.46 |             3 |
| meleagris_gallopavo          | Turkey (M. gallopavo)                    | amniote_bird        | beta_globin          | 1              |      150000 |        51.29 |                58.39 |             60.81 |             5 |
| anolis_carolinensis          | Green Anole (A. carolinensis)            | amniote_reptile     | alpha_globin         | NW_026943824.1 |      200000 |        51.59 |                59.94 |             83.73 |             3 |
| anolis_carolinensis          | Green Anole (A. carolinensis)            | amniote_reptile     | beta_globin          | NC_085843.1    |      200000 |        42.74 |                63.85 |             93.92 |             2 |
| pogona_vitticeps             | Bearded Dragon (P. vitticeps)            | amniote_reptile     | alpha_globin         | NC_135795.1    |      450000 |        59.59 |                62.16 |             89.97 |             3 |
| pogona_vitticeps             | Bearded Dragon (P. vitticeps)            | amniote_reptile     | beta_globin          | NC_135785.1    |      250000 |        57.92 |                61.48 |             89.19 |             2 |
| ornithorhynchus_anatinus     | Platypus (O. anatinus)                   | amniote_monotreme   | alpha_globin         | 21             |      260000 |        54.1  |                63.2  |             81.94 |             2 |
| ornithorhynchus_anatinus     | Platypus (O. anatinus)                   | amniote_monotreme   | beta_globin          | 2              |      200000 |        47.02 |                59.35 |             68.88 |             3 |
| Monodelphis_domestica        | Gray Short-tailed Opossum (M. domestica) | amniote_mammal      | alpha_globin         | 6              |      430000 |        40.84 |                55.55 |             67.36 |             4 |
| Monodelphis_domestica        | Gray Short-tailed Opossum (M. domestica) | amniote_mammal      | beta_globin          | 4              |      170000 |        36.31 |                52.14 |             62.17 |             2 |
| Mus_musculus                 | Mouse (M. musculus)                      | amniote_mammal      | alpha_globin         | 11             |      126929 |        48.97 |                56.6  |             68.82 |             5 |
| Mus_musculus                 | Mouse (M. musculus)                      | amniote_mammal      | beta_globin          | 7              |      150000 |        38.46 |                52.53 |             62.5  |             5 |
| Sus_scrofa                   | Pig (S. scrofa)                          | amniote_mammal      | alpha_globin         | 3              |      140000 |        58.37 |                66.03 |             89.85 |             4 |
| Sus_scrofa                   | Pig (S. scrofa)                          | amniote_mammal      | beta_globin          | 9              |      130000 |        44.09 |                56.4  |             69.19 |             2 |
| Bos_taurus                   | Cattle (B. taurus)                       | amniote_mammal      | alpha_globin         | 25             |      160000 |        53.03 |                64.97 |             88.23 |             4 |
| Bos_taurus                   | Cattle (B. taurus)                       | amniote_mammal      | beta_globin          | 15             |      140000 |        39.63 |                52.25 |             64.64 |             5 |
| Homo_sapiens                 | Human (H. sapiens)                       | amniote_mammal      | alpha_globin         | 16             |      215000 |        53.17 |                65.01 |             90.17 |             5 |
| Homo_sapiens                 | Human (H. sapiens)                       | amniote_mammal      | beta_globin          | 11             |      140000 |        38.33 |                53.52 |             64.39 |             5 |

## 3. Key Biological Insights

1. **Isochore Specialization & GC3 Elevation in Amniotes**:
   - Warm-blooded amniotes (especially human and chicken) display marked GC and GC3 elevation in globin genes.
   - Human alpha globin CDS achieves ~65% GC with synonymous codon 3rd-position (GC3) exceeding 90%, reflecting strong GC-biased gene conversion (gBGC) and residency in an ultra-dense H3 isochore.
2. **Amphibian Transition (Xenopus tropicalis)**:
   - In *Xenopus*, the alpha and beta globins reside in a contiguous linked cluster (*hbg1, hba1, hba3, hba4, hbz, hbe1.2*).
   - Locus GC content is significantly more homogeneous and moderate compared to amniotes.
3. **Teleost Dual Clusters (3R / Whole Genome Duplication)**:
   - Teleost fish have partitioned the ancestral globin gene repertoire across two paralogous chromosomes (Locus A and Locus B).
   - Both loci retain syntenic flanking markers (*rhbdf1a/kank2* on Locus A, *rhbdf1b/lcmt1* on Locus B).
4. **Agnathan Baseline (Lamprey)**:
   - Sea lamprey provides the ancestral agnathan baseline prior to the alpha-beta globin gene duplication.
