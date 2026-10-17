# Comparative Vertebrate Globin Loci Analysis Report

## 1. Executive Summary

This report summarizes the comparative genomic and composition analysis of vertebrate alpha (and linked beta) globin loci across three major evolutionary lineages:
- **Fish** (Cyclostomes and Teleost 3R duplicates: Zebrafish, Medaka, Fugu, Lamprey)
- **Amphibians** (*Xenopus tropicalis*)
- **Amniotes** (Mammals, Birds, Reptiles: Human, Cattle, Pig, Mouse, Opossum, Chicken, Turtle)

### Clade-Level Composition Summary

| clade               |   num_species |   num_loci |   mean_locus_length_bp |   mean_locus_gc |   mean_exonic_gc |   mean_intergenic_gc |   mean_globin_cds_gc |   mean_globin_gc3 |   total_globins_found |
|:--------------------|--------------:|-----------:|-----------------------:|----------------:|-----------------:|---------------------:|---------------------:|------------------:|----------------------:|
| amniote_bird        |             8 |         16 |       171188           |           50.38 |            54.67 |                49.9  |                58.34 |             60.36 |                    50 |
| amniote_crocodylian |             1 |          2 |            1.5e+06     |           48.44 |            50.39 |                48.28 |                53.1  |             67.94 |                    10 |
| amniote_mammal      |            14 |         28 |       147570           |           45.08 |            53.68 |                44.36 |                57.93 |             61.37 |                   126 |
| amniote_monotreme   |             1 |          2 |       230000           |           50.56 |            59.6  |                50.19 |                61.31 |             78.06 |                     5 |
| amniote_reptile     |             9 |         18 |       322501           |           51.02 |            55.88 |                50.34 |                60.22 |             58.24 |                    51 |
| amphibian           |             7 |         11 |            1.82455e+06 |           44.03 |            47.97 |                43.94 |                53.08 |             49.69 |                    91 |
| chondrichthyes      |             1 |          1 |       280001           |           49.67 |            42.47 |                49.75 |                47.53 |             46.4  |                     4 |
| cyclostome          |             3 |          3 |       420001           |           47.52 |            55.23 |                47.01 |                57.78 |             56.75 |                    15 |
| non_teleost_fish    |             1 |          1 |       250000           |           47.59 |            53.62 |                46.64 |                61.69 |             84.55 |                     4 |
| teleost_fish        |             5 |          9 |       114534           |           39.88 |            47.26 |                38.33 |                52.9  |             59.14 |                    59 |

## 2. Locus-by-Locus Genomic & GC Metrics

| species                      | common_name                              | clade               | locus_id             | chromosome     |   length_bp |   overall_gc |   globin_mean_cds_gc |   globin_mean_gc3 |   num_globins |
|:-----------------------------|:-----------------------------------------|:--------------------|:---------------------|:---------------|------------:|-------------:|---------------------:|------------------:|--------------:|
| petromyzon_marinus           | Sea Lamprey (P. marinus)                 | cyclostome          | proto_globin_cluster | GL478636       |      130000 |        44.1  |                62.27 |             91.03 |             2 |
| lepisosteus_oculatus         | Spotted Gar (L. oculatus)                | non_teleost_fish    | alpha_beta_cluster   | LG13           |      250000 |        47.59 |                61.69 |             84.55 |             4 |
| danio_rerio                  | Zebrafish (D. rerio)                     | teleost_fish        | alpha_beta_locus_A   | 3              |      327800 |        36.05 |                52.44 |             62.56 |            12 |
| danio_rerio                  | Zebrafish (D. rerio)                     | teleost_fish        | alpha_beta_locus_B   | 12             |      193000 |        34.73 |                51.28 |             62.64 |             4 |
| oryzias_latipes              | Medaka (O. latipes)                      | teleost_fish        | alpha_beta_locus_A   | 8              |      105000 |        39.75 |                50.06 |             58.93 |             8 |
| oryzias_latipes              | Medaka (O. latipes)                      | teleost_fish        | alpha_beta_locus_B   | 19             |       80000 |        39.11 |                55.14 |             72.2  |             2 |
| takifugu_rubripes            | Fugu (T. rubripes)                       | teleost_fish        | alpha_beta_locus_A   | 5              |       65000 |        46.94 |                55.1  |             72.79 |             1 |
| takifugu_rubripes            | Fugu (T. rubripes)                       | teleost_fish        | alpha_beta_locus_B   | 1              |       45000 |        43.71 |                52.88 |             64.81 |             4 |
| gasterosteus_aculeatus       | Stickleback (G. aculeatus)               | teleost_fish        | alpha_beta_cluster   | NC_135699.1    |       90001 |        45.12 |                56.19 |             48.22 |            11 |
| xenopus_tropicalis           | Clawed Frog (X. tropicalis)              | amphibian           | alpha_beta_cluster   | 9              |      950000 |        39.83 |                48.4  |             49.88 |             6 |
| bufo_bufo                    | Common Toad (B. bufo)                    | amphibian           | alpha_beta_cluster   | NC_053395.1    |     1720000 |        44.85 |                51.48 |             60.71 |            13 |
| Chrysemys_picta_bellii       | Western Painted Turtle (C. picta)        | amniote_reptile     | alpha_globin         | KK088206.1     |      450000 |        46.02 |                55.89 |             73.38 |             4 |
| Chrysemys_picta_bellii       | Western Painted Turtle (C. picta)        | amniote_reptile     | beta_globin          | KK083065.1     |       80000 |        45.34 |                46.11 |             55.11 |             2 |
| chelonia_mydas               | Green Sea Turtle (C. mydas)              | amniote_reptile     | alpha_globin         | NC_051250.2    |      180001 |        46.39 |                59.85 |             41.04 |             3 |
| chelonia_mydas               | Green Sea Turtle (C. mydas)              | amniote_reptile     | beta_globin          | NC_057849.1    |      150001 |        51.88 |                61.04 |             39.19 |             1 |
| pelodiscus_sinensis          | Softshell Turtle (P. sinensis)           | amniote_reptile     | alpha_globin         | NC_134726.1    |      300001 |        46.98 |                58.73 |             60.05 |             3 |
| pelodiscus_sinensis          | Softshell Turtle (P. sinensis)           | amniote_reptile     | beta_globin          | NC_134711.1    |       40001 |        55.77 |                64.26 |             61.04 |             3 |
| alligator_mississippiensis   | American Alligator (A. mississippiensis) | amniote_crocodylian | alpha_globin         | NC_081836.1    |     1700000 |        48.87 |                50.24 |             59.68 |             3 |
| alligator_mississippiensis   | American Alligator (A. mississippiensis) | amniote_crocodylian | beta_globin          | NC_081824.1    |     1300000 |        48.01 |                55.96 |             76.19 |             7 |
| Gallus_gallus                | Chicken (G. gallus)                      | amniote_bird        | alpha_globin         | 14             |      230000 |        48.9  |                56.87 |             70.31 |             3 |
| Gallus_gallus                | Chicken (G. gallus)                      | amniote_bird        | beta_globin          | 1              |      120000 |        50.64 |                59.28 |             77.8  |             2 |
| aquila_chrysaetos_chrysaetos | Golden Eagle (A. chrysaetos)             | amniote_bird        | alpha_globin         | 25             |      350000 |        50.24 |                60.36 |             61.25 |             2 |
| aquila_chrysaetos_chrysaetos | Golden Eagle (A. chrysaetos)             | amniote_bird        | beta_globin          | 19             |      150000 |        48.89 |                56.22 |             76.05 |             3 |
| anas_platyrhynchos           | Mallard Duck (A. platyrhynchos)          | amniote_bird        | alpha_globin         | 15             |      250000 |        47.48 |                55.45 |             68.81 |             2 |
| anas_platyrhynchos           | Mallard Duck (A. platyrhynchos)          | amniote_bird        | beta_globin          | 1              |       90000 |        54.94 |                59.49 |             77.73 |             3 |
| meleagris_gallopavo          | Turkey (M. gallopavo)                    | amniote_bird        | alpha_globin         | 16             |      170000 |        47.06 |                55.69 |             72.21 |             3 |
| meleagris_gallopavo          | Turkey (M. gallopavo)                    | amniote_bird        | beta_globin          | 1              |      150000 |        51.29 |                58.66 |             62.85 |             5 |
| struthio_camelus_australis   | African Ostrich (S. camelus)             | amniote_bird        | alpha_globin         | NW_009271137.1 |      190001 |        46.99 |                57.32 |             57.64 |             3 |
| struthio_camelus_australis   | African Ostrich (S. camelus)             | amniote_bird        | beta_globin          | NW_009270861.1 |       29001 |        49.42 |                59.64 |             51.34 |             4 |
| anolis_carolinensis          | Green Anole (A. carolinensis)            | amniote_reptile     | alpha_globin         | NW_026943824.1 |      200000 |        51.59 |                60.29 |             84.04 |             3 |
| anolis_carolinensis          | Green Anole (A. carolinensis)            | amniote_reptile     | beta_globin          | NC_085843.1    |      200000 |        42.74 |                63.85 |             93.92 |             2 |
| pogona_vitticeps             | Bearded Dragon (P. vitticeps)            | amniote_reptile     | alpha_globin         | NC_135795.1    |      450000 |        59.59 |                62.16 |             89.97 |             3 |
| pogona_vitticeps             | Bearded Dragon (P. vitticeps)            | amniote_reptile     | beta_globin          | NC_135785.1    |      250000 |        57.92 |                61.48 |             89.19 |             2 |
| podarcis_muralis             | Common Wall Lizard (P. muralis)          | amniote_reptile     | alpha_globin         | NC_135668.1    |      220001 |        53.84 |                67.37 |             53.84 |             2 |
| podarcis_muralis             | Common Wall Lizard (P. muralis)          | amniote_reptile     | beta_globin          | NC_135658.1    |      235001 |        49.1  |                59.78 |             38.72 |             3 |
| thamnophis_elegans           | Western Garter Snake (T. elegans)        | amniote_reptile     | alpha_globin         | NC_045554.1    |       85001 |        55.52 |                57.45 |             52.48 |             2 |
| thamnophis_elegans           | Western Garter Snake (T. elegans)        | amniote_reptile     | beta_globin          | NC_045546.1    |      345001 |        53.53 |                63.62 |             36.83 |             2 |
| ornithorhynchus_anatinus     | Platypus (O. anatinus)                   | amniote_monotreme   | alpha_globin         | 21             |      260000 |        54.1  |                63.16 |             83.84 |             2 |
| ornithorhynchus_anatinus     | Platypus (O. anatinus)                   | amniote_monotreme   | beta_globin          | 2              |      200000 |        47.02 |                59.46 |             72.27 |             3 |
| Monodelphis_domestica        | Gray Short-tailed Opossum (M. domestica) | amniote_mammal      | alpha_globin         | 6              |      430000 |        40.84 |                55.44 |             67.36 |             4 |
| Monodelphis_domestica        | Gray Short-tailed Opossum (M. domestica) | amniote_mammal      | beta_globin          | 4              |      170000 |        36.31 |                52.14 |             62.17 |             2 |
| Mus_musculus                 | Mouse (M. musculus)                      | amniote_mammal      | alpha_globin         | 11             |      126929 |        48.97 |                56.6  |             68.82 |             5 |
| Mus_musculus                 | Mouse (M. musculus)                      | amniote_mammal      | beta_globin          | 7              |      150000 |        38.46 |                51.94 |             62.3  |             5 |
| Sus_scrofa                   | Pig (S. scrofa)                          | amniote_mammal      | alpha_globin         | 3              |      140000 |        58.37 |                66.03 |             89.85 |             4 |
| Sus_scrofa                   | Pig (S. scrofa)                          | amniote_mammal      | beta_globin          | 9              |      130000 |        44.09 |                56.4  |             69.19 |             2 |
| Bos_taurus                   | Cattle (B. taurus)                       | amniote_mammal      | alpha_globin         | 25             |      160000 |        53.03 |                64.3  |             89.17 |             4 |
| Bos_taurus                   | Cattle (B. taurus)                       | amniote_mammal      | beta_globin          | 15             |      140000 |        39.63 |                52.36 |             64.86 |             5 |
| Homo_sapiens                 | Human (H. sapiens)                       | amniote_mammal      | alpha_globin         | 16             |      215000 |        53.17 |                65.87 |             90.48 |             5 |
| Homo_sapiens                 | Human (H. sapiens)                       | amniote_mammal      | beta_globin          | 11             |      140000 |        38.33 |                54.19 |             65.41 |             5 |
| equus_asinus                 | Donkey (E. asinus)                       | amniote_mammal      | alpha_globin         | NC_091803.1    |      160001 |        52.26 |                65.27 |             63.29 |             4 |
| equus_asinus                 | Donkey (E. asinus)                       | amniote_mammal      | beta_globin          | NC_091809.1    |      100001 |        38.43 |                53.4  |             60.38 |             4 |
| equus_caballus               | Horse (E. caballus)                      | amniote_mammal      | alpha_globin         | NC_091696.1    |      170001 |        51.8  |                64.26 |             62.01 |             3 |
| equus_caballus               | Horse (E. caballus)                      | amniote_mammal      | beta_globin          | NC_091690.1    |       80001 |        38.79 |                53.29 |             60.55 |             4 |
| miniopterus_natalensis       | Bat (M. natalensis)                      | amniote_mammal      | alpha_globin         | NC_142083.1    |      160001 |        48.98 |                64.18 |             62.44 |             6 |
| miniopterus_natalensis       | Bat (M. natalensis)                      | amniote_mammal      | beta_globin          | NC_142091.1    |       50001 |        38.82 |                52.78 |             59.01 |             3 |
| puma_concolor                | Puma (P. concolor)                       | amniote_mammal      | alpha_globin         | NC_141682.1    |      170001 |        51.54 |                65.31 |             62.2  |             6 |
| puma_concolor                | Puma (P. concolor)                       | amniote_mammal      | beta_globin          | NC_141676.1    |      120001 |        40.17 |                54.64 |             34.74 |             4 |
| sturnus_vulgaris             | Starling (S. vulgaris)                   | amniote_bird        | alpha_globin         | NC_142258.1    |      310001 |        49.97 |                59.91 |             38.93 |             3 |
| sturnus_vulgaris             | Starling (S. vulgaris)                   | amniote_bird        | beta_globin          | NC_142244.1    |       40001 |        56.37 |                62.16 |             59.46 |             3 |
| anser_cygnoides              | Swan Goose (A. cygnoides)                | amniote_bird        | alpha_globin         | NC_089887.1    |      360001 |        49.97 |                56.39 |             37.84 |             3 |
| anser_cygnoides              | Swan Goose (A. cygnoides)                | amniote_bird        | beta_globin          | NC_089873.1    |       30001 |        51.95 |                60.21 |             58.08 |             4 |
| coturnix_japonica            | Quail (C. japonica)                      | amniote_bird        | alpha_globin         | NC_140675.1    |      240001 |        49.72 |                56.54 |             38.31 |             3 |
| coturnix_japonica            | Quail (C. japonica)                      | amniote_bird        | beta_globin          | NC_140663.1    |       30001 |        52.21 |                59.29 |             57.09 |             4 |
| rattus_norvegicus            | Norway Rat (R. norvegicus)               | amniote_mammal      | alpha_globin         | NC_141543.1    |      160001 |        47.71 |                56.5  |             59.44 |             5 |
| rattus_norvegicus            | Norway Rat (R. norvegicus)               | amniote_mammal      | beta_globin          | NC_141534.1    |      320001 |        39.14 |                52.38 |             58.79 |             9 |
| pimephales_promelas          | Minnow (P. promelas)                     | teleost_fish        | alpha_beta_major     | NC_142771.1    |      100001 |        36.83 |                52.54 |             48.61 |            13 |
| pimephales_promelas          | Minnow (P. promelas)                     | teleost_fish        | alpha_beta_embryonic | NC_142785.1    |       25001 |        36.65 |                50.5  |             41.52 |             4 |
| eretmochelys_imbricata       | Hawksbill Turtle (E. imbricata)          | amniote_reptile     | alpha_globin         | NC_135581.1    |      960001 |        47.62 |                58.65 |             36.91 |             3 |
| eretmochelys_imbricata       | Hawksbill Turtle (E. imbricata)          | amniote_reptile     | beta_globin          | NC_135572.1    |       40001 |        54.57 |                59.53 |             38.06 |             3 |
| geotrypetes_seraphini        | Gabon Caecilian (G. seraphini)           | amphibian           | alpha_beta_linked    | NC_047094.1    |      650001 |        43.1  |                49.07 |             38.46 |             3 |
| geotrypetes_seraphini        | Gabon Caecilian (G. seraphini)           | amphibian           | beta_secondary       | NC_047088.1    |       30001 |        41.8  |                47.97 |             56.08 |             1 |
| rhinatrema_bivittatum        | Two-lined Caecilian (R. bivittatum)      | amphibian           | alpha_beta_linked    | NC_042628.1    |      510001 |        43.9  |                46.84 |             40.62 |             4 |
| microcaecilia_unicolor       | Tiny Caecilian (M. unicolor)             | amphibian           | alpha_beta_linked    | NC_044038.1    |      700001 |        43.58 |                50.63 |             42.06 |             4 |
| microcaecilia_unicolor       | Tiny Caecilian (M. unicolor)             | amphibian           | beta_secondary       | NC_044037.1    |       10001 |        42.54 |                50.68 |             56.08 |             1 |
| ambystoma_mexicanum          | Axolotl (A. mexicanum)                   | amphibian           | alpha_beta_linked    | NC_090924.1    |     8200001 |        46.41 |                59.15 |             45.14 |            17 |
| ambystoma_mexicanum          | Axolotl (A. mexicanum)                   | amphibian           | beta_secondary       | NC_090925.1    |     1150001 |        47.09 |                59.99 |             52.05 |            12 |
| triturus_cristatus           | Crested Newt (T. cristatus)              | amphibian           | alpha_beta_linked    | NC_140593.1    |     5100001 |        45.25 |                59.48 |             52.5  |            22 |
| triturus_cristatus           | Crested Newt (T. cristatus)              | amphibian           | beta_secondary       | NC_140590.1    |     1050001 |        45.96 |                60.15 |             53.06 |             8 |
| zootoca_vivipara             | Viviparous Lizard (Z. vivipara)          | amniote_reptile     | alpha_globin         | NC_142596.1    |      460001 |        52.9  |                62.48 |             50.06 |             6 |
| zootoca_vivipara             | Viviparous Lizard (Z. vivipara)          | amniote_reptile     | beta_globin          | NC_142587.1    |     1160001 |        47.14 |                61.48 |             54.56 |             4 |
| eptatretus_burgeri           | Inshore Hagfish (E. burgeri)             | cyclostome          | proto_globin         | NC_140977.1    |      230001 |        44.32 |                50.76 |             37.23 |             3 |
| lampetra_fluviatilis         | River Lamprey (L. fluviatilis)           | cyclostome          | proto_globin         | NC_136297.1    |      900001 |        54.14 |                60.3  |             41.98 |            10 |
| tursiops_truncatus           | Bottlenose Dolphin (T. truncatus)        | amniote_mammal      | alpha_globin         | NC_047048.1    |      120001 |        54.86 |                63.63 |             63.18 |             5 |
| tursiops_truncatus           | Bottlenose Dolphin (T. truncatus)        | amniote_mammal      | beta_globin          | NC_047041.1    |       60001 |        38.74 |                53.8  |             60.22 |             4 |
| zalophus_californianus       | California Sea Lion (Z. californianus)   | amniote_mammal      | alpha_globin         | NC_045604.1    |      120001 |        55    |                64.53 |             42.82 |             6 |
| zalophus_californianus       | California Sea Lion (Z. californianus)   | amniote_mammal      | beta_globin          | NC_045605.1    |       50001 |        41.22 |                54.95 |             61.22 |             5 |
| loxodonta_africana           | African Savanna Elephant (L. africana)   | amniote_mammal      | alpha_globin         | NC_087353.1    |      140001 |        49.87 |                62.09 |             42.46 |             6 |
| loxodonta_africana           | African Savanna Elephant (L. africana)   | amniote_mammal      | beta_globin          | NC_087348.1    |       50001 |        39.59 |                52.83 |             35.09 |             4 |
| scyliorhinus_canicula        | Small-spotted Catshark (S. canicula)     | chondrichthyes      | alpha_beta_cluster   | NC_052160.1    |      280001 |        49.67 |                47.53 |             46.4  |             4 |
| heterocephalus_glaber        | Naked Mole-Rat (H. glaber)               | amniote_mammal      | alpha_globin         | NW_004624913.1 |      200001 |        48.01 |                59.27 |             41.26 |             4 |
| heterocephalus_glaber        | Naked Mole-Rat (H. glaber)               | amniote_mammal      | beta_globin          | NW_004624817.1 |      100001 |        36.23 |                53.7  |             59.7  |             3 |

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
