"""Gene annotation and globin classification module."""

import re
from typing import Dict, List, Any, Optional
from globin_pipeline.data_manager import DataManager


class GeneAnnotator:
    """Extracts, categorizes, and coordinates protein-coding and globin genes within loci."""

    # Specific unannotated globin genes identified in lower vertebrates / marsupials
    KNOWN_UNANNOTATED_GLOBINS = {
        # Lamprey agnathous globins
        "ENSPMAG00000001587": {"category": "other_globin", "subtype": "cyclostome_hemoglobin", "name": "pm_globin1"},
        "ENSPMAG00000001592": {"category": "other_globin", "subtype": "cyclostome_hemoglobin", "name": "pm_globin2"},
        # Opossum uncharacterized alpha globins
        "ENSMODG00000042487": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "ENSMODG00000037989": {"category": "alpha_globin", "subtype": "HBA_like", "name": "HBA_like_2"},
        "ENSMODG00000016637": {"category": "alpha_globin", "subtype": "HBA_like", "name": "HBA_like_3"},
        "ENSMODG00000016635": {"category": "alpha_globin", "subtype": "HBZ_like", "name": "HBZ_like"},
        # Turtle uncharacterized alpha globins
        "ENSCPBG00000026394": {"category": "alpha_globin", "subtype": "alpha_globin_like", "name": "HBA_like"},
        "ENSCPBG00000026395": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA"},
        # Turtle uncharacterized beta globins (KK083065.1)
        "ENSCPBG00000021825": {"category": "beta_globin", "subtype": "beta_globin_like", "name": "HBB_like_1"},
        "ENSCPBG00000021826": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        # Pig uncharacterized alpha-like
        "ENSSSCG00000007978": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA"},
        # Duck (Anas platyrhynchos) uncharacterized globins
        "ENSAPLG00020018293": {"category": "alpha_globin", "subtype": "HBA_like", "name": "HBA_like"},
        "ENSAPLG00020006106": {"category": "beta_globin", "subtype": "HBB_like", "name": "HBB_like_2"},
        "ENSAPLG00020006112": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "ENSAPLG00020006120": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        # Turkey (Meleagris gallopavo) uncharacterized alpha globins
        "ENSMGAG00000008793": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "ENSMGAG00000008803": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "ENSMGAG00000008823": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        # Ostrich (Struthio camelus australis) uncharacterized alpha globins
        "ENSSCUG00000008763": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "ENSSCUG00000008765": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        # Green Anole (Anolis carolinensis) uncharacterized beta globin (scaffold GL344018.1)
        "ENSACAG00000012173": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Box turtle (Terrapene carolina triunguis) uncharacterized alpha globins
        "ENSTMTG00000002687": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "ENSTMTG00000002694": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_2"},
        # Softshell turtle (Pelodiscus sinensis) uncharacterized beta globin
        "ENSPSIG00000017186": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Bearded Dragon (Pogona vitticeps) uncharacterized alpha globins
        "ENSPVIG00000000601": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "ENSPVIG00000001216": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_2"},
        # Turkey (Meleagris gallopavo) uncharacterized beta globins (chr 1)
        "ENSMGAG00000014946": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "ENSMGAG00000014944": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like_1"},
        "ENSMGAG00000014945": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "ENSMGAG00000021241": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like_2"},
        "ENSMGAG00000014964": {"category": "beta_globin", "subtype": "HBB", "name": "HBBA_like"},
        # Box turtle (Terrapene carolina triunguis) uncharacterized beta globin (PRFB02001019.1)
        "ENSTMTG00000010826": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Medaka (Oryzias latipes) uncharacterized globins on chr 8 (Locus A)
        "ENSORLG00000025529": {"category": "globin_coding_locus_A", "subtype": "hbb", "name": "hbb_like_1"},
        "ENSORLG00000023522": {"category": "globin_coding_locus_A", "subtype": "hbb", "name": "hbb_like_2"},
        "ENSORLG00000005267": {"category": "globin_coding_locus_A", "subtype": "hba", "name": "hba_like_1"},
        "ENSORLG00000027445": {"category": "globin_coding_locus_A", "subtype": "hba", "name": "hba_like_2"},
        "ENSORLG00000005283": {"category": "globin_coding_locus_A", "subtype": "hbb", "name": "hbb_like_3"},
        # Platypus (Ornithorhynchus anatinus) uncharacterized beta globins (Chr 2)
        "ENSOANG00000043521": {"category": "beta_globin", "subtype": "HBE", "name": "HBE"},
        "ENSOANG00000049640": {"category": "beta_globin", "subtype": "HBB", "name": "HBB"},
        "ENSOANG00000038856": {"category": "beta_globin", "subtype": "HBE_like", "name": "HBE_like"},
        # American Alligator (Alligator mississippiensis) NCBI RefSeq (rAllMis1)
        "GeneID:102559155": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC102559155": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:102559387": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "LOC102559387": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "GeneID:102559610": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "LOC102559610": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "GeneID:102563321": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "LOC102563321": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "GeneID:106738221": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "LOC106738221": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "GeneID:132249606": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "LOC132249606": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "GeneID:102576565": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_3"},
        "LOC102576565": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_3"},
        "GeneID:132243239": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_4"},
        "LOC132243239": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_4"},
        "GeneID:102559263": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_5"},
        "LOC102559263": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_5"},
        "GeneID:102559498": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_6"},
        "LOC102559498": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_6"},
        # Green Anole (Anolis carolinensis) NCBI RefSeq (rAnoCar3.1.pri)
        "GeneID:134294182": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC134294182": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:100556620": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "LOC100556620": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "GeneID:103282340": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "LOC103282340": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "GeneID:134297303": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "LOC134297303": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "GeneID:100552694": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC100552694": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Central Bearded Dragon (Pogona vitticeps) NCBI RefSeq (PviZW2.1)
        "GeneID:140702511": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC140702511": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:110089136": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "LOC110089136": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "GeneID:110089112": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "LOC110089112": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "GeneID:110078124": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "LOC110078124": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "GeneID:140701279": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC140701279": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Golden Eagle (Aquila chrysaetos) Ensembl (bAquChr1.2)
        "ENSACCG00020013324": {"category": "alpha_globin", "subtype": "HBA1", "name": "HBA1_like"},
        "ENSACCG00020010446": {"category": "beta_globin", "subtype": "HBE", "name": "HBE1_like"},
        "ENSACCG00020010456": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "ENSACCG00020010465": {"category": "beta_globin", "subtype": "HBB", "name": "HBRHO_like"},
        # Common Toad (Bufo bufo) NCBI RefSeq (aBufBuf1.1)
        "GeneID:121008078": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA3_like"},
        "GeneID:121008055": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA3_like_2"},
        "GeneID:121008040": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA4_like"},
        "GeneID:121007733": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA5_like_1"},
        "GeneID:121007734": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA5_like_2"},
        "GeneID:121007735": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA3_like_3"},
        "GeneID:121007736": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "GeneID:121008069": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "GeneID:121008018": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "GeneID:121008046": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_3"},
        "GeneID:121008071": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_4"},
        "GeneID:121007942": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_5"},
        "GeneID:121008074": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_6"},
        "GeneID:121007737": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_7"},
        # Common Wall Lizard (Podarcis muralis) NCBI RefSeq (rPodMur119.hap1.1)
        "GeneID:114584330": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "LOC114584330": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "GeneID:114584329": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "LOC114584329": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "GeneID:114595222": {"category": "beta_globin", "subtype": "HBB", "name": "HBB1_like"},
        "LOC114595222": {"category": "beta_globin", "subtype": "HBB", "name": "HBB1_like"},
        "GeneID:144327399": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "LOC144327399": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "GeneID:114595212": {"category": "beta_globin", "subtype": "HBB", "name": "HBB2_like"},
        "LOC114595212": {"category": "beta_globin", "subtype": "HBB", "name": "HBB2_like"},
        # Western Garter Snake (Thamnophis elegans) NCBI RefSeq (rThaEle1.pri)
        "GeneID:116517500": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "LOC116517500": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "GeneID:116517505": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "LOC116517505": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "GeneID:116510062": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "LOC116510062": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "GeneID:116510063": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "LOC116510063": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        # Green Sea Turtle (Chelonia mydas) NCBI RefSeq (rCheMyd1.pri.v2)
        "GeneID:102945818": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC102945818": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:102945589": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "LOC102945589": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "GeneID:102946728": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "LOC102946728": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "GeneID:119563550": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC119563550": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        # Three-spined Stickleback (Gasterosteus aculeatus) NCBI RefSeq (fGasAcu3.hap1.1)
        "GeneID:120827429": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_1"},
        "LOC120827429": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_1"},
        "GeneID:120828497": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_1"},
        "LOC120828497": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_1"},
        "GeneID:120827725": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_2"},
        "LOC120827725": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_2"},
        "GeneID:120827721": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_2"},
        "LOC120827721": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_2"},
        "GeneID:120827724": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_3"},
        "LOC120827724": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_3"},
        "GeneID:120827723": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_3"},
        "LOC120827723": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_3"},
        "GeneID:120827726": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_4"},
        "LOC120827726": {"category": "alpha_globin", "subtype": "hbae", "name": "hbae_4"},
        "GeneID:120827722": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_4"},
        "LOC120827722": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_4"},
        "GeneID:120827720": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_5"},
        "LOC120827720": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_5"},
        "GeneID:120827727": {"category": "alpha_globin", "subtype": "hba", "name": "hba_1"},
        "LOC120827727": {"category": "alpha_globin", "subtype": "hba", "name": "hba_1"},
        "GeneID:100174873": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_6"},
        "LOC100174873": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_6"},
        # Fugu (Takifugu rubripes) Ensembl (fTakRub1.2)
        "ENSTRUG00000016923": {"category": "beta_globin", "subtype": "hbb", "name": "hbb_dup"},
        # Small-spotted Catshark (Scyliorhinus canicula) NCBI RefSeq (sScyCan1.1)
        "GeneID:119978771": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "LOC119978771": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_1"},
        "GeneID:119978759": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_2"},
        "LOC119978759": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like_2"},
        "GeneID:119978731": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "LOC119978731": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "GeneID:119978783": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "LOC119978783": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        # Common Ostrich (Struthio camelus australis) NCBI RefSeq (StrCam1 / ASM69896v1)
        "GeneID:104143901": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "LOC104143901": {"category": "alpha_globin", "subtype": "HBAA", "name": "HBAA_like"},
        "GeneID:104143902": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "LOC104143902": {"category": "alpha_globin", "subtype": "HBAD", "name": "HBAD_like"},
        "GeneID:104143903": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC104143903": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:104141513": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like_fragment"},
        "LOC104141513": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like_fragment"},
        "GeneID:104141507": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "LOC104141507": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "GeneID:104141505": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC104141505": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "GeneID:104141514": {"category": "beta_globin", "subtype": "HBB", "name": "HBBA_like"},
        "LOC104141514": {"category": "beta_globin", "subtype": "HBB", "name": "HBBA_like"},
        # Naked Mole-Rat (Heterocephalus glaber) NCBI RefSeq (HetGla_female_1.0)
        "GeneID:101705305": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC101705305": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:101704937": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA1_like"},
        "LOC101704937": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA1_like"},
        "GeneID:101704581": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "LOC101704581": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "GeneID:101724100": {"category": "alpha_globin", "subtype": "HBQ", "name": "HBQ_like"},
        "LOC101724100": {"category": "alpha_globin", "subtype": "HBQ", "name": "HBQ_like"},
        "GeneID:101696428": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC101696428": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "GeneID:101718185": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "LOC101718185": {"category": "beta_globin", "subtype": "HBE", "name": "HBE_like"},
        "GeneID:101718529": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "LOC101718529": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        # California Sea Lion (Zalophus californianus) NCBI RefSeq (mZalCal1.pri.v2)
        "GeneID:113933389": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ2_like"},
        "LOC113933389": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ2_like"},
        "GeneID:113933387": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "LOC113933387": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "GeneID:113914335": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "LOC113914335": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "GeneID:113914336": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
        "LOC113914336": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
        "GeneID:118355949": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "LOC118355949": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_1"},
        "GeneID:113914337": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        "LOC113914337": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like_2"},
        # Bottlenose Dolphin (Tursiops truncatus) NCBI RefSeq (mTurTru1.mat.Y)
        "GeneID:117308132": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA1_like"},
        "LOC117308132": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA1_like"},
        "GeneID:117308133": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "LOC117308133": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA2_like"},
        "GeneID:117313420": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "LOC117313420": {"category": "beta_globin", "subtype": "HBB", "name": "HBB_like"},
        "GeneID:101335834": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
        "LOC101335834": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
        "GeneID:101336124": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "LOC101336124": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "GeneID:101315568": {"category": "beta_globin", "subtype": "HBE", "name": "HBE1_like"},
        "LOC101315568": {"category": "beta_globin", "subtype": "HBE", "name": "HBE1_like"},
        # African Savanna Elephant (Loxodonta africana) NCBI RefSeq (mLoxAfr1.hap2)
        "GeneID:100657757": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "LOC100657757": {"category": "alpha_globin", "subtype": "HBZ", "name": "HBZ_like"},
        "GeneID:135232945": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "LOC135232945": {"category": "alpha_globin", "subtype": "HBA", "name": "HBA_like"},
        "GeneID:100677792": {"category": "beta_globin", "subtype": "HBE", "name": "HBE1_like"},
        "LOC100677792": {"category": "beta_globin", "subtype": "HBE", "name": "HBE1_like"},
        "GeneID:100677510": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "LOC100677510": {"category": "beta_globin", "subtype": "HBG", "name": "HBG_like"},
        "GeneID:135232176": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
        "LOC135232176": {"category": "beta_globin", "subtype": "HBD", "name": "HBD_like"},
    }

    FLANKING_MARKERS = {
        "rhbdf1", "rhbdf1a", "rhbdf1b",
        "mpg",
        "nprl3",
        "luc7l", "luc7l2",
        "kank2",
        "lcmt1",
        "dchs1",
        "stim1", "stim1a", "stim1b",
        "snrnp25", "polr3k", "gtf3c1",
        "itfg3", "c16orf35",
        "capn1", "capn11", "capn12",
        "trim68",
        "ifitm1", "ifitm2", "ifitm3",
        "folh1",
        "fam234a", "rgs11", "arhgdig", "pdia2", "axin1", "pgap6", "rrm1",
    }

    def __init__(self, data_manager: Optional[DataManager] = None):
        self.dm = data_manager or DataManager()

    def classify_gene(self, gene: Dict[str, Any]) -> Dict[str, str]:
        """Classifies a gene feature into globin subtypes or flanking markers."""
        gene_id = str(gene.get("id") or "")
        raw_name = gene.get("external_name")
        if raw_name is None or (isinstance(raw_name, float) and str(raw_name) == "nan"):
            gene_name = gene_id
        else:
            gene_name = str(raw_name).strip() or gene_id
        name_lower = gene_name.lower()

        raw_desc = gene.get("description")
        if raw_desc is None or (isinstance(raw_desc, float) and str(raw_desc) == "nan"):
            desc = ""
        else:
            desc = str(raw_desc).lower()

        raw_bio = gene.get("biotype")
        if raw_bio is None or (isinstance(raw_bio, float) and str(raw_bio) == "nan"):
            biotype = ""
        else:
            biotype = str(raw_bio).strip()

        # 1. Check known unannotated lookup table (by ID or Symbol)
        if gene_id in self.KNOWN_UNANNOTATED_GLOBINS:
            info = self.KNOWN_UNANNOTATED_GLOBINS[gene_id]
            return {
                "category": info["category"],
                "subtype": info["subtype"],
                "display_symbol": info["name"],
                "is_globin": "True",
            }
        if gene_name in self.KNOWN_UNANNOTATED_GLOBINS:
            info = self.KNOWN_UNANNOTATED_GLOBINS[gene_name]
            return {
                "category": info["category"],
                "subtype": info["subtype"],
                "display_symbol": info["name"],
                "is_globin": "True",
            }

        # 2. Check if pseudogene or non-coding RNA
        if "pseudogene" in biotype:
            # Check if globin pseudogene
            if "globin" in desc or "hemoglobin" in desc or "hba" in name_lower or "hbb" in name_lower or "hbz" in name_lower:
                return {"category": "globin_pseudogene", "subtype": "globin_pseudogene", "display_symbol": gene_name, "is_globin": "Pseudogene"}
            return {"category": "other_pseudogene", "subtype": "other_pseudogene", "display_symbol": gene_name, "is_globin": "False"}

        if biotype != "protein_coding":
            return {"category": "non_coding", "subtype": biotype, "display_symbol": gene_name, "is_globin": "False"}

        # 3. Check Flanking markers, non-globin growth factors, and non-erythroid intracellular globins
        if (
            name_lower.startswith("hbegf")
            or "heparin-binding egf" in desc
            or any(w in desc or w in name_lower for w in ["cytoglobin", "neuroglobin", "myoglobin"])
            or name_lower in self.FLANKING_MARKERS
            or any(name_lower.startswith(m) for m in ["rhbdf1", "kank2", "lcmt1", "nprl3", "luc7l", "mpg", "capn1", "trim68", "or51", "or52", "olfr", "fam234", "rgs11", "arhgdig", "pdia2", "axin1", "pgap6", "rrm1"])
        ):
            sub_type = "CYGB" if ("cyto" in desc or "cyto" in name_lower) else gene_name.upper()
            return {"category": "cluster_non_globin", "subtype": sub_type, "display_symbol": gene_name, "is_globin": "False"}

        # 4. Check Alpha Globins
        # HBA, HBA1, HBA2, HBAD, HBZ, HBM, HBQ, hbae1, hbae3, hbae4, hbaa1, pi, etc.
        alpha_patterns = [
            r"^hb[azmq][1-5]?$", r"^hba[1-4]$", r"^hbad.*", r"^hba-", r"^hbae[1-5]$", r"^hbaa[1-2]$", r"^ba1l$",
            r"hemoglobin.*alpha", r"hemoglobin.*zeta", r"hemoglobin.*mu", r"hemoglobin.*theta", r"hemoglobin.*pi"
        ]
        for pat in alpha_patterns:
            if re.search(pat, name_lower) or re.search(pat, desc):
                return {"category": "alpha_globin", "subtype": gene_name, "display_symbol": gene_name, "is_globin": "True"}

        # 5. Check Beta Globins
        # HBB, HBD, HBE1, HBG1, HBG2, Hbb-b1, Hbb-bs, Hbb-bt, Hbb-y, Hbb-bh1, HBBA, rho, etc.
        beta_patterns = [
            r"^hbb$", r"^hbb[1-4]$", r"^hbd$", r"^hbe[1-4]$", r"^hbg[1-2]$", r"^hbb-", r"^hbba[1-2]$", r"^hbbe[1-3]$", r"^gb-beta",
            r"hemoglobin.*beta", r"hemoglobin.*gamma", r"hemoglobin.*epsilon", r"hemoglobin.*delta", r"hemoglobin.*rho"
        ]
        for pat in beta_patterns:
            if re.search(pat, name_lower) or re.search(pat, desc):
                return {"category": "beta_globin", "subtype": gene_name, "display_symbol": gene_name, "is_globin": "True"}

        # 6. Check generic globin keyword
        if "globin" in desc or "globin" in name_lower or "hemoglobin" in desc:
            return {"category": "other_globin", "subtype": "globin_family", "display_symbol": gene_name, "is_globin": "True"}

        # 7. Check if uncharacterized novel prediction (no approved gene symbol)
        if (gene_name.startswith("ENS") and "G0" in gene_name) or (gene_name.startswith("LOC") and "novel" in desc):
            return {"category": "other_locus_feature", "subtype": "novel_uncharacterized", "display_symbol": gene_name, "is_globin": "False"}

        # 8. Otherwise: other protein-coding gene in cluster (olfactory receptors, etc.)
        return {"category": "cluster_non_globin", "subtype": "intervening_gene", "display_symbol": gene_name, "is_globin": "False"}

    @staticmethod
    def get_canonical_transcript(details: Dict[str, Any], is_globin: bool = False) -> Optional[Dict[str, Any]]:
        """Identifies and returns the authentic canonical/principal transcript model from gene details.

        Guards against Ensembl automated chimeric read-through transcripts where an automated
        prediction fuses distant exons (e.g. >10-50 kb away) and gets tagged as canonical
        due to an artificial long ORF, while an authentic compact biological transcript exists.
        """
        if not details or "Transcript" not in details:
            return None
        transcripts = details.get("Transcript", [])
        if not transcripts:
            return None
        if len(transcripts) == 1:
            return transcripts[0]

        # Detect globin gene
        desc_lower = str(details.get("description") or "").lower()
        name_lower = str(details.get("display_name") or details.get("external_name") or "").lower()
        if not is_globin:
            if "globin" in desc_lower or "hemoglobin" in desc_lower or any(
                name_lower.startswith(p) for p in ["hb", "hba", "hbb", "hbd", "hbe", "hbg", "hbz", "hbm", "hbq"]
            ):
                is_globin = True

        pc_transcripts = [t for t in transcripts if t.get("biotype") == "protein_coding"]
        if not pc_transcripts:
            pc_transcripts = transcripts

        def get_t_stats(t):
            span = t.get("end", 0) - t.get("start", 0) + 1
            trans = t.get("Translation", {})
            cds_aa = trans.get("length", 0) if trans else 0
            if not cds_aa and t.get("CDS"):
                cds_aa = sum(c.get("end", 0) - c.get("start", 0) + 1 for c in t["CDS"]) // 3
            ex_count = len(t.get("Exon", []))
            is_canon = t.get("is_canonical") == 1
            return span, cds_aa, ex_count, is_canon

        # BIOLOGICAL GLOBIN LOGIC
        if is_globin:
            # Vertebrate hemoglobins: standard polypeptide ~140-155 aa (CDS ~420-465 bp), 3 exons, span ~1-3.5 kb
            candidates = []
            for t in pc_transcripts:
                span, cds_aa, ex_count, is_canon = get_t_stats(t)
                if cds_aa >= 100:  # functional globin chain
                    # Heavily penalize chimeric readthrough spans (> 4,000 bp)
                    span_penalty = max(0, span - 4000) / 500.0
                    # Prefer classic globin CDS length (~141-147 aa)
                    cds_diff = abs(cds_aa - 144)
                    # Prefer standard 3-exon structure
                    exon_diff = abs(ex_count - 3)
                    score = span_penalty * 100 + cds_diff * 2 + exon_diff * 5
                    if is_canon and span <= 4500:
                        score -= 10
                    candidates.append((score, span, cds_aa, t))
            if candidates:
                candidates.sort(key=lambda x: x[0])
                return candidates[0][3]

        # GENERAL LOGIC FOR ALL GENES (filtering out chimeric readthroughs)
        # 1. Match explicit canonical_transcript identifier, provided it is not an extreme read-through
        canon_id = details.get("canonical_transcript")
        if canon_id:
            canon_clean = str(canon_id).split(".")[0]
            for t in pc_transcripts:
                if str(t.get("id", "")).split(".")[0] == canon_clean:
                    span, cds_aa, _, _ = get_t_stats(t)
                    compact = [
                        t2 for t2 in pc_transcripts
                        if get_t_stats(t2)[0] < span * 0.25 and get_t_stats(t2)[1] >= 0.85 * cds_aa and get_t_stats(t2)[1] >= 100
                    ]
                    if not compact:
                        return t
                    return min(compact, key=lambda t2: get_t_stats(t2)[0])

        # 2. Match is_canonical == 1 flag
        for t in pc_transcripts:
            if t.get("is_canonical") == 1:
                span, cds_aa, _, _ = get_t_stats(t)
                compact = [
                    t2 for t2 in pc_transcripts
                    if get_t_stats(t2)[0] < span * 0.25 and get_t_stats(t2)[1] >= 0.85 * cds_aa and get_t_stats(t2)[1] >= 100
                ]
                if not compact:
                    return t
                return min(compact, key=lambda t2: get_t_stats(t2)[0])

        # 3. Choose protein-coding transcript with longest translation, breaking ties with shortest span
        def t_score(t):
            span, cds_aa, _, _ = get_t_stats(t)
            return (cds_aa, -span)
        return max(pc_transcripts, key=t_score)

    def annotate_interval(self, interval: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Fetches all genes in interval and attaches classification, coordinates, and structure."""
        species = interval["species"]
        chrom = interval["chromosome"]
        start = interval["start"]
        end = interval["end"]
        locus_id = interval.get("locus_id", "alpha_globin")
        clade = interval.get("clade", "unknown")

        raw_genes = self.dm.fetch_region_genes(species, chrom, start, end)
        annotated_genes = []

        # Sort genes by coordinate
        sorted_genes = sorted(raw_genes, key=lambda g: g.get("start", 0))

        for g in sorted_genes:
            cls_info = self.classify_gene(g)
            raw_bio = str(g.get("biotype") or "").strip().lower()

            # Filter strictly for protein_coding genes; exclude non-coding features and pseudogenes
            is_known = (g.get("id") in self.KNOWN_UNANNOTATED_GLOBINS) or (g.get("external_name") in self.KNOWN_UNANNOTATED_GLOBINS)
            is_pc = (raw_bio == "protein_coding") or is_known
            if not is_pc or "pseudogene" in raw_bio or cls_info["category"] in ["globin_pseudogene", "other_pseudogene", "non_coding", "other_locus_feature"]:
                continue

            g_start = g.get("start", 0) - 1  # convert to 0-based
            g_end = g.get("end", 0)
            strand = "+" if g.get("strand") == 1 else ("-" if g.get("strand") == -1 else ".")

            # Check canonical transcript to obtain accurate isoform boundaries
            g_id = g.get("id")
            try:
                details = self.dm.fetch_gene_details(g_id)
            except Exception:
                details = None

            if details:
                target_t = self.get_canonical_transcript(details, is_globin=(cls_info["is_globin"] == "True"))
                if target_t and target_t.get("start") and target_t.get("end"):
                    g_start = target_t.get("start") - 1
                    g_end = target_t.get("end")
                    if target_t.get("strand"):
                        strand = "+" if target_t.get("strand") == 1 else ("-" if target_t.get("strand") == -1 else strand)

            entry = {
                "species": species,
                "clade": clade,
                "locus_id": locus_id,
                "chromosome": chrom,
                "gene_id": g_id,
                "gene_name": g.get("external_name") or g.get("id"),
                "display_symbol": cls_info["display_symbol"],
                "biotype": g.get("biotype"),
                "category": cls_info["category"],
                "subtype": cls_info["subtype"],
                "is_globin": cls_info["is_globin"],
                "start": g_start,
                "end": g_end,
                "span_bp": g_end - g_start,
                "strand": strand,
                "description": g.get("description", "") or "",
            }
            annotated_genes.append(entry)

        return annotated_genes
