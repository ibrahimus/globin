#!/usr/bin/env python3
"""
Update data/centromere_positions.json with recovered centromere coordinates
for all newly expanded species (minnow, rat, puma, horses, birds, sea turtles,
caecilians, bats, fugu, and giant urodeles).
"""

import json
import os

CEN_PATH = "data/centromere_positions.json"

NEW_CENTROMERES = {
    # Rat (Rattus norvegicus): both Chr 1 and Chr 10 are telocentric/acrocentric
    "rattus_norvegicus::NC_141534.1": 2500000,
    "rattus_norvegicus::NC_141543.1": 2500000,

    # Puma (Puma concolor): B4 is subtelomeric p-arm, E1 is small metacentric
    "puma_concolor::NC_141682.1": 12320000,
    "puma_concolor::NC_141676.1": 23500000,

    # Fathead Minnow (Pimephales promelas): Cyprinid pericentromeric heterochromatin
    "pimephales_promelas::NC_142771.1": 46000000,
    "pimephales_promelas::NC_142785.1": 10200000,

    # Horses (Equus caballus & Equus asinus): Chr 13 submetacentric, Chr 7 metacentric
    "equus_caballus::NC_091696.1": 18360000,
    "equus_caballus::NC_091690.1": 44660000,
    "equus_asinus::NC_091803.1": 18000000,
    "equus_asinus::NC_091809.1": 44000000,

    # Birds (Anser cygnoides, Coturnix japonica, Sturnus vulgaris): Chr 1 macro, Chr 14 micro
    "anser_cygnoides::NC_089873.1": 55000000,
    "anser_cygnoides::NC_089887.1": 1000000,
    "coturnix_japonica::NC_140663.1": 52000000,
    "coturnix_japonica::NC_140675.1": 1000000,
    "sturnus_vulgaris::NC_142244.1": 50000000,
    "sturnus_vulgaris::NC_142258.1": 1000000,

    # Hawksbill Sea Turtle (Eretmochelys imbricata): syntenic to green turtle Chr 1 and Chr 10
    "eretmochelys_imbricata::NC_135572.1": 150000000,
    "eretmochelys_imbricata::NC_135581.1": 35000000,

    # Caecilians (Geotrypetes, Microcaecilia, Rhinatrema): mapped via pericentromeric gene deserts
    "geotrypetes_seraphini::NC_047094.1": 43300000,
    "geotrypetes_seraphini::NC_047088.1": 73500000,
    "microcaecilia_unicolor::NC_044038.1": 137400000,
    "microcaecilia_unicolor::NC_044037.1": 202500000,
    "rhinatrema_bivittatum::NC_042628.1": 30600000,

    # Natal Long-Fingered Bat (Miniopterus natalensis): Chr 10 and Chr 2
    "miniopterus_natalensis::NC_142083.1": 56700000,
    "miniopterus_natalensis::NC_142091.1": 2240000,

    # Fugu (Takifugu rubripes): Chr 1 and Chr 5
    "takifugu_rubripes::1": 6300000,
    "takifugu_rubripes::5": 4040000,

    # Giant Urodeles (Ambystoma mexicanum & Triturus cristatus)
    "ambystoma_mexicanum::NC_090924.1": 1247000000,
    "ambystoma_mexicanum::NC_090925.1": 1069000000,
    "triturus_cristatus::NC_140593.1": 721000000,
    "triturus_cristatus::NC_140590.1": 1312000000,

    # Marine Divers (Tursiops truncatus & Zalophus californianus)
    "tursiops_truncatus::NC_047048.1": 35000000,
    "tursiops_truncatus::NC_047041.1": 50580000,
    "zalophus_californianus::NC_045604.1": 35000000,
    "zalophus_californianus::NC_045605.1": 50800000,

    # African Savanna Elephant (Loxodonta africana)
    "loxodonta_africana::NC_087353.1": 35700000,
    "loxodonta_africana::NC_087348.1": 43500000,

    # Small-spotted Catshark (Scyliorhinus canicula)
    "scyliorhinus_canicula::NC_052160.1": 68000000,

    # Common Ostrich (Struthio camelus australis)
    "struthio_camelus_australis::NW_009271137.1": 100000,
    "struthio_camelus_australis::NW_009270861.1": 50000,

    # Naked Mole-Rat (Heterocephalus glaber)
    "heterocephalus_glaber::NW_004624913.1": 825000,
    "heterocephalus_glaber::NW_004624817.1": 4480000,

    # Squamate Snakes (Vipera berus, Candoia aspera, Notechis scutatus, Pantherophis guttatus)
    # Alpha loci on acrocentric microchromosomes (1.0 Mb tip); Beta loci on submetacentric macrochromosomes (35-38 Mb)
    "vipera_berus::NC_138216.1": 1000000,
    "vipera_berus::NC_138206.1": 35000000,
    "candoia_aspera::NC_086166.1": 1000000,
    "candoia_aspera::NC_086157.1": 35000000,
    "notechis_scutatus::NC_142020.1": 1000000,
    "notechis_scutatus::NC_142013.1": 35000000,
    "pantherophis_guttatus::NC_142076.1": 1000000,
    "pantherophis_guttatus::NC_142067.1": 38000000,
}


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cen_paths = [
        os.path.join(base_dir, "data", "centromere_positions.json"),
        os.path.join(base_dir, "1_data_from_mike", "data", "centromere_positions.json"),
    ]

    for p in cen_paths:
        data = {}
        if os.path.exists(p):
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)

        n_before = len(data)
        for k, v in NEW_CENTROMERES.items():
            data[k] = v
        n_after = len(data)

        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

        print(f"Updated {p}: {n_before} -> {n_after} centromere definitions.")


if __name__ == "__main__":
    main()
