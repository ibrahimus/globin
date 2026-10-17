#!/usr/bin/env python3
"""
Update data/chromosome_lengths.json with exact chromosome lengths from NCBI assembly reports
for all 25 Phase 1 chromosomes.
"""

import json
import os
import sys
import urllib.request
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.compute_phase1_expansion_baselines import PHASE1_SPECIES

CHR_LEN_PATH = "data/chromosome_lengths.json"
CEN_PATH = "data/centromere_positions.json"


def update_lengths():
    with open(CHR_LEN_PATH, "r", encoding="utf-8") as f:
        chr_lengths = json.load(f)

    df = pd.read_csv("data/curated_intervals.csv")
    target_chroms = {}
    for idx, row in df.iterrows():
        target_chroms.setdefault(row["species"], set()).add(str(row["chromosome"]))

    added = 0
    for sp, cfg in PHASE1_SPECIES.items():
        gff_url = cfg["gff_url"]
        base_url = gff_url.rsplit("/", 1)[0]
        asm_name = base_url.split("/")[-1]
        rpt_url = f"{base_url}/{asm_name}_assembly_report.txt"

        needed = target_chroms.get(sp, set())
        req = urllib.request.Request(rpt_url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp:
            for line in resp:
                s = line.decode("utf-8").strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split("\t")
                refseq_acc = parts[6]
                seq_len = int(parts[8])
                if refseq_acc in needed:
                    key = f"{sp}::{refseq_acc}"
                    chr_lengths[key] = seq_len
                    added += 1
                    print(f"Added {key} -> {seq_len:,} bp ({seq_len/1e6:.2f} Mb)")

    with open(CHR_LEN_PATH, "w", encoding="utf-8") as f:
        json.dump(chr_lengths, f, indent=2)

    print(f"\nSuccessfully updated {CHR_LEN_PATH}! Added {added} chromosome lengths.")

    # Validation against centromeres
    with open(CEN_PATH, "r", encoding="utf-8") as f:
        centromeres = json.load(f)

    print("\n--- Validating all 74 intervals for centromere containment (0 < cen_bp < chr_len) ---")
    valid_count = 0
    invalid = []
    for idx, row in df.iterrows():
        sp = row["species"]
        chrom = str(row["chromosome"])
        key = f"{sp}::{chrom}"
        clen = chr_lengths.get(key)
        cen_bp = centromeres.get(key)
        if clen and cen_bp and 0 < cen_bp < clen:
            valid_count += 1
        else:
            invalid.append((key, clen, cen_bp))

    print(f"Valid centromeres: {valid_count}/{len(df)}")
    if invalid:
        print("INVALID intervals:")
        for inv in invalid:
            print(f"  {inv[0]}: chr_len={inv[1]}, cen_bp={inv[2]}")
    else:
        print("ALL 74 intervals have valid centromere containment (0 < cen_bp < chr_len)!")


if __name__ == "__main__":
    update_lengths()
