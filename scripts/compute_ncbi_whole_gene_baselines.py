"""Compute empirical whole-gene GC baselines for NCBI RefSeq vertebrate species.

Streams GFF to index exact genomic coordinates (start, end, chromosome) for all protein-coding genes,
then streams the genomic FASTA (fna.gz) chromosome-by-chromosome to compute true whole-gene GC
(including exons and introns). Merges with existing CDS GC / GC3 values, strictly eliminating
any whole_gene_gc = cds_gc fallback shortcut.
"""

import os
import re
import gzip
import json
import time
import urllib.request
import numpy as np
import pandas as pd
from typing import Dict, Tuple, Optional

CACHE_DIR = "data/genome_baselines"

NCBI_SPECIES = {
    "alligator_mississippiensis": {
        "assembly": "rAllMis1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_genomic.fna.gz",
        "genomic_mu": 44.0,
        "genomic_sigma": 4.1,
    },
    "anolis_carolinensis": {
        "assembly": "rAnoCar3.1.pri",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/594/765/GCF_035594765.1_rAnoCar3.1.pri/GCF_035594765.1_rAnoCar3.1.pri_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/594/765/GCF_035594765.1_rAnoCar3.1.pri/GCF_035594765.1_rAnoCar3.1.pri_genomic.fna.gz",
        "genomic_mu": 40.5,
        "genomic_sigma": 4.0,
    },
    "pogona_vitticeps": {
        "assembly": "PviZW2.1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/106/095/GCF_051106095.1_PviZW2.1/GCF_051106095.1_PviZW2.1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/106/095/GCF_051106095.1_PviZW2.1/GCF_051106095.1_PviZW2.1_genomic.fna.gz",
        "genomic_mu": 41.5,
        "genomic_sigma": 4.0,
    },
    "bufo_bufo": {
        "assembly": "aBufBuf1.1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/171/765/GCF_905171765.1_aBufBuf1.1/GCF_905171765.1_aBufBuf1.1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/171/765/GCF_905171765.1_aBufBuf1.1/GCF_905171765.1_aBufBuf1.1_genomic.fna.gz",
        "genomic_mu": 41.5,
        "genomic_sigma": 4.0,
    },
}


def calculate_gc(seq: str) -> Tuple[float, int]:
    if not seq:
        return 0.0, 0
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    valid = gc + seq.count("A") + seq.count("T")
    if valid == 0:
        return 0.0, 0
    return round((gc / valid) * 100.0, 2), valid


def process_ncbi_species(species: str, cfg: Dict, base_dir: str = CACHE_DIR):
    print(f"\n========================================================")
    print(f"Processing NCBI RefSeq: {species} ({cfg['assembly']})")
    print(f"========================================================")

    csv_path = os.path.join(base_dir, f"{species}_genes_gc.csv.gz")
    json_path = os.path.join(base_dir, f"{species}_stats.json")

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Existing baseline table not found: {csv_path}")

    # Load existing CDS GC / GC3 data
    df_existing = pd.read_csv(csv_path)
    print(f"[{species}] Loaded existing table with {len(df_existing):,} genes")

    # Step 1: Stream GFF to extract genomic coordinates for GeneID features
    t0 = time.time()
    gff_url = cfg["gff_url"]
    print(f"[{species}] Streaming GFF: {gff_url.split('/')[-1]}...")
    req_gff = urllib.request.Request(gff_url, headers={"User-Agent": "Mozilla/5.0"})

    genes_by_chrom: Dict[str, list] = {}
    gene_coords: Dict[str, Tuple[str, int, int]] = {}

    with urllib.request.urlopen(req_gff, timeout=60) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                s = line.decode("utf-8", errors="ignore")
                if s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) >= 9 and parts[2] == "gene":
                    attr = parts[8]
                    m = re.search(r"GeneID:(\d+)", attr)
                    if m:
                        gid = f"GeneID:{m.group(1)}"
                        chrom = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        # Keep the widest span if multiple gene records exist
                        if gid in gene_coords:
                            old_chrom, old_start, old_end = gene_coords[gid]
                            if old_chrom == chrom:
                                start = min(start, old_start)
                                end = max(end, old_end)
                        gene_coords[gid] = (chrom, start, end)

    for gid, (chrom, start, end) in gene_coords.items():
        genes_by_chrom.setdefault(chrom, []).append((gid, start, end))

    print(f"[{species}] Parsed GFF: {len(gene_coords):,} genes across {len(genes_by_chrom):,} contigs/chromosomes in {time.time()-t0:.2f}s")

    # Step 2: Stream genomic FASTA chromosome-by-chromosome to compute real whole-gene GC
    t1 = time.time()
    fna_url = cfg["fna_url"]
    print(f"[{species}] Streaming genomic FASTA: {fna_url.split('/')[-1]}...")
    req_fna = urllib.request.Request(fna_url, headers={"User-Agent": "Mozilla/5.0"})

    gene_whole_gc: Dict[str, Tuple[float, int]] = {}
    current_chrom = None
    seq_parts = []
    chrom_count = 0

    def process_chrom(chrom_name, parts):
        if not chrom_name or chrom_name not in genes_by_chrom:
            return
        chrom_seq = "".join(parts).upper()
        chr_len = len(chrom_seq)
        for gid, s, e in genes_by_chrom[chrom_name]:
            s_idx = max(0, s - 1)
            e_idx = min(chr_len, e)
            gseq = chrom_seq[s_idx:e_idx]
            gc_val, valid_len = calculate_gc(gseq)
            if valid_len > 0:
                gene_whole_gc[gid] = (gc_val, valid_len)

    with urllib.request.urlopen(req_fna, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                line_str = line.decode("utf-8", errors="ignore")
                if line_str.startswith(">"):
                    if current_chrom and seq_parts:
                        process_chrom(current_chrom, seq_parts)
                        chrom_count += 1
                        if chrom_count % 50 == 0:
                            print(f"  Processed {chrom_count} contigs, {len(gene_whole_gc):,} whole-gene GCs computed...")
                    current_chrom = line_str[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line_str.strip())

            # Final chromosome
            if current_chrom and seq_parts:
                process_chrom(current_chrom, seq_parts)
                chrom_count += 1

    print(f"[{species}] Completed FASTA streaming: {len(gene_whole_gc):,} whole-gene GCs computed across {chrom_count} contigs in {time.time()-t1:.2f}s")

    # Step 3: Update existing table with true whole-gene GC and accurate GFF coordinates
    updated_records = []
    whole_vals = []
    cds_vals = []
    gc3_vals = []

    for _, row in df_existing.iterrows():
        gid = row["gene_id"]
        c_gc = row.get("cds_gc")
        c_len = row.get("cds_length")
        gc3 = row.get("gc3")
        sym = row.get("gene_symbol")

        if gid in gene_whole_gc:
            w_gc, w_len = gene_whole_gc[gid]
            chrom, s, e = gene_coords[gid]
        else:
            # Gene not found in genomic FASTA / GFF
            w_gc = None
            w_len = None
            chrom = row.get("chromosome")
            s = row.get("start")
            e = row.get("end")

        updated_records.append({
            "species": species,
            "gene_id": gid,
            "gene_symbol": sym,
            "biotype": "protein_coding",
            "chromosome": chrom,
            "start": s,
            "end": e,
            "whole_gene_len": w_len,
            "whole_gene_gc": w_gc,
            "cds_length": c_len,
            "cds_gc": c_gc,
            "gc3": gc3,
        })

        if w_gc is not None and pd.notna(w_gc):
            whole_vals.append(w_gc)
        if c_gc is not None and pd.notna(c_gc):
            cds_vals.append(c_gc)
        if gc3 is not None and pd.notna(gc3):
            gc3_vals.append(gc3)

    df_updated = pd.DataFrame(updated_records)

    # Check difference between whole_gene_gc and cds_gc
    valid_both = df_updated.dropna(subset=["whole_gene_gc", "cds_gc"])
    mean_abs_diff = (valid_both["whole_gene_gc"] - valid_both["cds_gc"]).abs().mean()
    print(f"[{species}] Verification: {len(valid_both):,} genes with both whole-gene and CDS GC")
    print(f"[{species}] Mean absolute difference |whole_gene_gc - cds_gc| = {mean_abs_diff:.4f}%")

    if mean_abs_diff < 0.1:
        raise ValueError(f"ERROR: {species} whole_gene_gc and cds_gc are still nearly identical ({mean_abs_diff:.4f}%)! Fallback detected.")

    arr_w = np.array(whole_vals)
    arr_c = np.array(cds_vals)
    arr_3 = np.array(gc3_vals)

    stats = {
        "species": species,
        "assembly": cfg["assembly"],
        "num_genes": len(df_updated),
        "mean_genomic_gc": cfg["genomic_mu"],
        "std_genomic_gc": cfg["genomic_sigma"],
        # Whole-Gene Body Metrics (Empirical exons + introns)
        "num_whole_genes": len(whole_vals),
        "mean_whole_gene_gc": round(float(np.mean(arr_w)), 2),
        "std_whole_gene_gc": round(float(np.std(arr_w)), 2),
        "median_whole_gene_gc": round(float(np.median(arr_w)), 2),
        "q25_whole_gene_gc": round(float(np.percentile(arr_w, 25)), 2),
        "q75_whole_gene_gc": round(float(np.percentile(arr_w, 75)), 2),
        # Coding Sequence Metrics (Empirical CDS)
        "num_cds_genes": len(cds_vals),
        "mean_cds_gc": round(float(np.mean(arr_c)), 2),
        "std_cds_gc": round(float(np.std(arr_c)), 2),
        "median_cds_gc": round(float(np.median(arr_c)), 2),
        "q25_cds_gc": round(float(np.percentile(arr_c, 25)), 2),
        "q75_cds_gc": round(float(np.percentile(arr_c, 75)), 2),
        # Codon GC3 Metrics
        "num_gc3_genes": len(gc3_vals),
        "mean_gc3": round(float(np.mean(arr_3)), 2),
        "std_gc3": round(float(np.std(arr_3)), 2),
        "median_gc3": round(float(np.median(arr_3)), 2),
        "q25_gc3": round(float(np.percentile(arr_3, 25)), 2),
        "q75_gc3": round(float(np.percentile(arr_3, 75)), 2),
        # Exonic GC backwards compatibility
        "mean_exonic_gc": round(float(np.mean(arr_c)), 2),
        "std_exonic_gc": round(float(np.std(arr_c)), 2),
        "mean_gene_gc": round(float(np.mean(arr_w)), 2),
        "std_gene_gc": round(float(np.std(arr_w)), 2),
    }

    # Save to disk
    df_updated.to_csv(csv_path, index=False, compression="gzip")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[{species}] Successfully saved updated baseline table ({len(df_updated):,} genes) and stats.")
    print(f"  Whole-Gene GC: {stats['mean_whole_gene_gc']}% +/- {stats['std_whole_gene_gc']}% (median: {stats['median_whole_gene_gc']}%)")
    print(f"  Coding CDS GC: {stats['mean_cds_gc']}% +/- {stats['std_cds_gc']}% (median: {stats['median_cds_gc']}%)")
    print(f"  Codon GC3:     {stats['mean_gc3']}% +/- {stats['std_gc3']}% (median: {stats['median_gc3']}%)")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compute true whole-gene GC baselines for NCBI RefSeq genomes.")
    parser.add_argument("--species", type=str, default=None, help="Process single NCBI species (optional)")
    args = parser.parse_args()

    targets = [args.species] if args.species else list(NCBI_SPECIES.keys())
    for sp in targets:
        if sp not in NCBI_SPECIES:
            print(f"Species {sp} not recognized in NCBI_SPECIES")
            continue
        process_ncbi_species(sp, NCBI_SPECIES[sp])


if __name__ == "__main__":
    main()
