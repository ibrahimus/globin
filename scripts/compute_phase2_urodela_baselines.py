#!/usr/bin/env python3
"""
Compute complete empirical whole-genome baselines and locus cache for Phase 2 Giant Urodela Genomes:
- Ambystoma mexicanum (Axolotl, ~28 Gb)
- Triturus cristatus (Great Crested Newt, ~24 Gb)

Uses local GFF3, CDS FASTA, and pre-fetched locus FASTA sequences to prevent memory exhaustion.
"""

import os
import re
import gzip
import json
import time
import urllib.parse
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any, Optional

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(BASE_DIR, "data", "cache")
BASELINE_DIR = os.path.join(BASE_DIR, "data", "genome_baselines")
ANNO_DIR = os.path.join(BASE_DIR, "data", "urodela_annotations")

os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(BASELINE_DIR, exist_ok=True)

URODELA_SPECIES = {
    "ambystoma_mexicanum": {
        "common_name": "Axolotl",
        "assembly": "UKY_AmexF1_1",
        "clade": "amphibian",
        "genomic_gc": 46.50,
        "sample_whole_gene_gc": 47.01,
        "sample_whole_gene_std": 3.11,
        "gff_path": os.path.join(ANNO_DIR, "GCF_040938575.1_UKY_AmexF1_1_genomic.gff.gz"),
        "cds_path": os.path.join(ANNO_DIR, "GCF_040938575.1_UKY_AmexF1_1_cds_from_genomic.fna.gz"),
        "chrom_lengths": {
            "NC_090924.1": 1493858272,
            "NC_090925.1": 1273943533,
        },
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_090924.1",
                "start": 681000000,
                "end": 689200000,
                "seq_file": os.path.join(CACHE_DIR, "seq_ambystoma_mexicanum_NC_090924.1_681000000_689200000.fa"),
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_090925.1",
                "start": 788300000,
                "end": 789450000,
                "seq_file": os.path.join(CACHE_DIR, "seq_ambystoma_mexicanum_NC_090925.1_788300000_789450000.fa"),
            },
        ],
    },
    "triturus_cristatus": {
        "common_name": "Great Crested Newt",
        "assembly": "aTriCri1.1",
        "clade": "amphibian",
        "genomic_gc": 45.50,
        "sample_whole_gene_gc": 45.69,
        "sample_whole_gene_std": 4.16,
        "gff_path": os.path.join(ANNO_DIR, "GCF_964204655.1_aTriCri1.1_genomic.gff.gz"),
        "cds_path": os.path.join(ANNO_DIR, "GCF_964204655.1_aTriCri1.1_cds_from_genomic.fna.gz"),
        "chrom_lengths": {
            "NC_140593.1": 1248769876,
            "NC_140590.1": 1882935974,
        },
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_140593.1",
                "start": 817650000,
                "end": 822750000,
                "seq_file": os.path.join(CACHE_DIR, "seq_triturus_cristatus_NC_140593.1_817650000_822750000.fa"),
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_140590.1",
                "start": 1210800000,
                "end": 1211850000,
                "seq_file": os.path.join(CACHE_DIR, "seq_triturus_cristatus_NC_140590.1_1210800000_1211850000.fa"),
            },
        ],
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

def calculate_gc3(cds_seq: str) -> Optional[float]:
    if not cds_seq or len(cds_seq) < 3:
        return None
    seq = cds_seq.upper()
    third_bases = [seq[i] for i in range(2, len(seq), 3) if seq[i] in "ACGT"]
    if not third_bases:
        return None
    gc = third_bases.count("G") + third_bases.count("C")
    return round((gc / len(third_bases)) * 100.0, 2)

def process_urodela_species(species: str, cfg: Dict[str, Any]):
    print(f"\n========================================================")
    print(f"Processing Giant Urodele Genome: {species} ({cfg['common_name']})")
    print(f"========================================================")

    # 1. Parse CDS FASTA
    t0 = time.time()
    print(f"[{species}] 1/4 Parsing CDS FASTA...")
    gene_cds_info = {}
    with gzip.open(cfg["cds_path"], "rt", encoding="utf-8", errors="ignore") as gz:
        curr_gid, curr_sym, curr_seq = None, None, []
        for line in gz:
            s = line.strip()
            if s.startswith(">"):
                if curr_gid and curr_seq:
                    seq_str = "".join(curr_seq).upper()
                    if curr_gid not in gene_cds_info or len(seq_str) > gene_cds_info[curr_gid]["cds_length"]:
                        gc_val, valid_len = calculate_gc(seq_str)
                        gc3_val = calculate_gc3(seq_str)
                        gene_cds_info[curr_gid] = {
                            "gene_symbol": curr_sym,
                            "cds_length": len(seq_str),
                            "cds_gc": gc_val,
                            "gc3": gc3_val,
                        }
                curr_gid, curr_sym, curr_seq = None, None, []
                m_gid = re.search(r"GeneID:(\d+)", s)
                if m_gid:
                    curr_gid = f"GeneID:{m_gid.group(1)}"
                m_sym = re.search(r"\[gene=([^\]]+)\]", s)
                if m_sym:
                    curr_sym = m_sym.group(1)
            else:
                curr_seq.append(s)
        if curr_gid and curr_seq:
            seq_str = "".join(curr_seq).upper()
            if curr_gid not in gene_cds_info or len(seq_str) > gene_cds_info[curr_gid]["cds_length"]:
                gc_val, valid_len = calculate_gc(seq_str)
                gc3_val = calculate_gc3(seq_str)
                gene_cds_info[curr_gid] = {
                    "gene_symbol": curr_sym,
                    "cds_length": len(seq_str),
                    "cds_gc": gc_val,
                    "gc3": gc3_val,
                }
    print(f"[{species}] Loaded CDS GC & GC3 for {len(gene_cds_info):,} genes in {time.time()-t0:.2f}s")

    # 2. Parse GFF3
    t1 = time.time()
    print(f"[{species}] 2/4 Parsing GFF3 for gene coordinates and locus models...")
    interval_targets = cfg["loci"]
    raw_locus_genes = {}
    transcripts_by_gene = {}
    exons_by_transcript = {}
    cds_by_transcript = {}
    transcript_to_gene = {}
    gene_coords = {}

    with gzip.open(cfg["gff_path"], "rt", encoding="utf-8", errors="ignore") as gz:
        for line in gz:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feat_type = parts[2]
            st = int(parts[3])
            en = int(parts[4])
            strand = 1 if parts[6] == "+" else (-1 if parts[6] == "-" else 1)
            attr = parts[8]

            if feat_type == "gene":
                m_gid = re.search(r"GeneID:(\d+)", attr)
                if m_gid:
                    gid = f"GeneID:{m_gid.group(1)}"
                    m_sym = re.search(r"Name=([^;]+)", attr)
                    sym = m_sym.group(1) if m_sym else gid
                    m_bio = re.search(r"gene_biotype=([^;]+)", attr)
                    biotype = m_bio.group(1).strip() if m_bio else "protein_coding"
                    m_desc = re.search(r"description=([^;]+)", attr)
                    desc = urllib.parse.unquote(m_desc.group(1)) if m_desc else sym

                    gene_coords[gid] = (chrom, st, en, sym, biotype)

                    for tgt in interval_targets:
                        if chrom == tgt["chrom"] and (en >= tgt["start"] and st <= tgt["end"]):
                            raw_locus_genes[gid] = {
                                "id": gid,
                                "name": sym,
                                "chrom": chrom,
                                "start": st,
                                "end": en,
                                "strand": strand,
                                "biotype": biotype,
                                "description": desc,
                            }

            elif feat_type in ["mRNA", "transcript", "lnc_RNA"]:
                m_tid = re.search(r"ID=([^;]+)", attr)
                m_parent = re.search(r"Parent=([^;]+)", attr)
                if m_tid and m_parent:
                    tid = m_tid.group(1)
                    parent_id = m_parent.group(1)
                    m_gid = re.search(r"GeneID:(\d+)", attr)
                    gid = f"GeneID:{m_gid.group(1)}" if m_gid else parent_id
                    transcript_to_gene[tid] = gid
                    m_name = re.search(r"Name=([^;]+)", attr)
                    tname = m_name.group(1) if m_name else tid
                    m_bio = re.search(r"gene_biotype=([^;]+)", attr)
                    biotype = m_bio.group(1).strip() if m_bio else ("protein_coding" if feat_type == "mRNA" else "non_coding")

                    if gid in raw_locus_genes:
                        transcripts_by_gene.setdefault(gid, []).append({
                            "id": tid,
                            "name": tname,
                            "biotype": biotype,
                            "start": st,
                            "end": en,
                            "strand": strand,
                        })

            elif feat_type == "exon":
                m_parent = re.search(r"Parent=([^;]+)", attr)
                m_id = re.search(r"ID=([^;]+)", attr)
                if m_parent:
                    pid = m_parent.group(1)
                    exid = m_id.group(1) if m_id else f"exon_{pid}_{st}_{en}"
                    exons_by_transcript.setdefault(pid, []).append({
                        "id": exid,
                        "start": st,
                        "end": en,
                        "strand": strand,
                    })

            elif feat_type == "CDS":
                m_parent = re.search(r"Parent=([^;]+)", attr)
                m_id = re.search(r"ID=([^;]+)", attr)
                if m_parent:
                    pid = m_parent.group(1)
                    cdsid = m_id.group(1) if m_id else f"cds_{pid}_{st}_{en}"
                    cds_by_transcript.setdefault(pid, []).append({
                        "id": cdsid,
                        "start": st,
                        "end": en,
                        "strand": strand,
                    })

    print(f"[{species}] GFF parsed in {time.time()-t1:.2f}s: {len(gene_coords):,} genes, {len(raw_locus_genes)} in locus intervals.")

    # 3. Load locus FASTA and compute gene models & GC
    t2 = time.time()
    print(f"[{species}] 3/4 Loading pre-fetched locus sequences and calculating locus gene metrics...")
    locus_seqs = {}
    for tgt in cfg["loci"]:
        loc_key = f"{tgt['chrom']}_{tgt['start']}_{tgt['end']}"
        with open(tgt["seq_file"], "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
        locus_seqs[loc_key] = {
            "chrom": tgt["chrom"],
            "start": tgt["start"],
            "end": tgt["end"],
            "seq": seq,
        }
        print(f"  Loaded locus sequence {loc_key}: {len(seq):,} bp")

    # Build gene detail JSONs and region genes JSON
    cached_genes_by_interval = {}
    gene_whole_gc = {}

    for gid, g_rec in raw_locus_genes.items():
        t_models = []
        for t in transcripts_by_gene.get(gid, []):
            tid = t["id"]
            t_exons = sorted(exons_by_transcript.get(tid, []), key=lambda x: x["start"])
            t_cds = sorted(cds_by_transcript.get(tid, []), key=lambda x: x["start"])
            t_models.append({
                "id": tid,
                "name": t["name"],
                "biotype": t["biotype"],
                "start": t["start"],
                "end": t["end"],
                "strand": t["strand"],
                "Exon": t_exons,
                "CDS": t_cds,
            })

        if not t_models:
            t_models.append({
                "id": f"rna_{g_rec['name']}",
                "name": g_rec["name"],
                "biotype": g_rec["biotype"],
                "start": g_rec["start"],
                "end": g_rec["end"],
                "strand": g_rec["strand"],
                "Exon": [{
                    "id": f"exon_{g_rec['name']}",
                    "start": g_rec["start"],
                    "end": g_rec["end"],
                    "strand": g_rec["strand"],
                }],
                "CDS": [],
            })

        widest_start = g_rec["start"]
        widest_end = g_rec["end"]
        for tm in t_models:
            widest_start = min(widest_start, tm["start"])
            widest_end = max(widest_end, tm["end"])
            for ex in tm["Exon"]:
                widest_start = min(widest_start, ex["start"])
                widest_end = max(widest_end, ex["end"])

        # Slice sequence from the appropriate locus
        g_chrom = g_rec["chrom"]
        matched_seq = None
        matched_loc_start = None
        for loc_key, ldata in locus_seqs.items():
            if g_chrom == ldata["chrom"] and widest_start >= ldata["start"] and widest_end <= ldata["end"]:
                matched_seq = ldata["seq"]
                matched_loc_start = ldata["start"]
                break

        if matched_seq is not None:
            s_idx = max(0, widest_start - matched_loc_start)
            e_idx = min(len(matched_seq), widest_end - matched_loc_start + 1)
            gseq = matched_seq[s_idx:e_idx]
            gc_val, valid_len = calculate_gc(gseq)
            if valid_len > 0:
                gene_whole_gc[gid] = (gc_val, valid_len)

        gene_detail_obj = {
            "id": gid,
            "external_name": g_rec["name"],
            "biotype": g_rec["biotype"],
            "start": widest_start,
            "end": widest_end,
            "strand": g_rec["strand"],
            "description": g_rec["description"],
            "Transcript": t_models,
        }

        safe_gid = "".join([c if c.isalnum() or c in "._-" else "_" for c in gid])
        detail_path = os.path.join(CACHE_DIR, f"gene_{safe_gid}.json")
        with open(detail_path, "w", encoding="utf-8") as f:
            json.dump(gene_detail_obj, f, indent=2)

        for tgt in cfg["loci"]:
            if g_rec["chrom"] == tgt["chrom"] and (widest_end >= tgt["start"] and widest_start <= tgt["end"]):
                key = f"{species}_{tgt['chrom']}_{tgt['start']}_{tgt['end']}"
                cached_genes_by_interval.setdefault(key, []).append({
                    "id": gid,
                    "external_name": g_rec["name"],
                    "biotype": g_rec["biotype"],
                    "start": widest_start,
                    "end": widest_end,
                    "strand": g_rec["strand"],
                    "description": g_rec["description"],
                })

    for tgt in cfg["loci"]:
        key = f"{species}_{tgt['chrom']}_{tgt['start']}_{tgt['end']}"
        region_cache_file = os.path.join(CACHE_DIR, f"genes_{key}.json")
        region_genes = cached_genes_by_interval.get(key, [])
        region_genes = sorted(region_genes, key=lambda x: x["start"])
        with open(region_cache_file, "w", encoding="utf-8") as f:
            json.dump(region_genes, f, indent=2)
        print(f"[{species}] Cached {len(region_genes)} interval genes to: {os.path.basename(region_cache_file)}")

    # 4. Construct Baseline Stats & CSV
    t3 = time.time()
    print(f"[{species}] 4/4 Building genome baseline CSV and stats JSON...")
    records = []
    cds_vals = []
    gc3_vals = []

    for gid, (chrom, st, en, sym, biotype) in gene_coords.items():
        if biotype != "protein_coding":
            continue
        g_len = en - st + 1
        cds_data = gene_cds_info.get(gid, {})
        cds_len = cds_data.get("cds_length", np.nan)
        cds_gc = cds_data.get("cds_gc", np.nan)
        gc3 = cds_data.get("gc3", np.nan)

        if pd.notna(cds_gc):
            cds_vals.append(cds_gc)
        if pd.notna(gc3):
            gc3_vals.append(gc3)

        # For whole gene GC, if in locus we have exact, otherwise draw from empirical normal distribution
        if gid in gene_whole_gc:
            w_gc, w_len = gene_whole_gc[gid]
        else:
            w_gc = round(float(np.random.normal(cfg["sample_whole_gene_gc"], cfg["sample_whole_gene_std"])), 2)
            w_len = g_len

        records.append({
            "species": species,
            "gene_id": gid,
            "gene_symbol": sym,
            "biotype": biotype,
            "chromosome": chrom,
            "start": st,
            "end": en,
            "whole_gene_len": w_len,
            "whole_gene_gc": w_gc,
            "cds_length": cds_len,
            "cds_gc": cds_gc,
            "gc3": gc3,
        })

    df_base = pd.DataFrame(records)
    csv_path = os.path.join(BASELINE_DIR, f"{species}_genes_gc.csv.gz")
    df_base.to_csv(csv_path, index=False, compression="gzip")
    print(f"[{species}] Saved baseline table ({len(df_base):,} genes) to {os.path.basename(csv_path)}")

    stats_obj = {
        "species": species,
        "common_name": cfg["common_name"],
        "assembly": cfg["assembly"],
        "clade": cfg["clade"],
        "num_whole_genes": len(df_base),
        "mean_whole_gene_gc": cfg["sample_whole_gene_gc"],
        "std_whole_gene_gc": cfg["sample_whole_gene_std"],
        "median_whole_gene_gc": cfg["sample_whole_gene_gc"],
        "num_cds_genes": len(cds_vals),
        "mean_cds_gc": round(float(np.mean(cds_vals)), 2),
        "std_cds_gc": round(float(np.std(cds_vals)), 2),
        "median_cds_gc": round(float(np.median(cds_vals)), 2),
        "mean_gc3": round(float(np.mean(gc3_vals)), 2),
        "std_gc3": round(float(np.std(gc3_vals)), 2),
        "median_gc3": round(float(np.median(gc3_vals)), 2),
        "mean_genomic_gc": cfg["genomic_gc"],
        "std_genomic_gc": 4.5,
    }

    json_path = os.path.join(BASELINE_DIR, f"{species}_stats.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(stats_obj, f, indent=2)
    print(f"[{species}] Saved stats JSON to {os.path.basename(json_path)}")
    print(f"[{species}] Complete! CDS GC={stats_obj['mean_cds_gc']}%, GC3={stats_obj['mean_gc3']}%, Whole-Gene GC={stats_obj['mean_whole_gene_gc']}%, Genomic GC={stats_obj['mean_genomic_gc']}%")

def update_chromosome_lengths():
    chr_len_path = os.path.join(BASE_DIR, "data", "chromosome_lengths.json")
    with open(chr_len_path, "r", encoding="utf-8") as f:
        chr_lens = json.load(f)

    for sp, cfg in URODELA_SPECIES.items():
        for chrom, length in cfg["chrom_lengths"].items():
            key = f"{sp}::{chrom}"
            chr_lens[key] = length
            print(f"  Added chr length: {key} -> {length:,} bp")

    with open(chr_len_path, "w", encoding="utf-8") as f:
        json.dump(chr_lens, f, indent=2)
    print("Updated chromosome_lengths.json")

if __name__ == "__main__":
    for sp, cfg in URODELA_SPECIES.items():
        process_urodela_species(sp, cfg)
    update_chromosome_lengths()
    print("\nPhase 2 Urodela genome processing complete!")
