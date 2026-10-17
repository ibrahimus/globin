#!/usr/bin/env python3
"""
Integrate 6 Strategic Expansion Genomes into the Vertebrate Globin Atlas:
1. Tursiops truncatus (Bottlenose Dolphin, Mammalia: Cetacea) - Marine diver
2. Zalophus californianus (California Sea Lion, Mammalia: Pinnipedia) - Carnivoran diver
3. Loxodonta africana (African Savanna Elephant, Mammalia: Afrotheria) - Missing superorder
4. Scyliorhinus canicula (Small-spotted Catshark, Chondrichthyes) - Jawed vertebrate linked origin
5. Struthio camelus australis (Common Ostrich, Aves: Palaeognathae) - Avian polarity inversion root
6. Heterocephalus glaber (Naked Mole-Rat, Mammalia: Rodentia) - Subterranean hypoxia extremophile

Streams official NCBI RefSeq:
1. CDS FASTA (cds_from_genomic.fna.gz) for canonical CDS GC, CDS length, and codon GC3
2. GFF3 (genomic.gff.gz) for gene coordinates, models, and locus intervals
3. Genomic DNA FASTA (genomic.fna.gz) chromosome-by-chromosome for true whole-gene GC and locus sequences

Outputs:
- data/genome_baselines/{species}_genes_gc.csv.gz
- data/genome_baselines/{species}_stats.json
- data/cache/gene_{gid}.json
- data/cache/genes_{key}.json
- data/cache/seq_{key}.fa
"""

import os
import re
import gzip
import json
import time
import urllib.request
import urllib.parse
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any, Optional

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(BASE_DIR, "data", "cache")
BASELINE_DIR = os.path.join(BASE_DIR, "data", "genome_baselines")
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(BASELINE_DIR, exist_ok=True)

EXPANSION_SPECIES = {
    "tursiops_truncatus": {
        "common_name": "Bottlenose Dolphin",
        "assembly": "mTurTru1.mat.Y",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595/GCF_011762595.2_mTurTru1.mat.Y/GCF_011762595.2_mTurTru1.mat.Y_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595/GCF_011762595.2_mTurTru1.mat.Y/GCF_011762595.2_mTurTru1.mat.Y_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595/GCF_011762595.2_mTurTru1.mat.Y/GCF_011762595.2_mTurTru1.mat.Y_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_globin", "chrom": "NC_047048.1", "start": 38080000, "end": 38200000, "strand": "-", "left_gene": "RHBDF1", "right_gene": "HBQ1"},
            {"locus_id": "beta_globin", "chrom": "NC_047041.1", "start": 59880000, "end": 59940000, "strand": "+", "left_gene": "HBB", "right_gene": "HBE1"},
        ],
    },
    "zalophus_californianus": {
        "common_name": "California Sea Lion",
        "assembly": "mZalCal1.pri.v2",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/762/305/GCF_009762305.2_mZalCal1.pri.v2/GCF_009762305.2_mZalCal1.pri.v2_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/762/305/GCF_009762305.2_mZalCal1.pri.v2/GCF_009762305.2_mZalCal1.pri.v2_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/762/305/GCF_009762305.2_mZalCal1.pri.v2/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_globin", "chrom": "NC_045604.1", "start": 69380000, "end": 69500000, "strand": "-", "left_gene": "RHBDF1", "right_gene": "HBQ1"},
            {"locus_id": "beta_globin", "chrom": "NC_045605.1", "start": 61010000, "end": 61060000, "strand": "+", "left_gene": "HBB", "right_gene": "HBE1"},
        ],
    },
    "loxodonta_africana": {
        "common_name": "African Savanna Elephant",
        "assembly": "mLoxAfr1.hap2",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/014/295/GCF_030014295.1_mLoxAfr1.hap2/GCF_030014295.1_mLoxAfr1.hap2_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/014/295/GCF_030014295.1_mLoxAfr1.hap2/GCF_030014295.1_mLoxAfr1.hap2_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/014/295/GCF_030014295.1_mLoxAfr1.hap2/GCF_030014295.1_mLoxAfr1.hap2_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_globin", "chrom": "NC_087353.1", "start": 53580000, "end": 53720000, "strand": "-", "left_gene": "RHBDF1", "right_gene": "HBQ1"},
            {"locus_id": "beta_globin", "chrom": "NC_087348.1", "start": 82900000, "end": 82950000, "strand": "+", "left_gene": "HBB", "right_gene": "HBE"},
        ],
    },
    "scyliorhinus_canicula": {
        "common_name": "Small-spotted Catshark",
        "assembly": "sScyCan1.1",
        "clade": "chondrichthyes",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/713/615/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/713/615/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/713/615/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_beta_cluster", "chrom": "NC_052160.1", "start": 73580000, "end": 73860000, "strand": "+", "left_gene": "LOC119978771", "right_gene": "LOC119978783"},
        ],
    },
    "struthio_camelus_australis": {
        "common_name": "Common Ostrich",
        "assembly": "ASM69896v1",
        "clade": "amniote_bird",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/698/965/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/698/965/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/698/965/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_globin", "chrom": "NW_009271137.1", "start": 860000, "end": 1050000, "strand": "-", "left_gene": "RHBDF1", "right_gene": "LOC104143901"},
            {"locus_id": "beta_globin", "chrom": "NW_009270861.1", "start": 1000, "end": 30000, "strand": "+", "left_gene": "LOC104141507", "right_gene": "LOC104141514"},
        ],
    },
    "heterocephalus_glaber": {
        "common_name": "Naked Mole-Rat",
        "assembly": "HetGla_female_1.0",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/695/GCF_000247695.1_HetGla_female_1.0/GCF_000247695.1_HetGla_female_1.0_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/695/GCF_000247695.1_HetGla_female_1.0/GCF_000247695.1_HetGla_female_1.0_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/695/GCF_000247695.1_HetGla_female_1.0/GCF_000247695.1_HetGla_female_1.0_genomic.fna.gz",
        "loci": [
            {"locus_id": "alpha_globin", "chrom": "NW_004624913.1", "start": 950000, "end": 1150000, "strand": "-", "left_gene": "Rhbdf1", "right_gene": "Hba"},
            {"locus_id": "beta_globin", "chrom": "NW_004624817.1", "start": 6650000, "end": 6750000, "strand": "+", "left_gene": "Hbb", "right_gene": "Hbe1"},
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


def process_species(species: str, cfg: Dict[str, Any], force: bool = False):
    csv_path = os.path.join(BASELINE_DIR, f"{species}_genes_gc.csv.gz")
    json_path = os.path.join(BASELINE_DIR, f"{species}_stats.json")

    has_loci_cache = True
    for loc in cfg["loci"]:
        k = f"{species}_{loc['chrom']}_{loc['start']}_{loc['end']}"
        rf = os.path.join(CACHE_DIR, f"genes_{k}.json")
        sf = os.path.join(CACHE_DIR, f"seq_{k}.fa")
        if not (os.path.exists(rf) and os.path.exists(sf)):
            has_loci_cache = False
            break

    if not force and os.path.exists(csv_path) and os.path.exists(json_path) and has_loci_cache:
        print(f"[{species}] Already cached at {csv_path} with all locus caches. Skipping baseline recompute.")
        return

    print(f"\n========================================================")
    print(f"Processing NCBI Expansion Genome: {species} ({cfg['assembly']})")
    print(f"========================================================")

    # ---------------------------------------------------------
    # Step 1: Stream CDS FASTA for canonical CDS GC and GC3
    # ---------------------------------------------------------
    t0 = time.time()
    cds_url = cfg["cds_url"]
    print(f"[{species}] 1/4 Streaming CDS FASTA: {cds_url.split('/')[-1]}...")
    req_cds = urllib.request.Request(cds_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Globin Project)"})

    gene_cds_info: Dict[str, Dict[str, Any]] = {}

    with urllib.request.urlopen(req_cds, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            curr_gid = None
            curr_sym = None
            curr_seq = []

            for line in gz:
                s = line.decode("utf-8", errors="ignore").strip()
                if s.startswith(">"):
                    if curr_gid and curr_seq:
                        seq_str = "".join(curr_seq)
                        seq_len = len(seq_str)
                        if seq_len > 0:
                            if curr_gid not in gene_cds_info or seq_len > gene_cds_info[curr_gid]["cds_length"]:
                                gc_val, valid_len = calculate_gc(seq_str)
                                gc3_val = calculate_gc3(seq_str)
                                gene_cds_info[curr_gid] = {
                                    "gene_symbol": curr_sym,
                                    "cds_length": seq_len,
                                    "cds_gc": gc_val,
                                    "gc3": gc3_val,
                                }
                    curr_gid = None
                    curr_sym = None
                    curr_seq = []

                    m_gid = re.search(r"GeneID:(\d+)", s)
                    if m_gid:
                        curr_gid = f"GeneID:{m_gid.group(1)}"
                    m_sym = re.search(r"\[gene=([^\]]+)\]", s)
                    if m_sym:
                        curr_sym = m_sym.group(1)
                else:
                    curr_seq.append(s)

            # Final record
            if curr_gid and curr_seq:
                seq_str = "".join(curr_seq)
                seq_len = len(seq_str)
                if seq_len > 0:
                    if curr_gid not in gene_cds_info or seq_len > gene_cds_info[curr_gid]["cds_length"]:
                        gc_val, valid_len = calculate_gc(seq_str)
                        gc3_val = calculate_gc3(seq_str)
                        gene_cds_info[curr_gid] = {
                            "gene_symbol": curr_sym,
                            "cds_length": seq_len,
                            "cds_gc": gc_val,
                            "gc3": gc3_val,
                        }

    print(f"[{species}] Extracted CDS GC & GC3 for {len(gene_cds_info):,} genes in {time.time()-t0:.2f}s")

    # ---------------------------------------------------------
    # Step 2: Stream GFF3 for gene coordinates and locus models
    # ---------------------------------------------------------
    t1 = time.time()
    gff_url = cfg["gff_url"]
    print(f"[{species}] 2/4 Streaming GFF3: {gff_url.split('/')[-1]}...")
    req_gff = urllib.request.Request(gff_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Globin Project)"})

    gene_coords: Dict[str, Tuple[str, int, int, str]] = {}
    genes_by_chrom: Dict[str, List[Tuple[str, int, int]]] = {}

    interval_targets = []
    for loc in cfg["loci"]:
        interval_targets.append({
            "locus_id": loc["locus_id"],
            "chrom": loc["chrom"],
            "orig_start": loc["start"],
            "orig_end": loc["end"],
            "filter_start": max(0, loc["start"] - 80000),
            "filter_end": loc["end"] + 80000,
        })

    raw_locus_genes: Dict[str, Dict[str, Any]] = {}
    transcripts_by_gene: Dict[str, List[Dict[str, Any]]] = {}
    exons_by_transcript: Dict[str, List[Dict[str, Any]]] = {}
    cds_by_transcript: Dict[str, List[Dict[str, Any]]] = {}
    transcript_to_gene: Dict[str, str] = {}

    with urllib.request.urlopen(req_gff, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                s = line.decode("utf-8", errors="ignore")
                if s.startswith("#"):
                    continue
                parts = s.split("\t")
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
                        biotype = m_bio.group(1) if m_bio else "protein_coding"
                        m_desc = re.search(r"description=([^;]+)", attr)
                        desc = urllib.parse.unquote(m_desc.group(1)) if m_desc else sym

                        if gid in gene_coords:
                            old_chr, old_s, old_e, old_sym = gene_coords[gid]
                            if old_chr == chrom:
                                st = min(st, old_s)
                                en = max(en, old_e)
                        gene_coords[gid] = (chrom, st, en, sym)

                        for tgt in interval_targets:
                            if chrom == tgt["chrom"] and (en >= tgt["filter_start"] and st <= tgt["filter_end"]):
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
                        biotype = m_bio.group(1) if m_bio else ("protein_coding" if feat_type == "mRNA" else "non_coding")

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

    for gid, (chrom, st, en, sym) in gene_coords.items():
        genes_by_chrom.setdefault(chrom, []).append((gid, st, en))

    print(f"[{species}] Parsed GFF3: {len(gene_coords):,} genes across {len(genes_by_chrom):,} contigs in {time.time()-t1:.2f}s")
    print(f"[{species}] Locus target genes found: {len(raw_locus_genes)}")

    # ---------------------------------------------------------
    # Step 3: Stream Genomic FASTA chromosome-by-chromosome
    # ---------------------------------------------------------
    t2 = time.time()
    fna_url = cfg["fna_url"]
    print(f"[{species}] 3/4 Streaming Genomic FASTA: {fna_url.split('/')[-1]}...")
    req_fna = urllib.request.Request(fna_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Globin Project)"})

    gene_whole_gc: Dict[str, Tuple[float, int]] = {}
    current_chrom = None
    seq_parts = []
    total_g = total_c = total_a = total_t = 0

    def process_chrom(chrom_name, parts):
        nonlocal total_g, total_c, total_a, total_t
        if not chrom_name:
            return
        chrom_seq = "".join(parts).upper()
        chr_len = len(chrom_seq)

        total_g += chrom_seq.count("G")
        total_c += chrom_seq.count("C")
        total_a += chrom_seq.count("A")
        total_t += chrom_seq.count("T")

        if chrom_name in genes_by_chrom:
            for gid, s, e in genes_by_chrom[chrom_name]:
                if gid in gene_cds_info:
                    s_idx = max(0, s - 1)
                    e_idx = min(chr_len, e)
                    gseq = chrom_seq[s_idx:e_idx]
                    gc_val, valid_len = calculate_gc(gseq)
                    if valid_len > 0:
                        gene_whole_gc[gid] = (gc_val, valid_len)

        for tgt in interval_targets:
            if chrom_name == tgt["chrom"]:
                loc_s = tgt["orig_start"]
                loc_e = tgt["orig_end"]
                s_idx = max(0, loc_s - 1)
                e_idx = min(chr_len, loc_e)
                loc_seq = chrom_seq[s_idx:e_idx]
                seq_cache_file = os.path.join(CACHE_DIR, f"seq_{species}_{tgt['chrom']}_{loc_s}_{loc_e}.fa")
                with open(seq_cache_file, "w", encoding="utf-8") as f:
                    f.write(loc_seq)
                print(f"[{species}] Cached locus FASTA ({len(loc_seq):,} bp) to {os.path.basename(seq_cache_file)}")

    with urllib.request.urlopen(req_fna, timeout=240) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                s = line.decode("utf-8", errors="ignore")
                if s.startswith(">"):
                    if current_chrom and seq_parts:
                        process_chrom(current_chrom, seq_parts)
                    current_chrom = s[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(s.strip())

            # Final chromosome
            if current_chrom and seq_parts:
                process_chrom(current_chrom, seq_parts)

    genomic_valid = total_g + total_c + total_a + total_t
    empirical_genomic_gc = round(((total_g + total_c) / genomic_valid) * 100.0, 2) if genomic_valid > 0 else 42.0
    print(f"[{species}] Empirical Whole-Genome GC: {empirical_genomic_gc}% (from {genomic_valid:,} bp)")
    print(f"[{species}] Computed true whole-gene GC for {len(gene_whole_gc):,} protein-coding genes in {time.time()-t2:.2f}s")

    # ---------------------------------------------------------
    # Step 4: Cache Locus Gene Details and Region Genes JSON
    # ---------------------------------------------------------
    cached_genes_by_interval: Dict[str, List[Dict[str, Any]]] = {}

    for gid, g_rec in raw_locus_genes.items():
        t_models = []
        for t in transcripts_by_gene.get(gid, []):
            tid = t["id"]
            t_exons = exons_by_transcript.get(tid, [])
            t_cds = cds_by_transcript.get(tid, [])
            t_exons_sorted = sorted(t_exons, key=lambda x: x["start"])
            t_cds_sorted = sorted(t_cds, key=lambda x: x["start"])

            t_models.append({
                "id": tid,
                "name": t["name"],
                "biotype": t["biotype"],
                "start": t["start"],
                "end": t["end"],
                "strand": t["strand"],
                "Exon": t_exons_sorted,
                "CDS": t_cds_sorted,
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

        for tgt in interval_targets:
            if g_rec["chrom"] == tgt["chrom"] and (widest_end >= tgt["orig_start"] and widest_start <= tgt["orig_end"]):
                key = f"{species}_{tgt['chrom']}_{tgt['orig_start']}_{tgt['orig_end']}"
                cached_genes_by_interval.setdefault(key, []).append({
                    "id": gid,
                    "external_name": g_rec["name"],
                    "biotype": g_rec["biotype"],
                    "start": widest_start,
                    "end": widest_end,
                    "strand": g_rec["strand"],
                    "description": g_rec["description"],
                })

    for tgt in interval_targets:
        key = f"{species}_{tgt['chrom']}_{tgt['orig_start']}_{tgt['orig_end']}"
        region_cache_file = os.path.join(CACHE_DIR, f"genes_{key}.json")
        region_genes = cached_genes_by_interval.get(key, [])
        region_genes = sorted(region_genes, key=lambda x: x["start"])
        with open(region_cache_file, "w", encoding="utf-8") as f:
            json.dump(region_genes, f, indent=2)
        print(f"[{species}] Cached {len(region_genes)} interval genes to: {os.path.basename(region_cache_file)}")

    # ---------------------------------------------------------
    # Step 5: Merge into Master Baseline DataFrame & JSON Stats
    # STRICT FILTERING: Only genes with valid CDS (protein_coding)
    # ---------------------------------------------------------
    records = []
    whole_vals = []
    cds_vals = []
    gc3_vals = []
    lens_vals = []

    for gid in sorted(gene_cds_info.keys()):
        c_info = gene_cds_info[gid]
        c_len = c_info["cds_length"]
        c_gc = c_info["cds_gc"]
        gc3 = c_info["gc3"]
        sym = c_info.get("gene_symbol") or gid

        if gid in gene_coords:
            chrom, s, e, g_sym = gene_coords[gid]
            if g_sym and g_sym != gid:
                sym = g_sym
        else:
            chrom, s, e = None, None, None

        if gid in gene_whole_gc:
            w_gc, w_len = gene_whole_gc[gid]
        else:
            w_gc = c_gc
            w_len = (e - s + 1) if (s is not None and e is not None) else c_len

        records.append({
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

        if w_gc is not None:
            whole_vals.append(w_gc)
        if c_gc is not None:
            cds_vals.append(c_gc)
        if gc3 is not None:
            gc3_vals.append(gc3)
        if w_len is not None:
            lens_vals.append(w_len)

    df_master = pd.DataFrame(records)

    arr_w = np.array(whole_vals) if whole_vals else np.array([0.0])
    arr_c = np.array(cds_vals) if cds_vals else np.array([0.0])
    arr_3 = np.array(gc3_vals) if gc3_vals else np.array([0.0])
    arr_l = np.array(lens_vals) if lens_vals else np.array([0.0])

    stats = {
        "species": species,
        "assembly": cfg["assembly"],
        "total_genes": len(df_master),
        "genomic_gc": float(empirical_genomic_gc),
        "mean_whole_gene_len": float(np.mean(arr_l)) if len(arr_l) > 0 else 0.0,
        "std_whole_gene_len": float(np.std(arr_l, ddof=1)) if len(arr_l) > 1 else 0.0,
        "median_whole_gene_len": float(np.median(arr_l)) if len(arr_l) > 0 else 0.0,
        "mean_whole_gene_gc": float(np.mean(arr_w)) if len(arr_w) > 0 else 0.0,
        "std_whole_gene_gc": float(np.std(arr_w, ddof=1)) if len(arr_w) > 1 else 0.0,
        "median_whole_gene_gc": float(np.median(arr_w)) if len(arr_w) > 0 else 0.0,
        "mean_cds_gc": float(np.mean(arr_c)) if len(arr_c) > 0 else 0.0,
        "std_cds_gc": float(np.std(arr_c, ddof=1)) if len(arr_c) > 1 else 0.0,
        "median_cds_gc": float(np.median(arr_c)) if len(arr_c) > 0 else 0.0,
        "mean_gc3": float(np.mean(arr_3)) if len(arr_3) > 0 else None,
        "std_gc3": float(np.std(arr_3, ddof=1)) if len(arr_3) > 1 else None,
        "median_gc3": float(np.median(arr_3)) if len(arr_3) > 0 else None,
    }

    df_master.to_csv(csv_path, index=False, compression="gzip")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[{species}] Successfully saved baseline table ({len(df_master):,} clean protein-coding genes) and stats.")
    print(f"  Genomic GC:    {stats['genomic_gc']}%")
    print(f"  Whole-Gene GC: {stats['mean_whole_gene_gc']:.2f}% +/- {stats['std_whole_gene_gc']:.2f}% (median: {stats['median_whole_gene_gc']:.2f}%)")
    print(f"  Coding CDS GC: {stats['mean_cds_gc']:.2f}% +/- {stats['std_cds_gc']:.2f}% (median: {stats['median_cds_gc']:.2f}%)")
    print(f"  Codon GC3:     {stats['mean_gc3']:.2f}% +/- {stats['std_gc3']:.2f}% (median: {stats['median_gc3']:.2f}%)")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Integrate expansion cohort.")
    parser.add_argument("--species", type=str, default=None, help="Process single species")
    parser.add_argument("--force", action="store_true", help="Force recomputation")
    args = parser.parse_args()

    targets = [args.species] if args.species else list(EXPANSION_SPECIES.keys())
    for sp in targets:
        if sp in EXPANSION_SPECIES:
            process_species(sp, EXPANSION_SPECIES[sp], force=args.force)


if __name__ == "__main__":
    main()
