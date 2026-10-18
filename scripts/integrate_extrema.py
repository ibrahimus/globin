#!/usr/bin/env python3
"""
Integrate 10 Physiological Extrema & Evolutionary Transition Genomes:
1. Latimeria chalumnae (West Indian Ocean Coelacanth) - Sarcopterygian living fossil bridge to tetrapods
2. Polypterus senegalus (Senegal Bichir) - Basal ray-finned fish (Cladistia) with functional lungs
3. Camelus dromedarius (Dromedary Camel) - Desert heat & dehydration specialist (oval erythrocytes)
4. Vicugna pacos (Alpaca) - Andean high-altitude hypobaric hypoxia champion (>4,500 m)
5. Desmodus rotundus (Common Vampire Bat) - Sole mammalian obligate blood feeder (hematophagy)
6. Mirounga angustirostris (Northern Elephant Seal) - Mammalian champion deep diver (>1,700 m, 120 min)
7. Physeter catodon (Sperm Whale) - Champion oceanic deep diver (>2,000 m)
8. Tachyglossus aculeatus (Short-beaked Echidna) - Egg-laying monotreme, subterranean burrower
9. Calypte anna (Anna's Hummingbird) - Highest vertebrate mass-specific metabolic rate & torpor
10. Aptenodytes patagonicus (King Penguin) - Extreme Antarctic diving bird (>300 m in sub-zero water)

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

EXTREMA_SPECIES = {
    "latimeria_chalumnae": {
        "common_name": "West Indian Ocean Coelacanth",
        "assembly": "fLatCha1.pri",
        "clade": "sarcopterygii",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_beta_cluster",
                "chrom": "NC_088143.1",
                "start": 166200000,
                "end": 166550000,
                "strand": "+",
                "left_gene": "LOC102356551",
                "right_gene": "LOC102354658",
            },
        ],
    },
    "polypterus_senegalus": {
        "common_name": "Senegal Bichir",
        "assembly": "ASM1683550v1",
        "clade": "non_teleost_fish",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/835/505/GCF_016835505.1_ASM1683550v1/GCF_016835505.1_ASM1683550v1_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/835/505/GCF_016835505.1_ASM1683550v1/GCF_016835505.1_ASM1683550v1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/835/505/GCF_016835505.1_ASM1683550v1/GCF_016835505.1_ASM1683550v1_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_beta_cluster",
                "chrom": "NC_053166.1",
                "start": 157250000,
                "end": 157900000,
                "strand": "+",
                "left_gene": "LOC120542922",
                "right_gene": "LOC120542357",
            },
        ],
    },
    "camelus_dromedarius": {
        "common_name": "Dromedary Camel",
        "assembly": "mCamDro1.pat",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/321/535/GCF_036321535.1_mCamDro1.pat/GCF_036321535.1_mCamDro1.pat_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/321/535/GCF_036321535.1_mCamDro1.pat/GCF_036321535.1_mCamDro1.pat_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/321/535/GCF_036321535.1_mCamDro1.pat/GCF_036321535.1_mCamDro1.pat_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_087459.1",
                "start": 6850000,
                "end": 7030000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "FAM234A",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_087447.1",
                "start": 45400000,
                "end": 45450000,
                "strand": "-",
                "left_gene": "HBE1",
                "right_gene": "LOC105097953",
            },
        ],
    },
    "vicugna_pacos": {
        "common_name": "Alpaca",
        "assembly": "VicPac4",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/564/905/GCF_048564905.1_VicPac4/GCF_048564905.1_VicPac4_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/564/905/GCF_048564905.1_VicPac4/GCF_048564905.1_VicPac4_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/564/905/GCF_048564905.1_VicPac4/GCF_048564905.1_VicPac4_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_133004.1",
                "start": 18680000,
                "end": 18880000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "FAM234A",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_132996.1",
                "start": 29710000,
                "end": 29750000,
                "strand": "-",
                "left_gene": "LOC102525863",
                "right_gene": "HBE1",
            },
        ],
    },
    "desmodus_rotundus": {
        "common_name": "Common Vampire Bat",
        "assembly": "HLdesRot8A.1",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/682/495/GCF_022682495.2_HLdesRot8A.1/GCF_022682495.2_HLdesRot8A.1_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/682/495/GCF_022682495.2_HLdesRot8A.1/GCF_022682495.2_HLdesRot8A.1_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/682/495/GCF_022682495.2_HLdesRot8A.1/GCF_022682495.2_HLdesRot8A.1_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_071387.1",
                "start": 97650000,
                "end": 97850000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "FAM234A",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_071391.1",
                "start": 44400000,
                "end": 44500000,
                "strand": "-",
                "left_gene": "LOC112317580",
                "right_gene": "LOC128781113",
            },
        ],
    },
    "mirounga_angustirostris": {
        "common_name": "Northern Elephant Seal",
        "assembly": "ASM2128878v3",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/288/785/GCF_021288785.2_ASM2128878v3/GCF_021288785.2_ASM2128878v3_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/288/785/GCF_021288785.2_ASM2128878v3/GCF_021288785.2_ASM2128878v3_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/288/785/GCF_021288785.2_ASM2128878v3/GCF_021288785.2_ASM2128878v3_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_072358.1",
                "start": 128750000,
                "end": 128950000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "FAM234A",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_072366.1",
                "start": 47000000,
                "end": 48450000,
                "strand": "+",
                "left_gene": "HBE1",
                "right_gene": "CLPB",
            },
        ],
    },
    "physeter_catodon": {
        "common_name": "Sperm Whale",
        "assembly": "ASM283717v5",
        "clade": "amniote_mammal",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/837/175/GCF_002837175.3_ASM283717v5/GCF_002837175.3_ASM283717v5_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/837/175/GCF_002837175.3_ASM283717v5/GCF_002837175.3_ASM283717v5_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/837/175/GCF_002837175.3_ASM283717v5/GCF_002837175.3_ASM283717v5_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NW_021146805.1",
                "start": 1,
                "end": 25000,
                "strand": "-",
                "left_gene": "HBZ",
                "right_gene": "HBQ1",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_041229.1",
                "start": 60300000,
                "end": 61400000,
                "strand": "+",
                "left_gene": "LOC102992241",
                "right_gene": "UBQLN3",
            },
        ],
    },
    "tachyglossus_aculeatus": {
        "common_name": "Short-beaked Echidna",
        "assembly": "mTacAcu1.pri",
        "clade": "monotreme",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/852/505/GCF_015852505.1_mTacAcu1.pri/GCF_015852505.1_mTacAcu1.pri_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/852/505/GCF_015852505.1_mTacAcu1.pri/GCF_015852505.1_mTacAcu1.pri_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/852/505/GCF_015852505.1_mTacAcu1.pri/GCF_015852505.1_mTacAcu1.pri_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_052086.1",
                "start": 10520000,
                "end": 10730000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "FAM234A",
            },
            {
                "locus_id": "beta_globin",
                "chrom": "NC_052067.1",
                "start": 137800000,
                "end": 138000000,
                "strand": "+",
                "left_gene": "LOC119943796",
                "right_gene": "RRM1",
            },
        ],
    },
    "calypte_anna": {
        "common_name": "Anna's Hummingbird",
        "assembly": "bCalAnn1_v1.p",
        "clade": "amniote_bird",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_044260.1",
                "start": 14200000,
                "end": 14500000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "LUC7L",
            },
        ],
    },
    "aptenodytes_patagonicus": {
        "common_name": "King Penguin",
        "assembly": "bAptPat1.pri.cur",
        "clade": "amniote_bird",
        "cds_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/965/638/725/GCF_965638725.1_bAptPat1.pri.cur/GCF_965638725.1_bAptPat1.pri.cur_cds_from_genomic.fna.gz",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/965/638/725/GCF_965638725.1_bAptPat1.pri.cur/GCF_965638725.1_bAptPat1.pri.cur_genomic.gff.gz",
        "fna_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/965/638/725/GCF_965638725.1_bAptPat1.pri.cur/GCF_965638725.1_bAptPat1.pri.cur_genomic.fna.gz",
        "loci": [
            {
                "locus_id": "alpha_globin",
                "chrom": "NC_134961.1",
                "start": 140000,
                "end": 480000,
                "strand": "+",
                "left_gene": "RHBDF1",
                "right_gene": "LUC7L",
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
    print(f"Processing Genome: {species} ({cfg['assembly']})")
    print(f"========================================================")

    # ---------------------------------------------------------
    # Step 1: Stream CDS FASTA for canonical CDS GC and GC3
    # ---------------------------------------------------------
    t0 = time.time()
    cds_url = cfg["cds_url"]
    print(f"[{species}] 1/4 Streaming CDS FASTA: {cds_url.split('/')[-1]}...")
    req_cds = urllib.request.Request(cds_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Globin Project)"})

    gene_cds_info: Dict[str, Dict[str, Any]] = {}

    with urllib.request.urlopen(req_cds, timeout=180) as resp:
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

    with urllib.request.urlopen(req_gff, timeout=180) as resp:
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

    with urllib.request.urlopen(req_fna, timeout=300) as resp:
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
            "whole_gene_length": w_len,
            "whole_gene_gc": w_gc,
            "cds_length": c_len,
            "cds_gc": c_gc,
            "gc3": gc3,
        })
        whole_vals.append(w_gc)
        cds_vals.append(c_gc)
        if gc3 is not None:
            gc3_vals.append(gc3)
        lens_vals.append(w_len)

    df = pd.DataFrame(records)
    df.to_csv(csv_path, index=False, compression="gzip")
    print(f"[{species}] Saved {len(df):,} baseline genes to {csv_path}")

    stats_obj = {
        "species": species,
        "common_name": cfg["common_name"],
        "assembly": cfg["assembly"],
        "clade": cfg["clade"],
        "source": "NCBI_RefSeq",
        "empirical_genomic_gc": empirical_genomic_gc,
        "n_genes": len(df),
        "mean_gene_gc": round(float(np.mean(whole_vals)), 2) if whole_vals else 0.0,
        "std_gene_gc": round(float(np.std(whole_vals)), 2) if whole_vals else 0.0,
        "median_gene_gc": round(float(np.median(whole_vals)), 2) if whole_vals else 0.0,
        "mean_cds_gc": round(float(np.mean(cds_vals)), 2) if cds_vals else 0.0,
        "std_cds_gc": round(float(np.std(cds_vals)), 2) if cds_vals else 0.0,
        "median_cds_gc": round(float(np.median(cds_vals)), 2) if cds_vals else 0.0,
        "mean_gc3": round(float(np.mean(gc3_vals)), 2) if gc3_vals else 0.0,
        "std_gc3": round(float(np.std(gc3_vals)), 2) if gc3_vals else 0.0,
        "median_gc3": round(float(np.median(gc3_vals)), 2) if gc3_vals else 0.0,
        "mean_gene_length": round(float(np.mean(lens_vals)), 1) if lens_vals else 0.0,
        "median_gene_length": round(float(np.median(lens_vals)), 1) if lens_vals else 0.0,
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(stats_obj, f, indent=2)
    print(f"[{species}] Saved summary stats to {json_path}")


def main():
    import sys
    force = "--force" in sys.argv
    for sp, cfg in EXTREMA_SPECIES.items():
        process_species(sp, cfg, force=force)
    print("\nAll 10 physiological extrema & evolutionary transition species successfully processed!")


if __name__ == "__main__":
    main()
