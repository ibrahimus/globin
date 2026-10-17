"""Extract full gene models (transcripts, exons, CDS, introns) for NCBI RefSeq species.

Streams the official NCBI RefSeq GFF3 files for:
- Pogona vitticeps (PviZW2.1)
- Anolis carolinensis (rAnoCar3.1.pri)
- Alligator mississippiensis (rAllMis1)
- Bufo bufo (aBufBuf1.1)

Extracts all genes within the locus intervals specified in config/species_manifest.yaml,
generating fully structured gene details JSON files matching the Ensembl schema
(with Transcript, Exon, and CDS models) and caching full gene coordinates and CDS sequences.
"""

import os
import re
import gzip
import json
import time
import urllib.request
import yaml
from typing import Dict, List, Any, Tuple

CACHE_DIR = "data/cache"
os.makedirs(CACHE_DIR, exist_ok=True)

NCBI_SPECIES_GFF = {
    "alligator_mississippiensis": {
        "assembly": "rAllMis1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_genomic.gff.gz",
    },
    "anolis_carolinensis": {
        "assembly": "rAnoCar3.1.pri",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/594/765/GCF_035594765.1_rAnoCar3.1.pri/GCF_035594765.1_rAnoCar3.1.pri_genomic.gff.gz",
    },
    "pogona_vitticeps": {
        "assembly": "PviZW2.1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/106/095/GCF_051106095.1_PviZW2.1/GCF_051106095.1_PviZW2.1_genomic.gff.gz",
    },
    "bufo_bufo": {
        "assembly": "aBufBuf1.1",
        "gff_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/171/765/GCF_905171765.1_aBufBuf1.1/GCF_905171765.1_aBufBuf1.1_genomic.gff.gz",
    },
}


def load_manifest_intervals() -> Dict[str, List[Dict[str, Any]]]:
    manifest_path = "config/species_manifest.yaml"
    if not os.path.exists(manifest_path):
        manifest_path = "../config/species_manifest.yaml"
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f)

    species_intervals = {}
    for sp, cfg in manifest.get("species", {}).items():
        if sp in NCBI_SPECIES_GFF:
            species_intervals[sp] = cfg.get("loci", [])
    return species_intervals


def extract_species_gff(species: str, intervals: List[Dict[str, Any]]):
    cfg = NCBI_SPECIES_GFF[species]
    gff_url = cfg["gff_url"]
    print(f"\n========================================================")
    print(f"Extracting NCBI Gene Models for: {species} ({cfg['assembly']})")
    print(f"URL: {gff_url.split('/')[-1]}")
    print(f"========================================================")

    # Prepare interval filters with 100 kb padding to avoid clipping flankers
    interval_targets = []
    for loc in intervals:
        chrom = str(loc["chromosome"])
        start = int(loc["start"]) - 100000
        end = int(loc["end"]) + 100000
        interval_targets.append({
            "locus_id": loc["locus_id"],
            "chrom": chrom,
            "orig_start": int(loc["start"]),
            "orig_end": int(loc["end"]),
            "filter_start": max(0, start),
            "filter_end": end,
        })
        print(f"  Target locus {loc['locus_id']}: {chrom}:{loc['start']}-{loc['end']}")

    # Stream GFF
    t0 = time.time()
    req = urllib.request.Request(gff_url, headers={"User-Agent": "Mozilla/5.0 (Python globin-pipeline)"})

    raw_genes = {}        # gene_internal_id -> dict
    gene_by_name = {}     # Name -> gene_internal_id
    gene_by_geneid = {}   # GeneID:123 -> gene_internal_id
    transcripts = {}      # transcript_internal_id -> dict
    exons = []            # list of (parent_transcript, start, end, strand)
    cds_list = []         # list of (parent_transcript, start, end, strand)

    with urllib.request.urlopen(req, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                s = line.decode("utf-8", errors="ignore")
                if s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) < 9:
                    continue

                chrom = parts[0]
                feat = parts[2]
                s_pos = int(parts[3])
                e_pos = int(parts[4])
                strand_str = parts[6]
                strand_int = 1 if strand_str == "+" else (-1 if strand_str == "-" else 1)
                attr_str = parts[8].strip()

                # Check if this feature overlaps any of our interval targets
                in_target = False
                for tgt in interval_targets:
                    if chrom == tgt["chrom"] and not (e_pos < tgt["filter_start"] or s_pos > tgt["filter_end"]):
                        in_target = True
                        break
                if not in_target:
                    continue

                attrs = {}
                for item in attr_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k] = v

                if feat == "gene":
                    g_id = attrs.get("ID", "")
                    g_name = attrs.get("Name", "")
                    dbxref = attrs.get("Dbxref", "")
                    gene_id_tag = None
                    if "GeneID:" in dbxref:
                        gene_id_tag = "GeneID:" + dbxref.split("GeneID:")[1].split(",")[0]
                    elif "GeneID:" in g_id:
                        gene_id_tag = "GeneID:" + g_id.split("GeneID:")[1]
                    else:
                        gene_id_tag = f"GeneID:{g_name}"

                    biotype = attrs.get("gene_biotype", "protein_coding")
                    desc = attrs.get("description", g_name)

                    gene_record = {
                        "internal_id": g_id,
                        "gene_id": gene_id_tag,
                        "name": g_name,
                        "chrom": chrom,
                        "start": s_pos,
                        "end": e_pos,
                        "strand": strand_int,
                        "biotype": biotype,
                        "description": desc,
                    }
                    raw_genes[g_id] = gene_record
                    if g_name:
                        gene_by_name[g_name] = g_id
                    if gene_id_tag:
                        gene_by_geneid[gene_id_tag] = g_id

                elif feat in ("mRNA", "transcript", "lncRNA", "ncRNA"):
                    t_id = attrs.get("ID", "")
                    parent = attrs.get("Parent", "")
                    t_name = attrs.get("Name", t_id)
                    biotype = "protein_coding" if feat == "mRNA" else feat

                    transcripts[t_id] = {
                        "id": t_id,
                        "name": t_name,
                        "parent": parent,
                        "chrom": chrom,
                        "start": s_pos,
                        "end": e_pos,
                        "strand": strand_int,
                        "biotype": biotype,
                    }

                elif feat == "exon":
                    parent = attrs.get("Parent", "")
                    ex_id = attrs.get("ID", f"exon_{chrom}_{s_pos}_{e_pos}")
                    exons.append({
                        "id": ex_id,
                        "parent": parent,
                        "chrom": chrom,
                        "start": s_pos,
                        "end": e_pos,
                        "strand": strand_int,
                    })

                elif feat == "CDS":
                    parent = attrs.get("Parent", "")
                    cds_id = attrs.get("ID", f"cds_{chrom}_{s_pos}_{e_pos}")
                    cds_list.append({
                        "id": cds_id,
                        "parent": parent,
                        "chrom": chrom,
                        "start": s_pos,
                        "end": e_pos,
                        "strand": strand_int,
                    })

    print(f"[{species}] Streamed GFF in {time.time()-t0:.2f}s: found {len(raw_genes)} genes, {len(transcripts)} transcripts, {len(exons)} exons, {len(cds_list)} CDS segments")

    # Link exons and CDS to transcripts
    transcripts_by_gene: Dict[str, list] = {}
    for t_id, t_data in transcripts.items():
        parent_g = t_data["parent"]
        # Match parent gene
        g_internal = None
        if parent_g in raw_genes:
            g_internal = parent_g
        elif parent_g in gene_by_name:
            g_internal = gene_by_name[parent_g]
        else:
            # Try removing 'gene-' prefix
            clean_parent = re.sub(r"^gene-", "", parent_g)
            if clean_parent in raw_genes:
                g_internal = clean_parent
            elif clean_parent in gene_by_name:
                g_internal = gene_by_name[clean_parent]

        if g_internal:
            transcripts_by_gene.setdefault(g_internal, []).append(t_data)

    exons_by_transcript: Dict[str, list] = {}
    for ex in exons:
        exons_by_transcript.setdefault(ex["parent"], []).append(ex)

    cds_by_transcript: Dict[str, list] = {}
    for cds in cds_list:
        cds_by_transcript.setdefault(cds["parent"], []).append(cds)

    # Build and cache full Ensembl-compatible gene models for each gene
    cached_genes_for_interval: Dict[str, list] = {}

    for g_internal, g_rec in raw_genes.items():
        gid = g_rec["gene_id"]
        gene_transcripts = transcripts_by_gene.get(g_internal, [])

        t_models = []
        for t in gene_transcripts:
            tid = t["id"]
            t_exons = exons_by_transcript.get(tid, [])
            t_cds = cds_by_transcript.get(tid, [])

            # Sort exons and CDS by start position
            t_exons_sorted = sorted(t_exons, key=lambda x: x["start"])
            t_cds_sorted = sorted(t_cds, key=lambda x: x["start"])

            exon_records = [{
                "id": ex["id"],
                "start": ex["start"],
                "end": ex["end"],
                "strand": ex["strand"],
            } for ex in t_exons_sorted]

            cds_records = [{
                "id": c["id"],
                "start": c["start"],
                "end": c["end"],
                "strand": c["strand"],
            } for c in t_cds_sorted]

            t_models.append({
                "id": tid,
                "name": t["name"],
                "biotype": t["biotype"],
                "start": t["start"],
                "end": t["end"],
                "strand": t["strand"],
                "Exon": exon_records,
                "CDS": cds_records,
            })

        # If gene had no explicit transcripts in GFF, create dummy transcript spanning gene
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

        # Compute widest coordinates from transcripts/exons if wider than gene feature
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

        # Save to data/cache/gene_{gid}.json
        safe_gid = "".join([c if c.isalnum() or c in "._-" else "_" for c in gid])
        detail_path = os.path.join(CACHE_DIR, f"gene_{safe_gid}.json")
        with open(detail_path, "w", encoding="utf-8") as f:
            json.dump(gene_detail_obj, f, indent=2)

        # Check which target intervals this gene belongs to
        for tgt in interval_targets:
            if g_rec["chrom"] == tgt["chrom"] and (widest_end >= tgt["orig_start"] and widest_start <= tgt["orig_end"]):
                key = f"{species}_{tgt['chrom']}_{tgt['orig_start']}_{tgt['orig_end']}"
                cached_genes_for_interval.setdefault(key, []).append({
                    "id": gid,
                    "external_name": g_rec["name"],
                    "biotype": g_rec["biotype"],
                    "start": widest_start,
                    "end": widest_end,
                    "strand": g_rec["strand"],
                    "description": g_rec["description"],
                })

    # Cache region gene list files
    for tgt in interval_targets:
        key = f"{species}_{tgt['chrom']}_{tgt['orig_start']}_{tgt['orig_end']}"
        region_cache_file = os.path.join(CACHE_DIR, f"genes_{key}.json")
        region_genes = cached_genes_for_interval.get(key, [])
        # Sort by start coordinate
        region_genes = sorted(region_genes, key=lambda x: x["start"])
        with open(region_cache_file, "w", encoding="utf-8") as f:
            json.dump(region_genes, f, indent=2)
        print(f"[{species}] Cached {len(region_genes)} interval genes to: {os.path.basename(region_cache_file)}")


def main():
    intervals_by_sp = load_manifest_intervals()
    for sp in ["pogona_vitticeps", "anolis_carolinensis", "alligator_mississippiensis", "bufo_bufo"]:
        if sp in intervals_by_sp:
            extract_species_gff(sp, intervals_by_sp[sp])


if __name__ == "__main__":
    main()
