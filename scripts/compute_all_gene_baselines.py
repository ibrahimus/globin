"""Script to calculate empirical Whole-Gene GC, CDS GC, and GC3 for all protein-coding genes across vertebrate genomes."""

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
os.makedirs(CACHE_DIR, exist_ok=True)

SPECIES_ASSEMBLIES = {
    # Active Vertebrates (24 species)
    "Gallus_gallus": "bGalGal1.mat.broiler.GRCg7b",
    "Homo_sapiens": "GRCh38",
    "Mus_musculus": "GRCm39",
    "Bos_taurus": "ARS-UCD2.0",
    "Sus_scrofa": "Sscrofa11.1",
    "Monodelphis_domestica": "ASM229v1",
    "Chrysemys_picta_bellii": "Chrysemys_picta_bellii-3.0.3",
    "pelodiscus_sinensis": "PelSin_1.0",
    "gopherus_evgoodei": "rGopEvg1_v1.p",
    "terrapene_carolina_triunguis": "T_m_triunguis-2.0",
    "crocodylus_porosus": "CroPor_comp1",
    "taeniopygia_guttata": "bTaeGut1_v1.p",
    "anas_platyrhynchos": "ASM874695v1",
    "meleagris_gallopavo": "Turkey_5.1",
    "struthio_camelus_australis": "ASM69896v1",
    "anolis_carolinensis": "AnoCar2.0v2",
    "pogona_vitticeps": "pvi1.1",
    "xenopus_tropicalis": "UCB_Xtro_10.0",
    "danio_rerio": "GRCz11",
    "oryzias_latipes": "ASM223467v1",
    "takifugu_rubripes": "fTakRub1.2",
    "lepisosteus_oculatus": "LepOcu1",
    "petromyzon_marinus": "Pmarinus_7.0",
    "ornithorhynchus_anatinus": "mOrnAna1.p.v1",
    "aquila_chrysaetos_chrysaetos": "bAquChr1.2",
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


def resolve_urls(species: str, assembly: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    sp_lower = species.lower()
    sp_cap = species[0].upper() + species[1:]

    # 1. CDS URL
    cds_url = f"https://ftp.ensembl.org/pub/current_fasta/{sp_lower}/cds/{sp_cap}.{assembly}.cds.all.fa.gz"

    # 2. FASTA DNA URL
    fasta_dir = f"https://ftp.ensembl.org/pub/current_fasta/{sp_lower}/dna/"
    dna_url = None
    try:
        req = urllib.request.Request(fasta_dir, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=12) as resp:
            html = resp.read().decode("utf-8")
        primaries = re.findall(r'href=[\"\'](' + re.escape(sp_cap) + r'\.[^\"\']*\.dna\.primary_assembly\.fa\.gz)[\"\']', html)
        if not primaries:
            primaries = re.findall(r'href=[\"\'](' + re.escape(sp_cap) + r'\.[^\"\']*\.dna\.toplevel\.fa\.gz)[\"\']', html)
        if primaries:
            dna_url = fasta_dir + primaries[0]
    except Exception as e:
        print(f"Error accessing FASTA dir {fasta_dir}: {e}")

    # 3. GTF URL across releases 114, 113, 112, 111, 110
    gtf_url = None
    for rel in [114, 113, 112, 111, 110]:
        gtf_dir = f"https://ftp.ensembl.org/pub/release-{rel}/gtf/{sp_lower}/"
        try:
            req_gtf = urllib.request.Request(gtf_dir, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req_gtf, timeout=10) as resp:
                gtf_html = resp.read().decode("utf-8")
            gtf_files = re.findall(r'href=[\"\'](' + re.escape(sp_cap) + r'\.[^\"\']*\.' + str(rel) + r'\.gtf\.gz)[\"\']', gtf_html)
            if gtf_files:
                gtf_url = gtf_dir + gtf_files[0]
                break
        except:
            continue

    return cds_url, dna_url, gtf_url


def process_species(species: str, assembly: str, force: bool = False):
    sp_clean = species.lower()
    table_path = os.path.join(CACHE_DIR, f"{sp_clean}_genes_gc.csv.gz")
    stats_path = os.path.join(CACHE_DIR, f"{sp_clean}_stats.json")

    # Check if already processed with whole_gene_gc
    if not force and os.path.exists(table_path) and os.path.exists(stats_path):
        try:
            df_existing = pd.read_csv(table_path)
            if "whole_gene_gc" in df_existing.columns and df_existing["whole_gene_gc"].notna().sum() > 1000:
                print(f"[{species}] Already cached with whole_gene_gc ({len(df_existing):,} genes). Skipping.")
                return
        except Exception:
            pass

    print(f"\n==================================================")
    print(f"Processing: {species} ({assembly})")
    print(f"==================================================")

    cds_url, dna_url, gtf_url = resolve_urls(species, assembly)
    if not dna_url or not gtf_url:
        print(f"Error: Missing DNA or GTF URL for {species}. DNA={dna_url}, GTF={gtf_url}")
        return

    # Step 1: Stream GTF to get protein-coding gene genomic coordinates
    print(f"[{species}] Streaming GTF from: {gtf_url.split('/')[-1]}...")
    t0 = time.time()
    req_gtf = urllib.request.Request(gtf_url, headers={"User-Agent": "Mozilla/5.0"})
    genes_by_chrom = {}
    gene_meta = {}

    with urllib.request.urlopen(req_gtf, timeout=45) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for line in gz:
                line_str = line.decode("utf-8", errors="ignore")
                if line_str.startswith("#"):
                    continue
                parts = line_str.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "gene":
                    attr = parts[8]
                    if 'gene_biotype "protein_coding"' in attr:
                        gid = attr.split('gene_id "')[1].split('"')[0]
                        sym = gid
                        if 'gene_name "' in attr:
                            sym = attr.split('gene_name "')[1].split('"')[0]
                        chrom = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        genes_by_chrom.setdefault(chrom, []).append((gid, start, end))
                        gene_meta[gid] = {
                            "gene_id": gid,
                            "gene_symbol": sym,
                            "chromosome": chrom,
                            "start": start,
                            "end": end,
                        }

    print(f"[{species}] Extracted {len(gene_meta):,} protein-coding genes in {time.time()-t0:.2f}s")

    # Step 2: Stream FASTA DNA and compute whole-gene GC for each gene
    print(f"[{species}] Streaming assembly DNA from: {dna_url.split('/')[-1]}...")
    t1 = time.time()
    req_dna = urllib.request.Request(dna_url, headers={"User-Agent": "Mozilla/5.0"})
    gene_whole_gc = {}

    with urllib.request.urlopen(req_dna, timeout=90) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            current_chrom = None
            seq_parts = []
            for line in gz:
                line_str = line.decode("utf-8", errors="ignore")
                if line_str.startswith(">"):
                    if current_chrom and current_chrom in genes_by_chrom:
                        chrom_seq = "".join(seq_parts).upper()
                        chr_len = len(chrom_seq)
                        for gid, s, e in genes_by_chrom[current_chrom]:
                            s_idx = max(0, s - 1)
                            e_idx = min(chr_len, e)
                            gseq = chrom_seq[s_idx:e_idx]
                            gc_val, valid_len = calculate_gc(gseq)
                            if valid_len > 0:
                                gene_whole_gc[gid] = (gc_val, valid_len)
                    current_chrom = line_str[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line_str.strip())

            # Final chromosome
            if current_chrom and current_chrom in genes_by_chrom:
                chrom_seq = "".join(seq_parts).upper()
                chr_len = len(chrom_seq)
                for gid, s, e in genes_by_chrom[current_chrom]:
                    s_idx = max(0, s - 1)
                    e_idx = min(chr_len, e)
                    gseq = chrom_seq[s_idx:e_idx]
                    gc_val, valid_len = calculate_gc(gseq)
                    if valid_len > 0:
                        gene_whole_gc[gid] = (gc_val, valid_len)

    print(f"[{species}] Computed whole-gene GC for {len(gene_whole_gc):,} genes in {time.time()-t1:.2f}s")

    # Step 3: Stream CDS FASTA (or reuse existing cached CDS table)
    genes_cds = {}
    if os.path.exists(table_path):
        try:
            df_old = pd.read_csv(table_path)
            for _, r in df_old.iterrows():
                genes_cds[r["gene_id"]] = {
                    "cds_gc": r.get("cds_gc"),
                    "cds_len": r.get("cds_length", 0),
                    "gc3": r.get("gc3"),
                    "gene_symbol": r.get("gene_symbol"),
                }
        except Exception:
            pass

    if not genes_cds:
        print(f"[{species}] Streaming CDS FASTA from: {cds_url.split('/')[-1]}...")
        req_cds = urllib.request.Request(cds_url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req_cds, timeout=45) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                current_gid = None
                current_sym = None
                current_seq = []
                for line in gz:
                    line_str = line.decode("utf-8", errors="ignore")
                    if line_str.startswith(">"):
                        if current_gid and current_seq:
                            s_str = "".join(current_seq)
                            if current_gid not in genes_cds or len(s_str) > genes_cds[current_gid]["cds_len"]:
                                gc_val, valid_len = calculate_gc(s_str)
                                gc3_val = calculate_gc3(s_str)
                                genes_cds[current_gid] = {
                                    "cds_gc": gc_val,
                                    "cds_len": valid_len,
                                    "gc3": gc3_val,
                                    "gene_symbol": current_sym or current_gid,
                                }
                        current_gid = None
                        current_sym = None
                        current_seq = []
                        for p in line_str.strip().split():
                            if p.startswith("gene:"):
                                current_gid = p.split(":")[1].split(".")[0]
                            elif p.startswith("gene_symbol:"):
                                current_sym = p.split(":")[1]
                    else:
                        current_seq.append(line_str.strip())
                if current_gid and current_seq:
                    s_str = "".join(current_seq)
                    if current_gid not in genes_cds or len(s_str) > genes_cds[current_gid]["cds_len"]:
                        gc_val, valid_len = calculate_gc(s_str)
                        gc3_val = calculate_gc3(s_str)
                        genes_cds[current_gid] = {
                            "cds_gc": gc_val,
                            "cds_len": valid_len,
                            "gc3": gc3_val,
                            "gene_symbol": current_sym or current_gid,
                        }

    # Step 4: Assemble comprehensive all-genes DataFrame
    records = []
    whole_vals = []
    cds_vals = []
    gc3_vals = []

    for gid, meta in gene_meta.items():
        w_info = gene_whole_gc.get(gid)
        w_gc = w_info[0] if w_info else None
        w_len = w_info[1] if w_info else None

        c_info = genes_cds.get(gid, {})
        c_gc = c_info.get("cds_gc")
        c_len = c_info.get("cds_len", 0)
        gc3 = c_info.get("gc3")
        sym = c_info.get("gene_symbol") or meta["gene_symbol"]

        records.append({
            "species": species,
            "gene_id": gid,
            "gene_symbol": sym,
            "biotype": "protein_coding",
            "chromosome": meta["chromosome"],
            "start": meta["start"],
            "end": meta["end"],
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

    df = pd.DataFrame(records)

    arr_w = np.array(whole_vals) if whole_vals else np.array([0.0])
    arr_c = np.array(cds_vals) if cds_vals else np.array([0.0])
    arr_3 = np.array(gc3_vals) if gc3_vals else np.array([0.0])

    stats = {
        "species": species,
        "assembly": assembly,
        "num_genes": len(df),
        # Whole-Gene Body Metrics (ALL ~20,000 protein-coding genes)
        "num_whole_genes": len(whole_vals),
        "mean_whole_gene_gc": round(float(np.mean(arr_w)), 2),
        "std_whole_gene_gc": round(float(np.std(arr_w)), 2),
        "median_whole_gene_gc": round(float(np.median(arr_w)), 2),
        "q25_whole_gene_gc": round(float(np.percentile(arr_w, 25)), 2),
        "q75_whole_gene_gc": round(float(np.percentile(arr_w, 75)), 2),
        # Coding Sequence Metrics (ALL ~20,000 protein-coding genes)
        "num_cds_genes": len(cds_vals),
        "mean_cds_gc": round(float(np.mean(arr_c)), 2),
        "std_cds_gc": round(float(np.std(arr_c)), 2),
        "median_cds_gc": round(float(np.median(arr_c)), 2),
        "q25_cds_gc": round(float(np.percentile(arr_c, 25)), 2),
        "q75_cds_gc": round(float(np.percentile(arr_c, 75)), 2),
        # Codon GC3 Metrics (ALL ~20,000 protein-coding genes)
        "num_gc3_genes": len(gc3_vals),
        "mean_gc3": round(float(np.mean(arr_3)), 2),
        "std_gc3": round(float(np.std(arr_3)), 2),
        "median_gc3": round(float(np.median(arr_3)), 2),
        "q25_gc3": round(float(np.percentile(arr_3, 25)), 2),
        "q75_gc3": round(float(np.percentile(arr_3, 75)), 2),
    }

    # Save to disk
    df.to_csv(table_path, index=False, compression="gzip")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[{species}] Saved empirical baseline table ({len(df):,} genes) and stats.")
    print(f"  Whole-Gene GC: {stats['mean_whole_gene_gc']}% +/- {stats['std_whole_gene_gc']}% (median: {stats['median_whole_gene_gc']}%)")
    print(f"  Coding CDS GC: {stats['mean_cds_gc']}% +/- {stats['std_cds_gc']}% (median: {stats['median_cds_gc']}%)")
    print(f"  Codon GC3:     {stats['mean_gc3']}% +/- {stats['std_gc3']}% (median: {stats['median_gc3']}%)")


def main():
    import argparse
    from concurrent.futures import ProcessPoolExecutor, as_completed

    parser = argparse.ArgumentParser(description="Compute empirical whole-gene baselines across species.")
    parser.add_argument("--workers", type=int, default=3, help="Number of parallel species to process (default: 3)")
    parser.add_argument("--force", action="store_true", help="Force recomputation even if cached")
    parser.add_argument("--species", type=str, default=None, help="Process single species (optional)")
    args = parser.parse_args()

    if args.species:
        if args.species in SPECIES_ASSEMBLIES:
            process_species(args.species, SPECIES_ASSEMBLIES[args.species], force=args.force)
        else:
            print(f"Species {args.species} not found in SPECIES_ASSEMBLIES.")
        return

    print(f"Starting Empirical Baseline Computation for {len(SPECIES_ASSEMBLIES)} Vertebrate Species with {args.workers} workers...")
    
    # Filter out species already cached unless force is set
    to_process = []
    for sp, asm in SPECIES_ASSEMBLIES.items():
        sp_clean = sp.lower()
        table_path = os.path.join(CACHE_DIR, f"{sp_clean}_genes_gc.csv.gz")
        stats_path = os.path.join(CACHE_DIR, f"{sp_clean}_stats.json")
        already_done = False
        if not args.force and os.path.exists(table_path) and os.path.exists(stats_path):
            try:
                with open(stats_path, "r", encoding="utf-8") as f:
                    st = json.load(f)
                if st.get("mean_whole_gene_gc") is not None and st.get("num_whole_genes", 0) > 1000:
                    already_done = True
            except Exception:
                pass
        if already_done:
            print(f"[{sp}] Already cached with whole_gene_gc ({st.get('num_whole_genes', 0):,} genes). Skipping.")
        else:
            to_process.append((sp, asm))

    print(f"\nSpecies to compute: {len(to_process)} out of {len(SPECIES_ASSEMBLIES)}")
    if not to_process:
        print("All species are already computed and cached!")
        return

    start_time = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_sp = {
            executor.submit(process_species, sp, asm, args.force): sp
            for sp, asm in to_process
        }
        for future in as_completed(future_to_sp):
            sp = future_to_sp[future]
            try:
                future.result()
            except Exception as e:
                print(f"[{sp}] Failed with error: {e}")

    elapsed = time.time() - start_time
    print(f"\nAll species baseline processing completed in {elapsed/60:.2f} minutes!")


if __name__ == "__main__":
    main()

