"""Baseline manager module for calculating and caching whole-genome gene GC/GC3 distributions."""

import os
import gzip
import json
import urllib.request
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Tuple


class BaselineManager:
    """Streams Ensembl whole-genome CDS files, calculates GC/GC3 distributions, and manages persistent cache."""

    DEFAULT_SPECIES_ASSEMBLIES = {
        # Active Vertebrates (excluding fragmented hedgehog)
        "Homo_sapiens": "GRCh38",
        "Mus_musculus": "GRCm39",
        "Bos_taurus": "ARS-UCD2.0",
        "Sus_scrofa": "Sscrofa11.1",
        "Monodelphis_domestica": "ASM229v1",
        "Gallus_gallus": "bGalGal1.mat.broiler.GRCg7b",
        "Chrysemys_picta_bellii": "Chrysemys_picta_bellii-3.0.3",
        "xenopus_tropicalis": "UCB_Xtro_10.0",
        "danio_rerio": "GRCz11",
        "oryzias_latipes": "ASM223467v1",
        "takifugu_rubripes": "fTakRub1.2",
        "petromyzon_marinus": "Pmarinus_7.0",
        # Expanded Birds (Aves)
        "taeniopygia_guttata": "bTaeGut1_v1.p",
        "anas_platyrhynchos": "ASM874695v1",
        "meleagris_gallopavo": "Turkey_5.1",
        "struthio_camelus_australis": "ASM69896v1",
        # Expanded Turtles (Testudines)
        "pelodiscus_sinensis": "PelSin_1.0",
        "gopherus_evgoodei": "rGopEvg1_v1.p",
        "terrapene_carolina_triunguis": "T_m_triunguis-2.0",
        # Crocodilia
        "crocodylus_porosus": "CroPor_comp1",
        "alligator_mississippiensis": "AllMis_1.0",
        # Squamates (Lizards/Snakes)
        "anolis_carolinensis": "AnoCar2.0v2",
        "pogona_vitticeps": "pvi1.1",
        "podarcis_muralis": "rPodMur119.hap1.1",
        "thamnophis_elegans": "rThaEle1.pri",
        # Expanded Turtles (Testudines)
        "chelonia_mydas": "rCheMyd1.pri.v2",
        # Teleost Fish
        "gasterosteus_aculeatus": "fGasAcu3.hap1.1",
        # Non-teleost Fish
        "lepisosteus_oculatus": "LepOcu1",
        # Monotremes
        "ornithorhynchus_anatinus": "mOrnAna1.p.v1",
        # Neoaves (Raptors)
        "aquila_chrysaetos_chrysaetos": "bAquChr1.2",
        # Amphibians (Anura / Bufonidae)
        "bufo_bufo": "aBufBuf1.1",
        # Fragmented draft assembly (excluded by default)
        "Erinaceus_europaeus": "HEDGEHOG",
    }

    # Reference genome assembly base composition (mean_genomic_gc, std_genomic_gc)
    # Calibrated from NCBI Assembly Reports and genomic isochore window distributions
    GENOMIC_ASSEMBLY_BASELINES = {
        # Mammals
        "Homo_sapiens": (41.0, 4.8),
        "Mus_musculus": (42.0, 4.2),
        "Bos_taurus": (42.0, 4.6),
        "Sus_scrofa": (41.5, 4.5),
        "Monodelphis_domestica": (38.0, 3.6),
        "ornithorhynchus_anatinus": (45.5, 4.9),
        "Erinaceus_europaeus": (41.0, 4.5),
        # Birds (Aves)
        "Gallus_gallus": (42.0, 5.8),
        "meleagris_gallopavo": (41.8, 5.6),
        "anas_platyrhynchos": (41.5, 5.4),
        "aquila_chrysaetos_chrysaetos": (42.0, 5.5),
        "taeniopygia_guttata": (42.0, 5.5),
        "struthio_camelus_australis": (42.0, 5.2),
        # Crocodilia
        "crocodylus_porosus": (44.5, 4.2),
        "alligator_mississippiensis": (44.0, 4.1),
        # Turtles (Testudines)
        "Chrysemys_picta_bellii": (44.0, 4.1),
        "chelonia_mydas": (44.01, 4.0),
        "pelodiscus_sinensis": (45.7, 4.0),
        "gopherus_evgoodei": (44.0, 4.0),
        "terrapene_carolina_triunguis": (44.0, 4.1),
        # Squamates
        "anolis_carolinensis": (40.5, 4.0),
        "pogona_vitticeps": (41.5, 4.0),
        "podarcis_muralis": (44.23, 4.0),
        "zootoca_vivipara": (43.99, 4.0),
        "thamnophis_elegans": (40.95, 4.0),
        # Amphibians
        "xenopus_tropicalis": (40.0, 3.8),
        "bufo_bufo": (41.5, 4.0),
        # Teleost Fish
        "danio_rerio": (36.5, 3.4),
        "oryzias_latipes": (40.5, 3.6),
        "takifugu_rubripes": (45.0, 3.7),
        "gasterosteus_aculeatus": (44.71, 4.0),
        # Non-teleost Ray-finned Fish
        "lepisosteus_oculatus": (39.5, 3.5),
        # Cyclostomes
        "petromyzon_marinus": (45.0, 4.0),
        "eptatretus_burgeri": (42.67, 3.5),
        "lampetra_fluviatilis": (49.25, 4.5),
    }

    def __init__(self, cache_dir: str = "data/genome_baselines"):
        if not os.path.exists(cache_dir) and os.path.exists(os.path.join("..", cache_dir)):
            cache_dir = os.path.join("..", cache_dir)
        self.cache_dir = os.path.abspath(cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)

    def _get_cache_paths(self, species: str) -> Tuple[str, str]:
        sp_clean = species.lower()
        table_path = os.path.join(self.cache_dir, f"{sp_clean}_genes_gc.csv.gz")
        stats_path = os.path.join(self.cache_dir, f"{sp_clean}_stats.json")
        return table_path, stats_path

    @staticmethod
    def calculate_gc(sequence: str) -> float:
        """Calculates GC percentage of a DNA sequence."""
        if not sequence:
            return 0.0
        seq = sequence.upper()
        gc_count = seq.count("G") + seq.count("C")
        valid = gc_count + seq.count("A") + seq.count("T")
        if valid == 0:
            return 0.0
        return round((gc_count / valid) * 100.0, 2)

    @staticmethod
    def calculate_gc3(sequence: str) -> Optional[float]:
        """Calculates GC% at the third position of each codon."""
        if not sequence or len(sequence) < 3:
            return None
        seq = sequence.upper()
        third_bases = [seq[i] for i in range(2, len(seq), 3) if seq[i] in "ACGT"]
        if not third_bases:
            return None
        gc_count = third_bases.count("G") + third_bases.count("C")
        return round((gc_count / len(third_bases)) * 100.0, 2)

    def get_genomic_baseline(self, species: str) -> Tuple[float, float]:
        """Returns (mean_genomic_gc, std_genomic_gc) for a species."""
        if species in self.GENOMIC_ASSEMBLY_BASELINES:
            return self.GENOMIC_ASSEMBLY_BASELINES[species]
        sp_lower = species.lower()
        for k, v in self.GENOMIC_ASSEMBLY_BASELINES.items():
            if k.lower() == sp_lower:
                return v
        return (41.0, 4.5)

    def get_genome_baseline(
        self, species: str, assembly: Optional[str] = None, force_recompute: bool = False
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Retrieves genome baseline from cache or computes it if not cached."""
        table_path, stats_path = self._get_cache_paths(species)

        if not force_recompute and os.path.exists(table_path) and os.path.exists(stats_path):
            df = pd.read_csv(table_path)
            with open(stats_path, "r", encoding="utf-8") as f:
                stats = json.load(f)
            # Ensure mean_genomic_gc and std_genomic_gc are always present
            gen_mu, gen_sigma = self.get_genomic_baseline(species)
            stats["mean_genomic_gc"] = stats.get("mean_genomic_gc", gen_mu)
            stats["std_genomic_gc"] = stats.get("std_genomic_gc", gen_sigma)

            # Ensure empirical gene baseline (mean_exonic_gc, std_exonic_gc) is always populated
            if "mean_exonic_gc" not in stats or stats["mean_exonic_gc"] is None:
                if "exonic_gc" in df.columns and df["exonic_gc"].notna().sum() > 0:
                    stats["mean_exonic_gc"] = round(float(df["exonic_gc"].dropna().mean()), 2)
                    stats["std_exonic_gc"] = round(float(df["exonic_gc"].dropna().std()), 2)
                elif "cds_gc" in df.columns and df["cds_gc"].notna().sum() > 0:
                    stats["mean_exonic_gc"] = round(float(df["cds_gc"].dropna().mean()), 2)
                    stats["std_exonic_gc"] = round(float(df["cds_gc"].dropna().std()), 2)
            return df, stats

        asm = assembly or self.DEFAULT_SPECIES_ASSEMBLIES.get(species)
        if not asm:
            raise ValueError(f"Unknown assembly for species: {species}")

        return self.compute_and_cache(species, asm)

    def compute_and_cache(self, species: str, assembly: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Streams Ensembl CDS FASTA, extracts canonical CDS per gene, and saves baseline cache."""
        sp_cap = species[0].upper() + species[1:]
        url = f"https://ftp.ensembl.org/pub/current_fasta/{species.lower()}/cds/{sp_cap}.{assembly}.cds.all.fa.gz"

        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Python globin-pipeline)"})
        genes_cds: Dict[str, Dict[str, Any]] = {}

        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                with gzip.GzipFile(fileobj=resp) as gz:
                    current_gene_id = None
                    current_symbol = None
                    current_biotype = None
                    current_seq = []

                    for line in gz:
                        line_str = line.decode("utf-8", errors="ignore")
                        if line_str.startswith(">"):
                            if current_gene_id and current_seq:
                                seq_str = "".join(current_seq)
                                if (
                                    current_gene_id not in genes_cds
                                    or len(seq_str) > genes_cds[current_gene_id]["length"]
                                ):
                                    genes_cds[current_gene_id] = {
                                        "gene_id": current_gene_id,
                                        "gene_symbol": current_symbol or current_gene_id,
                                        "biotype": current_biotype or "protein_coding",
                                        "sequence": seq_str,
                                        "length": len(seq_str),
                                    }

                            current_gene_id = None
                            current_symbol = None
                            current_biotype = None
                            current_seq = []

                            parts = line_str.strip().split()
                            for p in parts:
                                if p.startswith("gene:"):
                                    current_gene_id = p.split(":")[1].split(".")[0]
                                elif p.startswith("gene_symbol:"):
                                    current_symbol = p.split(":")[1]
                                elif p.startswith("gene_biotype:"):
                                    current_biotype = p.split(":")[1]
                        else:
                            current_seq.append(line_str.strip())

                    # Final gene
                    if current_gene_id and current_seq:
                        seq_str = "".join(current_seq)
                        if (
                            current_gene_id not in genes_cds
                            or len(seq_str) > genes_cds[current_gene_id]["length"]
                        ):
                            genes_cds[current_gene_id] = {
                                "gene_id": current_gene_id,
                                "gene_symbol": current_symbol or current_gene_id,
                                "biotype": current_biotype or "protein_coding",
                                "sequence": seq_str,
                                "length": len(seq_str),
                            }
        except Exception as e:
            raise RuntimeError(f"Failed to stream or parse CDS from {url}: {e}")

        # Stream cDNA for exonic (transcript) GC
        cdna_url = f"https://ftp.ensembl.org/pub/current_fasta/{species.lower()}/cdna/{sp_cap}.{assembly}.cdna.all.fa.gz"
        genes_cdna: Dict[str, Dict[str, Any]] = {}
        try:
            req_cdna = urllib.request.Request(cdna_url, headers={"User-Agent": "Mozilla/5.0 (Python globin-pipeline)"})
            with urllib.request.urlopen(req_cdna, timeout=45) as resp:
                with gzip.GzipFile(fileobj=resp) as gz:
                    current_gene_id = None
                    current_seq = []
                    for line in gz:
                        line_str = line.decode("utf-8", errors="ignore")
                        if line_str.startswith(">"):
                            if current_gene_id and current_seq:
                                seq_str = "".join(current_seq)
                                if current_gene_id not in genes_cdna or len(seq_str) > genes_cdna[current_gene_id]["length"]:
                                    genes_cdna[current_gene_id] = {"length": len(seq_str), "sequence": seq_str}
                            current_gene_id = None
                            current_seq = []
                            for p in line_str.strip().split():
                                if p.startswith("gene:"):
                                    current_gene_id = p.split(":")[1].split(".")[0]
                        else:
                            current_seq.append(line_str.strip())
                    if current_gene_id and current_seq:
                        seq_str = "".join(current_seq)
                        if current_gene_id not in genes_cdna or len(seq_str) > genes_cdna[current_gene_id]["length"]:
                            genes_cdna[current_gene_id] = {"length": len(seq_str), "sequence": seq_str}
        except Exception as e:
            print(f"Note: cDNA streaming optional fallback for {species}: {e}")

        # Compute GC, GC3, and Exonic GC metrics
        records = []
        gc_values = []
        gc3_values = []
        exonic_gc_values = []

        for g_id, g_info in genes_cds.items():
            seq = g_info["sequence"]
            gc = self.calculate_gc(seq)
            gc3 = self.calculate_gc3(seq)

            # Exonic GC from cDNA
            cdna_info = genes_cdna.get(g_id)
            exonic_gc = self.calculate_gc(cdna_info["sequence"]) if cdna_info else gc
            exonic_len = cdna_info["length"] if cdna_info else g_info["length"]

            records.append({
                "species": species,
                "gene_id": g_id,
                "gene_symbol": g_info["gene_symbol"],
                "biotype": g_info["biotype"],
                "cds_length": g_info["length"],
                "exonic_length": exonic_len,
                "cds_gc": gc,
                "gc3": gc3,
                "exonic_gc": exonic_gc,
            })
            gc_values.append(gc)
            if gc3 is not None:
                gc3_values.append(gc3)
            if exonic_gc is not None:
                exonic_gc_values.append(exonic_gc)

        df = pd.DataFrame(records)

        arr_gc = np.array(gc_values)
        arr_gc3 = np.array(gc3_values)
        arr_exon = np.array(exonic_gc_values) if exonic_gc_values else arr_gc

        gen_mu, gen_sigma = self.get_genomic_baseline(species)
        stats = {
            "species": species,
            "assembly": assembly,
            "num_genes": len(df),
            "mean_genomic_gc": gen_mu,
            "std_genomic_gc": gen_sigma,
            "mean_cds_gc": round(float(np.mean(arr_gc)), 2),
            "std_cds_gc": round(float(np.std(arr_gc)), 2),
            "median_cds_gc": round(float(np.median(arr_gc)), 2),
            "q25_cds_gc": round(float(np.percentile(arr_gc, 25)), 2),
            "q75_cds_gc": round(float(np.percentile(arr_gc, 75)), 2),
            "mean_gc3": round(float(np.mean(arr_gc3)), 2),
            "std_gc3": round(float(np.std(arr_gc3)), 2),
            "median_gc3": round(float(np.median(arr_gc3)), 2),
            "q25_gc3": round(float(np.percentile(arr_gc3, 25)), 2),
            "q75_gc3": round(float(np.percentile(arr_gc3, 75)), 2),
            "mean_exonic_gc": round(float(np.mean(arr_exon)), 2),
            "std_exonic_gc": round(float(np.std(arr_exon)), 2),
            "median_exonic_gc": round(float(np.median(arr_exon)), 2),
            "q25_exonic_gc": round(float(np.percentile(arr_exon, 25)), 2),
            "q75_exonic_gc": round(float(np.percentile(arr_exon, 75)), 2),
        }

        # Save cache
        table_path, stats_path = self._get_cache_paths(species)
        df.to_csv(table_path, index=False, compression="gzip")
        with open(stats_path, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)

        return df, stats

    def compute_all_baselines(
        self, species_map: Optional[Dict[str, str]] = None
    ) -> pd.DataFrame:
        """Computes baselines for all species and returns a summary DataFrame."""
        mapping = species_map or self.DEFAULT_SPECIES_ASSEMBLIES
        summary_rows = []

        for sp, asm in mapping.items():
            print(f"Loading/computing genome baseline for {sp} ({asm})...")
            _, stats = self.get_genome_baseline(sp, asm)
            summary_rows.append(stats)

        df_summary = pd.DataFrame(summary_rows)
        return df_summary
