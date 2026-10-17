"""Unified NCBI RefSeq genome fetcher and streaming processor."""

import os
import re
import gzip
import time
import urllib.request
import urllib.error
import urllib.parse
from typing import Dict, List, Tuple, Any, Optional

from globin_pipeline.exceptions import (
    GenomeNotFoundError,
    AnnotationMissingError,
    ScaffoldOnlyError,
    BaselineComputationError,
)


def calculate_gc(seq: str) -> Tuple[float, int]:
    """Calculates GC percentage and count of valid nucleotides (A, C, G, T)."""
    if not seq:
        return 0.0, 0
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    valid = gc + seq.count("A") + seq.count("T")
    if valid == 0:
        return 0.0, 0
    return round((gc / valid) * 100.0, 2), valid


def calculate_gc3(cds_seq: str) -> Optional[float]:
    """Calculates codon GC3 percentage at 3rd codon positions."""
    if not cds_seq or len(cds_seq) < 3:
        return None
    seq = cds_seq.upper()
    third_bases = [seq[i] for i in range(2, len(seq), 3) if seq[i] in "ACGT"]
    if not third_bases:
        return None
    gc = third_bases.count("G") + third_bases.count("C")
    return round((gc / len(third_bases)) * 100.0, 2)


class NCBIRefSeqFetcher:
    """Streams and parses official NCBI RefSeq assembly reports, CDS, GFF, and genomic FASTA."""

    BASE_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

    def __init__(self, user_agent: str = "Mozilla/5.0 (Bioinformatics Globin Project)"):
        self.headers = {"User-Agent": user_agent}

    def resolve_ftp_urls(self, accession: str, timeout: int = 45) -> Dict[str, str]:
        """Resolves an NCBI accession (e.g. GCF_000001405.40) to direct FTP download URLs.
        
        Raises:
            GenomeNotFoundError: If accession directory is not found on NCBI FTP.
            AnnotationMissingError: If GFF3 or CDS files are missing from the release.
        """
        accession = accession.strip()
        parts = accession.split("_")
        if len(parts) < 2 or not parts[0].startswith("GC"):
            raise GenomeNotFoundError(accession, "Accession must start with 'GCF_' (RefSeq) or 'GCA_' (GenBank).")

        prefix = parts[0]
        num_part = parts[1].split(".")[0]
        if len(num_part) < 9:
            num_part = num_part.zfill(9)

        p1, p2, p3 = num_part[0:3], num_part[3:6], num_part[6:9]
        dir_url = f"{self.BASE_FTP}/{prefix}/{p1}/{p2}/{p3}/"

        try:
            req = urllib.request.Request(dir_url, headers=self.headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                html = resp.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise GenomeNotFoundError(accession, f"Directory not found on NCBI FTP ({dir_url}).")
            raise GenomeNotFoundError(accession, f"HTTP error {e.code}: {e.reason}")
        except Exception as e:
            raise GenomeNotFoundError(accession, f"Network connection error: {e}")

        # Locate assembly subfolder matching accession
        matches = re.findall(rf'href="({re.escape(accession)}_[^"]+)/?"', html)
        if not matches:
            raise GenomeNotFoundError(accession, f"No assembly subfolder matching '{accession}' found in {dir_url}.")

        subdir = matches[0].strip("/")
        base_file_url = f"{dir_url}{subdir}/{subdir}"

        urls = {
            "assembly_name": subdir.replace(f"{accession}_", ""),
            "assembly_report": f"{base_file_url}_assembly_report.txt",
            "cds": f"{base_file_url}_cds_from_genomic.fna.gz",
            "gff": f"{base_file_url}_genomic.gff.gz",
            "fna": f"{base_file_url}_genomic.fna.gz",
        }

        # Check for presence of annotation
        # RefSeq assemblies with active annotation contain cds_from_genomic.fna.gz
        if f"{subdir}_cds_from_genomic.fna.gz" not in html and f"{subdir}_genomic.gff.gz" not in html:
            # Let's probe the CDS url with a lightweight HEAD request
            try:
                head_req = urllib.request.Request(urls["cds"], headers=self.headers, method="HEAD")
                with urllib.request.urlopen(head_req, timeout=10):
                    pass
            except Exception:
                raise AnnotationMissingError(accession, urls["assembly_name"])

        return urls

    def fetch_assembly_report(self, report_url: str, timeout: int = 45) -> Dict[str, Any]:
        """Streams assembly_report.txt and extracts chromosome lengths, roles, and RefSeq accessions."""
        req = urllib.request.Request(report_url, headers=self.headers)
        chrom_map: Dict[str, Dict[str, Any]] = {}
        has_chromosomes = False

        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                for line in resp:
                    s = line.decode("utf-8", errors="ignore").strip()
                    if not s or s.startswith("#"):
                        continue
                    parts = s.split("\t")
                    if len(parts) >= 9:
                        chr_name = parts[0]
                        role = parts[1]
                        refseq_acc = parts[6]
                        seq_len = int(parts[8])

                        if role == "assembled-molecule":
                            has_chromosomes = True

                        chrom_info = {
                            "chr_name": chr_name,
                            "role": role,
                            "length": seq_len,
                            "refseq_acc": refseq_acc,
                        }
                        chrom_map[refseq_acc] = chrom_info
                        chrom_map[chr_name] = chrom_info
        except Exception as e:
            raise GenomeNotFoundError(report_url, f"Failed to stream assembly report: {e}")

        return {
            "has_chromosomes": has_chromosomes,
            "chromosomes": chrom_map,
        }

    def stream_cds_metrics(self, cds_url: str, timeout: int = 90) -> Dict[str, Dict[str, Any]]:
        """Streams CDS FASTA and calculates canonical CDS length, CDS GC, and codon GC3 per gene."""
        req = urllib.request.Request(cds_url, headers=self.headers)
        gene_cds_info: Dict[str, Dict[str, Any]] = {}

        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
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
        except Exception as e:
            raise BaselineComputationError(f"Failed to stream CDS FASTA from {cds_url}: {e}")

        if not gene_cds_info:
            raise BaselineComputationError(f"Zero CDS sequences found in {cds_url}.")

        return gene_cds_info

    def stream_gff(
        self,
        gff_url: str,
        valid_cds_gids: Optional[set] = None,
        target_interval: Optional[Tuple[str, int, int]] = None,
        timeout: int = 120,
    ) -> Dict[str, Any]:
        """Streams GFF3 and parses gene coordinates, symbols, biotypes, descriptions, and models.
        
        If valid_cds_gids is supplied, filters strictly for protein-coding genes with valid CDS.
        """
        req = urllib.request.Request(gff_url, headers=self.headers)
        gene_coords: Dict[str, Tuple[str, int, int, str, str, str]] = {}
        genes_by_chrom: Dict[str, List[Tuple[str, int, int]]] = {}
        all_genes_list: List[Dict[str, Any]] = []

        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                with gzip.GzipFile(fileobj=resp) as gz:
                    for line in gz:
                        s = line.decode("utf-8", errors="ignore")
                        if s.startswith("#"):
                            continue
                        parts = s.split("\t")
                        if len(parts) < 9 or parts[2] != "gene":
                            continue

                        chrom = parts[0]
                        st = int(parts[3])
                        en = int(parts[4])
                        strand = 1 if parts[6] == "+" else (-1 if parts[6] == "-" else 1)
                        attr = parts[8]

                        m_gid = re.search(r"GeneID:(\d+)", attr)
                        if not m_gid:
                            continue
                        gid = f"GeneID:{m_gid.group(1)}"

                        # If filter provided, skip non-coding RNAs lacking CDS
                        if valid_cds_gids is not None and gid not in valid_cds_gids:
                            continue

                        m_sym = re.search(r"Name=([^;]+)", attr)
                        sym = m_sym.group(1) if m_sym else gid
                        m_bio = re.search(r"gene_biotype=([^;]+)", attr)
                        biotype = m_bio.group(1).strip() if m_bio else "protein_coding"
                        m_desc = re.search(r"description=([^;]+)", attr)
                        desc = urllib.parse.unquote(m_desc.group(1)) if m_desc else sym

                        if gid in gene_coords:
                            old_chr, old_s, old_e, old_sym, old_bio, old_desc = gene_coords[gid]
                            if old_chr == chrom:
                                st = min(st, old_s)
                                en = max(en, old_e)

                        gene_coords[gid] = (chrom, st, en, sym, biotype, desc)
                        genes_by_chrom.setdefault(chrom, []).append((gid, st, en))
                        all_genes_list.append({
                            "gene_id": gid,
                            "symbol": sym,
                            "chromosome": chrom,
                            "start": st,
                            "end": en,
                            "strand": strand,
                            "biotype": biotype,
                            "description": desc,
                        })
        except Exception as e:
            raise BaselineComputationError(f"Failed to stream GFF3 from {gff_url}: {e}")

        return {
            "gene_coords": gene_coords,
            "genes_by_chrom": genes_by_chrom,
            "all_genes": all_genes_list,
        }
