"""Data manager module for fetching and caching sequences and annotations from Ensembl."""

import os
import json
import time
import urllib.request
import urllib.error
from typing import Dict, List, Any, Optional


class DataManager:
    """Handles communication with Ensembl REST API with local caching and rate limiting."""

    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, cache_dir: str = "data/cache"):
        self.cache_dir = os.path.abspath(cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)
        self.last_request_time = 0.0
        self.min_interval = 0.08  # ~12 req/sec to stay safely under Ensembl's 15 req/sec limit

    def _rate_limit(self):
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()

    def _fetch_url(self, endpoint: str, headers: Optional[Dict[str, str]] = None, max_retries: int = 3) -> Any:
        url = f"{self.BASE_URL}{endpoint}"
        req_headers = {"Content-Type": "application/json"}
        if headers:
            req_headers.update(headers)

        for attempt in range(max_retries):
            self._rate_limit()
            req = urllib.request.Request(url, headers=req_headers)
            try:
                with urllib.request.urlopen(req, timeout=45) as resp:
                    data = resp.read()
                    if req_headers.get("Content-Type") == "text/plain":
                        return data.decode("utf-8")
                    return json.loads(data.decode("utf-8"))
            except urllib.error.HTTPError as e:
                if e.code == 429:  # Rate limited
                    retry_after = float(e.headers.get("Retry-After", 2.0))
                    time.sleep(retry_after + 1.0)
                elif e.code == 404:
                    return None
                elif e.code in [500, 502, 503, 504]:
                    if attempt < max_retries - 1:
                        time.sleep(2 * (attempt + 1))
                    else:
                        print(f"Warning: Ensembl HTTP {e.code} for {endpoint}: returning None")
                        return None
                else:
                    return None
            except Exception as e:
                if attempt < max_retries - 1:
                    time.sleep(2 * (attempt + 1))
                else:
                    print(f"Warning: request error for {endpoint}: {e}")
                    return None
        return None

    def _get_cache_path(self, key: str, ext: str = "json") -> str:
        safe_key = "".join([c if c.isalnum() or c in "._-" else "_" for c in key])
        return os.path.join(self.cache_dir, f"{safe_key}.{ext}")

    def fetch_region_genes(self, species: str, chromosome: str, start: int, end: int) -> List[Dict[str, Any]]:
        """Fetches all genes overlapping a genomic region (1-based closed query to Ensembl)."""
        cache_key = f"genes_{species}_{chromosome}_{start}_{end}"
        cache_file = self._get_cache_path(cache_key, "json")

        if os.path.exists(cache_file):
            with open(cache_file, "r", encoding="utf-8") as f:
                cached = json.load(f)
                if cached:
                    return cached

        # Ensembl overlap uses 1-based start and end
        # Our start is 0-based, so query start+1 to end
        endpoint = f"/overlap/region/{species.lower()}/{chromosome}:{start + 1}-{end}?feature=gene"
        data = self._fetch_url(endpoint)
        if data:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
            return data

        # Fallback to local baseline cache if Ensembl returned nothing or for NCBI assemblies
        baseline_file = os.path.join(os.path.dirname(self.cache_dir), "genome_baselines", f"{species.lower()}_genes_gc.csv.gz")
        if not os.path.exists(baseline_file):
            baseline_file = os.path.join("data", "genome_baselines", f"{species.lower()}_genes_gc.csv.gz")
        if os.path.exists(baseline_file):
            import pandas as pd
            df_base = pd.read_csv(baseline_file)
            sub = df_base[(df_base["chromosome"] == str(chromosome)) & (df_base["end"] >= start) & (df_base["start"] <= end)]
            if not sub.empty:
                genes = []
                for _, row in sub.iterrows():
                    genes.append({
                        "id": str(row["gene_id"]),
                        "external_name": str(row["gene_symbol"]),
                        "biotype": str(row.get("biotype", "protein_coding")),
                        "start": int(row["start"]),
                        "end": int(row["end"]),
                        "strand": 1,
                        "description": str(row["gene_symbol"]),
                    })
                with open(cache_file, "w", encoding="utf-8") as f:
                    json.dump(genes, f, indent=2)
                return genes

        return []

    def fetch_gene_details(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """Fetches full gene details including transcripts and exons."""
        if not gene_id:
            return None

        cache_key = f"gene_{gene_id}"
        cache_file = self._get_cache_path(cache_key, "json")

        if os.path.exists(cache_file):
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)

        if not str(gene_id).startswith("ENS"):
            return None

        endpoint = f"/lookup/id/{gene_id}?expand=1"
        data = self._fetch_url(endpoint)
        if data:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        return data

    def fetch_region_sequence(self, species: str, chromosome: str, start: int, end: int) -> str:
        """Fetches genomic DNA sequence for a region (start+1 to end, 1-based Ensembl or NCBI)."""
        cache_key = f"seq_{species}_{chromosome}_{start}_{end}"
        cache_file = self._get_cache_path(cache_key, "fa")

        if os.path.exists(cache_file):
            with open(cache_file, "r", encoding="utf-8") as f:
                cached = f.read().strip()
                if cached:
                    return cached

        # Check if NCBI chromosome
        if str(chromosome).startswith("NC_") or str(chromosome).startswith("NW_"):
            try:
                ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chromosome}&seq_start={start + 1}&seq_stop={end}&rettype=fasta&retmode=text"
                req = urllib.request.Request(ncbi_url, headers={"User-Agent": "Mozilla/5.0 (Python globin-pipeline)"})
                with urllib.request.urlopen(req, timeout=45) as resp:
                    lines = resp.read().decode("utf-8").splitlines()
                    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
                    if seq:
                        with open(cache_file, "w", encoding="utf-8") as f:
                            f.write(seq)
                        return seq
            except Exception as e:
                print(f"Warning: NCBI efetch sequence error for {chromosome}:{start}-{end}: {e}")
                return ""
            return ""

        endpoint = f"/sequence/region/{species.lower()}/{chromosome}:{start + 1}..{end}:1"
        seq = self._fetch_url(endpoint, headers={"Content-Type": "text/plain"})
        if seq:
            seq = seq.strip().upper()
            with open(cache_file, "w", encoding="utf-8") as f:
                f.write(seq)
            return seq

        return ""

    def fetch_transcript_cds(self, transcript_id: str) -> str:
        """Fetches coding sequence (CDS) for a transcript."""
        if not transcript_id:
            return ""

        cache_key = f"cds_{transcript_id}"
        cache_file = self._get_cache_path(cache_key, "fa")

        if os.path.exists(cache_file):
            with open(cache_file, "r", encoding="utf-8") as f:
                cached = f.read().strip()
                if cached:
                    return cached

        if not str(transcript_id).startswith("ENS"):
            return ""

        endpoint = f"/sequence/id/{transcript_id}?type=cds"
        seq = self._fetch_url(endpoint, headers={"Content-Type": "text/plain"})
        if seq:
            seq = seq.strip().upper()
            with open(cache_file, "w", encoding="utf-8") as f:
                f.write(seq)
            return seq

        return ""

    def fetch_chromosome_info(self, species: str, chromosome: str) -> Optional[Dict[str, Any]]:
        """Fetches chromosome length and metadata from Ensembl."""
        if str(chromosome).startswith("NC_") or str(chromosome).startswith("NW_"):
            return None

        cache_key = f"chrom_info_{species}_{chromosome}"
        cache_file = self._get_cache_path(cache_key, "json")

        if os.path.exists(cache_file):
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)

        endpoint = f"/assembly/info/{species.lower()}/{chromosome}"
        data = self._fetch_url(endpoint)
        if data:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        return data
