"""High-level Python API for the globin comparative pipeline."""

import os
import json
import yaml
import importlib.resources
import pandas as pd
from typing import Dict, List, Any, Optional

from globin_pipeline.models import LocusInterval, CandidateLocus
from globin_pipeline.fetchers.ncbi import NCBIRefSeqFetcher
from globin_pipeline.synteny import GlobinSyntenyDetector
from globin_pipeline.baseline_manager import BaselineManager
from globin_pipeline.visualizer import Visualizer
from globin_pipeline.exceptions import (
    GenomeNotFoundError,
    AnnotationMissingError,
    LocusNotFoundError,
    InvalidCoordinateError,
)


class GlobinPipeline:
    """Primary interface for running comparative globin analysis and adding new genomes."""

    def __init__(self, data_dir: Optional[str] = None):
        # Locate package data or fallback to local directory
        self.data_dir = data_dir
        if not self.data_dir:
            pkg_data = os.path.join(os.path.dirname(__file__), "data")
            if os.path.exists(pkg_data):
                self.data_dir = pkg_data
            else:
                self.data_dir = "data"

        self.baseline_manager = BaselineManager(cache_dir="data/genome_baselines")
        self.fetcher = NCBIRefSeqFetcher()
        self.synteny_detector = GlobinSyntenyDetector()
        self.visualizer = Visualizer(self.baseline_manager)

    def get_curated_intervals(self) -> pd.DataFrame:
        """Returns the master DataFrame of all 80 curated intervals across 44 vertebrate species."""
        csv_path = os.path.join(self.data_dir, "curated_intervals.csv")
        if not os.path.exists(csv_path):
            csv_path = "data/curated_intervals.csv"
        return pd.read_csv(csv_path)

    def get_species_manifest(self) -> Dict[str, Any]:
        """Returns the full declarative species manifest dictionary."""
        yaml_path = os.path.join(self.data_dir, "species_manifest.yaml")
        if not os.path.exists(yaml_path):
            yaml_path = "config/species_manifest.yaml"
        with open(yaml_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)

    def get_gene_zscores(self, results_dir: str = "results") -> pd.DataFrame:
        """Loads standardized gene Z-scores for all 1,017 analyzed genes."""
        csv_path = os.path.join(results_dir, "gene_zscores.csv")
        if not os.path.exists(csv_path):
            csv_path = os.path.join("..", results_dir, "gene_zscores.csv")
        if os.path.exists(csv_path):
            return pd.read_csv(csv_path)
        raise FileNotFoundError(f"Gene Z-scores table not found at {csv_path}. Run pipeline first.")

    def get_cross_species_summary(self, results_dir: str = "results") -> pd.DataFrame:
        """Builds or loads the 44-species locus-level summary dataset."""
        df_genes = self.get_gene_zscores(results_dir)
        return self.visualizer.build_cross_species_summary(df_genes)

    def add_genome(
        self,
        accession: str,
        common_name: str,
        clade: str,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        locus_id: Optional[str] = None,
        auto_detect: bool = True,
        save_baseline: bool = True,
    ) -> Dict[str, Any]:
        """Imports an NCBI RefSeq genome by accession, streams annotations, detects or accepts globin loci,
        and computes clean empirical protein-coding baselines.
        
        Args:
            accession: RefSeq accession (e.g. 'GCF_965207765.1').
            common_name: User-facing common name (e.g. 'Viviparous Lizard').
            clade: Taxonomic clade (e.g. 'amniote_reptile', 'amniote_mammal', 'amniote_bird').
            chrom: Optional chromosome accession (e.g. 'NC_142596.1').
            start: Optional 0-based or 1-based start coordinate.
            end: Optional end coordinate.
            locus_id: Optional locus identifier ('alpha_globin', 'beta_globin', 'alpha_beta_linked').
            auto_detect: If True, uses GlobinSyntenyDetector to identify candidate loci if coordinates are omitted.
            save_baseline: If True, saves baseline CSV and JSON stats to data/genome_baselines/.
            
        Returns:
            Dictionary containing genome metadata, detected/curated loci, and baseline statistics.
        """
        urls = self.fetcher.resolve_ftp_urls(accession)
        assembly_name = urls["assembly_name"]
        
        # 1. Fetch assembly metadata
        asm_meta = self.fetcher.fetch_assembly_report(urls["assembly_report"])
        chrom_lens = {k: v["length"] for k, v in asm_meta["chromosomes"].items()}

        # 2. Stream CDS for protein-coding filter
        cds_metrics = self.fetcher.stream_cds_metrics(urls["cds"])
        valid_cds_gids = set(cds_metrics.keys())

        # 3. Stream GFF3 with strict protein-coding filtering
        gff_data = self.fetcher.stream_gff(urls["gff"], valid_cds_gids=valid_cds_gids)
        all_genes = gff_data["all_genes"]

        # 4. Resolve Locus Coordinates (Manual or Auto-Detect)
        chosen_loci: List[CandidateLocus] = []

        if chrom and start is not None and end is not None:
            if start >= end:
                raise InvalidCoordinateError(chrom, start, end)
            loc_id = locus_id or "alpha_globin"
            chosen_loci.append(
                CandidateLocus(
                    locus_id=loc_id,
                    chromosome=chrom,
                    start=start,
                    end=end,
                    span_bp=(end - start),
                    notes="Manually specified locus coordinates",
                )
            )
        elif auto_detect:
            detected = self.synteny_detector.detect_loci(all_genes, chrom_lengths=chrom_lens)
            chosen_loci.extend(detected)
        else:
            raise LocusNotFoundError(
                species=common_name,
                accession=accession,
                details="No coordinates provided and auto_detect=False.",
            )

        species_slug = common_name.lower().replace(" ", "_").replace("(", "").replace(")", "").replace(".", "")

        return {
            "accession": accession,
            "assembly_name": assembly_name,
            "species_slug": species_slug,
            "common_name": common_name,
            "clade": clade,
            "total_protein_coding_genes": len(valid_cds_gids),
            "chromosomes_available": len(chrom_lens),
            "candidate_loci": chosen_loci,
            "urls": urls,
        }
