"""Unit and integration tests for NCBI RefSeq Fetcher, Synteny Detector, and GlobinPipeline API."""

import pytest
from globin_pipeline import (
    GlobinPipeline,
    NCBIRefSeqFetcher,
    GlobinSyntenyDetector,
    GenomeNotFoundError,
    LocusNotFoundError,
    InvalidCoordinateError,
)


def test_api_bundled_data():
    """Verify GlobinPipeline initializes and accesses bundled 44-species reference data."""
    pipeline = GlobinPipeline()
    df_intervals = pipeline.get_curated_intervals()
    assert len(df_intervals) >= 80
    assert df_intervals["species"].nunique() >= 44

    manifest = pipeline.get_species_manifest()
    sp_dict = manifest.get("species", manifest)
    assert len(sp_dict) >= 44
    assert "zootoca_vivipara" in sp_dict
    assert "eptatretus_burgeri" in sp_dict
    assert "lampetra_fluviatilis" in sp_dict
    assert "pelodiscus_sinensis" in sp_dict


def test_ncbi_fetcher_url_resolution():
    """Verify resolve_ftp_urls resolves valid RefSeq accession and constructs file paths."""
    fetcher = NCBIRefSeqFetcher()
    urls = fetcher.resolve_ftp_urls("GCF_965207765.1")
    assert "rZooViv2.hap1.1" in urls["assembly_name"]
    assert urls["assembly_report"].endswith("_assembly_report.txt")
    assert urls["cds"].endswith("_cds_from_genomic.fna.gz")
    assert urls["gff"].endswith("_genomic.gff.gz")
    assert urls["fna"].endswith("_genomic.fna.gz")


def test_ncbi_fetcher_invalid_accession():
    """Verify resolve_ftp_urls raises GenomeNotFoundError on non-existent accession."""
    fetcher = NCBIRefSeqFetcher()
    with pytest.raises(GenomeNotFoundError):
        fetcher.resolve_ftp_urls("GCF_999999999.1")


def test_synteny_detector_synthetic():
    """Verify GlobinSyntenyDetector clusters globin homologs and synteny anchors."""
    detector = GlobinSyntenyDetector(cluster_distance_bp=100000, margin_padding_bp=10000)
    synthetic_genes = [
        {"gene_id": "1", "symbol": "RHBDF1", "description": "rhomboid", "chromosome": "chr1", "start": 100000, "end": 120000, "strand": 1},
        {"gene_id": "2", "symbol": "MPG", "description": "glycosylase", "chromosome": "chr1", "start": 130000, "end": 140000, "strand": -1},
        {"gene_id": "3", "symbol": "NPRL3", "description": "GATOR complex", "chromosome": "chr1", "start": 145000, "end": 160000, "strand": -1},
        {"gene_id": "4", "symbol": "HBA1", "description": "hemoglobin subunit alpha", "chromosome": "chr1", "start": 170000, "end": 172000, "strand": 1},
        {"gene_id": "5", "symbol": "HBA2", "description": "hemoglobin subunit alpha", "chromosome": "chr1", "start": 180000, "end": 182000, "strand": 1},
    ]

    candidates = detector.detect_loci(synthetic_genes)
    assert len(candidates) == 1
    locus = candidates[0]
    assert locus.suggested_type == "alpha_globin"
    assert locus.chromosome == "chr1"
    assert locus.start <= 100000
    assert locus.end >= 182000
    assert len(locus.globin_genes) == 2
    assert len(locus.anchor_genes) == 3


def test_synteny_detector_no_globins():
    """Verify GlobinSyntenyDetector raises LocusNotFoundError if no globins are present."""
    detector = GlobinSyntenyDetector()
    non_globin_genes = [
        {"gene_id": "1", "symbol": "ACTB", "description": "actin", "chromosome": "chr1", "start": 1000, "end": 2000, "strand": 1},
        {"gene_id": "2", "symbol": "GAPDH", "description": "dehydrogenase", "chromosome": "chr1", "start": 3000, "end": 4000, "strand": 1},
    ]
    with pytest.raises(LocusNotFoundError):
        detector.detect_loci(non_globin_genes)


def test_api_invalid_coordinates():
    """Verify GlobinPipeline raises InvalidCoordinateError if start >= end."""
    pipeline = GlobinPipeline()
    with pytest.raises(InvalidCoordinateError):
        pipeline.add_genome(
            "GCF_965207765.1",
            common_name="Test Lizard",
            clade="amniote_reptile",
            chrom="NC_142596.1",
            start=500000,
            end=400000,
        )
