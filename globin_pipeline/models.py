"""Data models and type definitions for the globin comparative pipeline."""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional


@dataclass
class GeneFeature:
    """Represents an individual gene model and its sequence metrics."""
    gene_id: str
    symbol: str
    chromosome: str
    start: int
    end: int
    strand: int
    biotype: str = "protein_coding"
    description: str = ""
    is_globin: bool = False
    comparative_group: str = "other_genes"
    whole_gene_len: Optional[int] = None
    whole_gene_gc: Optional[float] = None
    cds_len: Optional[int] = None
    cds_gc: Optional[float] = None
    gc3: Optional[float] = None
    z_score_whole: Optional[float] = None
    z_score_cds: Optional[float] = None
    z_score_gc3: Optional[float] = None


@dataclass
class LocusInterval:
    """Represents a curated globin locus interval."""
    species: str
    common_name: str
    clade: str
    locus_id: str
    assembly: str
    chromosome: str
    start: int
    end: int
    span_bp: int
    strand: str
    left_gene: str = ""
    right_gene: str = ""
    status: str = "active"
    notes: str = ""
    genes: List[GeneFeature] = field(default_factory=list)
    chrom_length: Optional[int] = None
    centromere_bp: Optional[int] = None
    dist_to_telomere_bp: Optional[int] = None
    rel_dist_to_telomere: Optional[float] = None

    @property
    def is_subtelomeric(self) -> bool:
        if self.rel_dist_to_telomere is not None:
            return self.rel_dist_to_telomere < 0.05
        return False


@dataclass
class CandidateLocus:
    """Represents an auto-detected candidate globin locus."""
    locus_id: str
    chromosome: str
    start: int
    end: int
    span_bp: int
    globin_genes: List[Dict[str, Any]] = field(default_factory=list)
    anchor_genes: List[Dict[str, Any]] = field(default_factory=list)
    confidence: float = 1.0
    suggested_type: str = "alpha_globin"  # alpha_globin, beta_globin, linked, or proto_globin
    notes: str = ""


@dataclass
class GenomeBaseline:
    """Represents whole-genome baseline statistics for a species."""
    species: str
    assembly: str
    total_genes: int
    genomic_gc: float
    mean_whole_gene_len: float
    std_whole_gene_len: float
    median_whole_gene_len: float
    mean_whole_gene_gc: float
    std_whole_gene_gc: float
    median_whole_gene_gc: float
    mean_cds_gc: float
    std_cds_gc: float
    median_cds_gc: float
    mean_gc3: Optional[float] = None
    std_gc3: Optional[float] = None
    median_gc3: Optional[float] = None
