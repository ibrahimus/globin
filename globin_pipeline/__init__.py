"""Globin Pipeline: Comparative genomics analysis and reference atlas for vertebrate globin loci."""

__version__ = "1.0.0"

from globin_pipeline.api import GlobinPipeline
from globin_pipeline.models import LocusInterval, GeneFeature, CandidateLocus, GenomeBaseline
from globin_pipeline.fetchers.ncbi import NCBIRefSeqFetcher
from globin_pipeline.synteny import GlobinSyntenyDetector
from globin_pipeline.exceptions import (
    GlobinPipelineError,
    GenomeNotFoundError,
    AnnotationMissingError,
    ScaffoldOnlyError,
    LocusNotFoundError,
    InvalidCoordinateError,
    BaselineComputationError,
)

__all__ = [
    "__version__",
    "GlobinPipeline",
    "LocusInterval",
    "GeneFeature",
    "CandidateLocus",
    "GenomeBaseline",
    "NCBIRefSeqFetcher",
    "GlobinSyntenyDetector",
    "GlobinPipelineError",
    "GenomeNotFoundError",
    "AnnotationMissingError",
    "ScaffoldOnlyError",
    "LocusNotFoundError",
    "InvalidCoordinateError",
    "BaselineComputationError",
]
