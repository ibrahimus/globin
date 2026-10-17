"""Custom domain exceptions for the globin comparative genomics pipeline."""


class GlobinPipelineError(Exception):
    """Base exception for all globin-pipeline errors."""
    pass


class GenomeNotFoundError(GlobinPipelineError):
    """Raised when an accession cannot be located on NCBI or FTP servers."""
    def __init__(self, accession: str, reason: str = ""):
        msg = f"Genome accession '{accession}' could not be resolved or found on NCBI."
        if reason:
            msg += f" Details: {reason}"
        super().__init__(msg)
        self.accession = accession


class AnnotationMissingError(GlobinPipelineError):
    """Raised when an assembly exists but lacks RefSeq/Ensembl gene annotation releases."""
    def __init__(self, accession: str, assembly_name: str = ""):
        msg = (
            f"Assembly '{accession}' ({assembly_name}) has no annotation release (GFF3 / CDS files missing).\n"
            f"Remediation: Ensure the accession is an annotated RefSeq build (prefixed with 'GCF_') "
            f"with an active NCBI Annotation Release, or supply your own curated GFF/CDS files."
        )
        super().__init__(msg)
        self.accession = accession


class ScaffoldOnlyError(GlobinPipelineError):
    """Raised when an assembly lacks chromosome-level anchoring and physical telomere distances cannot be measured."""
    def __init__(self, accession: str, scaffold_name: str = ""):
        msg = (
            f"Assembly '{accession}' is fragmented at contig/scaffold level (locus resides on '{scaffold_name}').\n"
            f"Physical and relative distance to chromosome telomeres cannot be calculated without chromosome-scale anchoring."
        )
        super().__init__(msg)
        self.accession = accession


class LocusNotFoundError(GlobinPipelineError):
    """Raised when automated synteny detection fails to locate globin clusters or anchors."""
    def __init__(self, species: str, accession: str = "", details: str = ""):
        msg = (
            f"Could not auto-detect globin locus in '{species}' ({accession}).\n"
            f"No recognizable globin genes (HBA, HBB, HBZ, HBD, or proto-globins) or syntenic anchors "
            f"(NPRL3, MPG, RHBDF1, LUC7L) were identified in the annotation.\n"
            f"Remediation: Specify locus boundaries manually via --chrom, --start, and --end flags."
        )
        if details:
            msg += f"\nDetails: {details}"
        super().__init__(msg)
        self.species = species


class InvalidCoordinateError(GlobinPipelineError):
    """Raised when supplied genomic coordinates are invalid or outside chromosome boundaries."""
    def __init__(self, chrom: str, start: int, end: int, chrom_len: int = 0):
        msg = f"Invalid coordinates {chrom}:{start}-{end}."
        if start >= end:
            msg += f" Start ({start}) must be strictly less than end ({end})."
        elif chrom_len > 0 and (start < 0 or end > chrom_len):
            msg += f" Coordinates exceed chromosome boundary (0 to {chrom_len:,} bp)."
        super().__init__(msg)


class BaselineComputationError(GlobinPipelineError):
    """Raised when baseline calculation fails due to corrupted files or missing CDS features."""
    pass
