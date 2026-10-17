"""Automated globin locus and synteny anchor detector."""

import re
from typing import List, Dict, Any, Optional

from globin_pipeline.models import CandidateLocus
from globin_pipeline.exceptions import LocusNotFoundError


ALPHA_REGEX = re.compile(
    r"\b(HBA|HBA1|HBA2|HBA-T1|HBA-T2|HBAD|HBAA|HBZ|HBM|HBQ|HBQ1|HBA_like|HBA-like)\b|"
    r"(hemoglobin.*alpha|alpha.*globin|subunit.*alpha|subunit.*pi|subunit.*zeta)",
    re.IGNORECASE,
)

BETA_REGEX = re.compile(
    r"\b(HBB|HBB1|HBB2|HBD|HBE|HBE1|HBG|HBG1|HBG2|HBBP1|HBB-like|HBB_like|rho|epsilon)\b|"
    r"(hemoglobin.*beta|beta.*globin|subunit.*beta|beta-1)",
    re.IGNORECASE,
)

PROTO_REGEX = re.compile(
    r"\b(globin|globin-like|globin-1|globin-3|globin-F1|myoglobin|cytoglobin)\b",
    re.IGNORECASE,
)

ALPHA_ANCHORS = re.compile(r"\b(NPRL3|MPG|RHBDF1|LUC7L|UBE2I|FBXL16|SNRNP25)\b", re.IGNORECASE)
BETA_ANCHORS = re.compile(r"\b(LUC7L|LUC7L2|OR51|OR52|IMPG2|TTC5|PPIL3|CDCA8)\b", re.IGNORECASE)


class GlobinSyntenyDetector:
    """Scans whole-genome gene annotations to auto-detect globin clusters and syntenic anchors."""

    def __init__(self, cluster_distance_bp: int = 500000, margin_padding_bp: int = 35000):
        self.cluster_distance_bp = cluster_distance_bp
        self.margin_padding_bp = margin_padding_bp

    def detect_loci(
        self,
        genes: List[Dict[str, Any]],
        chrom_lengths: Optional[Dict[str, int]] = None,
    ) -> List[CandidateLocus]:
        """Detects and clusters globin genes and syntenic anchors across the genome.
        
        Args:
            genes: List of gene dicts with keys: symbol, description, chromosome, start, end, strand.
            chrom_lengths: Optional mapping of chromosome RefSeq/names to lengths in bp.
            
        Returns:
            List of CandidateLocus objects sorted by confidence and globin gene count.
            
        Raises:
            LocusNotFoundError: If zero globin genes or syntenic anchors are found.
        """
        globlist = []
        anchorlist = []

        for g in genes:
            sym = g.get("symbol", "")
            desc = g.get("description", "")
            target_str = f"{sym} {desc}"

            is_alpha = bool(ALPHA_REGEX.search(target_str))
            is_beta = bool(BETA_REGEX.search(target_str))
            is_proto = bool(PROTO_REGEX.search(target_str))

            is_anc_a = bool(ALPHA_ANCHORS.search(target_str))
            is_anc_b = bool(BETA_ANCHORS.search(target_str))

            if is_alpha or is_beta or is_proto:
                g_type = "alpha" if is_alpha else ("beta" if is_beta else "proto")
                globlist.append({**g, "globin_type": g_type})
            elif is_anc_a or is_anc_b:
                anchorlist.append({**g, "anchor_type": "alpha" if is_anc_a else "beta"})

        if not globlist:
            raise LocusNotFoundError(
                species="Unknown",
                details=f"Scanned {len(genes):,} genes across all chromosomes, but zero globin homologs were found.",
            )

        # Group by chromosome
        globs_by_chrom: Dict[str, List[Dict[str, Any]]] = {}
        for g in globlist:
            globs_by_chrom.setdefault(g["chromosome"], []).append(g)

        anchors_by_chrom: Dict[str, List[Dict[str, Any]]] = {}
        for a in anchorlist:
            anchors_by_chrom.setdefault(a["chromosome"], []).append(a)

        candidates: List[CandidateLocus] = []

        for chrom, chr_globs in globs_by_chrom.items():
            chr_globs = sorted(chr_globs, key=lambda x: x["start"])
            chr_anchors = anchors_by_chrom.get(chrom, [])

            # Cluster globins within cluster_distance_bp
            clusters: List[List[Dict[str, Any]]] = []
            curr_cluster = [chr_globs[0]]

            for g in chr_globs[1:]:
                if g["start"] - curr_cluster[-1]["end"] <= self.cluster_distance_bp:
                    curr_cluster.append(g)
                else:
                    clusters.append(curr_cluster)
                    curr_cluster = [g]
            clusters.append(curr_cluster)

            for cl in clusters:
                min_s = min(g["start"] for g in cl)
                max_e = max(g["end"] for g in cl)

                # Find anchors within margin of this cluster
                cl_anchors = [
                    a for a in chr_anchors
                    if (a["end"] >= (min_s - self.cluster_distance_bp) and a["start"] <= (max_e + self.cluster_distance_bp))
                ]

                all_members = cl + cl_anchors
                locus_start = max(0, min(m["start"] for m in all_members) - self.margin_padding_bp)
                
                chr_len = chrom_lengths.get(chrom, 0) if chrom_lengths else 0
                raw_end = max(m["end"] for m in all_members) + self.margin_padding_bp
                locus_end = min(chr_len, raw_end) if chr_len > 0 else raw_end

                # Determine locus identity
                n_alpha = sum(1 for g in cl if g["globin_type"] == "alpha")
                n_beta = sum(1 for g in cl if g["globin_type"] == "beta")
                n_proto = sum(1 for g in cl if g["globin_type"] == "proto")

                has_nprl3 = any("nprl3" in a["symbol"].lower() for a in cl_anchors)
                has_rhbdf1 = any("rhbdf1" in a["symbol"].lower() for a in cl_anchors)
                has_or = any("or5" in a["symbol"].lower() for a in cl_anchors)

                if n_alpha > 0 and n_beta > 0:
                    locus_type = "alpha_beta_linked"
                    locus_id = "alpha_beta_linked"
                    confidence = 0.98
                elif n_alpha > 0 or has_nprl3 or has_rhbdf1:
                    locus_type = "alpha_globin"
                    locus_id = "alpha_globin"
                    confidence = 0.95 if (n_alpha > 0 and has_nprl3) else 0.85
                elif n_beta > 0 or has_or:
                    locus_type = "beta_globin"
                    locus_id = "beta_globin"
                    confidence = 0.95 if (n_beta > 0 and has_or) else 0.85
                else:
                    locus_type = "proto_globin"
                    locus_id = "proto_globin"
                    confidence = 0.80

                notes = f"Detected {len(cl)} globins ({n_alpha} alpha, {n_beta} beta, {n_proto} proto) and {len(cl_anchors)} synteny anchors"

                candidates.append(
                    CandidateLocus(
                        locus_id=locus_id,
                        chromosome=chrom,
                        start=locus_start,
                        end=locus_end,
                        span_bp=(locus_end - locus_start),
                        globin_genes=cl,
                        anchor_genes=cl_anchors,
                        confidence=confidence,
                        suggested_type=locus_type,
                        notes=notes,
                    )
                )

        # Sort by confidence and globin gene count
        candidates.sort(key=lambda c: (c.confidence, len(c.globin_genes)), reverse=True)
        return candidates
