"""Parser and auditing module for incoming data batches from collaborator."""

import os
import csv
from typing import Dict, List, Any, Optional
import yaml


class BatchParser:
    """Parses and audits a dated batch folder from Mike."""

    def __init__(self, batch_dir: str):
        if not os.path.exists(batch_dir):
            if os.path.exists(os.path.join("..", batch_dir)):
                batch_dir = os.path.join("..", batch_dir)
            elif os.path.exists(os.path.basename(batch_dir.rstrip("/"))):
                batch_dir = os.path.basename(batch_dir.rstrip("/"))
        self.batch_dir = os.path.abspath(batch_dir)
        self.batch_name = os.path.basename(self.batch_dir.rstrip("/"))
        self.config_path = os.path.join(self.batch_dir, "homologues_config_multi_gene.yaml")
        self.homologues_path = os.path.join(self.batch_dir, "Homo_sapiens_selected_flanking_homologues.csv")
        self.intervals_path = os.path.join(self.batch_dir, "alpha_globin_domain_intervals.csv")
        self.descriptions_path = os.path.join(self.batch_dir, "file_descriptions.txt")

    def load_config(self) -> Dict[str, Any]:
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(f"Config file not found: {self.config_path}")
        with open(self.config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)

    def load_homologues(self) -> List[Dict[str, Any]]:
        if not os.path.exists(self.homologues_path):
            raise FileNotFoundError(f"Homologues file not found: {self.homologues_path}")
        with open(self.homologues_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)

    def load_intervals(self) -> List[Dict[str, Any]]:
        if not os.path.exists(self.intervals_path):
            raise FileNotFoundError(f"Intervals file not found: {self.intervals_path}")
        with open(self.intervals_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)

    def load_descriptions(self) -> str:
        if os.path.exists(self.descriptions_path):
            with open(self.descriptions_path, "r", encoding="utf-8") as f:
                return f.read()
        return ""

    def audit(self) -> Dict[str, Any]:
        """Runs a comprehensive integrity and biological audit on the batch."""
        config = self.load_config()
        homologues = self.load_homologues()
        intervals = self.load_intervals()

        target_species_config = set(config.get("target_species", []))
        query_species = config.get("query_species", "Homo_sapiens")
        query_genes = config.get("query_genes", [])

        species_in_hom = set(r["target_species"] for r in homologues if r.get("target_species"))
        species_in_int = set(r["target_species"] for r in intervals if r.get("target_species"))

        missing_in_intervals = target_species_config - species_in_int
        missing_in_homologues = target_species_config - species_in_hom

        # Inspect each species in intervals
        interval_audits = []
        for inter in intervals:
            sp = inter["target_species"]
            chrom = str(inter["Chromosome"])
            try:
                start = int(inter["Start"])
                end = int(inter["End"])
            except ValueError:
                start, end = -1, -1

            strand = inter.get("Strand", "")
            left_gene = inter.get("left_gene_name", "")
            right_gene = inter.get("right_gene_name", "")

            issues = []
            flags = []

            # 1. Coordinate check
            if start < 0 or end < 0 or start >= end:
                issues.append(f"Invalid coordinate range: {start} to {end}")

            # 2. Fixed endpoint check
            if left_gene == "Fixed endpoint" or right_gene == "Fixed endpoint":
                flags.append("Contains 'Fixed endpoint' boundary (arbitrary window / missing flank)")

            # 3. Find homologues for this species
            sp_homs = [h for h in homologues if h.get("target_species") == sp]
            hba_homs = [h for h in sp_homs if h.get("query_gene") == "HBA2" and h.get("status") == "Found"]

            if not hba_homs:
                flags.append("HBA2 homologue not found in Ensembl Compara (e.g. cyclostome / divergent)")
            else:
                for hba in hba_homs:
                    hba_chr = str(hba.get("Chromosome"))
                    hba_name = hba.get("gene_name", "HBA2")
                    try:
                        hba_start = float(hba.get("Start", -1))
                        hba_end = float(hba.get("End", -1))
                    except ValueError:
                        hba_start, hba_end = -1, -1

                    # Check chromosome mismatch
                    if hba_chr != chrom:
                        issues.append(
                            f"Chromosome mismatch: HBA2 ortholog ({hba_name}) is on chr{hba_chr}, "
                            f"but interval is on chr{chrom} (indicates duplicate locus / 3R paralog!)"
                        )
                    else:
                        # Check coordinate inclusion
                        if hba_start < start or hba_end > end:
                            issues.append(
                                f"HBA2 ortholog ({hba_name} at {hba_start:.0f}-{hba_end:.0f}) "
                                f"is OUTSIDE the interval ({start}-{end})! Synteny assumption failed."
                            )

            interval_audits.append({
                "species": sp,
                "common_name": inter.get("species_common_name", ""),
                "assembly": inter.get("assembly_name", ""),
                "chromosome": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "left_gene": left_gene,
                "right_gene": right_gene,
                "issues": issues,
                "flags": flags,
            })

        return {
            "batch_name": self.batch_name,
            "query_species": query_species,
            "query_genes": query_genes,
            "target_species_count": len(target_species_config),
            "missing_in_intervals": list(missing_in_intervals),
            "missing_in_homologues": list(missing_in_homologues),
            "interval_audits": interval_audits,
        }

    def generate_markdown_report(self, audit_result: Optional[Dict[str, Any]] = None) -> str:
        """Formats the audit results into a detailed Markdown report."""
        if audit_result is None:
            audit_result = self.audit()

        lines = [
            f"# Batch Audit Report: `{audit_result['batch_name']}`",
            "",
            f"- **Query Species**: `{audit_result['query_species']}`",
            f"- **Target Species in Config**: {audit_result['target_species_count']}",
            f"- **Missing in Intervals Table**: {', '.join(audit_result['missing_in_intervals']) or 'None'}",
            "",
            "## Species Interval Audit Summary",
            "",
            "| Species | Common Name | Assembly | Chromosome | Span (bp) | Flanks | Status |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
        ]

        for item in audit_result["interval_audits"]:
            sp = item["species"]
            cname = item["common_name"]
            asm = item["assembly"]
            chrom = item["chromosome"]
            span = item["end"] - item["start"] if item["end"] > item["start"] else 0
            flanks = f"{item['left_gene']} ... {item['right_gene']}"

            status_badges = []
            if item["issues"]:
                status_badges.append(f"❌ **FAIL**: {'; '.join(item['issues'])}")
            elif item["flags"]:
                status_badges.append(f"⚠️ **WARN**: {'; '.join(item['flags'])}")
            else:
                status_badges.append("✅ **OK**")

            lines.append(f"| `{sp}` | {cname} | `{asm}` | `{chrom}` | {span:,} | {flanks} | {'<br>'.join(status_badges)} |")

        lines.append("")
        lines.append("## Critical Discrepancies Requiring Curation")
        lines.append("")
        lines.append("1. **Xenopus tropicalis Alpha Globins Excluded**:")
        lines.append("   - Left boundary was set to `rhbdf1` (40.71 Mb), but alpha globins (`hba1, hba3, hba4, hbz`) reside at 40.36-40.46 Mb.")
        lines.append("2. **Teleost Dual Loci (Danio rerio, Oryzias latipes, Takifugu rubripes)**:")
        lines.append("   - Ensembl Compara HBA2 orthologs are on Chr 12, 19, and 1 (Locus B), while the intervals are defined on Chr 3, 8, and 5 (Locus A). Both loci must be tracked.")
        lines.append("3. **Western Painted Turtle (Chrysemys picta bellii) Omission**:")
        lines.append("   - Missing from intervals table due to inverted flank coordinate calculation on scaffold KK088206.1.")
        lines.append("4. **Hedgehog (Erinaceus europaeus) Assembly Fragmentation**:")
        lines.append("   - GeneScaffold_103 contains NPRL3/RHBDF1 but no globins (globins are on unplaced scaffolds).")
        lines.append("5. **Human Reference Omission**:")
        lines.append("   - Homo sapiens (the reference amniote query species) is omitted from the intervals table.")
        lines.append("")

        return "\n".join(lines)
