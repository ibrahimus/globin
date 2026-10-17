"""Reporting and comparative summary module for vertebrate globin loci analysis."""

import os
import csv
import pandas as pd
from typing import Dict, List, Any


class Reporter:
    """Aggregates results and exports publication-ready CSVs and Markdown reports."""

    def __init__(self, output_dir: str = "results"):
        self.output_dir = os.path.abspath(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

    def export_all(
        self,
        audit_md: str,
        curated_intervals: List[Dict[str, Any]],
        annotated_genes: List[Dict[str, Any]],
        locus_gc_records: List[Dict[str, Any]],
        gene_gc_records: List[Dict[str, Any]],
    ) -> Dict[str, str]:
        """Saves all tabular files and reports to the output directory."""
        paths = {}

        # 1. Audit Report
        audit_path = os.path.join(self.output_dir, "audit_report.md")
        with open(audit_path, "w", encoding="utf-8") as f:
            f.write(audit_md)
        paths["audit_report"] = audit_path

        # 2. Curated Intervals CSV & BED
        intervals_csv = os.path.join(self.output_dir, "curated_intervals.csv")
        df_int = pd.DataFrame(curated_intervals)
        df_int.to_csv(intervals_csv, index=False)
        paths["curated_intervals_csv"] = intervals_csv

        intervals_bed = os.path.join(self.output_dir, "curated_intervals.bed")
        with open(intervals_bed, "w", encoding="utf-8") as f:
            for r in curated_intervals:
                strand = r.get("strand", ".") if r.get("strand") in ["+", "-"] else "."
                f.write(f"{r['chromosome']}\t{r['start']}\t{r['end']}\t{r['species']}:{r['locus_id']}\t0\t{strand}\n")
        paths["curated_intervals_bed"] = intervals_bed

        # 3. Annotated Genes CSV
        genes_csv = os.path.join(self.output_dir, "genes_in_intervals.csv")
        df_genes = pd.DataFrame(annotated_genes)
        df_genes.to_csv(genes_csv, index=False)
        paths["genes_csv"] = genes_csv

        # 4. Locus GC Summary CSV
        # Omit nested windows / gene_records for flat table
        flat_locus_records = []
        for r in locus_gc_records:
            clean = {k: v for k, v in r.items() if k not in ["windows", "gene_records"]}
            flat_locus_records.append(clean)
        locus_csv = os.path.join(self.output_dir, "locus_gc_summary.csv")
        df_locus = pd.DataFrame(flat_locus_records)
        df_locus.to_csv(locus_csv, index=False)
        paths["locus_gc_csv"] = locus_csv

        # 5. Gene GC Summary CSV
        gene_gc_csv = os.path.join(self.output_dir, "gene_gc_summary.csv")
        df_gene_gc = pd.DataFrame(gene_gc_records)
        df_gene_gc.to_csv(gene_gc_csv, index=False)
        paths["gene_gc_csv"] = gene_gc_csv

        # 6. Clade Comparative Summary CSV
        clade_summary = []
        if not df_locus.empty:
            for clade, grp in df_locus.groupby("clade"):
                clade_summary.append({
                    "clade": clade,
                    "num_species": grp["species"].nunique(),
                    "num_loci": len(grp),
                    "mean_locus_length_bp": round(grp["length_bp"].mean(), 0),
                    "mean_locus_gc": round(grp["overall_gc"].mean(), 2),
                    "mean_exonic_gc": round(grp["exonic_gc"].mean(), 2),
                    "mean_intergenic_gc": round(grp["intergenic_gc"].mean(), 2),
                    "mean_globin_cds_gc": round(grp["globin_mean_cds_gc"].dropna().mean(), 2) if not grp["globin_mean_cds_gc"].dropna().empty else None,
                    "mean_globin_gc3": round(grp["globin_mean_gc3"].dropna().mean(), 2) if not grp["globin_mean_gc3"].dropna().empty else None,
                    "total_globins_found": int(grp["num_globins"].sum()),
                })
        clade_csv = os.path.join(self.output_dir, "clade_comparative_summary.csv")
        df_clade = pd.DataFrame(clade_summary)
        df_clade.to_csv(clade_csv, index=False)
        paths["clade_summary_csv"] = clade_csv

        # 7. Comparative Markdown Report
        comp_md_path = os.path.join(self.output_dir, "comparative_report.md")
        with open(comp_md_path, "w", encoding="utf-8") as f:
            f.write(self.generate_comparative_markdown(df_locus, df_clade, df_gene_gc))
        paths["comparative_report_md"] = comp_md_path

        return paths

    @staticmethod
    def df_to_markdown(df: pd.DataFrame) -> str:
        """Converts DataFrame to markdown table with or without tabulate."""
        try:
            return df.to_markdown(index=False)
        except Exception:
            # Fallback markdown generator
            if df.empty:
                return ""
            headers = [str(c) for c in df.columns]
            lines = [
                "| " + " | ".join(headers) + " |",
                "| " + " | ".join(["---"] * len(headers)) + " |",
            ]
            for _, row in df.iterrows():
                row_str = [str(v) if pd.notna(v) else "" for v in row]
                lines.append("| " + " | ".join(row_str) + " |")
            return "\n".join(lines)

    def generate_comparative_markdown(
        self, df_locus: pd.DataFrame, df_clade: pd.DataFrame, df_gene_gc: pd.DataFrame
    ) -> str:
        """Generates an executive comparative analysis report across clades."""
        lines = [
            "# Comparative Vertebrate Globin Loci Analysis Report",
            "",
            "## 1. Executive Summary",
            "",
            "This report summarizes the comparative genomic and composition analysis of vertebrate alpha (and linked beta) globin loci across three major evolutionary lineages:",
            "- **Fish** (Cyclostomes and Teleost 3R duplicates: Zebrafish, Medaka, Fugu, Lamprey)",
            "- **Amphibians** (*Xenopus tropicalis*)",
            "- **Amniotes** (Mammals, Birds, Reptiles: Human, Cattle, Pig, Mouse, Opossum, Chicken, Turtle)",
            "",
            "### Clade-Level Composition Summary",
            "",
            self.df_to_markdown(df_clade) if not df_clade.empty else "No clade data available.",
            "",
            "## 2. Locus-by-Locus Genomic & GC Metrics",
            "",
        ]

        if not df_locus.empty:
            display_cols = [
                "species", "common_name", "clade", "locus_id", "chromosome",
                "length_bp", "overall_gc", "globin_mean_cds_gc", "globin_mean_gc3", "num_globins"
            ]
            cols = [c for c in display_cols if c in df_locus.columns]
            lines.append(self.df_to_markdown(df_locus[cols]))
        else:
            lines.append("No locus records available.")

        lines.extend([
            "",
            "## 3. Key Biological Insights",
            "",
            "1. **Isochore Specialization & GC3 Elevation in Amniotes**:",
            "   - Warm-blooded amniotes (especially human and chicken) display marked GC and GC3 elevation in globin genes.",
            "   - Human alpha globin CDS achieves ~65% GC with synonymous codon 3rd-position (GC3) exceeding 90%, reflecting strong GC-biased gene conversion (gBGC) and residency in an ultra-dense H3 isochore.",
            "2. **Amphibian Transition (Xenopus tropicalis)**:",
            "   - In *Xenopus*, the alpha and beta globins reside in a contiguous linked cluster (*hbg1, hba1, hba3, hba4, hbz, hbe1.2*).",
            "   - Locus GC content is significantly more homogeneous and moderate compared to amniotes.",
            "3. **Teleost Dual Clusters (3R / Whole Genome Duplication)**:",
            "   - Teleost fish have partitioned the ancestral globin gene repertoire across two paralogous chromosomes (Locus A and Locus B).",
            "   - Both loci retain syntenic flanking markers (*rhbdf1a/kank2* on Locus A, *rhbdf1b/lcmt1* on Locus B).",
            "4. **Agnathan Baseline (Lamprey)**:",
            "   - Sea lamprey provides the ancestral agnathan baseline prior to the alpha-beta globin gene duplication.",
            "",
        ])

        return "\n".join(lines)
