"""GC content and isochore analysis module."""

import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from globin_pipeline.data_manager import DataManager


class GCAnalyzer:
    """Calculates GC content metrics across compartments, codons (GC3), and isochore windows."""

    def __init__(self, data_manager: Optional[DataManager] = None):
        self.dm = data_manager or DataManager()

    @staticmethod
    def calculate_gc(sequence: str) -> float:
        """Calculates GC percentage of a DNA sequence, ignoring Ns and gaps."""
        if not sequence:
            return 0.0
        seq = sequence.upper()
        gc_count = seq.count("G") + seq.count("C")
        valid_bases = gc_count + seq.count("A") + seq.count("T")
        if valid_bases == 0:
            return 0.0
        return round((gc_count / valid_bases) * 100.0, 2)

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        if not sequence:
            return ""
        complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement)[::-1]

    @staticmethod
    def calculate_gc3(cds_sequence: str) -> Optional[float]:
        """Calculates GC3 (GC percentage at the third position of codons)."""
        if not cds_sequence or len(cds_sequence) < 3:
            return None
        seq = cds_sequence.upper()
        # Extract third codon positions (0-indexed: 2, 5, 8...)
        third_bases = [seq[i] for i in range(2, len(seq), 3) if seq[i] in "ACGT"]
        if not third_bases:
            return None
        gc_count = third_bases.count("G") + third_bases.count("C")
        return round((gc_count / len(third_bases)) * 100.0, 2)

    def analyze_sliding_windows(
        self, sequence: str, window_size: int = 10000, step_size: int = 2000
    ) -> Dict[str, Any]:
        """Calculates GC% across sliding windows to characterize isochores."""
        seq_len = len(sequence)
        if seq_len < window_size:
            overall = self.calculate_gc(sequence)
            return {
                "window_size": window_size,
                "step_size": step_size,
                "num_windows": 1,
                "mean_gc": overall,
                "std_gc": 0.0,
                "min_gc": overall,
                "max_gc": overall,
                "windows": [(0, seq_len, overall)],
            }

        windows = []
        gc_values = []
        for start in range(0, seq_len - window_size + 1, step_size):
            end = start + window_size
            win_seq = sequence[start:end]
            win_gc = self.calculate_gc(win_seq)
            windows.append((start, end, win_gc))
            gc_values.append(win_gc)

        arr = np.array(gc_values)
        return {
            "window_size": window_size,
            "step_size": step_size,
            "num_windows": len(windows),
            "mean_gc": round(float(np.mean(arr)), 2),
            "std_gc": round(float(np.std(arr)), 2),
            "min_gc": round(float(np.min(arr)), 2),
            "max_gc": round(float(np.max(arr)), 2),
            "windows": windows,
        }

    def analyze_locus(
        self, interval: Dict[str, Any], annotated_genes: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Performs full compartment and isochore GC analysis for a locus."""
        species = interval["species"]
        chrom = interval["chromosome"]
        start = interval["start"]
        end = interval["end"]
        locus_id = interval.get("locus_id", "alpha_globin")
        clade = interval.get("clade", "unknown")

        # 1. Fetch locus sequence
        locus_seq = self.dm.fetch_region_sequence(species, chrom, start, end)
        locus_length = len(locus_seq)
        overall_locus_gc = self.calculate_gc(locus_seq)

        # 2. Sliding windows
        sw_10kb = self.analyze_sliding_windows(locus_seq, window_size=10000, step_size=2000)
        sw_100kb = self.analyze_sliding_windows(locus_seq, window_size=100000, step_size=20000)

        # 3. Analyze individual genes (CDS GC and GC3)
        gene_gc_records = []
        all_cds_gc = []
        all_gc3 = []
        globin_cds_gc = []
        globin_gc3 = []

        # Create locus-relative mask for exons to calculate intergenic / intronic GC
        # Mask: 0 = intergenic, 1 = genic
        locus_mask = np.zeros(locus_length, dtype=np.uint8)

        for g in annotated_genes:
            g_id = g["gene_id"]
            try:
                details = self.dm.fetch_gene_details(g_id)
            except Exception as e:
                print(f"Warning: could not fetch gene details for {g_id}: {e}")
                details = None
            cds_gc = None
            gc3 = None
            exonic_gc = None
            intronic_gc = None
            cds_len = 0
            exonic_len = 0
            intronic_len = 0

            # 1. Determine canonical transcript and true gene start and end
            target_t = None
            if details and "Transcript" in details:
                from globin_pipeline.annotator import GeneAnnotator
                is_glo = str(g.get("is_globin", False)).lower() == "true"
                target_t = GeneAnnotator.get_canonical_transcript(details, is_globin=is_glo)

            g_s = g["start"]
            g_e = g["end"]
            if target_t and target_t.get("start") and target_t.get("end"):
                g_s = target_t.get("start") - 1
                g_e = target_t.get("end")
            elif details and details.get("start") and details.get("end"):
                g_s = details.get("start") - 1
                g_e = details.get("end")

            # Extract full genomic gene sequence (exons + introns + UTRs) for canonical isoform
            if g_s >= start and g_e <= end:
                gene_genomic_seq = locus_seq[g_s - start:g_e - start]
            else:
                gene_genomic_seq = self.dm.fetch_region_sequence(species, chrom, g_s, g_e)

            full_gene_gc = self.calculate_gc(gene_genomic_seq)
            full_gene_len = len(gene_genomic_seq)

            coding_span_gc = None
            coding_span_len = 0

            if target_t:
                # Coding span boundaries (start codon to stop codon, excluding 5' and 3' UTRs)
                cds_s = None
                cds_e = None
                if "CDS" in target_t and target_t["CDS"]:
                    sorted_cds_coords = sorted(target_t["CDS"], key=lambda c: c.get("start", 0))
                    cds_s = sorted_cds_coords[0].get("start", 0)
                    cds_e = sorted_cds_coords[-1].get("end", 0)
                elif "Translation" in target_t and target_t["Translation"]:
                    cds_s = target_t["Translation"].get("start", 0)
                    cds_e = target_t["Translation"].get("end", 0)

                if cds_s and cds_e and cds_e >= cds_s:
                    coding_span_len = cds_e - cds_s + 1
                    c_s_rel = cds_s - 1 - g_s
                    c_e_rel = cds_e - g_s
                    if 0 <= c_s_rel < len(gene_genomic_seq) and c_e_rel <= len(gene_genomic_seq):
                        coding_span_seq = gene_genomic_seq[c_s_rel:c_e_rel]
                    else:
                        coding_span_seq = self.dm.fetch_region_sequence(species, chrom, cds_s - 1, cds_e)
                    coding_span_gc = self.calculate_gc(coding_span_seq)
                else:
                    coding_span_gc = full_gene_gc
                    coding_span_len = full_gene_len

                # 1. CDS calculation (for protein coding)
                if target_t.get("biotype") == "protein_coding":
                    t_id = target_t.get("id")
                    cds_seq = self.dm.fetch_transcript_cds(t_id)

                    # If not fetched via API, splice CDS from genomic sequence using CDS feature coordinates
                    if not cds_seq and "CDS" in target_t and target_t["CDS"]:
                        sorted_cds = sorted(target_t["CDS"], key=lambda c: c.get("start", 0))
                        cds_parts = []
                        for c in sorted_cds:
                            c_s = c.get("start", 0) - 1 - g_s
                            c_e = c.get("end", 0) - g_s
                            if 0 <= c_s < len(gene_genomic_seq):
                                cds_parts.append(gene_genomic_seq[max(0, c_s):min(len(gene_genomic_seq), c_e)])
                        raw_cds = "".join(cds_parts)
                        if raw_cds:
                            if target_t.get("strand") == -1:
                                cds_seq = self.reverse_complement(raw_cds)
                            else:
                                cds_seq = raw_cds

                    if cds_seq:
                        cds_gc = self.calculate_gc(cds_seq)
                        gc3 = self.calculate_gc3(cds_seq)
                        cds_len = len(cds_seq)
                        all_cds_gc.append(cds_gc)
                        if gc3 is not None:
                            all_gc3.append(gc3)

                        if g["is_globin"] == "True":
                            globin_cds_gc.append(cds_gc)
                            if gc3 is not None:
                                globin_gc3.append(gc3)

                # 2. Exon and Intron calculation from target transcript
                exons = target_t.get("Exon", [])
                if exons:
                    sorted_exons = sorted(exons, key=lambda e: e.get("start", 0))
                    ex_seq_parts = []
                    intron_seq_parts = []

                    for i, ex in enumerate(sorted_exons):
                        e_s = ex.get("start", 0) - 1 - g_s
                        e_e = ex.get("end", 0) - g_s
                        if 0 <= e_s < len(gene_genomic_seq):
                            ex_seq_parts.append(gene_genomic_seq[max(0, e_s):min(len(gene_genomic_seq), e_e)])

                        # Intron between current and next exon
                        if i < len(sorted_exons) - 1:
                            next_s = sorted_exons[i + 1].get("start", 0) - 1 - g_s
                            if 0 <= e_e < len(gene_genomic_seq) and e_e < next_s:
                                intron_seq_parts.append(gene_genomic_seq[e_e:min(len(gene_genomic_seq), next_s)])

                    exonic_seq = "".join(ex_seq_parts)
                    if exonic_seq:
                        exonic_gc = self.calculate_gc(exonic_seq)
                        exonic_len = len(exonic_seq)

                    intronic_seq = "".join(intron_seq_parts)
                    if intronic_seq:
                        intronic_gc = self.calculate_gc(intronic_seq)
                        intronic_len = len(intronic_seq)
                    elif len(sorted_exons) <= 1:
                        intronic_gc = None
                        intronic_len = 0

            # Mark all exons on locus mask for overall locus partitioning
            if details and "Transcript" in details:
                for t in details.get("Transcript", []):
                    for ex in t.get("Exon", []):
                        ex_s = ex.get("start", 0) - 1 - start
                        ex_e = ex.get("end", 0) - start
                        if 0 <= ex_s < locus_length and 0 < ex_e <= locus_length:
                            locus_mask[max(0, ex_s):min(locus_length, ex_e)] = 1

            if coding_span_gc is None:
                coding_span_gc = full_gene_gc
                coding_span_len = full_gene_len

            if cds_gc is None:
                import os
                baseline_file = os.path.join("data", "genome_baselines", f"{species.lower()}_genes_gc.csv.gz")
                if not os.path.exists(baseline_file):
                    baseline_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "genome_baselines", f"{species.lower()}_genes_gc.csv.gz")
                if os.path.exists(baseline_file):
                    if not hasattr(self, "_baseline_cache"):
                        self._baseline_cache = {}
                    if species not in self._baseline_cache:
                        import pandas as pd
                        self._baseline_cache[species] = pd.read_csv(baseline_file).set_index("gene_id")
                    b_df = self._baseline_cache[species]
                    if g_id in b_df.index:
                        import pandas as pd
                        b_row = b_df.loc[g_id]
                        if isinstance(b_row, pd.DataFrame):
                            b_row = b_row.iloc[0]
                        if "cds_gc" in b_row and pd.notna(b_row["cds_gc"]):
                            cds_gc = float(b_row["cds_gc"])
                        if "gc3" in b_row and pd.notna(b_row["gc3"]):
                            gc3 = float(b_row["gc3"])
                        if "cds_length" in b_row and pd.notna(b_row["cds_length"]):
                            cds_len = int(b_row["cds_length"])
                        if (full_gene_gc is None or full_gene_gc == 0.0) and "whole_gene_gc" in b_row and pd.notna(b_row["whole_gene_gc"]):
                            full_gene_gc = float(b_row["whole_gene_gc"])
                            full_gene_len = int(b_row.get("whole_gene_len", cds_len) if hasattr(b_row, "get") else cds_len)
                        if (coding_span_gc is None or coding_span_gc == 0.0):
                            coding_span_gc = full_gene_gc
                            coding_span_len = full_gene_len

                        if cds_gc is not None:
                            all_cds_gc.append(cds_gc)
                            if gc3 is not None:
                                all_gc3.append(gc3)
                            if g["is_globin"] == "True":
                                globin_cds_gc.append(cds_gc)
                                if gc3 is not None:
                                    globin_gc3.append(gc3)

            log10_len = round(float(np.log10(coding_span_len)), 4) if coding_span_len > 0 else None

            gene_gc_records.append({
                "species": species,
                "clade": clade,
                "locus_id": locus_id,
                "gene_id": g_id,
                "gene_name": g["gene_name"],
                "display_symbol": g["display_symbol"],
                "category": g["category"],
                "subtype": g["subtype"],
                "is_globin": g["is_globin"],
                "start": g_s,
                "end": g_e,
                "coding_span_gc": coding_span_gc,
                "coding_span_len": coding_span_len,
                "log10_coding_span_len": log10_len,
                "full_gene_gc": full_gene_gc,
                "full_gene_len": full_gene_len,
                "exonic_gc": exonic_gc,
                "exonic_len": exonic_len,
                "intronic_gc": intronic_gc,
                "intronic_len": intronic_len,
                "cds_gc": cds_gc,
                "cds_len": cds_len,
                "gc3": gc3,
            })

        # 4. Calculate Exonic vs Intergenic GC
        exonic_seq = "".join([locus_seq[i] for i in range(locus_length) if locus_mask[i] == 1])
        intergenic_seq = "".join([locus_seq[i] for i in range(locus_length) if locus_mask[i] == 0])

        exonic_gc = self.calculate_gc(exonic_seq)
        intergenic_gc = self.calculate_gc(intergenic_seq)

        return {
            "species": species,
            "common_name": interval.get("common_name", ""),
            "clade": clade,
            "locus_id": locus_id,
            "chromosome": chrom,
            "start": start,
            "end": end,
            "length_bp": locus_length,
            "overall_gc": overall_locus_gc,
            "exonic_gc": exonic_gc,
            "intergenic_gc": intergenic_gc,
            "mean_cds_gc": round(float(np.mean(all_cds_gc)), 2) if all_cds_gc else None,
            "mean_gc3": round(float(np.mean(all_gc3)), 2) if all_gc3 else None,
            "globin_mean_cds_gc": round(float(np.mean(globin_cds_gc)), 2) if globin_cds_gc else None,
            "globin_mean_gc3": round(float(np.mean(globin_gc3)), 2) if globin_gc3 else None,
            "sw_10kb_mean": sw_10kb["mean_gc"],
            "sw_10kb_min": sw_10kb["min_gc"],
            "sw_10kb_max": sw_10kb["max_gc"],
            "sw_100kb_mean": sw_100kb["mean_gc"],
            "sw_100kb_min": sw_100kb["min_gc"],
            "sw_100kb_max": sw_100kb["max_gc"],
            "gene_records": gene_gc_records,
            "num_globins": len([g for g in annotated_genes if g["is_globin"] == "True"]),
            "num_total_genes": len(annotated_genes),
        }
