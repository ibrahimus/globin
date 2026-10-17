"""Curator module to reconcile and standardize vertebrate globin locus intervals."""

import csv
import os
from typing import Dict, List, Any, Optional


class LocusCurator:
    """Standardizes, corrects, and curates genomic intervals for comparative analysis."""

    CLADE_MAP = {
        "Homo_sapiens": "amniote_mammal",
        "Mus_musculus": "amniote_mammal",
        "Bos_taurus": "amniote_mammal",
        "Sus_scrofa": "amniote_mammal",
        "Monodelphis_domestica": "amniote_mammal",
        "Erinaceus_europaeus": "amniote_mammal",
        "Gallus_gallus": "amniote_bird",
        "Chrysemys_picta_bellii": "amniote_reptile",
        "xenopus_tropicalis": "amphibian",
        "danio_rerio": "teleost_fish",
        "oryzias_latipes": "teleost_fish",
        "takifugu_rubripes": "teleost_fish",
        "petromyzon_marinus": "cyclostome",
    }

    CURATED_RULES = {
        # Human reference (omitted in Mike's file)
        "Homo_sapiens": {
            "locus_id": "alpha_globin",
            "assembly_name": "GRCh38",
            "Chromosome": "16",
            "Start": 60000,
            "End": 275000,
            "Strand": "-",
            "left_gene_name": "RHBDF1",
            "right_gene_name": "LUC7L",
            "curation_note": "Human amniote reference locus (GRCh38)",
            "status": "curated_reference",
        },
        # Xenopus correction (Mike's interval missed all 4 alpha globins)
        "xenopus_tropicalis": {
            "locus_id": "alpha_beta_cluster",
            "assembly_name": "UCB_Xtro_10.0",
            "Chromosome": "9",
            "Start": 40315000,
            "End": 41265000,
            "Strand": "+",
            "left_gene_name": "snrnp25/hba1",
            "right_gene_name": "fam234a",
            "curation_note": "Corrected Xenopus inversion: spans from alpha globins (hba1/hbz) through rhbdf1 to fam234a",
            "status": "curated_correction",
        },
        # Western painted turtle recovery
        "Chrysemys_picta_bellii": {
            "locus_id": "alpha_globin",
            "assembly_name": "Chrysemys_picta_bellii-3.0.3",
            "Chromosome": "KK088206.1",
            "Start": 1415000,
            "End": 1865000,
            "Strand": "+",
            "left_gene_name": "ENSCPBG00000026394/HBA",
            "right_gene_name": "RHBDF1",
            "curation_note": "Recovered turtle locus: spans from HBA/HBAD/HBZ globins through NPRL3/MPG to RHBDF1",
            "status": "curated_recovered",
        },
        # Hedgehog annotation warning
        "Erinaceus_europaeus": {
            "locus_id": "alpha_flanking_only",
            "assembly_name": "HEDGEHOG",
            "Chromosome": "GeneScaffold_103",
            "Start": 0,
            "End": 129753,
            "Strand": "-",
            "left_gene_name": "RHBDF1",
            "right_gene_name": "Fixed endpoint",
            "curation_note": "Fragmented assembly: GeneScaffold_103 contains NPRL3/RHBDF1 but alpha globins reside on scaffold_235869",
            "status": "fragmented_assembly_no_globins",
        },
        # Teleost Dual Loci (Locus A = Major, Locus B = Minor 3R paralog)
        "danio_rerio_locus_b": {
            "target_species": "danio_rerio",
            "species_common_name": "zebrafish",
            "locus_id": "alpha_beta_locus_B",
            "assembly_name": "GRCz11",
            "Chromosome": "12",
            "Start": 20222000,
            "End": 20415000,
            "Strand": "-",
            "left_gene_name": "rhbdf1b",
            "right_gene_name": "lcmt1",
            "curation_note": "Teleost 3R paralogous minor cluster: rhbdf1b - hbae4/hbbe3/hbae5/hbbe2 - lcmt1",
            "status": "curated_teleost_3R_locus_b",
        },
        "oryzias_latipes_locus_b": {
            "target_species": "oryzias_latipes",
            "species_common_name": "Japanese medaka HdrR",
            "locus_id": "alpha_beta_locus_B",
            "assembly_name": "ASM223467v1",
            "Chromosome": "19",
            "Start": 2970000,
            "End": 3050000,
            "Strand": "-",
            "left_gene_name": "rhbdf1b",
            "right_gene_name": "lcmt1",
            "curation_note": "Teleost 3R paralogous minor cluster: rhbdf1b - hbae4 - lcmt1",
            "status": "curated_teleost_3R_locus_b",
        },
        "takifugu_rubripes_locus_b": {
            "target_species": "takifugu_rubripes",
            "species_common_name": "fugu",
            "locus_id": "alpha_beta_locus_B",
            "assembly_name": "fTakRub1.2",
            "Chromosome": "1",
            "Start": 10330000,
            "End": 10375000,
            "Strand": "-",
            "left_gene_name": "rhbdf1b",
            "right_gene_name": "lcmt1",
            "curation_note": "Teleost 3R paralogous minor cluster: rhbdf1b - hbae4 - lcmt1",
            "status": "curated_teleost_3R_locus_b",
        },
    }

    BETA_GLOBIN_RULES = {
        "Homo_sapiens": {
            "species": "Homo_sapiens",
            "common_name": "human",
            "clade": "amniote_mammal",
            "locus_id": "beta_globin",
            "assembly": "GRCh38",
            "chromosome": "11",
            "start": 5200000,
            "end": 5340000,
            "strand": "-",
            "left_gene": "HBE1",
            "right_gene": "HBB",
            "status": "curated_beta_locus",
            "curation_note": "Amniote beta-globin locus (chr11): HBE1 - HBG2 - HBG1 - HBD - HBB flanked by OR51 and CAPN1",
        },
        "Mus_musculus": {
            "species": "Mus_musculus",
            "common_name": "mouse",
            "clade": "amniote_mammal",
            "locus_id": "beta_globin",
            "assembly": "GRCm39",
            "chromosome": "7",
            "start": 103420000,
            "end": 103570000,
            "strand": "+",
            "left_gene": "Hbb-y",
            "right_gene": "Hbb-bs",
            "status": "curated_beta_locus",
            "curation_note": "Amniote beta-globin locus (chr7): Hbb-y - Hbb-bh1 - Hbb-bt - Hbb-bs",
        },
        "Bos_taurus": {
            "species": "Bos_taurus",
            "common_name": "cattle",
            "clade": "amniote_mammal",
            "locus_id": "beta_globin",
            "assembly": "ARS-UCD2.0",
            "chromosome": "15",
            "start": 48320000,
            "end": 48460000,
            "strand": "+",
            "left_gene": "HBE1",
            "right_gene": "HBB",
            "status": "curated_beta_locus",
            "curation_note": "Amniote beta-globin locus (chr15): HBE1 - HBE4 - HBB - HBE2",
        },
        "Sus_scrofa": {
            "species": "Sus_scrofa",
            "common_name": "pig",
            "clade": "amniote_mammal",
            "locus_id": "beta_globin",
            "assembly": "Sscrofa11.1",
            "chromosome": "9",
            "start": 4750000,
            "end": 4880000,
            "strand": "-",
            "left_gene": "HBE1",
            "right_gene": "HBB",
            "status": "curated_beta_locus",
            "curation_note": "Amniote beta-globin locus (chr9): HBE1 - HBB flanked by OR cluster",
        },
        "Monodelphis_domestica": {
            "species": "Monodelphis_domestica",
            "common_name": "gray short-tailed opossum",
            "clade": "amniote_mammal",
            "locus_id": "beta_globin",
            "assembly": "ASM229v1",
            "chromosome": "4",
            "start": 352150000,
            "end": 352320000,
            "strand": "+",
            "left_gene": "HBE1",
            "right_gene": "HBE1_dup",
            "status": "curated_beta_locus",
            "curation_note": "Marsupial beta-globin locus (chr4): HBE1 cluster",
        },
        "Gallus_gallus": {
            "species": "Gallus_gallus",
            "common_name": "chicken",
            "clade": "amniote_bird",
            "locus_id": "beta_globin",
            "assembly": "bGalGal1.mat.broiler.GRCg7b",
            "chromosome": "1",
            "start": 195940000,
            "end": 196060000,
            "strand": "+",
            "left_gene": "HBE1",
            "right_gene": "HBBA",
            "status": "curated_beta_locus",
            "curation_note": "Avian beta-globin locus (chr1): HBE1 - HBBA",
        },
        "Chrysemys_picta_bellii": {
            "species": "Chrysemys_picta_bellii",
            "common_name": "Western painted turtle",
            "clade": "amniote_reptile",
            "locus_id": "beta_globin",
            "assembly": "Chrysemys_picta_bellii-3.0.3",
            "chromosome": "KK083065.1",
            "start": 720000,
            "end": 800000,
            "strand": "+",
            "left_gene": "ENSCPBG00000021825",
            "right_gene": "ENSCPBG00000021826",
            "status": "curated_beta_locus",
            "curation_note": "Reptilian beta-globin locus (KK083065.1): beta-globin homologs syntenic with LUC7L2",
        },
    }

    @classmethod
    def load_from_manifest(
        cls, manifest_path: str = "config/species_manifest.yaml", status_filter: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """Loads curated locus intervals directly from the central species manifest."""
        import yaml
        if not os.path.exists(manifest_path):
            alt_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", manifest_path)
            if os.path.exists(alt_path):
                manifest_path = alt_path
            elif os.path.exists(os.path.join("..", manifest_path)):
                manifest_path = os.path.join("..", manifest_path)
            else:
                return []
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        species_dict = data.get("species", {})
        filter_set = set(status_filter) if status_filter else {"active"}

        curated = []
        for sp, sp_info in species_dict.items():
            status = sp_info.get("status", "active")
            if status not in filter_set or status == "excluded":
                continue
            common_name = sp_info.get("common_name", sp)
            clade = sp_info.get("clade", "unknown")
            assembly = sp_info.get("assembly", "")

            for loc in sp_info.get("loci", []):
                start = int(loc.get("start", 0))
                end = int(loc.get("end", 0))
                curated.append({
                    "species": sp,
                    "common_name": common_name,
                    "clade": clade,
                    "locus_id": loc.get("locus_id", "alpha_globin"),
                    "assembly": assembly,
                    "chromosome": str(loc.get("chromosome", "")),
                    "start": start,
                    "end": end,
                    "span_bp": end - start,
                    "strand": loc.get("strand", "+"),
                    "left_gene": loc.get("left_gene", ""),
                    "right_gene": loc.get("right_gene", ""),
                    "status": f"manifest_{status}",
                    "curation_note": loc.get("notes", ""),
                })
        return curated

    def __init__(self, raw_intervals: Optional[List[Dict[str, Any]]] = None):
        self.raw_intervals = raw_intervals or []

    def curate(self, manifest_path: Optional[str] = "config/species_manifest.yaml") -> List[Dict[str, Any]]:
        """Processes intervals from manifest if available, else raw intervals."""
        if manifest_path and os.path.exists(manifest_path):
            curated_manifest = self.load_from_manifest(manifest_path, status_filter=["active"])
            if curated_manifest:
                return curated_manifest

        curated: List[Dict[str, Any]] = []

        # 1. Add human reference first
        human_rule = self.CURATED_RULES["Homo_sapiens"]
        curated.append({
            "species": "Homo_sapiens",
            "common_name": "human",
            "clade": self.CLADE_MAP["Homo_sapiens"],
            "locus_id": human_rule["locus_id"],
            "assembly": human_rule["assembly_name"],
            "chromosome": human_rule["Chromosome"],
            "start": human_rule["Start"],
            "end": human_rule["End"],
            "span_bp": human_rule["End"] - human_rule["Start"],
            "strand": human_rule["Strand"],
            "left_gene": human_rule["left_gene_name"],
            "right_gene": human_rule["right_gene_name"],
            "status": human_rule["status"],
            "curation_note": human_rule["curation_note"],
        })

        seen_species = set()
        for raw in self.raw_intervals:
            sp = raw.get("target_species", "")
            seen_species.add(sp)

            # Check if this species has an explicit curated override
            if sp in self.CURATED_RULES:
                rule = self.CURATED_RULES[sp]
                curated.append({
                    "species": sp,
                    "common_name": raw.get("species_common_name", ""),
                    "clade": self.CLADE_MAP.get(sp, "unknown"),
                    "locus_id": rule.get("locus_id", "alpha_globin"),
                    "assembly": rule.get("assembly_name", raw.get("assembly_name", "")),
                    "chromosome": rule.get("Chromosome", str(raw.get("Chromosome", ""))),
                    "start": int(rule.get("Start", raw.get("Start", 0))),
                    "end": int(rule.get("End", raw.get("End", 0))),
                    "span_bp": int(rule.get("End", 0)) - int(rule.get("Start", 0)),
                    "strand": rule.get("Strand", raw.get("Strand", "")),
                    "left_gene": rule.get("left_gene_name", raw.get("left_gene_name", "")),
                    "right_gene": rule.get("right_gene_name", raw.get("right_gene_name", "")),
                    "status": rule.get("status", "curated_correction"),
                    "curation_note": rule.get("curation_note", ""),
                })
            else:
                # Retain raw with clade mapping
                locus_id = "alpha_beta_locus_A" if sp in ["danio_rerio", "oryzias_latipes", "takifugu_rubripes"] else "alpha_globin"
                start = int(raw.get("Start", 0))
                end = int(raw.get("End", 0))
                curated.append({
                    "species": sp,
                    "common_name": raw.get("species_common_name", ""),
                    "clade": self.CLADE_MAP.get(sp, "unknown"),
                    "locus_id": locus_id,
                    "assembly": raw.get("assembly_name", ""),
                    "chromosome": str(raw.get("Chromosome", "")),
                    "start": start,
                    "end": end,
                    "span_bp": end - start,
                    "strand": raw.get("Strand", ""),
                    "left_gene": raw.get("left_gene_name", ""),
                    "right_gene": raw.get("right_gene_name", ""),
                    "status": "valid_as_provided",
                    "curation_note": "Direct from Mike's batch",
                })

            # If teleost, add Locus B
            locus_b_key = f"{sp}_locus_b"
            if locus_b_key in self.CURATED_RULES:
                rule_b = self.CURATED_RULES[locus_b_key]
                curated.append({
                    "species": sp,
                    "common_name": raw.get("species_common_name", ""),
                    "clade": self.CLADE_MAP.get(sp, "unknown"),
                    "locus_id": rule_b["locus_id"],
                    "assembly": rule_b["assembly_name"],
                    "chromosome": rule_b["Chromosome"],
                    "start": rule_b["Start"],
                    "end": rule_b["End"],
                    "span_bp": rule_b["End"] - rule_b["Start"],
                    "strand": rule_b["Strand"],
                    "left_gene": rule_b["left_gene_name"],
                    "right_gene": rule_b["right_gene_name"],
                    "status": rule_b["status"],
                    "curation_note": rule_b["curation_note"],
                })

        # Check if turtle was completely omitted from raw_intervals
        if "Chrysemys_picta_bellii" not in seen_species:
            turtle_rule = self.CURATED_RULES["Chrysemys_picta_bellii"]
            curated.append({
                "species": "Chrysemys_picta_bellii",
                "common_name": "Western painted turtle",
                "clade": self.CLADE_MAP["Chrysemys_picta_bellii"],
                "locus_id": turtle_rule["locus_id"],
                "assembly": turtle_rule["assembly_name"],
                "chromosome": turtle_rule["Chromosome"],
                "start": turtle_rule["Start"],
                "end": turtle_rule["End"],
                "span_bp": turtle_rule["End"] - turtle_rule["Start"],
                "strand": turtle_rule["Strand"],
                "left_gene": turtle_rule["left_gene_name"],
                "right_gene": turtle_rule["right_gene_name"],
                "status": turtle_rule["status"],
                "curation_note": turtle_rule["curation_note"],
            })

        # Append dedicated amniote beta globin loci
        for sp_key, beta_rule in self.BETA_GLOBIN_RULES.items():
            curated.append({
                "species": beta_rule["species"],
                "common_name": beta_rule["common_name"],
                "clade": beta_rule["clade"],
                "locus_id": beta_rule["locus_id"],
                "assembly": beta_rule["assembly"],
                "chromosome": beta_rule["chromosome"],
                "start": beta_rule["start"],
                "end": beta_rule["end"],
                "span_bp": beta_rule["end"] - beta_rule["start"],
                "strand": beta_rule["strand"],
                "left_gene": beta_rule["left_gene"],
                "right_gene": beta_rule["right_gene"],
                "status": beta_rule["status"],
                "curation_note": beta_rule["curation_note"],
            })

        return curated

    @staticmethod
    def save_csv(intervals: List[Dict[str, Any]], output_path: str):
        """Saves curated intervals to CSV."""
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        if not intervals:
            return
        fieldnames = list(intervals[0].keys())
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(intervals)

    @staticmethod
    def save_bed(intervals: List[Dict[str, Any]], output_path: str):
        """Saves curated intervals to standard 0-based BED format."""
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            for row in intervals:
                chrom = row["chromosome"]
                start = row["start"]
                end = row["end"]
                name = f"{row['species']}:{row['locus_id']}"
                score = "0"
                strand = row["strand"] if row["strand"] in ["+", "-"] else "."
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
