"""Unit tests for globin_pipeline components."""

import pytest
from globin_pipeline.parser import BatchParser
from globin_pipeline.curator import LocusCurator
from globin_pipeline.annotator import GeneAnnotator
from globin_pipeline.gc_analyzer import GCAnalyzer


def test_batch_parser():
    parser = BatchParser("1_data_from_mike/20260930")
    audit_res = parser.audit()
    assert audit_res["query_species"] == "Homo_sapiens"
    assert audit_res["target_species_count"] == 12
    assert "Chrysemys_picta_bellii" in audit_res["missing_in_intervals"]

    # Verify Xenopus error detection
    xenopus_audit = [a for a in audit_res["interval_audits"] if a["species"] == "xenopus_tropicalis"][0]
    assert any("OUTSIDE the interval" in issue for issue in xenopus_audit["issues"])

    # Verify Teleost chromosome mismatch detection
    danio_audit = [a for a in audit_res["interval_audits"] if a["species"] == "danio_rerio"][0]
    assert any("Chromosome mismatch" in issue for issue in danio_audit["issues"])


def test_locus_curator():
    parser = BatchParser("1_data_from_mike/20260930")
    curator = LocusCurator(parser.load_intervals())
    curated = curator.curate()

    # Total intervals should be 16 (1 human + 10 raw + 3 teleost locus B + 1 turtle + 1 hedgehog)
    species_list = [c["species"] for c in curated]
    assert "Homo_sapiens" in species_list
    assert "Chrysemys_picta_bellii" in species_list
    assert species_list.count("danio_rerio") == 2
    assert species_list.count("oryzias_latipes") == 2
    assert species_list.count("takifugu_rubripes") == 2

    # Check Xenopus corrected coordinates
    xen_cur = [c for c in curated if c["species"] == "xenopus_tropicalis"][0]
    assert xen_cur["start"] == 40315000
    assert xen_cur["end"] == 41265000


def test_gene_classification():
    annotator = GeneAnnotator()

    # Alpha globin test
    g_hba2 = {"id": "ENSG00000188536", "external_name": "HBA2", "biotype": "protein_coding", "description": "hemoglobin subunit alpha 2"}
    cls_hba2 = annotator.classify_gene(g_hba2)
    assert cls_hba2["is_globin"] == "True"
    assert cls_hba2["category"] == "alpha_globin"

    # Beta globin test
    g_hbe = {"id": "ENSXETG00000045178", "external_name": "hbe1.2", "biotype": "protein_coding", "description": "hemoglobin subunit gamma 2"}
    cls_hbe = annotator.classify_gene(g_hbe)
    assert cls_hbe["is_globin"] == "True"
    assert cls_hbe["category"] == "beta_globin"

    # Flanking / cluster non-globin marker test
    g_rhb = {"id": "ENSG00000140798", "external_name": "RHBDF1", "biotype": "protein_coding", "description": "rhomboid 5 homolog 1"}
    cls_rhb = annotator.classify_gene(g_rhb)
    assert cls_rhb["is_globin"] == "False"
    assert cls_rhb["category"] == "cluster_non_globin"

    # Lamprey unannotated globin test
    g_lamp = {"id": "ENSPMAG00000001587", "external_name": "", "biotype": "protein_coding", "description": ""}
    cls_lamp = annotator.classify_gene(g_lamp)
    assert cls_lamp["is_globin"] == "True"
    assert cls_lamp["category"] == "other_globin"


def test_gc_calculations():
    # Test GC percentage calculation
    seq = "ATGC"
    assert GCAnalyzer.calculate_gc(seq) == 50.0

    seq_high_gc = "GGCC"
    assert GCAnalyzer.calculate_gc(seq_high_gc) == 100.0

    seq_low_gc = "AATT"
    assert GCAnalyzer.calculate_gc(seq_low_gc) == 0.0

    # Test GC3 calculation
    # Codons: ATG (G=3rd), GCC (C=3rd), TTT (T=3rd), CGA (A=3rd)
    # Third bases: G, C, T, A -> 2/4 = 50.0%
    cds = "ATGGCC TTTCGA".replace(" ", "")
    assert GCAnalyzer.calculate_gc3(cds) == 50.0

    # Codons with 100% GC3: ATG (G), GCG (G), CCC (C) -> G, G, C -> 100.0%
    cds_high_gc3 = "ATGGCGCCG"
    assert GCAnalyzer.calculate_gc3(cds_high_gc3) == 100.0
