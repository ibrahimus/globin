"""Command-line interface for globin-pipeline."""

import os
import sys
import click
import numpy as np
import pandas as pd
from globin_pipeline.parser import BatchParser
from globin_pipeline.curator import LocusCurator
from globin_pipeline.data_manager import DataManager
from globin_pipeline.annotator import GeneAnnotator
from globin_pipeline.gc_analyzer import GCAnalyzer
from globin_pipeline.baseline_manager import BaselineManager
from globin_pipeline.zscore_analyzer import ZScoreAnalyzer
from globin_pipeline.visualizer import Visualizer
from globin_pipeline.reporter import Reporter
from globin_pipeline import __version__


@click.group()
@click.version_option(version=__version__)
def main():
    """Globin Pipeline: Comparative genomics analysis for vertebrate globin loci."""
    pass


@main.command()
@click.argument("batch_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--output", "-o", default=None, help="Path to write audit Markdown report.")
def audit(batch_dir, output):
    """Audits an incoming batch directory from collaborator."""
    click.echo(f"Auditing batch directory: {batch_dir}...")
    parser = BatchParser(batch_dir)
    res = parser.audit()
    report_md = parser.generate_markdown_report(res)

    if output:
        with open(output, "w", encoding="utf-8") as f:
            f.write(report_md)
        click.echo(f"Audit report written to: {output}")
    else:
        click.echo(report_md)


@main.command()
@click.argument("batch_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--output", "-o", default="data/curated_intervals.csv", help="Path to output curated CSV.")
@click.option("--bed", default=None, help="Optional path to output curated BED.")
def curate(batch_dir, output, bed):
    """Standardizes and curates locus intervals from raw collaborator files."""
    click.echo(f"Curating intervals from: {batch_dir}...")
    parser = BatchParser(batch_dir)
    raw_intervals = parser.load_intervals()
    curator = LocusCurator(raw_intervals)
    curated = curator.curate()

    curator.save_csv(curated, output)
    click.echo(f"Curated intervals saved to: {output} ({len(curated)} intervals)")

    if bed:
        curator.save_bed(curated, bed)
        click.echo(f"Curated BED saved to: {bed}")


@main.command()
@click.option("--cache-dir", default="data/genome_baselines", help="Directory to cache genome baselines.")
@click.option("--output", "-o", default="results/genome_baseline_stats.csv", help="Output summary CSV.")
def baseline(cache_dir, output):
    """Computes and caches whole-genome coding gene GC/GC3 baselines for all species."""
    click.echo("Computing whole-genome baselines (streaming Ensembl CDS)...")
    bm = BaselineManager(cache_dir=cache_dir)
    df_summary = bm.compute_all_baselines()
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    df_summary.to_csv(output, index=False)
    click.echo(f"Genome baseline summary saved to: {output}")


@main.command()
@click.option("--genes", "-g", default="results/genes_in_intervals.csv", help="Path to annotated genes CSV.")
@click.option("--gene-gc", default="results/gene_gc_summary.csv", help="Path to gene GC summary CSV.")
@click.option("--output", "-o", default="results/gene_zscores.csv", help="Path to output Z-scores CSV.")
@click.option("--stats-output", default="results/statistical_comparisons.csv", help="Output stats CSV.")
def zscore(genes, gene_gc, output, stats_output):
    """Calculates standardized Z-scores and comparative hypothesis tests."""
    click.echo("Calculating standardized Z-scores and relative ranks...")
    df_genes = pd.read_csv(genes).to_dict(orient="records")
    df_gc = pd.read_csv(gene_gc).to_dict(orient="records")

    bm = BaselineManager()
    za = ZScoreAnalyzer(bm)
    df_res, df_stats = za.analyze_genes(df_genes, df_gc)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    df_res.to_csv(output, index=False)
    df_stats.to_csv(stats_output, index=False)

    click.echo(f"Gene Z-scores saved to: {output}")
    click.echo(f"Statistical comparisons saved to: {stats_output}")


@main.command(name="generate-plots")
@click.option("--species", "-s", default=None, help="Target species (e.g. 'Homo_sapiens'). Defaults to all.")
@click.option("--zscores", "-z", default=None, help="Path to gene Z-scores CSV (defaults to results/gene_zscores.csv or bundled dataset).")
@click.option("--intervals", "-i", default=None, help="Path to curated intervals CSV (defaults to results/curated_intervals.csv or bundled dataset).")
@click.option("--output-dir", "-o", default="results", help="Base output directory for generated galleries.")
@click.option("--plot-type", "-t", type=click.Choice(["all", "diagnostic", "chromosome", "bivariate", "bivariate-inset", "meta"]), default="all", help="Type of plot to generate.")
def generate_plots(species, zscores, intervals, output_dir, plot_type):
    """Generates publication plot galleries: diagnostic_plots, chromosome_plots, bivariate_length_gc_plots, bivariate_with_inset_plots."""
    import importlib.resources as pkg_resources

    # 1. Resolve Z-scores
    if zscores and os.path.exists(zscores):
        df_genes = pd.read_csv(zscores)
    elif os.path.exists("results/gene_zscores.csv"):
        df_genes = pd.read_csv("results/gene_zscores.csv")
    elif os.path.exists("../results/gene_zscores.csv"):
        df_genes = pd.read_csv("../results/gene_zscores.csv")
    else:
        ref_file = pkg_resources.files("globin_pipeline.data") / "gene_zscores.csv"
        df_genes = pd.read_csv(ref_file)

    # 2. Resolve Curated Intervals
    if intervals and os.path.exists(intervals):
        df_intervals = pd.read_csv(intervals)
    elif os.path.exists("results/curated_intervals.csv"):
        df_intervals = pd.read_csv("results/curated_intervals.csv")
    elif os.path.exists("../results/curated_intervals.csv"):
        df_intervals = pd.read_csv("../results/curated_intervals.csv")
    else:
        ref_int = pkg_resources.files("globin_pipeline.data") / "curated_intervals.csv"
        df_intervals = pd.read_csv(ref_int)

    viz = Visualizer()
    diag_dir = os.path.join(output_dir, "diagnostic_plots")
    chr_dir = os.path.join(output_dir, "chromosome_plots")
    biv_dir = os.path.join(output_dir, "bivariate_length_gc_plots")
    inset_dir = os.path.join(output_dir, "bivariate_with_inset_plots")

    if plot_type in ["all", "diagnostic"]:
        os.makedirs(diag_dir, exist_ok=True)
    if plot_type in ["all", "chromosome"]:
        os.makedirs(chr_dir, exist_ok=True)
    if plot_type in ["all", "bivariate"]:
        os.makedirs(biv_dir, exist_ok=True)
    if plot_type in ["all", "bivariate-inset"]:
        os.makedirs(inset_dir, exist_ok=True)

    if species and species.lower() != "all":
        target_species = [s for s in df_genes["species"].unique() if s.lower() == species.lower()]
        if not target_species:
            raise click.ClickException(f"Species '{species}' not found in dataset ({len(df_genes['species'].unique())} available).")
    else:
        target_species = sorted(df_genes["species"].unique())

    click.echo(f"=== Generating {plot_type} plots for {len(target_species)} species in '{output_dir}' ===")
    for idx, sp in enumerate(target_species, 1):
        click.echo(f"  [{idx:02d}/{len(target_species)}] {sp}...")
        if plot_type in ["all", "diagnostic"]:
            viz.plot_diagnostic(sp, df_genes, output_dir=diag_dir)
            viz.plot_diagnostic_interactive(sp, df_genes, output_dir=diag_dir)
        if plot_type in ["all", "chromosome"]:
            viz.plot_chromosome_linear(sp, df_genes, df_intervals=df_intervals, output_dir=chr_dir)
        if plot_type in ["all", "bivariate"]:
            viz.plot_length_vs_gc(sp, df_genes, output_dir=biv_dir)
        if plot_type in ["all", "bivariate-inset"]:
            viz.plot_bivariate_with_chromosome_inset(sp, df_genes, df_intervals=df_intervals, output_dir=inset_dir)

    if plot_type in ["all", "meta"] and (species is None or species.lower() == "all"):
        click.echo("  Generating cross-species meta-graphs...")
        meta_png = os.path.join(output_dir, "meta_comparative_graph.png")
        clean_png = os.path.join(output_dir, "meta_comparative_graph_clean.png")
        viz.plot_meta_comparative(df_genes, output_path=meta_png, clean_mode=False)
        viz.plot_meta_interactive(df_genes, output_path=meta_png.replace(".png", ".html"), clean_mode=False)
        viz.plot_meta_comparative(df_genes, output_path=clean_png, clean_mode=True)
        viz.plot_meta_interactive(df_genes, output_path=clean_png.replace(".png", ".html"), clean_mode=True)

    click.echo(f"✓ All requested plots successfully generated in '{output_dir}/'!")


@main.command(name="plot-diagnostic")
@click.option("--species", "-s", default=None, help="Target species (e.g. 'Homo_sapiens'). Defaults to all.")
@click.option("--zscores", "-z", default=None, help="Path to gene Z-scores CSV.")
@click.option("--intervals", "-i", default=None, help="Path to curated intervals CSV.")
@click.option("--output-dir", "-o", default="results", help="Directory for plots.")
@click.pass_context
def plot_diagnostic(ctx, species, zscores, intervals, output_dir):
    """Alias for generate-plots."""
    ctx.invoke(generate_plots, species=species, zscores=zscores, intervals=intervals, output_dir=output_dir, plot_type="all")


@main.command(name="plot-meta")
@click.option("--zscores", "-z", default="results/gene_zscores.csv", help="Path to gene Z-scores CSV.")
@click.option("--output", "-o", default="results/meta_comparative_graph.png", help="Path for meta-graph.")
def plot_meta(zscores, output):
    """Generates cross-species meta-graph comparing all vertebrates (comprehensive and clean versions)."""
    click.echo(f"Generating cross-species meta-graph: {output}...")
    df_genes = pd.read_csv(zscores)
    viz = Visualizer()
    out = viz.plot_meta_comparative(df_genes, output_path=output, clean_mode=False)
    html_out = output.replace(".png", ".html")
    viz.plot_meta_interactive(df_genes, output_path=html_out, clean_mode=False)
    click.echo(f"  ✓ Saved (comprehensive): {out}, {out.replace('.png', '.pdf')}, and {html_out}")

    clean_output = output.replace(".png", "_clean.png")
    clean_out = viz.plot_meta_comparative(df_genes, output_path=clean_output, clean_mode=True)
    clean_html_out = clean_output.replace(".png", ".html")
    viz.plot_meta_interactive(df_genes, output_path=clean_html_out, clean_mode=True)
    click.echo(f"  ✓ Saved (clean): {clean_out}, {clean_out.replace('.png', '.pdf')}, and {clean_html_out}")


@main.command(name="run-all")
@click.argument("batch_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", "-o", default="results", help="Directory to store all results.")
def run_all(batch_dir, output_dir):
    """Runs complete end-to-end pipeline: audit -> curate -> analyze -> baseline -> zscore -> plot -> report."""
    click.echo(f"=== Starting Comprehensive Globin Analysis Pipeline ===")
    click.echo(f"Input batch: {batch_dir}")
    click.echo(f"Output directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # 1. Audit
    click.echo("\n[1/7] Auditing batch files...")
    parser = BatchParser(batch_dir)
    audit_res = parser.audit()
    audit_md = parser.generate_markdown_report(audit_res)

    # 2. Curate
    click.echo("[2/7] Curating genomic intervals...")
    raw_intervals = parser.load_intervals()
    curator = LocusCurator(raw_intervals)
    curated_intervals = curator.curate()
    curator.save_csv(curated_intervals, os.path.join(output_dir, "curated_intervals.csv"))
    curator.save_bed(curated_intervals, os.path.join(output_dir, "curated_intervals.bed"))

    # 3. Locus & Gene Analysis
    click.echo("[3/7] Extracting genes and analyzing locus GC...")
    dm = DataManager()
    annotator = GeneAnnotator(dm)
    analyzer = GCAnalyzer(dm)
    reporter = Reporter(output_dir)

    all_annotated_genes = []
    locus_gc_records = []
    all_gene_gc_records = []

    for idx, inter in enumerate(curated_intervals, 1):
        sp = inter["species"]
        loc_id = inter.get("locus_id", "alpha_globin")
        click.echo(f"  [{idx:02d}/{len(curated_intervals):02d}] {sp} ({inter.get('clade')}): {loc_id}...")

        genes = annotator.annotate_interval(inter)
        all_annotated_genes.extend(genes)

        gc_res = analyzer.analyze_locus(inter, genes)
        locus_gc_records.append(gc_res)
        all_gene_gc_records.extend(gc_res.get("gene_records", []))

    # Initial reporting
    reporter.export_all(
        audit_md=audit_md,
        curated_intervals=curated_intervals,
        annotated_genes=all_annotated_genes,
        locus_gc_records=locus_gc_records,
        gene_gc_records=all_gene_gc_records,
    )

    # 4. Genome Baselines
    click.echo("\n[4/7] Ensuring whole-genome CDS baselines are computed and cached...")
    bm = BaselineManager(cache_dir="data/genome_baselines")
    unique_species = set(i["species"] for i in curated_intervals)
    sp_asm_map = {i["species"]: i["assembly"] for i in curated_intervals}
    df_baselines = bm.compute_all_baselines(sp_asm_map)
    df_baselines.to_csv(os.path.join(output_dir, "genome_baseline_stats.csv"), index=False)

    # 5. Z-Scores and Statistical Testing
    click.echo("\n[5/7] Computing standardized Z-scores and hypothesis tests...")
    za = ZScoreAnalyzer(bm)
    df_zscores, df_comp_stats = za.analyze_genes(all_annotated_genes, all_gene_gc_records)
    df_zscores.to_csv(os.path.join(output_dir, "gene_zscores.csv"), index=False)
    df_comp_stats.to_csv(os.path.join(output_dir, "statistical_comparisons.csv"), index=False)

    # 6. Per-species Diagnostic Plots
    diag_dir = os.path.join(output_dir, "diagnostic_plots")
    click.echo(f"\n[6/7] Generating 3-panel diagnostic figures in: {diag_dir}...")
    viz = Visualizer(bm)
    for sp in sorted(unique_species):
        click.echo(f"  Rendering diagnostic plot (PNG & HTML) for {sp}...")
        viz.plot_diagnostic(sp, df_zscores, output_dir=diag_dir)
        viz.plot_diagnostic_interactive(sp, df_zscores, output_dir=diag_dir)

    # 7. Cross-Species Meta-Graph
    meta_png = os.path.join(output_dir, "meta_comparative_graph.png")
    meta_html = os.path.join(output_dir, "meta_comparative_graph.html")
    meta_clean_png = os.path.join(output_dir, "meta_comparative_graph_clean.png")
    meta_clean_html = os.path.join(output_dir, "meta_comparative_graph_clean.html")
    click.echo(f"\n[7/7] Generating cross-species meta-graphs (comprehensive & clean): {meta_png}...")
    viz.plot_meta_comparative(df_zscores, output_path=meta_png, clean_mode=False)
    viz.plot_meta_interactive(df_zscores, output_path=meta_html, clean_mode=False)
    viz.plot_meta_comparative(df_zscores, output_path=meta_clean_png, clean_mode=True)
    viz.plot_meta_interactive(df_zscores, output_path=meta_clean_html, clean_mode=True)

@main.command(name="run-manifest")
@click.option("--manifest", default="config/species_manifest.yaml", help="Path to species manifest.")
@click.option("--output-dir", "-o", default="results", help="Directory to store all results.")
def run_manifest(manifest, output_dir):
    """Runs end-to-end comparative analysis from declarative species manifest."""
    click.echo("=== Starting Manifest-Driven Vertebrate Globin Pipeline ===")
    if not os.path.exists(manifest) and os.path.exists(os.path.join("..", manifest)):
        manifest = os.path.join("..", manifest)
    if output_dir == "results" and os.path.exists("../results"):
        output_dir = "../results"
    click.echo(f"Manifest: {manifest}")
    click.echo(f"Output directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load intervals from manifest
    click.echo("\n[1/6] Loading curated intervals from manifest...")
    curated_intervals = LocusCurator.load_from_manifest(manifest, status_filter=["active"])
    # Exclude hedgehog
    curated_intervals = [i for i in curated_intervals if i["species"] != "Erinaceus_europaeus"]
    LocusCurator.save_csv(curated_intervals, os.path.join(output_dir, "curated_intervals.csv"))
    LocusCurator.save_bed(curated_intervals, os.path.join(output_dir, "curated_intervals.bed"))
    click.echo(f"  Curated intervals: {len(curated_intervals)} across {len(set(i['species'] for i in curated_intervals))} species")

    # 2. Extract genes and analyze locus GC
    click.echo("\n[2/6] Extracting genes and analyzing locus GC...")
    dm = DataManager()
    annotator = GeneAnnotator(dm)
    analyzer = GCAnalyzer(dm)
    reporter = Reporter(output_dir)

    all_annotated_genes = []
    locus_gc_records = []
    all_gene_gc_records = []

    for idx, inter in enumerate(curated_intervals, 1):
        sp = inter["species"]
        loc_id = inter.get("locus_id", "alpha_globin")
        click.echo(f"  [{idx:02d}/{len(curated_intervals):02d}] {sp} ({inter.get('clade')}): {loc_id}...")

        genes = annotator.annotate_interval(inter)
        all_annotated_genes.extend(genes)

        gc_res = analyzer.analyze_locus(inter, genes)
        locus_gc_records.append(gc_res)
        all_gene_gc_records.extend(gc_res.get("gene_records", []))

    # Export intermediate results
    reporter.export_all(
        audit_md="",
        curated_intervals=curated_intervals,
        annotated_genes=all_annotated_genes,
        locus_gc_records=locus_gc_records,
        gene_gc_records=all_gene_gc_records,
    )

    # 3. Whole-genome baselines
    click.echo("\n[3/6] Ensuring whole-genome baselines are computed and cached...")
    bm = BaselineManager(cache_dir="data/genome_baselines")
    unique_species = sorted(set(i["species"] for i in curated_intervals))
    sp_asm_map = {i["species"]: i["assembly"] for i in curated_intervals}
    df_baselines = bm.compute_all_baselines(sp_asm_map)
    df_baselines.to_csv(os.path.join(output_dir, "genome_baseline_stats.csv"), index=False)

    # 4. Standardized Z-Scores
    click.echo("\n[4/6] Computing standardized Z-scores and hypothesis tests...")
    za = ZScoreAnalyzer(bm)
    df_zscores, df_comp_stats = za.analyze_genes(all_annotated_genes, all_gene_gc_records)
    df_zscores.to_csv(os.path.join(output_dir, "gene_zscores.csv"), index=False)
    df_comp_stats.to_csv(os.path.join(output_dir, "statistical_comparisons.csv"), index=False)

    # 5. Diagnostic Plots
    diag_dir = os.path.join(output_dir, "diagnostic_plots")
    click.echo(f"\n[5/6] Generating diagnostic figures in: {diag_dir}...")
    viz = Visualizer(bm)
    for sp in unique_species:
        click.echo(f"  Rendering diagnostic plot for {sp}...")
        try:
            viz.plot_diagnostic(sp, df_zscores, output_dir=diag_dir)
            viz.plot_diagnostic_interactive(sp, df_zscores, output_dir=diag_dir)
        except Exception as e:
            click.echo(f"    Warning: plot error for {sp}: {e}")

    # 6. Cross-Species Meta-Graph
    meta_png = os.path.join(output_dir, "meta_comparative_graph.png")
    meta_html = os.path.join(output_dir, "meta_comparative_graph.html")
    meta_clean_png = os.path.join(output_dir, "meta_comparative_graph_clean.png")
    meta_clean_html = os.path.join(output_dir, "meta_comparative_graph_clean.html")
    click.echo(f"\n[6/6] Generating cross-species meta-graphs (comprehensive & clean): {meta_png}...")
    viz.plot_meta_comparative(df_zscores, output_path=meta_png, clean_mode=False)
    viz.plot_meta_interactive(df_zscores, output_path=meta_html, clean_mode=False)
    viz.plot_meta_comparative(df_zscores, output_path=meta_clean_png, clean_mode=True)
    viz.plot_meta_interactive(df_zscores, output_path=meta_clean_html, clean_mode=True)

    click.echo("\n=== Manifest Pipeline Completed Successfully! ===")


@main.command(name="add-genome")
@click.argument("accession")
@click.option("--common-name", "-n", required=True, help="Common name (e.g. 'Burmese Python').")
@click.option("--clade", "-c", required=True, help="Taxonomic clade (e.g. 'amniote_reptile', 'amniote_bird').")
@click.option("--chrom", default=None, help="Chromosome accession if manually specifying coordinates.")
@click.option("--start", default=None, type=int, help="Start coordinate in bp.")
@click.option("--end", default=None, type=int, help="End coordinate in bp.")
@click.option("--locus-id", default=None, help="Locus ID ('alpha_globin', 'beta_globin', 'alpha_beta_linked').")
@click.option("--auto-detect/--no-auto-detect", default=True, help="Auto-detect globin synteny clusters.")
def add_genome_cmd(accession, common_name, clade, chrom, start, end, locus_id, auto_detect):
    """Imports an NCBI RefSeq genome by accession, streams annotations, detects/accepts globin loci,
    and computes clean empirical protein-coding baselines.
    """
    from globin_pipeline.api import GlobinPipeline
    from globin_pipeline.exceptions import GlobinPipelineError

    click.echo(f"=== Importing NCBI Genome: {accession} ({common_name}) ===")
    pipeline = GlobinPipeline()

    try:
        res = pipeline.add_genome(
            accession=accession,
            common_name=common_name,
            clade=clade,
            chrom=chrom,
            start=start,
            end=end,
            locus_id=locus_id,
            auto_detect=auto_detect,
        )
        click.secho(f"\n[SUCCESS] Genome resolved: {res['assembly_name']}", fg="green", bold=True)
        click.echo(f"  Clean protein-coding genes: {res['total_protein_coding_genes']:,}")
        click.echo(f"  Chromosomes/Scaffolds:      {res['chromosomes_available']:,}")
        
        click.echo(f"\n  Detected/Curated Globin Loci ({len(res['candidate_loci'])}):")
        for idx, loc in enumerate(res['candidate_loci'], 1):
            click.secho(f"    [{idx}] {loc.locus_id} ({loc.suggested_type}): {loc.chromosome}:{loc.start:,}-{loc.end:,} ({loc.span_bp:,} bp)", fg="cyan")
            click.echo(f"        Confidence: {loc.confidence*100:.0f}% | {loc.notes}")
            for g in loc.globin_genes[:6]:
                click.echo(f"          - Globin: {g.get('symbol')} ({g.get('start'):,}-{g.get('end'):,})")
            for a in loc.anchor_genes[:4]:
                click.echo(f"          - Anchor: {a.get('symbol')} ({a.get('start'):,}-{a.get('end'):,})")

    except GlobinPipelineError as e:
        click.secho(f"\n[ERROR] Pipeline Error: {e}", fg="red", err=True)
        raise SystemExit(1)
    except Exception as e:
        click.secho(f"\n[FATAL] Unexpected error importing {accession}: {e}", fg="red", err=True)
        raise SystemExit(1)


@main.command(name="export-tables")
@click.option("--output-dir", "-o", default="results", help="Directory to save supplementary tables.")
def export_tables_cmd(output_dir):
    """Exports Supplementary Tables S1-S4 and cross-species summary tables."""
    import subprocess
    click.echo("=== Building and Exporting Data Release Tables ===")
    script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "scripts", "build_supplementary_data_package.py")
    if os.path.exists(script):
        subprocess.run(["python3", script], check=True)
    else:
        click.secho(f"Export script not found: {script}", fg="red")


if __name__ == "__main__":
    main()

