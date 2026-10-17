"""Interactive Streamlit Dashboard for Vertebrate Globin Comparative Genomics."""

import os
import sys
import importlib.resources as pkg_resources
from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from globin_pipeline.visualizer import Visualizer


st.set_page_config(
    page_title="Vertebrate Globin Genomics Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for scientific typography and clean card styling
st.markdown("""
<style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1e293b;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.05rem;
        color: #64748b;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 1rem;
        margin-bottom: 0.8rem;
    }
    .metric-val {
        font-size: 1.6rem;
        font-weight: 700;
        color: #0f172a;
    }
    .metric-lbl {
        font-size: 0.85rem;
        color: #64748b;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .species-tag {
        display: inline-block;
        padding: 0.2rem 0.6rem;
        border-radius: 12px;
        font-size: 0.8rem;
        font-weight: 600;
        background-color: #e0e7ff;
        color: #3730a3;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_data
def load_datasets():
    """Load reference datasets with fallback to bundled package data."""
    def resolve_csv(filename):
        candidates = [
            PROJECT_ROOT / "results" / filename,
            PROJECT_ROOT / "1_data_from_mike" / "results" / filename,
        ]
        for c in candidates:
            if c.exists():
                return pd.read_csv(c)
        try:
            ref = pkg_resources.files("globin_pipeline.data") / filename
            return pd.read_csv(ref)
        except Exception:
            return pd.DataFrame()

    df_genes = resolve_csv("gene_zscores.csv")
    if df_genes.empty:
        df_genes = resolve_csv("Supplementary_Table_S3_All_1017_Genes_Annotations_ZScores.csv")

    df_intervals = resolve_csv("curated_intervals.csv")
    if df_intervals.empty:
        df_intervals = resolve_csv("Supplementary_Table_S2_Globin_Loci_Summary.csv")

    df_baselines = resolve_csv("genome_baseline_stats.csv")
    if df_baselines.empty:
        df_baselines = resolve_csv("Supplementary_Table_S1_Species_Baselines.csv")

    df_stats = resolve_csv("Supplementary_Table_S4_Statistical_Contrasts.csv")
    if df_stats.empty:
        df_stats = resolve_csv("statistical_comparisons.csv")

    manifest_species = {}
    manifest_candidates = [
        PROJECT_ROOT / "config" / "species_manifest.yaml",
        PROJECT_ROOT / "globin_pipeline" / "data" / "species_manifest.yaml",
    ]
    for mc in manifest_candidates:
        if mc.exists():
            try:
                with open(mc) as f:
                    manifest_species = yaml.safe_load(f).get("species", {})
                    break
            except Exception:
                pass

    return df_genes, df_intervals, df_baselines, df_stats, manifest_species


df_genes, df_intervals, df_baselines, df_stats, manifest_species = load_datasets()

# Initialize Visualizer
@st.cache_resource
def get_visualizer():
    return Visualizer()

viz = get_visualizer()


def get_or_render_plot(species: str, plot_type: str) -> str:
    """Returns the path to a requested plot, generating it on-the-fly if missing."""
    results_dir = PROJECT_ROOT / "results"
    
    if plot_type == "chromosome":
        target = results_dir / "chromosome_plots" / f"chromosome_linear_{species}.png"
        if not target.exists():
            viz.plot_chromosome_linear(species, df_genes, df_intervals=df_intervals, output_dir=str(results_dir / "chromosome_plots"))
        return str(target) if target.exists() else None

    elif plot_type == "inset":
        target = results_dir / "bivariate_with_inset_plots" / f"bivariate_with_inset_{species}.png"
        if not target.exists():
            viz.plot_bivariate_with_chromosome_inset(species, df_genes, df_intervals=df_intervals, output_dir=str(results_dir / "bivariate_with_inset_plots"))
        return str(target) if target.exists() else None

    elif plot_type == "diagnostic":
        target = results_dir / "diagnostic_plots" / f"diagnostic_{species}.png"
        if not target.exists():
            viz.plot_diagnostic(species, df_genes, output_dir=str(results_dir / "diagnostic_plots"))
        return str(target) if target.exists() else None

    elif plot_type == "bivariate":
        target = results_dir / "bivariate_length_gc_plots" / f"bivariate_length_gc_{species}.png"
        if not target.exists():
            viz.plot_length_vs_gc(species, df_genes, output_dir=str(results_dir / "bivariate_length_gc_plots"))
        return str(target) if target.exists() else None

    return None


# Sidebar Navigation
st.sidebar.markdown("## 🧬 Globin Explorer")
mode = st.sidebar.radio(
    "Navigation Mode",
    [
        "⚖️ 1-on-1 Head-to-Head Comparison",
        "🔬 Single Species Deep Dive",
        "🌐 Cross-Species Meta Explorer",
        "📋 Locus & Gene Data Table",
        "📊 Statistical Contrasts & Hypotheses",
    ],
    index=0,
)

species_list = sorted(df_genes["species"].unique()) if not df_genes.empty else []

# Helper: format species display name
def format_sp_name(sp):
    if not df_intervals.empty and sp in df_intervals["species"].values:
        common = df_intervals[df_intervals["species"] == sp]["common_name"].iloc[0]
        clade = df_intervals[df_intervals["species"] == sp]["clade"].iloc[0]
        return f"{sp} ({common} • {clade})"
    return sp


# ==============================================================================
# VIEW 1: 1-on-1 HEAD-TO-HEAD COMPARISON
# ==============================================================================
if mode == "⚖️ 1-on-1 Head-to-Head Comparison":
    st.markdown('<div class="main-title">⚖️ 1-on-1 Head-to-Head Species Comparison</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Directly contrast globin genomic architectures, isochore coupling, and base composition between any two vertebrate genomes.</div>', unsafe_allow_html=True)

    col_sel_a, col_sel_b = st.columns(2)
    with col_sel_a:
        default_a = species_list.index("Homo_sapiens") if "Homo_sapiens" in species_list else 0
        sp_a = st.selectbox("Select Species A", species_list, index=default_a, format_func=format_sp_name, key="sp_a")
    with col_sel_b:
        default_b = species_list.index("Gallus_gallus") if "Gallus_gallus" in species_list else (1 if len(species_list) > 1 else 0)
        sp_b = st.selectbox("Select Species B", species_list, index=default_b, format_func=format_sp_name, key="sp_b")

    st.markdown("---")

    # Comparative Metrics Summary Scorecard
    st.markdown("### 📊 Architectural & Compositional Contrast")
    
    def get_species_metrics(sp):
        m = {}
        # Baselines
        base_row = df_baselines[df_baselines["species"] == sp] if not df_baselines.empty else pd.DataFrame()
        b = base_row.iloc[0] if not base_row.empty else {}

        # Clade
        b_clade = str(b.get("clade", "N/A"))
        m["clade"] = b_clade if b_clade not in ["nan", "None", "N/A", ""] else "N/A"

        # Assembly resolution with manifest fallback
        asm = str(b.get("assembly", "")).strip()
        if not asm or asm in ["nan", "None", "N/A"] or asm.startswith("NC_"):
            m_asm = manifest_species.get(sp, {}).get("assembly")
            asm = str(m_asm).strip() if m_asm else "N/A"
        if not asm or asm in ["nan", "None"]:
            asm = "N/A"
        m["assembly"] = asm

        m["mean_cds_gc"] = float(b.get("mean_cds_gc", 0.0)) if pd.notna(b.get("mean_cds_gc")) else 0.0
        m["std_cds_gc"] = float(b.get("std_cds_gc", 0.0)) if pd.notna(b.get("std_cds_gc")) else 0.0

        # Gene count resolution across all possible columns
        gene_count = "N/A"
        for col in ["num_genes", "total_genes", "num_cds_genes", "num_whole_genes"]:
            val = b.get(col)
            if pd.notna(val) and float(val) > 0:
                gene_count = int(float(val))
                break
        m["gene_count"] = gene_count

        # Loci & Intervals
        sp_genes = df_genes[df_genes["species"] == sp]
        sp_intervals = df_intervals[df_intervals["species"] == sp]
        m["intervals_count"] = len(sp_intervals)
        if m["clade"] == "N/A" and not sp_intervals.empty:
            m["clade"] = str(sp_intervals["clade"].iloc[0])
        
        globins = sp_genes[sp_genes["is_globin"] == True]

        def _is_alpha(r):
            grp = str(r.get("comparative_group", "")).lower()
            sym = str(r.get("display_symbol", "")).lower()
            name = str(r.get("gene_name", "")).lower()
            if "alpha" in grp: return True
            if "beta" in grp: return False
            return any(sym.startswith(p) for p in ["hba", "hbaa", "hbae", "hbz", "hbm", "hbq"]) or "hba" in sym or "hba" in name

        def _is_beta(r):
            grp = str(r.get("comparative_group", "")).lower()
            sym = str(r.get("display_symbol", "")).lower()
            name = str(r.get("gene_name", "")).lower()
            if "beta" in grp: return True
            if "alpha" in grp: return False
            return any(sym.startswith(p) for p in ["hbb", "hbba", "hbbe", "hbd", "hbg", "ba1l", "gb-beta", "rho"]) or "hbb" in sym or "hbb" in name

        alpha_genes = globins[globins.apply(_is_alpha, axis=1)]
        beta_genes = globins[globins.apply(_is_beta, axis=1)]
        z_col = "z_score_full_gene_gc" if "z_score_full_gene_gc" in sp_genes.columns else ("zscore_full_gene" if "zscore_full_gene" in sp_genes.columns else None)

        m["alpha_gc"] = float(alpha_genes["cds_gc"].mean()) if not alpha_genes.empty else np.nan
        m["alpha_z"] = float(alpha_genes[z_col].mean()) if (z_col and not alpha_genes.empty) else np.nan
        m["beta_gc"] = float(beta_genes["cds_gc"].mean()) if not beta_genes.empty else np.nan
        m["beta_z"] = float(beta_genes[z_col].mean()) if (z_col and not beta_genes.empty) else np.nan

        # Proto-globin for cyclostomes
        m["proto_gc"] = float(globins["cds_gc"].mean()) if (np.isnan(m["alpha_gc"]) and np.isnan(m["beta_gc"]) and not globins.empty) else np.nan
        m["proto_z"] = float(globins[z_col].mean()) if (np.isnan(m["alpha_z"]) and np.isnan(m["beta_z"]) and z_col and not globins.empty) else np.nan

        # Polarity
        if not np.isnan(m["alpha_z"]) and not np.isnan(m["beta_z"]):
            m["polarity"] = m["alpha_z"] - m["beta_z"]
        else:
            m["polarity"] = np.nan

        return m

    m_a = get_species_metrics(sp_a)
    m_b = get_species_metrics(sp_b)

    def _format_alpha(m_dict):
        if not np.isnan(m_dict["alpha_gc"]):
            return f"{m_dict['alpha_gc']:.1f}% (Z = {m_dict['alpha_z']:+.2f}σ)"
        if not np.isnan(m_dict["proto_gc"]):
            return f"Proto: {m_dict['proto_gc']:.1f}% (Z = {m_dict['proto_z']:+.2f}σ)"
        return "N/A"

    def _format_beta(m_dict):
        if not np.isnan(m_dict["beta_gc"]):
            return f"{m_dict['beta_gc']:.1f}% (Z = {m_dict['beta_z']:+.2f}σ)"
        if not np.isnan(m_dict["proto_gc"]):
            return "N/A (Ancestral Proto-Globin)"
        return "N/A"

    def _format_polarity(m_dict):
        if not np.isnan(m_dict["polarity"]):
            return f"{m_dict['polarity']:+.2f}σ"
        if not np.isnan(m_dict["proto_gc"]):
            return "N/A (Monomeric Proto-Globin)"
        return "N/A"

    def _format_polarity_diff(m1, m2):
        if not np.isnan(m1["polarity"]) and not np.isnan(m2["polarity"]):
            if (m1["polarity"] > 0 and m2["polarity"] < 0) or (m1["polarity"] < 0 and m2["polarity"] > 0):
                return "⚡ Polarity Inversion"
            return "Conserved Polarity"
        return "N/A"

    scorecard_data = [
        {"Metric": "Phylogenetic Clade", sp_a: m_a["clade"], sp_b: m_b["clade"], "Contrast / Difference": "Same Clade" if m_a["clade"] == m_b["clade"] else "Cross-Clade"},
        {"Metric": "Genome Assembly Accession", sp_a: m_a["assembly"], sp_b: m_b["assembly"], "Contrast / Difference": "—"},
        {"Metric": "Genome Baseline CDS GC (%)", sp_a: f"{m_a['mean_cds_gc']:.2f}% (±{m_a['std_cds_gc']:.2f}%)", sp_b: f"{m_b['mean_cds_gc']:.2f}% (±{m_b['std_cds_gc']:.2f}%)", "Contrast / Difference": f"{m_a['mean_cds_gc'] - m_b['mean_cds_gc']:+.2f}%"},
        {"Metric": "Curated Globin Intervals", sp_a: f"{m_a['intervals_count']} loci", sp_b: f"{m_b['intervals_count']} loci", "Contrast / Difference": "—"},
        {"Metric": "α-Globin Coding GC / Whole-Gene Z", sp_a: _format_alpha(m_a), sp_b: _format_alpha(m_b), "Contrast / Difference": f"ΔZ = {m_a['alpha_z'] - m_b['alpha_z']:+.2f}σ" if not np.isnan(m_a['alpha_z']) and not np.isnan(m_b['alpha_z']) else "N/A"},
        {"Metric": "β-Globin Coding GC / Whole-Gene Z", sp_a: _format_beta(m_a), sp_b: _format_beta(m_b), "Contrast / Difference": f"ΔZ = {m_a['beta_z'] - m_b['beta_z']:+.2f}σ" if not np.isnan(m_a['beta_z']) and not np.isnan(m_b['beta_z']) else "N/A"},
        {"Metric": "Subunit Polarity (ΔZ = α - β)", sp_a: _format_polarity(m_a), sp_b: _format_polarity(m_b), "Contrast / Difference": _format_polarity_diff(m_a, m_b)},
    ]
    st.table(pd.DataFrame(scorecard_data).set_index("Metric"))

    st.markdown("---")

    # Synchronized Visualizer Comparison Tabs
    st.markdown("### 🖼️ Synchronized Visual Comparison")
    comp_plot_type = st.radio(
        "Select Comparison Plot View",
        [
            "📍 Chromosome Linear Tracks (Centromeres & Telomeres)",
            "🎯 Bivariate with Chromosome Mini-Inset",
            "📈 Whole-Genome Diagnostic GC Distributions",
            "📐 Standard Bivariate (Length vs. GC)",
        ],
        horizontal=True,
    )

    plot_code_map = {
        "📍 Chromosome Linear Tracks (Centromeres & Telomeres)": "chromosome",
        "🎯 Bivariate with Chromosome Mini-Inset": "inset",
        "📈 Whole-Genome Diagnostic GC Distributions": "diagnostic",
        "📐 Standard Bivariate (Length vs. GC)": "bivariate",
    }
    p_code = plot_code_map[comp_plot_type]

    col_img_a, col_img_b = st.columns(2)
    with col_img_a:
        st.markdown(f"#### {format_sp_name(sp_a)}")
        p_path_a = get_or_render_plot(sp_a, p_code)
        if p_path_a and os.path.exists(p_path_a):
            st.image(p_path_a, use_container_width=True)
        else:
            st.info(f"Plot not available for {sp_a}.")

    with col_img_b:
        st.markdown(f"#### {format_sp_name(sp_b)}")
        p_path_b = get_or_render_plot(sp_b, p_code)
        if p_path_b and os.path.exists(p_path_b):
            st.image(p_path_b, use_container_width=True)
        else:
            st.info(f"Plot not available for {sp_b}.")


# ==============================================================================
# VIEW 2: SINGLE SPECIES DEEP DIVE
# ==============================================================================
elif mode == "🔬 Single Species Deep Dive":
    st.markdown('<div class="main-title">🔬 Single Species Deep Dive</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Inspect full chromosome architectures, telomeric margins, GC distributions, and individual gene models.</div>', unsafe_allow_html=True)

    col1, col2 = st.columns([1, 3])
    with col1:
        clade_filter = st.selectbox("Filter by Clade (optional)", ["All"] + sorted(df_genes["clade"].dropna().unique().tolist()))
        if clade_filter != "All":
            filtered_species = sorted(df_genes[df_genes["clade"] == clade_filter]["species"].unique())
        else:
            filtered_species = species_list

        selected_species = st.selectbox("Select Target Species", filtered_species, format_func=format_sp_name)

        # Quick stats cards
        sp_int = df_intervals[df_intervals["species"] == selected_species]
        sp_base = df_baselines[df_baselines["species"] == selected_species] if not df_baselines.empty else pd.DataFrame()
        
        st.markdown("#### ℹ️ Species Profile")
        sm = get_species_metrics(selected_species)
        common_name = df_intervals[df_intervals["species"] == selected_species]["common_name"].iloc[0] if not sp_int.empty else selected_species
        st.markdown(f"**Common Name**: {common_name}")
        st.markdown(f"**Clade**: `{sm['clade']}`")
        st.markdown(f"**Assembly**: `{sm['assembly']}`")
        st.markdown(f"**Baseline CDS GC**: `{sm['mean_cds_gc']:.2f}%` (± `{sm['std_cds_gc']:.2f}%`)")
        if sm["gene_count"] != "N/A":
            st.markdown(f"**Total Coding Genes**: `{sm['gene_count']:,}`")
        else:
            st.markdown(f"**Total Coding Genes**: `N/A`")
        st.markdown(f"**Curated Loci**: `{len(sp_int)} intervals`")
        
        # Globin Subunit Metrics
        if not np.isnan(sm["alpha_gc"]):
            st.markdown(f"**α-Globin GC**: `{sm['alpha_gc']:.1f}%` (Z = `{sm['alpha_z']:+.2f}σ`)")
            st.markdown(f"**β-Globin GC**: `{sm['beta_gc']:.1f}%` (Z = `{sm['beta_z']:+.2f}σ`)")
            st.markdown(f"**Subunit Polarity**: `{sm['polarity']:+.2f}σ`")
        elif not np.isnan(sm["proto_gc"]):
            st.markdown(f"**Proto-Globin GC**: `{sm['proto_gc']:.1f}%` (Z = `{sm['proto_z']:+.2f}σ`)")

    with col2:
        tab_chr, tab_inset, tab_diag, tab_biv, tab_table = st.tabs([
            "📍 Chromosome Linear Track",
            "🎯 Bivariate with Mini-Inset",
            "📈 Diagnostic GC Distribution",
            "📐 Bivariate Length vs. GC",
            "📋 Interval Genes Table",
        ])

        with tab_chr:
            st.markdown("##### Linear Chromosome Arm Track (Centromeres & Telomeres)")
            p = get_or_render_plot(selected_species, "chromosome")
            if p and os.path.exists(p):
                st.image(p, use_container_width=True)
            else:
                st.info("Chromosome plot rendering...")

        with tab_inset:
            st.markdown("##### Bivariate log10(Coding Length) vs. GC with Mini Chromosome Inset")
            p = get_or_render_plot(selected_species, "inset")
            if p and os.path.exists(p):
                st.image(p, use_container_width=True)
            else:
                st.info("Bivariate with inset rendering...")

        with tab_diag:
            st.markdown("##### Whole-Genome GC Distribution & Individual Globin Percentiles")
            p = get_or_render_plot(selected_species, "diagnostic")
            if p and os.path.exists(p):
                st.image(p, use_container_width=True)
            else:
                st.info("Diagnostic plot rendering...")

        with tab_biv:
            st.markdown("##### Standard Bivariate log10(Length) vs. GC Scatter")
            p = get_or_render_plot(selected_species, "bivariate")
            if p and os.path.exists(p):
                st.image(p, use_container_width=True)
            else:
                st.info("Bivariate plot rendering...")

        with tab_table:
            st.markdown(f"##### Annotated Genes in {selected_species} Globin Intervals")
            sp_genes = df_genes[df_genes["species"] == selected_species].copy()
            show_cols = [c for c in ["display_symbol", "comparative_group", "chromosome", "start", "end", "strand", "exon_count", "cds_gc", "full_gene_gc", "z_score_full_gene_gc", "z_score_gc", "zscore_full_gene", "zscore_cds"] if c in sp_genes.columns]
            st.dataframe(sp_genes[show_cols], use_container_width=True)


# ==============================================================================
# VIEW 3: CROSS-SPECIES META EXPLORER
# ==============================================================================
elif mode == "🌐 Cross-Species Meta Explorer":
    st.markdown('<div class="main-title">🌐 Cross-Species Meta-Comparative Explorer</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Explore global macroevolutionary patterns across all 44 vertebrate species.</div>', unsafe_allow_html=True)

    tab_meta_scatter, tab_divergence, tab_drivers, tab_heatmap = st.tabs([
        "🌌 Meta-Comparative Standardized Graph",
        "📊 Summary 1: Locus Divergence",
        "🔬 Summary 2: Mechanistic Drivers (gBGC & Flankers)",
        "🗺️ Summary 3: Clade Composition Heatmap",
    ])

    results_dir = PROJECT_ROOT / "results"

    with tab_meta_scatter:
        clean_img = results_dir / "meta_comparative_graph_clean.png"
        full_img = results_dir / "meta_comparative_graph.png"
        show_clean = st.checkbox("Show Clean Filtered Version (Globins & Syntenic Anchors Only)", value=True)
        img_to_show = clean_img if show_clean and clean_img.exists() else full_img
        if img_to_show.exists():
            st.image(str(img_to_show), use_container_width=True)
        else:
            st.info("Meta graph not found. Run `globin-pipeline generate-plots` to render.")

    with tab_divergence:
        f1 = results_dir / "summary_figure1_locus_divergence.png"
        if f1.exists():
            st.image(str(f1), use_container_width=True)
        else:
            st.info("Summary figure 1 not found.")

    with tab_drivers:
        f2 = results_dir / "summary_figure2_mechanistic_drivers.png"
        if f2.exists():
            st.image(str(f2), use_container_width=True)
        else:
            st.info("Summary figure 2 not found.")

    with tab_heatmap:
        f3 = results_dir / "summary_figure3_clade_heatmap.png"
        if f3.exists():
            st.image(str(f3), use_container_width=True)
        else:
            st.info("Summary figure 3 not found.")


# ==============================================================================
# VIEW 4: LOCUS & GENE DATA TABLE
# ==============================================================================
elif mode == "📋 Locus & Gene Data Table":
    st.markdown('<div class="main-title">📋 Locus & Gene Repertoire Table</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Interactive query engine for all 1,017 genes across 80 intervals in 44 species.</div>', unsafe_allow_html=True)

    col_f1, col_f2, col_f3 = st.columns(3)
    with col_f1:
        c_sel = st.multiselect("Filter by Clade", sorted(df_genes["clade"].dropna().unique().tolist()))
    with col_f2:
        grp_sel = st.multiselect("Filter by Functional Group", sorted(df_genes["comparative_group"].dropna().unique().tolist()))
    with col_f3:
        search_kw = st.text_input("Search Gene Symbol or Species", "")

    filtered_df = df_genes.copy()
    if c_sel:
        filtered_df = filtered_df[filtered_df["clade"].isin(c_sel)]
    if grp_sel:
        filtered_df = filtered_df[filtered_df["comparative_group"].isin(grp_sel)]
    if search_kw:
        mask = filtered_df["display_symbol"].str.contains(search_kw, case=False, na=False) | filtered_df["species"].str.contains(search_kw, case=False, na=False)
        filtered_df = filtered_df[mask]

    st.markdown(f"**Displaying {len(filtered_df):,} of {len(df_genes):,} genes**")
    st.dataframe(filtered_df, use_container_width=True)

    csv_data = filtered_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Filtered Data as CSV",
        data=csv_data,
        file_name="globin_genes_filtered.csv",
        mime="text/csv",
    )


# ==============================================================================
# VIEW 5: STATISTICAL CONTRASTS & HYPOTHESES
# ==============================================================================
elif mode == "📊 Statistical Contrasts & Hypotheses":
    st.markdown('<div class="main-title">📊 Statistical Contrasts & Hypothesis Testing</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Formal hypothesis testing results evaluating Natural Selection vs. gBGC vs. Isochore Inertia.</div>', unsafe_allow_html=True)

    if not df_stats.empty:
        st.dataframe(df_stats, use_container_width=True)
    else:
        st.info("Statistical contrasts table not found.")
