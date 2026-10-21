"""Generates a standalone, interactive HTML dashboard (results/dashboard.html)
allowing instant exploration and 1-on-1 head-to-head comparison across all 44 vertebrate species.
Requires zero server or Python runtime to view—can be double-clicked in any browser.
"""

import json
import os
from pathlib import Path
import yaml
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent

def main():
    print("Loading datasets for static dashboard...")
    df_genes = pd.read_csv(PROJECT_ROOT / "results" / "gene_zscores.csv")
    df_intervals = pd.read_csv(PROJECT_ROOT / "results" / "curated_intervals.csv")
    df_baselines = pd.read_csv(PROJECT_ROOT / "results" / "genome_baseline_stats.csv")
    df_stats = pd.read_csv(PROJECT_ROOT / "results" / "Supplementary_Table_S4_Statistical_Contrasts.csv")

    manifest_path = PROJECT_ROOT / "config" / "species_manifest.yaml"
    manifest_species = {}
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest_species = yaml.safe_load(f).get("species", {})

    species_list = sorted(df_genes["species"].unique())
    print(f"Loaded {len(species_list)} species.")

    # Build species metadata dictionary
    species_data = {}
    for sp in species_list:
        sp_int = df_intervals[df_intervals["species"] == sp]
        common = sp_int["common_name"].iloc[0] if not sp_int.empty else sp
        clade = sp_int["clade"].iloc[0] if not sp_int.empty else "unknown"

        # Baseline stats
        base_row = df_baselines[df_baselines["species"] == sp]
        b = base_row.iloc[0] if not base_row.empty else {}

        # Assembly resolution
        assembly = str(b.get("assembly", "")).strip()
        if not assembly or assembly in ["nan", "None", "N/A"] or assembly.startswith("NC_"):
            m_asm = manifest_species.get(sp, {}).get("assembly")
            if m_asm:
                assembly = str(m_asm).strip()
        if not assembly or assembly in ["nan", "None"]:
            assembly = "N/A"

        mean_cds = float(b.get("mean_cds_gc", 0.0)) if pd.notna(b.get("mean_cds_gc")) else 0.0
        std_cds = float(b.get("std_cds_gc", 0.0)) if pd.notna(b.get("std_cds_gc")) else 0.0

        # Gene count resolution across all available columns
        num_genes = "N/A"
        for col in ["num_genes", "total_genes", "num_cds_genes", "num_whole_genes"]:
            val = b.get(col)
            if pd.notna(val) and float(val) > 0:
                num_genes = int(float(val))
                break

        # Globin locus metrics
        sp_genes = df_genes[df_genes["species"] == sp]
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
        alpha_gc = float(alpha_genes["cds_gc"].mean()) if not alpha_genes.empty else None
        alpha_z = float(alpha_genes[z_col].mean()) if (z_col and not alpha_genes.empty) else None
        beta_gc = float(beta_genes["cds_gc"].mean()) if not beta_genes.empty else None
        beta_z = float(beta_genes[z_col].mean()) if (z_col and not beta_genes.empty) else None

        # Handle cyclostomes (ancestral proto-globins)
        proto_gc = float(globins["cds_gc"].mean()) if (alpha_gc is None and beta_gc is None and not globins.empty) else None
        proto_z = float(globins[z_col].mean()) if (alpha_z is None and beta_z is None and z_col and not globins.empty) else None

        polarity = (alpha_z - beta_z) if (alpha_z is not None and beta_z is not None) else None

        species_data[sp] = {
            "species": sp,
            "common_name": common,
            "clade": clade,
            "assembly": assembly,
            "mean_cds_gc": mean_cds,
            "std_cds_gc": std_cds,
            "num_genes": num_genes,
            "loci_count": len(sp_int),
            "alpha_gc": alpha_gc,
            "alpha_z": alpha_z,
            "beta_gc": beta_gc,
            "beta_z": beta_z,
            "proto_gc": proto_gc,
            "proto_z": proto_z,
            "polarity": polarity,
        }

    # Prepare gene records for search table (compact JSON)
    genes_data = []
    for _, g in df_genes.iterrows():
        z_g = g.get("z_score_full_gene_gc") if pd.notna(g.get("z_score_full_gene_gc")) else g.get("zscore_full_gene", 0.0)
        z_c = g.get("z_score_gc") if pd.notna(g.get("z_score_gc")) else g.get("zscore_cds", 0.0)
        genes_data.append({
            "species": str(g.get("species", "")),
            "symbol": str(g.get("display_symbol", "")),
            "group": str(g.get("comparative_group", "")),
            "clade": str(g.get("clade", "")),
            "chrom": str(g.get("chromosome", "")),
            "start": int(g.get("start", 0)),
            "end": int(g.get("end", 0)),
            "cds_gc": round(float(g.get("cds_gc", 0.0)), 1) if pd.notna(g.get("cds_gc")) else 0.0,
            "gene_gc": round(float(g.get("full_gene_gc", 0.0)), 1) if pd.notna(g.get("full_gene_gc")) else 0.0,
            "z_gene": round(float(z_g), 2) if pd.notna(z_g) else 0.0,
            "z_cds": round(float(z_c), 2) if pd.notna(z_c) else 0.0,
        })

    # Convert stats table to list of dicts
    stats_data = df_stats.to_dict(orient="records")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Vertebrate Globin Genomics Explorer</title>
    <style>
        :root {{
            --primary: #2563eb;
            --primary-hover: #1d4ed8;
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --border: #e2e8f0;
            --text-main: #0f172a;
            --text-muted: #64748b;
            --accent: #4f46e5;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg);
            color: var(--text-main);
            line-height: 1.45;
            padding: 0.75rem 1.25rem;
        }}
        header {{
            background: linear-gradient(135deg, #1e293b, #0f172a);
            color: white;
            padding: 0.6rem 1.25rem;
            border-radius: 8px;
            margin-bottom: 0.65rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.08);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 0.5rem;
        }}
        .header-left h1 {{
            font-size: 1.25rem;
            font-weight: 700;
            margin: 0;
            line-height: 1.2;
            color: #f8fafc;
        }}
        .header-left p {{
            color: #94a3b8;
            font-size: 0.8rem;
            margin: 0.15rem 0 0 0;
        }}
        .badge {{
            display: inline-block;
            background: #334155;
            color: #38bdf8;
            padding: 0.2rem 0.55rem;
            border-radius: 9999px;
            font-size: 0.74rem;
            font-weight: 600;
            white-space: nowrap;
        }}

        nav.tabs {{
            display: flex;
            gap: 0.25rem;
            border-bottom: 2px solid var(--border);
            margin-bottom: 0.65rem;
            overflow-x: auto;
        }}
        nav.tabs button {{
            background: none;
            border: none;
            padding: 0.45rem 0.85rem;
            font-size: 0.88rem;
            font-weight: 600;
            color: var(--text-muted);
            cursor: pointer;
            border-bottom: 3px solid transparent;
            transition: all 0.2s;
            white-space: nowrap;
        }}
        nav.tabs button:hover {{ color: var(--primary); }}
        nav.tabs button.active {{
            color: var(--primary);
            border-bottom-color: var(--primary);
        }}

        .tab-pane {{ display: none; }}
        .tab-pane.active {{ display: block; }}

        .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1.25rem; }}
        @media (max-width: 900px) {{ .grid-2 {{ grid-template-columns: 1fr; }} }}

        .card {{
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 0.85rem 1.1rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
            margin-bottom: 0.85rem;
        }}
        .card h3 {{ font-size: 1.05rem; font-weight: 600; margin-bottom: 0.6rem; color: #1e293b; }}

        select, input {{
            width: 100%;
            padding: 0.6rem 0.8rem;
            border: 1px solid var(--border);
            border-radius: 6px;
            font-size: 0.95rem;
            background: white;
            color: var(--text-main);
            margin-bottom: 0.8rem;
        }}
        select:focus, input:focus {{ outline: none; border-color: var(--primary); }}

        /* COMBOBOX DYNAMIC SEARCH */
        .combobox-wrapper {{
            position: relative;
            width: 100%;
            margin-bottom: 0.8rem;
        }}
        .combobox-trigger {{
            width: 100%;
            padding: 0.55rem 0.8rem;
            background: white;
            border: 1px solid var(--border);
            border-radius: 8px;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: space-between;
            transition: all 0.2s;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04);
            user-select: none;
        }}
        .combobox-trigger:hover, .combobox-trigger:focus {{
            border-color: var(--primary);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.12);
            outline: none;
        }}
        .combobox-selection {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
            overflow: hidden;
            white-space: nowrap;
            text-overflow: ellipsis;
            font-size: 0.92rem;
        }}
        .combobox-name {{
            font-weight: 700;
            color: #0f172a;
        }}
        .combobox-common {{
            color: #475569;
            font-size: 0.86rem;
        }}
        .combobox-clade {{
            background: #e2e8f0;
            color: #334155;
            padding: 0.15rem 0.45rem;
            border-radius: 4px;
            font-size: 0.72rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.03em;
        }}
        .combobox-arrow {{
            font-size: 0.8rem;
            color: #94a3b8;
            transition: transform 0.2s;
            margin-left: 0.5rem;
        }}
        .combobox-wrapper.open .combobox-arrow {{
            transform: rotate(180deg);
            color: var(--primary);
        }}
        .combobox-dropdown {{
            display: none;
            position: absolute;
            top: calc(100% + 4px);
            left: 0;
            width: 100%;
            background: white;
            border: 1px solid var(--border);
            border-radius: 8px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.15), 0 8px 10px -6px rgba(0, 0, 0, 0.1);
            z-index: 1000;
            overflow: hidden;
            animation: fadeInDown 0.15s ease-out;
        }}
        @keyframes fadeInDown {{
            from {{ opacity: 0; transform: translateY(-6px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}
        .combobox-wrapper.open .combobox-dropdown {{
            display: block;
        }}
        .combobox-search-box {{
            padding: 0.6rem;
            border-bottom: 1px solid var(--border);
            background: #f8fafc;
        }}
        .combobox-search-box input {{
            width: 100%;
            padding: 0.5rem 0.75rem;
            border: 1px solid var(--border);
            border-radius: 6px;
            font-size: 0.9rem;
            background: white;
            margin-bottom: 0;
        }}
        .combobox-search-box input:focus {{
            outline: none;
            border-color: var(--primary);
            box-shadow: 0 0 0 2px rgba(37, 99, 235, 0.15);
        }}
        .combobox-options {{
            max-height: 280px;
            overflow-y: auto;
            padding: 0.3rem 0;
        }}
        .combobox-item {{
            padding: 0.55rem 0.8rem;
            display: flex;
            align-items: center;
            justify-content: space-between;
            cursor: pointer;
            transition: background 0.1s;
            font-size: 0.88rem;
        }}
        .combobox-item:hover, .combobox-item.highlighted {{
            background: #eff6ff;
        }}
        .combobox-item.selected {{
            background: #dbeafe;
            font-weight: 600;
        }}
        .combobox-item-left {{
            display: flex;
            flex-direction: column;
            gap: 0.1rem;
        }}
        .combobox-item-sp {{
            font-weight: 600;
            color: #0f172a;
        }}
        .combobox-item-common {{
            font-size: 0.8rem;
            color: #64748b;
        }}
        .combobox-item-clade {{
            font-size: 0.72rem;
            padding: 0.15rem 0.45rem;
            border-radius: 4px;
            background: #f1f5f9;
            color: #475569;
            font-weight: 600;
            white-space: nowrap;
        }}
        .combobox-no-results {{
            padding: 1.5rem;
            text-align: center;
            color: #94a3b8;
            font-size: 0.88rem;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.88rem;
            margin-top: 0.5rem;
        }}
        th, td {{
            padding: 0.6rem 0.75rem;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{ background: #f1f5f9; font-weight: 600; color: #334155; }}
        tr:hover td {{ background: #f8fafc; }}

        .btn-group {{
            display: flex;
            gap: 0.35rem;
            flex-wrap: wrap;
            margin-bottom: 0.5rem;
        }}
        .btn {{
            background: #f1f5f9;
            border: 1px solid var(--border);
            padding: 0.32rem 0.65rem;
            border-radius: 6px;
            font-size: 0.81rem;
            font-weight: 600;
            color: #475569;
            cursor: pointer;
            transition: all 0.15s;
            white-space: nowrap;
        }}
        .btn:hover {{ background: #e2e8f0; color: #0f172a; }}
        .btn.active {{ background: var(--primary); color: white; border-color: var(--primary); }}
        .btn-sm {{
            padding: 0.2rem 0.5rem;
            font-size: 0.76rem;
        }}

        .img-container {{
            background: white;
            border: 1px solid var(--border);
            border-radius: 8px;
            overflow: hidden;
            cursor: zoom-in;
            display: flex;
            flex-direction: column;
            transition: all 0.2s ease;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
        }}
        .img-container:hover {{
            border-color: #3b82f6;
            box-shadow: 0 6px 18px rgba(37, 99, 235, 0.15);
        }}
        .img-container img {{
            width: 100%;
            height: auto;
            display: block;
            border-radius: 7px 7px 0 0;
            transition: transform 0.2s;
            pointer-events: none;
        }}
        .img-container:hover img {{
            transform: scale(1.006);
        }}
        .img-caption-bar {{
            padding: 0.32rem 0.75rem;
            background: #f8fafc;
            border-top: 1px solid var(--border);
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.78rem;
            color: #64748b;
        }}
        .caption-text {{
            font-weight: 600;
            color: #334155;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        .caption-action {{
            color: #2563eb;
            font-weight: 600;
            font-size: 0.78rem;
            display: flex;
            align-items: center;
            gap: 0.2rem;
            white-space: nowrap;
        }}
        .img-container:hover .caption-action {{
            color: #1d4ed8;
            text-decoration: underline;
        }}

        /* FULLSCREEN INTERACTIVE LIGHTBOX */
        .lightbox-modal {{
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            width: 100vw;
            height: 100vh;
            background: rgba(10, 15, 30, 0.96);
            backdrop-filter: blur(10px);
            z-index: 999999;
            flex-direction: column;
            overflow: hidden !important;
            user-select: none;
        }}
        .lightbox-modal.active {{
            display: flex;
        }}
        .lightbox-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 0.6rem 1.25rem;
            background: #0f172a;
            border-bottom: 1px solid rgba(255, 255, 255, 0.12);
            color: white;
            flex-shrink: 0;
            z-index: 1000;
            position: relative;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.4);
            gap: 0.75rem;
            flex-wrap: wrap;
        }}
        .lightbox-title {{
            font-weight: 700;
            font-size: 0.95rem;
            color: #f8fafc;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        .lightbox-controls {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
            flex-wrap: wrap;
        }}
        .subpanel-nav {{
            display: flex;
            align-items: center;
            gap: 0.25rem;
            background: rgba(30, 41, 59, 0.75);
            padding: 0.2rem 0.35rem;
            border-radius: 8px;
            border: 1px solid rgba(255, 255, 255, 0.15);
        }}
        .subpanel-lbl {{
            font-size: 0.76rem;
            color: #94a3b8;
            font-weight: 600;
            padding: 0 0.25rem;
            white-space: nowrap;
        }}
        .subpanel-btn {{
            background: rgba(255, 255, 255, 0.08);
            border: 1px solid rgba(255, 255, 255, 0.15);
            color: #e2e8f0;
            padding: 0.25rem 0.45rem;
            border-radius: 5px;
            font-size: 0.75rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
            white-space: nowrap;
        }}
        .subpanel-btn:hover {{
            background: #2563eb;
            color: white;
            border-color: #60a5fa;
        }}
        .subpanel-btn.active {{
            background: #2563eb;
            color: white;
            border-color: #93c5fd;
            box-shadow: 0 0 10px rgba(37, 99, 235, 0.7);
        }}
        .zoom-btn-group {{
            display: flex;
            align-items: center;
            gap: 0.3rem;
        }}
        .tool-btn {{
            background: rgba(255, 255, 255, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.2);
            color: white;
            padding: 0.3rem 0.6rem;
            border-radius: 6px;
            font-size: 0.82rem;
            font-weight: 700;
            cursor: pointer;
            transition: all 0.15s;
        }}
        .tool-btn:hover {{
            background: rgba(255, 255, 255, 0.25);
        }}
        .tool-btn.close-btn {{
            background: #ef4444;
            border-color: #f87171;
            margin-left: 0.3rem;
        }}
        .tool-btn.close-btn:hover {{
            background: #dc2626;
        }}
        .zoom-indicator {{
            font-size: 0.8rem;
            font-weight: 600;
            color: #38bdf8;
            min-width: 52px;
            text-align: center;
        }}
        .lightbox-body {{
            flex: 1 1 0;
            overflow: hidden !important;
            position: relative;
            display: flex;
            align-items: center;
            justify-content: center;
            cursor: grab;
            width: 100%;
            height: 100%;
        }}
        .lightbox-body.grabbing {{
            cursor: grabbing;
        }}
        .lightbox-img-wrapper {{
            position: absolute;
            transform-origin: center center;
            transition: transform 0.28s cubic-bezier(0.25, 1, 0.5, 1);
            will-change: transform;
            display: flex;
            align-items: center;
            justify-content: center;
        }}
        .lightbox-img-wrapper.no-transition {{
            transition: none !important;
        }}
        .lightbox-img-wrapper img {{
            max-width: none;
            max-height: none;
            display: block;
            pointer-events: auto;
            cursor: crosshair;
            box-shadow: 0 10px 40px rgba(0, 0, 0, 0.8);
            border-radius: 4px;
        }}
        .lightbox-footer {{
            padding: 0.45rem 1.25rem;
            background: #0f172a;
            border-top: 1px solid rgba(255, 255, 255, 0.12);
            color: #94a3b8;
            font-size: 0.78rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-shrink: 0;
            z-index: 1000;
            position: relative;
        }}
        .lightbox-footer code {{
            background: rgba(255, 255, 255, 0.12);
            padding: 0.15rem 0.4rem;
            border-radius: 4px;
            color: #e2e8f0;
        }}

        .search-box {{
            display: flex;
            gap: 1rem;
            margin-bottom: 1rem;
        }}
    </style>
</head>
<body>

    <header>
        <div class="header-left">
            <h1>Vertebrate Globin Genomics Explorer</h1>
            <p>Interactive Comparative Atlas, Synteny Tracker &amp; 1-on-1 Architectural Visualizer</p>
        </div>
        <div class="header-right" style="display: flex; align-items: center; gap: 0.65rem;">
            <a href="overview_guide.html" target="_blank" class="btn" style="background: rgba(56, 189, 248, 0.15); border: 1px solid rgba(56, 189, 248, 0.4); color: #38bdf8; font-weight: 600; padding: 0.35rem 0.75rem; font-size: 0.82rem; border-radius: 6px; text-decoration: none; display: inline-flex; align-items: center; gap: 0.35rem; transition: all 0.2s ease;">Overview &amp; Guide ↗</a>
            <span class="badge">68 Vertebrate Genomes • {len(df_intervals)} Curated Loci • {len(df_genes):,} Verified Gene Models • 10 Evolutionary Clades</span>
        </div>
    </header>

    <nav class="tabs">
        <button class="active" onclick="showTab('meta')">Cross-Species Meta Explorer</button>
        <button onclick="showTab('comparison')">1-on-1 Head-to-Head Comparison</button>
        <button onclick="showTab('deepdive')">Single Species Deep Dive</button>
        <button onclick="showTab('genes')">Interval Genes Database</button>
        <button onclick="showTab('stats')">Statistical Contrasts</button>
    </nav>

    <!-- TAB 1: CROSS-SPECIES META EXPLORER (DEFAULT LANDING PAGE) -->
    <div id="tab-meta" class="tab-pane active">
        <div class="card" style="padding: 0.65rem 0.95rem;">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.45rem; flex-wrap: wrap; gap: 0.5rem;">
                <h3 style="margin: 0; font-size: 0.98rem; font-weight: 700;">Global Macroevolutionary Figures</h3>
                <div class="btn-group" id="meta-view-btns" style="margin-bottom: 0;">
                    <button class="btn" onclick="setMetaView('fig0')">Overview &amp; Guide</button>
                    <button class="btn active" onclick="setMetaView('fig1a')">Fig 1A: Curated Summary</button>
                    <button class="btn" onclick="setMetaView('fig1b')">Fig 1B: All 68 Species</button>
                    <button class="btn" onclick="setMetaView('fig2')">Fig 2: Mechanistic Drivers</button>
                    <button class="btn" onclick="setMetaView('fig3')">Fig 3: Synteny &amp; Codon Map</button>
                    <button class="btn" onclick="setMetaView('fig4')">Fig 4: Bernardi Compaction</button>
                </div>
            </div>
            <!-- Subpanel Navigator for Overview Figures -->
            <div id="fig0-subpanel-bar" style="display: none; align-items: center; gap: 0.35rem; padding: 0.35rem 0.6rem; background: #0f172a; border: 1px solid #334155; border-radius: 6px; margin-bottom: 0.45rem; flex-wrap: wrap;">
                <span style="font-size: 0.75rem; font-weight: 700; color: #94a3b8; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 0.25rem;">Overview Figures:</span>
                <button class="btn btn-sm active" id="subbtn-fig0" onclick="setMetaView('fig0')">Fig 0: Dashboard Guide</button>
                <button class="btn btn-sm" id="subbtn-fig0a" onclick="setMetaView('fig0a')">Fig 0A: Human Paradigm</button>
                <button class="btn btn-sm" id="subbtn-fig0b" onclick="setMetaView('fig0b')">Fig 0B: Synteny Pipeline</button>
                <button class="btn btn-sm" id="subbtn-fig0c" onclick="setMetaView('fig0c')">Fig 0C: Metrics &amp; Codon Dynamics</button>
                <a href="overview_guide.html" target="_blank" class="btn btn-sm" style="margin-left: auto; background: rgba(56, 189, 248, 0.18); color: #38bdf8; border: 1px solid rgba(56, 189, 248, 0.4); text-decoration: none; font-weight: 600;">Interactive Slide Navigator ↗</a>
            </div>
            <div class="img-container" onclick="openLightbox(this)">
                <img id="meta-img" src="summary_figure1a_locus_divergence_summary.png" alt="Fig 1A: Curated Macroevolutionary Locus Divergence Across Key Clades &amp; Extrema">
                <div class="img-caption-bar">
                    <span class="caption-text" id="caption-meta">Fig 1A: Curated Macroevolutionary Locus Divergence Across Key Clades &amp; Extrema</span>
                    <span class="caption-action">Inspect &amp; Zoom →</span>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 2: 1-on-1 HEAD-TO-HEAD COMPARISON -->
    <div id="tab-comparison" class="tab-pane">
        <div class="card">
            <h3>Configure 1-on-1 Species Contrast</h3>
            <div class="grid-2">
                <div>
                    <label style="font-size: 0.82rem; font-weight: 600; color: #475569; display: block; margin-bottom: 0.25rem;">SPECIES A</label>
                    <div class="combobox-wrapper" id="wrapper-comp-sp-a">
                        <div class="combobox-trigger" onclick="toggleCombobox('comp-sp-a')" tabindex="0" onkeydown="handleComboboxTriggerKey(event, 'comp-sp-a')">
                            <div class="combobox-selection" id="display-comp-sp-a">
                                <span>Loading species...</span>
                            </div>
                            <div class="combobox-arrow">▼</div>
                        </div>
                        <div class="combobox-dropdown" id="dropdown-comp-sp-a">
                            <div class="combobox-search-box">
                                <input type="text" id="search-comp-sp-a" placeholder="Search species, common name, clade, short name..." oninput="filterCombobox('comp-sp-a', this.value)" onkeydown="handleSearchInputKey(event, 'comp-sp-a')">
                            </div>
                            <div class="combobox-options" id="options-comp-sp-a"></div>
                        </div>
                        <select id="comp-sp-a" style="display: none;" onchange="updateComparison()"></select>
                    </div>
                </div>
                <div>
                    <label style="font-size: 0.82rem; font-weight: 600; color: #475569; display: block; margin-bottom: 0.25rem;">SPECIES B</label>
                    <div class="combobox-wrapper" id="wrapper-comp-sp-b">
                        <div class="combobox-trigger" onclick="toggleCombobox('comp-sp-b')" tabindex="0" onkeydown="handleComboboxTriggerKey(event, 'comp-sp-b')">
                            <div class="combobox-selection" id="display-comp-sp-b">
                                <span>Loading species...</span>
                            </div>
                            <div class="combobox-arrow">▼</div>
                        </div>
                        <div class="combobox-dropdown" id="dropdown-comp-sp-b">
                            <div class="combobox-search-box">
                                <input type="text" id="search-comp-sp-b" placeholder="Search species, common name, clade, short name..." oninput="filterCombobox('comp-sp-b', this.value)" onkeydown="handleSearchInputKey(event, 'comp-sp-b')">
                            </div>
                            <div class="combobox-options" id="options-comp-sp-b"></div>
                        </div>
                        <select id="comp-sp-b" style="display: none;" onchange="updateComparison()"></select>
                    </div>
                </div>
            </div>
        </div>

        <div class="card">
            <h3>Comparative Scorecard &amp; Isochore Polarity</h3>
            <div style="overflow-x: auto;">
                <table id="comp-table">
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th id="th-sp-a">Species A</th>
                            <th id="th-sp-b">Species B</th>
                            <th>Contrast / Difference</th>
                        </tr>
                    </thead>
                    <tbody id="comp-tbody"></tbody>
                </table>
            </div>
        </div>

        <div class="card">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem; flex-wrap: wrap; gap: 0.5rem;">
                <h3 style="margin: 0;">Synchronized Visual Comparison</h3>
                <div class="btn-group" id="comp-view-btns" style="margin-bottom: 0;">
                    <button class="btn active" onclick="setCompView('chromosome')">Linear Chromosome Tracks</button>
                    <button class="btn" onclick="setCompView('inset')">Bivariate with Mini Inset</button>
                    <button class="btn" onclick="setCompView('diagnostic')">Diagnostic GC Distribution</button>
                    <button class="btn" onclick="setCompView('bivariate')">Standard Bivariate</button>
                </div>
            </div>

            <div class="grid-2">
                <div>
                    <h4 id="lbl-img-a" style="margin-bottom: 0.35rem; color: #1e293b; font-size: 0.95rem;">Species A</h4>
                    <div class="img-container" onclick="openLightbox(this)">
                        <img id="img-comp-a" src="" alt="Plot A">
                        <div class="img-caption-bar">
                            <span class="caption-text" id="caption-comp-a">Species A Plot</span>
                            <span class="caption-action">Inspect &amp; Zoom →</span>
                        </div>
                    </div>
                </div>
                <div>
                    <h4 id="lbl-img-b" style="margin-bottom: 0.35rem; color: #1e293b; font-size: 0.95rem;">Species B</h4>
                    <div class="img-container" onclick="openLightbox(this)">
                        <img id="img-comp-b" src="" alt="Plot B">
                        <div class="img-caption-bar">
                            <span class="caption-text" id="caption-comp-b">Species B Plot</span>
                            <span class="caption-action">Inspect &amp; Zoom →</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 3: SINGLE SPECIES DEEP DIVE -->
    <div id="tab-deepdive" class="tab-pane">
        <div class="grid-2" style="grid-template-columns: 320px 1fr;">
            <div>
                <div class="card">
                    <h3>Select Species</h3>
                    <div class="combobox-wrapper" id="wrapper-deep-sp" style="margin-top: 0.35rem;">
                        <div class="combobox-trigger" onclick="toggleCombobox('deep-sp')" tabindex="0" onkeydown="handleComboboxTriggerKey(event, 'deep-sp')">
                            <div class="combobox-selection" id="display-deep-sp">
                                <span>Loading species...</span>
                            </div>
                            <div class="combobox-arrow">▼</div>
                        </div>
                        <div class="combobox-dropdown" id="dropdown-deep-sp">
                            <div class="combobox-search-box">
                                <input type="text" id="search-deep-sp" placeholder="Search species, common name, clade, short name..." oninput="filterCombobox('deep-sp', this.value)" onkeydown="handleSearchInputKey(event, 'deep-sp')">
                            </div>
                            <div class="combobox-options" id="options-deep-sp"></div>
                        </div>
                        <select id="deep-sp" style="display: none;" onchange="updateDeepDive()"></select>
                    </div>
                    
                    <h4 style="margin-top: 0.8rem; margin-bottom: 0.35rem; font-size: 0.9rem; color: #334155;">Species Profile</h4>
                    <div id="deep-profile" style="font-size: 0.85rem; line-height: 1.6; color: #475569;"></div>
                </div>
            </div>
            <div>
                <div class="card">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem; flex-wrap: wrap; gap: 0.5rem;">
                        <h3 id="deep-title" style="margin: 0;">Visualizations</h3>
                        <div class="btn-group" id="deep-view-btns" style="margin-bottom: 0;">
                            <button class="btn active" onclick="setDeepView('chromosome')">Linear Chromosome</button>
                            <button class="btn" onclick="setDeepView('inset')">Bivariate with Inset</button>
                            <button class="btn" onclick="setDeepView('diagnostic')">Diagnostic GC</button>
                            <button class="btn" onclick="setDeepView('bivariate')">Standard Bivariate</button>
                        </div>
                    </div>
                    <div class="img-container" onclick="openLightbox(this)">
                        <img id="deep-img" src="" alt="Deep dive plot">
                        <div class="img-caption-bar">
                            <span class="caption-text" id="caption-deep">Architecture &amp; Distribution Plot</span>
                            <span class="caption-action">Inspect &amp; Zoom →</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 4: INTERVAL GENES DATABASE -->
    <div id="tab-genes" class="tab-pane">
        <div class="card">
            <h3>Annotated Genes Search Database ({len(df_genes):,} Genes)</h3>
            <div class="search-box">
                <input type="text" id="gene-search" placeholder="Search by Gene Symbol, Species, or Clade (e.g. HBA, human, mammal)..." oninput="filterGenes()">
            </div>
            <div style="max-height: 600px; overflow-y: auto;">
                <table id="genes-table">
                    <thead>
                        <tr>
                            <th>Species</th>
                            <th>Gene Symbol</th>
                            <th>Role / Group</th>
                            <th>Clade</th>
                            <th>Chr</th>
                            <th>CDS GC%</th>
                            <th>Gene GC%</th>
                            <th>Z (Gene)</th>
                            <th>Z (CDS)</th>
                        </tr>
                    </thead>
                    <tbody id="genes-tbody"></tbody>
                </table>
            </div>
            <p id="genes-count" style="font-size: 0.8rem; color: var(--text-muted); margin-top: 0.5rem;"></p>
        </div>
    </div>

    <!-- TAB 5: STATISTICAL CONTRASTS -->
    <div id="tab-stats" class="tab-pane">
        <div class="card">
            <h3>Formal Macroevolutionary Hypothesis Tests (Supplementary Table S4)</h3>
            <div style="overflow-x: auto;">
                <table>
                    <thead>
                        <tr>
                            <th>Hypothesis / Contrast</th>
                            <th>Group 1</th>
                            <th>Group 2</th>
                            <th>Delta / Effect</th>
                            <th>Statistical Test</th>
                            <th>p-value</th>
                            <th>Biological Conclusion</th>
                        </tr>
                    </thead>
                    <tbody>
                        {"".join(f'''<tr>
                            <td style="font-weight: 600;">{s.get("hypothesis_or_contrast", "")}</td>
                            <td>{s.get("group_1", "")} (mean={s.get("mean_group_1", "")})</td>
                            <td>{s.get("group_2", "")} (mean={s.get("mean_group_2", "")})</td>
                            <td><b>{s.get("difference_or_delta", "")}</b></td>
                            <td>{s.get("statistical_test", "")} (t/U={s.get("test_statistic", "")})</td>
                            <td><span style="color: #dc2626; font-weight: 700;">{s.get("p_value_formatted", "")}</span></td>
                            <td style="font-size: 0.85rem;">{s.get("biological_conclusion", "")}</td>
                        </tr>''' for s in stats_data)}
                    </tbody>
                </table>
            </div>
        </div>
    </div>

    <!-- FULLSCREEN INTERACTIVE LIGHTBOX MODAL -->
    <div id="lightbox-modal" class="lightbox-modal">
        <div class="lightbox-header">
            <div class="lightbox-title" id="lightbox-title">High-Resolution Figure Inspection</div>
            <div class="lightbox-controls">
                <div id="lightbox-subpanels" class="subpanel-nav" style="display: none;">
                    <span class="subpanel-lbl">Diagnostic Panels:</span>
                    <button class="subpanel-btn active" onclick="zoomToPanel('all')" title="Full 6-panel overview">Full View</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('A')" title="Panel A: Whole-Gene GC% vs Genome Baseline">[A] Gene GC</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('B')" title="Panel B: Coding CDS GC% vs Genome Baseline">[B] CDS GC</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('C')" title="Panel C: Bivariate Gene Length vs GC%">[C] Length vs GC</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('D')" title="Panel D: Whole-Gene GC% vs Coding CDS GC%">[D] Gene vs CDS</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('E')" title="Panel E: Coding CDS GC% vs Exon GC3%">[E] CDS vs GC3</button>
                    <button class="subpanel-btn" onclick="zoomToPanel('F')" title="Panel F: Chromosomal GC Compartments & Synteny">[F] Synteny GC</button>
                </div>
                <div class="zoom-btn-group">
                    <button class="tool-btn" onclick="changeZoom(-0.35)" title="Zoom Out (-)">−</button>
                    <span id="zoom-level-text" class="zoom-indicator">100%</span>
                    <button class="tool-btn" onclick="changeZoom(0.35)" title="Zoom In (+)">+</button>
                    <button class="tool-btn" onclick="fitToScreen()" title="Fit to Screen (Key: 0)">Fit</button>
                    <button class="tool-btn" onclick="setZoom100()" title="100% Native Pixels">1:1</button>
                    <button class="tool-btn close-btn" onclick="closeLightbox()" title="Close (Esc)">&times; Close</button>
                </div>
            </div>
        </div>
        <div class="lightbox-body" id="lightbox-body" onwheel="handleWheel(event)" onmousedown="startDrag(event)" onmousemove="doDrag(event)" onmouseup="endDrag(event)" onmouseleave="endDrag(event)">
            <div class="lightbox-img-wrapper" id="lightbox-wrapper">
                <img id="lightbox-img" src="" alt="Zoomed Plot" draggable="false" onclick="handleImageClick(event)">
            </div>
        </div>
        <div class="lightbox-footer">
            <span id="lightbox-mouse-help"><b>Interactive Zoom:</b> Click to zoom in/out • Scroll wheel to zoom • Drag to pan</span>
            <span id="lightbox-shortcuts"><b>Shortcuts:</b> <code>Esc</code> Close • <code>0</code> Fit • <code>+</code>/<code>-</code> Zoom</span>
        </div>
    </div>

    <script>
        const speciesData = {json.dumps(species_data)};
        const genesData = {json.dumps(genes_data)};
        let currentCompView = 'chromosome';
        let currentDeepView = 'chromosome';

        function showTab(tabId) {{
            document.querySelectorAll('.tab-pane').forEach(el => el.classList.remove('active'));
            document.querySelectorAll('nav.tabs button').forEach(el => el.classList.remove('active'));
            const targetPane = document.getElementById('tab-' + tabId);
            if (targetPane) targetPane.classList.add('active');
            const clickedBtn = (window.event && window.event.target) ? window.event.target.closest('button') : null;
            if (clickedBtn && clickedBtn.closest('nav.tabs')) {{
                clickedBtn.classList.add('active');
            }} else {{
                document.querySelectorAll('nav.tabs button').forEach(btn => {{
                    const oc = btn.getAttribute('onclick') || '';
                    if (oc.indexOf(tabId) !== -1) btn.classList.add('active');
                }});
            }}
        }}

        /* COMBOBOX CONTROLLER & DYNAMIC SEARCH */
        function toggleCombobox(id) {{
            const wrapper = document.getElementById('wrapper-' + id);
            if (!wrapper) return;
            const isOpen = wrapper.classList.contains('open');
            closeAllComboboxes();
            if (!isOpen) {{
                wrapper.classList.add('open');
                const searchInput = document.getElementById('search-' + id);
                if (searchInput) {{
                    searchInput.value = '';
                    filterCombobox(id, '');
                    setTimeout(function() {{
                        searchInput.focus();
                    }}, 40);
                }}
            }}
        }}

        function closeAllComboboxes() {{
            document.querySelectorAll('.combobox-wrapper').forEach(function(w) {{
                w.classList.remove('open');
            }});
        }}

        document.addEventListener('click', function(e) {{
            if (!e.target.closest('.combobox-wrapper')) {{
                closeAllComboboxes();
            }}
        }});

        function filterCombobox(id, query) {{
            const optionsContainer = document.getElementById('options-' + id);
            if (!optionsContainer) return;
            const currentVal = document.getElementById(id).value;
            const q = (query || '').toLowerCase().trim();
            const words = q ? q.split(/\s+/) : [];
            const spKeys = Object.keys(speciesData).sort();

            let matched = [];
            spKeys.forEach(function(sp) {{
                const d = speciesData[sp];
                const sciName = sp.replace(/_/g, ' ');
                const parts = sp.split('_');
                const shortName = parts.length > 1 ? (parts[0][0] + '. ' + parts.slice(1).join(' ')) : sp;
                const searchCorpus = (sp + ' ' + sciName + ' ' + shortName + ' ' + (d.common_name || '') + ' ' + (d.clade || '')).toLowerCase();

                let isMatch = true;
                for (let i = 0; i < words.length; i++) {{
                    if (!searchCorpus.includes(words[i])) {{
                        isMatch = false;
                        break;
                    }}
                }}
                if (isMatch) {{
                    matched.push({{ sp: sp, d: d, shortName: shortName, sciName: sciName }});
                }}
            }});

            optionsContainer.innerHTML = '';

            if (matched.length === 0) {{
                const noRes = document.createElement('div');
                noRes.className = 'combobox-no-results';
                noRes.textContent = 'No species match "' + query + '"';
                optionsContainer.appendChild(noRes);
                return;
            }}

            matched.forEach(function(item, idx) {{
                const itemEl = document.createElement('div');
                itemEl.className = 'combobox-item' + (item.sp === currentVal ? ' selected' : '') + (idx === 0 ? ' highlighted' : '');
                itemEl.setAttribute('data-value', item.sp);
                itemEl.onclick = function() {{
                    selectSpecies(id, item.sp);
                }};

                const leftDiv = document.createElement('div');
                leftDiv.className = 'combobox-item-left';

                const spSpan = document.createElement('span');
                spSpan.className = 'combobox-item-sp';
                spSpan.textContent = item.sciName + ' ';

                const shortSpan = document.createElement('span');
                shortSpan.style.fontWeight = 'normal';
                shortSpan.style.color = '#94a3b8';
                shortSpan.style.fontSize = '0.78rem';
                shortSpan.textContent = '(' + item.shortName + ')';
                spSpan.appendChild(shortSpan);

                const commonSpan = document.createElement('span');
                commonSpan.className = 'combobox-item-common';
                commonSpan.textContent = item.d.common_name || item.sp;

                leftDiv.appendChild(spSpan);
                leftDiv.appendChild(commonSpan);

                const cladeSpan = document.createElement('span');
                cladeSpan.className = 'combobox-item-clade';
                cladeSpan.textContent = item.d.clade || '';

                itemEl.appendChild(leftDiv);
                itemEl.appendChild(cladeSpan);
                optionsContainer.appendChild(itemEl);
            }});
        }}

        function selectSpecies(id, sp, triggerUpdate) {{
            if (triggerUpdate === undefined) triggerUpdate = true;
            const selectEl = document.getElementById(id);
            if (!selectEl) return;
            selectEl.value = sp;

            const d = speciesData[sp];
            const displayEl = document.getElementById('display-' + id);
            if (displayEl && d) {{
                displayEl.innerHTML = '';
                const sciName = sp.replace(/_/g, ' ');

                const nameSpan = document.createElement('span');
                nameSpan.className = 'combobox-name';
                nameSpan.textContent = sciName;

                const commonSpan = document.createElement('span');
                commonSpan.className = 'combobox-common';
                commonSpan.textContent = '(' + (d.common_name || sp) + ')';

                const cladeSpan = document.createElement('span');
                cladeSpan.className = 'combobox-clade';
                cladeSpan.textContent = d.clade || '';

                displayEl.appendChild(nameSpan);
                displayEl.appendChild(commonSpan);
                displayEl.appendChild(cladeSpan);
            }}

            closeAllComboboxes();

            if (triggerUpdate) {{
                if (id === 'comp-sp-a' || id === 'comp-sp-b') {{
                    updateComparison();
                }} else if (id === 'deep-sp') {{
                    updateDeepDive();
                }}
            }}
        }}

        function handleComboboxTriggerKey(e, id) {{
            if (e.key === 'Enter' || e.key === ' ' || e.key === 'ArrowDown') {{
                e.preventDefault();
                toggleCombobox(id);
            }}
        }}

        function handleSearchInputKey(e, id) {{
            const optionsContainer = document.getElementById('options-' + id);
            if (!optionsContainer) return;
            const items = optionsContainer.querySelectorAll('.combobox-item');
            if (!items.length) return;

            let curIndex = -1;
            items.forEach(function(item, idx) {{
                if (item.classList.contains('highlighted')) curIndex = idx;
            }});

            if (e.key === 'ArrowDown') {{
                e.preventDefault();
                if (curIndex >= 0) items[curIndex].classList.remove('highlighted');
                curIndex = (curIndex + 1) % items.length;
                items[curIndex].classList.add('highlighted');
                items[curIndex].scrollIntoView({{ block: 'nearest' }});
            }} else if (e.key === 'ArrowUp') {{
                e.preventDefault();
                if (curIndex >= 0) items[curIndex].classList.remove('highlighted');
                curIndex = (curIndex - 1 + items.length) % items.length;
                items[curIndex].classList.add('highlighted');
                items[curIndex].scrollIntoView({{ block: 'nearest' }});
            }} else if (e.key === 'Enter') {{
                e.preventDefault();
                const target = curIndex >= 0 ? items[curIndex] : items[0];
                if (target) {{
                    const val = target.getAttribute('data-value');
                    selectSpecies(id, val);
                }}
            }} else if (e.key === 'Escape') {{
                e.preventDefault();
                closeAllComboboxes();
                const trigger = document.getElementById('wrapper-' + id).querySelector('.combobox-trigger');
                if (trigger) trigger.focus();
            }}
        }}

        function populateSpeciesDropdowns() {{
            const spKeys = Object.keys(speciesData).sort();
            const selA = document.getElementById('comp-sp-a');
            const selB = document.getElementById('comp-sp-b');
            const selDeep = document.getElementById('deep-sp');

            spKeys.forEach(sp => {{
                const d = speciesData[sp];
                const optA = new Option(`${{sp}} (${{d.common_name}} • ${{d.clade}})`, sp);
                const optB = new Option(`${{sp}} (${{d.common_name}} • ${{d.clade}})`, sp);
                const optDeep = new Option(`${{sp}} (${{d.common_name}} • ${{d.clade}})`, sp);
                selA.add(optA);
                selB.add(optB);
                selDeep.add(optDeep);
            }});

            selectSpecies('comp-sp-a', 'Homo_sapiens', false);
            selectSpecies('comp-sp-b', 'Gallus_gallus', false);
            selectSpecies('deep-sp', 'Homo_sapiens', false);

            filterCombobox('comp-sp-a', '');
            filterCombobox('comp-sp-b', '');
            filterCombobox('deep-sp', '');

            updateComparison();
            updateDeepDive();
            renderGenesTable(genesData.slice(0, 100));
        }}

        function getViewDisplayName(viewType) {{
            if (viewType === 'chromosome') return 'Linear Chromosome Tracks';
            if (viewType === 'inset') return 'Bivariate with Mini Inset';
            if (viewType === 'diagnostic') return '6-Panel Diagnostic GC Distribution';
            if (viewType === 'bivariate') return 'Standard Bivariate Plot';
            return 'Visual Plot';
        }}

        function getPlotPath(species, plotType) {{
            if (plotType === 'chromosome') return `chromosome_plots/chromosome_linear_${{species}}.png`;
            if (plotType === 'inset') return `bivariate_with_inset_plots/bivariate_with_inset_${{species}}.png`;
            if (plotType === 'diagnostic') return `diagnostic_plots/diagnostic_${{species}}.png`;
            if (plotType === 'bivariate') return `bivariate_length_gc_plots/bivariate_length_gc_${{species}}.png`;
            return '';
        }}

        function setCompView(viewType) {{
            currentCompView = viewType;
            document.querySelectorAll('#comp-view-btns button').forEach(b => b.classList.remove('active'));
            const clickedBtn = (window.event && window.event.target) ? window.event.target.closest('button') : null;
            if (clickedBtn && clickedBtn.closest('#comp-view-btns')) {{
                clickedBtn.classList.add('active');
            }} else {{
                document.querySelectorAll('#comp-view-btns button').forEach(btn => {{
                    const oc = btn.getAttribute('onclick') || '';
                    if (oc.indexOf(viewType) !== -1) btn.classList.add('active');
                }});
            }}
            updateComparisonImages();
        }}

        function updateComparisonImages() {{
            const spA = document.getElementById('comp-sp-a').value;
            const spB = document.getElementById('comp-sp-b').value;
            document.getElementById('lbl-img-a').innerText = `${{spA}} (${{speciesData[spA]?.common_name || ''}})`;
            document.getElementById('lbl-img-b').innerText = `${{spB}} (${{speciesData[spB]?.common_name || ''}})`;
            document.getElementById('img-comp-a').src = getPlotPath(spA, currentCompView);
            document.getElementById('img-comp-b').src = getPlotPath(spB, currentCompView);

            const vName = getViewDisplayName(currentCompView);
            const capA = document.getElementById('caption-comp-a');
            const capB = document.getElementById('caption-comp-b');
            if (capA) capA.innerText = `${{spA.replace(/_/g, ' ')}} • ${{vName}}`;
            if (capB) capB.innerText = `${{spB.replace(/_/g, ' ')}} • ${{vName}}`;
        }}

        function updateComparison() {{
            const spA = document.getElementById('comp-sp-a').value;
            const spB = document.getElementById('comp-sp-b').value;
            const dA = speciesData[spA];
            const dB = speciesData[spB];

            document.getElementById('th-sp-a').innerText = `${{spA}} (${{dA.common_name}})`;
            document.getElementById('th-sp-b').innerText = `${{spB}} (${{dB.common_name}})`;

            const rows = [
                ['Phylogenetic Clade', dA.clade, dB.clade, dA.clade === dB.clade ? 'Same Clade' : 'Cross-Clade'],
                ['Assembly Accession', dA.assembly, dB.assembly, '—'],
                ['Baseline Coding GC (%)', `${{dA.mean_cds_gc.toFixed(2)}}% (±${{dA.std_cds_gc.toFixed(2)}}%)`, `${{dB.mean_cds_gc.toFixed(2)}}% (±${{dB.std_cds_gc.toFixed(2)}}%)`, `${{(dA.mean_cds_gc - dB.mean_cds_gc).toFixed(2)}}%`],
                ['Curated Globin Loci', `${{dA.loci_count}} loci`, `${{dB.loci_count}} loci`, '—'],
                ['α-Globin GC / Whole-Gene Z', dA.alpha_gc ? `${{dA.alpha_gc.toFixed(1)}}% (Z=${{dA.alpha_z > 0 ? '+' : ''}}${{dA.alpha_z.toFixed(2)}}σ)` : (dA.proto_gc ? `Proto: ${{dA.proto_gc.toFixed(1)}}% (Z=${{dA.proto_z > 0 ? '+' : ''}}${{dA.proto_z.toFixed(2)}}σ)` : 'N/A'), dB.alpha_gc ? `${{dB.alpha_gc.toFixed(1)}}% (Z=${{dB.alpha_z > 0 ? '+' : ''}}${{dB.alpha_z.toFixed(2)}}σ)` : (dB.proto_gc ? `Proto: ${{dB.proto_gc.toFixed(1)}}% (Z=${{dB.proto_z > 0 ? '+' : ''}}${{dB.proto_z.toFixed(2)}}σ)` : 'N/A'), (dA.alpha_z && dB.alpha_z) ? `ΔZ = ${{ (dA.alpha_z - dB.alpha_z).toFixed(2) }}σ` : 'N/A'],
                ['β-Globin GC / Whole-Gene Z', dA.beta_gc ? `${{dA.beta_gc.toFixed(1)}}% (Z=${{dA.beta_z > 0 ? '+' : ''}}${{dA.beta_z.toFixed(2)}}σ)` : (dA.proto_gc ? 'N/A (Ancestral Proto-Globin)' : 'N/A'), dB.beta_gc ? `${{dB.beta_gc.toFixed(1)}}% (Z=${{dB.beta_z > 0 ? '+' : ''}}${{dB.beta_z.toFixed(2)}}σ)` : (dB.proto_gc ? 'N/A (Ancestral Proto-Globin)' : 'N/A'), (dA.beta_z && dB.beta_z) ? `ΔZ = ${{ (dA.beta_z - dB.beta_z).toFixed(2) }}σ` : 'N/A'],
                ['Subunit Polarity (ΔZ = α - β)', dA.polarity ? `${{dA.polarity > 0 ? '+' : ''}}${{dA.polarity.toFixed(2)}}σ` : (dA.proto_gc ? 'N/A (Monomeric Proto-Globin)' : 'N/A'), dB.polarity ? `${{dB.polarity > 0 ? '+' : ''}}${{dB.polarity.toFixed(2)}}σ` : (dB.proto_gc ? 'N/A (Monomeric Proto-Globin)' : 'N/A'), (dA.polarity && dB.polarity) ? ((dA.polarity > 0 && dB.polarity < 0) || (dA.polarity < 0 && dB.polarity > 0) ? 'Polarity Inversion' : 'Conserved Polarity') : 'N/A']
            ];

            const tbody = document.getElementById('comp-tbody');
            tbody.innerHTML = '';
            rows.forEach(r => {{
                const tr = document.createElement('tr');
                tr.innerHTML = `<td style="font-weight: 600;">${{r[0]}}</td><td>${{r[1]}}</td><td>${{r[2]}}</td><td style="font-weight: 600; color: #2563eb;">${{r[3]}}</td>`;
                tbody.appendChild(tr);
            }});

            updateComparisonImages();
        }}

        function setDeepView(viewType) {{
            currentDeepView = viewType;
            document.querySelectorAll('#deep-view-btns button').forEach(b => b.classList.remove('active'));
            const clickedBtn = (window.event && window.event.target) ? window.event.target.closest('button') : null;
            if (clickedBtn && clickedBtn.closest('#deep-view-btns')) {{
                clickedBtn.classList.add('active');
            }} else {{
                document.querySelectorAll('#deep-view-btns button').forEach(btn => {{
                    const oc = btn.getAttribute('onclick') || '';
                    if (oc.indexOf(viewType) !== -1) btn.classList.add('active');
                }});
            }}
            const sp = document.getElementById('deep-sp').value;
            document.getElementById('deep-img').src = getPlotPath(sp, currentDeepView);
            const cap = document.getElementById('caption-deep');
            if (cap) cap.innerText = `${{sp.replace(/_/g, ' ')}} • ${{getViewDisplayName(currentDeepView)}}`;
        }}

        function updateDeepDive() {{
            const sp = document.getElementById('deep-sp').value;
            const d = speciesData[sp];
            document.getElementById('deep-title').innerText = `${{sp}} (${{d.common_name}}) Architecture & Distribution`;
            document.getElementById('deep-img').src = getPlotPath(sp, currentDeepView);
            const cap = document.getElementById('caption-deep');
            if (cap) cap.innerText = `${{sp.replace(/_/g, ' ')}} • ${{getViewDisplayName(currentDeepView)}}`;

            let globinHtml = '';
            if (d.alpha_gc) {{
                globinHtml = `
                    <div><b>α-Globin GC:</b> ${{d.alpha_gc.toFixed(1)}}% (Z=${{d.alpha_z > 0 ? '+' : ''}}${{d.alpha_z.toFixed(2)}}σ)</div>
                    <div><b>β-Globin GC:</b> ${{d.beta_gc.toFixed(1)}}% (Z=${{d.beta_z > 0 ? '+' : ''}}${{d.beta_z.toFixed(2)}}σ)</div>
                    <div><b>Subunit Polarity:</b> ${{d.polarity > 0 ? '+' : ''}}${{d.polarity.toFixed(2)}}σ</div>
                `;
            }} else if (d.proto_gc) {{
                globinHtml = `<div><b>Proto-Globin GC:</b> ${{d.proto_gc.toFixed(1)}}% (Z=${{d.proto_z > 0 ? '+' : ''}}${{d.proto_z.toFixed(2)}}σ)</div>`;
            }}

            document.getElementById('deep-profile').innerHTML = `
                <div><b>Common Name:</b> ${{d.common_name}}</div>
                <div><b>Clade:</b> <code>${{d.clade}}</code></div>
                <div><b>Assembly:</b> <code>${{d.assembly}}</code></div>
                <div><b>Baseline CDS GC:</b> ${{d.mean_cds_gc.toFixed(2)}}% (±${{d.std_cds_gc.toFixed(2)}}%)</div>
                <div><b>Total Coding Genes:</b> ${{typeof d.num_genes === 'number' ? d.num_genes.toLocaleString() : d.num_genes}}</div>
                <div><b>Curated Loci:</b> ${{d.loci_count}} intervals</div>
                ${{globinHtml}}
            `;
        }}

        function setMetaView(figType) {{
            const subBar = document.getElementById('fig0-subpanel-bar');
            if (figType.startsWith('fig0') || figType === 'overview') {{
                if (subBar) subBar.style.display = 'flex';
                if (subBar) {{
                    subBar.querySelectorAll('button').forEach(b => b.classList.remove('active'));
                    const activeSubBtn = document.getElementById('subbtn-' + figType) || document.getElementById('subbtn-fig0');
                    if (activeSubBtn) activeSubBtn.classList.add('active');
                }}
                document.querySelectorAll('#meta-view-btns button').forEach(b => {{
                    const oc = b.getAttribute('onclick') || '';
                    if (oc.indexOf('fig0') !== -1) b.classList.add('active');
                    else b.classList.remove('active');
                }});
            }} else {{
                if (subBar) subBar.style.display = 'none';
                document.querySelectorAll('#meta-view-btns button').forEach(b => b.classList.remove('active'));
                const clickedBtn = (window.event && window.event.target) ? window.event.target.closest('button') : null;
                if (clickedBtn && clickedBtn.closest('#meta-view-btns')) {{
                    clickedBtn.classList.add('active');
                }} else {{
                    document.querySelectorAll('#meta-view-btns button').forEach(btn => {{
                        const oc = btn.getAttribute('onclick') || '';
                        if (oc.indexOf(figType) !== -1) btn.classList.add('active');
                    }});
                }}
            }}

            const img = document.getElementById('meta-img');
            let figName = 'Fig 1A: Curated Macroevolutionary Locus Divergence Across Key Clades & Extrema';
            if (figType === 'fig0' || figType === 'overview' || figType === 'fig0d') {{
                img.src = 'summary_figure0_dashboard_guide.png';
                figName = 'Figure 0: Interactive Dashboard Architecture & Analytical User Guide';
            }} else if (figType === 'fig0a') {{
                img.src = 'summary_figure0a_human_paradigm.png';
                figName = 'Figure 0A: The Human Globin Paradigm & Evolutionary Inquiry';
            }} else if (figType === 'fig0b') {{
                img.src = 'summary_figure0b_synteny_pipeline.png';
                figName = 'Figure 0B: Synteny-Anchored Locus Retrieval Pipeline & 4 Archetypes';
            }} else if (figType === 'fig0c') {{
                img.src = 'summary_figure0c_metric_anatomy_and_gc3.png';
                figName = 'Figure 0C: Anatomy of Calculated Metrics: Gene Architecture & Codon Dynamics';
            }} else if (figType === 'fig1a' || figType === 'fig1') {{
                img.src = 'summary_figure1a_locus_divergence_summary.png';
                figName = 'Fig 1A: Curated Macroevolutionary Locus Divergence Across Key Clades & Extrema';
            }} else if (figType === 'fig1b') {{
                img.src = 'summary_figure1b_locus_divergence_full.png';
                figName = 'Fig 1B: Complete Macroevolutionary Locus Divergence Across All 68 Vertebrates';
            }} else if (figType === 'fig2') {{
                img.src = 'summary_figure2_mechanistic_drivers.png';
                figName = 'Fig 2: Mechanistic Drivers of GC Polarity';
            }} else if (figType === 'fig3') {{
                img.src = 'summary_figure3_phylogenetic_synteny.png';
                figName = 'Fig 3: Macroevolutionary Synteny Architecture and Codon Divergence Across Vertebrates';
            }} else if (figType === 'fig4' || figType === 'bernardi') {{
                img.src = 'summary_figure_bernardi_compaction_law.png';
                figName = "Fig 4: Empirical Proof of Bernardi's Isochore Compaction Law Across Vertebrates";
            }}
            const cap = document.getElementById('caption-meta');
            if (cap) cap.innerText = figName;
        }}

        function renderGenesTable(data) {{
            const tbody = document.getElementById('genes-tbody');
            tbody.innerHTML = '';
            data.forEach(g => {{
                const tr = document.createElement('tr');
                tr.innerHTML = `
                    <td>${{g.species}}</td>
                    <td style="font-weight: 700;">${{g.symbol}}</td>
                    <td><span class="badge" style="background:#e0e7ff; color:#3730a3;">${{g.group}}</span></td>
                    <td>${{g.clade}}</td>
                    <td>${{g.chrom}}</td>
                    <td>${{g.cds_gc}}%</td>
                    <td>${{g.gene_gc}}%</td>
                    <td style="font-weight: 600;">${{g.z_gene > 0 ? '+' : ''}}${{g.z_gene}}σ</td>
                    <td>${{g.z_cds > 0 ? '+' : ''}}${{g.z_cds}}σ</td>
                `;
                tbody.appendChild(tr);
            }});
            document.getElementById('genes-count').innerText = `Showing ${{data.length}} of ${{genesData.length}} genes`;
        }}

        function filterGenes() {{
            const query = document.getElementById('gene-search').value.toLowerCase();
            const filtered = genesData.filter(g => 
                g.symbol.toLowerCase().includes(query) || 
                g.species.toLowerCase().includes(query) || 
                g.clade.toLowerCase().includes(query) ||
                g.group.toLowerCase().includes(query)
            );
            renderGenesTable(filtered.slice(0, 100));
        }}

        /* --- INTERACTIVE PAN & ZOOM LIGHTBOX LOGIC --- */
        let zoomScale = 1.0;
        let baseWidth = 800;
        let baseHeight = 800;
        let panX = 0;
        let panY = 0;
        let isDragging = false;
        let startDragX = 0;
        let startDragY = 0;
        let dragDistance = 0;
        let activePanel = 'all';
        let currentLightboxType = '';

        function openLightbox(container) {{
            const img = container.querySelector('img');
            if (!img || !img.src) return;

            const modal = document.getElementById('lightbox-modal');
            const lbImg = document.getElementById('lightbox-img');
            const titleEl = document.getElementById('lightbox-title');
            const subpanelNav = document.getElementById('lightbox-subpanels');
            const mouseHelp = document.getElementById('lightbox-mouse-help');
            const shortcuts = document.getElementById('lightbox-shortcuts');

            // Find which species and view this image represents
            let sp = '';
            if (img.id === 'img-comp-a') {{
                sp = document.getElementById('comp-sp-a')?.value || '';
            }} else if (img.id === 'img-comp-b') {{
                sp = document.getElementById('comp-sp-b')?.value || '';
            }} else if (img.id === 'deep-img') {{
                sp = document.getElementById('deep-sp')?.value || '';
            }}

            const isDiagnostic = img.src.includes('diagnostic_') || img.src.includes('/diagnostic/');
            currentLightboxType = isDiagnostic ? 'diagnostic' : 'standard';

            let titleHtml = '';
            if (sp && speciesData[sp]) {{
                const d = speciesData[sp];
                const sciName = sp.replace(/_/g, ' ');
                const viewName = isDiagnostic ? '6-Panel Diagnostic GC Distribution' : (img.alt || 'High-Resolution Visual');
                titleHtml = '<span><b>' + sciName + '</b> <span style="color: #94a3b8; font-weight: normal;">(' + d.common_name + ' • ' + d.clade + ')</span> — <span style="color: #38bdf8; font-weight: 600;">' + viewName + '</span></span>';
            }} else if (img.id === 'meta-img') {{
                const cap = document.getElementById('caption-meta');
                const metaTitle = cap ? cap.innerText : 'Macroevolutionary Meta Analysis';
                titleHtml = '<span><b>Macroevolutionary Meta Analysis</b> — <span style="color: #38bdf8; font-weight: 600;">' + metaTitle + '</span></span>';
            }} else {{
                titleHtml = '<span><b>High-Resolution Figure Inspection</b></span>';
            }}

            titleEl.innerHTML = titleHtml;
            subpanelNav.style.display = isDiagnostic ? 'flex' : 'none';

            if (isDiagnostic) {{
                if (mouseHelp) mouseHelp.innerHTML = '<b>Diagnostic View:</b> Click any panel to focus • Click again to return to full view • Drag to pan';
                if (shortcuts) shortcuts.innerHTML = '<b>Shortcuts:</b> <code>Esc</code> Close • <code>0</code> Fit • <code>+</code>/<code>-</code> Zoom • <code>1</code>–<code>6</code> Panels A–F';
            }} else {{
                if (mouseHelp) mouseHelp.innerHTML = '<b>Interactive Zoom:</b> Click to zoom in/out • Scroll wheel to zoom • Drag to pan';
                if (shortcuts) shortcuts.innerHTML = '<b>Shortcuts:</b> <code>Esc</code> Close • <code>0</code> Fit • <code>+</code>/<code>-</code> Zoom';
            }}

            modal.classList.add('active');
            document.body.style.overflow = 'hidden';

            lbImg.src = img.src;
            activePanel = 'all';
            updateSubpanelButtons();

            const triggerFit = function() {{
                fitToScreen();
            }};
            lbImg.onload = triggerFit;
            if (window.requestAnimationFrame) {{
                window.requestAnimationFrame(triggerFit);
            }}
            setTimeout(triggerFit, 50);
            setTimeout(triggerFit, 150);
        }}

        function closeLightbox() {{
            const modal = document.getElementById('lightbox-modal');
            modal.classList.remove('active');
            document.body.style.overflow = '';
        }}

        function applyTransform(withTransition = true) {{
            const wrapper = document.getElementById('lightbox-wrapper');
            if (!wrapper) return;
            if (withTransition) {{
                wrapper.classList.remove('no-transition');
            }} else {{
                wrapper.classList.add('no-transition');
            }}
            wrapper.style.transform = `translate(${{panX}}px, ${{panY}}px) scale(${{zoomScale}})`;

            let zoomText = '';
            if (activePanel && activePanel !== 'all') {{
                zoomText = `Panel ${{activePanel}} (${{zoomScale.toFixed(1)}}x)`;
            }} else if (Math.abs(zoomScale - 1.0) < 0.05) {{
                zoomText = 'Fit (1.0x)';
            }} else {{
                zoomText = `${{zoomScale.toFixed(1)}}x`;
            }}
            document.getElementById('zoom-level-text').innerText = zoomText;
        }}

        function fitToScreen() {{
            const body = document.getElementById('lightbox-body');
            const img = document.getElementById('lightbox-img');
            if (!body || !img) return;

            const bw = Math.max(300, (body.clientWidth > 50 ? body.clientWidth : window.innerWidth) - 50);
            const bh = Math.max(300, (body.clientHeight > 50 ? body.clientHeight : (window.innerHeight - 130)) - 50);
            const iw = img.naturalWidth || 6600;
            const ih = img.naturalHeight || 6600;

            const ratio = Math.min(bw / iw, bh / ih);
            baseWidth = Math.max(100, Math.round(iw * ratio));
            baseHeight = Math.max(100, Math.round(ih * ratio));

            // Setting explicit CSS width/height prevents GPU texture tiling artifacts & dropped corners
            img.style.width = baseWidth + 'px';
            img.style.height = baseHeight + 'px';

            zoomScale = 1.0;
            panX = 0;
            panY = 0;
            activePanel = 'all';
            updateSubpanelButtons();
            applyTransform(true);
        }}

        function setZoom100() {{
            const img = document.getElementById('lightbox-img');
            const iw = img?.naturalWidth || 6600;
            zoomScale = baseWidth > 0 ? (iw / baseWidth) : 2.5;
            panX = 0;
            panY = 0;
            activePanel = '';
            updateSubpanelButtons();
            applyTransform(true);
        }}

        function changeZoom(delta) {{
            const factor = delta > 0 ? 1.25 : 0.8;
            zoomScale = Math.max(0.6, Math.min(8.0, zoomScale * factor));
            applyTransform(true);
        }}

        function zoomToPanel(panel) {{
            if (panel === 'all') {{
                fitToScreen();
                return;
            }}

            const body = document.getElementById('lightbox-body');
            const img = document.getElementById('lightbox-img');
            if (!body || !img) return;

            const bw = Math.max(300, (body.clientWidth > 50 ? body.clientWidth : window.innerWidth) - 50);
            const bh = Math.max(300, (body.clientHeight > 50 ? body.clientHeight : (window.innerHeight - 130)) - 50);

            // Accurate subpanel empirical centers
            let cx = 0.5, cy = 0.5;
            if (panel === 'A') {{ cx = 0.262; cy = 0.164; }}
            else if (panel === 'B') {{ cx = 0.753; cy = 0.164; }}
            else if (panel === 'C') {{ cx = 0.262; cy = 0.466; }}
            else if (panel === 'D') {{ cx = 0.753; cy = 0.466; }}
            else if (panel === 'E') {{ cx = 0.262; cy = 0.795; }}
            else if (panel === 'F') {{ cx = 0.753; cy = 0.795; }}

            // Each panel occupies ~44% of baseWidth and ~33% of baseHeight
            // Target zoom comfortably fits the chosen panel to ~88% of screen
            const targetZoomX = (bw * 0.88) / (baseWidth * 0.44);
            const targetZoomY = (bh * 0.88) / (baseHeight * 0.33);
            zoomScale = Math.min(targetZoomX, targetZoomY);

            panX = (0.5 - cx) * baseWidth * zoomScale;
            panY = (0.5 - cy) * baseHeight * zoomScale;
            activePanel = panel;
            updateSubpanelButtons();
            applyTransform(true);
        }}

        function updateSubpanelButtons() {{
            document.querySelectorAll('.subpanel-btn').forEach(btn => {{
                btn.classList.remove('active');
                if (btn.getAttribute('onclick')?.includes(`'${{activePanel}}'`)) {{
                    btn.classList.add('active');
                }}
            }});
        }}

        function handleWheel(e) {{
            e.preventDefault();
            const zoomFactor = e.deltaY < 0 ? 1.15 : 0.87;
            const newScale = Math.max(0.6, Math.min(10.0, zoomScale * zoomFactor));

            const rect = document.getElementById('lightbox-body').getBoundingClientRect();
            const mouseX = e.clientX - rect.left - rect.width / 2;
            const mouseY = e.clientY - rect.top - rect.height / 2;

            panX = mouseX - (mouseX - panX) * (newScale / zoomScale);
            panY = mouseY - (mouseY - panY) * (newScale / zoomScale);
            zoomScale = newScale;

            applyTransform(false);
        }}

        function startDrag(e) {{
            if (e.button !== 0) return;
            isDragging = true;
            dragDistance = 0;
            startDragX = e.clientX - panX;
            startDragY = e.clientY - panY;
            document.getElementById('lightbox-body').classList.add('grabbing');
        }}

        function doDrag(e) {{
            if (!isDragging) return;
            dragDistance += Math.abs(e.movementX) + Math.abs(e.movementY);
            panX = e.clientX - startDragX;
            panY = e.clientY - startDragY;
            applyTransform(false);
        }}

        function endDrag(e) {{
            if (isDragging) {{
                isDragging = false;
                document.getElementById('lightbox-body').classList.remove('grabbing');
            }}
        }}

        function handleImageClick(e) {{
            if (dragDistance > 6) return; // User was dragging, not clicking

            const rect = e.target.getBoundingClientRect();
            const clickX = (e.clientX - rect.left) / rect.width;
            const clickY = (e.clientY - rect.top) / rect.height;

            if (currentLightboxType === 'diagnostic') {{
                // If already zoomed into a subpanel, click toggles back to full view
                if (activePanel !== 'all') {{
                    zoomToPanel('all');
                    return;
                }}

                // In full view, click zooms into clicked subpanel
                let targetP = 'A';
                if (clickX < 0.50) {{
                    if (clickY < 0.30) targetP = 'A';
                    else if (clickY < 0.63) targetP = 'C';
                    else targetP = 'E';
                }} else {{
                    if (clickY < 0.30) targetP = 'B';
                    else if (clickY < 0.63) targetP = 'D';
                    else targetP = 'F';
                }}
                zoomToPanel(targetP);
                return;
            }}

            // Standard plot zoom toggle (linear chromosome, bivariate, meta)
            if (zoomScale <= 1.3) {{
                zoomScale = 2.2;
                panX = (0.5 - clickX) * baseWidth * zoomScale;
                panY = (0.5 - clickY) * baseHeight * zoomScale;
                applyTransform(true);
            }} else {{
                fitToScreen();
            }}
        }}

        document.addEventListener('keydown', function(e) {{
            const modal = document.getElementById('lightbox-modal');
            if (!modal || !modal.classList.contains('active')) return;
            if (e.key === 'Escape') closeLightbox();
            if (e.key === '+' || e.key === '=') changeZoom(0.35);
            if (e.key === '-' || e.key === '_') changeZoom(-0.35);
            if (e.key === '0') fitToScreen();

            // 1-6 shortcuts only active for diagnostic plots
            if (currentLightboxType === 'diagnostic') {{
                if (e.key === '1') zoomToPanel('A');
                if (e.key === '2') zoomToPanel('B');
                if (e.key === '3') zoomToPanel('C');
                if (e.key === '4') zoomToPanel('D');
                if (e.key === '5') zoomToPanel('E');
                if (e.key === '6') zoomToPanel('F');
            }}
        }});

        window.onload = populateSpeciesDropdowns;
    </script>
</body>
</html>
"""

    out_file = PROJECT_ROOT / "results" / "dashboard.html"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    alt_file = PROJECT_ROOT / "1_data_from_mike" / "results" / "dashboard.html"
    if alt_file.parent.exists():
        with open(alt_file, "w", encoding="utf-8") as f:
            f.write(html_content)

    print(f"✓ Standalone HTML dashboard successfully created at: {out_file}")

if __name__ == "__main__":
    main()
