"""Generates a standalone, interactive HTML dashboard (results/dashboard.html)
allowing instant exploration and 1-on-1 head-to-head comparison across all 44 vertebrate species.
Requires zero server or Python runtime to view—can be double-clicked in any browser.
"""

import json
import os
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent

def main():
    print("Loading datasets for static dashboard...")
    df_genes = pd.read_csv(PROJECT_ROOT / "results" / "gene_zscores.csv")
    df_intervals = pd.read_csv(PROJECT_ROOT / "results" / "curated_intervals.csv")
    df_baselines = pd.read_csv(PROJECT_ROOT / "results" / "genome_baseline_stats.csv")
    df_stats = pd.read_csv(PROJECT_ROOT / "results" / "Supplementary_Table_S4_Statistical_Contrasts.csv")

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
        if not base_row.empty:
            b = base_row.iloc[0]
            assembly = str(b.get("assembly", "N/A"))
            mean_cds = float(b.get("mean_cds_gc", 0.0))
            std_cds = float(b.get("std_cds_gc", 0.0))
            num_genes = int(b.get("num_genes", 0)) if not pd.isna(b.get("num_genes")) else "N/A"
        else:
            assembly = "N/A"
            mean_cds = 0.0
            std_cds = 0.0
            num_genes = "N/A"

        # Globin locus metrics
        sp_genes = df_genes[df_genes["species"] == sp]
        alpha_genes = sp_genes[sp_genes["comparative_group"].str.contains("alpha_globin", case=False, na=False)]
        beta_genes = sp_genes[sp_genes["comparative_group"].str.contains("beta_globin", case=False, na=False)]

        z_col = "z_score_full_gene_gc" if "z_score_full_gene_gc" in sp_genes.columns else ("zscore_full_gene" if "zscore_full_gene" in sp_genes.columns else None)
        alpha_gc = float(alpha_genes["cds_gc"].mean()) if not alpha_genes.empty else None
        alpha_z = float(alpha_genes[z_col].mean()) if (z_col and not alpha_genes.empty) else None
        beta_gc = float(beta_genes["cds_gc"].mean()) if not beta_genes.empty else None
        beta_z = float(beta_genes[z_col].mean()) if (z_col and not beta_genes.empty) else None

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
            line-height: 1.5;
            padding: 1.5rem;
        }}
        header {{
            background: linear-gradient(135deg, #1e293b, #0f172a);
            color: white;
            padding: 1.8rem 2rem;
            border-radius: 12px;
            margin-bottom: 1.5rem;
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
        }}
        header h1 {{ font-size: 1.8rem; font-weight: 700; margin-bottom: 0.3rem; display: flex; align-items: center; gap: 0.5rem; }}
        header p {{ color: #94a3b8; font-size: 0.95rem; }}
        .badge {{
            display: inline-block;
            background: #334155;
            color: #38bdf8;
            padding: 0.2rem 0.6rem;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 600;
            margin-top: 0.5rem;
        }}

        nav.tabs {{
            display: flex;
            gap: 0.5rem;
            border-bottom: 2px solid var(--border);
            margin-bottom: 1.5rem;
            overflow-x: auto;
        }}
        nav.tabs button {{
            background: none;
            border: none;
            padding: 0.75rem 1.25rem;
            font-size: 0.95rem;
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

        .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1.5rem; }}
        @media (max-width: 900px) {{ .grid-2 {{ grid-template-columns: 1fr; }} }}

        .card {{
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 1.25rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
            margin-bottom: 1.25rem;
        }}
        .card h3 {{ font-size: 1.15rem; font-weight: 600; margin-bottom: 0.8rem; color: #1e293b; }}

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
            gap: 0.4rem;
            flex-wrap: wrap;
            margin-bottom: 1rem;
        }}
        .btn {{
            background: #f1f5f9;
            border: 1px solid var(--border);
            padding: 0.45rem 0.8rem;
            border-radius: 6px;
            font-size: 0.85rem;
            font-weight: 600;
            color: #475569;
            cursor: pointer;
            transition: all 0.15s;
        }}
        .btn:hover {{ background: #e2e8f0; color: #0f172a; }}
        .btn.active {{ background: var(--primary); color: white; border-color: var(--primary); }}

        .img-container {{
            background: white;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 0.5rem;
            text-align: center;
            min-height: 380px;
            display: flex;
            align-items: center;
            justify-content: center;
            overflow: hidden;
        }}
        .img-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
            transition: transform 0.2s;
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
        <h1>🧬 Vertebrate Globin Genomics Explorer</h1>
        <p>Interactive Comparative Atlas, Synteny Tracker & 1-on-1 Head-to-Head Architectural Visualizer</p>
        <span class="badge">44 Vertebrate Genomes • 80 Curated Loci • 1,017 Verified Gene Models • 9 Evolutionary Clades</span>
    </header>

    <nav class="tabs">
        <button class="active" onclick="showTab('comparison')">⚖️ 1-on-1 Head-to-Head Comparison</button>
        <button onclick="showTab('deepdive')">🔬 Single Species Deep Dive</button>
        <button onclick="showTab('meta')">🌐 Cross-Species Meta Explorer</button>
        <button onclick="showTab('genes')">📋 Interval Genes Database</button>
        <button onclick="showTab('stats')">📊 Statistical Contrasts</button>
    </nav>

    <!-- TAB 1: 1-on-1 HEAD-TO-HEAD COMPARISON -->
    <div id="tab-comparison" class="tab-pane active">
        <div class="card">
            <h3>⚖️ Configure 1-on-1 Species Contrast</h3>
            <div class="grid-2">
                <div>
                    <label style="font-size: 0.85rem; font-weight: 600; color: #475569;">SPECIES A</label>
                    <select id="comp-sp-a" onchange="updateComparison()"></select>
                </div>
                <div>
                    <label style="font-size: 0.85rem; font-weight: 600; color: #475569;">SPECIES B</label>
                    <select id="comp-sp-b" onchange="updateComparison()"></select>
                </div>
            </div>
        </div>

        <div class="card">
            <h3>📊 Comparative Scorecard & Isochore Polarity</h3>
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
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.8rem;">
                <h3>🖼️ Synchronized Visual Comparison</h3>
                <div class="btn-group" id="comp-view-btns">
                    <button class="btn active" onclick="setCompView('chromosome')">📍 Linear Chromosome Tracks</button>
                    <button class="btn" onclick="setCompView('inset')">🎯 Bivariate with Mini Inset</button>
                    <button class="btn" onclick="setCompView('diagnostic')">📈 Diagnostic GC Distribution</button>
                    <button class="btn" onclick="setCompView('bivariate')">📐 Standard Bivariate</button>
                </div>
            </div>

            <div class="grid-2">
                <div>
                    <h4 id="lbl-img-a" style="margin-bottom: 0.5rem; color: #1e293b;">Species A</h4>
                    <div class="img-container">
                        <img id="img-comp-a" src="" alt="Plot A">
                    </div>
                </div>
                <div>
                    <h4 id="lbl-img-b" style="margin-bottom: 0.5rem; color: #1e293b;">Species B</h4>
                    <div class="img-container">
                        <img id="img-comp-b" src="" alt="Plot B">
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 2: SINGLE SPECIES DEEP DIVE -->
    <div id="tab-deepdive" class="tab-pane">
        <div class="grid-2" style="grid-template-columns: 320px 1fr;">
            <div>
                <div class="card">
                    <h3>Select Species</h3>
                    <select id="deep-sp" onchange="updateDeepDive()"></select>
                    
                    <h4 style="margin-top: 1rem; margin-bottom: 0.5rem; font-size: 0.95rem; color: #334155;">Species Profile</h4>
                    <div id="deep-profile" style="font-size: 0.88rem; line-height: 1.7; color: #475569;"></div>
                </div>
            </div>
            <div>
                <div class="card">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.8rem;">
                        <h3 id="deep-title">Visualizations</h3>
                        <div class="btn-group" id="deep-view-btns">
                            <button class="btn active" onclick="setDeepView('chromosome')">📍 Linear Chromosome</button>
                            <button class="btn" onclick="setDeepView('inset')">🎯 Bivariate with Inset</button>
                            <button class="btn" onclick="setDeepView('diagnostic')">📈 Diagnostic GC</button>
                            <button class="btn" onclick="setDeepView('bivariate')">📐 Standard Bivariate</button>
                        </div>
                    </div>
                    <div class="img-container">
                        <img id="deep-img" src="" alt="Deep dive plot">
                    </div>
                </div>
            </div>
        </div>
    </div>

    <!-- TAB 3: CROSS-SPECIES META EXPLORER -->
    <div id="tab-meta" class="tab-pane">
        <div class="card">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.8rem;">
                <h3>Global Macroevolutionary Figures</h3>
                <div class="btn-group" id="meta-view-btns">
                    <button class="btn active" onclick="setMetaView('clean')">🌌 Standardized Meta Graph</button>
                    <button class="btn" onclick="setMetaView('fig1')">📊 Fig 1: Locus Divergence</button>
                    <button class="btn" onclick="setMetaView('fig2')">🔬 Fig 2: Mechanistic Drivers</button>
                    <button class="btn" onclick="setMetaView('fig3')">🗺️ Fig 3: Clade Heatmap</button>
                </div>
            </div>
            <div class="img-container">
                <img id="meta-img" src="meta_comparative_graph_clean.png" alt="Meta figure">
            </div>
        </div>
    </div>

    <!-- TAB 4: INTERVAL GENES DATABASE -->
    <div id="tab-genes" class="tab-pane">
        <div class="card">
            <h3>Annotated Genes Search Database (1,017 Genes)</h3>
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

    <script>
        const speciesData = {json.dumps(species_data)};
        const genesData = {json.dumps(genes_data)};
        let currentCompView = 'chromosome';
        let currentDeepView = 'chromosome';

        function showTab(tabId) {{
            document.querySelectorAll('.tab-pane').forEach(el => el.classList.remove('active'));
            document.querySelectorAll('nav.tabs button').forEach(el => el.classList.remove('active'));
            document.getElementById('tab-' + tabId).classList.add('active');
            event.target.classList.add('active');
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

            selA.value = 'Homo_sapiens';
            selB.value = 'Gallus_gallus';
            selDeep.value = 'Homo_sapiens';

            updateComparison();
            updateDeepDive();
            renderGenesTable(genesData.slice(0, 100));
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
            event.target.classList.add('active');
            updateComparisonImages();
        }}

        function updateComparisonImages() {{
            const spA = document.getElementById('comp-sp-a').value;
            const spB = document.getElementById('comp-sp-b').value;
            document.getElementById('lbl-img-a').innerText = `${{spA}} (${{speciesData[spA]?.common_name || ''}})`;
            document.getElementById('lbl-img-b').innerText = `${{spB}} (${{speciesData[spB]?.common_name || ''}})`;
            document.getElementById('img-comp-a').src = getPlotPath(spA, currentCompView);
            document.getElementById('img-comp-b').src = getPlotPath(spB, currentCompView);
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
                ['α-Globin GC / Whole-Gene Z', dA.alpha_gc ? `${{dA.alpha_gc.toFixed(1)}}% (Z=${{dA.alpha_z > 0 ? '+' : ''}}${{dA.alpha_z.toFixed(2)}}σ)` : 'N/A', dB.alpha_gc ? `${{dB.alpha_gc.toFixed(1)}}% (Z=${{dB.alpha_z > 0 ? '+' : ''}}${{dB.alpha_z.toFixed(2)}}σ)` : 'N/A', (dA.alpha_z && dB.alpha_z) ? `ΔZ = ${{ (dA.alpha_z - dB.alpha_z).toFixed(2) }}σ` : 'N/A'],
                ['β-Globin GC / Whole-Gene Z', dA.beta_gc ? `${{dA.beta_gc.toFixed(1)}}% (Z=${{dA.beta_z > 0 ? '+' : ''}}${{dA.beta_z.toFixed(2)}}σ)` : 'N/A', dB.beta_gc ? `${{dB.beta_gc.toFixed(1)}}% (Z=${{dB.beta_z > 0 ? '+' : ''}}${{dB.beta_z.toFixed(2)}}σ)` : 'N/A', (dA.beta_z && dB.beta_z) ? `ΔZ = ${{ (dA.beta_z - dB.beta_z).toFixed(2) }}σ` : 'N/A'],
                ['Subunit Polarity (ΔZ = α - β)', dA.polarity ? `${{dA.polarity > 0 ? '+' : ''}}${{dA.polarity.toFixed(2)}}σ` : 'N/A', dB.polarity ? `${{dB.polarity > 0 ? '+' : ''}}${{dB.polarity.toFixed(2)}}σ` : 'N/A', (dA.polarity && dB.polarity) ? ((dA.polarity > 0 && dB.polarity < 0) || (dA.polarity < 0 && dB.polarity > 0) ? '⚡ Polarity Inversion' : 'Conserved Polarity') : 'N/A']
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
            event.target.classList.add('active');
            const sp = document.getElementById('deep-sp').value;
            document.getElementById('deep-img').src = getPlotPath(sp, currentDeepView);
        }}

        function updateDeepDive() {{
            const sp = document.getElementById('deep-sp').value;
            const d = speciesData[sp];
            document.getElementById('deep-title').innerText = `${{sp}} (${{d.common_name}}) Architecture & Distribution`;
            document.getElementById('deep-img').src = getPlotPath(sp, currentDeepView);

            document.getElementById('deep-profile').innerHTML = `
                <div><b>Common Name:</b> ${{d.common_name}}</div>
                <div><b>Clade:</b> <code>${{d.clade}}</code></div>
                <div><b>Assembly:</b> <code>${{d.assembly}}</code></div>
                <div><b>Baseline CDS GC:</b> ${{d.mean_cds_gc.toFixed(2)}}% (±${{d.std_cds_gc.toFixed(2)}}%)</div>
                <div><b>Total Coding Genes:</b> ${{typeof d.num_genes === 'number' ? d.num_genes.toLocaleString() : d.num_genes}}</div>
                <div><b>Curated Loci:</b> ${{d.loci_count}} intervals</div>
            `;
        }}

        function setMetaView(figType) {{
            document.querySelectorAll('#meta-view-btns button').forEach(b => b.classList.remove('active'));
            event.target.classList.add('active');
            const img = document.getElementById('meta-img');
            if (figType === 'clean') img.src = 'meta_comparative_graph_clean.png';
            if (figType === 'fig1') img.src = 'summary_figure1_locus_divergence.png';
            if (figType === 'fig2') img.src = 'summary_figure2_mechanistic_drivers.png';
            if (figType === 'fig3') img.src = 'summary_figure3_clade_heatmap.png';
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

        window.onload = populateSpeciesDropdowns;
    </script>
</body>
</html>
"""

    out_file = PROJECT_ROOT / "results" / "dashboard.html"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"✓ Standalone HTML dashboard successfully created at: {out_file}")

if __name__ == "__main__":
    main()
