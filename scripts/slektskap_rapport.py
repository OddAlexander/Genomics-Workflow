#!/usr/bin/env python3
"""
slektskap_rapport.py -- Generer utbruddsrapport frå slektskapsanalyse

Lese inn:
  - SNP-matrise (snp-dists)
  - Parvise SNP-avstandar med vurdering
  - IQ-TREE newick-fil
  - snippy-core statistikk
  - Lokal HierCC cluster types (chewBBACA)
  - Enterobase cluster types (valfritt)
  - Per-prøve: MLST, AMRFinder, GAMBIT (frå hovudpipeline)

Produserer:
  - utbrudd_rapport.txt  (maskinleseleg oppsummering)
  - utbrudd_rapport.json (fullt datasett)
  - utbrudd_rapport.html (interaktiv HTML-rapport)
"""

import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path


# =============================================================================
# Argumentparsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Generer utbruddsrapport")
    p.add_argument("--outbreak-name",  required=True)
    p.add_argument("--species",        default="")
    p.add_argument("--samples",        required=True, help="Kommaseparert liste")
    p.add_argument("--results-dir",    required=True, help="Hovudpipeline results/")
    p.add_argument("--outbreak-dir",   required=True)
    p.add_argument("--ref",            default="")
    p.add_argument("--threshold",      type=int, default=20)
    p.add_argument("--matrix",         required=True)
    p.add_argument("--pairs",          required=True)
    p.add_argument("--tree",           required=True)
    p.add_argument("--core-stats",     required=True)
    p.add_argument("--hc",             default="")
    p.add_argument("--ent-ct",         default="")
    p.add_argument("--mst-meta",       default="")
    p.add_argument("--use-gubbins",    default="false")
    p.add_argument("--output-txt",     required=True)
    p.add_argument("--output-json",    required=True)
    p.add_argument("--output-html",    required=True)
    return p.parse_args()


# =============================================================================
# Hjelpefunksjonar -- lesing av pipeline-output
# =============================================================================

def safe_read(path):
    """Les fil eller returner None."""
    if not path or not os.path.exists(path):
        return None
    try:
        return open(path, encoding="utf-8", errors="replace").read()
    except Exception:
        return None


def parse_tsv(path):
    """Les TSV til liste av dicts."""
    txt = safe_read(path)
    if not txt:
        return []
    return list(csv.DictReader(txt.splitlines(), delimiter="\t"))


def parse_snp_matrix(path):
    """Les snp-dists matrise -> {sample: {sample: dist}}."""
    txt = safe_read(path)
    if not txt:
        return {}, []
    lines = [l for l in txt.splitlines() if l.strip()]
    reader = csv.reader(lines, delimiter="\t")
    header = next(reader)[1:]
    matrix = {}
    for row in reader:
        s1 = row[0]
        matrix[s1] = {}
        for i, val in enumerate(row[1:]):
            try:
                matrix[s1][header[i]] = int(val)
            except ValueError:
                matrix[s1][header[i]] = None
    return matrix, header


def parse_core_stats(path):
    """Les snippy-core stats (core.txt)."""
    txt = safe_read(path)
    if not txt:
        return {}
    stats = {}
    for line in txt.splitlines():
        if "\t" in line:
            k, _, v = line.partition("\t")
            stats[k.strip()] = v.strip()
    return stats


def parse_mlst(path):
    """Les mlst.tsv -> {'ST': '..', 'scheme': '..', 'alleles': {}}."""
    txt = safe_read(path)
    if not txt:
        return {}
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 2:
        return {}
    cols = lines[1].split("\t")
    # Format: FILE  SCHEME  ST  allele1  allele2 ...
    if len(cols) < 3:
        return {}
    scheme  = cols[1] if len(cols) > 1 else "-"
    st      = cols[2] if len(cols) > 2 else "-"
    headers = lines[0].split("\t")
    alleles = {}
    for i, h in enumerate(headers[3:], 3):
        alleles[h] = cols[i] if i < len(cols) else "-"
    return {"ST": st, "scheme": scheme, "alleles": alleles}


def parse_amrfinder(path):
    """Les AMRFinder TSV -> liste av {gene, type, subtype, subclass, contig, coverage, identity}."""
    rows = parse_tsv(path)
    genes = []
    for r in rows:
        gene_type = r.get("Type", r.get("Element type", "")).upper()
        if gene_type not in ("AMR", "VIRULENCE", "STRESS"):
            continue
        genes.append({
            "gene":     r.get("Element symbol", r.get("Gene symbol", "-")),
            "name":     r.get("Element name",   r.get("Sequence name", "-")),
            "type":     gene_type,
            "subtype":  r.get("Subtype", "-"),
            "subclass": r.get("Subclass", "-"),
            "scope":    r.get("Scope", "-"),
            "contig":   r.get("Contig id", "-"),
            "coverage": r.get("% Coverage of reference sequence", "-"),
            "identity": r.get("% Identity to reference sequence", "-"),
        })
    return genes


def parse_gambit(path):
    """Les GAMBIT CSV -> {'name': str, 'rank': str, 'distance': float}."""
    txt = safe_read(path)
    if not txt:
        return {}
    rows = list(csv.DictReader(txt.splitlines()))
    if not rows:
        return {}
    r = rows[0]
    def _f(k):
        try:
            return float(r.get(k, "") or "")
        except ValueError:
            return None
    return {
        "name":     r.get("predicted.name", "-") or "-",
        "rank":     r.get("predicted.rank", "-") or "-",
        "distance": _f("closest.distance"),
    }


def read_newick(path):
    txt = safe_read(path)
    return txt.strip() if txt else ""


def load_per_sample_data(samples, results_dir):
    """
    Hent per-prøve data frå hovudpipelinen.
    Returner dict: {sample -> {mlst, amr, gambit}}
    """
    data = {}
    for s in samples:
        base = Path(results_dir) / s

        # MLST -- prøv begge vanlege stiar
        mlst = parse_mlst(base / "MLST" / "mlst.tsv")

        # AMRFinder -- hovudpipeline skriv til AMRFinder/amrfinder.tsv
        amr = parse_amrfinder(base / "AMRFinder" / "amrfinder.tsv")
        if not amr:
            # Eldre pipeline skreiv til AMR/amrfinder_amr.tsv
            amr = parse_amrfinder(base / "AMR" / "amrfinder_amr.tsv")

        # GAMBIT
        gambit = parse_gambit(base / "GAMBIT" / "gambit.csv")

        data[s] = {"mlst": mlst, "amr": amr, "gambit": gambit}
    return data


# =============================================================================
# Byggelogikk
# =============================================================================

def build_report_data(args, samples):
    """Sett saman all rapport-data til eitt stort dict."""
    matrix, mat_samples = parse_snp_matrix(args.matrix)
    pairs       = parse_tsv(args.pairs)
    core_stats  = parse_core_stats(args.core_stats)
    newick      = read_newick(args.tree)
    per_sample  = load_per_sample_data(samples, args.results_dir)

    # HC cluster types -- lokal (prioritert) eller Enterobase
    local_hc = {}
    ent_hc   = {}
    if args.hc:
        for row in parse_tsv(args.hc):
            local_hc[row.get("sample", "")] = row
    if args.ent_ct:
        for row in parse_tsv(args.ent_ct):
            ent_hc[row.get("sample", "")] = row

    use_gubbins = args.use_gubbins.lower() in ("true", "1", "yes")

    # Utbruddssamandrag
    cluster_pairs = [p for p in pairs if p.get("vurdering") == "SANNSYNLIG UTBRUDD"]
    possible_pairs = [p for p in pairs if p.get("vurdering") == "MULIG TILKNYTNING"]

    return {
        "meta": {
            "outbreak_name": args.outbreak_name,
            "species":       args.species,
            "samples":       samples,
            "ref":           args.ref,
            "threshold":     args.threshold,
            "use_gubbins":   use_gubbins,
            "generated":     datetime.now().isoformat(timespec="seconds"),
            "n_samples":     len(samples),
        },
        "core_stats":     core_stats,
        "snp_matrix":     matrix,
        "snp_pairs":      pairs,
        "cluster_pairs":  cluster_pairs,
        "possible_pairs": possible_pairs,
        "local_hc":       local_hc,
        "enterobase_hc":  ent_hc,
        "per_sample":     per_sample,
        "newick":         newick,
    }


# =============================================================================
# TXT-rapport
# =============================================================================

def write_txt(data, path):
    lines = []
    m = data["meta"]
    lines += [
        "=" * 72,
        f"UTBRUDDSANALYSE: {m['outbreak_name']}",
        "=" * 72,
        f"Generert   : {m['generated']}",
        f"Art        : {m['species'] or 'ikkje angitt'}",
        f"Prøver     : {m['n_samples']}",
        f"Referanse  : {m['ref']}",
        f"SNP-terskel: {m['threshold']}",
        f"Gubbins    : {'ja' if m['use_gubbins'] else 'nei'}",
        "",
    ]

    # Core-statistikk
    cs = data["core_stats"]
    if cs:
        lines += ["--- Kjerne-alignment ---"]
        for k, v in cs.items():
            lines.append(f"  {k:<30} {v}")
        lines.append("")

    # Utbruddssamandrag
    cp = data["cluster_pairs"]
    pp = data["possible_pairs"]
    lines += [
        "--- Utbruddsvurdering ---",
        f"  Sannsynleg utbrudd  ({m['threshold']} SNP) : {len(cp)} par",
        f"  Muleg tilknytning  ({m['threshold']*2} SNP) : {len(pp)} par",
        "",
    ]
    if cp:
        lines.append("  Sannsynlege par:")
        for p in cp:
            lines.append(f"    {p['prøve_1']} -- {p['prøve_2']}: {p['snp_avstand']} SNP")
        lines.append("")

    # SNP-matrise
    matrix = data["snp_matrix"]
    samples = m["samples"]
    if matrix:
        lines.append("--- SNP-matrise ---")
        col_w = max(len(s) for s in samples) + 2
        hdr = " " * col_w + "".join(f"{s:>{col_w}}" for s in samples)
        lines.append(hdr)
        for s1 in samples:
            row_vals = ""
            for s2 in samples:
                v = matrix.get(s1, {}).get(s2)
                row_vals += f"{str(v) if v is not None else '-':>{col_w}}"
            lines.append(f"{s1:<{col_w}}{row_vals}")
        lines.append("")

    # Per-prøve oppsummering
    lines.append("--- Per prøve ---")
    for s in samples:
        ps  = data["per_sample"].get(s, {})
        mlst = ps.get("mlst", {})
        amr  = ps.get("amr", [])
        hc_row = (data["enterobase_hc"] or data["local_hc"]).get(s, {})
        hc_src = "Enterobase" if data["enterobase_hc"] else "lokal"
        amr_genes = [g["gene"] for g in amr if g["type"] == "AMR"]
        lines += [
            f"  {s}",
            f"    ST      : {mlst.get('ST', '-')} ({mlst.get('scheme', '-')})",
            f"    HC10    : {hc_row.get('HC10', '-')} [{hc_src}]",
            f"    HC50    : {hc_row.get('HC50', '-')} [{hc_src}]",
            f"    AMR-gen : {', '.join(amr_genes) if amr_genes else 'ingen'}",
            "",
        ]

    # Newick
    lines += ["--- Fylogenetisk tre (Newick) ---", data["newick"], ""]

    os.makedirs(Path(path).parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# =============================================================================
# JSON-rapport
# =============================================================================

def write_json(data, path):
    os.makedirs(Path(path).parent, exist_ok=True)

    # Gjer matrix JSON-serierbar (int / None)
    matrix_out = {}
    for s1, row in data["snp_matrix"].items():
        matrix_out[s1] = {s2: v for s2, v in row.items()}

    out = {
        "meta":          data["meta"],
        "core_stats":    data["core_stats"],
        "snp_matrix":    matrix_out,
        "snp_pairs":     data["snp_pairs"],
        "cluster_pairs": data["cluster_pairs"],
        "possible_pairs":data["possible_pairs"],
        "local_hc":      data["local_hc"],
        "enterobase_hc": data["enterobase_hc"],
        "per_sample": {
            s: {
                "mlst":   v["mlst"],
                "amr":    v["amr"],
                "gambit": v["gambit"],
            }
            for s, v in data["per_sample"].items()
        },
        "newick": data["newick"],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)


# =============================================================================
# HTML-rapport
# =============================================================================

_CSS = """
:root {
  --green: #27ae60; --red: #e74c3c; --orange: #e67e22;
  --blue: #2980b9;  --gray: #95a5a6;
  --bg: #f4f6f9;    --card: #ffffff;
  --border: #dee2e6;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: 'Segoe UI', Arial, sans-serif;
  background: var(--bg); color: #2c3e50; font-size: 14px;
}
header {
  background: #1a252f; color: #fff;
  padding: 20px 32px; display: flex;
  align-items: center; justify-content: space-between;
}
header h1 { font-size: 1.4rem; font-weight: 600; }
header .meta { font-size: 0.8rem; opacity: 0.7; }
.container { max-width: 1400px; margin: 0 auto; padding: 24px 32px; }
.section { margin-bottom: 32px; }
.section h2 {
  font-size: 1rem; font-weight: 600; text-transform: uppercase;
  letter-spacing: .05em; color: #1a252f;
  border-bottom: 2px solid var(--border);
  padding-bottom: 6px; margin-bottom: 14px;
}
.card {
  background: var(--card); border-radius: 8px;
  box-shadow: 0 1px 4px rgba(0,0,0,.07);
  padding: 18px 22px;
}
.card-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 14px; margin-bottom: 24px;
}
.stat-card {
  background: var(--card); border-radius: 8px;
  box-shadow: 0 1px 4px rgba(0,0,0,.07);
  padding: 16px 18px; text-align: center;
}
.stat-card .val {
  font-size: 2rem; font-weight: 700; line-height: 1;
}
.stat-card .lbl { font-size: 0.75rem; color: #666; margin-top: 4px; }
.val.green { color: var(--green); }
.val.red   { color: var(--red);   }
.val.orange{ color: var(--orange);}
table {
  width: 100%; border-collapse: collapse; font-size: 0.84rem;
}
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--border); }
th { background: #f0f3f7; font-weight: 600; position: sticky; top: 0; }
tr:hover td { background: #f8f9fb; }
.badge {
  display: inline-block; padding: 2px 7px; border-radius: 3px;
  font-size: 0.73rem; font-weight: 600; text-transform: uppercase;
}
.badge-green  { background: #d5f5e3; color: #1e8449; }
.badge-red    { background: #fadbd8; color: #922b21; }
.badge-orange { background: #fdebd0; color: #935116; }
.badge-gray   { background: #eaecee; color: #555; }
.verdict-outbreak { color: var(--red);    font-weight: 700; }
.verdict-possible { color: var(--orange); font-weight: 600; }
.verdict-unrelated{ color: var(--gray); }
.matrix-table td, .matrix-table th {
  text-align: right; min-width: 52px; font-size: 0.78rem; font-family: monospace;
}
.matrix-table td.diag { background: #edf2f7; color: #aaa; }
.matrix-table td.cluster { background: #fde8e8; font-weight: 700; color: var(--red); }
.matrix-table td.possible { background: #fef3e2; color: var(--orange); }
.tree-box {
  font-family: monospace; font-size: 0.76rem; white-space: pre-wrap;
  word-break: break-all; max-height: 160px; overflow-y: auto;
  background: #f8f9fb; border: 1px solid var(--border);
  border-radius: 4px; padding: 10px 14px;
}
.hc-table td { font-family: monospace; }
.amr-chip {
  display: inline-block; background: #fde8e8; color: #922b21;
  border-radius: 3px; padding: 1px 6px; margin: 1px 2px;
  font-size: 0.75rem;
}
.vir-chip {
  display: inline-block; background: #fef3e2; color: #935116;
  border-radius: 3px; padding: 1px 6px; margin: 1px 2px;
  font-size: 0.75rem;
}
code { font-family: monospace; font-size: 0.85em; }
"""

_VERDICT_CLASS = {
    "SANNSYNLIG UTBRUDD": "verdict-outbreak",
    "MULIG TILKNYTNING":  "verdict-possible",
    "IKKE RELATERT":      "verdict-unrelated",
}

_VERDICT_BADGE = {
    "SANNSYNLIG UTBRUDD": "badge-red",
    "MULIG TILKNYTNING":  "badge-orange",
    "IKKE RELATERT":      "badge-gray",
}


def _h(text):
    """HTML-escape."""
    return (str(text)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;"))


def build_html(data):
    m         = data["meta"]
    samples   = m["samples"]
    matrix    = data["snp_matrix"]
    pairs     = data["snp_pairs"]
    cluster   = data["cluster_pairs"]
    possible  = data["possible_pairs"]
    local_hc  = data["local_hc"]
    ent_hc    = data["enterobase_hc"]
    ps        = data["per_sample"]
    cs        = data["core_stats"]
    newick    = data["newick"]
    threshold = m["threshold"]

    hc_data   = ent_hc if ent_hc else local_hc
    hc_src    = "Enterobase" if ent_hc else "Lokal HierCC"
    hc_levels = [2, 5, 10, 20, 50, 100, 200]

    # ---- Samandragskort
    n_cluster  = len(cluster)
    n_possible = len(possible)
    if n_cluster > 0:
        outbreak_badge = f'<span class="badge badge-red">AKTIVT UTBRUDD</span>'
    elif n_possible > 0:
        outbreak_badge = f'<span class="badge badge-orange">MULEG TILKNYTNING</span>'
    else:
        outbreak_badge = f'<span class="badge badge-green">INGEN KLYNGE</span>'

    stat_cards = "".join([
        f'<div class="stat-card"><div class="val {"red" if n_cluster > 0 else "green"}">{n_cluster}</div>'
        f'<div class="lbl">Par ≤ {threshold} SNP<br>(sannsynleg utbrudd)</div></div>',
        f'<div class="stat-card"><div class="val {"orange" if n_possible > 0 else "green"}">{n_possible}</div>'
        f'<div class="lbl">Par ≤ {threshold*2} SNP<br>(muleg tilknyting)</div></div>',
        f'<div class="stat-card"><div class="val">{len(samples)}</div>'
        f'<div class="lbl">Prøver</div></div>',
        f'<div class="stat-card"><div class="val">{threshold}</div>'
        f'<div class="lbl">SNP-terskel</div></div>',
    ])

    # ---- Core statistikk
    core_rows = ""
    for k, v in cs.items():
        core_rows += f"<tr><td>{_h(k)}</td><td><code>{_h(v)}</code></td></tr>"
    core_section = ""
    if core_rows:
        core_section = (
            '<div class="section">'
            '<h2>Kjerne-alignment</h2>'
            '<div class="card"><table>'
            f'<thead><tr><th>Parameter</th><th>Verdi</th></tr></thead>'
            f'<tbody>{core_rows}</tbody>'
            '</table></div></div>'
        )

    # ---- SNP-matrise
    def cell_class(s1, s2, d):
        if s1 == s2:
            return "diag"
        if d is not None:
            if d <= threshold:
                return "cluster"
            if d <= threshold * 2:
                return "possible"
        return ""

    mat_header = "<th></th>" + "".join(f"<th>{_h(s)}</th>" for s in samples)
    mat_rows = ""
    for s1 in samples:
        mat_rows += f"<tr><th>{_h(s1)}</th>"
        for s2 in samples:
            d = matrix.get(s1, {}).get(s2)
            cls = cell_class(s1, s2, d)
            val = str(d) if d is not None else "-"
            mat_rows += f'<td class="{cls}">{_h(val)}</td>'
        mat_rows += "</tr>"

    matrix_section = (
        '<div class="section">'
        '<h2>SNP-matrise</h2>'
        '<div class="card" style="overflow-x:auto">'
        '<table class="matrix-table">'
        f'<thead><tr>{mat_header}</tr></thead>'
        f'<tbody>{mat_rows}</tbody>'
        '</table></div></div>'
    )

    # ---- Parvise avstandar
    pair_rows = ""
    for p in sorted(pairs, key=lambda x: int(x.get("snp_avstand", 9999) or 9999)):
        v   = p.get("vurdering", "")
        cls = _VERDICT_CLASS.get(v, "")
        bls = _VERDICT_BADGE.get(v, "badge-gray")
        pair_rows += (
            f'<tr>'
            f'<td><code>{_h(p.get("prøve_1",""))}</code></td>'
            f'<td><code>{_h(p.get("prøve_2",""))}</code></td>'
            f'<td style="text-align:right">{_h(p.get("snp_avstand",""))}</td>'
            f'<td>{_h(p.get("terskel",""))}</td>'
            f'<td><span class="badge {bls}">{_h(v)}</span></td>'
            f'</tr>'
        )
    pairs_section = (
        '<div class="section">'
        '<h2>Parvise SNP-avstandar</h2>'
        '<div class="card" style="overflow-x:auto">'
        '<table>'
        '<thead><tr><th>Prøve 1</th><th>Prøve 2</th>'
        '<th style="text-align:right">SNP</th><th>Terskel</th>'
        '<th>Vurdering</th></tr></thead>'
        f'<tbody>{pair_rows}</tbody>'
        '</table></div></div>'
    )

    # ---- HC cluster types
    hc_section = ""
    if hc_data:
        hc_hdr = (
            "<tr><th>Prøve</th>"
            + "".join(f"<th>HC{t}</th>" for t in hc_levels)
            + f'<th>Kjelde</th></tr>'
        )
        hc_rows_html = ""
        for s in samples:
            row = hc_data.get(s, {})
            hc_rows_html += f'<tr><td><code>{_h(s)}</code></td>'
            for t in hc_levels:
                hc_rows_html += f'<td>{_h(row.get(f"HC{t}", "-"))}</td>'
            hc_rows_html += f'<td>{_h(row.get("source", hc_src))}</td></tr>'
        hc_section = (
            '<div class="section">'
            f'<h2>cgMLST HierCC Cluster Types ({_h(hc_src)})</h2>'
            '<div class="card" style="overflow-x:auto">'
            '<table class="hc-table">'
            f'<thead>{hc_hdr}</thead>'
            f'<tbody>{hc_rows_html}</tbody>'
            '</table></div></div>'
        )

    # ---- Per-prøve oversikt
    def amr_chips(genes):
        out = ""
        for g in genes:
            if g["type"] == "AMR":
                out += f'<span class="amr-chip">{_h(g["gene"])}</span>'
            elif g["type"] == "VIRULENCE":
                out += f'<span class="vir-chip">{_h(g["gene"])}</span>'
        return out or "<span style='color:#aaa'>ingen</span>"

    sample_hdr = (
        "<tr><th>Prøve</th><th>Art (GAMBIT)</th><th>ST</th>"
        "<th>HC10</th><th>HC50</th>"
        "<th>AMR / Virulens</th></tr>"
    )
    sample_rows = ""
    for s in samples:
        d       = ps.get(s, {})
        mlst    = d.get("mlst", {})
        amr     = d.get("amr",  [])
        gambit  = d.get("gambit", {})
        hc_row  = hc_data.get(s, {})
        sample_rows += (
            f'<tr>'
            f'<td><code>{_h(s)}</code></td>'
            f'<td>{_h(gambit.get("name", "-"))}</td>'
            f'<td>{_h(mlst.get("ST", "-"))}</td>'
            f'<td>{_h(hc_row.get("HC10", "-"))}</td>'
            f'<td>{_h(hc_row.get("HC50", "-"))}</td>'
            f'<td>{amr_chips(amr)}</td>'
            f'</tr>'
        )

    sample_section = (
        '<div class="section">'
        '<h2>Prøveoversikt</h2>'
        '<div class="card" style="overflow-x:auto">'
        '<table>'
        f'<thead>{sample_hdr}</thead>'
        f'<tbody>{sample_rows}</tbody>'
        '</table></div></div>'
    )

    # ---- AMR-detaljtabell
    amr_detail_rows = ""
    for s in samples:
        for g in ps.get(s, {}).get("amr", []):
            type_cls = {
                "AMR":       "badge-red",
                "VIRULENCE": "badge-orange",
                "STRESS":    "badge-gray",
            }.get(g["type"], "badge-gray")
            amr_detail_rows += (
                f'<tr>'
                f'<td><code>{_h(s)}</code></td>'
                f'<td><span class="badge {type_cls}">{_h(g["type"])}</span></td>'
                f'<td><code>{_h(g["gene"])}</code></td>'
                f'<td>{_h(g["name"])}</td>'
                f'<td>{_h(g["subclass"])}</td>'
                f'<td style="text-align:right">{_h(g["coverage"])}</td>'
                f'<td style="text-align:right">{_h(g["identity"])}</td>'
                f'<td><code style="font-size:0.75em">{_h(g["contig"])}</code></td>'
                f'</tr>'
            )
    amr_detail_section = ""
    if amr_detail_rows:
        amr_detail_section = (
            '<div class="section">'
            '<h2>AMR / Virulens-detaljar</h2>'
            '<div class="card" style="overflow-x:auto">'
            '<table>'
            '<thead><tr>'
            '<th>Prøve</th><th>Type</th><th>Gen</th><th>Namn</th>'
            '<th>Subklasse</th><th>%Cov</th><th>%Id</th><th>Contig</th>'
            '</tr></thead>'
            f'<tbody>{amr_detail_rows}</tbody>'
            '</table></div></div>'
        )

    # ---- Newick-tre
    tree_section = (
        '<div class="section">'
        '<h2>Fylogenetisk tre (Newick)</h2>'
        '<div class="card">'
        f'<div class="tree-box">{_h(newick)}</div>'
        '</div></div>'
    )

    # ---- Metadata
    gubbins_str = "ja" if m["use_gubbins"] else "nei"
    meta_html = (
        f'<p>Art: <strong>{_h(m["species"] or "ikkje angitt")}</strong> &nbsp;|&nbsp; '
        f'Referanse: <code>{_h(Path(m["ref"]).name if m["ref"] else "-")}</code> &nbsp;|&nbsp; '
        f'Gubbins: {gubbins_str}</p>'
    )

    html = f"""<!DOCTYPE html>
<html lang="no">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Utbruddsrapport — {_h(m['outbreak_name'])}</title>
<style>{_CSS}</style>
</head>
<body>
<header>
  <div>
    <h1>Utbruddsanalyse &mdash; {_h(m['outbreak_name'])}</h1>
    <div class="meta">Generert: {_h(m['generated'])} &nbsp;|&nbsp; {_h(m['n_samples'])} prøver</div>
  </div>
  <div>{outbreak_badge}</div>
</header>
<div class="container">

<div class="section">
  <h2>Samandrag</h2>
  <div class="card-grid">{stat_cards}</div>
  {meta_html}
</div>

{core_section}
{matrix_section}
{pairs_section}
{hc_section}
{sample_section}
{amr_detail_section}
{tree_section}

</div>
</body>
</html>
"""
    return html


def write_html(data, path):
    os.makedirs(Path(path).parent, exist_ok=True)
    html = build_html(data)
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)


# =============================================================================
# Inngangspunkt
# =============================================================================

def main():
    args    = parse_args()
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    print(f"[slektskap_rapport] {args.outbreak_name} — {len(samples)} prøver", flush=True)

    data = build_report_data(args, samples)

    write_txt(data,  args.output_txt)
    write_json(data, args.output_json)
    write_html(data, args.output_html)

    n_cluster  = len(data["cluster_pairs"])
    n_possible = len(data["possible_pairs"])
    print(f"[slektskap_rapport] Ferdig. "
          f"Utbruddspar: {n_cluster}, Muleg tilknyting: {n_possible}", flush=True)
    print(f"  TXT : {args.output_txt}")
    print(f"  JSON: {args.output_json}")
    print(f"  HTML: {args.output_html}")


if __name__ == "__main__":
    main()
