#!/usr/bin/env python3
"""Generate HTML report for the cgMLST pipeline run."""
import argparse, csv, json, math, re, subprocess
from datetime import datetime
from pathlib import Path


def nonempty(path):
    if not path:
        return False
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_read(path):
    return Path(path).read_text().strip() if nonempty(path) else None


def _is_valid_allele(v):
    """chewBBACA writes integers for resolved alleles; tags like LNF/NIPH/INF-N for missing."""
    try:
        return int(v) > 0
    except (TypeError, ValueError):
        return False


def parse_allele_profiles(path):
    """Parse chewBBACA results_alleles.tsv → {n_loci, sample_stats}."""
    if not nonempty(path):
        return {"n_loci": 0, "sample_stats": {}}
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        return {"n_loci": 0, "sample_stats": {}}
    n_loci = len(lines[0].split("\t")) - 1
    stats  = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        parts    = line.split("\t")
        n_called = sum(1 for v in parts[1:] if _is_valid_allele(v))
        stats[parts[0]] = {"n_called": n_called, "n_missing": n_loci - n_called}
    return {"n_loci": n_loci, "sample_stats": stats}


def parse_dists(path):
    """Parse cgmlst-dists TSV → {samples, matrix}."""
    if not nonempty(path):
        return {"samples": [], "matrix": []}
    rows = list(csv.reader(Path(path).read_text().splitlines(), delimiter="\t"))
    if not rows:
        return {"samples": [], "matrix": []}
    matrix = []
    for row in rows[1:]:
        try:
            matrix.append([int(v) for v in row[1:]])
        except (IndexError, ValueError):
            pass
    return {"samples": rows[0][1:], "matrix": matrix}


def parse_partitions(path):
    """Parse cluster TSV → {headers, rows}."""
    if not nonempty(path):
        return {"headers": [], "rows": []}
    reader = list(csv.DictReader(Path(path).read_text().splitlines(), delimiter="\t"))
    if not reader:
        return {"headers": [], "rows": []}
    return {"headers": list(reader[0].keys()), "rows": [dict(r) for r in reader]}


def get_sample_species(results_dir, sample):
    p = Path(results_dir) / sample / "species.txt"
    return p.read_text().strip().replace("_", " ") if p.exists() else "-"


def parse_fastp(json_path):
    """Read totals from a fastp JSON; returns dict with reads/bases/q30/dup_pct or Nones."""
    if not nonempty(json_path):
        return {"reads": None, "bases": None, "q30_pct": None, "dup_pct": None}
    summary = json.loads(Path(json_path).read_text()).get("summary", {})
    before  = summary.get("before_filtering", {})
    after   = summary.get("after_filtering",  {})
    return {
        "reads":   after.get("total_reads"),
        "bases":   after.get("total_bases"),
        "q30_pct": round(after.get("q30_rate",          0) * 100, 1),
        "dup_pct": round(before.get("duplication_rate", 0) * 100, 1),
    }


def parse_quast(quast_dir):
    """Return assembly stats from QUAST's report.tsv: {total_length, contigs, N50}."""
    tsv = Path(quast_dir) / "report.tsv"
    out = {"total_length": None, "contigs": None, "N50": None}
    if not tsv.exists():
        return out
    rows = [l.split("\t") for l in tsv.read_text().splitlines() if "\t" in l]
    data = {r[0].strip(): r[1].strip().replace(",", "") for r in rows if len(r) >= 2}
    try:
        out["total_length"] = int(data.get("Total length") or 0) or None
    except ValueError: pass
    try:
        out["N50"] = int(data.get("N50") or 0) or None
    except ValueError: pass
    out["contigs"] = data.get("# contigs") or None
    return out


def parse_bracken_primary_pct(path):
    """Top-species fraction from Bracken's species output, as a percentage."""
    if not nonempty(path):
        return None
    rows = list(csv.DictReader(Path(path).read_text().splitlines(), delimiter="\t"))
    if not rows:
        return None
    top = max(rows, key=lambda r: float(r.get("fraction_total_reads") or 0))
    try:
        return round(float(top["fraction_total_reads"]) * 100, 1)
    except (ValueError, KeyError):
        return None


def parse_quast_total_length(quast_dir):
    """Return Total length (assembly size in bp) from QUAST's report.tsv, or None."""
    tsv = Path(quast_dir) / "report.tsv"
    if not tsv.exists():
        return None
    for line in tsv.read_text().splitlines():
        if line.startswith("Total length\t"):
            try:
                return int(line.split("\t")[1].replace(",", ""))
            except (IndexError, ValueError):
                return None
    return None


def estimate_depth(bases, genome_size):
    """Approximate mean coverage = trimmed bases / assembled length."""
    if not bases or not genome_size:
        return None
    return round(bases / genome_size, 1)


def compute_mst(n, matrix):
    """Prim's algorithm on a dense distance matrix → list of {from, to, dist} edges."""
    if n < 2:
        return []
    visited    = [False] * n
    visited[0] = True
    edges      = []
    for _ in range(n - 1):
        best_d, best_i, best_j = float("inf"), -1, -1
        for i in range(n):
            if not visited[i]:
                continue
            for j in range(n):
                if visited[j]:
                    continue
                if matrix[i][j] < best_d:
                    best_d, best_i, best_j = matrix[i][j], i, j
        if best_j == -1:
            break
        visited[best_j] = True
        edges.append({"from": best_i, "to": best_j, "dist": best_d})
    return edges


def force_layout(n, edges, width=800, height=500, iters=500):
    """Spring-embedder layout → list of {x, y} positions (one per node)."""
    if n == 0:
        return []
    if n == 1:
        return [{"x": width / 2, "y": height / 2}]
    TAU = 2 * math.pi
    pos = [
        [width  / 2 + (width  / 3) * math.cos(TAU * i / n),
         height / 2 + (height / 3) * math.sin(TAU * i / n)]
        for i in range(n)
    ]
    max_d    = max((e["dist"] for e in edges), default=1) or 1
    ideal    = min(width, height) / max(n ** 0.55, 2)
    k_spring = 0.06
    k_rep    = ideal * ideal * 1.8

    for it in range(iters):
        fx = [0.0] * n
        fy = [0.0] * n
        for e in edges:
            i, j   = e["from"], e["to"]
            dx     = pos[j][0] - pos[i][0]
            dy     = pos[j][1] - pos[i][1]
            d      = math.sqrt(dx * dx + dy * dy) or 0.001
            target = ideal * (1 + e["dist"] / max_d * 1.5)
            f      = k_spring * (d - target) / d
            fx[i] += f * dx;  fy[i] += f * dy
            fx[j] -= f * dx;  fy[j] -= f * dy
        for i in range(n):
            for j in range(i + 1, n):
                dx = pos[j][0] - pos[i][0]
                dy = pos[j][1] - pos[i][1]
                d2 = dx * dx + dy * dy or 0.001
                d  = math.sqrt(d2)
                f  = k_rep / (d2 * d)
                fx[i] -= f * dx;  fy[i] -= f * dy
                fx[j] += f * dx;  fy[j] += f * dy
        cool = max(0.05, 1 - it / iters)
        for i in range(n):
            pos[i][0] = max(60, min(width  - 60, pos[i][0] + fx[i] * cool))
            pos[i][1] = max(40, min(height - 40, pos[i][1] + fy[i] * cool))

    # Rescale to fill the viewport with a uniform margin.
    xs = [p[0] for p in pos]
    ys = [p[1] for p in pos]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    rx     = (x1 - x0) or 1
    ry     = (y1 - y0) or 1
    pad    = 60
    for p in pos:
        p[0] = pad + (p[0] - x0) / rx * (width  - 2 * pad)
        p[1] = pad + (p[1] - y0) / ry * (height - 2 * pad)

    return [{"x": round(p[0], 1), "y": round(p[1], 1)} for p in pos]


def collect_tool_versions():
    def ver(cmd):
        try:
            r     = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            first = ((r.stdout + r.stderr).strip().splitlines() or [""])[0]
            m     = re.search(r"(\d+\.\d+(?:\.\d+)*)", first)
            return m.group(1) if m else (first or "-")
        except Exception:
            return "-"
    return {
        "chewBBACA":    ver(["pixi", "run", "--environment", "cgmlst",    "chewBBACA.py", "--version"]),
        "cgmlst-dists": ver(["pixi", "run", "--environment", "cgmlst",    "cgmlst-dists", "--version"]),
        "GrapeTree":    ver(["pixi", "run", "--environment", "grapetree", "grapetree",    "--version"]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alleles",     required=True)
    ap.add_argument("--dists",       required=True)
    ap.add_argument("--grapetree",   required=True)
    ap.add_argument("--partitions",  required=True)
    ap.add_argument("--template",    required=True)
    ap.add_argument("--samples",     required=True)
    ap.add_argument("--schema",      required=True)
    ap.add_argument("--thresholds",  default="7", help="Comma-separated allele thresholds (e.g. '7' or '5,10,24')")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--output",      required=True)
    args = ap.parse_args()

    samples    = args.samples.split(",")
    thresholds = [int(t.strip()) for t in args.thresholds.split(",") if t.strip()]
    profiles   = parse_allele_profiles(args.alleles)
    dists      = parse_dists(args.dists)
    parts      = parse_partitions(args.partitions)
    n_loci     = profiles["n_loci"]

    # chewBBACA, cgmlst-dists and the cluster TSV all key samples by the flat
    # underscore form (19-03-2026_005a). Map back to the leaf name ("005a") so
    # every table in the report shows the same short ID.
    to_leaf = {s.replace("/", "_"): s.split("/")[-1] for s in samples}

    # ReporTree emits one cluster column per threshold (e.g. "single-7x1.0",
    # "single-24x1.0"). Map each ReporTree column back to its integer threshold
    # so the report can label columns as cluster_<N>.
    threshold_columns = []        # list of (threshold:int, reportree_col:str)
    if len(parts["headers"]) > 1:
        id_col = parts["headers"][0]
        for col in parts["headers"][1:]:
            m = re.search(r"(\d+)", col)
            if m:
                threshold_columns.append((int(m.group(1)), col))

    # cluster_maps[threshold][underscore_id] -> cluster label
    cluster_maps = {t: {} for t, _ in threshold_columns}
    if threshold_columns:
        id_col = parts["headers"][0]
        for row in parts["rows"]:
            for t, col in threshold_columns:
                cluster_maps[t][row[id_col]] = row.get(col, "-")

    # Normalise distance-matrix and partition-table IDs to leaf form.
    dists["samples"] = [to_leaf.get(s, s) for s in dists["samples"]]
    if parts["headers"]:
        id_col = parts["headers"][0]
        for row in parts["rows"]:
            row[id_col] = to_leaf.get(row[id_col], row[id_col])

    # Rename ReporTree's per-threshold columns to a cleaner cluster_<N> form
    # so the report rendering and template stay readable.
    rename_map = {col: f"cluster_{t}" for t, col in threshold_columns}
    if rename_map:
        parts["headers"] = [rename_map.get(h, h) for h in parts["headers"]]
        for row in parts["rows"]:
            for old, new in rename_map.items():
                if old in row:
                    row[new] = row.pop(old)

    results_dir = Path(args.results_dir)
    sample_data = []
    for s in samples:
        underscore_id = s.replace("/", "_")
        st            = profiles["sample_stats"].get(underscore_id, {})
        n_called      = st.get("n_called")
        n_missing     = st.get("n_missing")
        pct           = round(n_called / n_loci * 100, 1) if n_loci and n_called is not None else None

        fp          = parse_fastp(results_dir / s / "QC/fastp.json")
        quast_stats = parse_quast(results_dir / s / "QUAST")
        genome_size = quast_stats["total_length"]
        est_depth   = estimate_depth(fp["bases"], genome_size)
        kraken_pct  = parse_bracken_primary_pct(results_dir / s / "ID_Kraken2/bracken_species.txt")

        # Cluster per threshold + primary cluster (the smallest/tightest threshold).
        clusters = {t: cluster_maps[t].get(underscore_id, "-") for t in cluster_maps}
        primary  = clusters[min(clusters)] if clusters else "-"

        sample_data.append({
            "id":                 s.split("/")[-1],
            "species":            get_sample_species(args.results_dir, s),
            "n_loci":             n_loci,
            "n_called":           n_called,
            "n_missing":          n_missing,
            "pct_called":         pct,
            "clusters":           clusters,
            "cluster":            primary,
            "reads":              fp["reads"],
            "q30_pct":            fp["q30_pct"],
            "est_depth":          est_depth,
            "genome_size":        genome_size,
            "contigs":            quast_stats["contigs"],
            "N50":                quast_stats["N50"],
            "kraken_primary_pct": kraken_pct,
        })

    # MST positions, pre-computed so the HTML renderer just draws them.
    mst_edges = compute_mst(len(dists["samples"]), dists["matrix"])
    mst_nodes = force_layout(len(dists["samples"]), mst_edges)

    data = {
        "run_date":      datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_samples":     len(samples),
        "schema":        args.schema,
        "thresholds":    thresholds,
        "threshold":     thresholds[0] if thresholds else None,   # back-compat: primary
        "n_loci":        n_loci,
        "samples":       sample_data,
        "dists":         dists,
        "partitions":    parts,
        "grapetree":     safe_read(args.grapetree) or "",
        "mst_nodes":     mst_nodes,
        "mst_edges":     mst_edges,
        "tool_versions": collect_tool_versions(),
    }

    Path(args.output).write_text(
        Path(args.template).read_text().replace(
            "__CGMLST_JSON_DATA__", json.dumps(data, ensure_ascii=False)
        )
    )


if __name__ == "__main__":
    main()
