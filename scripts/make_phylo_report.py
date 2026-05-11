#!/usr/bin/env python3
"""Generate HTML report for phylogenomics pipeline run."""
import argparse, json, csv, subprocess, re
from datetime import datetime
from pathlib import Path


def nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_read(path):
    return Path(path).read_text().strip() if nonempty(path) else None


def parse_quast(report_html):
    tsv = Path(report_html).parent / "report.tsv"
    if not tsv.exists():
        return {"contigs": "-", "total_length": 0, "N50": 0}
    rows = [l.split("\t") for l in tsv.read_text().splitlines()]
    data = {r[0].strip(): r[1].strip().replace(",", "") for r in rows if len(r) >= 2}
    return {
        "contigs":      data.get("# contigs", "-"),
        "total_length": int(data.get("Total length") or 0),
        "N50":          int(data.get("N50") or 0),
    }


def parse_core_txt(path):
    if not nonempty(path):
        return []
    rows = []
    for line in Path(path).read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 7 and parts[0] != "ID":
            try:
                rows.append({
                    "id": parts[0], "length": int(parts[1]),
                    "aligned": int(parts[2]), "unaligned": int(parts[3]),
                    "variant": int(parts[4]), "het": int(parts[5]),
                    "lowcov": int(parts[7]) if len(parts) > 7 else 0,
                })
            except ValueError:
                pass
    return rows


def count_core_snps(path):
    return max(0, len(Path(path).read_text().splitlines()) - 1) if nonempty(path) else None


def count_gubbins_snps(path):
    if not nonempty(path):
        return None
    for line in Path(path).read_text().splitlines():
        if not line.startswith(">"):
            return len(line.strip())
    return None


def collect_tool_versions():
    def ver(cmd):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            out = (r.stdout + r.stderr).strip()
            first = out.splitlines()[0] if out else ""
            m = re.search(r"(\d+\.\d+(?:\.\d+)*)", first)
            return m.group(1) if m else (first or "-")
        except Exception:
            return "-"

    return {
        "fastp":     ver(["pixi", "run", "fastp", "--version"]),
        "FastQC":    ver(["pixi", "run", "fastqc", "--version"]),
        "Shovill":   ver(["pixi", "run", "shovill", "--version"]),
        "Prokka":    ver(["pixi", "run", "--environment", "prokka", "prokka", "--version"]),
        "QUAST":     ver(["pixi", "run", "quast", "--version"]),
        "Skani":     ver(["pixi", "run", "--environment", "identification", "skani", "--version"]),
        "Snippy":    ver(["pixi", "run", "snippy", "--version"]),
        "Gubbins":   ver(["pixi", "run", "--environment", "gubbins", "run_gubbins.py", "--version"]),
        "snp-dists": ver(["pixi", "run", "snp-dists", "-v"]),
        "IQ-TREE":   ver(["pixi", "run", "iqtree", "--version"]),
        "MultiQC":   ver(["pixi", "run", "multiqc", "--version"]),
    }


def parse_gubbins_gff(path):
    if not nonempty(path):
        return None
    regions, bases = 0, 0
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 5:
            try:
                bases += int(parts[4]) - int(parts[3]) + 1
                regions += 1
            except ValueError:
                pass
    return {"regions": regions, "bases_masked": bases}


def parse_snp_dists(path):
    txt = safe_read(path)
    if not txt:
        return {"samples": [], "matrix": []}
    rows = list(csv.reader(txt.splitlines(), delimiter="\t"))
    matrix = []
    for row in rows[1:]:
        try:
            matrix.append([int(v) for v in row[1:]])
        except (ValueError, IndexError):
            pass
    return {"samples": rows[0][1:] if rows else [], "matrix": matrix}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--phylo-dir",         required=True)
    ap.add_argument("--samples",           required=True)
    ap.add_argument("--template",          required=True)
    ap.add_argument("--output",            required=True)
    ap.add_argument("--ref",               default="-")
    ap.add_argument("--mincov",            type=int,   default=10)
    ap.add_argument("--minfrac",           type=float, default=0.9)
    ap.add_argument("--snp-dists-gubbins", default=None)
    ap.add_argument("--core-txt",          default=None)
    ap.add_argument("--core-tab",          default=None)
    ap.add_argument("--gubbins-fasta",     default=None)
    ap.add_argument("--gubbins-gff",       default=None)
    ap.add_argument("--snp-threshold",     type=int, default=0)
    args = ap.parse_args()

    pd      = Path(args.phylo_dir)
    samples = args.samples.split(",")

    ref_path  = Path(args.ref)
    ref_auto  = ref_path.name == "reference.gbk"
    ref_label = ref_path.parent.parent.name if ref_auto else ref_path.stem

    sample_data = [
        {"id": s,
         "species": (safe_read(pd / s / "ID_Skani/species.txt") or "-").replace("_", " "),
         **parse_quast(pd / s / "QUAST/report.html")}
        for s in samples
    ]

    data = {
        "run_date":          datetime.now().strftime("%Y-%m-%d %H:%M"),
        "phylo_dir":         str(pd),
        "reference":         ref_label,
        "ref_path":          str(args.ref),
        "ref_auto":          ref_auto,
        "n_samples":         len(samples),
        "snippy_mincov":     args.mincov,
        "snippy_minfrac":    args.minfrac,
        "core_coverage":     parse_core_txt(args.core_txt) if args.core_txt else [],
        "core_snp_count":    count_core_snps(args.core_tab) if args.core_tab else None,
        "gubbins_snp_count": count_gubbins_snps(args.gubbins_fasta) if args.gubbins_fasta else None,
        "samples":           sample_data,
        "snp_dists":         parse_snp_dists(pd / "SNP_Dists/snp_dists.tsv"),
        "snp_dists_gubbins": parse_snp_dists(args.snp_dists_gubbins) if args.snp_dists_gubbins else {"samples": [], "matrix": []},
        "iqtree":            safe_read(pd / "IQtree/iqtree.treefile") or "",
        "tool_versions":     collect_tool_versions(),
        "gubbins_recomb":    parse_gubbins_gff(args.gubbins_gff) if args.gubbins_gff else None,
        "snp_threshold":     args.snp_threshold or None,
    }

    Path(args.output).write_text(
        Path(args.template).read_text().replace("__PHYLO_JSON_DATA__", json.dumps(data, ensure_ascii=False))
    )
