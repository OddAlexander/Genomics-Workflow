#!/usr/bin/env python3
"""Generate HTML report for phylogenomics pipeline run."""
import argparse, json, csv
from datetime import datetime
from pathlib import Path


def safe_read(path):
    p = Path(path)
    return p.read_text().strip() if p.exists() and p.stat().st_size > 0 else None


def parse_quast(report_html):
    tsv = Path(report_html).parent / "report.tsv"
    if not tsv.exists():
        return {"contigs": "-", "total_length": 0, "N50": 0}
    data = {}
    for line in tsv.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            data[parts[0].strip()] = parts[1].strip().replace(",", "")
    return {
        "contigs":      data.get("# contigs", "-"),
        "total_length": int(data.get("Total length", 0) or 0),
        "N50":          int(data.get("N50", 0) or 0),
    }


def parse_snp_dists(path):
    txt = safe_read(path)
    if not txt:
        return {"samples": [], "matrix": []}
    rows = list(csv.reader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return {"samples": [], "matrix": []}
    header = rows[0][1:]
    matrix = []
    for row in rows[1:]:
        if row:
            try:
                matrix.append([int(v) for v in row[1:]])
            except (ValueError, IndexError):
                pass
    return {"samples": header, "matrix": matrix}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--phylo-dir",          required=True)
    ap.add_argument("--samples",            required=True)
    ap.add_argument("--template",           required=True)
    ap.add_argument("--output",             required=True)
    ap.add_argument("--ref",                default="-")
    ap.add_argument("--snp-dists-gubbins",  default=None)
    args = ap.parse_args()

    pd      = Path(args.phylo_dir)
    samples = args.samples.split(",")

    sample_data = []
    for s in samples:
        sp_dir  = pd / s
        species = safe_read(sp_dir / "ID_Skani/species.txt") or "-"
        quast   = parse_quast(sp_dir / "QUAST/report.html")
        sample_data.append({
            "id":      s,
            "species": species.replace("_", " "),
            **quast,
        })

    snp_data         = parse_snp_dists(pd / "SNP_Dists/snp_dists.tsv")
    snp_data_gubbins = parse_snp_dists(args.snp_dists_gubbins) if args.snp_dists_gubbins else {"samples": [], "matrix": []}
    iqtree           = safe_read(pd / "IQtree/iqtree.treefile") or ""

    ref_path = Path(args.ref)
    if ref_path.name == "reference.gbk":
        ref_label = ref_path.parent.parent.name
        ref_auto  = True
    else:
        ref_label = ref_path.stem
        ref_auto  = False

    data = {
        "run_date":   datetime.now().strftime("%Y-%m-%d %H:%M"),
        "phylo_dir":  str(pd),
        "reference":  ref_label,
        "ref_path":   str(args.ref),
        "ref_auto":   ref_auto,
        "n_samples":  len(samples),
        "samples":   sample_data,
        "snp_dists":        snp_data,
        "snp_dists_gubbins": snp_data_gubbins,
        "iqtree":    iqtree,
    }

    template = Path(args.template).read_text()
    html     = template.replace("__PHYLO_JSON_DATA__", json.dumps(data, ensure_ascii=False))
    Path(args.output).write_text(html)
