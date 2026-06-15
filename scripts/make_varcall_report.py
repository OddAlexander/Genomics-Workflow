#!/usr/bin/env python3
"""Generate per-sample HTML report for the variant-calling pipeline.

Reference picked from Mash screen, reads mapped with bowtie2, variants called
by bcftools, optionally annotated with bcftools csq. Pulls "real" depth/breadth
from samtools coverage and mapping rate from samtools flagstat -- the figures
the user explicitly wanted to surface.
"""
import argparse, csv, gzip, json, re, subprocess
from datetime import datetime
from pathlib import Path


def nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_read(path):
    return Path(path).read_text().strip() if nonempty(path) else None


def parse_fastp(json_path):
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


def parse_samtools_coverage(path):
    """samtools coverage TSV -> weighted mean depth + breadth % over all contigs."""
    out = {"mean_depth": None, "breadth_pct": None, "covered_bases": 0, "total_bases": 0,
           "contigs": []}
    if not nonempty(path):
        return out
    total_len = total_weighted = total_covered = 0
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 7:
            continue
        try:
            length     = int(p[2]) - int(p[1]) + 1
            covered    = int(p[4])
            mean_depth = float(p[6])
        except (IndexError, ValueError):
            continue
        total_len      += length
        total_weighted += mean_depth * length
        total_covered  += covered
        out["contigs"].append({
            "name":       p[0],
            "length":     length,
            "covered":    covered,
            "mean_depth": round(mean_depth, 1),
            "breadth":    round(covered / length * 100, 2) if length else 0,
        })
    if total_len:
        out.update({
            "mean_depth":     round(total_weighted / total_len, 1),
            "breadth_pct":    round(total_covered / total_len * 100, 2),
            "covered_bases":  total_covered,
            "total_bases":    total_len,
        })
    return out


def parse_flagstat(path):
    out = {"total_reads": None, "mapped": None, "mapped_pct": None,
           "duplicates": None, "properly_paired_pct": None}
    if not nonempty(path):
        return out
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if "in total" in line:
            try:
                out["total_reads"] = int(line.split()[0])
            except ValueError:
                pass
        elif " duplicates" in line and "primary" not in line:
            try:
                out["duplicates"] = int(line.split()[0])
            except ValueError:
                pass
        elif " mapped (" in line and "primary" not in line and "mate" not in line:
            try:
                out["mapped"] = int(line.split()[0])
            except ValueError:
                pass
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                out["mapped_pct"] = float(m.group(1))
        elif "properly paired" in line:
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                out["properly_paired_pct"] = float(m.group(1))
    return out


_SPECIES_RE = re.compile(r"\b([A-Z][a-z]+)\s+([a-z]+)\b")


def parse_mash_screen(path, top_n=10):
    """Top-N Mash screen hits, each with parsed Genus species."""
    rows = []
    if not nonempty(path):
        return rows
    for line in Path(path).read_text().splitlines()[:top_n]:
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        try:
            identity = float(parts[0])
        except ValueError:
            continue
        comment = parts[-1].strip()
        m = _SPECIES_RE.search(comment)
        rows.append({
            "identity":      round(identity, 4),
            "shared_hashes": parts[1] if len(parts) > 1 else "-",
            "median_mult":   parts[2] if len(parts) > 2 else "-",
            "p_value":       parts[3] if len(parts) > 3 else "-",
            "query":         parts[4] if len(parts) > 4 else "-",
            "species":       f"{m.group(1)} {m.group(2)}" if m else "-",
            "comment":       comment,
        })
    return rows


def parse_vcf_summary(vcf_path):
    """bcftools stats summary -> {snps, indels, n_records, titv}."""
    out = {"snps": None, "indels": None, "n_records": None, "titv": None,
           "n_pass": None, "n_filtered": None}
    try:
        r = subprocess.run(
            ["pixi", "run", "bcftools", "stats", str(vcf_path)],
            capture_output=True, text=True, timeout=120,
        )
        text = r.stdout
    except Exception:
        return out
    for line in text.splitlines():
        if line.startswith("SN") and "number of SNPs:" in line:
            try: out["snps"] = int(line.split(":")[-1].strip())
            except ValueError: pass
        elif line.startswith("SN") and "number of indels:" in line:
            try: out["indels"] = int(line.split(":")[-1].strip())
            except ValueError: pass
        elif line.startswith("SN") and "number of records:" in line:
            try: out["n_records"] = int(line.split(":")[-1].strip())
            except ValueError: pass
        elif line.startswith("TSTV"):
            parts = line.split("\t")
            try: out["titv"] = float(parts[4])
            except (IndexError, ValueError): pass
    return out


def _open_vcf(path):
    """Return a text-mode reader for a .vcf.gz or plain .vcf."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return p.open()


def parse_vcf_records(vcf_path, limit=300):
    """Extract variant records (PASS or unfiltered) with optional BCSQ consequence."""
    records = []
    n_pass = n_filt = 0
    try:
        reader = _open_vcf(vcf_path)
    except Exception:
        return records, n_pass, n_filt
    with reader as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            filt = parts[6]
            if filt in (".", "PASS"):
                n_pass += 1
            else:
                n_filt += 1
            if len(records) >= limit:
                continue
            info = parts[7]
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                else:
                    info_d[kv] = True
            csq_str = info_d.get("BCSQ", "")
            cons, gene, change = "-", "-", "-"
            if csq_str and csq_str != "-" and isinstance(csq_str, str):
                first = csq_str.split(",")[0]
                cf = first.split("|")
                if len(cf) >= 2:
                    cons = cf[0] or "-"
                    gene = cf[1] or "-"
                if len(cf) >= 7:
                    change = cf[6] or "-"   # e.g. p.Ser83Leu
            records.append({
                "chrom":  parts[0],
                "pos":    parts[1],
                "ref":    parts[3],
                "alt":    parts[4],
                "qual":   parts[5],
                "filter": filt,
                "depth":  info_d.get("DP", "-"),
                "csq":    cons,
                "gene":   gene,
                "change": change,
            })
    return records, n_pass, n_filt


def parse_bracken_primary_pct(path):
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


def parse_reference_tsv(path):
    """reference.tsv has a header row + one data row: species, fasta, gff."""
    if not nonempty(path):
        return {"species": "-", "fasta": "-", "gff": None}
    lines = Path(path).read_text().strip().splitlines()
    if len(lines) < 2:
        return {"species": "-", "fasta": "-", "gff": None}
    sp, fasta, gff = (lines[1].split("\t") + ["-", "-", "-"])[:3]
    return {"species": sp, "fasta": fasta, "gff": (gff if gff and gff != "-" else None)}


_TOOL_DISPLAY = {
    "fastp":    "fastp",
    "fastqc":   "FastQC",
    "mash":     "Mash",
    "bowtie2":  "Bowtie2",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "multiqc":  "MultiQC",
    "snakemake": "Snakemake",
}

_TOOL_ENV = {
    "fastp":    "default",
    "fastqc":   "default",
    "mash":     "identification",
    "bowtie2":  "default",
    "samtools": "default",
    "bcftools": "default",
    "multiqc":  "default",
    "snakemake": "default",
}


def collect_versions():
    """Read tool versions from pixi conda-meta filenames. Returns list of
    {tool, version} dicts for the tools actually installed."""
    pixi_envs = Path(__file__).parent.parent / ".pixi" / "envs"
    versions = []
    seen = set()
    for pkg, env in _TOOL_ENV.items():
        meta_dir = pixi_envs / env / "conda-meta"
        if not meta_dir.is_dir():
            continue
        prefix = pkg + "-"
        for f in meta_dir.iterdir():
            if f.name.startswith(prefix) and f.suffix == ".json":
                parts = f.stem.split("-")
                pkg_parts = pkg.split("-")
                version = parts[len(pkg_parts)] if len(parts) > len(pkg_parts) else "?"
                display = _TOOL_DISPLAY.get(pkg, pkg)
                if display not in seen:
                    versions.append({"tool": display, "version": version})
                    seen.add(display)
                break
    return versions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample",      required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--template",    required=True)
    ap.add_argument("--ref-tsv",     required=True)
    ap.add_argument("--mash-screen", required=True)
    ap.add_argument("--coverage",    required=True)
    ap.add_argument("--flagstat",    required=True)
    ap.add_argument("--stats",       required=True)
    ap.add_argument("--vcf",         required=True)
    ap.add_argument("--fastp",       required=True)
    ap.add_argument("--bracken",     required=True)
    ap.add_argument("--output",      required=True)
    args = ap.parse_args()

    ref          = parse_reference_tsv(args.ref_tsv)
    fastp_d      = parse_fastp(args.fastp)
    cov          = parse_samtools_coverage(args.coverage)
    flag         = parse_flagstat(args.flagstat)
    mash_top     = parse_mash_screen(args.mash_screen)
    vcf_summary  = parse_vcf_summary(args.vcf)
    vcf_records, n_pass, n_filt = parse_vcf_records(args.vcf)
    bracken_pct  = parse_bracken_primary_pct(args.bracken)
    vcf_summary.update({"n_pass": n_pass, "n_filtered": n_filt})

    data = {
        "sample":        args.sample,
        "completed":     datetime.now().strftime("%Y-%m-%d %H:%M"),
        "reference":     {**ref, "species_display": ref["species"].replace("_", " ")},
        "mash":          mash_top,
        "fastp":         fastp_d,
        "coverage":      cov,
        "flagstat":      flag,
        "vcf_summary":   vcf_summary,
        "vcf_records":   vcf_records,
        "bracken_pct":   bracken_pct,
        "versions": collect_versions(),
    }

    Path(args.output).write_text(
        Path(args.template).read_text()
            .replace("__VARCALL_JSON_DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__SAMPLE_NAME__", args.sample)
    )


if __name__ == "__main__":
    main()
