#!/usr/bin/env python3
"""Generate HTML report for the phylogenomics pipeline run."""
import argparse, csv, json, re, subprocess
from datetime import datetime
from pathlib import Path


def nonempty(path):
    if not path:
        return False
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_read(path):
    return Path(path).read_text().strip() if nonempty(path) else None


def parse_checkm(path):
    """Read CheckM lineage_wf --tab_table; return completeness/contamination/strain_het or Nones."""
    txt = safe_read(path)
    empty = {"completeness": None, "contamination": None, "strain_het": None}
    if not txt:
        return empty
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return empty
    r = rows[0]
    def _f(k):
        try: return float(r.get(k, "") or "")
        except (TypeError, ValueError): return None
    return {
        "completeness":  _f("Completeness"),
        "contamination": _f("Contamination"),
        "strain_het":    _f("Strain heterogeneity"),
    }


def parse_quast(quast_dir):
    tsv = Path(quast_dir) / "report.tsv"
    if not tsv.exists():
        return {"contigs": "-", "total_length": 0, "N50": 0, "largest_contig": 0, "gc_pct": None}
    rows = [l.split("\t") for l in tsv.read_text().splitlines() if "\t" in l]
    d    = {r[0].strip(): r[1].strip().replace(",", "") for r in rows}
    gc   = d.get("GC (%)", "")
    return {
        "contigs":        d.get("# contigs", "-"),
        "total_length":   int(d.get("Total length")  or 0),
        "N50":            int(d.get("N50")            or 0),
        "largest_contig": int(d.get("Largest contig") or 0),
        "gc_pct":         float(gc) if gc else None,
    }


def parse_fastp(json_path):
    if not nonempty(json_path):
        return {"reads_raw": None, "reads_trimmed": None, "q30_pct": None, "dup_pct": None}
    summary = json.loads(Path(json_path).read_text()).get("summary", {})
    before  = summary.get("before_filtering", {})
    after   = summary.get("after_filtering",  {})
    return {
        "reads_raw":     before.get("total_reads"),
        "reads_trimmed": after.get("total_reads"),
        "q30_pct":       round(after.get("q30_rate",          0) * 100, 1),
        "dup_pct":       round(before.get("duplication_rate", 0) * 100, 1),
    }


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


def parse_skani_ani(tsv_path):
    if not nonempty(tsv_path):
        return None
    lines = Path(tsv_path).read_text().splitlines()
    try:
        return round(float(lines[1].split("\t")[2]), 2)
    except (IndexError, ValueError):
        return None


def get_mean_depth(bam_path):
    """Return (mean depth, breadth %) from samtools coverage, or (None, None)."""
    if not Path(bam_path).exists():
        return None, None
    try:
        r = subprocess.run(
            ["pixi", "run", "samtools", "coverage", str(bam_path)],
            capture_output=True, text=True, timeout=120,
        )
    except Exception:
        return None, None
    total_len = total_weighted = total_covered = 0
    for line in r.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 7:
            continue
        try:
            length          = int(p[2]) - int(p[1]) + 1
            total_len      += length
            total_weighted += float(p[6]) * length
            total_covered  += int(p[4])
        except (IndexError, ValueError):
            pass
    if total_len == 0:
        return None, None
    return (round(total_weighted / total_len, 1),
            round(total_covered  / total_len * 100, 1))


def parse_core_txt(path):
    """Per-sample stats from snippy-core's core.txt."""
    if not nonempty(path):
        return []
    rows = []
    for line in Path(path).read_text().splitlines():
        p = line.split("\t")
        if len(p) < 7 or p[0] == "ID":
            continue
        try:
            rows.append({
                "id":        p[0],
                "length":    int(p[1]),
                "aligned":   int(p[2]),
                "unaligned": int(p[3]),
                "variant":   int(p[4]),
                "het":       int(p[5]),
                "lowcov":    int(p[7]) if len(p) > 7 else 0,
            })
        except ValueError:
            pass
    return rows


def count_core_snps(path):
    return max(0, len(Path(path).read_text().splitlines()) - 1) if nonempty(path) else None


def count_gubbins_snps(path):
    """First non-header line of the Gubbins polymorphic FASTA = sample's SNP-only sequence."""
    if not nonempty(path):
        return None
    for line in Path(path).read_text().splitlines():
        if not line.startswith(">"):
            return len(line.strip())
    return None


def count_parsnp_snps(vcf_path):
    if not nonempty(vcf_path):
        return None
    return sum(1 for l in Path(vcf_path).read_text().splitlines()
               if l and not l.startswith("#"))


def parse_gubbins_gff(path):
    if not nonempty(path):
        return None
    regions = bases = 0
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) >= 5:
            try:
                bases   += int(p[4]) - int(p[3]) + 1
                regions += 1
            except ValueError:
                pass
    return {"regions": regions, "bases_masked": bases}


def parse_snp_dists(path):
    """Parse snp-dists TSV → {samples, matrix}."""
    txt = safe_read(path)
    if not txt:
        return {"samples": [], "matrix": []}
    rows   = list(csv.reader(txt.splitlines(), delimiter="\t"))
    matrix = []
    for row in rows[1:]:
        try:
            matrix.append([int(v) for v in row[1:]])
        except (IndexError, ValueError):
            pass
    return {"samples": rows[0][1:] if rows else [], "matrix": matrix}


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
        "fastp":     ver(["pixi", "run", "fastp",    "--version"]),
        "FastQC":    ver(["pixi", "run", "fastqc",   "--version"]),
        "Shovill":   ver(["pixi", "run", "shovill",  "--version"]),
        "Prokka":    ver(["pixi", "run", "--environment", "prokka",         "prokka",         "--version"]),
        "QUAST":     ver(["pixi", "run", "quast",    "--version"]),
        "Skani":     ver(["pixi", "run", "--environment", "identification", "skani",          "--version"]),
        "Snippy":    ver(["pixi", "run", "snippy",   "--version"]),
        "Gubbins":   ver(["pixi", "run", "--environment", "gubbins",        "run_gubbins.py", "--version"]),
        "snp-dists": ver(["pixi", "run", "snp-dists", "-v"]),
        "IQ-TREE":   ver(["pixi", "run", "iqtree",   "--version"]),
        "Parsnp":    ver(["pixi", "run", "--environment", "parsnp",         "parsnp",         "--version"]),
        "MultiQC":   ver(["pixi", "run", "multiqc",  "--version"]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phylo-dir",         required=True)
    ap.add_argument("--samples",           required=True)
    ap.add_argument("--template",          required=True)
    ap.add_argument("--output",            required=True)
    ap.add_argument("--ref",               default="-")
    ap.add_argument("--mincov",            type=int,   default=10)
    ap.add_argument("--minfrac",           type=float, default=0.9)
    ap.add_argument("--snp-threshold",     type=int,   default=0)
    ap.add_argument("--snp-dists-gubbins", default=None)
    ap.add_argument("--snp-dists-parsnp",  default=None)
    ap.add_argument("--core-txt",          default=None)
    ap.add_argument("--core-tab",          default=None)
    ap.add_argument("--gubbins-fasta",     default=None)
    ap.add_argument("--gubbins-gff",       default=None)
    ap.add_argument("--parsnp-tree",       default=None)
    ap.add_argument("--parsnp-vcf",        default=None)
    args = ap.parse_args()

    pd      = Path(args.phylo_dir)
    samples = args.samples.split(",")

    # Reference label: auto-annotated Prokka refs live at <dir>/Prokka/reference.gbk;
    # external refs come in as <name>.gbk.
    ref_path  = Path(args.ref)
    ref_auto  = ref_path.name == "reference.gbk"
    ref_label = ref_path.parent.parent.name if ref_auto else ref_path.stem

    # Both snippy-core (via flat symlinks) and Parsnp (via copied FASTAs) emit
    # underscore-form labels (19-03-2026_005a). Map them back to the leaf name
    # so every table in the report displays the same short ID ("005a").
    to_leaf = {s.replace("/", "_"): s.split("/")[-1] for s in samples}

    def normalise(dists):
        return {**dists, "samples": [to_leaf.get(s, s) for s in dists["samples"]]}

    # Kraken2/Bracken is not part of the phylo pipeline; fall back to the
    # main pipeline's results/ dir (sibling of results_phylo/) if available.
    main_results = pd.parent.parent / "results"

    def kraken_pct(sample):
        for base in (pd, main_results):
            p = base / sample / "ID_Kraken2" / "bracken_species.txt"
            if p.exists():
                return parse_bracken_primary_pct(p)
        return None

    sample_data = []
    for s in samples:
        depth, breadth = get_mean_depth(pd / "Snippy" / s / "snps.bam")
        species_raw    = safe_read(pd / s / "ID_Skani/species.txt") or "-"
        sample_data.append({
            "id":                 s.split("/")[-1],
            "species":            species_raw.replace("_", " "),
            "ani":                parse_skani_ani(pd / s / "ID_Skani/skani.tsv"),
            "depth":              depth,
            "breadth":            breadth,
            "kraken_primary_pct": kraken_pct(s),
            **parse_quast(pd / s / "QUAST"),
            **parse_fastp(pd / s / "QC/fastp.json"),
            **parse_checkm(pd / s / "CheckM/quality.tsv"),
        })

    core_coverage = parse_core_txt(args.core_txt)
    for row in core_coverage:
        row["id"] = to_leaf.get(row["id"], row["id"])

    empty_matrix = {"samples": [], "matrix": []}
    data = {
        "run_date":          datetime.now().strftime("%Y-%m-%d %H:%M"),
        "phylo_dir":         str(pd),
        "reference":         ref_label,
        "ref_path":          str(args.ref),
        "ref_auto":          ref_auto,
        "n_samples":         len(samples),
        "snippy_mincov":     args.mincov,
        "snippy_minfrac":    args.minfrac,
        "snp_threshold":     args.snp_threshold or None,
        "core_snp_count":    count_core_snps(args.core_tab),
        "gubbins_snp_count": count_gubbins_snps(args.gubbins_fasta),
        "gubbins_recomb":    parse_gubbins_gff(args.gubbins_gff),
        "parsnp_snp_count":  count_parsnp_snps(args.parsnp_vcf),
        "core_coverage":     core_coverage,
        "samples":           sample_data,
        "snp_dists":         normalise(parse_snp_dists(pd / "SNP_Dists/snp_dists.tsv")),
        "snp_dists_gubbins": normalise(parse_snp_dists(args.snp_dists_gubbins)) if args.snp_dists_gubbins else empty_matrix,
        "snp_dists_parsnp":  normalise(parse_snp_dists(args.snp_dists_parsnp))  if args.snp_dists_parsnp  else empty_matrix,
        "iqtree":            safe_read(pd / "IQtree/iqtree.treefile") or "",
        "parsnp_tree":       safe_read(args.parsnp_tree) if args.parsnp_tree else None,
        "tool_versions":     collect_tool_versions(),
    }

    Path(args.output).write_text(
        Path(args.template).read_text().replace(
            "__PHYLO_JSON_DATA__", json.dumps(data, ensure_ascii=False)
        )
    )


if __name__ == "__main__":
    main()
