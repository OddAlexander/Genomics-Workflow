#!/usr/bin/env python3
"""Generate HTML report for a single sample from pipeline outputs."""
import argparse, csv, fcntl, gzip, json, re
from datetime import datetime
from pathlib import Path


def nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_read(path):
    return Path(path).read_text().strip() if nonempty(path) else None


def tsv_first_row(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    return rows[0] if rows else None


def parse_mlst(path):
    txt = safe_read(path)
    if not txt:
        return {"ST": "-", "scheme": "-"}
    parts = txt.split("\t")
    return {"ST":     parts[2] if len(parts) > 2 else "-",
            "scheme": parts[1] if len(parts) > 1 else "-"}


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


def parse_quast(report_html):
    tsv = Path(report_html).parent / "report.tsv"
    if not tsv.exists():
        return {"contigs": "-", "total_length": 0, "N50": 0, "largest_contig": 0}
    rows = [l.split("\t") for l in tsv.read_text().splitlines()]
    data = {r[0].strip(): r[1].strip().replace(",", "") for r in rows if len(r) >= 2}
    return {
        "contigs":        data.get("# contigs", "-"),
        "total_length":   int(data.get("Total length") or 0),
        "N50":            int(data.get("N50") or 0),
        "largest_contig": int(data.get("Largest contig") or 0),
    }


def parse_fastp(json_path):
    """Read totals from a fastp JSON; returns reads/bases/q30/dup_pct or Nones."""
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


def estimate_depth(bases, genome_size):
    """Approximate mean coverage = trimmed bases / assembled length."""
    if not bases or not genome_size:
        return None
    return round(bases / genome_size, 1)


def parse_bracken(path):
    txt = safe_read(path)
    empty = {"primary_species": "-", "primary_pct": 0.0, "secondary_species": "-", "secondary_pct": 0.0}
    if not txt:
        return empty
    frac = lambda r: float(r.get("fraction_total_reads") or 0)
    rows = sorted(csv.DictReader(txt.splitlines(), delimiter="\t"), key=frac, reverse=True)
    if not rows:
        return empty
    sp = lambda r: r.get("name", "-").replace(" ", "_")
    return {
        "primary_species":   sp(rows[0]),
        "primary_pct":       frac(rows[0]) * 100,
        "secondary_species": sp(rows[1]) if len(rows) > 1 else "-",
        "secondary_pct":     frac(rows[1]) * 100 if len(rows) > 1 else 0.0,
    }


def parse_kraken2_unclassified(path):
    txt = safe_read(path)
    if not txt:
        return 0.0
    for line in txt.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 6 and parts[3].strip() == "U":
            try:
                return float(parts[0].strip())
            except ValueError:
                pass
    return 0.0


def parse_amrfinder(path):
    txt = safe_read(path)
    if not txt:
        return {"total_amr_genes": 0, "total_virulence_genes": 0, "total_stress_genes": 0, "_rows": []}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    return {
        "total_amr_genes":       sum(1 for r in rows if r.get("Type") == "AMR"),
        "total_virulence_genes": sum(1 for r in rows if r.get("Type") == "VIRULENCE"),
        "total_stress_genes":    sum(1 for r in rows if r.get("Type") == "STRESS"),
        "_rows": rows,
    }


def parse_mobsuite(path):
    txt = safe_read(path)
    if not txt:
        return {"count": 0, "conjugative": 0, "mobilizable": 0,
                "non_mobilizable": 0, "details": [], "_contig_mob": {}}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    mob  = lambda r: r.get("predicted_mobility", "").upper()
    contig_mob = {
        r.get("file_id", ""): {
            "mobility": mob(r),
            "replicon": r.get("rep_type(s)", "-"),
            "cluster":  r.get("cluster_id", "-"),
        }
        for r in rows
    }
    return {
        "count":           len(rows),
        "conjugative":     sum(1 for r in rows if "CONJUGATIVE"     in mob(r)),
        "mobilizable":     sum(1 for r in rows if mob(r) == "MOBILIZABLE"),
        "non_mobilizable": sum(1 for r in rows if "NON_MOBILIZABLE" in mob(r)),
        "details": [{
            "file_id":    r.get("file_id", "-"),
            "replicon":   r.get("rep_type(s)", "-"),
            "relaxase":   r.get("relaxase_type(s)", "-"),
            "mobility":   mob(r),
            "host_range": r.get("mash_nearest_neighbor", "-"),
            "cluster":    r.get("cluster_id", "-"),
            "size":       r.get("size", 0),
        } for r in rows],
        "_contig_mob": contig_mob,
    }


def build_amr_details(amr_rows, contig_mob, contig_mge, contig_plasmid):
    details, high_risk, medium_risk = [], [], []
    for g in amr_rows:
        contig       = g.get("Contig id", "-")
        mob_info     = contig_mob.get(contig, {})
        mobility     = mob_info.get("mobility", "-")
        gene         = g.get("Element symbol", "-")
        element_type = g.get("Type", "-")
        plasmids     = contig_plasmid.get(contig, [])
        if element_type == "AMR":
            if mobility == "CONJUGATIVE" or plasmids:
                high_risk.append(gene)
            elif mobility == "MOBILIZABLE":
                medium_risk.append(gene)
        details.append({
            "gene":      gene,
            "gene_name": g.get("Element name", "-"),
            "type":      element_type,
            "subtype":   g.get("Subtype", "-"),
            "class":     g.get("Class", "-"),
            "subclass":  g.get("Subclass", "-"),
            "scope":     g.get("Scope", "-"),
            "location":  contig,
            "mobility":  mobility,
            "replicon":  mob_info.get("replicon", "-"),
            "cluster":   mob_info.get("cluster", "-"),
            "mge":       contig_mge.get(contig, []),
            "plasmids":  plasmids,
        })
    return details, high_risk, medium_risk


def parse_staphscope(path):
    r = tsv_first_row(path)
    if not r:
        return None
    return {
        "spa_type":    r.get("spa_type",    "-") or "-",
        "sccmec_type": r.get("sccmec_type", "-") or "-",
        "mrsa":        r.get("mrsa_status", "").upper() == "MRSA",
    }


def parse_kleborate(path):
    r = tsv_first_row(path)
    if not r:
        return None

    def _col(*suffixes):
        for suffix in suffixes:
            for k, v in r.items():
                if k.endswith(suffix):
                    return v or "-"
        return "-"

    def _st():
        for k, v in r.items():
            if k.endswith("__mlst__ST") or k.endswith("__mlst_achtman__ST"):
                return (v or "-").lstrip("ST") if v and v not in ("-", "") else "-"
        return "-"

    result = {
        "ST":              _st(),
        "K_type":          _col("__K_locus", "__K_type", "__wzi"),
        "O_type":          _col("__O_locus", "__O_type", "__ectyper__O-type"),
        "virulence_score": _col("__virulence_score"),
        "Bla_Carb":        _col("__Carbapenem", "__Bla_Carb_acquired"),
        "Bla_ESBL":        _col("__Cephalosporin", "__Bla_ESBL_acquired"),
    }
    return None if all(v in ("-", "", None) for v in result.values()) else result


def parse_ectyper(path):
    r = tsv_first_row(path)
    return {"serotype": f"{r.get('O-type', '-')}:{r.get('H-type', '-')}"} if r else None


def parse_emmtyper(path):
    txt = safe_read(path)
    if not txt:
        return None
    for line in txt.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            return {"emm_type": parts[1] or "-"}
    return None


def parse_skani(path):
    txt  = safe_read(path)
    empty = {"status": "skipped", "top_hit": "-", "ani": None, "af_query": None, "af_ref": None}
    if not txt:
        return empty
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows or rows[0].get("status") == "skipped":
        return empty
    r   = rows[0]
    ani = r.get("ANI")
    return {
        "status":   "ok",
        "top_hit":  r.get("Ref_name", r.get("Ref_file", "-").split("/")[-1].replace(".fna", "").replace(".fa", "")),
        "ani":      float(ani) if ani else None,
        "af_query": r.get("Align_fraction_query", "-"),
        "af_ref":   r.get("Align_fraction_ref", "-"),
    }


def parse_seqsero2(path):
    r = tsv_first_row(path)
    if not r:
        return None
    return {
        "serotype":          r.get("Predicted serotype",          "-") or "-",
        "antigenic_profile": r.get("Predicted antigenic profile",  "-") or "-",
        "O_antigen":         r.get("O antigen prediction",         "-") or "-",
        "comment":           r.get("Note", "") or "",
    }


def parse_hicap(path):
    r = tsv_first_row(path)
    if not r:
        return None
    return {
        "serotype":     r.get("predicted_serotype", "-") or "-",
        "genes":        r.get("genes_identified",   "-") or "-",
        "completeness": r.get("percent_coverage", r.get("completeness", "-")) or "-",
    }


def parse_gambit(path):
    txt = safe_read(path)
    empty = {"status": "skipped", "predicted_name": "-", "closest_description": "-", "closest_distance": None}
    if not txt:
        return empty
    rows = list(csv.DictReader(txt.splitlines()))
    if not rows:
        return empty
    r    = rows[0]
    name = (r.get("predicted.name") or "").strip()
    dist = r.get("closest.distance", "")
    return {
        "status":              "ok" if name else "no_match",
        "predicted_name":      name or "-",
        "closest_description": (r.get("closest.description") or "-").strip(),
        "closest_distance":    float(dist) if dist else None,
    }


def parse_gbsserotyper(path):
    """GBS-SBG TSV (Name / Serotype / Uncertainty) -> serotype + uncertainty."""
    txt = safe_read(path)
    if not txt:
        return None
    for line in txt.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            serotype    = parts[1].strip() or "NT"
            uncertainty = parts[2].strip() if len(parts) > 2 else ""
            return {"serotype": serotype, "uncertainty": uncertainty}
    return None


def parse_lrefinder(path):
    """LRE-Finder .res (KMA TSV) -> detected linezolid resistance markers.

    LRE-Finder reports a row for any template with non-zero k-mer overlap,
    including: (a) the wrong-species 23S probe (always low identity), and
    (b) acquired-gene templates that share a few k-mers with unrelated reads.
    Without filtering, those fragmentary hits would flip the report banner
    from green to red. So we apply standard AMR-calling thresholds here:
      - acquired genes (cfr/optrA/poxtA/...): require coverage >= 90% AND identity >= 98%
      - 23S references: require coverage >= 50% AND identity >= 80% (drops cross-species probe)

    For 23S, the *clinical* number is `mut_pct = 100 - identity` -- the fraction
    of 23S rRNA copies carrying a resistance mutation (mosaicism)."""
    txt = safe_read(path)
    if not txt:
        return None
    acquired, ribosomal = [], []
    for line in txt.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        try:
            template, identity, coverage, depth = p[0], float(p[4]), float(p[5]), float(p[8])
        except (IndexError, ValueError):
            continue
        is_23s = "23S" in template or template.lower().startswith("rrl")
        if is_23s:
            if coverage < 50 or identity < 80:
                continue   # cross-species probe, not informative
        else:
            if coverage < 90 or identity < 98:
                continue   # fragmentary k-mer overlap, not a real acquired-gene call
        row = {"template": template, "identity": round(identity, 2),
               "coverage": round(coverage, 2), "depth": round(depth, 1),
               "mut_pct":  round(100.0 - identity, 2) if is_23s else None}
        (ribosomal if is_23s else acquired).append(row)
    if not acquired and not ribosomal:
        return None
    return {
        "acquired_genes": acquired,
        "ribosomal_hits": ribosomal,
        "n_acquired":     len(acquired),
        # Max mutated-copy % across 23S references is the clinically actionable signal.
        "max_23s_mut_pct": max((r["mut_pct"] for r in ribosomal if r["mut_pct"] is not None),
                                default=None),
    }


def parse_plasmidfinder(path):
    txt = safe_read(path)
    if not txt:
        return {"hits": [], "_contig_plasmid": {}}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    contig_plasmid = {}
    hits = []
    for r in rows:
        # Shovill appends metadata to contig name — keep only the contigXXXXX prefix
        contig  = r.get("Contig", "-").split(" ")[0] or "-"
        plasmid = r.get("Plasmid", "-")
        hits.append({
            "plasmid":   plasmid,
            "identity":  r.get("Identity", "-"),
            "contig":    contig,
            "database":  r.get("Database", "-"),
            "accession": r.get("Accession number", "-"),
        })
        contig_plasmid.setdefault(contig, []).append(plasmid)
    return {"hits": hits, "_contig_plasmid": contig_plasmid}


def parse_mefinder(path):
    txt = safe_read(path)
    if not txt:
        return {"count": 0, "elements": []}
    data_lines = [l for l in txt.splitlines() if not l.startswith("#") and l.strip()]
    rows = list(csv.DictReader(data_lines))
    def _g(r, *keys):
        for k in keys:
            if r.get(k): return r[k]
        return "-"
    elements = [{
        "nr":     i,
        "name":   _g(r, "name",   "Name"),
        "type":   _g(r, "type",   "Type"),
        "contig": _g(r, "contig", "contig_id", "sequence_id"),
        "start":  _g(r, "start",  "Start"),
        "end":    _g(r, "end",    "End"),
    } for i, r in enumerate(rows, 1)]
    return {"count": len(elements), "elements": elements}


def _aggregate_samtools_coverage(path):
    """Walk a `samtools coverage` TSV and return length-weighted (mean_depth,
    breadth_pct) across all contigs, or (None, None) if the file is missing
    or has no usable rows. Shared by parse_self_coverage and parse_varcall."""
    if not nonempty(path):
        return (None, None)
    total_len = total_dep = total_cov = 0
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 7:
            continue
        try:
            length = int(p[2]) - int(p[1]) + 1
            total_len += length
            total_dep += float(p[6]) * length
            total_cov += int(p[4])
        except (IndexError, ValueError):
            pass
    if not total_len:
        return (None, None)
    return (round(total_dep / total_len, 1),
            round(total_cov / total_len * 100, 2))


def parse_varcall(sp_dir):
    """Opportunistically surface the varcall pipeline's outputs if it has been
    run for this sample (results/<sample>/VarCall/...). Returns None when
    varcall hasn't run, so the rest of the report renders cleanly without it.

    Provides reference-mapped depth/breadth (complementary to self-mapped depth
    from parse_self_coverage), chosen reference species, variant PASS/filtered
    counts, and a link to the full varcall_report.html."""
    vc_dir          = Path(sp_dir) / "VarCall"
    mean_depth, br  = _aggregate_samtools_coverage(vc_dir / "QC" / "coverage.tsv")
    if mean_depth is None:
        return None

    mapped_pct = None
    flagstat = vc_dir / "QC" / "flagstat.txt"
    if nonempty(flagstat):
        for line in flagstat.read_text().splitlines():
            if " mapped (" in line and "primary" not in line and "mate" not in line:
                m = re.search(r"\(([0-9.]+)%", line)
                if m:
                    mapped_pct = float(m.group(1))
                    break

    ref_species = "-"
    ref_tsv = vc_dir / "Mash" / "reference.tsv"
    if nonempty(ref_tsv):
        lines = ref_tsv.read_text().splitlines()
        if len(lines) >= 2:
            ref_species = lines[1].split("\t")[0]

    n_pass = n_filt = 0
    vcf = vc_dir / "VCF" / "variants.annotated.vcf.gz"
    if not vcf.exists():
        vcf = vc_dir / "VCF" / "variants.vcf.gz"
    if vcf.exists():
        try:
            with gzip.open(vcf, "rt") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 7:
                        continue
                    if parts[6] in ("PASS", "."):
                        n_pass += 1
                    else:
                        n_filt += 1
        except OSError:
            pass

    return {
        "reference_species": ref_species,
        "mean_depth":        mean_depth,
        "breadth_pct":       br,
        "mapping_rate_pct":  mapped_pct,
        "variants_pass":     n_pass,
        "variants_filtered": n_filt,
        "report_link":       "VarCall/varcall_report.html",
    }


def parse_self_coverage(path):
    """samtools coverage on reads self-mapped to the Shovill assembly.
    Returns length-weighted mean depth + breadth, SeqSphere-equivalent metric.
    Shares the aggregation helper with parse_varcall."""
    mean_depth, br = _aggregate_samtools_coverage(path)
    return {"mean_depth": mean_depth, "breadth_pct": br}


def parse_pasty(path):
    """Pasty P. aeruginosa O-antigen serotyping -> serotype string."""
    r = tsv_first_row(path)
    if not r:
        return None
    stype = r.get("type") or r.get("serogroup") or "-"
    return {"serotype": stype or "-"} if stype and stype not in ("-", "") else None


def derive_purity(primary_pct, species_sp, bracken_sp):
    if primary_pct >= 80 and species_sp == bracken_sp:
        return "PURE"
    if primary_pct >= 60:
        return "LIKELY_PURE"
    return "MIXED"


def contamination_check(primary_pct, secondary_pct, secondary_sp, uncl_pct):
    warnings = []
    if secondary_pct >= 10 and secondary_sp not in ("-", "unknown", ""):
        warnings.append(f"Secondary species {secondary_sp.replace('_', ' ')} accounts for {secondary_pct:.1f}%")
    if uncl_pct >= 20:
        warnings.append(f"{uncl_pct:.1f}% unclassified reads")
    if primary_pct < 80:
        warnings.append(f"Primary species below 80% ({primary_pct:.1f}%)")
    return warnings


def concordance_score(species_sp, bracken_sp):
    return 2 if species_sp == bracken_sp and species_sp not in ("unknown", "-", "") else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample",      required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--template",    required=True)
    ap.add_argument("--output",      required=True)
    ap.add_argument("--kraken2-db",  default="-")
    ap.add_argument("--log-path",    default=None)
    args = ap.parse_args()

    rd     = Path(args.results_dir)
    sp_dir = rd / args.sample

    mlst_data    = parse_mlst(sp_dir / "MLST/mlst.tsv")
    quast_data   = parse_quast(sp_dir / "QUAST/report.html")
    quast_data.update(parse_checkm(sp_dir / "CheckM/quality.tsv"))
    fastp_data   = parse_fastp(sp_dir / "QC/fastp.json")
    est_depth    = estimate_depth(fastp_data["bases"], quast_data["total_length"])
    bracken_data = parse_bracken(sp_dir / "ID_Kraken2/bracken_species.txt")
    uncl_pct     = parse_kraken2_unclassified(sp_dir / "ID_Kraken2/kraken2_report.txt")
    amr_raw      = parse_amrfinder(sp_dir / "AMRFinder/amrfinder.tsv")
    mob_data     = parse_mobsuite(sp_dir / "MOBSuite/mobtyper.tsv")
    mge_data     = parse_mefinder(sp_dir / "MEfinder/mefinder.tsv")
    pf_data      = parse_plasmidfinder(sp_dir / "PlasmidFinder/results_tab.tsv")
    skani_data   = parse_skani(sp_dir / "ID_Skani/skani.tsv")
    gambit_data  = parse_gambit(sp_dir / "ID_Gambit/gambit.csv")
    species_sp   = (safe_read(sp_dir / "species.txt") or "unknown").strip()
    self_cov     = parse_self_coverage(sp_dir / "QC/self_coverage.tsv")
    varcall_data = parse_varcall(sp_dir)   # None when varcall hasn't been run for this sample

    contig_mge = {}
    for el in mge_data.get("elements", []):
        c = el.get("contig", "-")
        if c and c != "-":
            contig_mge.setdefault(c, []).append({"name": el["name"], "type": el["type"]})

    amr_details, high_risk, medium_risk = build_amr_details(
        amr_raw.get("_rows", []), mob_data.get("_contig_mob", {}),
        contig_mge, pf_data.get("_contig_plasmid", {}))

    purity  = derive_purity(bracken_data["primary_pct"], species_sp, bracken_data["primary_species"])
    score   = concordance_score(species_sp, bracken_data["primary_species"])
    species = species_sp if species_sp not in ("unknown", "") else bracken_data["primary_species"]
    contam_warnings = contamination_check(
        bracken_data["primary_pct"], bracken_data["secondary_pct"],
        bracken_data["secondary_species"], uncl_pct)

    data = {
        "sample":    args.sample,
        "species":   species,
        "completed": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "mlst":      mlst_data,
        "assembly":  quast_data,
        "qc": {
            "reads":         fastp_data["reads"],
            "q30_pct":       fastp_data["q30_pct"],
            "dup_pct":       fastp_data["dup_pct"],
            "est_depth":     est_depth,                # fastp bases / QUAST length -- upper bound
            "mapped_depth":  self_cov["mean_depth"],   # self-mapped (bwa-mem2 vs Shovill assembly)
            "breadth_pct":   self_cov["breadth_pct"],
        },
        "amr": {
            "total_amr_genes":       amr_raw["total_amr_genes"],
            "total_virulence_genes": amr_raw["total_virulence_genes"],
            "total_stress_genes":    amr_raw["total_stress_genes"],
            "high_risk_genes":       high_risk,
            "medium_risk_genes":     medium_risk,
            "details":               amr_details,
        },
        "plasmids":      {k: v for k, v in mob_data.items() if not k.startswith("_")},
        "plasmidfinder": {"hits": pf_data.get("hits", [])},
        "mge":           mge_data,
        "species_id": {
            "purity":                 purity,
            "concordance_score":      score,
            "contamination_warnings": contam_warnings,
            "bracken_species":        bracken_data["primary_species"],
            "gambit_species":          species_sp,
            "id_method":              "Bracken + GAMBIT + SKANI",
            "primary_species":        bracken_data["primary_species"],
            "primary_pct":            bracken_data["primary_pct"],
            "secondary_species":      bracken_data["secondary_species"],
            "secondary_pct":          bracken_data["secondary_pct"],
            "unclassified_pct":       uncl_pct,
            "skani":                  skani_data,
            "gambit":                 gambit_data,
            "kraken2_db":             args.kraken2_db,
        },
    }

    for parser, key, fname in [
        (parse_kleborate, "kleborate", "Kleborate/kleborate_output.tsv"),
        (parse_ectyper,   "ectyper",   "ECTyper/ectyper.tsv"),
        (parse_emmtyper,  "emmtyper",  "EmmTyper/emmtyper.txt"),
        (parse_seqsero2,  "seqsero2",  "SeqSero2/seqsero2.tsv"),
        (parse_hicap,     "hicap",     "Hicap/hicap.tsv"),
        (parse_pasty,        "pasty",        "Pasty/pasty.tsv"),
        (parse_lrefinder,    "lrefinder",    "LRE-Finder/lre.res"),
        (parse_gbsserotyper, "gbsserotyper", "GBSSeroTyper/serotype.tsv"),
    ]:
        result = parser(sp_dir / fname)
        if result:
            data[key] = result

    ss = parse_staphscope(sp_dir / "Staphscope/Staphscope_final_report/staphscope_comprehensive_report.tsv")
    if ss:
        data["staphscope"] = ss
        data["spatyper"]   = {"spa_type": ss["spa_type"], "CC": "-"}
        data["sccmec"]     = {"type": ss["sccmec_type"], "mrsa": ss["mrsa"]}

    if varcall_data:
        data["varcall"] = varcall_data

    Path(args.output).write_text(
        Path(args.template).read_text()
            .replace("__PIPELINE_JSON_DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__SAMPLE_NAME__", args.sample)
    )

    # --- Cumulative run log ---
    # Static, single TSV — re-running a sample upserts (replaces) its row.
    log_path = Path(args.log_path) if args.log_path else \
               Path(args.results_dir) / "run_log.tsv"
    header = "\t".join(["#", "Sample", "Species (GAMBIT)", "GAMBIT hit", "Skani ANI",
                         "Bracken top", "ST", "MLST tool", "Species-specific typing",
                         "Depth (mapped)", "N50", "Q30", "Date"])
    skani  = data["species_id"]["skani"]
    gambit = data["species_id"]["gambit"]

    typing_parts = []
    if data.get("kleborate"):
        k     = data["kleborate"]
        parts = ["Kleborate:"]
        if k.get("ST")       and k["ST"]       != "-": parts.append(f"ST{k['ST']}")
        if k.get("K_type")   and k["K_type"]   != "-": parts.append(f"K={k['K_type']}")
        if k.get("O_type")   and k["O_type"]   != "-": parts.append(f"O={k['O_type']}")
        if k.get("Bla_Carb") and k["Bla_Carb"] != "-": parts.append(f"Carb={k['Bla_Carb']}")
        if k.get("Bla_ESBL") and k["Bla_ESBL"] != "-": parts.append(f"ESBL={k['Bla_ESBL']}")
        if len(parts) > 1:
            typing_parts.append(" ".join(parts))
    if data.get("ectyper"):
        typing_parts.append(f"ECTyper: {data['ectyper'].get('serotype', '-')}")
    if data.get("staphscope"):
        ss = data["staphscope"]
        typing_parts.append(f"StaphScope: spa={ss.get('spa_type','-')} SCCmec={ss.get('sccmec_type','-')} "
                            f"({'MRSA' if ss.get('mrsa') else 'MSSA'})")
    if data.get("emmtyper"):
        typing_parts.append(f"emmtyper: {data['emmtyper'].get('emm_type', '-')}")
    if data.get("seqsero2"):
        typing_parts.append(f"SeqSero2: {data['seqsero2'].get('serotype', '-')}")
    if data.get("hicap"):
        typing_parts.append(f"hicap: {data['hicap'].get('serotype', '-')}")
    if data.get("pasty"):
        typing_parts.append(f"Pasty: {data['pasty'].get('serotype', '-')}")
    if data.get("gbsserotyper"):
        typing_parts.append(f"GBS-SBG: {data['gbsserotyper'].get('serotype', '-')}")

    row = "\t".join([
        "{nr}",
        args.sample,
        species_sp.replace("_", " ") if species_sp not in ("unknown", "") else "-",
        gambit.get("closest_description", "-") if gambit.get("status") == "ok" else "-",
        f"{skani['ani']:.2f}" if skani.get("status") == "ok" and skani.get("ani") else "-",
        bracken_data["primary_species"].replace("_", " "),
        mlst_data["ST"],
        f"mlst ({mlst_data['scheme']})" if mlst_data["scheme"] != "-" else "mlst",
        "; ".join(typing_parts) if typing_parts else "-",
        f"{self_cov['mean_depth']:.1f}" if self_cov.get("mean_depth") else "-",
        str(quast_data.get("N50") or "-"),
        f"{fastp_data['q30_pct']:.1f}" if fastp_data.get("q30_pct") is not None else "-",
        datetime.now().strftime("%Y-%m-%d %H:%M"),
    ])

    # Parallel report jobs all write to the same TSV; flock serialises the
    # read-modify-write so rows don't clobber each other. Re-runs of the same
    # sample overwrite the prior row (matched by column 1 = sample ID).
    log_path.touch()
    with open(log_path, "r+") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            existing = f.read().splitlines()
            body = existing[1:] if existing and existing[0].startswith("#") else existing
            # Drop any prior row for this sample, then append the fresh one.
            data = [l for l in body
                    if l.strip() and (l.split("\t") + [""])[1] != args.sample]
            data.append(row)
            # Renumber column 0 (#) so the file stays 1..N after deletions.
            renumbered = []
            for i, l in enumerate(data, 1):
                parts    = l.split("\t")
                parts[0] = str(i)
                renumbered.append("\t".join(parts))
            f.seek(0)
            f.truncate()
            f.write("\n".join([header] + renumbered) + "\n")
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
