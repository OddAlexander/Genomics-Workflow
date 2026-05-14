#!/usr/bin/env python3
"""Generate HTML report for a single sample from pipeline outputs."""
import argparse, csv, json
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


def parse_gambit_full(path):
    txt = safe_read(path)
    if not txt:
        return {}
    rows = list(csv.DictReader(txt.splitlines()))
    if not rows:
        return {}
    r = rows[0]
    def _f(key):
        try:
            return float(r.get(key, "") or "")
        except ValueError:
            return None
    return {
        "predicted_name":      r.get("predicted.name",      "-") or "-",
        "predicted_rank":      r.get("predicted.rank",      "-") or "-",
        "predicted_threshold": _f("predicted.threshold"),
        "closest_distance":    _f("closest.distance"),
        "closest_description": r.get("closest.description", "-") or "-",
        "next_name":           r.get("next.name",           "-") or "-",
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


def derive_purity(primary_pct, gambit_sp, bracken_sp):
    if primary_pct >= 80 and gambit_sp == bracken_sp:
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


def concordance_score(gambit_sp, bracken_sp):
    return 2 if gambit_sp == bracken_sp and gambit_sp not in ("unknown", "-", "") else 0


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
    fastp_data   = parse_fastp(sp_dir / "QC/fastp.json")
    est_depth    = estimate_depth(fastp_data["bases"], quast_data["total_length"])
    bracken_data = parse_bracken(sp_dir / "ID_Kraken2/bracken_species.txt")
    uncl_pct     = parse_kraken2_unclassified(sp_dir / "ID_Kraken2/kraken2_report.txt")
    amr_raw      = parse_amrfinder(sp_dir / "AMRFinder/amrfinder.tsv")
    mob_data     = parse_mobsuite(sp_dir / "MOBSuite/mobtyper.tsv")
    mge_data     = parse_mefinder(sp_dir / "MEfinder/mefinder.tsv")
    pf_data      = parse_plasmidfinder(sp_dir / "PlasmidFinder/results_tab.tsv")
    skani_data   = parse_skani(sp_dir / "ID_Skani/skani.tsv")
    gambit_full  = parse_gambit_full(sp_dir / "ID_GAMBIT/gambit.csv")
    gambit_sp    = (safe_read(sp_dir / "species.txt") or "unknown").strip()

    contig_mge = {}
    for el in mge_data.get("elements", []):
        c = el.get("contig", "-")
        if c and c != "-":
            contig_mge.setdefault(c, []).append({"name": el["name"], "type": el["type"]})

    amr_details, high_risk, medium_risk = build_amr_details(
        amr_raw.get("_rows", []), mob_data.get("_contig_mob", {}),
        contig_mge, pf_data.get("_contig_plasmid", {}))

    purity  = derive_purity(bracken_data["primary_pct"], gambit_sp, bracken_data["primary_species"])
    score   = concordance_score(gambit_sp, bracken_data["primary_species"])
    species = gambit_sp if gambit_sp not in ("unknown", "") else bracken_data["primary_species"]
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
            "reads":     fastp_data["reads"],
            "q30_pct":   fastp_data["q30_pct"],
            "dup_pct":   fastp_data["dup_pct"],
            "est_depth": est_depth,
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
            "purity":                     purity,
            "concordance_score":          score,
            "contamination_warnings":     contam_warnings,
            "bracken_species":            bracken_data["primary_species"],
            "gambit_species":             gambit_sp,
            "gambit_predicted_name":      gambit_full.get("predicted_name",      "-"),
            "gambit_predicted_rank":      gambit_full.get("predicted_rank",      "-"),
            "gambit_predicted_threshold": gambit_full.get("predicted_threshold"),
            "gambit_distance":            gambit_full.get("closest_distance"),
            "gambit_closest_description": gambit_full.get("closest_description", "-"),
            "gambit_next_name":           gambit_full.get("next_name",           "-"),
            "id_method":                  "Bracken + GAMBIT",
            "primary_species":            bracken_data["primary_species"],
            "primary_pct":                bracken_data["primary_pct"],
            "secondary_species":          bracken_data["secondary_species"],
            "secondary_pct":              bracken_data["secondary_pct"],
            "unclassified_pct":           uncl_pct,
            "skani":                      skani_data,
            "kraken2_db":                 args.kraken2_db,
        },
    }

    for parser, key, fname in [
        (parse_kleborate, "kleborate", "Kleborate/kleborate_output.tsv"),
        (parse_ectyper,   "ectyper",   "ECTyper/ectyper.tsv"),
        (parse_emmtyper,  "emmtyper",  "EmmTyper/emmtyper.txt"),
        (parse_seqsero2,  "seqsero2",  "SeqSero2/seqsero2.tsv"),
        (parse_hicap,     "hicap",     "Hicap/hicap.tsv"),
    ]:
        result = parser(sp_dir / fname)
        if result:
            data[key] = result

    ss = parse_staphscope(sp_dir / "Staphscope/Staphscope_final_report/staphscope_comprehensive_report.tsv")
    if ss:
        data["staphscope"] = ss
        data["spatyper"]   = {"spa_type": ss["spa_type"], "CC": "-"}
        data["sccmec"]     = {"type": ss["sccmec_type"], "mrsa": ss["mrsa"]}

    Path(args.output).write_text(
        Path(args.template).read_text()
            .replace("__PIPELINE_JSON_DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__SAMPLE_NAME__", args.sample)
    )

    # --- Cumulative run log ---
    log_path = Path(args.log_path) if args.log_path else \
               Path(args.results_dir) / f"pipeline_run_{datetime.now().strftime('%Y%m%d_%H%M')}.tsv"
    header = "\t".join(["#", "Sample", "GAMBIT", "Bracken", "Skani hit", "Skani ANI",
                         "ST", "MLST tool", "Species-specific typing", "Date"])
    skani  = data["species_id"]["skani"]

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

    pasty_txt = safe_read(sp_dir / "Pasty/pasty.tsv")
    if pasty_txt:
        pasty_rows = list(csv.DictReader(pasty_txt.splitlines(), delimiter="\t"))
        if pasty_rows and pasty_rows[0].get("type") and pasty_rows[0]["type"] != "-":
            typing_parts.append(f"Pasty: {pasty_rows[0]['type']}")

    row = "\t".join([
        "{nr}",
        args.sample,
        gambit_full.get("predicted_name", "-"),
        bracken_data["primary_species"].replace("_", " "),
        skani.get("top_hit", "-") if skani.get("status") == "ok" else "-",
        f"{skani['ani']:.2f}" if skani.get("status") == "ok" and skani.get("ani") else "-",
        mlst_data["ST"],
        f"mlst ({mlst_data['scheme']})" if mlst_data["scheme"] != "-" else "mlst",
        "; ".join(typing_parts) if typing_parts else "-",
        datetime.now().strftime("%Y-%m-%d %H:%M"),
    ])

    if log_path.exists():
        lines = log_path.read_text().splitlines()
        nr    = len([l for l in lines if l and not l.startswith("#")])
    else:
        lines = [header]
        nr    = 0
    lines.append(row.format(nr=nr + 1))
    log_path.write_text("\n".join(lines) + "\n")
