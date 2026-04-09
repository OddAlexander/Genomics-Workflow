#!/usr/bin/env python3
"""Generate HTML report for a single sample from pipeline outputs."""
import argparse, csv, json, os
from datetime import datetime
from pathlib import Path


def safe_read(path):
    p = Path(path)

    return p.read_text().strip() if p.exists() and p.stat().st_size > 0 else None


def parse_mlst(path):
    txt = safe_read(path)
    if not txt:
        return {"ST": "-", "scheme": "-"}
    parts = txt.split("\t")
    return {"ST": parts[2] if len(parts) > 2 else "-",
            "scheme": parts[1] if len(parts) > 1 else "-"}


def parse_quast(report_html):
    tsv = Path(report_html).parent / "report.tsv"
    if not tsv.exists():
        return {"contigs": "-", "total_length": 0, "N50": 0, "largest_contig": 0}
    data = {}
    for line in tsv.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            data[parts[0].strip()] = parts[1].strip().replace(",", "")
    return {
        "contigs":        data.get("# contigs", "-"),
        "total_length":   int(data.get("Total length", 0) or 0),
        "N50":            int(data.get("N50", 0) or 0),
        "largest_contig": int(data.get("Largest contig", 0) or 0),
    }


def parse_bracken(path):
    # kreport2: %  cov_reads  taxon_reads  rank  taxID  name
    txt = safe_read(path)
    result = {"primary_species": "-", "primary_pct": 0.0,
              "secondary_species": "-", "secondary_pct": 0.0}
    if not txt:
        return result
    rows = []
    for line in txt.splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        try:
            rows.append((float(parts[0].strip()), parts[5].strip()))
        except ValueError:
            continue
    rows.sort(reverse=True)
    if rows:
        result["primary_pct"]     = rows[0][0]
        result["primary_species"] = rows[0][1].replace(" ", "_")
    if len(rows) > 1:
        result["secondary_pct"]     = rows[1][0]
        result["secondary_species"] = rows[1][1].replace(" ", "_")
    return result


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
        return {"total_amr_genes": 0, "total_virulence_genes": 0, "_rows": []}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    return {
        "total_amr_genes":      sum(1 for r in rows if "AMR" in r.get("Element type", "")),
        "total_virulence_genes": sum(1 for r in rows if "VIRULENCE" in r.get("Element type", "")),
        "_rows": rows,
    }


def parse_mobsuite(path):
    txt = safe_read(path)
    if not txt:
        return {"count": 0, "conjugative": 0, "mobilizable": 0,
                "non_mobilizable": 0, "details": [], "_contig_mob": {}}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    contig_mob = {
        r.get("file_id", ""): {
            "mobility": r.get("predicted_mobility", "-").upper(),
            "replicon": r.get("rep_type(s)", "-"),
            "cluster":  r.get("cluster_id", "-"),
        }
        for r in rows
    }
    return {
        "count":           len(rows),
        "conjugative":     sum(1 for r in rows if "CONJUGATIVE"    in r.get("predicted_mobility", "").upper()),
        "mobilizable":     sum(1 for r in rows if r.get("predicted_mobility", "").upper() == "MOBILIZABLE"),
        "non_mobilizable": sum(1 for r in rows if "NON_MOBILIZABLE" in r.get("predicted_mobility", "").upper()),
        "details": [{
            "file_id":    r.get("file_id", "-"),
            "replicon":   r.get("rep_type(s)", "-"),
            "relaxase":   r.get("relaxase_type(s)", "-"),
            "mobility":   r.get("predicted_mobility", "-").upper(),
            "host_range": r.get("mash_nearest_neighbor", "-"),
            "cluster":    r.get("cluster_id", "-"),
            "size":       r.get("size", 0),
        } for r in rows],
        "_contig_mob": contig_mob,
    }


def build_amr_details(amr_rows, contig_mob):
    details, high_risk, medium_risk = [], [], []
    for g in amr_rows:
        contig   = g.get("Contig id", g.get("Protein identifier", "-"))
        mob_info = contig_mob.get(contig, {})
        mobility = mob_info.get("mobility", "-")
        gene     = g.get("Gene symbol", g.get("Element symbol", "-"))
        if "AMR" in g.get("Element type", ""):
            if mobility == "CONJUGATIVE":
                high_risk.append(gene)
            elif mobility == "MOBILIZABLE":
                medium_risk.append(gene)
        details.append({
            "gene":      gene,
            "gene_name": g.get("Sequence name", g.get("Element name", "-")),
            "class":     g.get("Class", "-"),
            "location":  contig,
            "mobility":  mobility,
            "replicon":  mob_info.get("replicon", "-"),
            "cluster":   mob_info.get("cluster", "-"),
        })
    return details, high_risk, medium_risk


def parse_sccmec(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return None
    r = rows[0]
    t = r.get("type", r.get("Type", "-")) or "-"
    return {"type": t, "mrsa": t not in ("-", "none", "")}


def parse_spatyper(path):
    txt = safe_read(path)
    if not txt:
        return None
    for line in txt.splitlines():
        if line.startswith("Type") or not line.strip():
            continue
        parts = line.strip().split("\t")
        return {"spa_type": parts[0] or "-", "CC": parts[1] if len(parts) > 1 else "-"}
    return {"spa_type": "-", "CC": "-"}


def parse_kleborate(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return None
    r = rows[0]
    return {
        "K_type":          r.get("K_locus", r.get("K_type", "-")),
        "O_type":          r.get("O_locus", r.get("O_type", "-")),
        "virulence_score": r.get("virulence_score", "-"),
        "Bla_Carb":        r.get("Bla_Carb_acquired", "-"),
        "Bla_ESBL":        r.get("Bla_ESBL_acquired", "-"),
    }


def parse_ectyper(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return None
    r = rows[0]
    return {"serotype": f"{r.get('O-type','-')}:{r.get('H-type','-')}"}


def parse_emmtyper(path):
    txt = safe_read(path)
    if not txt:
        return None
    for line in txt.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            return {"emm_type": parts[1] or "-"}
    return None


def parse_fastani(path):
    txt = safe_read(path)
    if not txt:
        return {"status": "hoppet_over", "top_hit": "-", "ani": None, "af": None}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows or rows[0].get("status") == "hoppet_over":
        return {"status": "hoppet_over", "top_hit": "-", "ani": None, "af": None}
    # FastANI output: query, reference, ANI, bidirectional_frags, total_frags
    r = rows[0]
    ref = r.get("reference", list(r.values())[1] if len(r) > 1 else "-")
    ani = r.get("ANI", list(r.values())[2] if len(r) > 2 else None)
    return {
        "status":  "ok",
        "top_hit": ref.split("/")[-1].replace(".fna", "").replace(".fa", ""),
        "ani":     float(ani) if ani else None,
        "af":      r.get("bidirectional_fragment_count", "-"),
    }


def parse_gambit_distance(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines()))
    if not rows:
        return None
    try:
        return float(rows[0].get("closest.distance", "") or "")
    except ValueError:
        return None


def parse_seqsero2(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return None
    r = rows[0]
    return {
        "serotype":        r.get("Predicted serotype", "-") or "-",
        "antigenic_profile": r.get("Predicted antigenic profile", "-") or "-",
        "O_antigen":       r.get("O antigen prediction", "-") or "-",
        "comment":         r.get("Note", "") or "",
    }


def parse_hicap(path):
    txt = safe_read(path)
    if not txt:
        return None
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if not rows:
        return None
    r = rows[0]
    return {
        "serotype":     r.get("predicted_serotype", "-") or "-",
        "genes":        r.get("genes_identified", "-") or "-",
        "completeness": r.get("percent_coverage", r.get("completeness", "-")) or "-",
    }


def parse_mefinder(path):
    txt = safe_read(path)
    if not txt:
        return {"count": 0, "elements": []}
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    elements = [{
        "nr":     i,
        "name":   r.get("name", r.get("Name", "-")),
        "type":   r.get("type", r.get("Type", "-")),
        "contig": r.get("contig_id", r.get("sequence_id", "-")),
        "start":  r.get("start", r.get("Start", "-")),
        "end":    r.get("end",   r.get("End",   "-")),
    } for i, r in enumerate(rows, 1)]
    return {"count": len(elements), "elements": elements}


def derive_purity(primary_pct, gambit_sp, bracken_sp):
    concordant = (gambit_sp == bracken_sp)
    if primary_pct >= 80 and concordant:
        return "RENKULTUR"
    elif primary_pct >= 60:
        return "MULIG RENKULTUR"
    else:
        return "BLANDING"


def concordance_score(gambit_sp, bracken_sp):
    if gambit_sp == bracken_sp and gambit_sp not in ("unknown", "-", ""):
        return 2  # begge enige
    return 0  # ingen samsvar


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample",      required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--template",    required=True)
    ap.add_argument("--output",      required=True)
    ap.add_argument("--kraken2-db",  default="-")
    args = ap.parse_args()

    rd     = Path(args.results_dir)
    sp_dir = rd / args.sample

    mlst_data    = parse_mlst(sp_dir / "MLST/mlst.tsv")
    quast_data   = parse_quast(sp_dir / "QUAST/report.html")
    bracken_data = parse_bracken(sp_dir / "ID_Kraken2/bracken_species.txt")
    uncl_pct     = parse_kraken2_unclassified(sp_dir / "ID_Kraken2/kraken2_report.txt")
    amr_raw      = parse_amrfinder(sp_dir / "AMRFinder/amrfinder.tsv")
    mob_data     = parse_mobsuite(sp_dir / "MOBSuite/mobtyper.tsv")
    mge_data     = parse_mefinder(sp_dir / "MEfinder/mefinder.tsv")
    fastani_data = parse_fastani(sp_dir / "ID_FastANI/fastani.tsv")
    gambit_dist  = parse_gambit_distance(sp_dir / "ID_GAMBIT/gambit.csv")

    gambit_sp = (safe_read(sp_dir / "species.txt") or "unknown").strip()

    amr_details, high_risk, medium_risk = build_amr_details(
        amr_raw.get("_rows", []), mob_data.get("_contig_mob", {}))

    purity  = derive_purity(bracken_data["primary_pct"], gambit_sp, bracken_data["primary_species"])
    score   = concordance_score(gambit_sp, bracken_data["primary_species"])
    species = gambit_sp if gambit_sp not in ("unknown", "") else bracken_data["primary_species"]

    data = {
        "sample":   args.sample,
        "species":  species,
        "completed": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "mlst":     mlst_data,
        "assembly": quast_data,
        "amr": {
            "total_amr_genes":       amr_raw["total_amr_genes"],
            "total_virulence_genes": amr_raw["total_virulence_genes"],
            "high_risk_genes":       high_risk,
            "medium_risk_genes":     medium_risk,
            "details":               amr_details,
        },
        "plasmids": {k: v for k, v in mob_data.items() if not k.startswith("_")},
        "mge":      mge_data,
        "species_id": {
            "purity":             purity,
            "concordance_score":  score,
            "bracken_species":    bracken_data["primary_species"],
            "gambit_species":     gambit_sp,
            "gambit_distance":    gambit_dist,
            "id_method":          "Bracken + GAMBIT",
            "primary_species":    bracken_data["primary_species"],
            "primary_pct":        bracken_data["primary_pct"],
            "secondary_species":  bracken_data["secondary_species"],
            "secondary_pct":      bracken_data["secondary_pct"],
            "unclassified_pct":   uncl_pct,
            "fastani":            fastani_data,
            "kraken2_db":         args.kraken2_db,
        },
    }

    # Species-specific — only add if file exists
    for parser, key, fname in [
        (parse_sccmec,   "sccmec",    "SCCmec/sccmec.txt"),
        (parse_spatyper, "spatyper",  "SpaTyper/spatyper.txt"),
        (parse_kleborate,"kleborate", "Kleborate/kleborate_output.tsv"),
        (parse_ectyper,  "ectyper",   "ECTyper/ectyper.tsv"),
        (parse_emmtyper, "emmtyper",  "EmmTyper/emmtyper.txt"),
        (parse_seqsero2, "seqsero2",  "SeqSero2/seqsero2.tsv"),
        (parse_hicap,    "hicap",     "Hicap/hicap.tsv"),
    ]:
        result = parser(sp_dir / fname)
        if result:
            data[key] = result

    template = Path(args.template).read_text()
    html = (template
            .replace("__PIPELINE_JSON_DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__SAMPLE_NAME__", args.sample))
    Path(args.output).write_text(html)
