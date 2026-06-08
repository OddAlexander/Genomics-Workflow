#!/usr/bin/env bash
# =============================================================================
# pixi_link-amr-plasmid.sh
# Links AMRFinder resistance genes to MOB-suite plasmid typing results
# by matching contig names
#
# Usage:
#   ./pixi_link-amr-plasmid.sh <sample_dir>
#
# Requirements:
#   - AMR/amrfinder_amr.tsv        (from run-amr.sh)
#   - Plasmids/contig_report.txt   (from run-mobsuite.sh)
#   - Plasmids/plasmid_summary.tsv (from run-mobsuite.sh)
#
# Output:
#   <sample_dir>/AMR/amr_plasmid_links.tsv   Combined report
# =============================================================================

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()     { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
success() { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${RESET}"; }
warn()    { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠ $*${RESET}"; }
die()     { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ $*${RESET}" >&2; exit 1; }
header()  { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${RESET}"
            echo -e "${BOLD}${CYAN}  $*${RESET}"
            echo -e "${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; }

if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "     A M R   —   P L A S M I D   L I N K E R               "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir>"
    echo ""
    echo -e "${BOLD}Requirements:${RESET}"
    echo "  Run run-amr.sh and run-mobsuite.sh first"
    echo ""
    echo -e "${BOLD}Output:${RESET}"
    echo "  <sample_dir>/AMR/amr_plasmid_links.tsv"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  $0 results/001k/"
    echo "============================================================"
    exit 0
fi

SAMPLE_DIR="${1%/}"
[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

SAMPLE=$(basename "$SAMPLE_DIR")

AMR_TSV="$SAMPLE_DIR/AMR/amrfinder_amr.tsv"
CONTIG_REPORT="$SAMPLE_DIR/Plasmids/contig_report.txt"
PLASMID_SUMMARY="$SAMPLE_DIR/Plasmids/plasmid_summary.tsv"
OUTPUT_TSV="$SAMPLE_DIR/AMR/amr_plasmid_links.tsv"

[[ -f "$AMR_TSV" ]]         || die "AMR results not found: $AMR_TSV\nRun run-amr.sh first"
[[ -f "$CONTIG_REPORT" ]]   || die "MOB-suite contig report not found: $CONTIG_REPORT\nRun run-mobsuite.sh first"
[[ -f "$PLASMID_SUMMARY" ]] || die "MOB-suite plasmid summary not found: $PLASMID_SUMMARY\nRun run-mobsuite.sh first"

log "Sample:  $SAMPLE"
log "AMR:     $AMR_TSV"
log "Plasmid: $PLASMID_SUMMARY"

header "Linking AMR genes to plasmid typing results"

python3 - "$AMR_TSV" "$CONTIG_REPORT" "$PLASMID_SUMMARY" "$OUTPUT_TSV" "$SAMPLE" << 'PYEOF'
import sys, csv

amr_file      = sys.argv[1]
contig_file   = sys.argv[2]
plasmid_file  = sys.argv[3]
out_file      = sys.argv[4]
sample        = sys.argv[5]

RED    = '\033[0;31m'
YELLOW = '\033[1;33m'
GREEN  = '\033[0;32m'
BOLD   = '\033[1m'
RESET  = '\033[0m'

# ── Load MOB-suite contig report ───────────────────────────────────────────────
# contig_report.txt columns vary by version but always have contig and molecule_type
# Typical: file_id, contig_id, molecule_type, ...
contig_to_molecule = {}  # contig_id → chromosome/plasmid/unclassified
contig_to_file     = {}  # contig_id → plasmid file name (e.g. plasmid_1.fasta)

with open(contig_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # Handle different column names across MOB-suite versions
        contig = (row.get('contig_id') or row.get('contig') or
                  row.get('sequence_id') or '').strip()
        mol    = (row.get('molecule_type') or row.get('type') or '').strip().lower()
        fid    = (row.get('file_id') or row.get('primary_cluster_id') or '').strip()
        if contig:
            contig_to_molecule[contig] = mol
            contig_to_file[contig]     = fid

# ── Load MOB-suite plasmid summary ────────────────────────────────────────────
# keyed by file_id
plasmid_info = {}  # file_id → dict of typing info

with open(plasmid_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        fid = row.get('file_id', '').strip()
        if fid:
            plasmid_info[fid] = row

# ── Load AMRFinder results ─────────────────────────────────────────────────────
amr_rows = []
with open(amr_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        amr_rows.append(row)

if not amr_rows:
    print(f"  {YELLOW}No AMR genes found in {amr_file}{RESET}")
    sys.exit(0)

# ── Build linked results ───────────────────────────────────────────────────────
# AMRFinder v4 column names (exact)
linked = []
for row in amr_rows:
    contig       = row.get('Contig id', '').strip()
    gene         = row.get('Element symbol', '').strip()
    gene_name    = row.get('Element name', '').strip()
    el_type      = row.get('Type', '').strip()
    subtype      = row.get('Subtype', '').strip()
    el_class     = row.get('Class', '').strip()
    el_subclass  = row.get('Subclass', '').strip()
    pct_cov      = row.get('% Coverage of reference', '').strip()
    pct_id       = row.get('% Identity to reference', '').strip()
    scope        = row.get('Scope', '').strip()

    mol  = contig_to_molecule.get(contig, 'unknown')
    fid  = contig_to_file.get(contig, '-')

    # Get plasmid details if this contig is on a plasmid
    pinfo     = plasmid_info.get(fid, {})
    rep_type  = pinfo.get('rep_type(s)', '-') or '-'
    relaxase  = pinfo.get('relaxase_type(s)', '-') or '-'
    mobility  = pinfo.get('predicted_mobility', '-') or '-'
    cluster   = pinfo.get('primary_cluster_id', '-') or '-'
    host_range= pinfo.get('predicted_host_range_overall_rank', '-') or '-'
    p_size    = pinfo.get('size', '-') or '-'

    linked.append({
        'sample':      sample,
        'gene':        gene,
        'gene_name':   gene_name,
        'type':        el_type,
        'subtype':     subtype,
        'class':       el_class,
        'subclass':    el_subclass,
        'pct_coverage':pct_cov,
        'pct_identity':pct_id,
        'scope':       scope,
        'contig':      contig,
        'location':    mol,
        'plasmid_file':fid,
        'replicon':    rep_type,
        'relaxase':    relaxase,
        'mobility':    mobility,
        'host_range':  host_range,
        'cluster':     cluster,
        'plasmid_size':p_size,
    })

# ── Write TSV ──────────────────────────────────────────────────────────────────
fieldnames = ['sample','gene','gene_name','type','subtype','class','subclass',
              'pct_coverage','pct_identity','scope','contig','location',
              'plasmid_file','replicon','relaxase','mobility',
              'host_range','cluster','plasmid_size']

with open(out_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(linked)

# ── Print summary ──────────────────────────────────────────────────────────────
plasmid_genes = []
chrom_genes   = []

for r in linked:
    mob = r['mobility']
    loc = r['location']

    if 'plasmid' in loc.lower():
        plasmid_genes.append(r)
    elif 'chromosome' in loc.lower():
        chrom_genes.append(r)

# Print plasmid-borne genes with full details
if plasmid_genes:
    print(f"\n  {BOLD}Plasmid-borne resistance/virulence genes:{RESET}")
    print(f"  {'─'*90}")
    print(f"  {'Gene':<16} {'Name':<30} {'Class':<14} {'Mobility':<16} {'Replicon':<12} {'Cluster'}")
    print(f"  {'─'*90}")
    for r in plasmid_genes:
        mob = r['mobility']
        if mob == 'CONJUGATIVE':
            mob_str = f"{RED}{mob:<16}{RESET}"
        elif mob == 'MOBILIZABLE':
            mob_str = f"{YELLOW}{mob:<16}{RESET}"
        else:
            mob_str = f"{GREEN}{mob:<16}{RESET}"
        name = r['gene_name'][:28] if len(r['gene_name']) > 28 else r['gene_name']
        print(f"  {r['gene']:<16} {name:<30} {r['class']:<14} {mob_str} {r['replicon']:<12} {r['cluster']}")
    print(f"  {'─'*90}")

if chrom_genes:
    print(f"\n  {BOLD}Chromosomal resistance/virulence genes:{RESET}")
    print(f"  {'─'*70}")
    print(f"  {'Gene':<16} {'Name':<30} {'Class':<14} {'Subclass'}")
    print(f"  {'─'*70}")
    for r in chrom_genes:
        name = r['gene_name'][:28] if len(r['gene_name']) > 28 else r['gene_name']
        print(f"  {r['gene']:<16} {name:<30} {r['class']:<14} {r['subclass']}")
    print(f"  {'─'*70}")

print(f"\n  Plasmid-borne genes:  {len(plasmid_genes)}")
print(f"  Chromosomal genes:    {len(chrom_genes)}")
print(f"  Unknown location:     {len(linked) - len(plasmid_genes) - len(chrom_genes)}")

# Clinical risk highlight
conj_plasmid_genes = [r for r in plasmid_genes if r['mobility'] == 'CONJUGATIVE']
if conj_plasmid_genes:
    print(f"\n  {RED}{BOLD}⚠ HIGH RISK — genes on CONJUGATIVE plasmids:{RESET}")
    for r in conj_plasmid_genes:
        print(f"    {RED}• {r['gene']} ({r['gene_name']}) — {r['class']} — {r['replicon']} — cluster {r['cluster']}{RESET}")

mob_plasmid_genes = [r for r in plasmid_genes if r['mobility'] == 'MOBILIZABLE']
if mob_plasmid_genes:
    print(f"\n  {YELLOW}⚠ MODERATE RISK — genes on MOBILIZABLE plasmids:{RESET}")
    for r in mob_plasmid_genes:
        print(f"    {YELLOW}• {r['gene']} ({r['gene_name']}) — {r['class']} — {r['replicon']} — cluster {r['cluster']}{RESET}")
PYEOF

success "Linking done"
echo ""
echo "────────────────────────────────────────────────────────"
echo "  Output: $OUTPUT_TSV"
echo "────────────────────────────────────────────────────────"
