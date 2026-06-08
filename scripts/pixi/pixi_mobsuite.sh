#!/usr/bin/env bash
# =============================================================================
# run-mobsuite.sh
# Runs MOB-suite plasmid typing on an assembly
#
# Usage:
#   ./run-mobsuite.sh <sample_dir>
#
# Tools run:
#   mob_typer  — replicon, relaxase, conjugation prediction, host-range
#   mob_recon  — reconstruct individual plasmid sequences from draft assembly
#
# Output:
#   <sample_dir>/Plasmids/
#     plasmid_summary.tsv       Clean summary of key columns
#     mobtyper_results.txt      Full MOB-typer raw output
#     chromosome.fasta          Predicted chromosomal contigs
#     plasmid_*.fasta           Reconstructed plasmid sequences
#     contig_report.txt         Per-contig classification
#
# Optional env vars:
#   THREADS — number of threads (default: 8)
#
# First run: databases are auto-downloaded (~1.5 GB) via mob_init
# Manual init: pixi run --environment mobsuite -- mob_init
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

# ── No arguments — show help ──────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "         M O B - S U I T E   P L A S M I D   T Y P I N G  "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir>"
    echo ""
    echo -e "${BOLD}What it does:${RESET}"
    echo "  mob_typer  — predicts replicon, relaxase type, conjugation"
    echo "               potential and host-range for each plasmid"
    echo "  mob_recon  — reconstructs individual plasmid sequences from"
    echo "               the draft assembly contigs"
    echo ""
    echo -e "${BOLD}Mobility classification:${RESET}"
    echo "  CONJUGATIVE     — has relaxase + MPF genes, can self-transfer"
    echo "  MOBILIZABLE     — has relaxase, needs helper plasmid to transfer"
    echo "  NON-MOBILIZABLE — cannot transfer by conjugation"
    echo ""
    echo -e "${BOLD}Output:${RESET}"
    echo "  <sample_dir>/Plasmids/"
    echo "    plasmid_summary.tsv    Clean summary (key columns only)"
    echo "    mobtyper_results.txt   Full raw MOB-typer output"
    echo "    contig_report.txt      Per-contig chromosome vs plasmid"
    echo "    chromosome.fasta       Chromosomal contigs"
    echo "    plasmid_*.fasta        Reconstructed plasmid sequences"
    echo ""
    echo -e "${BOLD}First run:${RESET}"
    echo "  Databases auto-download (~1.5 GB) on first use"
    echo "  Or initialise manually:"
    echo "  pixi run --environment mobsuite -- mob_init"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  $0 results/001k/"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
THREADS="${THREADS:-8}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

# ── Find assembly ─────────────────────────────────────────────────────────────
ASSEMBLY_DIR="$SAMPLE_DIR/Assembly"
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly/ not found in $SAMPLE_DIR\nRun bacterial-pipeline.sh first"

if [[ -f "$ASSEMBLY_DIR/scaffolds.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/scaffolds.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fa" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fa"
else
    die "No assembly file found in $ASSEMBLY_DIR"
fi

SAMPLE=$(basename "$SAMPLE_DIR")
log "Sample:   $SAMPLE"
log "Assembly: $(basename "$ASSEMBLY_FA")"

OUTPUT_DIR="$SAMPLE_DIR/Plasmids"
mkdir -p "$OUTPUT_DIR"
log "Output:   $OUTPUT_DIR"

STEP_START=0
start_timer() { STEP_START=$SECONDS; }
elapsed()     { echo $(( SECONDS - STEP_START ))s; }

MOBTYPER_OUT="$OUTPUT_DIR/mobtyper_results.txt"
SUMMARY_TSV="$OUTPUT_DIR/plasmid_summary.tsv"

# =============================================================================
header "Step 1/2 — MOB-recon: reconstruct plasmid sequences"
# =============================================================================
start_timer

pixi run --environment mobsuite -- mob_recon \
    --infile "$ASSEMBLY_FA" \
    --outdir "$OUTPUT_DIR" \
    --num_threads "$THREADS" \
    --force

success "MOB-recon done ($(elapsed))"

PLASMID_COUNT=$(find "$OUTPUT_DIR" -name "plasmid_*.fasta" 2>/dev/null | wc -l)
log "Plasmid sequences reconstructed: $PLASMID_COUNT"

# =============================================================================
header "Step 2/2 — MOB-typer: replicon, relaxase, conjugation typing"
# =============================================================================
start_timer

pixi run --environment mobsuite -- mob_typer \
    --infile "$ASSEMBLY_FA" \
    --out_file "$MOBTYPER_OUT" \
    --num_threads "$THREADS"

success "MOB-typer done ($(elapsed))"

# =============================================================================
header "Results"
# =============================================================================
echo ""

# ── Contig classification summary ─────────────────────────────────────────────
if [[ -f "$OUTPUT_DIR/contig_report.txt" ]]; then
    CHROM_COUNT=$(grep -c "chromosome" "$OUTPUT_DIR/contig_report.txt" 2>/dev/null || echo 0)
    PLASMID_CONTIGS=$(grep -c "plasmid" "$OUTPUT_DIR/contig_report.txt" 2>/dev/null || echo 0)
    UNCLASSIFIED=$(grep -c "unclassified" "$OUTPUT_DIR/contig_report.txt" 2>/dev/null || echo 0)
    echo "  Contig classification:"
    echo "    Chromosomal:   $CHROM_COUNT"
    echo "    Plasmid:       $PLASMID_CONTIGS"
    echo "    Unclassified:  $UNCLASSIFIED"
    echo ""
fi

# ── Parse MOB-typer into clean summary ────────────────────────────────────────
if [[ -f "$MOBTYPER_OUT" ]]; then
    PLASMID_LINES=$(( $(wc -l < "$MOBTYPER_OUT") - 1 ))

    if [[ "$PLASMID_LINES" -gt 0 ]]; then
        python3 - "$MOBTYPER_OUT" "$SUMMARY_TSV" "$SAMPLE" << 'PYEOF'
import csv, sys

infile  = sys.argv[1]
outfile = sys.argv[2]
sample  = sys.argv[3]

RED    = '\033[0;31m'
YELLOW = '\033[1;33m'
GREEN  = '\033[0;32m'
RESET  = '\033[0m'

KEY_COLS = [
    "file_id",
    "rep_type(s)",
    "relaxase_type(s)",
    "predicted_mobility",
    "mge_type(s)",
    "primary_cluster_id",
    "secondary_cluster_id",
    "predicted_host_range_overall_rank",
    "observed_mobility",
    "size",
    "gc",
    "num_contigs",
]

rows = []
with open(infile) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rows.append(row)

# Write clean TSV
with open(outfile, 'w') as out:
    out.write('sample\t' + '\t'.join(KEY_COLS) + '\n')
    for row in rows:
        vals = [sample] + [row.get(c, '-') or '-' for c in KEY_COLS]
        out.write('\t'.join(vals) + '\n')

# Print human-readable summary
for i, row in enumerate(rows, 1):
    mob      = row.get("predicted_mobility", "-") or "-"
    rep      = row.get("rep_type(s)", "-") or "-"
    relax    = row.get("relaxase_type(s)", "-") or "-"
    cluster  = row.get("primary_cluster_id", "-") or "-"
    hostrank = row.get("predicted_host_range_overall_rank", "-") or "-"
    size     = row.get("size", "-") or "-"
    gc       = row.get("gc", "-") or "-"
    fid      = row.get("file_id", f"plasmid_{i}")
    contigs  = row.get("num_contigs", "-") or "-"

    if mob == "CONJUGATIVE":
        mob_str = f"{RED}{mob}{RESET}"
        risk    = f"{RED}HIGH — can self-transfer to other bacteria{RESET}"
    elif mob == "MOBILIZABLE":
        mob_str = f"{YELLOW}{mob}{RESET}"
        risk    = f"{YELLOW}MODERATE — needs helper plasmid to transfer{RESET}"
    else:
        mob_str = f"{GREEN}{mob}{RESET}"
        risk    = f"{GREEN}LOW — cannot transfer by conjugation{RESET}"

    print(f"  {'─'*54}")
    print(f"  Plasmid {i}: {fid}")
    print(f"    Size:        {size} bp  |  GC: {gc}%  |  Contigs: {contigs}")
    print(f"    Replicon:    {rep}")
    print(f"    Relaxase:    {relax}")
    print(f"    Mobility:    {mob_str}")
    print(f"    Risk:        {risk}")
    print(f"    Host range:  {hostrank}")
    print(f"    Cluster:     {cluster}")
print(f"  {'─'*54}")
PYEOF

    else
        warn "No plasmids detected in this assembly"
        printf "sample\tfile_id\trep_type(s)\trelaxase_type(s)\tpredicted_mobility\tcluster\n" \
            > "$SUMMARY_TSV"
    fi
fi

echo ""
echo "────────────────────────────────────────────────────────"
echo "  Output directory:    $OUTPUT_DIR"
echo "  Clean summary:       $SUMMARY_TSV"
echo "  Full raw output:     $MOBTYPER_OUT"
echo "  Contig report:       $OUTPUT_DIR/contig_report.txt"
[[ "$PLASMID_COUNT" -gt 0 ]] && \
echo "  Plasmid sequences:   $OUTPUT_DIR/plasmid_*.fasta ($PLASMID_COUNT files)"
echo "  Chromosome:          $OUTPUT_DIR/chromosome.fasta"
echo "────────────────────────────────────────────────────────"
