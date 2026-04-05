#!/usr/bin/env bash
# =============================================================================
# assembly-qc.sh
# Runs: QUAST → CheckM → MultiQC
# MultiQC scans the whole sample folder covering fastp, FastQC, QUAST, CheckM
#
# Usage:
#   ./assembly-qc.sh <sample_dir>
#
# Run after bacterial-pipeline.sh has completed.
# Assembly/ subfolder is detected automatically.
#
# Output:
#   <sample_dir>/AssemblyQC/QUAST/
#   <sample_dir>/AssemblyQC/CheckM/
#   <sample_dir>/QC/MultiQC/multiqc_report.html   ← combined report
#
# Optional env vars:
#   THREADS — number of threads (default: 8)
# =============================================================================

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()     { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
success() { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${RESET}"; }
die()     { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ $*${RESET}" >&2; exit 1; }
header()  { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${RESET}"
            echo -e "${BOLD}${CYAN}  $*${RESET}"
            echo -e "${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; }

[[ $# -lt 1 ]] && die "Usage: $0 <sample_dir>"

SAMPLE_DIR="${1%/}"
THREADS="${THREADS:-8}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

# ── Auto-detect assembly file ──────────────────────────────────────────────────
ASSEMBLY_DIR="$SAMPLE_DIR/Assembly"
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly/ not found in $SAMPLE_DIR\nRun bacterial-pipeline.sh first"

if [[ -f "$ASSEMBLY_DIR/scaffolds.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/scaffolds.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fasta"
    log "Using contigs.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fa" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fa"
    log "Using contigs.fa"
else
    die "No assembly file found in $ASSEMBLY_DIR\nExpected: scaffolds.fasta, contigs.fasta, or contigs.fa"
fi

log "Sample:   $(basename "$SAMPLE_DIR")"
log "Assembly: $(basename "$ASSEMBLY_FA")"

# ── Output folders ─────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/AssemblyQC"
QUAST_DIR="$OUTPUT_DIR/QUAST"
CHECKM_DIR="$OUTPUT_DIR/CheckM"
MULTIQC_DIR="$SAMPLE_DIR/QC/MultiQC"

mkdir -p "$QUAST_DIR" "$CHECKM_DIR" "$MULTIQC_DIR"

STEP_START=0
start_timer() { STEP_START=$SECONDS; }
elapsed()     { echo $(( SECONDS - STEP_START ))s; }

# =============================================================================
header "Step 1/3 — QUAST: assembly statistics"
# =============================================================================
start_timer

pixi run -- quast.py "$ASSEMBLY_FA" -o "$QUAST_DIR"

N50=$(grep "N50" "$QUAST_DIR/report.tsv" 2>/dev/null | awk '{print $2}' || echo "unknown")
LARGEST=$(grep "Largest contig" "$QUAST_DIR/report.tsv" 2>/dev/null | awk '{print $3}' || echo "unknown")
CONTIGS=$(grep "# contigs " "$QUAST_DIR/report.tsv" 2>/dev/null | awk '{print $NF}' || echo "unknown")
TOTAL=$(grep "Total length " "$QUAST_DIR/report.tsv" 2>/dev/null | awk '{print $NF}' || echo "unknown")

success "QUAST done ($(elapsed))"
echo ""
echo "  Contigs:        $CONTIGS"
echo "  Total length:   $TOTAL bp"
echo "  N50:            $N50 bp"
echo "  Largest contig: $LARGEST bp"

# =============================================================================
header "Step 2/3 — CheckM: completeness and contamination"
# =============================================================================
start_timer

pixi run --environment checkm -- checkm lineage_wf \
  "$ASSEMBLY_DIR" \
  "$CHECKM_DIR" \
  -t "$THREADS" \
  -x fasta

COMPLETENESS=$(grep -oP "'Completeness': \K[0-9.]+" \
  "$CHECKM_DIR/storage/bin_stats_ext.tsv" 2>/dev/null | head -1 || echo "unknown")
CONTAMINATION=$(grep -oP "'Contamination': \K[0-9.]+" \
  "$CHECKM_DIR/storage/bin_stats_ext.tsv" 2>/dev/null | head -1 || echo "unknown")

success "CheckM done ($(elapsed))"
echo ""
echo "  Completeness:   ${COMPLETENESS}%"
echo "  Contamination:  ${CONTAMINATION}%"

pixi run --environment checkm -- checkm qa \
  "$CHECKM_DIR/lineage.ms" \
  "$CHECKM_DIR" \
  -o 1 \
  > "$CHECKM_DIR/checkm_summary.txt"

cat "$CHECKM_DIR/checkm_summary.txt"

# =============================================================================
header "Step 3/3 — MultiQC: combined QC report"
# =============================================================================
# Scans whole sample dir — picks up fastp, FastQC, QUAST and CheckM
start_timer

pixi run -- multiqc \
  "$SAMPLE_DIR" \
  --outdir "$MULTIQC_DIR" \
  --force

success "MultiQC done ($(elapsed))"

# =============================================================================
header "Done"
# =============================================================================
echo ""
echo "  Contigs:         $CONTIGS"
echo "  Total length:    $TOTAL bp"
echo "  N50:             $N50 bp"
echo "  Largest contig:  $LARGEST bp"
echo "  Completeness:    ${COMPLETENESS}%"
echo "  Contamination:   ${CONTAMINATION}%"
echo ""
echo "  QUAST report:    $QUAST_DIR/report.html"
echo "  CheckM output:   $CHECKM_DIR"
echo "  MultiQC report:  $MULTIQC_DIR/multiqc_report.html"