#!/usr/bin/env bash
# =============================================================================
# run-snippy.sh
# Runs Snippy — SNP calling against a reference genome
#
# Usage:
#   ./run-snippy.sh <sample_dir> <reference>
#
# Arguments:
#   sample_dir  — sample folder containing a Trimmed/ subfolder
#   reference   — reference genome (.gbk, .gb, .gbff, or .fasta/.fa/.fna)
#
# Output:
#   <sample_dir>/SNPs/
#     snps.vcf        All variants (VCF format)
#     snps.tab        Variants as tab-delimited table (easier to read)
#     snps.html       Visual HTML report
#     snps.txt        Summary statistics
#     snps.aligned.fa Masked consensus — use for phylogeny
#
# Optional env vars:
#   THREADS   — number of threads (default: 8)
#   MIN_COV   — minimum coverage to call a SNP (default: 10)
#   MIN_FRAC  — minimum allele frequency (default: 0.9)
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
if [[ $# -lt 2 ]]; then
    echo "============================================================"
    echo "           S N P   C A L L I N G   —   S N I P P Y         "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> <reference>"
    echo ""
    echo -e "${BOLD}Arguments:${RESET}"
    echo "  sample_dir  Sample folder containing Trimmed/ subfolder"
    echo "  reference   Reference genome file:"
    echo "                .gbk / .gb / .gbff  — GenBank format (recommended,"
    echo "                                       includes gene annotations)"
    echo "                .fasta / .fa / .fna  — FASTA format"
    echo ""
    echo -e "${BOLD}Optional env vars:${RESET}"
    echo "  THREADS   Number of threads (default: 8)"
    echo "  MIN_COV   Minimum read depth to call SNP (default: 10)"
    echo "  MIN_FRAC  Minimum allele frequency 0-1 (default: 0.9)"
    echo ""
    echo -e "${BOLD}Output:${RESET}"
    echo "  <sample_dir>/SNPs/"
    echo "    snps.tab        Variant table — easiest to read"
    echo "    snps.vcf        VCF format"
    echo "    snps.html       Visual report"
    echo "    snps.aligned.fa Masked consensus for phylogeny"
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/ reference.gbk"
    echo "  $0 results/001k/ /databases/refs/Kpneumoniae_HS11286.gbff"
    echo "  MIN_COV=20 $0 results/001k/ reference.gbk"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
REFERENCE="$2"
THREADS="${THREADS:-8}"
MIN_COV="${MIN_COV:-10}"
MIN_FRAC="${MIN_FRAC:-0.9}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -f "$REFERENCE" ]]  || die "Reference file not found: $REFERENCE"

# ── Find trimmed reads ────────────────────────────────────────────────────────
TRIM_DIR="$SAMPLE_DIR/Trimmed"
[[ -d "$TRIM_DIR" ]] || die "Trimmed directory not found: $TRIM_DIR\nRun bacterial-pipeline.sh first"

R1_FILE=$(find "$TRIM_DIR" -maxdepth 1 -name "*_R1*.fastq.gz" | sort | head -1)
[[ -n "$R1_FILE" ]] || die "No R1 fastq.gz found in $TRIM_DIR"
R2_FILE="${R1_FILE/_R1/_R2}"
[[ -f "$R2_FILE" ]] || die "R2 not found: $R2_FILE"

SAMPLE=$(basename "$SAMPLE_DIR")
log "Sample:    $SAMPLE"
log "R1:        $(basename "$R1_FILE")"
log "R2:        $(basename "$R2_FILE")"
log "Reference: $(basename "$REFERENCE")"
log "Min depth: $MIN_COV x"
log "Min freq:  $MIN_FRAC"

# ── Output folder ─────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/SNPs"
mkdir -p "$OUTPUT_DIR"
log "Output:    $OUTPUT_DIR"

STEP_START=$SECONDS

# =============================================================================
header "Snippy — SNP calling vs $(basename "$REFERENCE")"
# =============================================================================

pixi run -- snippy \
  --cpus "$THREADS" \
  --outdir "$OUTPUT_DIR" \
  --ref "$REFERENCE" \
  --R1 "$R1_FILE" \
  --R2 "$R2_FILE" \
  --mincov "$MIN_COV" \
  --minfrac "$MIN_FRAC" \
  --force

# ── Summary ───────────────────────────────────────────────────────────────────
ELAPSED=$(( SECONDS - STEP_START ))s

SNP_COUNT=$(grep -v "^#" "$OUTPUT_DIR/snps.vcf" 2>/dev/null \
  | awk '$5!="." && $4!="-"' | wc -l || echo "0")
INDEL_COUNT=$(grep -v "^#" "$OUTPUT_DIR/snps.vcf" 2>/dev/null \
  | awk 'length($4)!=length($5) && $5!="."' | wc -l || echo "0")

success "Snippy done ($ELAPSED)"
echo ""

# Print the snps.txt summary if it exists
if [[ -f "$OUTPUT_DIR/snps.txt" ]]; then
    cat "$OUTPUT_DIR/snps.txt"
    echo ""
fi

if [[ "$SNP_COUNT" -gt 0 ]]; then
    echo ""
    log "Top variants (snps.tab):"
    echo ""
    head -6 "$OUTPUT_DIR/snps.tab" | column -t -s $'\t'
fi

echo ""
echo "────────────────────────────────────────"
echo "  Reference:    $(basename "$REFERENCE")"
echo "  SNPs:         $SNP_COUNT"
echo "  Indels:       $INDEL_COUNT"
echo ""
echo "  Full table:   $OUTPUT_DIR/snps.tab"
echo "  VCF:          $OUTPUT_DIR/snps.vcf"
echo "  HTML report:  $OUTPUT_DIR/snps.html"
echo "  For phylogeny: $OUTPUT_DIR/snps.aligned.fa"
echo "────────────────────────────────────────"
