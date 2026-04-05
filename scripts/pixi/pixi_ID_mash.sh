#!/usr/bin/env bash
# =============================================================================
# pixi_ID_mash.sh
# Runs Mash screen for rapid species identification from trimmed reads
#
# Usage:
#   ./pixi_ID_mash.sh <sample_dir>
#
# The sample_dir should contain a Trimmed/ subfolder with *_R1*.fastq.gz.
# Results are written to <sample_dir>/ID_Mash/
#
# Required env vars:
#   MASH_DB — path to Mash reference sketch (.msh file)
#             Download: wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
#
# Optional env vars:
#   THREADS — number of threads (default: 8)
#   MASH_MIN_IDENTITY — minimum identity to report (default: 0.90)
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

# ── No arguments — show help ──────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "           M A S H   I D E N T I F I C A T I O N           "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir>"
    echo ""
    echo -e "${BOLD}Required env vars:${RESET}"
    echo "  MASH_DB   — path to .msh reference sketch"
    echo "              export MASH_DB=/path/to/refseq.msh"
    echo ""
    echo -e "${BOLD}Optional env vars:${RESET}"
    echo "  THREADS           — number of threads (default: 8)"
    echo "  MASH_MIN_IDENTITY — minimum identity to report (default: 0.90)"
    echo ""
    echo -e "${BOLD}Download Mash database (once, ~2GB):${RESET}"
    echo "  wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh \\"
    echo "       -O /path/to/refseq.msh"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  export MASH_DB=/databases/refseq.msh"
    echo "  $0 results/001k/"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
THREADS="${THREADS:-8}"
MASH_MIN_IDENTITY="${MASH_MIN_IDENTITY:-0.90}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -n "${MASH_DB:-}" ]] || die "MASH_DB is not set.\nExport it first:\n  export MASH_DB=/path/to/refseq.msh"
[[ -f "$MASH_DB" ]]     || die "MASH_DB file not found: $MASH_DB"

# ── Find trimmed R1 reads ─────────────────────────────────────────────────────
TRIM_DIR="$SAMPLE_DIR/Trimmed"
[[ -d "$TRIM_DIR" ]] || die "Trimmed directory not found: $TRIM_DIR"

R1_FILE=$(find "$TRIM_DIR" -maxdepth 1 -name "*_R1*.fastq.gz" | sort | head -1)
[[ -n "$R1_FILE" ]] || die "No R1 fastq.gz found in $TRIM_DIR"
[[ -f "$R1_FILE" ]]  || die "R1 file not found: $R1_FILE"

log "Reads:    $(basename "$R1_FILE")"
log "Database: $(basename "$MASH_DB")"
log "Min identity: $MASH_MIN_IDENTITY"

# ── Output folder ─────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/ID_Mash"
OUTPUT_RAW="$OUTPUT_DIR/mash_screen_raw.tsv"
OUTPUT_TSV="$OUTPUT_DIR/mash_screen.tsv"
mkdir -p "$OUTPUT_DIR"
log "Output:   $OUTPUT_DIR"

# =============================================================================
header "Mash screen — rapid species identification"
# =============================================================================

pixi run --environment identification -- mash screen \
  -w \
  -p "$THREADS" \
  -i "$MASH_MIN_IDENTITY" \
  "$MASH_DB" \
  "$R1_FILE" \
  > "$OUTPUT_RAW"

# ── Sort by identity (best matches first) and add header ─────────────────────
{
  printf "identity\tshared_hashes\tmedian_multiplicity\tp_value\tquery_id\tquery_comment\n"
  sort -rn "$OUTPUT_RAW"
} > "$OUTPUT_TSV"

rm "$OUTPUT_RAW"

# ── Show top 10 hits ──────────────────────────────────────────────────────────
success "Mash done"
echo ""
log "Top 10 matches:"
echo ""
head -11 "$OUTPUT_TSV" | column -t -s $'\t'
echo ""

TOTAL_HITS=$(( $(wc -l < "$OUTPUT_TSV") - 1 ))
log "Total matches: $TOTAL_HITS  (identity >= $MASH_MIN_IDENTITY)"
log "Results saved to: $OUTPUT_TSV"
