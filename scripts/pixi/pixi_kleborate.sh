#!/usr/bin/env bash
# =============================================================================
# pixi_kleborate.sh
# Runs Kleborate on an assembly
#
# Usage:
#   ./pixi_kleborate.sh <sample_dir> [preset]
#
# Arguments:
#   sample_dir  — sample folder containing an Assembly/ subfolder
#   preset      — Kleborate preset (default: kpsc)
#
# Output:
#   results/001k/Kleborate/kleborate_kpsc.tsv
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
    echo "         K L E B O R A T E   P I P E L I N E               "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> [preset]"
    echo ""
    echo -e "${BOLD}Available presets:${RESET}"
    echo "  kpsc        K. pneumoniae species complex — AMR, virulence, K/O typing (default)"
    echo "  kosc        K. oxytoca species complex — MLST, AMR"
    echo "  escherichia Escherichia / Shigella — MLST, pathotyping"
    echo "  species     Species identification only (fastest)"
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/"
    echo "  $0 results/001k/ kosc"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
PRESET="${2:-kpsc}"
THREADS="${THREADS:-8}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

# ── Validate preset ───────────────────────────────────────────────────────────
case "$PRESET" in
    kpsc|kosc|escherichia|species) ;;
    *) die "Unknown preset: $PRESET\nValid options: kpsc, kosc, escherichia, species" ;;
esac

# ── Find assembly file ────────────────────────────────────────────────────────
ASSEMBLY_DIR="$SAMPLE_DIR/Assembly"
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly directory not found: $ASSEMBLY_DIR"

if [[ -f "$ASSEMBLY_DIR/scaffolds.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/scaffolds.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fa" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fa"
else
    die "No assembly file found in $ASSEMBLY_DIR\nExpected: scaffolds.fasta, contigs.fasta, or contigs.fa"
fi

log "Assembly: $(basename "$ASSEMBLY_FA")"
log "Preset:   $PRESET"

# ── Output paths ──────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/Kleborate"
TMP_DIR="$OUTPUT_DIR/.tmp_$(date +%s)"
RESULT_TSV="$OUTPUT_DIR/kleborate_${PRESET}.tsv"
mkdir -p "$TMP_DIR"
log "Output:   $RESULT_TSV"

# =============================================================================
header "Kleborate — typing, virulence and AMR ($PRESET)"
# =============================================================================

pixi run --environment identification -- kleborate \
  --preset "$PRESET" \
  --trim_headers \
  --threads "$THREADS" \
  -a "$ASSEMBLY_FA" \
  -o "$TMP_DIR"

# ── Merge output files and clean up headers ───────────────────────────────────
[[ -d "$TMP_DIR" ]] || die "Kleborate produced no output"

first=true
for f in "$TMP_DIR"/*.txt; do
    [[ -f "$f" ]] || die "No result .txt files found in output"
    if [[ "$first" == true ]]; then
        cat "$f" > "$RESULT_TSV"
        first=false
    else
        # Append without header for additional files
        tail -n +2 "$f" >> "$RESULT_TSV"
    fi
done

# Clean up verbose v3 column headers
sed -i '1s/[^[:space:]]*mlst__ST/ST/g'                          "$RESULT_TSV"
sed -i '1s/[^[:space:]]*species_complex__species/species/g'     "$RESULT_TSV"
sed -i '1s/[^[:space:]]*__/\t/g'                                "$RESULT_TSV"

rm -rf "$TMP_DIR"

# ── Show results ──────────────────────────────────────────────────────────────
success "Kleborate done"
echo ""
column -t -s $'\t' "$RESULT_TSV"
echo ""
log "Results saved to: $RESULT_TSV"
