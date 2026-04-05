#!/usr/bin/env bash
# =============================================================================
# pixi_kraken2.sh
# Runs Kraken2 (taxonomic classification) + Bracken (species abundance)
#
# Usage:
#   ./pixi_kraken2.sh <sample_dir>
#
# The sample_dir should contain a Trimmed/ subfolder with *_R1*.fastq.gz.
# Results are written to <sample_dir>/ID_Kraken2/
#
# Kraken2:  assigns each read to a taxon — good for species ID + contamination
# Bracken:  re-estimates species abundance from Kraken2 results
#
# Required env vars:
#   KRAKEN2_DB — path to Kraken2 database directory
#                Download: https://benlangmead.github.io/aws-indexes/k2
#                Recommended: PlusPF (standard + protozoa + fungi)
#
# Optional env vars:
#   THREADS         — number of threads (default: 8)
#   BRACKEN_LENGTH  — read length for Bracken (default: 150)
#   MIN_HIT_GROUPS  — min hit groups for classification (default: 3)
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
    echo "     K R A K E N 2   +   B R A C K E N                     "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir>"
    echo ""
    echo -e "${BOLD}Required env vars:${RESET}"
    echo "  KRAKEN2_DB  — path to Kraken2 database"
    echo "                export KRAKEN2_DB=/databases/kraken2_db"
    echo ""
    echo -e "${BOLD}Optional env vars:${RESET}"
    echo "  THREADS         Number of threads (default: 8)"
    echo "  BRACKEN_LENGTH  Read length for Bracken (default: 150)"
    echo "  MIN_HIT_GROUPS  Min hit groups (default: 3)"
    echo ""
    echo -e "${BOLD}Download database (once):${RESET}"
    echo "  # PlusPF ~50GB — bacteria, archaea, virus, plasmid, fungi, protozoa"
    echo "  wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz"
    echo "  mkdir -p /databases/kraken2_db"
    echo "  tar -xzf k2_pluspf_20240605.tar.gz -C /databases/kraken2_db"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  export KRAKEN2_DB=/databases/kraken2_db"
    echo "  $0 results/001k/"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
THREADS="${THREADS:-8}"
BRACKEN_LENGTH="${BRACKEN_LENGTH:-150}"
MIN_HIT_GROUPS="${MIN_HIT_GROUPS:-3}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -n "${KRAKEN2_DB:-}" ]] || die "KRAKEN2_DB is not set.\nExport it first:\n  export KRAKEN2_DB=/databases/kraken2_db"
[[ -d "$KRAKEN2_DB" ]]     || die "KRAKEN2_DB directory not found: $KRAKEN2_DB"

# ── Find trimmed R1/R2 reads ──────────────────────────────────────────────────
TRIM_DIR="$SAMPLE_DIR/Trimmed"
[[ -d "$TRIM_DIR" ]] || die "Trimmed directory not found: $TRIM_DIR\nRun bacterial-pipeline.sh first"

R1_FILE=$(find "$TRIM_DIR" -maxdepth 1 -name "*_R1*.fastq.gz" | sort | head -1)
[[ -n "$R1_FILE" ]] || die "No R1 fastq.gz found in $TRIM_DIR"
R2_FILE="${R1_FILE/_R1/_R2}"
[[ -f "$R2_FILE" ]] || die "R2 not found: $R2_FILE"

SAMPLE=$(basename "$SAMPLE_DIR")
log "Sample:   $SAMPLE"
log "R1:       $(basename "$R1_FILE")"
log "R2:       $(basename "$R2_FILE")"
log "Database: $KRAKEN2_DB"

# ── Output folder ─────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/ID_Kraken2"
mkdir -p "$OUTPUT_DIR"
log "Output:   $OUTPUT_DIR"

STEP_START=0
start_timer() { STEP_START=$SECONDS; }
elapsed()     { echo $(( SECONDS - STEP_START ))s; }

# =============================================================================
header "Step 1/2 — Kraken2: taxonomic classification"
# =============================================================================
start_timer

KRAKEN_REPORT="$OUTPUT_DIR/kraken2_report.txt"
KRAKEN_OUTPUT="$OUTPUT_DIR/kraken2_output.txt"

pixi run --environment identification -- kraken2 \
  --db "$KRAKEN2_DB" \
  --threads "$THREADS" \
  --paired \
  --gzip-compressed \
  --report "$KRAKEN_REPORT" \
  --minimum-hit-groups "$MIN_HIT_GROUPS" \
  --report-minimizer-data \
  "$R1_FILE" "$R2_FILE" \
  > "$KRAKEN_OUTPUT"

# Parse summary stats from report
CLASSIFIED=$(grep -P "^\s+\d+\s+\d+\s+\d+\s+U\s+0\s+unclassified" "$KRAKEN_REPORT" \
  | awk '{print $1}' || echo "unknown")
UNCLASSIFIED=$(grep "unclassified" "$KRAKEN_REPORT" \
  | awk '{print $1}' | head -1 || echo "unknown")

success "Kraken2 done ($(elapsed))"
echo ""
echo "  Classified:   $(grep -v "unclassified" "$KRAKEN_REPORT" | awk 'NR==1{print $1}')%"
echo "  Unclassified: $(grep "unclassified" "$KRAKEN_REPORT" | awk '{print $1}' | head -1)%"
echo ""

# Show top species hits
log "Top species hits:"
echo ""
grep -P "\s+S\s+" "$KRAKEN_REPORT" \
  | sort -rn -k1 \
  | head -10 \
  | awk '{printf "  %6.2f%%  %s\n", $1, $NF}' \
  | sed 's/_/ /g'
echo ""

# Contamination check — warn if top hit is not dominant
TOP_PCT=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 | head -1 | awk '{print $1}')
SECOND_PCT=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 | sed -n '2p' | awk '{print $1}')

if (( $(echo "$SECOND_PCT > 5" | bc -l) )); then
    warn "Possible contamination — second species at ${SECOND_PCT}%"
fi

# =============================================================================
header "Step 2/2 — Bracken: species abundance re-estimation"
# =============================================================================
start_timer

BRACKEN_OUTPUT="$OUTPUT_DIR/bracken_species.txt"
BRACKEN_REPORT="$OUTPUT_DIR/bracken_report.txt"

# Check if Bracken database files exist for this read length
BRACKEN_DB="$KRAKEN2_DB/database${BRACKEN_LENGTH}mers.kmer_distrib"
if [[ ! -f "$BRACKEN_DB" ]]; then
    warn "Bracken database not found for ${BRACKEN_LENGTH}bp reads: $BRACKEN_DB"
    warn "Available kmer_distrib files:"
    ls "$KRAKEN2_DB"/*.kmer_distrib 2>/dev/null | xargs -I{} basename {} || echo "  none found"
    warn "Skipping Bracken — run with correct BRACKEN_LENGTH"
else
    pixi run --environment identification -- bracken \
      -d "$KRAKEN2_DB" \
      -i "$KRAKEN_REPORT" \
      -o "$BRACKEN_OUTPUT" \
      -w "$BRACKEN_REPORT" \
      -r "$BRACKEN_LENGTH" \
      -l S \
      -t 10

    success "Bracken done ($(elapsed))"
    echo ""
    log "Top species (Bracken re-estimated):"
    echo ""
    # Sort by fraction_total_reads (col 7), skip header
    awk -F'\t' 'NR>1{print}' "$BRACKEN_OUTPUT" \
      | sort -t$'\t' -k7 -rn \
      | head -10 \
      | awk -F'\t' '{printf "  %6.2f%%  %s\n", $7*100, $1}'
    echo ""
fi

# =============================================================================
# Contamination assessment
# =============================================================================
header "Contamination assessment"

# Use Bracken output if available, otherwise fall back to Kraken2 report
if [[ -f "${BRACKEN_OUTPUT:-}" ]]; then
    TOP_SPECIES=$(awk -F'\t' 'NR>1{print}' "$BRACKEN_OUTPUT" \
      | sort -t$'\t' -k7 -rn | head -1 | cut -f1)
    TOP_PCT=$(awk -F'\t' 'NR>1{print}' "$BRACKEN_OUTPUT" \
      | sort -t$'\t' -k7 -rn | head -1 | awk -F'\t' '{printf "%.1f", $7*100}')
    SECOND_SPECIES=$(awk -F'\t' 'NR>1{print}' "$BRACKEN_OUTPUT" \
      | sort -t$'\t' -k7 -rn | sed -n '2p' | cut -f1)
    SECOND_PCT=$(awk -F'\t' 'NR>1{print}' "$BRACKEN_OUTPUT" \
      | sort -t$'\t' -k7 -rn | sed -n '2p' | awk -F'\t' '{printf "%.1f", $7*100}')
else
    TOP_SPECIES=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 \
      | head -1 | awk '{print $NF}' | sed 's/_/ /g')
    TOP_PCT=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 \
      | head -1 | awk '{print $1}')
    SECOND_SPECIES=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 \
      | sed -n '2p' | awk '{print $NF}' | sed 's/_/ /g')
    SECOND_PCT=$(grep -P "\s+S\s+" "$KRAKEN_REPORT" | sort -rn -k1 \
      | sed -n '2p' | awk '{print $1}')
fi

UNCLASSIFIED_PCT=$(grep "unclassified" "$KRAKEN_REPORT" \
  | awk '{print $1}' | head -1)
HUMAN_PCT=$(grep -P "Homo sapiens" "$KRAKEN_REPORT" \
  | awk '{print $1}' | head -1 || echo "0")
HUMAN_PCT="${HUMAN_PCT:-0}"

echo ""
echo "  Primary species:   ${TOP_SPECIES} (${TOP_PCT}%)"
echo "  Secondary species: ${SECOND_SPECIES} (${SECOND_PCT}%)"
echo "  Unclassified:      ${UNCLASSIFIED_PCT}%"
echo "  Human reads:       ${HUMAN_PCT}%"
echo ""

# Verdict
if (( $(echo "$TOP_PCT >= 95" | bc -l) )); then
    echo -e "  Purity:  ${GREEN}CLEAN — dominant species ≥95%${RESET}"
    PURITY_VERDICT="CLEAN"
elif (( $(echo "$TOP_PCT >= 80" | bc -l) )); then
    echo -e "  Purity:  ${GREEN}LIKELY CLEAN — dominant species ≥80%${RESET}"
    PURITY_VERDICT="LIKELY CLEAN"
elif (( $(echo "$TOP_PCT >= 60" | bc -l) )); then
    echo -e "  Purity:  ${YELLOW}⚠ POSSIBLE CONTAMINATION — dominant species ${TOP_PCT}%${RESET}"
    PURITY_VERDICT="POSSIBLE CONTAMINATION"
else
    echo -e "  Purity:  ${RED}✘ LIKELY CONTAMINATED OR MIXED CULTURE — dominant species only ${TOP_PCT}%${RESET}"
    PURITY_VERDICT="LIKELY CONTAMINATED"
fi

if (( $(echo "$HUMAN_PCT > 5" | bc -l) )); then
    echo -e "  ${YELLOW}⚠ High human reads (${HUMAN_PCT}%) — possible DNA extraction contamination${RESET}"
fi

echo ""
echo "  Notes:"
echo "  Cross-classification noise: closely related species (e.g. K. variicola vs"
echo "  K. pneumoniae, ~97% ANI) share most k-mers. Reads from one species may be"
echo "  assigned to the other. Minor signals (<20%) within the same genus are"
echo "  usually noise, not contamination. Use FastANI for definitive ID."
echo ""
echo "  Rules of thumb:"
echo "  ┌─────────────────────────────────────────────────────────┐"
echo "  │ >95% one species        →  Clean isolate                │"
echo "  │ 80-95%, rest same genus →  Likely clean (noise)         │"
echo "  │ >10% unrelated species  →  Possible contamination       │"
echo "  │ >5% human reads         →  Extraction contamination     │"
echo "  │ Multiple species >5%    →  Mixed culture/contamination  │"
echo "  └─────────────────────────────────────────────────────────┘"

# Save contamination summary to TSV
CONTAM_TSV="$OUTPUT_DIR/contamination_summary.tsv"
{
echo -e "metric\tvalue"
echo -e "sample\t$SAMPLE"
echo -e "primary_species\t$TOP_SPECIES"
echo -e "primary_pct\t$TOP_PCT"
echo -e "secondary_species\t$SECOND_SPECIES"
echo -e "secondary_pct\t$SECOND_PCT"
echo -e "unclassified_pct\t$UNCLASSIFIED_PCT"
echo -e "human_pct\t$HUMAN_PCT"
echo -e "purity_verdict\t$PURITY_VERDICT"
} > "$CONTAM_TSV"

echo ""
log "Contamination summary saved to: $CONTAM_TSV"