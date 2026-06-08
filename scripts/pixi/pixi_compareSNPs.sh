#!/usr/bin/env bash
# =============================================================================
# pixi_compareSNP.sh
# Compares two isolates using SNP distance, genome size and contig stats
# to predict if they are the same clone
#
# Usage:
#   ./pixi_compareSNP.sh <sample1_dir> <sample2_dir> [reference]
#
# Arguments:
#   sample1_dir  — sample folder containing SNPs/ and Assembly/
#   sample2_dir  — sample folder containing SNPs/ and Assembly/
#   reference    — (optional) reference genome used for snippy
#                  Default: /databases/snippy_refs/latest
#
# Requirements:
#   Both samples must have SNPs/ folder from run-snippy.sh,
#   run against the same reference.
#
# Output:
#   comparison_<sample1>_vs_<sample2>.tsv
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
    echo "       I S O L A T E   C O M P A R I S O N                 "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample1_dir> <sample2_dir> [reference]"
    echo ""
    echo -e "${BOLD}Arguments:${RESET}"
    echo "  sample1_dir  Sample folder with SNPs/ and Assembly/"
    echo "  sample2_dir  Sample folder with SNPs/ and Assembly/"
    echo "  reference    (optional) reference genome"
    echo "               Default: /databases/snippy_refs/latest"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  $0 results/001k/ results/002k/"
    echo "  $0 results/001k/ results/002k/ reference.gbk"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE1_DIR="${1%/}"
SAMPLE2_DIR="${2%/}"
REFERENCE="${3:-/databases/snippy_refs/latest}"

# Convert to absolute paths
SAMPLE1_DIR=$(realpath "$SAMPLE1_DIR")
SAMPLE2_DIR=$(realpath "$SAMPLE2_DIR")
REFERENCE=$(realpath "$REFERENCE")

SAMPLE1=$(basename "$SAMPLE1_DIR")
SAMPLE2=$(basename "$SAMPLE2_DIR")

[[ -d "$SAMPLE1_DIR" ]]      || die "Sample 1 not found: $SAMPLE1_DIR"
[[ -d "$SAMPLE2_DIR" ]]      || die "Sample 2 not found: $SAMPLE2_DIR"
[[ -f "$REFERENCE" ]]        || die "Reference not found: $REFERENCE\nDefault path: /databases/snippy_refs/latest"
[[ -d "$SAMPLE1_DIR/SNPs" ]] || die "SNPs/ not found in $SAMPLE1_DIR — run run-snippy.sh first"
[[ -d "$SAMPLE2_DIR/SNPs" ]] || die "SNPs/ not found in $SAMPLE2_DIR — run run-snippy.sh first"
[[ -d "$SAMPLE1_DIR/Assembly" ]] || die "Assembly/ not found in $SAMPLE1_DIR"
[[ -d "$SAMPLE2_DIR/Assembly" ]] || die "Assembly/ not found in $SAMPLE2_DIR"

REFERENCE=$(realpath "$REFERENCE")
OUTPUT_TSV="comparison_${SAMPLE1}_vs_${SAMPLE2}.tsv"
CORE_PREFIX="/tmp/core_${SAMPLE1}_${SAMPLE2}"
SYMDIR="/tmp/snippy_compare_${SAMPLE1}_${SAMPLE2}"

log "Sample 1:  $SAMPLE1"
log "Sample 2:  $SAMPLE2"
log "Reference: $(basename "$REFERENCE")"

# =============================================================================
header "Step 1/3 — SNP distance via snippy-core + snp-dists"
# =============================================================================

# Create symlinks named after samples so snippy-core labels them correctly
mkdir -p "$SYMDIR"
ln -sf "$(realpath "$SAMPLE1_DIR/SNPs")" "$SYMDIR/$SAMPLE1"
ln -sf "$(realpath "$SAMPLE2_DIR/SNPs")" "$SYMDIR/$SAMPLE2"

pixi run -- snippy-core \
  --ref "$REFERENCE" \
  --prefix "$CORE_PREFIX" \
  "$SYMDIR/$SAMPLE1" \
  "$SYMDIR/$SAMPLE2"

SNP_DISTANCE=$(pixi run -- snp-dists "${CORE_PREFIX}.full.aln" 2>/dev/null \
  | awk 'NR==2{print $3}')

rm -rf "$SYMDIR" "${CORE_PREFIX}"*

log "SNP distance: $SNP_DISTANCE SNPs"

# =============================================================================
header "Step 2/3 — Assembly statistics"
# =============================================================================

get_assembly_stats() {
    local dir="$1"
    local fa=""

    if [[ -f "$dir/Assembly/scaffolds.fasta" ]]; then
        fa="$dir/Assembly/scaffolds.fasta"
    elif [[ -f "$dir/Assembly/contigs.fasta" ]]; then
        fa="$dir/Assembly/contigs.fasta"
    elif [[ -f "$dir/Assembly/contigs.fa" ]]; then
        fa="$dir/Assembly/contigs.fa"
    else
        echo "0	0	0	0"
        return
    fi

    python3 - "$fa" << 'PYEOF'
import sys

fasta = sys.argv[1]
lengths = []
current = 0
with open(fasta) as f:
    for line in f:
        if line.startswith(">"):
            if current > 0:
                lengths.append(current)
            current = 0
        else:
            current += len(line.strip())
    if current > 0:
        lengths.append(current)

lengths.sort(reverse=True)
total   = sum(lengths)
n_ctg   = len(lengths)
largest = lengths[0] if lengths else 0

cumsum = 0
n50 = 0
for l in lengths:
    cumsum += l
    if cumsum >= total / 2:
        n50 = l
        break

print(f"{total}\t{n_ctg}\t{largest}\t{n50}")
PYEOF
}

STATS1=$(get_assembly_stats "$SAMPLE1_DIR")
STATS2=$(get_assembly_stats "$SAMPLE2_DIR")

TOTAL1=$(echo "$STATS1" | cut -f1)
NCTG1=$(echo  "$STATS1" | cut -f2)
LARGE1=$(echo "$STATS1" | cut -f3)
N50_1=$(echo  "$STATS1" | cut -f4)

TOTAL2=$(echo "$STATS2" | cut -f1)
NCTG2=$(echo  "$STATS2" | cut -f2)
LARGE2=$(echo "$STATS2" | cut -f3)
N50_2=$(echo  "$STATS2" | cut -f4)

SIZE_DIFF=$(( TOTAL1 - TOTAL2 ))
SIZE_DIFF_ABS=${SIZE_DIFF#-}
CTGDIFF=$(( NCTG1 - NCTG2 ))
CTGDIFF_ABS=${CTGDIFF#-}

log "Assembly stats collected"

# =============================================================================
header "Step 3/3 — Prediction"
# =============================================================================

SCORE=0
NOTES=()

# SNP distance — most important criterion
if   [[ "$SNP_DISTANCE" -le 10 ]]; then
    SCORE=$(( SCORE + 3 ))
    SNP_VERDICT="SAME CLONE"
    SNP_COLOR="$GREEN"
    NOTES+=("SNP distance $SNP_DISTANCE ≤ 10 — strong evidence same clone")
elif [[ "$SNP_DISTANCE" -le 50 ]]; then
    SCORE=$(( SCORE + 1 ))
    SNP_VERDICT="POSSIBLY RELATED"
    SNP_COLOR="$YELLOW"
    NOTES+=("SNP distance $SNP_DISTANCE ≤ 50 — possibly related")
else
    SCORE=$(( SCORE - 2 ))
    SNP_VERDICT="DIFFERENT CLONES"
    SNP_COLOR="$RED"
    NOTES+=("SNP distance $SNP_DISTANCE > 50 — likely different clones")
fi

# Genome size difference
if   [[ "$SIZE_DIFF_ABS" -lt 50000 ]]; then
    SCORE=$(( SCORE + 2 ))
    SIZE_VERDICT="SIMILAR"
    SIZE_COLOR="$GREEN"
    NOTES+=("Genome size difference ${SIZE_DIFF_ABS} bp < 50kb — similar genomes")
elif [[ "$SIZE_DIFF_ABS" -lt 200000 ]]; then
    SIZE_VERDICT="MODERATE DIFFERENCE"
    SIZE_COLOR="$YELLOW"
    NOTES+=("Genome size difference ${SIZE_DIFF_ABS} bp — possible plasmid acquisition/loss")
else
    SCORE=$(( SCORE - 1 ))
    SIZE_VERDICT="LARGE DIFFERENCE"
    SIZE_COLOR="$RED"
    NOTES+=("Genome size difference ${SIZE_DIFF_ABS} bp > 200kb — significant extra genetic material")
fi

# Contig count difference
if [[ "$CTGDIFF_ABS" -lt 50 ]]; then
    SCORE=$(( SCORE + 1 ))
    CTG_VERDICT="SIMILAR"
    CTG_COLOR="$GREEN"
    NOTES+=("Contig count difference $CTGDIFF_ABS < 50 — similar assembly complexity")
else
    CTG_VERDICT="DIFFERENT"
    CTG_COLOR="$YELLOW"
    NOTES+=("Contig count difference $CTGDIFF_ABS — different assembly complexity")
fi

# Final verdict
if   [[ "$SCORE" -ge 5 ]]; then
    VERDICT="SAME CLONE"
    VERDICT_COLOR="$GREEN"
elif [[ "$SCORE" -ge 2 ]]; then
    VERDICT="LIKELY RELATED"
    VERDICT_COLOR="$YELLOW"
else
    VERDICT="DIFFERENT CLONES"
    VERDICT_COLOR="$RED"
fi

# ── Print table ───────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}╔══════════════════════════════════════════════════════════════╗${RESET}"
echo -e "${BOLD}║           ISOLATE COMPARISON REPORT                         ║${RESET}"
echo -e "${BOLD}╚══════════════════════════════════════════════════════════════╝${RESET}"
echo ""
printf "%-28s %-20s %-20s\n" "Metric" "$SAMPLE1" "$SAMPLE2"
echo "────────────────────────────────────────────────────────────────"
printf "%-28s %-20s %-20s\n" "Genome size (bp)"    "$TOTAL1"  "$TOTAL2"
printf "%-28s %-20s %-20s\n" "Contigs"             "$NCTG1"   "$NCTG2"
printf "%-28s %-20s %-20s\n" "Largest contig (bp)" "$LARGE1"  "$LARGE2"
printf "%-28s %-20s %-20s\n" "N50 (bp)"            "$N50_1"   "$N50_2"
echo "────────────────────────────────────────────────────────────────"
printf "%-28s %-20s\n" "SNP distance"      "$SNP_DISTANCE SNPs"
printf "%-28s %-20s\n" "Genome size diff"  "${SIZE_DIFF_ABS} bp"
printf "%-28s %-20s\n" "Contig count diff" "$CTGDIFF_ABS"
echo "────────────────────────────────────────────────────────────────"
echo ""
echo -e "  SNP verdict:    ${SNP_COLOR}${SNP_VERDICT}${RESET}"
echo -e "  Genome size:    ${SIZE_COLOR}${SIZE_VERDICT}${RESET}"
echo -e "  Contig count:   ${CTG_COLOR}${CTG_VERDICT}${RESET}"
echo ""
echo "────────────────────────────────────────────────────────────────"
echo -e "  PREDICTION:  ${VERDICT_COLOR}${BOLD}${VERDICT}${RESET}"
echo "────────────────────────────────────────────────────────────────"
echo ""
echo "Notes:"
for note in "${NOTES[@]}"; do
    echo "  • $note"
done
echo ""
warn "Always interpret genomic results with epidemiological context"
warn "(ward, collection dates, patient data)"

# ── Save TSV ──────────────────────────────────────────────────────────────────
{
echo -e "metric\t${SAMPLE1}\t${SAMPLE2}\tnotes"
echo -e "genome_size_bp\t${TOTAL1}\t${TOTAL2}\tdiff=${SIZE_DIFF_ABS}bp"
echo -e "contigs\t${NCTG1}\t${NCTG2}\tdiff=${CTGDIFF_ABS}"
echo -e "largest_contig_bp\t${LARGE1}\t${LARGE2}\t"
echo -e "N50_bp\t${N50_1}\t${N50_2}\t"
echo -e "SNP_distance\t${SNP_DISTANCE}\t-\t${SNP_VERDICT}"
echo -e "genome_size_verdict\t${SIZE_VERDICT}\t-\t"
echo -e "prediction\t${VERDICT}\t-\tscore=${SCORE}"
} > "$OUTPUT_TSV"

echo ""
success "Report saved to: $OUTPUT_TSV"
