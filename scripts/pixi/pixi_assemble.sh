#!/usr/bin/env bash
# =============================================================================
# Assembly pipeline
# Runs: fastp → FastQC → SPAdes
#
# Usage:
#   ./pixi_assemble.sh <input_dir> <output_base_dir>
#
# Supported R1/R2 filename patterns:
#   sampleA_R1.fastq.gz     / sampleA_R2.fastq.gz
#   sampleA_R1_001.fastq.gz / sampleA_R2_001.fastq.gz
#   sampleA_1.fastq.gz      / sampleA_2.fastq.gz
#   ERR123456_1.fastq.gz    / ERR123456_2.fastq.gz
#
# Optional env vars:
#   THREADS — number of threads (default: 8)
# =============================================================================

set -euo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()     { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
success() { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${RESET}"; }
die()     { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ $*${RESET}" >&2; exit 1; }
header()  { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${RESET}"
            echo -e "${BOLD}${CYAN}  $*${RESET}"
            echo -e "${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; }

# ── Arguments ─────────────────────────────────────────────────────────────────
[[ $# -lt 2 ]] && die "Usage: $0 <input_dir> <output_base_dir>"

INPUT_DIR="${1%/}"
OUTPUT_BASE="${2%/}"
THREADS="${THREADS:-8}"

[[ -d "$INPUT_DIR" ]] || die "Input directory not found: $INPUT_DIR"

# ── Auto-detect R1/R2 ─────────────────────────────────────────────────────────
log "Scanning $INPUT_DIR for paired fastq.gz files..."

R1_FILE=""
R2_FILE=""

R1_FILE=$(find "$INPUT_DIR" -maxdepth 1 -name "*_R1*.fastq.gz" | sort | head -1)
if [[ -n "$R1_FILE" ]]; then
    R2_FILE="${R1_FILE/_R1/_R2}"
fi

if [[ -z "$R1_FILE" ]] || [[ ! -f "$R2_FILE" ]]; then
    R1_FILE=$(find "$INPUT_DIR" -maxdepth 1 -name "*_1.fastq.gz" | sort | head -1)
    if [[ -n "$R1_FILE" ]]; then
        R2_FILE="${R1_FILE/_1.fastq.gz/_2.fastq.gz}"
    fi
fi

[[ -n "$R1_FILE" ]] || die "No R1 fastq.gz found in $INPUT_DIR\nExpected: *_R1*.fastq.gz or *_1.fastq.gz"
[[ -f "$R1_FILE" ]] || die "R1 not found: $R1_FILE"
[[ -f "$R2_FILE" ]] || die "R2 not found: $R2_FILE"

log "Found R1: $(basename "$R1_FILE")"
log "Found R2: $(basename "$R2_FILE")"

# ── Sample name from filename ─────────────────────────────────────────────────
R1_BASENAME=$(basename "$R1_FILE")
SAMPLE="${R1_BASENAME%%_R1*}"
if [[ "$SAMPLE" == "$R1_BASENAME" ]]; then
    SAMPLE="${R1_BASENAME%%_1.fastq*}"
fi
log "Sample name: $SAMPLE"

# ── Output folders ────────────────────────────────────────────────────────────
OUTPUT_DIR="$OUTPUT_BASE/$SAMPLE"
TRIM_DIR="$OUTPUT_DIR/Trimmed"
FASTQC_DIR="$OUTPUT_DIR/QC/FastQC"
MULTIQC_DIR="$OUTPUT_DIR/QC/MultiQC"
ASSEMBLY_DIR="$OUTPUT_DIR/Assembly"

mkdir -p "$TRIM_DIR" "$FASTQC_DIR" "$MULTIQC_DIR" "$ASSEMBLY_DIR"
log "Output directory: $OUTPUT_DIR"

# ── Timing ────────────────────────────────────────────────────────────────────
STEP_START=0
start_timer() { STEP_START=$SECONDS; }
elapsed()     { echo $(( SECONDS - STEP_START ))s; }

# ── Summary ───────────────────────────────────────────────────────────────────
SUMMARY="$OUTPUT_DIR/pipeline_summary.txt"
{
echo "Pipeline Summary"
echo "================"
echo "Sample:     $SAMPLE"
echo "R1:         $(basename "$R1_FILE")"
echo "R2:         $(basename "$R2_FILE")"
echo "Output:     $OUTPUT_DIR"
echo "Threads:    $THREADS"
echo "Started:    $(date)"
echo ""
} > "$SUMMARY"

# =============================================================================
header "Step 1/3 — fastp: adapter trimming"
# =============================================================================
start_timer

pixi run -- fastp \
  --thread "$THREADS" \
  -i "$R1_FILE" \
  -I "$R2_FILE" \
  -o "$TRIM_DIR/${SAMPLE}_R1.fastq.gz" \
  -O "$TRIM_DIR/${SAMPLE}_R2.fastq.gz" \
  -h "$TRIM_DIR/${SAMPLE}_fastp.html"

success "fastp done ($(elapsed))"
echo "fastp:          $(elapsed)" >> "$SUMMARY"

# =============================================================================
header "Step 2/3 — FastQC: quality reports on trimmed reads"
# =============================================================================
start_timer

pixi run -- fastqc \
  --threads "$THREADS" \
  "$TRIM_DIR/${SAMPLE}_R1.fastq.gz" \
  "$TRIM_DIR/${SAMPLE}_R2.fastq.gz" \
  -o "$FASTQC_DIR"

success "FastQC done ($(elapsed))"
echo "FastQC:         $(elapsed)" >> "$SUMMARY"

# =============================================================================
header "Step 3/3 — SPAdes: genome assembly"
# =============================================================================
start_timer

pixi run -- spades.py \
  --pe1-1 "$TRIM_DIR/${SAMPLE}_R1.fastq.gz" \
  --pe1-2 "$TRIM_DIR/${SAMPLE}_R2.fastq.gz" \
  -o "$ASSEMBLY_DIR" \
  --threads "$THREADS"

ASSEMBLY_FA="$ASSEMBLY_DIR/scaffolds.fasta"
[[ -f "$ASSEMBLY_FA" ]] || die "SPAdes failed — scaffolds.fasta not found"

CONTIGS=$(grep -c "^>" "$ASSEMBLY_FA" || true)
TOTAL_BP=$(grep -v "^>" "$ASSEMBLY_FA" | tr -d '\n' | wc -c || true)

success "SPAdes done ($(elapsed)) — $CONTIGS contigs, ~${TOTAL_BP} bp"
echo "SPAdes:         $(elapsed) — $CONTIGS contigs, ${TOTAL_BP} bp" >> "$SUMMARY"

# =============================================================================
header "Pipeline complete"
# =============================================================================
{
echo ""
echo "Finished: $(date)"
echo ""
echo "Output structure:"
echo "  $OUTPUT_DIR/"
echo "  ├── Trimmed/"
echo "  │   ├── ${SAMPLE}_R1.fastq.gz"
echo "  │   ├── ${SAMPLE}_R2.fastq.gz"
echo "  │   └── ${SAMPLE}_fastp.html"
echo "  ├── QC/"
echo "  │   ├── FastQC/"
echo "  │   └── MultiQC/multiqc_report.html"
echo "  ├── Assembly/"
echo "  │   ├── scaffolds.fasta     ($CONTIGS contigs, ~${TOTAL_BP} bp)"
echo "  │   ├── contigs.fasta"
echo "  │   └── assembly_graph.fastg"
echo "  └── pipeline_summary.txt"
} | tee -a "$SUMMARY"

log "Summary saved to: $SUMMARY"
