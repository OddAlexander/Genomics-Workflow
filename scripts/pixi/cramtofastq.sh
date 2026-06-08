#!/usr/bin/env bash
# ─────────────────────────────────────────────
#  cram_to_fastq.sh
#  Recursively find CRAM files and convert to
#  paired-end FASTQ.gz using samtools via pixi.
# ─────────────────────────────────────────────
 
set -euo pipefail
 
# ── Colours ───────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'
 
# ── Usage ─────────────────────────────────────
usage() {
  echo -e "Usage: $0 <input_folder> [--reference <ref.fa>]"
  echo -e "  <input_folder>        Root folder to search for .cram files"
  echo -e "  --reference <ref.fa>  Optional reference FASTA (recommended)"
  exit 1
}
 
[[ $# -lt 1 ]] && usage
 
INPUT_DIR="$1"; shift
REFERENCE=""
 
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference|-r) REFERENCE="$2"; shift 2 ;;
    *) echo -e "${RED}Unknown argument: $1${RESET}"; usage ;;
  esac
done
 
if [[ ! -d "$INPUT_DIR" ]]; then
  echo -e "${RED}Error: '$INPUT_DIR' is not a directory.${RESET}"
  exit 1
fi
 
# ── Check pixi + samtools ──────────────────────
if ! command -v pixi &>/dev/null; then
  echo -e "${RED}Error: pixi not found. Install from https://pixi.sh${RESET}"
  exit 1
fi
 
# ── Collect CRAM files ────────────────────────
mapfile -t CRAM_FILES < <(find "$INPUT_DIR" -type f -name "*.cram" | sort)
 
TOTAL=${#CRAM_FILES[@]}
if [[ $TOTAL -eq 0 ]]; then
  echo -e "${YELLOW}No .cram files found in '$INPUT_DIR'.${RESET}"
  exit 0
fi
 
echo -e "${BOLD}${CYAN}Found $TOTAL CRAM file(s) to process.${RESET}\n"
 
# ── Report arrays ─────────────────────────────
declare -a SUCCESS_LIST=()
declare -a SKIPPED_LIST=()
declare -a ERROR_LIST=()
declare -A ERROR_MSG=()
 
# ── Process each CRAM ─────────────────────────
for CRAM in "${CRAM_FILES[@]}"; do
  DIR=$(dirname "$CRAM")
  BASE=$(basename "$CRAM" .cram)
  R1="${DIR}/${BASE}_R1.fastq.gz"
  R2="${DIR}/${BASE}_R2.fastq.gz"
 
  echo -e "${BOLD}Processing:${RESET} $CRAM"
 
  # Skip if outputs already exist
  if [[ -f "$R1" && -f "$R2" ]]; then
    echo -e "  ${YELLOW}⚠ Skipping — output files already exist.${RESET}"
    SKIPPED_LIST+=("$CRAM")
    continue
  fi
 
  # Build samtools command
  REF_FLAG=""
  [[ -n "$REFERENCE" ]] && REF_FLAG="--reference $REFERENCE"
 
  CMD="samtools collate -u -O --no-PG \"$CRAM\" | \
    samtools fastq $REF_FLAG \
      -1 \"$R1\" \
      -2 \"$R2\" \
      -0 /dev/null \
      -s /dev/null \
      -n"
 
  # Run via pixi and capture stderr
  STDERR_TMP=$(mktemp)
  EXIT_CODE=0
  pixi run bash -c "$CMD" 2>"$STDERR_TMP" || EXIT_CODE=$?
 
  STDERR_CONTENT=$(cat "$STDERR_TMP")
  rm -f "$STDERR_TMP"
 
  if [[ $EXIT_CODE -ne 0 ]]; then
    # Detect common error types
    if echo "$STDERR_CONTENT" | grep -qi "reference"; then
      ERR_TYPE="Missing or mismatched reference genome"
    elif echo "$STDERR_CONTENT" | grep -qi "truncat"; then
      ERR_TYPE="Truncated or corrupt CRAM file"
    elif echo "$STDERR_CONTENT" | grep -qi "permission"; then
      ERR_TYPE="Permission denied"
    else
      ERR_TYPE="Unknown error"
    fi
 
    ERROR_LIST+=("$CRAM")
    ERROR_MSG["$CRAM"]="${ERR_TYPE}: ${STDERR_CONTENT:0:300}"
 
    echo -e "  ${RED}✗ Failed — ${ERR_TYPE}${RESET}"
    echo -e "  ${RED}$(echo "$STDERR_CONTENT" | head -3)${RESET}"
 
    # Clean up partial outputs
    [[ -f "$R1" ]] && rm -f "$R1"
    [[ -f "$R2" ]] && rm -f "$R2"
  else
    SUCCESS_LIST+=("$CRAM")
    echo -e "  ${GREEN}✓ Done → $(basename "$R1"), $(basename "$R2")${RESET}"
  fi
 
  echo ""
done
 
# ── Final Report ──────────────────────────────
REPORT_FILE="${INPUT_DIR}/cram_conversion_report.txt"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
 
{
  echo "═══════════════════════════════════════════════════"
  echo "  CRAM → FASTQ Conversion Report"
  echo "  Generated: $TIMESTAMP"
  echo "═══════════════════════════════════════════════════"
  echo ""
  echo "Input folder : $INPUT_DIR"
  [[ -n "$REFERENCE" ]] && echo "Reference    : $REFERENCE" || echo "Reference    : None provided"
  echo ""
  echo "Summary"
  echo "───────────────────────────────────────────────────"
  echo "  Total found  : $TOTAL"
  echo "  Succeeded    : ${#SUCCESS_LIST[@]}"
  echo "  Skipped      : ${#SKIPPED_LIST[@]}"
  echo "  Failed       : ${#ERROR_LIST[@]}"
  echo ""
 
  if [[ ${#SUCCESS_LIST[@]} -gt 0 ]]; then
    echo "✓ Successful conversions"
    echo "───────────────────────────────────────────────────"
    for F in "${SUCCESS_LIST[@]}"; do
      DIR=$(dirname "$F"); BASE=$(basename "$F" .cram)
      echo "  $F"
      echo "    → ${DIR}/${BASE}_R1.fastq.gz"
      echo "    → ${DIR}/${BASE}_R2.fastq.gz"
    done
    echo ""
  fi
 
  if [[ ${#SKIPPED_LIST[@]} -gt 0 ]]; then
    echo "⚠ Skipped (output already exists)"
    echo "───────────────────────────────────────────────────"
    for F in "${SKIPPED_LIST[@]}"; do echo "  $F"; done
    echo ""
  fi
 
  if [[ ${#ERROR_LIST[@]} -gt 0 ]]; then
    echo "✗ Errors"
    echo "───────────────────────────────────────────────────"
    for F in "${ERROR_LIST[@]}"; do
      echo "  $F"
      echo "    Error: ${ERROR_MSG[$F]}"
      echo ""
    done
 
    echo "Troubleshooting tips"
    echo "───────────────────────────────────────────────────"
    echo "  Missing reference:"
    echo "    Re-run with: --reference /path/to/reference.fa"
    echo "    Check CRAM header: samtools view -H <file.cram> | grep '^@SQ'"
    echo "    Look for 'UR:' tag to find the original reference URL."
    echo ""
    echo "  Truncated CRAM:"
    echo "    Validate with: samtools quickcheck <file.cram>"
    echo "    Try re-downloading the file if possible."
  fi
 
  echo "═══════════════════════════════════════════════════"
} | tee "$REPORT_FILE"
 
echo -e "\n${BOLD}Report saved to:${RESET} $REPORT_FILE"
