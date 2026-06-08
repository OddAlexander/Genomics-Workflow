#!/usr/bin/env bash
# =============================================================================
# run-spatyper.sh
# Runs spaTyper — S. aureus protein A (spa) gene typing from WGS assembly
#
# Usage:
#   ./run-spatyper.sh <sample_dir>
#
# spa typing is ONLY meaningful for Staphylococcus aureus.
# Run this after confirming species identity with Kraken2 or FastANI.
#
# Output:
#   <sample_dir>/SpaType/
#     spatyper.tsv     Clean TSV with sample, spa type, repeat pattern
#     spatyper_raw.txt Full spaTyper stdout
#
# spa type interpretation:
#   t008 → CC8  USA300 (CA-MRSA)
#   t032 → CC22 EMRSA-15 (dominant European HA-MRSA)
#   t003 → CC8  Common HA-MRSA
#   t002 → CC5  Paediatric MRSA / MSSA
#   t041 → CC8  EMRSA-16
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
    echo "         S P A   T Y P I N G   ( S . a u r e u s )         "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir>"
    echo ""
    echo -e "${BOLD}Note:${RESET}"
    echo "  spa typing is only meaningful for Staphylococcus aureus."
    echo "  Confirm species with Kraken2 or FastANI first."
    echo ""
    echo -e "${BOLD}Common spa types:${RESET}"
    echo "  t008 → CC8  USA300 — dominant CA-MRSA"
    echo "  t032 → CC22 EMRSA-15 — dominant European HA-MRSA"
    echo "  t003 → CC8  Common hospital MRSA"
    echo "  t002 → CC5  Paediatric MRSA / MSSA"
    echo "  t041 → CC8  EMRSA-16 (UK)"
    echo "  t127 → CC1  Community MSSA/MRSA"
    echo ""
    echo -e "${BOLD}Output:${RESET}"
    echo "  <sample_dir>/SpaType/"
    echo "    spatyper.tsv      Clean TSV: sample, spa_type, repeats"
    echo "    spatyper_raw.txt  Full raw spaTyper output"
    echo ""
    echo -e "${BOLD}Example:${RESET}"
    echo "  $0 results/001s/"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"

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
OUTPUT_DIR="$SAMPLE_DIR/SpaType"
RAW_OUT="$OUTPUT_DIR/spatyper_raw.txt"
TSV_OUT="$OUTPUT_DIR/spatyper.tsv"

mkdir -p "$OUTPUT_DIR"

log "Sample:   $SAMPLE"
log "Assembly: $(basename "$ASSEMBLY_FA")"
log "Output:   $OUTPUT_DIR"

# =============================================================================
header "spa typing — Staphylococcus aureus protein A repeat region"
# =============================================================================

pixi run -- spaTyper -f "$ASSEMBLY_FA" > "$RAW_OUT" 2>&1

# =============================================================================
header "Results"
# =============================================================================

# Parse raw output — extract spa type and repeat pattern
python3 - "$RAW_OUT" "$TSV_OUT" "$SAMPLE" << 'PYEOF'
import sys

raw_file   = sys.argv[1]
tsv_file   = sys.argv[2]
sample     = sys.argv[3]

GREEN  = '\033[0;32m'
YELLOW = '\033[1;33m'
RED    = '\033[0;31m'
BOLD   = '\033[1m'
RESET  = '\033[0m'

# Known high-risk spa types and their clonal complexes
SPA_CC = {
    't008': ('CC8',  'USA300 — dominant CA-MRSA'),
    't032': ('CC22', 'EMRSA-15 — dominant European HA-MRSA'),
    't003': ('CC8',  'Common HA-MRSA'),
    't002': ('CC5',  'Paediatric MRSA / MSSA'),
    't041': ('CC8',  'EMRSA-16 UK lineage'),
    't127': ('CC1',  'Community MSSA/MRSA'),
    't037': ('CC239','Asian HA-MRSA'),
    't045': ('CC45', 'Community MSSA, occasional MRSA'),
    't044': ('CC8',  'USA300 variant'),
    't026': ('CC30', 'SW-MRSA — Southwest Pacific MRSA'),
}

results = []
with open(raw_file) as f:
    for line in f:
        line = line.strip()
        # Tab-separated results line (skip header and info lines)
        cols = line.split("\t")
        if len(cols) == 3 and not line.startswith("Sequence name"):
            results.append({
                'contig':   cols[0].strip(),
                'repeats':  cols[1].strip(),
                'spa_type': cols[2].strip(),
            })

# Write TSV
with open(tsv_file, 'w') as out:
    out.write('sample\tcontig\tspa_type\trepeats\tclonal_complex\tnotes\n')
    for r in results:
        spa  = r['spa_type']
        cc, note = SPA_CC.get(spa, ('-', '-'))
        out.write(f"{sample}\t{r['contig']}\t{spa}\t{r['repeats']}\t{cc}\t{note}\n")

# Print summary
if not results:
    print(f"  {YELLOW}No spa type detected — spa gene may be absent, fragmented,")
    print(f"  or this is not S. aureus.{RESET}")
else:
    # Group by spa type (multiple contigs may give same result)
    seen = {}
    for r in results:
        spa = r['spa_type']
        if spa not in seen:
            seen[spa] = r

    for spa, r in seen.items():
        cc, note = SPA_CC.get(spa, ('-', 'Not in common type database'))
        if spa.startswith('t') and spa != '-':
            color = GREEN
        else:
            color = YELLOW

        print(f"  {BOLD}spa type:{RESET}  {color}{spa}{RESET}")
        print(f"  Repeats:   {r['repeats']}")
        print(f"  CC:        {cc}")
        print(f"  Notes:     {note}")
        print()
PYEOF

success "spa typing done"
echo ""
echo "────────────────────────────────────────────────────────"
echo "  Clean TSV:  $TSV_OUT"
echo "  Raw output: $RAW_OUT"
echo ""
warn "spa typing is only valid for S. aureus — verify species first"
echo "────────────────────────────────────────────────────────"
