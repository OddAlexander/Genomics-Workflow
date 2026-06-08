#!/usr/bin/env bash
# =============================================================================
# pixi_mlst.sh
# Runs MLST — multilocus sequence typing
#
# Usage:
#   ./pixi_mlst.sh <sample_dir> [scheme]
#
# Arguments:
#   sample_dir  — sample folder containing an Assembly/ subfolder
#   scheme      — (optional) force a specific MLST scheme
#                 Default: auto-detect
#
# Output:
#   <sample_dir>/MLST/mlst.tsv   (clean TSV with header)
#
# Optional env vars:
#   THREADS — number of threads (default: 8)
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

if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "              M L S T   T Y P I N G                        "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> [scheme]"
    echo ""
    echo -e "${BOLD}Common schemes:${RESET}"
    echo "  efaecium         Enterococcus faecium"
    echo "  efaecalis        Enterococcus faecalis"
    echo "  ecoli            Escherichia coli"
    echo "  ecoli_achtman_4  Escherichia coli (Achtman 4)"
    echo "  klebsiella       Klebsiella pneumoniae"
    echo "  koxytoca         Klebsiella oxytoca"
    echo "  saureus          Staphylococcus aureus"
    echo "  spneumoniae      Streptococcus pneumoniae"
    echo "  spyogenes        Streptococcus pyogenes"
    echo "  salmonella       Salmonella enterica"
    echo "  campylobacter    Campylobacter jejuni/coli"
    echo "  abaumannii       Acinetobacter baumannii"
    echo "  paeruginosa      Pseudomonas aeruginosa"
    echo "  neisseria        Neisseria spp."
    echo "  hinfluenzae      Haemophilus influenzae"
    echo "  cdifficile       Clostridioides difficile"
    echo "  ecloacae         Enterobacter cloacae"
    echo "  smaltophilia     Stenotrophomonas maltophilia"
    echo ""
    echo -e "${BOLD}List ALL schemes:${RESET}"
    echo "  pixi run -- bash -c \"mlst --list\""
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/"
    echo "  $0 results/001k/ efaecium"
    echo "============================================================"
    exit 0
fi

SAMPLE_DIR="${1%/}"
SCHEME="${2:-}"
THREADS="${THREADS:-8}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

ASSEMBLY_DIR="$SAMPLE_DIR/Assembly"
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly directory not found: $ASSEMBLY_DIR"

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
[[ -n "$SCHEME" ]] && log "Scheme:   $SCHEME" || log "Scheme:   auto-detect"

OUTPUT_DIR="$SAMPLE_DIR/MLST"
mkdir -p "$OUTPUT_DIR"
OUTPUT_TSV="$OUTPUT_DIR/mlst.tsv"
RAW="$OUTPUT_DIR/.mlst_raw.tsv"
log "Output:   $OUTPUT_TSV"

# =============================================================================
header "MLST — multilocus sequence typing"
# =============================================================================

if [[ -n "$SCHEME" ]]; then
    pixi run -- bash -c "mlst --scheme '$SCHEME' '$ASSEMBLY_FA'" > "$RAW"
else
    pixi run -- bash -c "mlst '$ASSEMBLY_FA'" > "$RAW"
fi

# ── Reformat into clean TSV with proper header ────────────────────────────────
python3 << PYEOF
import sys

raw_path    = "$RAW"
out_path    = "$OUTPUT_TSV"
sample_name = "$SAMPLE"

with open(raw_path) as f:
    line = f.readline().strip()

cols = line.split('\t')
scheme  = cols[1] if len(cols) > 1 else '-'
st      = cols[2] if len(cols) > 2 else '-'
alleles = cols[3:] if len(cols) > 3 else []

# Parse locus(allele) into separate columns
loci = []
for a in alleles:
    if '(' in a and a.endswith(')'):
        locus  = a[:a.index('(')]
        allele = a[a.index('(')+1:-1]
        loci.append((locus, allele))
    else:
        loci.append((a, '-'))

with open(out_path, 'w') as out:
    headers = ['sample', 'scheme', 'ST'] + [l for l, _ in loci]
    out.write('\t'.join(headers) + '\n')
    values  = [sample_name, scheme, st] + [a for _, a in loci]
    out.write('\t'.join(values) + '\n')
PYEOF

rm -f "$RAW"

# ── Parse for display ─────────────────────────────────────────────────────────
SCHEME_FOUND=$(awk -F'\t' 'NR==2{print $2}' "$OUTPUT_TSV")
ST=$(awk -F'\t' 'NR==2{print $3}' "$OUTPUT_TSV")

success "MLST done"
echo ""
echo "────────────────────────────────────────"
echo ""

# Print as clean aligned table
column -t -s $'\t' "$OUTPUT_TSV"
echo ""

# ST interpretation
if [[ "$ST" == "-" ]]; then
    warn "No ST assigned — novel allele combination or wrong scheme"
    echo "  Tip: check https://pubmlst.org or try a different scheme"
elif [[ "$ST" == "~"* ]]; then
    warn "Closest ST: $ST — inexact match, possible novel allele"
else
    echo -e "  Result: ${GREEN}${SCHEME_FOUND} ST${ST}${RESET}"
    echo ""

    case "$SCHEME_FOUND:$ST" in
        efaecium:17)         echo -e "  ${YELLOW}⚠ ST17 — high-risk hospital-adapted VRE lineage${RESET}" ;;
        efaecium:18)         echo -e "  ${YELLOW}⚠ ST18 — high-risk hospital-adapted VRE lineage${RESET}" ;;
        efaecium:80)         echo -e "  ${YELLOW}⚠ ST80 — high-risk hospital-adapted VRE lineage${RESET}" ;;
        efaecium:117)        echo -e "  ${YELLOW}⚠ ST117 — high-risk hospital-adapted VRE lineage${RESET}" ;;
        klebsiella:258)      echo -e "  ${YELLOW}⚠ ST258 — dominant KPC-producing Klebsiella lineage${RESET}" ;;
        klebsiella:11)       echo -e "  ${YELLOW}⚠ ST11 — high-risk carbapenem-resistant Klebsiella${RESET}" ;;
        klebsiella:14)       echo -e "  ${YELLOW}⚠ ST14 — associated with NDM-producing Klebsiella${RESET}" ;;
        ecoli:131)           echo -e "  ${YELLOW}⚠ ST131 — dominant MDR E. coli, ESBL/FQ-resistant${RESET}" ;;
        ecoli_achtman_4:131) echo -e "  ${YELLOW}⚠ ST131 — dominant MDR E. coli, ESBL/FQ-resistant${RESET}" ;;
        ecoli:648)           echo -e "  ${YELLOW}⚠ ST648 — high-risk MDR E. coli lineage${RESET}" ;;       
        abaumannii:2)        echo -e "  ${YELLOW}⚠ ST2 — dominant carbapenem-resistant Acinetobacter${RESET}" ;;
        paeruginosa:235)     echo -e "  ${YELLOW}⚠ ST235 — high-risk MDR Pseudomonas lineage${RESET}" ;;
        
    esac
fi

echo ""
echo "────────────────────────────────────────"
log "Full results: $OUTPUT_TSV"
