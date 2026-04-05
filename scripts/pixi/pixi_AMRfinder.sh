#!/usr/bin/env bash
# =============================================================================
# pixi_AMRfinder.sh.sh
# Runs AMRFinder v4 (AMR + virulence) and Abricate PlasmidFinder
# Generates a summary report linking plasmid replicons to their gene cargo
#
# Usage:
#   ./pixi_AMRfinder.sh <sample_dir> [organism]
#
# Arguments:
#   sample_dir  — sample folder containing an Assembly/ subfolder
#   organism    — (optional) organism for point mutation screening
#                 Run with no args to see full list
#
# Output:
#   <sample_dir>/AMR/
#     amrfinder_all.tsv         All AMR + virulence + stress hits
#     amrfinder_amr.tsv         AMR genes only
#     amrfinder_virulence.tsv   Virulence genes only
#     amrfinder_stress.tsv      Stress response genes only
#     plasmidfinder.tsv         Plasmid replicons
#     plasmid_cargo.tsv         Genes found on same contig as each plasmid
#     report.txt                Human-readable summary report
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

# ── No arguments — show help ──────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "     A M R   +   V I R U L E N C E   +   P L A S M I D     "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> [organism]"
    echo ""
    echo -e "${BOLD}Organism options (enables point mutation screening):${RESET}"
    echo "  Acinetobacter_baumannii    Burkholderia_cepacia"
    echo "  Burkholderia_pseudomallei  Campylobacter"
    echo "  Clostridioides_difficile   Enterococcus_faecalis"
    echo "  Enterococcus_faecium       Escherichia"
    echo "  Klebsiella_oxytoca         Klebsiella_pneumoniae"
    echo "  Neisseria_gonorrhoeae      Neisseria_meningitidis"
    echo "  Pseudomonas_aeruginosa     Salmonella"
    echo "  Serratia_marcescens        Staphylococcus_aureus"
    echo "  Staphylococcus_pseudintermedius"
    echo "  Streptococcus_agalactiae   Streptococcus_pneumoniae"
    echo "  Streptococcus_pyogenes     Vibrio_cholerae"
    echo "  Vibrio_parahaemolyticus    Vibrio_vulnificus"
    echo ""
    echo -e "${BOLD}Optional env vars:${RESET}"
    echo "  THREADS — number of threads (default: 8)"
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/"
    echo "  $0 results/001k/ Klebsiella_pneumoniae"
    echo "  $0 results/001k/ Escherichia"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
ORGANISM="${2:-}"
THREADS="${THREADS:-8}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"

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
    die "No assembly file found in $ASSEMBLY_DIR"
fi

SAMPLE=$(basename "$SAMPLE_DIR")
log "Sample:   $SAMPLE"
log "Assembly: $(basename "$ASSEMBLY_FA")"
if [[ -n "$ORGANISM" ]]; then
    log "Organism: $ORGANISM (point mutation screening enabled)"
else
    warn "No organism set — point mutation screening disabled"
fi

# ── Output folder ─────────────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/AMR"
mkdir -p "$OUTPUT_DIR"
log "Output:   $OUTPUT_DIR"

STEP_START=0
start_timer() { STEP_START=$SECONDS; }
elapsed()     { echo $(( SECONDS - STEP_START ))s; }

# =============================================================================
header "Step 1/3 — AMRFinder v4: AMR + virulence + stress genes"
# =============================================================================
start_timer

AMR_ALL="$OUTPUT_DIR/amrfinder_all.tsv"
ORG_FLAG=""
[[ -n "$ORGANISM" ]] && ORG_FLAG="--organism $ORGANISM"

pixi run --environment amrfinder4 -- amrfinder \
  --nucleotide "$ASSEMBLY_FA" \
  --output "$AMR_ALL" \
  --threads "$THREADS" \
  --plus \
  $ORG_FLAG

HEAD=$(head -1 "$AMR_ALL")

for TYPE in AMR VIRULENCE STRESS; do
    {
        echo "$HEAD"
        awk -F'\t' -v t="$TYPE" 'NR>1 && $9==t' "$AMR_ALL"
    } > "$OUTPUT_DIR/amrfinder_${TYPE,,}.tsv"
done

AMR_COUNT=$(( $(wc -l < "$OUTPUT_DIR/amrfinder_amr.tsv") - 1 ))
VIR_COUNT=$(( $(wc -l < "$OUTPUT_DIR/amrfinder_virulence.tsv") - 1 ))
STR_COUNT=$(( $(wc -l < "$OUTPUT_DIR/amrfinder_stress.tsv") - 1 ))

success "AMRFinder done ($(elapsed)) — AMR: $AMR_COUNT  Virulence: $VIR_COUNT  Stress: $STR_COUNT"

# =============================================================================
header "Step 2/3 — Abricate PlasmidFinder: plasmid replicons"
# =============================================================================
start_timer

PLASMID_TSV="$OUTPUT_DIR/plasmidfinder.tsv"

pixi run -- abricate \
  --db plasmidfinder \
  --threads "$THREADS" \
  "$ASSEMBLY_FA" \
  > "$PLASMID_TSV"

PLASMID_COUNT=$(( $(wc -l < "$PLASMID_TSV") - 1 ))

success "PlasmidFinder done ($(elapsed)) — $PLASMID_COUNT replicon(s) found"

# =============================================================================
header "Step 3/3 — Linking plasmid replicons to gene cargo"
# =============================================================================

CARGO_TSV="$OUTPUT_DIR/plasmid_cargo.tsv"

python3 - "$PLASMID_TSV" "$AMR_ALL" "$CARGO_TSV" << 'PYEOF'
import sys

plasmid_file = sys.argv[1]
amr_file     = sys.argv[2]
cargo_file   = sys.argv[3]

# Load all AMR/virulence/stress genes with their contig
# AMRFinder col 2 = contig/sequence name, col 6 = gene, col 9 = type, col 12 = subclass
amr_by_contig = {}
with open(amr_file) as f:
    next(f)  # skip header
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 12:
            continue
        contig   = cols[1]
        gene     = cols[5]
        etype    = cols[8]
        subclass = cols[11] if len(cols) > 11 else ''
        if contig not in amr_by_contig:
            amr_by_contig[contig] = []
        amr_by_contig[contig].append((gene, etype, subclass))

# Load plasmid replicons
# Abricate col 2 = contig (SEQUENCE), col 6 = gene (replicon name)
# col 9 = %COVERAGE, col 10 = %IDENTITY
rows = []
with open(plasmid_file) as f:
    next(f)  # skip header
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 10:
            continue
        contig   = cols[1]
        replicon = cols[5]
        coverage = cols[8]
        identity = cols[9]

        # Find genes on same contig
        cargo = amr_by_contig.get(contig, [])
        if cargo:
            cargo_str = '; '.join(
                f"{g} [{t}]" + (f" ({s})" if s else '')
                for g, t, s in cargo
            )
        else:
            cargo_str = 'none detected'

        rows.append((replicon, contig, coverage, identity, cargo_str))

with open(cargo_file, 'w') as out:
    out.write('replicon\tcontig\tcoverage\tidentity\tgenes_on_contig\n')
    for row in rows:
        out.write('\t'.join(row) + '\n')

print(f"Linked {len(rows)} plasmid replicon(s) to gene cargo")
PYEOF

success "Plasmid cargo analysis done"

# =============================================================================
header "Generating report"
# =============================================================================

REPORT="$OUTPUT_DIR/report.txt"

{
echo "======================================================"
echo "  AMR / VIRULENCE / PLASMID REPORT"
echo "======================================================"
echo "Sample:    $SAMPLE"
echo "Assembly:  $(basename "$ASSEMBLY_FA")"
[[ -n "$ORGANISM" ]] && echo "Organism:  $ORGANISM" || echo "Organism:  not specified"
echo "Date:      $(date)"
echo ""

# ── AMR genes ──────────────────────────────────────────────────────────────
echo "------------------------------------------------------"
echo "AMR GENES ($AMR_COUNT)"
echo "------------------------------------------------------"
if [[ "$AMR_COUNT" -gt 0 ]]; then
    awk -F'\t' 'NR>1{
        printf "  %-22s %-18s %s%%id  %s%%cov\n", $6, $12, $16, $17
    }' "$OUTPUT_DIR/amrfinder_amr.tsv"
else
    echo "  None detected"
fi
echo ""

# ── Virulence genes ────────────────────────────────────────────────────────
echo "------------------------------------------------------"
echo "VIRULENCE GENES ($VIR_COUNT)"
echo "------------------------------------------------------"
if [[ "$VIR_COUNT" -gt 0 ]]; then
    awk -F'\t' 'NR>1{
        printf "  %-22s %s\n", $6, $7
    }' "$OUTPUT_DIR/amrfinder_virulence.tsv"
else
    echo "  None detected"
fi
echo ""

# ── Stress genes ───────────────────────────────────────────────────────────
echo "------------------------------------------------------"
echo "STRESS RESPONSE GENES ($STR_COUNT)"
echo "------------------------------------------------------"
if [[ "$STR_COUNT" -gt 0 ]]; then
    awk -F'\t' 'NR>1{
        printf "  %-22s %s\n", $6, $7
    }' "$OUTPUT_DIR/amrfinder_stress.tsv"
else
    echo "  None detected"
fi
echo ""

# ── Plasmids + cargo ───────────────────────────────────────────────────────
echo "------------------------------------------------------"
echo "PLASMID REPLICONS ($PLASMID_COUNT)"
echo "------------------------------------------------------"
if [[ "$PLASMID_COUNT" -gt 0 ]]; then
    awk -F'\t' 'NR>1{
        printf "  Replicon:  %s\n", $1
        printf "  Contig:    %s\n", $2
        printf "  Coverage:  %s%%   Identity: %s%%\n", $3, $4
        printf "  Cargo:     %s\n", $5
        printf "\n"
    }' "$OUTPUT_DIR/plasmid_cargo.tsv"
else
    echo "  None detected"
fi

echo "------------------------------------------------------"
echo "OUTPUT FILES"
echo "------------------------------------------------------"
echo "  amrfinder_all.tsv       All raw AMRFinder results"
echo "  amrfinder_amr.tsv       AMR genes"
echo "  amrfinder_virulence.tsv Virulence genes"
echo "  amrfinder_stress.tsv    Stress response genes"
echo "  plasmidfinder.tsv       Plasmid replicons (Abricate)"
echo "  plasmid_cargo.tsv       Plasmid-to-gene linkage"
echo "  report.txt              This report"
echo "======================================================"
} > "$REPORT"

# Print report to terminal
cat "$REPORT"

echo ""
log "Report saved to: $REPORT"
