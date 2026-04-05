#!/usr/bin/env bash
# =============================================================================
# pixi_ID_fastANI.sh
# Runs FastANI to compute ANI between your assembly and reference genomes
#
# Usage:
#   ./pixi_ID_fastANI.sh <sample_dir> [references]
#
# Arguments:
#   sample_dir   — sample folder containing an Assembly/ subfolder
#   references   — (optional) either:
#                    a) path to a directory of reference .fasta/.fna/.fa files
#                    b) path to a reference list text file (one path per line)
#                  Default: /databases/fastANI_refs/reference_list.txt
#
# Species names are resolved from:
#   /databases/fastANI_refs/ncbi_dataset/data/assembly_data_report.jsonl
#
# Output:
#   <sample_dir>/ID_FastANI/fastani.tsv   (sorted by ANI, with species names)
#
# Optional env vars:
#   THREADS        — number of threads (default: 8)
#   ANI_THRESHOLD  — minimum ANI % to report (default: 80)
#   JSONL_DB       — path to assembly_data_report.jsonl
#                    Default: /databases/fastANI_refs/ncbi_dataset/data/assembly_data_report.jsonl
# =============================================================================

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
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
    echo "        F A S T A N I   I D E N T I F I C A T I O N        "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> [references]"
    echo ""
    echo -e "${BOLD}Arguments:${RESET}"
    echo "  sample_dir   Sample folder containing Assembly/ subfolder"
    echo "  references   (optional) Directory of reference genomes or"
    echo "               a text file with one genome path per line"
    echo "               Default: /databases/fastANI_refs/reference_list.txt"
    echo ""
    echo -e "${BOLD}Optional env vars:${RESET}"
    echo "  THREADS        Number of threads (default: 8)"
    echo "  ANI_THRESHOLD  Minimum ANI % to report (default: 80)"
    echo ""
    echo -e "${BOLD}Interpretation:${RESET}"
    echo "  ANI >= 95%   Same species"
    echo "  ANI 85-95%   Same genus, different species"
    echo "  ANI < 85%    Different genus"
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/"
    echo "  $0 results/001k/ /databases/fastANI_refs/enterobacteriaceae/"
    echo "  $0 results/001k/ /databases/fastANI_refs/reference_list.txt"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
REFERENCES="${2:-/databases/fastANI_refs/reference_list.txt}"
THREADS="${THREADS:-8}"
ANI_THRESHOLD="${ANI_THRESHOLD:-80}"
JSONL_DB="${JSONL_DB:-/databases/fastANI_refs/ncbi_dataset/data/assembly_data_report.jsonl}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -e "$REFERENCES" ]] || die "References not found: $REFERENCES"

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

log "Query:     $(basename "$ASSEMBLY_FA")"

# ── Build reference list ──────────────────────────────────────────────────────
OUTPUT_DIR="$SAMPLE_DIR/ID_FastANI"
mkdir -p "$OUTPUT_DIR"

REF_LIST="$OUTPUT_DIR/reference_list.txt"

if [[ -f "$REFERENCES" ]]; then
    cp "$REFERENCES" "$REF_LIST"
    REF_COUNT=$(wc -l < "$REF_LIST")
    log "References: $REFERENCES ($REF_COUNT genomes)"
elif [[ -d "$REFERENCES" ]]; then
    find "$REFERENCES" -maxdepth 1 \
        \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \
           -o -name "*.fasta.gz" -o -name "*.fna.gz" \) \
        | sort > "$REF_LIST"
    REF_COUNT=$(wc -l < "$REF_LIST")
    [[ "$REF_COUNT" -gt 0 ]] || die "No genome files found in $REFERENCES"
    log "References: $REFERENCES ($REF_COUNT genomes)"
fi

log "Threshold: ANI >= ${ANI_THRESHOLD}%"
log "Output:    $OUTPUT_DIR"

# ── Output paths ──────────────────────────────────────────────────────────────
RAW_OUT="$OUTPUT_DIR/fastani_raw.txt"
RESULT_TSV="$OUTPUT_DIR/fastani.tsv"

# =============================================================================
header "FastANI — Average Nucleotide Identity"
# =============================================================================

pixi run --environment identification -- fastANI \
  -q "$ASSEMBLY_FA" \
  --rl "$REF_LIST" \
  -t "$THREADS" \
  -o "$RAW_OUT"

# ── Build accession → species name lookup from JSONL ─────────────────────────
LOOKUP_FILE="$OUTPUT_DIR/accession_lookup.txt"

if [[ -f "$JSONL_DB" ]]; then
    log "Building species name lookup from $(basename "$JSONL_DB")..."
    python3 - "$JSONL_DB" "$LOOKUP_FILE" << 'PYEOF'
import json, sys

jsonl_path = sys.argv[1]
lookup_path = sys.argv[2]

with open(jsonl_path) as f, open(lookup_path, 'w') as out:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
            accession = d.get('accession', '')
            organism  = d.get('organism', {}).get('organismName', '')
            if accession and organism:
                out.write(f"{accession}\t{organism}\n")
        except json.JSONDecodeError:
            continue
PYEOF
    log "Lookup built: $(wc -l < "$LOOKUP_FILE") entries"
else
    warn "JSONL not found at $JSONL_DB — species names will not be added"
    touch "$LOOKUP_FILE"
fi

# ── Add header, species name, filter by threshold, sort by ANI ───────────────
python3 - "$RAW_OUT" "$LOOKUP_FILE" "$RESULT_TSV" "$ANI_THRESHOLD" << 'PYEOF'
import sys, os

raw_path      = sys.argv[1]
lookup_path   = sys.argv[2]
result_path   = sys.argv[3]
threshold     = float(sys.argv[4])

# Load accession → species name lookup
lookup = {}
with open(lookup_path) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            lookup[parts[0]] = parts[1]

def accession_from_path(path):
    """Extract GCF_XXXXXX.X from a file path."""
    basename = os.path.basename(path)
    # GCF_000005845.2_ASM584v2_genomic.fna -> GCF_000005845.2
    parts = basename.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return basename

rows = []
with open(raw_path) as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 5:
            continue
        query, ref, ani, aligned, total = cols[:5]
        ani_f = float(ani)
        if ani_f < threshold:
            continue
        accession    = accession_from_path(ref)
        species_name = lookup.get(accession, 'Unknown')
        genome_frac  = f"{(int(aligned)/int(total))*100:.1f}%"
        rows.append((query, ref, ani_f, aligned, total, genome_frac, species_name, accession))

# Sort by ANI descending
rows.sort(key=lambda x: x[2], reverse=True)

with open(result_path, 'w') as out:
    out.write("query\tspecies\tANI\tgenome_fraction\taligned_fragments\ttotal_fragments\taccession\treference_path\n")
    for query, ref, ani_f, aligned, total, genome_frac, species, accession in rows:
        out.write(f"{query}\t{species}\t{ani_f:.4f}\t{genome_frac}\t{aligned}\t{total}\t{accession}\t{ref}\n")
PYEOF

rm -f "$RAW_OUT"

# ── Show results ──────────────────────────────────────────────────────────────
HITS=$(( $(wc -l < "$RESULT_TSV") - 1 ))

success "FastANI done — $HITS hits above ${ANI_THRESHOLD}% ANI"
echo ""

if [[ "$HITS" -gt 0 ]]; then
    log "Top 10 hits:"
    echo ""
    # Print only the readable columns (skip reference_path)
    head -11 "$RESULT_TSV" | cut -f1-7 | column -t -s $'\t'
    echo ""

    BEST_ANI=$(awk -F'\t' 'NR==2{print $3}' "$RESULT_TSV")
    BEST_SPECIES=$(awk -F'\t' 'NR==2{print $2}' "$RESULT_TSV")
    BEST_FRAC=$(awk -F'\t' 'NR==2{print $4}' "$RESULT_TSV")
    BEST_ACC=$(awk -F'\t' 'NR==2{print $7}' "$RESULT_TSV")

    echo "────────────────────────────────────────"
    if (( $(echo "$BEST_ANI >= 95" | bc -l) )); then
        echo -e "  Species:  ${GREEN}${BEST_SPECIES}${RESET}"
        echo -e "  ANI:      ${GREEN}${BEST_ANI}% ✔ Same species (>=95%)${RESET}"
    elif (( $(echo "$BEST_ANI >= 85" | bc -l) )); then
        echo -e "  Species:  ${YELLOW}${BEST_SPECIES}${RESET}"
        echo -e "  ANI:      ${YELLOW}${BEST_ANI}% — Same genus, different species${RESET}"
    else
        echo -e "  Species:  ${RED}${BEST_SPECIES}${RESET}"
        echo -e "  ANI:      ${RED}${BEST_ANI}% — Different genus${RESET}"
    fi
    echo "  Accession:        $BEST_ACC"
    echo "  Genome fraction:  $BEST_FRAC"
    echo "────────────────────────────────────────"
else
    log "No hits found above ${ANI_THRESHOLD}% ANI"
    log "Try:  ANI_THRESHOLD=70 $0 $SAMPLE_DIR"
fi

echo ""
log "Full results: $RESULT_TSV"