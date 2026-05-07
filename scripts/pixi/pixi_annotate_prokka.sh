#!/usr/bin/env bash
# Annotate one or more assembled FASTA files with Prokka to produce reference genomes.
# Output: GenBank (.gbk) and GFF3 (.gff) files suitable for use with Snippy, etc.
#
# Usage:
#   bash scripts/pixi/pixi_annotate_prokka.sh 19-03-2026              # alle samples fra en dato
#   bash scripts/pixi/pixi_annotate_prokka.sh 19-03-2026/005a         # ett enkelt sample
#   bash scripts/pixi/pixi_annotate_prokka.sh 19-03-2026 --genus Staphylococcus --species aureus
#
# Input resolves mot results/ i prosjektroten. Fulle stier fungerer også.
# Output lands in --outdir (default: /databases/prokka_refs/).
# Requires: pixi (used internally to run prokka in the prokka environment)

set -euo pipefail

PIXI_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT=""
OUTDIR="/databases/prokka_refs"
GENUS=""
SPECIES=""
STRAIN=""
THREADS="${THREADS:-8}"
RESULTS_DIR="$PIXI_ROOT/results"

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--outdir)   OUTDIR="$2";      shift 2 ;;
        --genus)       GENUS="$2";       shift 2 ;;
        --species)     SPECIES="$2";     shift 2 ;;
        --strain)      STRAIN="$2";      shift 2 ;;
        -t|--threads)  THREADS="$2";     shift 2 ;;
        --results-dir) RESULTS_DIR="$2"; shift 2 ;;
        -*) echo "Ukjent flagg: $1" >&2; exit 1 ;;
        *)  INPUT="$1"; shift ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "Bruk: $(basename "$0") <dato|sample-sti> [--genus X --species Y ...]" >&2
    echo "Eksempel: $(basename "$0") 19-03-2026/005a --genus Staphylococcus --species aureus" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Build optional prokka flags ────────────────────────────────────────────────
TAXON_FLAGS=""
[[ -n "$GENUS"   ]] && TAXON_FLAGS+=" --genus $GENUS"
[[ -n "$SPECIES" ]] && TAXON_FLAGS+=" --species $SPECIES"
[[ -n "$STRAIN"  ]] && TAXON_FLAGS+=" --strain $STRAIN"

# ── Resolve input → list of contigs.fa paths ──────────────────────────────────
# Accepts: full file path, full dir path, or a date/sample shorthand
# resolved against $RESULTS_DIR (e.g. "19-03-2026" or "19-03-2026/005a")
declare -a FASTAS=()

resolve_dir() {
    local dir="$1"
    while IFS= read -r -d '' fa; do
        FASTAS+=("$fa")
    done < <(find "$dir" -path "*/Assembly/contigs.fa" -print0 | sort -z)
}

if [[ -f "$INPUT" ]]; then
    FASTAS+=("$INPUT")
elif [[ -d "$INPUT" ]]; then
    resolve_dir "$INPUT"
else
    CANDIDATE="$RESULTS_DIR/$INPUT"
    if [[ -f "$CANDIDATE" ]]; then
        FASTAS+=("$CANDIDATE")
    elif [[ -d "$CANDIDATE" ]]; then
        resolve_dir "$CANDIDATE"
    else
        echo "Feil: finner ikke '$INPUT' som fil, katalog, eller under $RESULTS_DIR/" >&2
        exit 1
    fi
fi

if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "Feil: ingen Assembly/contigs.fa funnet for '$INPUT'" >&2
    exit 1
fi

echo "Annoterer ${#FASTAS[@]} assembly(ies) → $OUTDIR"
echo ""

# ── Run prokka ─────────────────────────────────────────────────────────────────
DONE=0
FAIL=0

for FA in "${FASTAS[@]}"; do
    # results/19-03-2026/005a/Assembly/contigs.fa → prefix 19-03-2026_005a
    SAMPLE_DIR=$(dirname "$(dirname "$FA")")   # strip /Assembly/contigs.fa
    DATE_PART=$(basename "$(dirname "$SAMPLE_DIR")")
    SAMPLE_PART=$(basename "$SAMPLE_DIR")
    PREFIX=$(echo "${DATE_PART}_${SAMPLE_PART}" | tr -cd '[:alnum:]_.-')
    SAMPLE_OUT="$OUTDIR/$PREFIX"

    if [[ -f "$SAMPLE_OUT/${PREFIX}.gbk" ]]; then
        echo "[$PREFIX] Allerede annotert — hopper over (slett $SAMPLE_OUT for å kjøre på nytt)"
        (( DONE++ )) || true
        continue
    fi

    echo "[$PREFIX] Annoterer $FA ..."
    if pixi run --environment prokka -- \
            prokka \
            --outdir "$SAMPLE_OUT" \
            --prefix "$PREFIX" \
            --cpus   "$THREADS" \
            --rfam \
            --compliant \
            $TAXON_FLAGS \
            "$FA" 2>&1; then
        echo "[$PREFIX] Ferdig: $SAMPLE_OUT/${PREFIX}.gbk"
        (( DONE++ )) || true
    else
        echo "[$PREFIX] FEIL — sjekk prokka-logg i $SAMPLE_OUT/" >&2
        (( FAIL++ )) || true
    fi
    echo ""
done

echo "─────────────────────────────────────────"
echo "Fullført: $DONE   Feilet: $FAIL"
echo "Referansefiler (.gbk): $OUTDIR/"
