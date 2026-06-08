#!/usr/bin/env bash
# Build a rapid Neighbor-Joining tree from assembled genomes using MashTree.
# Useful for a quick overview before committing to a full SNP analysis.
#
# Usage:
#   bash scripts/pixi/pixi_mashtree.sh 19-03-2026
#   bash scripts/pixi/pixi_mashtree.sh 19-03-2026 26-03-2026
#   bash scripts/pixi/pixi_mashtree.sh results/19-03-2026/005a/Assembly/contigs.fa ...
#
# Each positional argument is resolved against results/ if it is not an
# existing path — same shorthand as pixi_annotate_prokka.sh.
# Output goes to --outdir (default: results_phylo/MashTree/).

set -euo pipefail

PIXI_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESULTS_DIR="$PIXI_ROOT/results"
OUTDIR="$PIXI_ROOT/results_phylo/MashTree"
THREADS="${THREADS:-8}"
declare -a INPUTS=()

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--outdir)   OUTDIR="$2";   shift 2 ;;
        -t|--threads)  THREADS="$2";  shift 2 ;;
        -*) echo "Ukjent flagg: $1" >&2; exit 1 ;;
        *)  INPUTS+=("$1"); shift ;;
    esac
done

if [[ ${#INPUTS[@]} -eq 0 ]]; then
    echo "Bruk: $(basename "$0") <dato|sti> [<dato|sti> ...] [--outdir <dir>]" >&2
    echo "Eksempel: $(basename "$0") 19-03-2026 26-03-2026" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Resolve inputs → list of contigs.fa paths ─────────────────────────────────
declare -a FASTAS=()

resolve() {
    local arg="$1"
    if [[ -f "$arg" ]]; then
        FASTAS+=("$(realpath "$arg")")
    elif [[ -d "$arg" ]]; then
        while IFS= read -r -d '' fa; do
            FASTAS+=("$(realpath "$fa")")
        done < <(find "$arg" -path "*/Assembly/contigs.fa" -print0 | sort -z)
    else
        local candidate="$RESULTS_DIR/$arg"
        if [[ -f "$candidate" ]]; then
            FASTAS+=("$(realpath "$candidate")")
        elif [[ -d "$candidate" ]]; then
            while IFS= read -r -d '' fa; do
                FASTAS+=("$(realpath "$fa")")
            done < <(find "$candidate" -path "*/Assembly/contigs.fa" -print0 | sort -z)
        else
            echo "Advarsel: finner ikke '$arg' — hopper over" >&2
        fi
    fi
}

for INPUT in "${INPUTS[@]}"; do
    resolve "$INPUT"
done

if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "Feil: ingen Assembly/contigs.fa funnet for gitt input" >&2
    exit 1
fi

if [[ ${#FASTAS[@]} -lt 3 ]]; then
    echo "Advarsel: MashTree krever minst 3 genomer for et meningsfullt tre (fant ${#FASTAS[@]})" >&2
fi

echo "MashTree: ${#FASTAS[@]} genomer → $OUTDIR"
echo ""

# ── Run MashTree ───────────────────────────────────────────────────────────────
TREE_FILE="$OUTDIR/mashtree.dnd"
MATRIX_FILE="$OUTDIR/mashtree_distances.tsv"

pixi run --environment identification -- \
    mashtree \
        --numcpus "$THREADS" \
        --outmatrix "$MATRIX_FILE" \
        "${FASTAS[@]}" \
    > "$TREE_FILE"

echo "Tre:      $TREE_FILE"
echo "Matrise:  $MATRIX_FILE"
echo ""

# ── Print tree to terminal ─────────────────────────────────────────────────────
echo "Newick-tre:"
cat "$TREE_FILE"
