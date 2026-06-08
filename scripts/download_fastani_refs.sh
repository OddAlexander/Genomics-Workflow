#!/usr/bin/env bash
# Download RefSeq bacteria reference genomes for use as FastANI reference database.
# Uses the official NCBI datasets CLI with dehydrated download for robustness on large sets.
#
# Scope: complete + chromosome-level, RefSeq reference genomes (~6 900 genomes)
# Output: /databases/fastani_db/genomes/  + /databases/fastani_db/reference_list.txt
#
# Usage:
#   pixi run --environment identification bash scripts/download_fastani_refs.sh
#   pixi run --environment identification bash scripts/download_fastani_refs.sh --all
#
# Then set in config.yaml:
#   fastani_db: /databases/fastani_db/reference_list.txt

set -euo pipefail

OUT_DIR="/databases/fastani_db"
ZIP="$OUT_DIR/bacteria_dehydrated.zip"
DATA_DIR="$OUT_DIR/ncbi_dataset"
LIST_FILE="$OUT_DIR/reference_list.txt"
THREADS="${THREADS:-16}"
ALL_LEVELS=false

for arg in "$@"; do
    [[ "$arg" == "--all" ]] && ALL_LEVELS=true
done

if $ALL_LEVELS; then
    LEVELS="complete,chromosome,scaffold,contig"
    echo "Mode: all assembly levels"
else
    LEVELS="complete,chromosome"
    echo "Mode: complete + chromosome assemblies only (recommended)"
fi

sudo mkdir -p "$OUT_DIR"
sudo chown -R "$(id -u):$(id -g)" "$OUT_DIR"

# ── Step 1: Download dehydrated package (catalog only, no sequences yet) ──────
echo "Downloading genome catalog..."
datasets download genome taxon bacteria \
    --assembly-source refseq \
    --assembly-level "$LEVELS" \
    --reference \
    --include genome \
    --dehydrated \
    --filename "$ZIP"

unzip -q -o "$ZIP" -d "$DATA_DIR"

# ── Step 2: Rehydrate (download actual genome files, resumable) ───────────────
echo "Downloading genome sequences (resumable)..."
datasets rehydrate --directory "$DATA_DIR" --max-workers "$THREADS"

# ── Step 3: Build reference list ──────────────────────────────────────────────
echo "Building reference list..."
find "$DATA_DIR" -name "*.fna" | sort > "$LIST_FILE"

COUNT=$(wc -l < "$LIST_FILE")
echo "Done: $COUNT genomes written to $LIST_FILE"
