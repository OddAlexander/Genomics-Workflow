#!/usr/bin/env bash
# Install LRE-Finder (linezolid resistance detection for Enterococcus) into
# ./.lrefinder/ at the repo root. Idempotent: skips steps that are already done.
#
# Four things have to happen for LRE-Finder to be runnable, all handled below:
#   1. Clone the repo
#   2. Extract elmDB.tar.gz (ships only the FASTA + metadata, not the KMA index)
#   3. Compile getGene + sibling C utilities via `make`
#   4. Build the KMA binary index from elm.fsa (elm.length.b / elm.seq.b / ...)
#
# Override the target dir or upstream URL via env vars:
#   LREFINDER_DIR=/some/path  scripts/fetch_lrefinder.sh
#   LREFINDER_REPO=https://github.com/genomicepidemiology/LRE-Finder.git scripts/fetch_lrefinder.sh
#
# Requires the lrefinder pixi env:  pixi install
set -euo pipefail

cd "$(dirname "$0")/.."

DIR="${LREFINDER_DIR:-.lrefinder}"
REPO="${LREFINDER_REPO:-https://bitbucket.org/genomicepidemiology/lre-finder.git}"
PIXI_RUN=(pixi run --environment lrefinder)

# 1. Clone the repo (skip if already there)
if [ -d "$DIR/.git" ]; then
    echo "[1/4] repo already cloned at $DIR"
else
    echo "[1/4] cloning $REPO -> $DIR"
    rm -rf "$DIR"
    "${PIXI_RUN[@]}" git clone --depth 1 --recurse-submodules "$REPO" "$DIR"
fi

# 2. Extract the bundled database tarball (FASTA + metadata only -- the KMA
# binary index is built in step 4 below).
if [ -f "$DIR/elmDB/elm.fsa" ]; then
    echo "[2/4] elmDB already extracted"
elif [ -f "$DIR/elmDB.tar.gz" ]; then
    echo "[2/4] extracting elmDB.tar.gz"
    "${PIXI_RUN[@]}" tar -xzf "$DIR/elmDB.tar.gz" -C "$DIR"
else
    echo "[2/4] WARNING: no elmDB.tar.gz found; LRE-Finder won't have a database to query" >&2
fi

# 3. Compile getGene (and sibling C utilities) -- LRE-Finder calls getGene at runtime
if [ -x "$DIR/getGene" ]; then
    echo "[3/4] getGene already built"
elif [ -f "$DIR/Makefile" ]; then
    echo "[3/4] running make"
    (cd "$DIR" && "${PIXI_RUN[@]}" make)
else
    echo "[3/4] WARNING: no Makefile -- skipping C build (rule may fail at runtime)" >&2
fi

# 4. Build the KMA binary index from elm.fsa (the .tar.gz ships only the source)
if [ -f "$DIR/elmDB/elm.length.b" ]; then
    echo "[4/4] KMA index already built"
elif [ -f "$DIR/elmDB/elm.fsa" ]; then
    echo "[4/4] building KMA index"
    "${PIXI_RUN[@]}" kma_index -i "$DIR/elmDB/elm.fsa" -o "$DIR/elmDB/elm"
else
    echo "[4/4] WARNING: no elmDB/elm.fsa -- can't build the KMA index" >&2
fi

# Final sanity check
missing=0
for f in LRE-Finder.py elmDB/elm.length.b elmDB/elm.mutdist getGene; do
    if [ ! -e "$DIR/$f" ]; then
        echo "  missing: $DIR/$f" >&2
        missing=$((missing + 1))
    fi
done

if [ "$missing" -gt 0 ]; then
    echo "Setup incomplete: $missing required file(s) missing -- LRE-Finder rule will fail" >&2
    exit 2
fi
echo "OK -- LRE-Finder ready at $DIR/"
