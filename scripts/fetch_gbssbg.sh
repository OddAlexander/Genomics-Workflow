#!/usr/bin/env bash
# Install GBS-SBG (S. agalactiae serotyping by genome) into ./.gbssbg/
# at the repo root. Idempotent: skips steps already done.
#
# Override target dir or upstream URL via env vars:
#   GBSSBG_DIR=/some/path  scripts/fetch_gbssbg.sh
#
# Requires the identification pixi env: pixi install
set -euo pipefail

cd "$(dirname "$0")/.."

DIR="${GBSSBG_DIR:-.gbssbg}"
REPO="${GBSSBG_REPO:-https://github.com/swainechen/GBS-SBG.git}"
PIXI_RUN=(pixi run --environment identification)

if [ -d "$DIR/.git" ]; then
    echo "[1/2] repo already cloned at $DIR"
else
    echo "[1/2] cloning $REPO -> $DIR"
    rm -rf "$DIR"
    "${PIXI_RUN[@]}" git clone --depth 1 "$REPO" "$DIR"
fi

if [ ! -f "$DIR/GBS-SBG.pl" ]; then
    echo "Setup incomplete: GBS-SBG.pl not found in $DIR" >&2
    exit 2
fi

echo "[2/2] verifying BLAST is available"
"${PIXI_RUN[@]}" blastn -version 2>&1 | head -1

echo "OK -- GBS-SBG ready at $DIR/"
