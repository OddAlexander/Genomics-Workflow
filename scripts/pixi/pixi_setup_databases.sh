#!/usr/bin/env bash
# =============================================================================
# setup-databases.sh
# Sets up database environment variables for the bacterial genomics pipeline
#
# Usage:
#   ./setup-databases.sh
#
# This script adds database paths to ~/.bashrc so they persist across sessions.
# Run once on a new system after databases have been downloaded.
#
# Default database paths:
#   /databases/kraken2_db
#   /databases/mash_db/Bacteria_Archaea_type_assembly_set.msh
#   /databases/bakta_db/db
#   /databases/fastANI_refs/reference_list.txt
#   /databases/snippy_refs/latest
# =============================================================================

set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()     { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
success() { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${RESET}"; }
warn()    { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠ $*${RESET}"; }
skip()    { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⊘ $*${RESET}"; }

echo ""
echo -e "${BOLD}╔══════════════════════════════════════════════════════════════╗${RESET}"
echo -e "${BOLD}║       D A T A B A S E   S E T U P                          ║${RESET}"
echo -e "${BOLD}╚══════════════════════════════════════════════════════════════╝${RESET}"
echo ""

BASHRC="$HOME/.bashrc"
ADDED=0
SKIPPED=0

# Function to add export if not already present
add_export() {
    local varname="$1"
    local value="$2"
    local description="$3"
    local check_path="$4"

    # Check if already in bashrc
    if grep -q "export ${varname}=" "$BASHRC" 2>/dev/null; then
        skip "$varname already set in ~/.bashrc — skipping"
        (( SKIPPED++ )) || true
        return
    fi

    # Check if path exists
    if [[ ! -e "$check_path" ]]; then
        warn "$varname: path not found: $check_path"
        warn "  Add manually after downloading: export ${varname}=${value}"
        echo ""
    fi

    echo "export ${varname}=${value}" >> "$BASHRC"
    success "Added $varname=$value"
    echo "         ($description)"
    echo ""
    (( ADDED++ )) || true
}

log "Adding database paths to $BASHRC..."
echo ""

# ── Kraken2 ───────────────────────────────────────────────────────────────────
add_export \
    "KRAKEN2_DB" \
    "/databases/kraken2_db" \
    "Kraken2 + Bracken database directory" \
    "/databases/kraken2_db"

# ── Mash ──────────────────────────────────────────────────────────────────────
add_export \
    "MASH_DB" \
    "/databases/mash_db/Bacteria_Archaea_type_assembly_set.msh" \
    "Mash type strain sketch database" \
    "/databases/mash_db/refseq.msh"

# ── Bakta ─────────────────────────────────────────────────────────────────────
add_export \
    "BAKTA_DB" \
    "/databases/bakta_db/db" \
    "Bakta annotation database" \
    "/databases/bakta_db/db"

# ── FastANI ───────────────────────────────────────────────────────────────────
add_export \
    "FASTANI_REFS" \
    "/databases/fastANI_refs/reference_list.txt" \
    "FastANI reference genome list" \
    "/databases/fastANI_refs/reference_list.txt"

# ── Snippy reference ──────────────────────────────────────────────────────────
add_export \
    "SNIPPY_REF" \
    "/databases/snippy_refs/latest" \
    "Default Snippy reference genome" \
    "/databases/snippy_refs"

# ── CheckM ────────────────────────────────────────────────────────────────────
add_export \
    "CHECKM_DATA_PATH" \
    "$HOME/.checkm" \
    "CheckM marker gene database" \
    "$HOME/.checkm"

# =============================================================================
echo "────────────────────────────────────────────────────────────────"
echo ""
echo "  Added:   $ADDED variable(s)"
echo "  Skipped: $SKIPPED variable(s) (already set)"
echo ""

if [[ "$ADDED" -gt 0 ]]; then
    log "Reloading ~/.bashrc..."
    # shellcheck disable=SC1090
    source "$BASHRC"
    success "Done — database paths are now active"
else
    log "No changes made"
fi

echo ""
echo -e "${BOLD}Current database settings:${RESET}"
echo ""
for var in KRAKEN2_DB MASH_DB BAKTA_DB FASTANI_REFS SNIPPY_REF CHECKM_DATA_PATH; do
    val="${!var:-NOT SET}"
    if [[ "$val" == "NOT SET" ]]; then
        echo -e "  ${YELLOW}$var = NOT SET${RESET}"
    elif [[ -e "$val" ]]; then
        echo -e "  ${GREEN}$var = $val ✔${RESET}"
    else
        echo -e "  ${YELLOW}$var = $val (path not found)${RESET}"
    fi
done

echo ""
echo -e "${BOLD}Download databases:${RESET}"
echo ""
echo "  Kraken2 PlusPF (80GB):"
echo "    wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz"
echo "    mkdir -p /databases/kraken2_db"
echo "    tar -xzf k2_pluspf_20260226.tar.gz -C /databases/kraken2_db"
echo ""
echo "  Mash type strains (~500MB):"
echo "    mkdir -p /databases/mash_db"
echo "    wget https://figshare.com/ndownloader/files/24753609 \\"
echo "         -O /databases/mash_db/Bacteria_Archaea_type_assembly_set.msh"
echo ""
echo "  FastANI references:"
echo "    datasets download genome taxon enterobacteriaceae \\"
echo "      --assembly-level complete --reference \\"
echo "      --filename /databases/fastANI_refs/refs.zip"
echo "    unzip /databases/fastANI_refs/refs.zip -d /databases/fastANI_refs/"
echo "    find /databases/fastANI_refs/ -name '*.fna' \\"
echo "      > /databases/fastANI_refs/reference_list.txt"
echo ""
echo "  CheckM database:"
echo "    pixi run --environment checkm -- checkm data setRoot ~/.checkm"
echo ""
