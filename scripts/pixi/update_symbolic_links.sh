#!/usr/bin/env bash
# =============================================================================
# install-scripts.sh
# Creates symbolic links in ~/bin/ for all pipeline scripts
# so they are runnable from anywhere
#
# Usage:
#   ./install-scripts.sh [scripts_dir]
#
# Arguments:
#   scripts_dir — path to your scripts folder
#                 Default: ~/genomics/scripts/pixi
#
# Run once, or re-run after adding new scripts.
# =============================================================================

set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()     { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
success() { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${RESET}"; }
warn()    { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠ $*${RESET}"; }
skip()    { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⊘ $*${RESET}"; }
die()     { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ $*${RESET}" >&2; exit 1; }

SCRIPTS_DIR="${1:-$HOME/genomics/scripts/pixi}"
BIN_DIR="$HOME/bin"

[[ -d "$SCRIPTS_DIR" ]] || die "Scripts directory not found: $SCRIPTS_DIR"

# Create ~/bin if it doesn't exist
mkdir -p "$BIN_DIR"

# Ensure ~/bin is in PATH
if ! echo "$PATH" | grep -q "$BIN_DIR"; then
    warn "~/bin is not in your PATH"
    echo "  Add this to ~/.bashrc:"
    echo "  export PATH=\"\$HOME/bin:\$PATH\""
    echo ""
fi

echo ""
echo -e "${BOLD}Scripts directory: $SCRIPTS_DIR${RESET}"
echo -e "${BOLD}Symlink target:    $BIN_DIR${RESET}"
echo ""

LINKED=0
UPDATED=0
SKIPPED=0

for script in "$SCRIPTS_DIR"/*.sh; do
    [[ -f "$script" ]] || continue

    name=$(basename "$script")
    target="$BIN_DIR/$name"

    # Make sure script is executable
    chmod +x "$script"

    if [[ -L "$target" ]]; then
        # Already a symlink — check if it points to the right place
        current=$(readlink "$target")
        if [[ "$current" == "$script" ]]; then
            skip "$name — already linked"
            (( SKIPPED++ )) || true
        else
            # Update stale symlink
            ln -sf "$script" "$target"
            warn "$name — updated (was: $current)"
            (( UPDATED++ )) || true
        fi
    elif [[ -f "$target" ]]; then
        warn "$name — file already exists at $target (not a symlink, skipping)"
        (( SKIPPED++ )) || true
    else
        ln -s "$script" "$target"
        success "$name → $target"
        (( LINKED++ )) || true
    fi
done

echo ""
echo "────────────────────────────────────────"
echo "  Linked:  $LINKED"
echo "  Updated: $UPDATED"
echo "  Skipped: $SKIPPED"
echo "────────────────────────────────────────"
echo ""

if [[ $(( LINKED + UPDATED )) -gt 0 ]]; then
    log "Done — scripts are now available from anywhere"
    echo ""
    echo "  Example:"
    echo "    pixi_QC.sh results/001k/"
    echo "    pixi_assemble.sh data/001k/ results/"
fi

# Check if ~/bin is in PATH and add if missing
if ! grep -q 'HOME/bin' "$HOME/.bashrc" 2>/dev/null; then
    echo ""
    warn "Adding ~/bin to PATH in ~/.bashrc..."
    echo '' >> "$HOME/.bashrc"
    echo '# User scripts' >> "$HOME/.bashrc"
    echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
    success "Added to ~/.bashrc — run: source ~/.bashrc"
fi
