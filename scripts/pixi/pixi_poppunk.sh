#!/usr/bin/env bash
# =============================================================================
# run-poppunk.sh
# Assigns PopPUNK cluster type (CT) to a bacterial assembly
#
# Usage:
#   ./run-poppunk.sh <sample_dir> [database]
#
# Arguments:
#   sample_dir  — sample folder containing an Assembly/ subfolder
#   database    — (optional) path to PopPUNK database directory
#                 Default: /databases/poppunk/Enterococcus_faecium_v2_full
#
# Output:
#   <sample_dir>/PopPUNK/
#     poppunk_clusters.csv   Cluster assignment
#     poppunk_summary.tsv    Clean summary: sample, cluster, database
#
# Pre-built databases: https://www.bacpop.org/poppunk/
# Common databases:
#   /databases/poppunk/Enterococcus_faecium_v2_full   (default)
#   /databases/poppunk/Klebsiella_pneumoniae_v2_refs
#   /databases/poppunk/Staphylococcus_aureus_v2_refs
#
# Note: PopPUNK unwords.py requires a manual patch after installation.
# Run this once after pixi install --environment poppunk:
#   FILE=$(find ~/genomics/.pixi/envs/poppunk -name "unwords.py" | head -1)
#   # Then replace pkg_resources with importlib.resources (see docs)
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

DEFAULT_DB="/databases/poppunk/Enterococcus_faecium_v2_full"

# ── No arguments — show help ──────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "============================================================"
    echo "   P O P P U N K   C L U S T E R   T Y P I N G            "
    echo "============================================================"
    echo -e "${BOLD}Usage:${RESET} $0 <sample_dir> [database]"
    echo ""
    echo -e "${BOLD}Arguments:${RESET}"
    echo "  sample_dir  Sample folder containing Assembly/ subfolder"
    echo "  database    (optional) PopPUNK database path"
    echo "              Default: $DEFAULT_DB"
    echo ""
    echo -e "${BOLD}Available databases:${RESET}"
    echo "  /databases/poppunk/Enterococcus_faecium_v2_full    ← default"
    echo "  /databases/poppunk/Klebsiella_pneumoniae_v2_refs"
    echo "  /databases/poppunk/Staphylococcus_aureus_v2_refs"
    echo ""
    echo -e "${BOLD}Download databases from:${RESET}"
    echo "  https://www.bacpop.org/poppunk/"
    echo ""
    echo -e "${BOLD}Examples:${RESET}"
    echo "  $0 results/001k/"
    echo "  $0 results/001k/ /databases/poppunk/Klebsiella_pneumoniae_v2_refs"
    echo "============================================================"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
SAMPLE_DIR="${1%/}"
DATABASE="${2:-$DEFAULT_DB}"

[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -d "$DATABASE" ]]   || die "PopPUNK database not found: $DATABASE\nDownload from: https://www.bacpop.org/poppunk/"

# ── Find assembly ─────────────────────────────────────────────────────────────
ASSEMBLY_DIR="$SAMPLE_DIR/Assembly"
[[ -d "$ASSEMBLY_DIR" ]] || die "Assembly/ not found in $SAMPLE_DIR\nRun bacterial-pipeline.sh first"

if [[ -f "$ASSEMBLY_DIR/scaffolds.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/scaffolds.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fasta" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fasta"
elif [[ -f "$ASSEMBLY_DIR/contigs.fa" ]]; then
    ASSEMBLY_FA="$ASSEMBLY_DIR/contigs.fa"
else
    die "No assembly file found in $ASSEMBLY_DIR"
fi

# Convert to absolute paths — PopPUNK requires this
ASSEMBLY_FA=$(realpath "$ASSEMBLY_FA")
DATABASE=$(realpath "$DATABASE")

SAMPLE=$(basename "$SAMPLE_DIR")
OUTPUT_DIR="$SAMPLE_DIR/PopPUNK"
QUERY_FILE="$OUTPUT_DIR/query.txt"
SUMMARY_TSV="$OUTPUT_DIR/poppunk_summary.tsv"

mkdir -p "$OUTPUT_DIR"

log "Sample:   $SAMPLE"
log "Assembly: $(basename "$ASSEMBLY_FA")"
log "Database: $(basename "$DATABASE")"
log "Output:   $OUTPUT_DIR"

# =============================================================================
header "PopPUNK — cluster type assignment"
# =============================================================================

# Create query file (tab-separated: name <tab> assembly_path)
echo -e "${SAMPLE}\t${ASSEMBLY_FA}" > "$QUERY_FILE"
log "Query file: $QUERY_FILE"

STEP_START=$SECONDS

pixi run --environment poppunk -- poppunk_assign \
    --db "$DATABASE" \
    --query "$QUERY_FILE" \
    --output "$OUTPUT_DIR" \
    --threads 8 \


ELAPSED=$(( SECONDS - STEP_START ))

success "PopPUNK done (${ELAPSED}s)"

# =============================================================================
header "Results"
# =============================================================================

# Find the cluster output file
CLUSTER_FILE=$(find "$OUTPUT_DIR" -name "*_clusters.csv" 2>/dev/null | head -1)

if [[ -z "$CLUSTER_FILE" ]]; then
    warn "Cluster output file not found — check $OUTPUT_DIR for output"
    exit 1
fi

# Parse cluster assignment
CLUSTER=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $2}' "$CLUSTER_FILE" \
          || awk -F',' 'NR==2{print $2}' "$CLUSTER_FILE")

DB_NAME=$(basename "$DATABASE")

# Write summary TSV
{
    echo -e "sample\tcluster\tdatabase"
    echo -e "${SAMPLE}\t${CLUSTER}\t${DB_NAME}"
} > "$SUMMARY_TSV"

# Print result
echo ""
echo "  ┌─────────────────────────────────────────────┐"
printf "  │  %-12s  %-28s │\n" "Sample:" "$SAMPLE"
printf "  │  %-12s  %-28s │\n" "Cluster:" "$CLUSTER"
printf "  │  %-12s  %-28s │\n" "Database:" "${DB_NAME:0:28}"
echo "  └─────────────────────────────────────────────┘"
echo ""

# E. faecium clade interpretation
if [[ "$DB_NAME" == *"faecium"* ]] || [[ "$DB_NAME" == *"Efaecium"* ]]; then
    echo "  Database context: Enterococcus faecium"
    echo ""
    echo "  PopPUNK clusters broadly correspond to:"
    echo "    Clade A1 — hospital-adapted, high-risk VRE lineage"
    echo "    Clade A2 — transitional"
    echo "    Clade B  — community/animal, lower hospital risk"
    echo ""
    warn "Cross-reference cluster with MLST ST for full clinical context"
    warn "Clade assignment requires comparison against published cluster-clade tables"
fi

echo ""
echo "────────────────────────────────────────────────────────"
echo "  Summary TSV:   $SUMMARY_TSV"
echo "  Full output:   $CLUSTER_FILE"
echo "────────────────────────────────────────────────────────"
