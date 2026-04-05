#!/bin/bash

# --- Configuration ---
ENV_NAME="kleborate_env"
DEFAULT_PRESET="kpsc"
CPU_THREADS=4

# --- Stylish Help Menu ---
usage() {
    echo "============================================================"
    echo "            K L E B O R A T E   P I P E L I N E             "
    echo "============================================================"
    echo "Usage: ./kleborate_script.sh <input_file> <output_dir> [preset]"
    echo ""
    echo "ARGUMENTS:"
    echo "  1. input_file  : Path to a single assembly (e.g., contigs.fasta)"
    echo "  2. output_dir  : Folder to store results"
    echo "  3. preset      : Analysis target (Default: $DEFAULT_PRESET)"
    echo ""
    echo "AVAILABLE PRESETS:"
    echo "  kpsc        : K. pneumoniae species complex (AMR, Vir, K/O)"
    echo "  kosc        : K. oxytoca species complex (MLST, AMR)"
    echo "  escherichia : Escherichia genus (MLST, Pathotyping)"
    echo "  species     : Species identification ONLY (Fastest)"
    echo ""
    echo "EXAMPLES:"
    echo "  ./kleborate_script.sh contigs.fasta ./results"
    echo "  ./kleborate_script.sh isolate_01.fna ./out escherichia"
    echo "============================================================"
    exit 1
}

# Show menu if arguments are missing
if [ $# -lt 2 ]; then
    usage
fi

# --- Variables ---
INPUT_FILE=$1
FINAL_OUT_DIR=$2
PRESET=${3:-$DEFAULT_PRESET}

TEMP_KLEB_DIR="${FINAL_OUT_DIR}/raw_tmp_$(date +%s)"
FINAL_TSV="${FINAL_OUT_DIR}/kleborate_results.tsv"

# --- Validation ---
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found!"
    exit 1
fi

mkdir -p "$FINAL_OUT_DIR"

# --- Execution ---
echo "--> Environment: $ENV_NAME"
echo "--> Analyzing:   $(basename "$INPUT_FILE")"
echo "--> Preset:      $PRESET"
echo "--> Please wait..."

conda run -n "$ENV_NAME" kleborate \
    --preset "$PRESET" \
    --trim_headers \
    -a "$INPUT_FILE" \
    -o "$TEMP_KLEB_DIR" \
    --threads "$CPU_THREADS"

# --- Post-Processing ---
if [ -d "$TEMP_KLEB_DIR" ]; then
    # Merge result files
    first_file=true
    for f in "$TEMP_KLEB_DIR"/*.txt; do
        if [ "$first_file" = true ]; then
            cat "$f" > "$FINAL_TSV"
            first_file=false
        else
            tail -n +2 "$f" >> "$FINAL_TSV"
        fi
    done

    # Clean up the v3 long headers
    sed -i '1s/[^[:space:]]*mlst__ST/ST/g' "$FINAL_TSV"
    sed -i '1s/[^[:space:]]*species_complex__species/species/g' "$FINAL_TSV"

    # Cleanup temp folder
    rm -rf "$TEMP_KLEB_DIR"

    echo "------------------------------------------------------------"
    echo "SUCCESS: Results saved to $FINAL_TSV"
    echo "------------------------------------------------------------"
    
    # Final Visual Output
    column -t -s $'\t' "$FINAL_TSV"
else
    echo "ERROR: Kleborate did not produce any results. Check your input file."
    exit 1
fi
