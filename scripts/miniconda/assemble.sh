#!/bin/bash

# --- 1. CONFIGURATION ---
ENV_NAME="genomics"
THREADS=4
BAKTA_DB="/databases/bakta_db/db"

# --- 2. ARGUMENT & CONDA SETUP ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

IN_DIR=$1
OUT_DIR=$2

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# --- 3. PREPARATION ---
mkdir -p "$OUT_DIR/qc" "$OUT_DIR/assembly" "$OUT_DIR/annotation" "$OUT_DIR/amr" "$OUT_DIR/summary_report"

# --- 4. THE WORKFLOW ---

echo "--- [1/5] Quality Control (fastp) ---"
fastp -i "$IN_DIR"/*_R1.fastq.gz -I "$IN_DIR"/*_R2.fastq.gz \
      -o "$OUT_DIR/qc/clean_R1.fq.gz" -O "$OUT_DIR/qc/clean_R2.fq.gz" \
      --html "$OUT_DIR/qc/fastp_report.html" --thread $THREADS
FP_STATUS=$?

echo "--- [2/5] Assembly (SPAdes) ---"
spades.py -1 "$OUT_DIR/qc/clean_R1.fq.gz" -2 "$OUT_DIR/qc/clean_R2.fq.gz" \
          --isolate -o "$OUT_DIR/assembly/" -t $THREADS
SP_STATUS=$?

echo "--- [3/5] Annotation (Bakta) ---"
bakta --db "$BAKTA_DB" --force --keep-contig-headers --threads $THREADS \
      --output "$OUT_DIR/annotation/" "$OUT_DIR/assembly/scaffolds.fasta"
BK_STATUS=$?

echo "--- [4/5] AMR & Virulence Detection (AMRFinderPlus) ---"
# -n runs on nucleotide assembly, -O specifies the organism if known (optional)
amrfinder -n "$OUT_DIR/assembly/scaffolds.fasta" \
          --threads $THREADS \
          -o "$OUT_DIR/amr/amr_results.tsv"
AMR_STATUS=$?

echo "--- [5/5] MultiQC (Final Summary) ---"
multiqc "$OUT_DIR" -o "$OUT_DIR/summary_report/" --title "Bacterial Genomics Report"
MQ_STATUS=$?

conda deactivate

# --- 5. FINAL SUMMARY ---
echo -e "\n===================================================="
echo "                WORKFLOW SUMMARY"
echo "===================================================="

if [ $FP_STATUS -eq 0 ] && [ $SP_STATUS -eq 0 ] && [ $BK_STATUS -eq 0 ] && [ $AMR_STATUS -eq 0 ]; then
    echo "OVERALL STATUS: SUCCESS ✅"
else
    echo "OVERALL STATUS: COMPLETED WITH ERRORS ⚠️"
    echo "  - fastp: $FP_STATUS | SPAdes: $SP_STATUS | Bakta: $BK_STATUS | AMR: $AMR_STATUS"
fi

echo -e "\nREPORT LOCATIONS:"
echo "QC:      $OUT_DIR/qc/fastp_report.html"
echo "AMR:     $OUT_DIR/amr/amr_results.tsv"
echo "MASTER:  $OUT_DIR/summary_report/multiqc_report.html"
echo "===================================================="
