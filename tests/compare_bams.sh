#!/bin/bash
set -euo pipefail

# Script to compare two preprocessed BAM files
# Usage: ./compare_bams.sh <bam_old_path> <bam_new_path>

CONDA_ENV="bam-compare-env"
ENV_NAME="${CONDA_ENV}"

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <bam_old_path> <bam_new_path>"
    exit 1
fi

BAM_OLD="$1"
BAM_NEW="$2"

# Check if BAM files exist
if [ ! -f "$BAM_OLD" ]; then
    echo "Error: BAM file not found: $BAM_OLD"
    exit 1
fi

if [ ! -f "$BAM_NEW" ]; then
    echo "Error: BAM file not found: $BAM_NEW"
    exit 1
fi

echo "================================================"
echo "BAM File Comparison Script"
echo "================================================"
echo "BAM OLD: $BAM_OLD"
echo "BAM NEW: $BAM_NEW"
echo ""

# Ensure conda exists
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda command not found. Please install Conda/Miniconda first."
    exit 1
fi

# Create conda environment if it doesn't exist
echo "Setting up conda environment: $ENV_NAME"
if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Conda environment $ENV_NAME already exists, using it..."
else
    echo "Creating new conda environment with required tools..."
    if ! conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda --strict-channel-priority samtools=1.19; then
        echo "Pinned solve failed, retrying with unpinned samtools..."
        conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda --strict-channel-priority samtools
    fi
fi

# Use samtools through conda without relying on shell activation
SAMTOOLS_CMD=(conda run --no-capture-output -n "$ENV_NAME" samtools)

# Create output directory for comparison results
OUTPUT_DIR="bam_comparison_results"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "================================================"
echo "1. FLAGSTAT COMPARISON (Read Statistics)"
echo "================================================"
echo "File OLD: $BAM_OLD"
"${SAMTOOLS_CMD[@]}" flagstat "$BAM_OLD" | tee "$OUTPUT_DIR/bam_old.flagstat"
echo ""
echo "File NEW: $BAM_NEW"
"${SAMTOOLS_CMD[@]}" flagstat "$BAM_NEW" | tee "$OUTPUT_DIR/bam_new.flagstat"
echo ""
echo "Differences:"
diff "$OUTPUT_DIR/bam_old.flagstat" "$OUTPUT_DIR/bam_new.flagstat" || echo "Flagstat statistics differ between files"

echo ""
echo "================================================"
echo "2. INDEX STATISTICS (Per-Chromosome)"
echo "================================================"
echo "File OLD: $BAM_OLD"
"${SAMTOOLS_CMD[@]}" idxstats "$BAM_OLD" | tee "$OUTPUT_DIR/bam_old.idxstats"
echo ""
echo "File NEW: $BAM_NEW"
"${SAMTOOLS_CMD[@]}" idxstats "$BAM_NEW" | tee "$OUTPUT_DIR/bam_new.idxstats"
echo ""
echo "Differences:"
diff "$OUTPUT_DIR/bam_old.idxstats" "$OUTPUT_DIR/bam_new.idxstats" || echo "Index statistics differ between files"

echo ""
echo "================================================"
echo "3. SAM HEADER COMPARISON"
echo "================================================"
"${SAMTOOLS_CMD[@]}" view -H "$BAM_OLD" > "$OUTPUT_DIR/bam_old.header.sam"
"${SAMTOOLS_CMD[@]}" view -H "$BAM_NEW" > "$OUTPUT_DIR/bam_new.header.sam"
echo "File OLD header lines:"
wc -l "$OUTPUT_DIR/bam_old.header.sam"
echo "File NEW header lines:"
wc -l "$OUTPUT_DIR/bam_new.header.sam"
echo ""
echo "Header differences:"
diff "$OUTPUT_DIR/bam_old.header.sam" "$OUTPUT_DIR/bam_new.header.sam" || echo "Headers differ between files"

echo ""
echo "================================================"
echo "4. READ COUNT AND MAPPING STATISTICS"
echo "================================================"
BAM_OLD_TOTAL=$("${SAMTOOLS_CMD[@]}" view -c "$BAM_OLD")
BAM_NEW_TOTAL=$("${SAMTOOLS_CMD[@]}" view -c "$BAM_NEW")
BAM_OLD_MAPPED=$("${SAMTOOLS_CMD[@]}" view -c -F 4 "$BAM_OLD")
BAM_NEW_MAPPED=$("${SAMTOOLS_CMD[@]}" view -c -F 4 "$BAM_NEW")

echo "File OLD total reads: $BAM_OLD_TOTAL"
echo "File OLD mapped reads: $BAM_OLD_MAPPED"
echo ""
echo "File NEW total reads: $BAM_NEW_TOTAL"
echo "File NEW mapped reads: $BAM_NEW_MAPPED"
echo ""

if [ "$BAM_OLD_TOTAL" -eq "$BAM_NEW_TOTAL" ]; then
    echo "OK Total read counts match"
else
    DIFF=$((BAM_OLD_TOTAL - BAM_NEW_TOTAL))
    echo "DIFF Total read counts differ by: $DIFF"
fi

if [ "$BAM_OLD_MAPPED" -eq "$BAM_NEW_MAPPED" ]; then
    echo "OK Mapped read counts match"
else
    DIFF=$((BAM_OLD_MAPPED - BAM_NEW_MAPPED))
    echo "DIFF Mapped read counts differ by: $DIFF"
fi

echo ""
echo "================================================"
echo "5. DEPTH STATISTICS (Average Coverage)"
echo "================================================"
echo "File OLD: $BAM_OLD"
DEPTH_OLD=$("${SAMTOOLS_CMD[@]}" depth "$BAM_OLD" | awk '{sum+=$3; count++} END {if (count > 0) printf "%.2f\n", sum/count; else print "0"}')
echo "Average depth: $DEPTH_OLD"

echo ""
echo "File NEW: $BAM_NEW"
DEPTH_NEW=$("${SAMTOOLS_CMD[@]}" depth "$BAM_NEW" | awk '{sum+=$3; count++} END {if (count > 0) printf "%.2f\n", sum/count; else print "0"}')
echo "Average depth: $DEPTH_NEW"

echo ""
echo "================================================"
echo "SUMMARY"
echo "================================================"
echo "Detailed results saved in: $OUTPUT_DIR/"
echo "  - bam_old.flagstat"
echo "  - bam_new.flagstat"
echo "  - bam_old.idxstats"
echo "  - bam_new.idxstats"
echo "  - bam_old.header.sam"
echo "  - bam_new.header.sam"
echo ""
echo "Comparison complete!"
