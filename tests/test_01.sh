#!/bin/bash


source tests/assert.sh
output=output/test1
nextflow main.nf -profile test,conda --output $output

test -s $output/sample1/sample1.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample1/sample1.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample1/software_versions.PREPARE_BAM.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.INDEX_BAM.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.MARK_DUPLICATES.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.HS_METRICS.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.METRICS.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.COVERAGE_ANALYSIS.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.REALIGNMENT_AROUND_INDELS.txt || { echo "Missing software versions file!"; exit 1; }
test -s $output/sample1/software_versions.BQSR.txt || { echo "Missing software versions file!"; exit 1; }

test -s $output/sample2/sample2.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }