#!/bin/bash


source tests/assert.sh
output=output/test5
nextflow main.nf -profile test,conda --output $output --skip_deduplication --skip_bqsr --skip_metrics --known_indels1 false --known_indels2 false

test -s $output/sample1/TESTX_S1_L001.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample1/TESTX_S1_L001.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample2/TESTX_S1_L002.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample2/TESTX_S1_L002.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }