#!/bin/bash


source tests/assert.sh
output=output/test6
nextflow main.nf -profile test,conda --output $output --intervals false --skip_deduplication --skip_bqsr --skip_realignment

test -s $output/sample1/sample1.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample1/sample1.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }