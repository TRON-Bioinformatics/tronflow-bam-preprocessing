#!/bin/bash


source tests/assert.sh
output=output/test4
nextflow main.nf -profile test,conda,ci --output $output --skip_deduplication

test -s $output/sample1/sample1.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample1/sample1.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }