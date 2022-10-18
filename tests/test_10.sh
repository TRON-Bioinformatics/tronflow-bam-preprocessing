#!/bin/bash


source tests/assert.sh
output=output/test10
nextflow main.nf -profile test,conda --output $output --skip_realignment --split_cigarn

test -s $output/sample1/sample1.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample1/sample1.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample1/metrics/flagstat/sample1.flagstat.csv || { echo "Missing output metrics file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bai || { echo "Missing output BAI file!"; exit 1; }
test -s $output/sample2/metrics/flagstat/sample2.flagstat.csv || { echo "Missing output metrics file!"; exit 1; }