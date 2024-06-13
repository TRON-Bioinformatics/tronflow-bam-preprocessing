#!/bin/bash


source tests/assert.sh
output=output/test9
nextflow main.nf -profile test,conda,ci --output $output --skip_deduplication --skip_bqsr --skip_realignment  \
--input_files false --input_bam test_data/TESTX_S1_L001.bam

test -s $output/TESTX_S1_L001/TESTX_S1_L001.preprocessed.bam || { echo "Missing test 9 output file!"; exit 1; }
test -s $output/TESTX_S1_L001/TESTX_S1_L001.preprocessed.bai || { echo "Missing test 9 output file!"; exit 1; }
