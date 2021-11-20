#!/bin/bash


source tests/assert.sh
output=output/test8
nextflow main.nf -profile test,conda --output $output --collect_hs_metrics_min_base_quality 10 \
--collect_hs_metrics_min_mapping_quality 10 --remove_duplicates false --skip_bqsr --skip_realignment

test -s $output/sample1/TESTX_S1_L001.preprocessed.bam || { echo "Missing test 8 output file!"; exit 1; }
test -s $output/sample1/TESTX_S1_L001.preprocessed.bai || { echo "Missing test 8 output file!"; exit 1; }
test -s $output/sample1/metrics/TESTX_S1_L001.prepared.dedup.hs_metrics.txt || { echo "Missing test 8 output file!"; exit 1; }
test -s $output/sample1/metrics/TESTX_S1_L001.prepared.dedup_metrics.txt || { echo "Missing test 8 output file!"; exit 1; }
test -s $output/sample2/TESTX_S1_L002.preprocessed.bam || { echo "Missing test 8 output file!"; exit 1; }
test -s $output/sample2/TESTX_S1_L002.preprocessed.bai || { echo "Missing test 8 output file!"; exit 1; }