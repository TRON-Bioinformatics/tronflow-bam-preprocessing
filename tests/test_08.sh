#!/bin/bash


source tests/assert.sh
output=output/test8
nextflow main.nf -profile test,conda --output $output --collect_hs_metrics_min_base_quality 10 \
--collect_hs_metrics_min_mapping_quality 10 --remove_duplicates false --skip_bqsr --skip_realignment

test -s $output/sample1/sample1.preprocessed.bam || { echo "Missing BAM file!"; exit 1; }
test -s $output/sample1/sample1.preprocessed.bai || { echo "Missing BAI file!"; exit 1; }
test -s $output/sample1/metrics/hs_metrics/sample1.hs_metrics.txt || { echo "Missing HS metrics file!"; exit 1; }
test -s $output/sample1/metrics/flagstat/sample1.flagstat.csv || { echo "Missing dedup metrics file!"; exit 1; }
test -s $output/sample1/metrics/coverage/sample1.coverage.tsv || { echo "Missing horizontal coverage file!"; exit 1; }
test -s $output/sample1/metrics/coverage/sample1.depth.tsv || { echo "Missing depth of coverage file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bam || { echo "Missing BAM file!"; exit 1; }
test -s $output/sample2/sample2.preprocessed.bai || { echo "Missing BAI file!"; exit 1; }
