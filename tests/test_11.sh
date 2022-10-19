#!/bin/bash


source tests/assert.sh
output=output/test11
nextflow main.nf -profile test,conda --output $output --reference `pwd`/test_data/ucsc.hg19.minimal.without_indices.fasta

test -s `pwd`/test_data/ucsc.hg19.minimal.without_indices.fasta.fai || { echo "Missing output FAI index!"; exit 1; }
test -s `pwd`/test_data/ucsc.hg19.minimal.without_indices.dict || { echo "Missing output dict index!"; exit 1; }
