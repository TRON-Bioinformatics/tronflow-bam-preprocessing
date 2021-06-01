all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --output output/test1
	nextflow main.nf -profile test,conda --skip_bqsr --output output/test2
	nextflow main.nf -profile test,conda --skip_realignment --output output/test3
	nextflow main.nf -profile test,conda --skip_deduplication --output output/test4
	nextflow main.nf -profile test,conda --output output/test5 --skip_metrics --known_indels1 false --known_indels2 false
	nextflow main.nf -profile test,conda --output output/test6 --intervals false
	nextflow main.nf -profile test,conda --output output/test7 --hs_metrics_target_coverage target_coverage.txt --hs_metrics_per_base_coverage per_base_coverage.txt
	nextflow main.nf -profile test,conda --output output/test8 --hs_metrics_target_coverage target_coverage.txt --hs_metrics_per_base_coverage per_base_coverage.txt --collect_hs_metrics_min_base_quality 10 --collect_hs_metrics_min_mapping_quality 10 --remove_duplicates false

check:
	test -s output/test1/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test3/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test3/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test3/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test3/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test4/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test4/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test4/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test4/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test5/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test5/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test5/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test5/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test6/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test6/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test6/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test6/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test7/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test7/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test7/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test7/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test8/sample1/sample1.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test8/sample1/sample1.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test8/sample2/sample2.preprocessed.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test8/sample2/sample2.preprocessed.bai || { echo "Missing test 1 output file!"; exit 1; }