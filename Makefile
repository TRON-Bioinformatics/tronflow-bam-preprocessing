clean:
	rm -rf output
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	nextflow main.nf -profile test,conda --output output/test1
	#nextflow main.nf -profile test,conda --skip_bqsr --output output/test2
	#nextflow main.nf -profile test,conda --skip_realignment --output output/test3
	#nextflow main.nf -profile test,conda --skip_deduplication --output output/test4
	#nextflow main.nf -profile test,conda --output output/test5 --skip_metrics
	nextflow main.nf -profile test,conda --output output/test6 --intervals false
	nextflow main.nf -profile test,conda --output output/test7 --hs_metrics_target_coverage target_coverage.txt --hs_metrics_per_base_coverage per_base_coverage.txt
	#nextflow main.nf -profile test,conda --output output/test8 --hs_metrics_target_coverage target_coverage.txt --hs_metrics_per_base_coverage per_base_coverage.txt --collect_hs_metrics_min_base_quality 10 --collect_hs_metrics_min_mapping_quality 10
