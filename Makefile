clean:
	rm -rf output
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	#nextflow main.nf -profile test,conda --output output/test1
	#nextflow main.nf -profile test,conda --skip_bqsr --output output/test2
	#nextflow main.nf -profile test,conda --skip_realignment --output output/test3
	#nextflow main.nf -profile test,conda --skip_deduplication --output output/test4
	nextflow main.nf -profile test,conda --output output/test5 --skip_metrics
