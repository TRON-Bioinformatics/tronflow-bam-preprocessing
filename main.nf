#!/usr/bin/env nextflow

publish_dir = 'output'
params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"
params.dbsnp = "/projects/data/gatk_bundle/hg19/dbsnp_138.hg19.vcf"
params.known_indels1 = "/projects/data/gatk_bundle/hg19/1000G_phase1.indels.hg19.sites.vcf"
params.known_indels2 = "/projects/data/gatk_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.sorted.vcf"
params.intervals = false
params.hs_metrics_target_coverage = false
params.hs_metrics_per_base_coverage = false
params.skip_bqsr = false
params.skip_realignment = false
params.skip_deduplication = false
params.skip_metrics = false
params.output = false
params.platform = "ILLUMINA"

params.prepare_bam_cpus = 3
params.prepare_bam_memory = "8g"
params.mark_duplicates_cpus = 16
params.mark_duplicates_memory = "64g"
params.realignment_around_indels_cpus = 2
params.realignment_around_indels_memory = "32g"
params.bqsr_cpus = 3
params.bqsr_memory = "4g"



def helpMessage() {
    log.info"""
Usage:
    main.nf --input_files input_files

Input:
    * --input_files: the path to a tab-separated values file containing in each row the sample name, sample type (eg: tumor or normal) and path to the BAM file
    Sample type will be added to the BAM header @SN sample name
    The input file does not have header!
    Example input file:
    name1       tumor   tumor.1.bam
    name1       normal  normal.1.bam
    name2       tumor   tumor.2.bam

Optional input:
    * --reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    * --dbsnp: path to the dbSNP VCF
    * --known_indels1: path to a VCF of known indels
    * --known_indels2: path to a second VCF of known indels
    **NOTE**: if any of the above parameters is not provided, default hg19 resources under
    /projects/data/gatk_bundle/hg19/ will be used

    * --intervals: path to an intervals file to collect HS metrics from, this can be built with Picard's BedToIntervalList (default: None)
    * --hs_metrics_target_coverage: name of output file for target HS metrics (default: None)
    * --hs_metrics_per_base_coverage: name of output file for per base HS metrics (default: None)
    * --skip_bqsr: optionally skip BQSR (default: false)
    * --skip_realignment: optionally skip realignment (default: false)
    * --skip_deduplication: optionally skip deduplication (default: false)
    * --skip_metrics: optionally skip metrics (default: false)
    * --output: the folder where to publish output (default: ./output)
    * --platform: the platform to be added to the BAM header. Valid values: [ILLUMINA, SOLID, LS454, HELICOS and PACBIO] (default: ILLUMINA)

Computational resources:
    * --prepare_bam_cpus: (default: 3)
    * --prepare_bam_memory: (default: 8g)
    * --mark_duplicates_cpus: (default: 16)
    * --mark_duplicates_memory: (default: 64g)
    * --realignment_around_indels_cpus: (default: 2)
    * --realignment_around_indels_memory: (default: 32g)
    * --bqsr_cpus: (default: 3)
    * --bqsr_memory: (default: 4g)

 Output:
    * Preprocessed and indexed BAMs
    * Tab-separated values file with the absolute paths to the preprocessed BAMs, preprocessed_bams.txt

Optional output:
    * Recalibration report
    * Realignment intervals
    * Metrics
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if (params.output) {
  publish_dir = params.output
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'type', 'bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.type, file(row.bam)) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

/*
This step sets MAPQ to 0 for all unmapped reads + avoids soft clipping beyond the end of the reference genome
This step reorders chromosomes in the BAM file according to the provided reference (this step is required for GATK)
Adds the required read groups fields to the BAM file. The provided type is added to the BAM sample name.
*/
process prepareBam {
    cpus "${params.prepare_bam_cpus}"
    memory "${params.prepare_bam_memory}"
    tag "${name}"

    input:
    	set name, type, file(bam) from input_files

    output:
      set val(name),
        val("${bam.baseName}"),
        val(type), file("${bam.baseName}.prepared.bam"),
        file("${bam.baseName}.prepared.bai")  into prepared_bams, prepared_bams_for_metrics, prepared_bams_for_hs_metrics

    """
    mkdir tmp

    gatk CleanSam \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=tmp' \
    --INPUT ${bam} \
    --OUTPUT /dev/stdout | \
    gatk ReorderSam \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=tmp' \
    --INPUT /dev/stdin \
    --OUTPUT /dev/stdout \
    --SEQUENCE_DICTIONARY ${params.reference} | \
    gatk AddOrReplaceReadGroups \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=tmp' \
    --VALIDATION_STRINGENCY SILENT \
    --INPUT /dev/stdin \
    --OUTPUT ${bam.baseName}.prepared.bam \
    --REFERENCE_SEQUENCE ${params.reference} \
    --RGPU 1 \
    --RGID 1 \
    --RGSM ${type} \
    --RGLB 1 \
    --RGPL ${params.platform} \
    --SORT_ORDER coordinate \
    --CREATE_INDEX true
    """
}

/*
Adds the appropriate read groups to the BAM file.
The provided type is added to the BAM sample name.
*/
if (!params.skip_deduplication) {
	process markDuplicates {
	    cpus "${params.mark_duplicates_cpus}"
        memory "${params.mark_duplicates_memory}"
	    tag "${name}"
	    publishDir "${publish_dir}/${name}/metrics", mode: "copy", pattern: "*.dedup_metrics"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from prepared_bams

	    output:
	    	set val(name), val(bam_name), val(type),
	    	    file("${bam.baseName}.dedup.bam"), file("${bam.baseName}.dedup.bam.bai") into deduplicated_bams
	    	file("${bam.baseName}.dedup_metrics") optional true into deduplication_metrics

        script:
        dedup_metrics = params.skip_metrics ? "--metrics-file ${bam.baseName}.dedup_metrics" : ""
	    """
	    mkdir tmp

        gatk MarkDuplicatesSpark \
        --java-options '-Xmx${params.mark_duplicates_memory}  -Djava.io.tmpdir=tmp' \
        --input  ${bam} \
        --output ${bam.baseName}.dedup.bam \
        --conf 'spark.executor.cores=${task.cpus}' \
        ${dedup_metrics}
	    """
	}
}
else {
    deduplicated_bams = prepared_bams
}

if (! params.skip_metrics) {

    if (params.intervals) {

        process hsMetrics {
            cpus 1
            memory "2g"
            tag "${name}"
            publishDir "${publish_dir}/${name}/metrics", mode: "copy"

            input:
                set name, bam_name, type, file(bam), file(bai) from prepared_bams_for_hs_metrics

            output:
                file("*_metrics") optional true into txt_hs_metrics
                file("*.pdf") optional true into pdf_hs_metrics
                file(params.hs_metrics_target_coverage) optional true into target_hs_metrics
                file(params.hs_metrics_per_base_coverage) optional true into per_base_hs_metrics

            script:
            hs_metrics_target_coverage= params.hs_metrics_target_coverage ?
                "--PER_TARGET_COVERAGE ${params.hs_metrics_target_coverage} --REFERENCE_SEQUENCE ${params.reference}" :
                ""
            hs_metrics_per_base_coverage= params.hs_metrics_per_base_coverage ?
                "--PER_BASE_COVERAGE ${params.hs_metrics_per_base_coverage}" :
                ""
            """
            mkdir tmp

            gatk CollectHsMetrics \
            --java-options '-Xmx2g  -Djava.io.tmpdir=tmp' \
            --INPUT  ${bam} \
            --OUTPUT ${bam.baseName} \
            --TARGET_INTERVALS ${params.intervals} \
            --BAIT_INTERVALS ${params.intervals} \
            ${hs_metrics_target_coverage} ${hs_metrics_per_base_coverage}
            """
        }
    }

    process metrics {
	    cpus 1
        memory "2g"
	    tag "${name}"
	    publishDir "${publish_dir}/${name}/metrics", mode: "copy"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from prepared_bams_for_metrics

	    output:
	        file("*_metrics") optional true into txt_metrics
	        file("*.pdf") optional true into pdf_metrics

	    """
	    mkdir tmp

	    gatk CollectMultipleMetrics \
        --java-options '-Xmx2g  -Djava.io.tmpdir=tmp' \
        --INPUT  ${bam} \
        --OUTPUT ${bam.baseName} \
        --REFERENCE_SEQUENCE ${params.reference} \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectBaseDistributionByCycle \
        --PROGRAM CollectGcBiasMetrics \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM CollectSequencingArtifactMetrics \
        --PROGRAM CollectSequencingArtifactMetrics
	    """
	}
}

if (!params.skip_realignment) {
	process realignmentAroundindels {
	    cpus "${params.realignment_around_indels_cpus}"
        memory "${params.realignment_around_indels_memory}"
	    tag "${name}"
	    publishDir "${publish_dir}/${name}", mode: "copy", pattern: "*.RA.intervals"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from deduplicated_bams

	    output:
	    	set val(name), val(bam_name), val(type), file("${bam.baseName}.realigned.bam"), file("${bam.baseName}.realigned.bai") into realigned_bams
	    	file("${bam.baseName}.RA.intervals") into realignment_intervals

	    """
	    mkdir tmp

	    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T RealignerTargetCreator \
	    --input_file ${bam} \
	    --out ${bam.baseName}.RA.intervals \
	    --reference_sequence ${params.reference} \
	    --known ${params.known_indels1} \
	    --known ${params.known_indels2}

	    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T IndelRealigner \
	    --input_file ${bam} \
	    --out ${bam.baseName}.realigned.bam \
	    --reference_sequence ${params.reference} \
	    --targetIntervals ${bam.baseName}.RA.intervals \
	    --knownAlleles ${params.known_indels1} \
	    --knownAlleles ${params.known_indels2} \
	    --consensusDeterminationModel USE_SW \
	    --LODThresholdForCleaning 0.4 \
	    --maxReadsInMemory 600000
	    """
	}
}
else {
    realigned_bams = deduplicated_bams
}

if (!params.skip_bqsr) {
	process baseQualityScoreRecalibration {
	    cpus "${params.bqsr_cpus}"
        memory "${params.bqsr_memory}"
	    publishDir "${publish_dir}/${name}", mode: "copy"
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from realigned_bams

	    output:
	    	set val("${name}"), val("${type}"), val("${publish_dir}/${name}/${bam_name}.preprocessed.bam") into recalibrated_bams
            file "${bam_name}.recalibration_report.grp" into recalibration_report
            file "${bam_name}.preprocessed.bam" into recalibrated_bam
            file "${bam_name}.preprocessed.bai" into recalibrated_bai

	    """
	    mkdir tmp

	    gatk BaseRecalibrator \
	    --java-options '-Xmx${params.bqsr_memory} -Djava.io.tmpdir=tmp' \
	    --input ${bam} \
	    --output ${bam_name}.recalibration_report.grp \
	    --reference ${params.reference} \
	    --known-sites ${params.dbsnp}

	    gatk ApplyBQSR \
	    --java-options '-Xmx${params.bqsr_memory} -Djava.io.tmpdir=tmp' \
	    --input ${bam} \
	    --output ${bam_name}.preprocessed.bam \
	    --reference ${params.reference} \
	    --bqsr-recal-file ${bam_name}.recalibration_report.grp
	    """
	}
}
else {
	process createOutput {
	    cpus 1
	    memory '1g'
	    publishDir "${publish_dir}/${name}", mode: "copy"
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from realigned_bams

	    output:
	    	set val("${name}"), val("${type}"), val("${publish_dir}/${name}/${bam_name}.preprocessed.bam") into recalibrated_bams
			file "${bam_name}.preprocessed.bam" into recalibrated_bam
			file "${bam_name}.preprocessed.bai" into recalibrated_bai

		"""
		cp ${bam} ${bam_name}.preprocessed.bam
		cp ${bai} ${bam_name}.preprocessed.bai
		"""
	}
}

recalibrated_bams
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/preprocessed_bams.txt", newLine: true)
