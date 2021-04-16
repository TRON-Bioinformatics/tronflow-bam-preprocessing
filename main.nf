#!/usr/bin/env nextflow

publish_dir = 'output'
params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"
params.dbsnp = "/projects/data/gatk_bundle/hg19/dbsnp_138.hg19.vcf"
params.known_indels1 = "/projects/data/gatk_bundle/hg19/1000G_phase1.indels.hg19.sites.vcf"
params.known_indels2 = "/projects/data/gatk_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.sorted.vcf"
params.skip_bqsr = false
params.skip_realignment = false
params.skip_deduplication = false
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
    bam_preprocessing.nf --input_files input_files --reference reference.fasta

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name, sample type (tumor or normal) and path to the BAM file
    Sample type will be added to the BAM header @SN sample name
    The input file does not have header!
    Example input file:
    name1	tumor	tumor.1.bam
    name1	normal	normal.1.bam
    name2	tumor	tumor.2.bam

Optional input:
    * reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    * dbsnp: path to the dbSNP VCF
    * known_indels1: path to a VCF of known indels
    * known_indels2: path to a second VCF of known indels
    * NOTE: if any of the above parameters is not provided, default hg19 resources will be used
    * skip_bqsr: optionally skip BQSR
    * skip_realignment: optionally skip realignment
    * skip_deduplication: optionally skip deduplication
    * output: the folder where to publish output
    * platform: the platform to be added to the BAM header. Valid values: [ILLUMINA, SOLID, LS454, HELICOS and PACBIO] (default: ILLUMINA)
    * prepare_bam_cpus: default 3
    * prepare_bam_memory: default 8g
    * mark_duplicates_cpus: default 16
    * mark_duplicates_memory: default 64g
    * realignment_around_indels_cpus: default 2
    * realignment_around_indels_memory: default 32g
    * bqsr_cpus: default 3
    * bqsr_memory: default 4g

Output:
    * Preprocessed and indexed BAM
    * Tab-separated values file with the absolute paths to the preprocessed BAMs, preprocessed_bams.txt

Optional output:
    * Recalibration report
    * Realignment intervals
    * Duplication metrics
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
      set val(name), val("${bam.baseName}"), val(type),
        file("${bam.baseName}.prepared.bam"), file("${bam.baseName}.prepared.bai")  into prepared_bams

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

    rm -rf tmp
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
	    publishDir "${publish_dir}/${name}", mode: "copy", pattern: "*.dedup_metrics.txt"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from prepared_bams

	    output:
	    	set val(name), val(bam_name), val(type),
	    	    file("${bam.baseName}.dedup.bam"), file("${bam.baseName}.dedup.bam.bai") into deduplicated_bams
	    	file("${bam.baseName}.dedup_metrics.txt") into deduplication_metrics

	    """
	    mkdir tmp

        gatk MarkDuplicatesSpark \
        --java-options '-Xmx${params.mark_duplicates_memory}  -Djava.io.tmpdir=tmp' \
        --input  ${bam} \
        --output ${bam.baseName}.dedup.bam \
        --conf 'spark.executor.cores=${task.cpus}' \
        --metrics-file ${bam.baseName}.dedup_metrics.txt

        rm -rf tmp
	    """
	}
}
else {
    deduplicated_bams = prepared_bams
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

	    rm -rf tmp
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

	    rm -rf tmp
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
