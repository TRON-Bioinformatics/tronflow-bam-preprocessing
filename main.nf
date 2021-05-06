#!/usr/bin/env nextflow

publish_dir = 'output'
params.help= false
params.input_files = false
params.reference = false
params.dbsnp = false
params.known_indels1 = false
params.known_indels2 = false
params.intervals = false
params.hs_metrics_target_coverage = false
params.hs_metrics_per_base_coverage = false
params.skip_bqsr = false
params.skip_realignment = false
params.skip_deduplication = false
params.skip_metrics = false
params.output = false
params.platform = "ILLUMINA"
params.collect_hs_metrics_min_base_quality = false
params.collect_hs_metrics_min_mapping_quality = false

// computational resources
params.prepare_bam_cpus = 3
params.prepare_bam_memory = "8g"
params.mark_duplicates_cpus = 16
params.mark_duplicates_memory = "64g"
params.realignment_around_indels_cpus = 2
params.realignment_around_indels_memory = "32g"
params.bqsr_cpus = 3
params.bqsr_memory = "4g"



def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.reference) {
    log.error "--reference is required"
    exit 1
}

if (!params.skip_bqsr && !params.dbsnp) {
    log.error "--dbsnp is required to perform BQSR"
    exit 1
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
        dedup_metrics = params.skip_metrics ? "": "--metrics-file ${bam.baseName}.dedup_metrics"
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
            minimum_base_quality = params.collect_hs_metrics_min_base_quality ?
                "--MINIMUM_BASE_QUALITY ${params.collect_hs_metrics_min_base_quality}" : ""
            minimum_mapping_quality = params.collect_hs_metrics_min_mapping_quality ?
                "--MINIMUM_MAPPING_QUALITY ${params.collect_hs_metrics_min_mapping_quality}" : ""
            """
            mkdir tmp

            gatk CollectHsMetrics \
            --java-options '-Xmx2g  -Djava.io.tmpdir=tmp' \
            --INPUT  ${bam} \
            --OUTPUT ${bam.baseName} \
            --TARGET_INTERVALS ${params.intervals} \
            --BAIT_INTERVALS ${params.intervals} \
            ${hs_metrics_target_coverage} ${hs_metrics_per_base_coverage} ${minimum_base_quality} ${minimum_mapping_quality}
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

        script:
        known_indels = "" + params.known_indels1 ? " --known ${params.known_indels1}" : "" +
            params.known_indels2 ? " --known ${params.known_indels2}" : ""
        known_alleles = "" + params.known_indels1 ? " --knownAlleles ${params.known_indels1}" : "" +
            params.known_indels2 ? " --knownAlleles ${params.known_indels2}" : ""
	    """
	    mkdir tmp

	    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T RealignerTargetCreator \
	    --input_file ${bam} \
	    --out ${bam.baseName}.RA.intervals \
	    --reference_sequence ${params.reference} \
	    ${known_indels}

	    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T IndelRealigner \
	    --input_file ${bam} \
	    --out ${bam.baseName}.realigned.bam \
	    --reference_sequence ${params.reference} \
	    --targetIntervals ${bam.baseName}.RA.intervals \
	    --consensusDeterminationModel USE_SW \
	    --LODThresholdForCleaning 0.4 \
	    --maxReadsInMemory 600000 \
	    ${known_alleles}
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
