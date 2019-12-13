#!/usr/bin/env nextflow

picard_jar = "/code/picard/2.21.2/picard.jar"
// NOTE: we need GATK 3 for the realignment around indels as it has been discontinued
gatk3_jar = "/code/gatk/3.8.1.0/GenomeAnalysisTK.jar"
gatk4_jar = "/code/gatk/4.1.3.0/gatk-package-4.1.3.0-local.jar"
publish_dir = 'output'

params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"						// TODO: remove this hard coded bit
params.dbsnp = "/projects/data/gatk_bundle/hg19/dbsnp_138.hg19.vcf" 							// TODO: remove this hard coded bit
params.known_indels1 = "/projects/data/gatk_bundle/hg19/1000G_phase1.indels.hg19.sites.vcf"			// TODO: remove this hard coded bit
params.known_indels2 = "/projects/data/gatk_bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.sorted.vcf"	// TODO: remove this hard coded bit
params.skip_bqsr = false
params.skip_realignment = false
params.skip_deduplication = false
params.output = false
params.platform = "ILLUMINA"

def helpMessage() {
    log.info"""
Usage:
    bam_preprocessing.nf --input_files input_files --reference reference.fasta

This workflow is based on the implementation at /code/iCaM/scripts/mutect.sh

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

Output:
    * Preprocessed and indexed BAM
    * Tab-separated values file with the absolute paths to the preprocessed BAMs, preprocessed_bams.txt

Optional output:
    * Recalibration report
    * Realignment intervals
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
    cpus 3
    memory '8g'
    module 'java/1.8.0'
    module 'bioinf/samtools/1.9'
    tag "${name}"

    input:
    	set name, type, file(bam) from input_files

    output:
      set val(name), val("${bam.baseName}"), val(type), file("${bam.baseName}.prepared.bam")  into prepared_bams

    """
    java -Xmx8g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${picard_jar} \
    CleanSam \
    I=${bam} \
    O=${bam.baseName}.cleaned.bam \
    TMP_DIR=`pwd`/scratch/tmp

    samtools index ${bam.baseName}.cleaned.bam

    java -Xmx8g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${picard_jar} \
    ReorderSam \
    I=${bam.baseName}.cleaned.bam \
    O=${bam.baseName}.reordered.bam \
    SEQUENCE_DICTIONARY=${params.reference} \
    TMP_DIR=`pwd`/scratch/tmp

    rm -f ${bam.baseName}.cleaned.bam

    java -Xmx8g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${picard_jar} \
    AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=SILENT \
    I=${bam.baseName}.reordered.bam \
    O=${bam.baseName}.prepared.bam \
    R=${params.reference} \
    PU=1 \
    ID=1 \
    SM=${type} \
    LB=1 \
    PL=${params.platform} \
    SO=coordinate \
    TMP_DIR=`pwd`/scratch/tmp

    rm -f ${bam.baseName}.reordered.bam
    """
}

/*
Adds the appropriate read groups to the BAM file.
The provided type is added to the BAM sample name.
*/
if (!params.skip_deduplication) {
	process markDuplicates {
	    cpus 8
	    memory '64g'
	    module 'java/1.8.0'
	    module 'bioinf/samtools/1.9'
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam) from prepared_bams

	    output:
	    	set val(name), val(bam_name), val(type), file("${bam.baseName}.dedup.bam"), file("${bam.baseName}.dedup.bam.bai") into deduplicated_bams

	    """
	    # --create-output-bam-index false is required due to https://github.com/broadinstitute/gatk/issues/5919
      mkdir -p `pwd`/scratch/tmp
	    java -Xmx64g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk4_jar} \
      MarkDuplicatesSpark \
      --input  ${bam} \
      --output ${bam.baseName}.dedup.bam \
      --conf 'spark.executor.cores=8' \
      --create-output-bam-index false \
      -M ${bam.baseName}.dedup_metrics.txt

	    samtools index ${bam.baseName}.dedup.bam

      mv ${bam.baseName}.dedup_metrics.txt ${publish_dir}
	    """
	}
}
else {
	process skipMarkDuplicates {
	    cpus 1
	    memory '4g'
	    module 'bioinf/samtools/1.9'
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam) from prepared_bams

	    output:
	    	set val(name), val(bam_name), val(type), file("${bam}"), file("${bam.baseName}.bam.bai") into deduplicated_bams

	    """
	    samtools index ${bam}
	    """
	}
}

if (!params.skip_realignment) {
	process realignmentAroundindels {
	    cpus 2
	    memory '32g'
	    module 'java/1.8.0'		// GATK requires Java 8
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from deduplicated_bams

	    output:
	    	set val(name), val(bam_name), val(type), file("${bam.baseName}.realigned.bam"), file("${bam.baseName}.realigned.bai") into realigned_bams

	    """
      mkdir -p `pwd`/scratch/tmp
	    java -Xmx32g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk3_jar} \
      -T RealignerTargetCreator \
	    --input_file ${bam} \
	    --out ${bam.baseName}.RA.intervals \
	    --reference_sequence ${params.reference} \
	    --known ${params.known_indels1} \
	    --known ${params.known_indels2}

	    java -Xmx32g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk3_jar} \
      -T IndelRealigner \
	    --input_file ${bam} \
	    --out ${bam.baseName}.realigned.bam \
	    --reference_sequence ${params.reference} \
	    --targetIntervals ${bam.baseName}.RA.intervals \
	    --knownAlleles ${params.known_indels1} \
	    --knownAlleles ${params.known_indels2} \
	    --consensusDeterminationModel USE_SW \
	    --LODThresholdForCleaning 0.4 \
	    --maxReadsInMemory 600000

      mv ${bam.baseName}.RA.intervals ${publish_dir}
	    """
	}
}
else {
	realigned_bams = deduplicated_bams
}

if (!params.skip_bqsr) {
	process baseQualityScoreRecalibration {
	    cpus 3
	    memory '4g'
	    module 'java/1.8.0'		// GATK requires Java 8
	    publishDir "${publish_dir}", mode: "move"
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from realigned_bams

	    output:
	    	set val("${name}"), val("${type}"), val("${publish_dir}/${bam_name}.preprocessed.bam") into recalibrated_bams
        file "${bam_name}.recalibration_report.grp" into recalibration_report
        file "${bam_name}.preprocessed.bam" into recalibrated_bam
        file "${bam_name}.preprocessed.bai" into recalibrated_bai

	    """
      mkdir -p `pwd`/scratch/tmp
	    java -Xmx4g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk4_jar} \
      BaseRecalibrator \
	    --input ${bam} \
	    --output ${bam_name}.recalibration_report.grp \
	    --reference ${params.reference} \
	    --known-sites:VCF ${params.dbsnp}

	    java -Xmx4g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk4_jar} \
      ApplyBQSR \
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
	    publishDir "${publish_dir}", mode: "move"
	    tag "${name}"

	    input:
	    	set name, bam_name, type, file(bam), file(bai) from realigned_bams

	    output:
	    	set val("${name}"), val("${type}"), val("${publish_dir}/${bam_name}.preprocessed.bam") into recalibrated_bams
		file "${bam_name}.preprocessed.bam" into recalibrated_bam
		file "${bam_name}.preprocessed.bai" into recalibrated_bai

	    """
	    mv ${bam} ${bam_name}.preprocessed.bam
	    mv ${bai} ${bam_name}.preprocessed.bai
	    """
	}
}

recalibrated_bams
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/preprocessed_bams.txt", newLine: true)
