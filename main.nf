#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_BAM; SORT_AND_INDEX_BAM } from './modules/01_prepare_bam'
include { MARK_DUPLICATES; SPLIT_CIGAR_N_READS } from './modules/02_mark_duplicates'
include { METRICS; HS_METRICS; COVERAGE_ANALYSIS; FLAGSTAT } from './modules/03_metrics'
include { REALIGNMENT_AROUND_INDELS } from './modules/04_realignment_around_indels'
include { BQSR; CREATE_OUTPUT } from './modules/05_bqsr'
include { CREATE_FAIDX; CREATE_DICT } from './modules/00_reference_indices'


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

if (! params.input_files && ! params.input_bam) {
  exit 1, "Neither --input_files or --input_bam are provided!"
}
else if (params.input_files && params.input_bam) {
  exit 1, "Both --input_files and --input_bam are provided! Please, provide only one."
}
else if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'type', 'bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.type, file(row.bam)) }
    .set { input_files }
} else if (params.input_bam && params.input_name) {
  input_bam = file(params.input_bam)
  Channel
    .fromList([tuple(input_bam.name.take(input_bam.name.lastIndexOf('.')), params.input_name, input_bam)])
    .set { input_files }
}

workflow CHECK_REFERENCE {

    take:
        reference
    main:
        reference_file = file(reference)

        if (!reference_file.exists()) {
            error "--reference points to a non existing file"
        }

        if (!file("${reference}.fai").exists()) {
            CREATE_FAIDX(reference)
        }

        dict_file = reference.toString().replaceFirst(/\.(fa|fasta)$/, ".dict")

        if (!file(dict_file).exists()) {
            CREATE_DICT(reference)
        }
    emit:
        checked_reference = reference
}


workflow {

    CHECK_REFERENCE(params.reference)

    // PREPARE_BAM(input_files, CHECK_REFERENCE.out.checked_reference)

    if (!params.skip_deduplication) {
        MARK_DUPLICATES(input_files)
        deduplicated_bams = MARK_DUPLICATES.out.deduplicated_bams
    }
    else {
        SORT_AND_INDEX_BAM(input_files)
        deduplicated_bams = SORT_AND_INDEX_BAM.out.sorted_and_indexed_bams
    }

    if (params.split_cigarn) {
        SPLIT_CIGAR_N_READS(deduplicated_bams, CHECK_REFERENCE.out.checked_reference)
        deduplicated_bams = SPLIT_CIGAR_N_READS.out.split_cigarn_bams
    }

    if (! params.skip_metrics) {
        if (params.intervals) {
            HS_METRICS(deduplicated_bams)
        }
        METRICS(deduplicated_bams, CHECK_REFERENCE.out.checked_reference)
        COVERAGE_ANALYSIS(deduplicated_bams)
        FLAGSTAT(deduplicated_bams)
    }

    if (!params.skip_realignment) {
        REALIGNMENT_AROUND_INDELS(deduplicated_bams, CHECK_REFERENCE.out.checked_reference)
        realigned_bams = REALIGNMENT_AROUND_INDELS.out.realigned_bams
    }
    else {
        realigned_bams = deduplicated_bams
    }

    if (!params.skip_bqsr) {
        BQSR(realigned_bams, CHECK_REFERENCE.out.checked_reference)
        preprocessed_bams = BQSR.out.recalibrated_bams
    }
    else {
        CREATE_OUTPUT(realigned_bams)
        preprocessed_bams = CREATE_OUTPUT.out.recalibrated_bams
    }

    preprocessed_bams
        .map {it.join("\t")}
        .collectFile(name: "${params.output}/preprocessed_bams.txt", newLine: true)
}
