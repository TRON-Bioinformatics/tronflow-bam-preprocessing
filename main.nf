#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_BAM; INDEX_BAM } from './modules/01_prepare_bam'
include { MARK_DUPLICATES; SPLIT_CIGAR_N_READS } from './modules/02_mark_duplicates'
include { METRICS; HS_METRICS; COVERAGE_ANALYSIS; FLAGSTAT } from './modules/03_metrics'
include { REALIGNMENT_AROUND_INDELS } from './modules/04_realignment_around_indels'
include { BQSR; CREATE_OUTPUT } from './modules/05_bqsr'

params.help= false
params.input_files = false
params.input_name = "normal"
params.input_bam = false
params.reference = false
params.dbsnp = false
params.known_indels1 = false
params.known_indels2 = false
params.intervals = false
params.skip_bqsr = false
params.skip_realignment = false
params.skip_deduplication = false
params.remove_duplicates = true
params.skip_metrics = false
params.output = 'output'
params.platform = "ILLUMINA"
params.collect_hs_metrics_min_base_quality = false
params.collect_hs_metrics_min_mapping_quality = false
params.split_cigarn = false

// computational resources
params.prepare_bam_cpus = 3
params.prepare_bam_memory = "8g"
params.mark_duplicates_cpus = 2
params.mark_duplicates_memory = "16g"
params.realignment_around_indels_cpus = 2
params.realignment_around_indels_memory = "31g"
params.bqsr_cpus = 3
params.bqsr_memory = "4g"
params.metrics_cpus = 1
params.metrics_memory = "8g"
params.index_cpus = 1
params.index_memory = "8g"



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


workflow {

    PREPARE_BAM(input_files, params.reference)

    if (!params.skip_deduplication) {
        MARK_DUPLICATES(PREPARE_BAM.out.prepared_bams)
        deduplicated_bams = MARK_DUPLICATES.out.deduplicated_bams
    }
    else {
        INDEX_BAM(PREPARE_BAM.out.prepared_bams)
        deduplicated_bams = INDEX_BAM.out.indexed_bams
    }

    if (params.split_cigarn) {
        SPLIT_CIGAR_N_READS(deduplicated_bams, params.reference)
        deduplicated_bams = SPLIT_CIGAR_N_READS.out.split_cigarn_bams
    }

    if (! params.skip_metrics) {
        if (params.intervals) {
            HS_METRICS(deduplicated_bams)
        }
        METRICS(deduplicated_bams, params.reference)
        COVERAGE_ANALYSIS(deduplicated_bams)
        FLAGSTAT(deduplicated_bams)
    }

    if (!params.skip_realignment) {
        REALIGNMENT_AROUND_INDELS(deduplicated_bams, params.reference)
        realigned_bams = REALIGNMENT_AROUND_INDELS.out.realigned_bams
    }
    else {
        realigned_bams = deduplicated_bams
    }

    if (!params.skip_bqsr) {
        BQSR(realigned_bams, params.reference)
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
