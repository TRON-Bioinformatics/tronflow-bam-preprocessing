params.metrics_cpus = 1
params.metrics_memory = "8g"
params.collect_hs_metrics_min_base_quality = false
params.collect_hs_metrics_min_mapping_quality = false
params.output = 'output'
params.intervals = false


process HS_METRICS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/hs_metrics", mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true
    file("${name}.hs_metrics.txt")
    file("software_versions.${task.process}.txt")

    script:
    minimum_base_quality = params.collect_hs_metrics_min_base_quality ?
        "--MINIMUM_BASE_QUALITY ${params.collect_hs_metrics_min_base_quality}" : ""
    minimum_mapping_quality = params.collect_hs_metrics_min_mapping_quality ?
        "--MINIMUM_MAPPING_QUALITY ${params.collect_hs_metrics_min_mapping_quality}" : ""
    """
    mkdir tmp

    gatk BedToIntervalList \
    --INPUT ${params.intervals} \
    --OUTPUT my.intervals \
    --SEQUENCE_DICTIONARY ${bam}

    gatk CollectHsMetrics \
    --java-options '-Xmx${params.metrics_memory}  -Djava.io.tmpdir=./tmp' \
    --INPUT  ${bam} \
    --OUTPUT ${name}.hs_metrics.txt \
    --TARGET_INTERVALS my.intervals \
    --BAIT_INTERVALS my.intervals \
    ${minimum_base_quality} ${minimum_mapping_quality}

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    """
}

process METRICS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/gatk_multiple_metrics", mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    // NOTE: the method CollectMultipleMetrics has a hidden dependency to R for making plots
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0 r::r=3.6.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)
    val(reference)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true
    file("software_versions.${task.process}.txt")

    """
    mkdir tmp

    gatk CollectMultipleMetrics \
    --java-options '-Xmx${params.metrics_memory}  -Djava.io.tmpdir=./tmp' \
    --INPUT  ${bam} \
    --OUTPUT ${name} \
    --REFERENCE_SEQUENCE ${reference} \
    --PROGRAM QualityScoreDistribution \
    --PROGRAM MeanQualityByCycle \
    --PROGRAM CollectAlignmentSummaryMetrics \
    --PROGRAM CollectBaseDistributionByCycle \
    --PROGRAM CollectGcBiasMetrics \
    --PROGRAM CollectInsertSizeMetrics \
    --PROGRAM CollectSequencingArtifactMetrics \
    --PROGRAM CollectSequencingArtifactMetrics

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    """
}

process COVERAGE_ANALYSIS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/coverage", mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
        tuple val(name), val(type), file(bam), file(bai)

    output:
        file("${name}.coverage.tsv")
        file("${name}.depth.tsv")
        file("software_versions.${task.process}.txt")

    script:
    minimum_base_quality = params.collect_hs_metrics_min_base_quality ?
        "--min-BQ ${params.collect_hs_metrics_min_base_quality}" : ""
    minimum_mapping_quality = params.collect_hs_metrics_min_mapping_quality ?
        "--min-MQ ${params.collect_hs_metrics_min_mapping_quality}" : ""
    intervals = params.intervals ? "-b ${params.intervals}" : ""
    """
    samtools coverage ${minimum_base_quality} ${minimum_mapping_quality} ${bam} > ${name}.coverage.tsv
    samtools depth -s -d 0 -H ${intervals} ${bam} > ${name}.depth.tsv

    echo ${params.manifest} >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}
