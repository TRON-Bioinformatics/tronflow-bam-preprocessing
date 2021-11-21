params.metrics_cpus = 1
params.metrics_memory = "8g"
params.collect_hs_metrics_min_base_quality = false
params.collect_hs_metrics_min_mapping_quality = false
params.reference = false
params.output = 'output'
params.intervals_bed = false
params.intervals = false


process HS_METRICS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/hs_metrics", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true
    file("${name}.hs_metrics.txt")

    script:
    minimum_base_quality = params.collect_hs_metrics_min_base_quality ?
        "--MINIMUM_BASE_QUALITY ${params.collect_hs_metrics_min_base_quality}" : ""
    minimum_mapping_quality = params.collect_hs_metrics_min_mapping_quality ?
        "--MINIMUM_MAPPING_QUALITY ${params.collect_hs_metrics_min_mapping_quality}" : ""
    """
    mkdir tmp

    gatk CollectHsMetrics \
    --java-options '-Xmx${params.metrics_memory}  -Djava.io.tmpdir=tmp' \
    --INPUT  ${bam} \
    --OUTPUT ${name}.hs_metrics.txt \
    --TARGET_INTERVALS ${params.intervals} \
    --BAIT_INTERVALS ${params.intervals} \
    ${minimum_base_quality} ${minimum_mapping_quality}
    """
}

process METRICS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/gatk_multiple_metrics", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true

    """
    mkdir tmp

    gatk CollectMultipleMetrics \
    --java-options '-Xmx${params.metrics_memory}  -Djava.io.tmpdir=tmp' \
    --INPUT  ${bam} \
    --OUTPUT ${name} \
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

process COVERAGE_ANALYSIS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/coverage", mode: "copy"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
        tuple val(name), val(type), file(bam), file(bai)

    output:
        file("${name}.coverage.tsv")
        file("${name}.depth.tsv")

    script:
    minimum_base_quality = params.collect_hs_metrics_min_base_quality ?
        "--min-BQ ${params.collect_hs_metrics_min_base_quality}" : ""
    minimum_mapping_quality = params.collect_hs_metrics_min_mapping_quality ?
        "--min-MQ ${params.collect_hs_metrics_min_mapping_quality}" : ""
    intervals = params.intervals_bed ? "-b ${params.intervals_bed}" : ""
    """
    samtools coverage ${minimum_base_quality} ${minimum_mapping_quality} ${bam} > ${name}.coverage.tsv
    samtools depth -s -d 0 -H ${intervals} ${bam} > ${name}.depth.tsv
    """
}
