params.metrics_cpus = 1
params.metrics_memory = "8g"
params.collect_hs_metrics_min_base_quality = false
params.collect_hs_metrics_min_mapping_quality = false
params.reference = false
params.output = 'output'


process HS_METRICS {
    cpus "${params.metrics_cpus}"
    memory "${params.metrics_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics", mode: "copy"

    input:
    tuple val(name), val(bam_name), val(type), file(bam), file(bai)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true
    file("${bam.baseName}.hs_metrics.txt")

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
    --OUTPUT ${bam.baseName}.hs_metrics.txt \
    --TARGET_INTERVALS ${params.intervals} \
    --BAIT_INTERVALS ${params.intervals} \
    ${minimum_base_quality} ${minimum_mapping_quality}
    """
}

process METRICS {
    cpus "${params.metrics_cpus}"
    memory "${params.metrics_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics", mode: "copy"

    input:
    tuple val(name), val(bam_name), val(type), file(bam), file(bai)

    output:
    file("*_metrics") optional true
    file("*.pdf") optional true

    """
    mkdir tmp

    gatk CollectMultipleMetrics \
    --java-options '-Xmx${params.metrics_memory}  -Djava.io.tmpdir=tmp' \
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