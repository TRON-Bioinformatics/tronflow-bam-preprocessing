process HS_METRICS {
    cpus params.metrics_cpus
    memory params.metrics_memory
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/hs_metrics", mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=${params.gatk4_version}" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)

    output:
    path("*_metrics", optional: true)
    path("*.pdf", optional: true)
    path("${name}.hs_metrics.txt")
    path("software_versions.${task.process}.txt")

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
    conda (params.enable_conda ? "bioconda::gatk4=${params.gatk4_version} r::r=3.6.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)
    val(reference)

    output:
    path("*_metrics", optional: true)
    path("*.pdf", optional: true)
    path("software_versions.${task.process}.txt")

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

    conda (params.enable_conda ? "bioconda::samtools=${params.samtools_version}" : null)

    input:
        tuple val(name), val(type), file(bam), file(bai)

    output:
        path("${name}.coverage.tsv")
        path("${name}.depth.tsv")
        path("software_versions.${task.process}.txt")

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

process FLAGSTAT {
    cpus "${params.metrics_cpus}"
    memory "${params.metrics_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/flagstat", mode: "copy", pattern: "*.flagstat.csv"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::sambamba=${params.sambamba_version}" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)

    output:
    path("${name}.flagstat.csv")
    path("software_versions.${task.process}.txt")

    script:
    """
    sambamba flagstat \
        --nthreads=${task.cpus} \
        --tabular \
        ${bam} > ${name}.flagstat.csv

    echo ${params.manifest} >> software_versions.${task.process}.txt
    sambamba --version >> software_versions.${task.process}.txt
    """
}
