params.mark_duplicates_cpus = 2
params.mark_duplicates_memory = "16g"
params.remove_duplicates = true
params.skip_metrics = false
params.output = 'output'


process MARK_DUPLICATES {
    cpus "${params.mark_duplicates_cpus}"
    memory "${params.mark_duplicates_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/mark_duplicates", mode: "copy", pattern: "*.dedup_metrics.txt"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(type), file(bam)

    output:
    tuple val(name), val(type), file("${name}.dedup.bam"), file("${name}.dedup.bam.bai"), emit: deduplicated_bams
    file("${name}.dedup_metrics.txt") optional true
    file("software_versions.${task.process}.txt")

    script:
    dedup_metrics = params.skip_metrics ? "": "--METRICS_FILE ${name}.dedup_metrics.txt"
    remove_duplicates = params.remove_duplicates ? "--REMOVE_DUPLICATES true" : "--REMOVE_DUPLICATES false"
    """
    mkdir tmp

    gatk SortSam \
    --INPUT ${bam} \
    --OUTPUT ${name}.sorted.bam \
    --SORT_ORDER coordinate

    gatk MarkDuplicates \
    --java-options '-Xmx${params.mark_duplicates_memory}  -Djava.io.tmpdir=./tmp' \
    --INPUT  ${name}.sorted.bam \
    --OUTPUT ${name}.dedup.bam \
    --ASSUME_SORT_ORDER coordinate \
    --CREATE_INDEX true ${remove_duplicates} ${dedup_metrics}

    cp ${name}.dedup.bai ${name}.dedup.bam.bai

    rm -f ${name}.sorted.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    """
}
