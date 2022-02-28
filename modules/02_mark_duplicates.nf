params.mark_duplicates_cpus = 16
params.mark_duplicates_memory = "64g"
params.remove_duplicates = true
params.skip_metrics = false
params.output = 'output'


process MARK_DUPLICATES {
    cpus "${params.mark_duplicates_cpus}"
    memory "${params.mark_duplicates_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/mark_duplicates", mode: "copy", pattern: "*.dedup_metrics.txt"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(type), file(bam)

    output:
    tuple val(name), val(type), file("${name}.dedup.bam"), file("${name}.dedup.bam.bai"), emit: deduplicated_bams
    file("${name}.dedup_metrics.txt") optional true

    script:
    dedup_metrics = params.skip_metrics ? "": "--metrics-file ${name}.dedup_metrics.txt"
    remove_duplicates = params.remove_duplicates ? "--remove-all-duplicates true" : "--remove-all-duplicates false"
    """
    mkdir tmp

    gatk MarkDuplicatesSpark \
    --java-options '-Xmx${params.mark_duplicates_memory}  -Djava.io.tmpdir=tmp' \
    --input  ${bam} \
    --output ${name}.dedup.bam \
    --conf 'spark.executor.cores=${task.cpus}' ${remove_duplicates} ${dedup_metrics}
    """
}
