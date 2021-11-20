params.mark_duplicates_cpus = 16
params.mark_duplicates_memory = "64g"
params.remove_duplicates = true
params.skip_metrics = false
params.output = 'output'


process MARK_DUPLICATES {
    cpus "${params.mark_duplicates_cpus}"
    memory "${params.mark_duplicates_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics", mode: "copy", pattern: "*.dedup_metrics.txt"

    input:
    tuple val(name), val(bam_name), val(type), file(bam)

    output:
    tuple val(name), val(bam_name), val(type),
        file("${bam.baseName}.dedup.bam"), file("${bam.baseName}.dedup.bam.bai"), emit: deduplicated_bams
    file("${bam.baseName}.dedup_metrics.txt") optional true

    script:
    dedup_metrics = params.skip_metrics ? "": "--metrics-file ${bam.baseName}.dedup_metrics.txt"
    remove_duplicates = params.remove_duplicates ? "--remove-all-duplicates true" : "--remove-all-duplicates false"
    """
    mkdir tmp

    gatk MarkDuplicatesSpark \
    --java-options '-Xmx${params.mark_duplicates_memory}  -Djava.io.tmpdir=tmp' \
    --input  ${bam} \
    --output ${bam.baseName}.dedup.bam \
    --conf 'spark.executor.cores=${task.cpus}' ${remove_duplicates} ${dedup_metrics}
    """
}