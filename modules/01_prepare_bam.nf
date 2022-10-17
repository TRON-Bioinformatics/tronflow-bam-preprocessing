params.prepare_bam_cpus = 3
params.prepare_bam_memory = "8g"
params.index_cpus = 1
params.index_memory = "8g"
params.platform = "ILLUMINA"
params.output = 'output'

/*
This step sets MAPQ to 0 for all unmapped reads + avoids soft clipping beyond the end of the reference genome
This step reorders chromosomes in the BAM file according to the provided reference (this step is required for GATK)
Adds the required read groups fields to the BAM file. The provided type is added to the BAM sample name.
*/
process PREPARE_BAM {
    cpus "${params.prepare_bam_cpus}"
    memory "${params.prepare_bam_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0 bioconda::samtools=1.12" : null)

    input:
    tuple val(name), val(type), file(bam)
    val(reference)

    output:
    tuple val(name), val(type), file("${name}.prepared.bam"), emit: prepared_bams
    file("software_versions.${task.process}.txt")

    script:
    """
    mkdir tmp

    gatk AddOrReplaceReadGroups \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=./tmp' \
    --VALIDATION_STRINGENCY SILENT \
    --INPUT ${bam} \
    --OUTPUT /dev/stdout \
    --REFERENCE_SEQUENCE ${reference} \
    --RGPU 1 \
    --RGID 1 \
    --RGSM ${type} \
    --RGLB 1 \
    --RGPL ${params.platform} | \
    gatk CleanSam \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=./tmp' \
    --INPUT /dev/stdin \
    --OUTPUT /dev/stdout | \
    gatk ReorderSam \
    --java-options '-Xmx${params.prepare_bam_memory} -Djava.io.tmpdir=./tmp' \
    --INPUT /dev/stdin \
    --OUTPUT ${name}.prepared.bam \
    --SEQUENCE_DICTIONARY ${reference}

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}

process INDEX_BAM {
    cpus "${params.index_cpus}"
    memory "${params.index_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(type), file(bam)

    output:
    tuple val(name), val(type), file("${bam}"), file("${bam.baseName}.bai"), emit: indexed_bams
    file("software_versions.${task.process}.txt")

    script:
    """
    mkdir tmp

    gatk BuildBamIndex \
    --java-options '-Xmx8g  -Djava.io.tmpdir=./tmp' \
    --INPUT  ${bam}

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    """
}
