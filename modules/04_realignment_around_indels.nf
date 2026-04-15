process REALIGNMENT_AROUND_INDELS {
    cpus "${params.realignment_around_indels_cpus}"
    memory "${params.realignment_around_indels_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/metrics/realignment", mode: "copy", pattern: "*.RA.intervals"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    // NOTE: this dependency is fixed to GATK 3 as the realignment around indels is not anymore maintained in GATK 4
    // but still for some reason for GATK 3 to work the dependency to GATK 4.2.0.0 is needed
    conda (params.enable_conda ? "bioconda::gatk4=${params.gatk4_realignment_version} bioconda::gatk=${params.gatk3_version}" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)
    val(reference)

    output:
    tuple val(name), val(type), file("${name}.realigned.bam"), file("${name}.realigned.bai"), emit: realigned_bams
    path("${name}.RA.intervals")
    path("software_versions.${task.process}.txt")

    script:
    known_indels1 = params.known_indels1 ? " --known ${params.known_indels1}" : ""
    known_indels2 = params.known_indels2 ? " --known ${params.known_indels2}" : ""
    known_alleles1 = params.known_indels1 ? " --knownAlleles ${params.known_indels1}" : ""
    known_alleles2 = params.known_indels2 ? " --knownAlleles ${params.known_indels2}" : ""
    """
    mkdir tmp

    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=./tmp -T RealignerTargetCreator \
    --input_file ${bam} \
    --out ${name}.RA.intervals \
    --reference_sequence ${reference} ${known_indels1} ${known_indels2}

    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=./tmp -T IndelRealigner \
    --input_file ${bam} \
    --out ${name}.realigned.bam \
    --reference_sequence ${reference} \
    --targetIntervals ${name}.RA.intervals \
    --consensusDeterminationModel USE_SW \
    --LODThresholdForCleaning 0.4 \
    --maxReadsInMemory 600000 ${known_alleles1} ${known_alleles2}

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk3 --version >> software_versions.${task.process}.txt
    """
}
