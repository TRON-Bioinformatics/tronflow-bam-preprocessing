params.realignment_around_indels_cpus = 2
params.realignment_around_indels_memory = "31g"
params.known_indels1 = false
params.known_indels2 = false
params.reference = false
params.output = 'output'


process REALIGNMENT_AROUND_INDELS {
    cpus "${params.realignment_around_indels_cpus}"
    memory "${params.realignment_around_indels_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "*.RA.intervals"

    // NOTE: this dependency is fixed to GATK 3 as the realignment around indels is not anymore maintained in GATK 4
    // but still for some reason for GATK 3 to work the dependency to GATK 4 is needed
    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0 bioconda::gatk=3.8" : null)

    input:
    tuple val(name), val(bam_name), val(type), file(bam), file(bai)

    output:
    tuple val(name), val(bam_name), val(type), file("${bam.baseName}.realigned.bam"), file("${bam.baseName}.realigned.bai"), emit: realigned_bams
    file("${bam.baseName}.RA.intervals")

    script:
    known_indels1 = params.known_indels1 ? " --known ${params.known_indels1}" : ""
    known_indels2 = params.known_indels2 ? " --known ${params.known_indels2}" : ""
    known_alleles1 = params.known_indels1 ? " --knownAlleles ${params.known_indels1}" : ""
    known_alleles2 = params.known_indels2 ? " --knownAlleles ${params.known_indels2}" : ""
    """
    mkdir tmp

    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T RealignerTargetCreator \
    --input_file ${bam} \
    --out ${bam.baseName}.RA.intervals \
    --reference_sequence ${params.reference} ${known_indels1} ${known_indels2}

    gatk3 -Xmx${params.realignment_around_indels_memory} -Djava.io.tmpdir=tmp -T IndelRealigner \
    --input_file ${bam} \
    --out ${bam.baseName}.realigned.bam \
    --reference_sequence ${params.reference} \
    --targetIntervals ${bam.baseName}.RA.intervals \
    --consensusDeterminationModel USE_SW \
    --LODThresholdForCleaning 0.4 \
    --maxReadsInMemory 600000 ${known_alleles1} ${known_alleles2}
    """
}
