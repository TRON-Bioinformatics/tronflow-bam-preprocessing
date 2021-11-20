params.bqsr_cpus = 3
params.bqsr_memory = "4g"
params.reference = false
params.dbsnp = false
params.output = 'output'


process BQSR {
    cpus "${params.bqsr_cpus}"
    memory "${params.bqsr_memory}"
    publishDir "${params.output}/${name}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), val(bam_name), val(type), file(bam), file(bai)

    output:
    tuple val("${name}"), val("${type}"), val("${params.output}/${name}/${bam_name}.preprocessed.bam"), emit: recalibrated_bams
    file "${bam_name}.recalibration_report.grp"
    file "${bam_name}.preprocessed.bam"
    file "${bam_name}.preprocessed.bai"

    """
    mkdir tmp

    gatk BaseRecalibrator \
    --java-options '-Xmx${params.bqsr_memory} -Djava.io.tmpdir=tmp' \
    --input ${bam} \
    --output ${bam_name}.recalibration_report.grp \
    --reference ${params.reference} \
    --known-sites ${params.dbsnp}

    gatk ApplyBQSR \
    --java-options '-Xmx${params.bqsr_memory} -Djava.io.tmpdir=tmp' \
    --input ${bam} \
    --output ${bam_name}.preprocessed.bam \
    --reference ${params.reference} \
    --bqsr-recal-file ${bam_name}.recalibration_report.grp
    """
}


process CREATE_OUTPUT {
    cpus 1
    memory '1g'
    publishDir "${params.output}/${name}", mode: "copy"
    tag "${name}"

    input:
    tuple val(name), val(bam_name), val(type), file(bam), file(bai)

    output:
    tuple val("${name}"), val("${type}"), val("${params.output}/${name}/${bam_name}.preprocessed.bam"), emit: recalibrated_bams
    file "${bam_name}.preprocessed.bam"
    file "${bam_name}.preprocessed.bai"

    """
    cp ${bam} ${bam_name}.preprocessed.bam
    cp ${bai} ${bam_name}.preprocessed.bai
    """
}