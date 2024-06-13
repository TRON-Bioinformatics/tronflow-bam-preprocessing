params.mark_duplicates_cpus = 2
params.mark_duplicates_memory = "16g"
params.split_reads_cpus = 2
params.split_reads_memory = "4g"
params.remove_duplicates = true
params.output = 'output'
params.split_cigarn_args = ""


process MARK_DUPLICATES {
    cpus "${params.mark_duplicates_cpus}"
    memory "${params.mark_duplicates_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::sambamba=0.8.2" : null)

    input:
    tuple val(name), val(type), file(bam)

    output:
    tuple val(name), val(type), file("${name}.dedup.bam"), file("${name}.dedup.bam.bai"), emit: deduplicated_bams
    file("software_versions.${task.process}.txt")

    script:
    remove_duplicates_param = params.remove_duplicates ? "--remove-duplicates" : ""
    """
    mkdir tmp

    # sort
    sambamba sort \
        --nthreads=${task.cpus} \
        --tmpdir=./tmp \
        --out=${name}.sorted.bam \
        ${bam}

    # removes duplicates (sorted from the alignment process)
    sambamba markdup ${remove_duplicates_param} \
        --nthreads=${task.cpus} \
        --tmpdir=./tmp \
        ${name}.sorted.bam ${name}.dedup.bam

    rm -f ${name}.sorted.bam

    # indexes the output BAM file
    sambamba index \
        --nthreads=${task.cpus} \
        ${name}.dedup.bam ${name}.dedup.bam.bai

    echo ${params.manifest} >> software_versions.${task.process}.txt
    sambamba --version >> software_versions.${task.process}.txt
    """
}

process SPLIT_CIGAR_N_READS {
    cpus "${params.split_reads_cpus}"
    memory "${params.split_reads_memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(type), file(bam), file(bai)
    val(reference)

    output:
    tuple val(name), val(type), file("${name}.split_cigarn.bam"), file("${name}.split_cigarn.bam.bai"), emit: split_cigarn_bams
    file("software_versions.${task.process}.txt")

    script:
    """
    mkdir tmp

    gatk SplitNCigarReads \
    --java-options '-Xmx${params.split_reads_memory}  -Djava.io.tmpdir=./tmp' \
    --input ${bam} \
    --output ${name}.split_cigarn.bam \
    --create-output-bam-index true \
    --reference ${reference} \
    ${params.split_cigarn_args}

    cp ${name}.split_cigarn.bai ${name}.split_cigarn.bam.bai

    echo ${params.manifest} >> software_versions.${task.process}.txt
    gatk --version >> software_versions.${task.process}.txt
    """
}
