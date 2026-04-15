
process CREATE_FAIDX {
    cpus "1"
    memory "4g"
    tag "${reference}"

    conda (params.enable_conda ? "bioconda::samtools=${params.samtools_version}" : null)

    input:
    val(reference)

    script:
    """
    samtools faidx ${reference}
    """
}

process CREATE_DICT {
    cpus "1"
    memory "4g"
    tag "${reference}"

    conda (params.enable_conda ? "bioconda::gatk4=${params.gatk4_version}" : null)

    input:
    val(reference)

    script:
    """
    gatk CreateSequenceDictionary --REFERENCE ${reference}
    """
}