
process CREATE_FAIDX {
    cpus "1"
    memory "4g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
    val(reference)

    """
    samtools faidx ${reference}
    """
}

process CREATE_DICT {
    cpus "1"
    memory "4g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    val(reference)

    """
    gatk CreateSequenceDictionary --REFERENCE ${reference}
    """
}