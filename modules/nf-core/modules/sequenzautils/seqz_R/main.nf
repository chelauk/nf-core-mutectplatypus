process SEQUENZAUTILS_RSEQZ {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::r-sequenza=3.0.0" : null)
    container 'r-sequenza.sif'

    input:
    tuple val(patient), val(id), val(gender), path(seqz_bin)
    val  ploidy
    val  ccf

    output:
    tuple val(patient), val(id), path("${id}_${ploidy}"), emit: rseqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}_${ploidy}"
    """
    zcat ${seqz_bin} > $seqz_in
    analyse_cn_sequenza.R ${seqz_in} ${prefix} ${gender} ${ploidy} ${ccf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}_${ploidy}"
    """
    echo "analyse_cn_sequenza.R ${seqz_in} ${prefix} ${gender} ${ploidy} ${ccf}"
    mkdir ${id}_${ploidy}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
