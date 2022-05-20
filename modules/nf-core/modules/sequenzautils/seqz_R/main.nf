process SEQUENZAUTILS_RSEQZ {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::r-sequenza=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sequenza:3.0.0--r41h3342da4_4' :
        'quay.io/biocontainers/r-sequenza:3.0.0--r41h3342da4_4' }"

    input:
    tuple val(patient), val(id), path(seqz_bin)
    val  gender
    val  ploidy

    output:
    tuple val(patient), val(id), path(id), emit: rseqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}_${ploidy}"
    """
    analyse_cn_sequenza.R ${seqz_bin} ${prefix} ${gender} ${ploidy}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}_${ploidy}"
    """
    echo "analyse_cn_sequenza.R ${seqz_bin} ${prefix} ${gender} ${ploidy}"
    mkdir $id
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
