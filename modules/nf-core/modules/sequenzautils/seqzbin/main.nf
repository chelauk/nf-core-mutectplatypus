process SEQUENZAUTILS_BINNING {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h67e14b5_5' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39h67e14b5_5' }"

    input:
    tuple val(patient), val(id), path(concat_seqz)
    val(bin)

    output:
    tuple val(patient), val(id), path("*bin${bin}.seqz.gz"), emit: seqz_bin
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    sequenza-utils \\
        seqz_binning \\
        $args \\
        -w $concat_seqz \\
        -o - |\\
        gzip > ${id}_bin${bin}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    cat << 'EOF'
    sequenza-utils \\
        seqz_binning \\
        $args \\
        -w $concat_seqz \\
        -o - |\\
        gzip > ${id}_bin${bin}.seqz.gz
    EOF
    touch ${id}_bin${bin}.seqz.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
