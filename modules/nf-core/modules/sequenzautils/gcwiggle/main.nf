process SEQUENZAUTILS_GCWIGGLE {

    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"

    input:
    path(fasta)

    output:
    path("*.wig.gz"), emit: wiggle
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.toString().replaceFirst(~/\.[^\.]+$/, '')}"
    """
    sequenza-utils \\
        gc_wiggle \\
        $args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.toString().replaceFirst(~/\.[^\.]+$/, '')}"
    """
    echo -e " sequenza-utils \\
        gc_wiggle \\
        $args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz
        "

    touch ${prefix}.wig.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
