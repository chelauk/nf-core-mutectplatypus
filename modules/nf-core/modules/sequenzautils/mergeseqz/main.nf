process SEQUENZAUTILS_MERGESEQZ {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"
    
    input:
    tuple val(patient), val(id), path(seqz)

    output:
    tuple val(patient), val(id), path("*concat.seqz.gz"), emit: concat_seqz
    tuple val(patient), val(id), path("*concat.seqz.gz.tbi"), emit: concat_seqz_tbi
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    zcat $seqz | \
    awk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' | bgzip > \
    ${prefix}.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 ${prefix}_concat.seqz.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    cat << 'EOF' 
    zcat $seqz | \
    awk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' | bgzip > \
    ${prefix}_concat.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 ${prefix}_concat.seqz.gz
    EOF
    touch ${prefix}_concat.seqz.gz
    touch ${prefix}_concat.seqz.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
