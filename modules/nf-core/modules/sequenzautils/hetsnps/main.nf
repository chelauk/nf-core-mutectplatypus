process SEQUENZAUTILS_HETSNPS {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h67e14b5_5' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39h67e14b5_5' }"

    input:
    tuple val(patient), val(id), val(gender), path(concat_seqz)

    output:
    tuple val(patient), val(id), val(gender), path("*het.seqz"), emit: het_seqz
    path "versions.yml"                         , emit: versionsi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
	zcat ${concat_seqz} > temp
    head -1 temp > header 
	sed -n '/het/p' temp > temp2
	cat header temp2 > ${prefix}_het.seqz
    rm temp*

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
    zcat $concat_seqz | grep -ae "chromo|het" > ${prefix}_het.seqz
    EOF
    touch ${prefix}_het.seqz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
