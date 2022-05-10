process PLATYPUS_ZIP {

    tag "${patient}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/samtools:1.15--h3843a85_0' :
        'quay.io/biocontainers/samtools:1.15--0' }"
    input:
    tuple val(patient), path(vcf)

    output:
    tuple val(patient), path("*.vcf.gz"), path("*vcf.gz.tbi"), emit: vcf
    path "*versions.yml",                                      emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bgzip $vcf
    tabix -p vcf "${vcf}".gz
    touch versions.yml
	"""
    stub:
    def args = task.ext.args ?: ''
    """
    echo -e "bgzip $vcf
    tabix -p vcf "${vcf}".gz"
    touch "${vcf}".gz
    touch "${vcf}".gz.tbi
    touch versions.yml
    """
}
