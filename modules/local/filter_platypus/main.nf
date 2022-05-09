process PLATYPUS_FILTER {

    tag "$patient"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(patient), path(vcf), val(norm)

    output:
    tuple val(patient), path("*filtered.vcf"), emit: vcf
    tuple val(patient), path("*removed.vcf"),  optional:true, emit: rejected
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${patient}_platypus"
    """
    filter_platypus.py $vcf ${norm[0]} $prefix
    touch versions.yml
	"""
    stub:
    """
    echo "filter_platypus.py $vcf ${norm[0]}"
    touch "$patient"_platypus_filtered.vcf
    touch versions.yml
    """
}
