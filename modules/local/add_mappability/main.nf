process ADD_MAPPABILITY {
    cache false
    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    tuple path(mappability), val(patient), val(id), path(maf)

    output:
    tuple val(patient), val(id), path("*adj.maf"), emit: fixed_maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    add_mappability.py $mappability $maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    """
    echo "add_mappability.py $mappability $maf"
    touch "$id".adj.maf
    touch versions.yml
    """
}
