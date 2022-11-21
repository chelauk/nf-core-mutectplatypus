process PLATYPUS_FILTER {

    tag "$patient"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    tuple val(patient), path(vcf), val(norm)
    val (tef)

    output:
    tuple val(patient), path("*filtered.vcf"), emit: vcf
    tuple val(patient), path("*removed.vcf"),  optional:true, emit: rejected
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${patient}_platypus"
    def my_vcf = "${vcf.toString().minus(".gz")}"
    """
    gunzip $vcf
    filter_platypus.py $my_vcf ${norm[0]} $tef
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    """
    echo "filter_platypus.py $vcf ${norm[0]} $tef"
    touch "$patient"_platypus_filtered.vcf
    touch versions.yml
    """
}
