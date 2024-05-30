process PLATYPUS_FILTER {

    tag "${meta_control.patient}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta_control), val(meta_tumour), path(vcf)
    val (tef)

    output:
    tuple val(meta_control), val(meta_tumour), path("*filtered.vcf"), emit: vcf
    tuple val(meta_control), path("*removed.vcf"),  optional:true, emit: rejected
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta_control.patient}_platypus"
    def my_vcf = "${vcf.toString().minus(".gz")}"
    """
    gunzip $vcf
    filter_platypus.py $my_vcf ${meta_control.id} $tef
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    """
    echo "filter_platypus.py $vcf ${meta_control.id} $tef"
    touch "${meta_control.patient}"_platypus_filtered.vcf
    touch versions.yml
    """
}
