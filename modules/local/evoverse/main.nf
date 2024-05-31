process EVOVERSE_CNAQC {
    debug true
    errorStrategy 'ignore'
    tag "${meta.id}_${caller}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container 'r-evoverse.sif'

    input:
    tuple val(meta), path(vcf), val(tissue), path(segments)
    val(ploidy)
    path(drivers)
    val(coverage)
    val(caller)

    output:
    tuple val(meta), path("*pdf"), path("*rds"), emit: evoverse_pdf
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def vcf_file = "${vcf[0]}"
    def prefix = task.ext.prefix ?: "${segments}"
    """
    evoverse.R $prefix ${meta.id} $vcf_file $drivers $coverage $caller
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    def prefix = task.ext.prefix ?: "${segments}"
    """
    echo "evoverse $prefix ${meta.id} $vcf $drivers $coverage $caller"
    touch ${meta.id}_${ploidy}_${caller}.pdf
    touch ${meta.id}_${ploidy}.${caller}.rds
    touch versions.yml
    """
}
