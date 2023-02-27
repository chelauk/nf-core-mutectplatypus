process EVOVERSE_CNAQC {
    errorStrategy 'ignore'
    tag "$patient"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container 'r-evoverse.sif'

    input:
    tuple val(patient), val(id), path(segments), path(vcf), path(tbi)
    val(ploidy)
    path(drivers)

    output:
    tuple val(patient), val(id), path("*pdf"), path("*rds"), emit: evoverse_pdf
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${patient}_${id}"
    """
    evoverse.R $prefix $id $vcf $drivers
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    def prefix = task.ext.prefix ?: "${patient}_${id}"
    """
    echo "evoverse $prefix $id $vcf $drivers"
    touch ${id}_${ploidy}.pdf
    touch ${id}_${ploidy}.rds
    touch versions.yml
    """
}
