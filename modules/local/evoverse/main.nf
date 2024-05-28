process EVOVERSE_CNAQC {
    errorStrategy 'ignore'
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container 'r-evoverse.sif'

    input:
    tuple val(meta), path(vcf), val(tissue), path(segments)
    val(ploidy)
    path(drivers)

    output:
    tuple val(meta), path("*pdf"), path("*rds"), emit: evoverse_pdf
    path "*versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def vcf_file = "${vcf[0]}"
    def prefix = task.ext.prefix ?: "${segments}"
    """
    zgrep -P "^#|PASS" $vcf_file | gzip > temp.vcf.gz
    evoverse.R $prefix ${meta.id} temp.vcf.gz $drivers
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    def prefix = task.ext.prefix ?: "${segments}"
    """
    echo "evoverse $prefix ${meta.id} $vcf $drivers"
    touch ${meta.id}_${ploidy}.pdf
    touch ${meta.id}_${ploidy}.rds
    touch versions.yml
    """
}
