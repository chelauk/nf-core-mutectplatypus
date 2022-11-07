process RUN_BEAGLEBERG {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container 'beagleberg.img'

    input:
    tuple val(patient), val(id), path(het_seqz), path(lrr_rds), path(phased_vcfs)
    path(chr_arm_boundaries)

    output:
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    run_beagleberg.R
    """

	stub:
	def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo "run_beagleberg.R  ${id} ${chr_arm_boundaries}"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preprocess_het_snps: 3.0.0
    END_VERSIONS
    """
}


