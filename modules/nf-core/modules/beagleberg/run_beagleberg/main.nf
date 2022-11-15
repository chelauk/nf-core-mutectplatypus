process RUN_BEAGLEBERG {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container 'beagleberg.sif'

    input:
    tuple val(patient), val(id), path(het_seqz), path(lrr_rds), path(phased_vcfs)
    path(chr_arm_boundaries)

    output:
	path "beagleberg/*"                         , emit: beagleberg_out
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    run_beagleberg.R \
	${baseDir}/bin/run_battenberg_using_mimmal_src.R \
	beagleberg \
    ${lrr_rds} \
	${chr_arm_boundaries} \
	${id} \
	${het_seqz} 
    
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagleberg: 1.0.0
    END_VERSIONS
    """

	stub:
	def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo "run_beagleberg.R  ${id} ${chr_arm_boundaries}"
    mkdir beagleberg
    cd beagleberg
    for i in {1..22}
        do
        touch ${id}_subclones_chr"\$i".png
        done

    cd ..
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preprocess_het_snps: 3.0.0
    END_VERSIONS
    """
}


