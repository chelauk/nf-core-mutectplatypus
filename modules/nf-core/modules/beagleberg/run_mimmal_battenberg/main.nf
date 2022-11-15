process MIMMAL_BATTENBERG {
    tag "$id"
    label 'process_medium'

    container 'beagleberg.sif'

    input:
    tuple val(patient), val(id), val(gender), path(baf), path(lrr)
    path(chr_arm_boundaries)

    output:
    path ("./mimmal_battenberg/*")                    , emit: mimmal_out
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    run_battenberg_using_mimmal.R \
	${baseDir}/bin/run_battenberg_using_mimmal_src.R mimmal_battenberg  \
	${baf} \
	${lrr} \
	${chr_arm_boundaries} \
	${id}
    
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      beagleberg: 1.0.0
    END_VERSIONS
    """

	stub:
	def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo runBattenbergUsingMiMMAl.R ${patient} ${id} ${gender} ${baf} ${lrr} ${chr_arm_boundaries}
    mkdir mimmal_battenberg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preprocess_het_snps: 3.0.0
    END_VERSIONS
    """
}


