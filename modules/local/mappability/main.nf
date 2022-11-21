process MAPPABILITY {

    tag "$patient"

    label 'process_low'
    
    container 'mappability.sif'

    input:
        tuple val(patient),  path(vcf), path(tbi)
        path mappability_bw
        val pan

    output:
        path "*out"

    script:
    """
    R --slave --no-restore --file=${baseDir}/scripts/assess_mutation_mappability.R \
	--args $patient $vcf $mappability_bw $pan
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
    """

    stub:
    """
    echo "assess_mutation_mappability.R $patient $vcf $mappability_bw $pan "
    touch ${patient}.out
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
    """
}
