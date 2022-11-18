process MAPPABILITY {
    errorStrategy 'ignore'
    tag "$patient"

    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcfr" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-vcfr:1.8.0--r36h0357c0b_3' :
        'quay.io/biocontainers/r-vcfr:1.8.0--r36h0357c0b_3' }"

    input:
        vcf
        mappability_bigwig
        variable
    output:
        path "*out"

    script:
    """
    assess_mutation_mappability.R
    """


    script:
    """
    evoverse.R $id $ploidy $vcf $drivers
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        1
    END_VERSIONS
	"""
    stub:
    """
    echo "evoverse $id $ploidy $segments $vcf $drivers"
    touch ${id}_${ploidy}.pdf
    touch ${id}_${ploidy}.rds
    touch versions.yml
    """
}
