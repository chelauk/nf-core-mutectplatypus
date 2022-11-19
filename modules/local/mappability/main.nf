process MAPPABILITY {

    tag "$patient"

    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcfr" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-vcfr:1.8.0--r36h0357c0b_3' :
        'quay.io/biocontainers/r-vcfr:1.8.0--r36h0357c0b_3' }"

    input:
        tuple val(patient),  path(vcf), path(tbi)
        path mappability_bw
        val pan

    output:
        path "*out"

    script:
    """
    assess_mutation_mappability.R $patient $vcf $mappability_bw $pan
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
