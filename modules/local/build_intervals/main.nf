process BUILD_INTERVALS {
    tag "$fai"
    label 'process_low'

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gawk:5.1.0"
    } else {
        container "quay.io/biocontainers/gawk:5.1.0"
    }

    input:
    path fai

    output:
    path "${fai.baseName}.bed"

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fai} > ${fai.baseName}.bed
    """
    stub:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fai} > ${fai.baseName}.bed
    """
}
