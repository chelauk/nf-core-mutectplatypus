process GATK4_CALCULATECONTAMINATION {
    tag "${patient}_${sample}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(patient), val(sample), path(table)

    output:
    tuple val(patient), val(sample), path('*.contamination.table'), emit: contamination
    tuple val(patient), val(sample), path('*.segmentation.table') , emit: segmentation, optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}_${sample}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK CalculateContamination] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CalculateContamination \\
        --input $table \\
        -segments ${prefix}.segmentation.table \\
        --output ${prefix}.contamination.table \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}_${sample}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK CalculateContamination] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    echo -e "gatk --java-options "-Xmx${avail_mem}g" CalculateContamination \\
        --input $table \\
        -segments ${prefix}.segmentation.table \\
        --output ${prefix}.contamination.table \\
        --tmp-dir . \\
        $args"

    touch ${prefix}.contamination.table
    touch ${prefix}.segmentation.table
    touch versions.yml 
    """
}
