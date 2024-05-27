process GATK4_GATHERPILEUPSUMMARIES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"


    input:
    tuple val(meta), path(table)
    path  dict

    output:
    tuple val(meta), path("*.pileupsummaries.table"), emit: table
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = table.collect{ "--I $it" }.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GatherPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" GatherPileupSummaries \\
        $input_list \\
        --O ${prefix}.unsorted \\
        --sequence-dictionary $dict \\
        --tmp-dir . \\
        $args
    
    head -2 ${prefix}.unsorted > header
    tail -n +3 ${prefix}.unsorted | sort -k1,1 -k2,2n > temporary_file
    cat header temporary_file > ${prefix}.pileupsummaries.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def input_list = table.collect{ "--I $it" }.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GatherPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    echo -e "gatk --java-options "-Xmx${avail_mem}g" GatherPileupSummaries \\
        $input_list \\
        --O ${prefix}.pileupsummaries.table \\
        --sequence-dictionary $dict \\
        --tmp-dir . \\
        $args"
    
    touch ${prefix}.pileupsummaries.table 
    touch versions.yml
    """
}
