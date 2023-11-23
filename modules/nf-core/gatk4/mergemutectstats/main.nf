process GATK4_MERGEMUTECTSTATS {
    tag "$patient"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(patient), path(stats)

    output:
    tuple val(patient), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${patient}"
    def input_list = stats.collect{ "--stats ${it}"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MergeMutectStats] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MergeMutectStats \\
        $input_list \\
        --output ${prefix}_mutect2_concatenated.vcf.gz.stats \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${patient}"
    def input_list = stats.collect{ "--stats ${it}"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MergeMutectStats] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    echo -e "gatk --java-options "-Xmx${avail_mem}g" MergeMutectStats \\
        $input_list \\
        --output ${prefix}.vcf.gz.stats \\
        --tmp-dir . \\
        $args"
    
    touch ${prefix}.vcf.gz.stats
    touch versions.yml
    """
}
