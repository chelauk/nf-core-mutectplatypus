process GATK4_GETPILEUPSUMMARIES {
    tag "$id_intervals"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), val(id_intervals), path(bams), path(intervals)
    path  fasta
    path  fai
    path  dict
    path  variants
    path  variants_tbi

    output:
    tuple val(meta), val(id_intervals), path('*.pileups.table'), emit: table
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id_intervals}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def reference_command = fasta ? "--reference $fasta" : ''

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GetPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" GetPileupSummaries \\
        --input ${bams[0]} \\
        --variant $variants \\
        --output ${prefix}.pileups.table \\
        $reference_command \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id_intervals}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def reference_command = fasta ? "--reference $fasta" : ''

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GetPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    echo -e "gatk --java-options "-Xmx${avail_mem}g" GetPileupSummaries \\
        --input ${bams[0]} \\
        --variant $variants \\
        --output ${prefix}.pileups.table \\
        $reference_command \\
        $interval_command \\
        --tmp-dir . \\
        $args"
    
    touch ${prefix}.pileups.table
    touch versions.yml
    """
}
