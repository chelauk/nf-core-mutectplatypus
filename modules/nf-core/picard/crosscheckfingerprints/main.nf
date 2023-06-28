process PICARD_CROSSCHECKFINGERPRINTS {
    //tag "${task.name}_${sample_tumour}_${patient}"
    label 'process_medium'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(patient), path(control_bam), path(control_bai), val(sample_tumour), path(tumour_bam), path(tumour_bai)
    path haplotype_map

    output:
    tuple val(patient), val(sample_tumour), path("*.crosscheck_metrics.txt"), emit: crosscheck_metrics
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CrosscheckFingerprints] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CrosscheckFingerprints \\
        $args \\
        --NUM_THREADS ${task.cpus} \\
        --INPUT ${control_bam} \\
        --SECOND_INPUT ${tumour_bam} \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${sample_tumour}_vs_control.crosscheck_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard CrosscheckFingerprints --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d: )
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CrosscheckFingerprints] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    echo -e "picard \\
        -Xmx${avail_mem}M \\
        CrosscheckFingerprints \\
        $args \\
        --NUM_THREADS ${task.cpus} \\
        --INPUT ${control_bam} \\
        --SECOND_INPUT ${tumour_bam} \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${sample_tumour}_vs_control.crosscheck_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard CrosscheckFingerprints --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d: )
    END_VERSIONS"
    touch ${sample_tumour}_vs_control.crosscheck_metrics.txt
    touch versions.yml
    """
}
