// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_GETPILEUPSUMMARIES {
    tag "$patient_interval"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(patient), val(id), val(id_intervals), val(status), path(bam), path(bai), path(intervals)
    path germline_resource
    path germline_resource_idx

    output:
    tuple val(meta), path('*.pileups.table'), emit: table
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${id_intervals}"
    """
    gatk GetPileupSummaries \\
        -I $bam \\
        -V $germline_resource \\
        -L $intervals \\
        -O ${prefix}.pileups.table \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${id_intervals}"
    """
    echo -e "gatk GetPileupSummaries \\\n
        -I $bam \\\n
        -V $germline_resource \\\n
        -L $intervals \\\n
        -O ${prefix}.pileups.table \\\n
        $options.args"

    touch versions.yml
    """
}
