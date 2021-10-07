// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_GATHERPILEUPSUMMARIES {
    tag "$id"
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
    tuple val(patient), val(id), val(id_intervals), val(status), path(pileup_tables)
    path dict

    output:
	tuple val(patient), val(id),  val(status), path('*.pileups_gathered.table'), emit: gathered_table
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${id}"
    def allPileups = pileup_tables.collect{ "-I ${it} " }.join(' ')

    """
    gatk GatherPileupSummaries \\
        --sequence-dict $dict \\
        ${allPileups} \\
		-O ${prefix}.pileups_gathered.table \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${id}"
    def allPileups = pileup_tables.collect{ "-I ${it} " }.join(' ')
    """
    echo -n "gatk GatherPileupSummaries \\
        --sequence-dict $dict \\
        ${allPileups} \\
        -O ${prefix}.pileups_gathered.table \\
        $options.args"
		touch ${prefix}.pileups_gathered.table
		touch versions.yml
   """
}
