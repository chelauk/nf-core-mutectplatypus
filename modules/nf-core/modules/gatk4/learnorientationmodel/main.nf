// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_LEARNORIENTATION {
    tag "$patient"
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
	tuple val(patient), path(f1r2)

    output:
    tuple val(patient), path("*tar.gz")         , emit: orientation_model
    path "versions.yml"                                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def input  = f1r2.collect{ "-I ${it} " }.join(' ')
    """
    gatk LearnReadOrientationModel \\
        ${input} \\
        -O ${prefix}.read-orientation-model.tar.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
	def input  = f1r2.collect{ "-I ${it} " }.join(' ')
    """
    echo -e  "gatk LearnReadOrientationModel \\
        ${input} \\
        -O ${prefix}.read-orientation-model.tar.gz \\
        $options.args"

    touch ${prefix}.read-orientation-model.tar.gz
	touch versions.yml
	"""
}

