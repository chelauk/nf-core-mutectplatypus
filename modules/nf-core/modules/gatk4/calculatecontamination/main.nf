// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_CALCULATECONTAMINATION {
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
	tuple val(patient), val(id), val(status), path(table)

    output:
    tuple val(patient), val(id), val(status), path('*.contamination.table')               , emit: contamination
    tuple val(patient), val(id), val(status), path('*.segmentation.table') , optional:true, emit: segmentation
    path "versions.yml"                                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${id}"
    """
	awk -v FS='\t' -v OFS='\t' 'NR<3{print \$0;next}{print \$0| "sort -k1,1 -k2,2n"}' $table > temp && mv temp  $table 
    gatk CalculateContamination \\
        -I ${table} \\
        -segments ${prefix}.segmentation.table \\
        -O ${prefix}.contamination.table \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${id}"
    """
    echo -e  "gatk CalculateContamination \\
        -I ${table} \\
        -segments ${prefix}.segmentation.table \\
        -O ${prefix}.contamination.table \\
        $options.args"

    touch ${prefix}.segmentation.table
	touch ${prefix}.contamination.table
	touch versions.yml
	"""
}

