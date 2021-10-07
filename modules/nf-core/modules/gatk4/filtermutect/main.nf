// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_FILTERMUTECT {
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
    tuple val(patient), path(vcf) 
	tuple val(patient_con), val(id_con), val(status_con), path(contamination_table)    
    tuple val(patient_seg), val(id_seg), val(status_seg), path(segmentation_table) 

    output:
    tuple val(patient), path("*.vcf.gz")      , emit: vcf
    path "versions.yml"                       , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def allsegs = segmentation_table.collect{ "--segmentation-table ${it} " }.join(' ')
    def allconts = contamination_table.collect{ "--contamination-table ${it} " }.join(' ')
    """
    gatk FilterMutectCalls \
    -V ${vcf} \
    ${allsegs} \
    ${allconts} \
    --ob-priors priors.tar.gz \
    -O ${prefix}.mutect2.filtered.vcf
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${id}"
    """
    echo -e  "gatk 
	touch versions.yml
	"""
}

