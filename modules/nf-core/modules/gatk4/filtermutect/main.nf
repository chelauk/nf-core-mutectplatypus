// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_FILTERMUTECT {
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
    tuple val(patient), path(vcf), path(orientation_model), path(contamination_table), path(segmentation_table)
    path fasta
	path fasta_fai

    output:
    tuple val(patient), path("*.filtered.vcf.gz")      , emit: vcf
    path "versions.yml"                       , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def allsegs = segmentation_table.collect{ "--segmentation-table ${it} " }.join(' ')
    def allconts = contamination_table.collect{ "--contamination-table ${it} " }.join(' ')
    """
    gatk FilterMutectCalls \\
	-R ${fasta} \\
    -V ${vcf} \\
    ${allsegs} \\
    ${allconts} \\
    --ob-priors $orientation_model \\\
    -O ${prefix}.mutect2.filtered.vcf.gz
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def allsegs = segmentation_table.collect{ "--segmentation-table ${it} " }.join(' ')
    def allconts = contamination_table.collect{ "--contamination-table ${it} " }.join(' ')
    """
    echo -e  "gatk FilterMutectCalls \\
	-R ${fasta} \\
    -V ${vcf} \\
    ${allsegs} \\
    ${allconts} \\
    --ob-priors $orientation_model \\
    -O ${prefix}.mutect2.filtered.vcf.gz"
    touch ${prefix}.mutect2.filtered.vcf.gz
    touch versions.yml
	"""
}

