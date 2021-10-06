// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MUTECT2 {
    tag "$patient"
    label 'process_medium'
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
    tuple val(patient), val(which_tumour), val(which_control), path(bam)
    path fasta
    path fasta_fai
    path dict
    path germline_resource
    path germline_resource_idx

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def inputsList = []
    def normalsList = []
    def panelsCommand = ''
    def inputsCommand = ''

    bam.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ')
    which_norm.each() {a -> normalsList.add(" -normal " + a ) }
    normalsCommand = normalsList.join( ' ')
    panelsCommand = " --germline-resource $germline_resource --f1r2-tar-gz ${prefix}.f1r2.tar.gz"

    """
    gatk Mutect2 \\
        -R ${fasta} \\
        ${inputsCommand} \\
        ${normalsCommand} \\
        ${panelsCommand} \\
        -O ${prefix}.vcf.gz \\
        $options.args
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${patient}"
    def inputsList = []
    def normalsList = []
    def panelsCommand = ''
    def inputsCommand = ''

    bam.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ')
    which_norm.each() {a -> normalsList.add(" -normal " + a ) }
    normalsCommand = normalsList.join( ' ')
    panelsCommand = " --germline-resource $germline_resource --f1r2-tar-gz ${prefix}.f1r2.tar.gz"

    """
    echo -e "gatk gatk Mutect2 \\\n
             -R ${fasta} \\\n
        ${inputsCommand} \\\n
        ${normalsCommand} \\\n
        ${panelsCommand} \\\n
        -O ${prefix}.vcf.gz \\\n
        $options.args        \n"
    touch ${patient}.vcf.gz
    touch ${patient}.vcf.gz.tbi"
    touch ${patient}.stats"
    touch ${patient}.f1r2.tar.gz"
    touch "gatk.versions.yml"
    """
}
