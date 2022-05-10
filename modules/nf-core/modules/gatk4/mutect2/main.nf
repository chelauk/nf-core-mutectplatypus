process GATK4_MUTECT2 {
    tag "$patient"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0"
    }

    input:
    tuple val(patient), val(interval_patient), val(which_tumour), val(which_norm), path(bam), path(bai), path(intervals)
    path fasta
    path fasta_fai
    path dict
    path pon
    path pon_idx
    path germline_resource
    path germline_resource_idx

    output:
    tuple val(patient), path("*.vcf")         , emit: vcf
//    tuple val(patient), path("*.tbi")         , emit: tbi
    tuple val(patient), path("*.stats")       , emit: stats
    tuple val(patient), path("*.f1r2.tar.gz") , optional:true, emit: f1r2
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${interval_patient}"
    def inputsList = []
    def normalsList = []
    def panelsCommand = ''
    def inputsCommand = ''

    bam.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ')
    which_norm.each() {a -> normalsList.add(" -normal " + a ) }
    normalsCommand = normalsList.join( ' ')
    if (pon) {
        pon_command = "--panel-of-normals $pon"
        }

    """
    gatk Mutect2 \\
        -R ${fasta} \\
        ${inputsCommand} \\
        ${normalsCommand} \\
        --germline-resource $germline_resource \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        $args \\
        ${pon_command} \\
        -L $intervals \\
        -O ${prefix}.vcf \\
        --af-of-alleles-not-in-resource 0.0000025 \\
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 
    
    echo "GATK4.2.0" > versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${interval_patient}"
    def inputsList = []
    def normalsList = []
    def panelsCommand = ''
    def inputsCommand = ''

    bam.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ')
    which_norm.each() {a -> normalsList.add(" -normal " + a ) }
    normalsCommand = normalsList.join( ' ')
    if (pon) {
        pon_command = "--panel-of-normals $pon"
        }
    """
    echo -e "gatk gatk Mutect2 \\
            -R ${fasta} \\
            -L ${intervals} \\
            ${inputsCommand} \\
            ${normalsCommand} \\
            --germline-resource $germline_resource \\
            --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
            ${pon_command} \\
            -O ${prefix}.vcf.gz \\
            $args        \n"
    touch ${prefix}.vcf
    touch ${prefix}.vcf.tbi
    touch ${prefix}.stats
    touch ${prefix}.f1r2.tar.gz
    touch versions.yml
    """
}