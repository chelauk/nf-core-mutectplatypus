process GATK4_MUTECT2 {
    tag "$patient"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0"
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
    val af_resource
    path germline_resource
    path germline_resource_idx

    output:
    tuple val(patient), val(interval_patient), path("*.vcf")         , emit: vcf
//    tuple val(patient), path("*.tbi")         , emit: tbi
    tuple val(patient), val(interval_patient), path("*.stats")       , emit: stats
    tuple val(patient), val(interval_patient), path("*.f1r2.tar.gz") , optional:true, emit: f1r2
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${interval_patient}"
    def inputs_list = []
    def normals_list = []
    def panel_of_normals_command = ''
    def inputs_command = ''

    bam.each() {a -> inputs_list.add(" -I " + a ) }
    inputs_command = inputs_list.join( ' ')
    which_norm.each() {a -> normals_list.add(" -normal " + a ) }
    normals_command = normals_list.join( ' ')
    if (pon) {
        panel_of_normals_command = "--panel-of-normals $pon"
        } else {
            panel_of_normals_command = ""
    }

    """
    if [ ! -d tmpdir ]
    then
      mkdir -p tmpdir
    fi
    gatk Mutect2 \\
        -R ${fasta} \\
        ${inputs_command} \\
        ${normals_command} \\
        --tmp-dir ./tmpdir \\
        --germline-resource $germline_resource \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        $args $task.cpus \\
        ${panel_of_normals_command} \\
        -L $intervals \\
        -O ${prefix}.vcf \\
        --af-of-alleles-not-in-resource $af_resource \\
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 
    
    echo "GATK4.2.0" > versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${interval_patient}"
    def inputs_list = []
    def normals_list = []
    def panel_of_normals_command = ''
    def inputs_command = ''

    bam.each() {a -> inputs_list.add(" -I " + a ) }
    inputs_command = inputs_list.join( ' ')
    which_norm.each() {a -> normals_list.add(" -normal " + a ) }
    normals_command = normals_list.join( ' ')
    if (pon) {
        panel_of_normals_command = "--panel-of-normals $pon"
        } else {
            panel_of_normals_command = ""
        }
    """
    echo -e "gatk gatk Mutect2 \\
            -R ${fasta} \\
            -L ${intervals} \\
            ${inputs_command} \\
            ${normals_command} \\
            --germline-resource $germline_resource \\
            --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
            ${panel_of_normals_command} \\
            -O ${prefix}.vcf.gz \\
            $args        \n"
    touch ${prefix}.vcf
    touch ${prefix}.vcf.tbi
    touch ${prefix}.stats
    touch ${prefix}.f1r2.tar.gz
    touch versions.yml
    """
}
