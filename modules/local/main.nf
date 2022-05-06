process PLATYPUS_VARIANT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platypus-variant:0.8.1.2--py27hb763d49_0' :
        'quay.io/biocontainers/platypus-variant:0.8.1.2--py27hb763d49_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = input.collect{ "--input $it"}.join(" ")
    def interval_command = intervals ? "--intervals $intervals" : ""
    def pon_command = panel_of_normals ? "--panel-of-normals $panel_of_normals" : ""
    def gr_command = germline_resource ? "--germline-resource $germline_resource" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    if [[ -f ${intervalBed} ]]; then
        awk 'BEGIN{OFS=""}{print \$1,":",\$2,"-",\$3}' ${intervalBed} > ${intervalBed}.txt
    fi

    platypus callVariants \
        --refFile=${fasta} --bamFiles=${bams.join(',')} \
        --output=${intervalBed.baseName}_${idPatient}.vcf \
        --source=${mutect2Vcf.join(',')} \
        --filterReadPairsWithSmallInserts=0 \
        --maxReads=100000000 \
        --maxVariants=100 \
        --minPosterior=0 \
        --nCPU=${task.cpus} \
        --getVariantsFromBAMs=0 \
        ${intervalsOptions} \
        --logFileName ${idPatient}_${intervalBed.baseName}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
