process PLATYPUS_CALLVARIANTS {
    tag "${meta.patient}_${meta.interval}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::platypus-variant:0.8.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platypus-variant:0.8.1.2--py27hb763d49_0' :
        'quay.io/biocontainers/platypus-variant:0.8.1.2--py27hb763d49_0' }"

    input:
    tuple val(meta), path(vcf), val(which_tumour), val(which_norm), path(bam), path(bai), path(intervals)
    path fasta
    path fasta_fai

    output:
    tuple val("${meta.patient}"), path("*platypus.vcf")     , emit: vcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.interval}"
    def interval_command = intervals ? "--regions=${intervals}.txt" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[Platypus callvariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    if [[ -f ${intervals} ]]; then
        awk 'BEGIN{OFS=""}{print \$1,":",\$2,"-",\$3}' ${intervals} > ${intervals}.txt
    fi

    if [[ ! -f ${vcf}.gz ]]; then
        bgzip ${vcf}
        tabix -p vcf ${vcf}.gz
    fi 

    platypus callVariants \
        --refFile=${fasta} --bamFiles=${bam.join(',')} \
        --output=${prefix}.platypus.vcf \
        --source=${vcf}.gz \
        --filterReadPairsWithSmallInserts=0 \
        --maxReads=100000000 \
        --maxVariants=100 \
        --minPosterior=0 \
        --nCPU=${task.cpus} \
        --getVariantsFromBAMs=0 \
        ${interval_command} \
        --logFileName ${prefix}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.interval}"
    def interval_command = intervals ? "--regions=${intervals}.txt" : ""
    """
    if [[ -f ${intervals} ]]; then
        awk 'BEGIN{OFS=""}{print \$1,":",\$2,"-",\$3}' ${intervals} > ${intervals}.txt
    fi
    echo "
    platypus callVariants \
        --refFile=${fasta} --bamFiles=${bam.join(',')} \
        --output=${prefix}.platypus.vcf \
        --source=${vcf.join(',')} \
        --filterReadPairsWithSmallInserts=0 \
        --maxReads=100000000 \
        --maxVariants=100 \
        --minPosterior=0 \
        --nCPU=${task.cpus} \
        --getVariantsFromBAMs=0 \
        ${interval_command} \
        --logFileName ${prefix}.log"

    touch ${prefix}.platypus.vcf
    touch versions.yml
    """
}
