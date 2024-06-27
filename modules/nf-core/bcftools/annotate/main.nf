process BCFTOOLS_MAPPABILITY {
    tag "$meta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    path(mappability) 
    path(mappability_tbi)

    output:
    tuple val(meta), path("*mappability.vcf")     , emit: vcf

    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output = "${meta}.mutect.mappability.vcf"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    echo '##INFO=<ID=MAPPABILITY,Number=1,Type=Float,Description="k100 mappability">' > mappability.h
    bcftools \\
    annotate -a $mappability \\
    -h mappability.h \\
    -c CHROM,FROM,TO,MAPPABILITY \\
    $vcf > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def output = "${meta}.mutect.mappability.vcf"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    echo '##INFO=<ID=MAPPABILITY,Number=1,Type=Float,Description="k100 mappability">' > mappability.h
    touch $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub 
    END_VERSIONS
    """
}
