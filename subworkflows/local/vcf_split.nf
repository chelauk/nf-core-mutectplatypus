process VCF_SPLIT{
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf)

    ouput:
    tuple val(meta), path("*mono.vcf") , emit: vcf
    path "versions.yml"           , emit: versions

    script:
    """
    sed -n -e '/tumor_sample/s/##tumor_sample=//p' $vcf > temp
    while read -r sample
      do
      bcftools view $vcf -s \$sample > "\$sample"_mutect2.mono.vcf
      done<temp
    """


}