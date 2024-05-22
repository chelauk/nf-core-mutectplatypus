process VCF_SPLIT{
    tag "$meta_tumour.patient"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(patient), val(meta_control), val(meta_tumour), path(vcf)

    output:
    tuple val(meta_control), val(meta_tumour), path("*mono.vcf") , emit: vcf
    path "versions.yml"           , emit: versions

    script:
    """
    bcftools view $vcf -s $meta_control.id,$meta_tumour.id > ${meta_tumour.id}_mutect2.mono.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*(bedtools) v//; s/ .*\$//')
    END_VERSIONS
    """


}
