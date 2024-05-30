process VCF_SPLIT{
    debug true
    tag "${meta_tumour.id}"
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
    def suffix = task.ext.args ?: "mutect2"
    def control = task.ext.control ?: "${meta_control.id}" ?: "${meta_control.sample}" 
    def tumour = task.ext.tumour ?: "${meta_tumour.id}" ?: "${meta_tumour.sample}" 
    """
    bcftools view $vcf -s $control,$tumour > ${tumour}_${suffix}.mono.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*(bedtools) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def suffix = task.ext.args ?: "mutect2"
    def control = task.ext.control ?: "${meta_control.id}" ?: "${meta_control.sample}" 
    def tumour = task.ext.tumour ?: "${meta_tumour.id}" ?: "${meta_tumour.sample}" 
    """
    touch ${tumour}_${suffix}.mono.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub 
    END_VERSIONS
    """


}
