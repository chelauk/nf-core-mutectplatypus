process VCF2MAF {
    tag "${patient}_${id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcf2maf" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'vcf2maf.sif' : null }"

    input:
    tuple val(patient), val(id), path(vcf)
    fasta

    output:
    tuple val(patient), val(id), path("*maf"), emit: maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    perl /usr/bin/vcf2maf/vcf2maf.pl --input-vcf $vcf --output-maf ${id}.maf  --tumor-id ${id}  --ref-fasta 
    """
    
    stub:
    """
    touch $id.maf
    """
}
