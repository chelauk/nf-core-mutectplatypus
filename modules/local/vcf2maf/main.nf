process VCF2MAF {
    tag "${patient}_${id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcf2maf" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'vcf2maf.1.6.21--hdfd78af_0.sif' : null }"

    input:
    tuple val(patient), val(id), path(vcf)
    path fasta

    output:
    tuple val(patient), val(id), path("*maf"), emit: maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf2maf.pl \
    --input-vcf $vcf \
    --output-maf ${id}.maf \
    --ref-fasta $fasta \
    --inhibit-vep \
    --tumor-id $id \
    --ncbi-build GRCh38
    """
    
    stub:
    """
    touch $id.maf
    """
}
