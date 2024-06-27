process VCF2MAF {
    tag "${meta.tumour}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcf2maf" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'vcf2maf.1.6.21--hdfd78af_0.sif' : null }"

    input:
    tuple val(meta), path(vcf)
    path fasta

    output:
    tuple val(meta), path("*maf"), emit: maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf2maf.pl \\
    --input-vcf $vcf \\
    --retain-info MAPPABILITY \\
    --output-maf ${meta.tumour}.maf \\
    --ref-fasta $fasta \\
    --tumor-id ${meta.tumour} \\
    --normal-id ${meta.control} \\
    --vcf-tumor-id ${meta.tumour} \\
    --vcf-normal-id ${meta.control} \\
    --inhibit-vep \\
    --ncbi-build GRCh38
    """
    
    stub:
    """
    touch ${meta.tumour}.maf
    """
}
