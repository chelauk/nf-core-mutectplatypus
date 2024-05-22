process VCF2MAF {
    tag "$patient"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcf2maf" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'vcf2maf.1.6.21--hdfd78af_0.sif' : null }"

    input:
    tuple val(patient), val(tumour_id), val(control_id), path(vcf)
    path fasta

    output:
    tuple val(patient), path("*maf"), emit: maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf2maf.pl \\
    --input-vcf $vcf \\
    --retain-info MAPPABILITY \\
    --output-maf ${tumour_id}.maf \\
    --ref-fasta $fasta \\
    --tumor-id ${tumour_id} \\
    --normal-id ${control_id} \\
    --vcf-tumor-id ${tumour_id} \\
    --vcf-normal-id ${control_id} \\
    --inhibit-vep \\
    --ncbi-build GRCh38
    """
    
    stub:
    """
    touch $id.maf
    """
}
