process CONCAT_VCF {
    tag "$patient"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::htslib=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0"
    } else {
        container "quay.io/biocontainers/htslib:1.12--hd3b49d5_0"
    }


    input:
    tuple val(patient), val(intervals), path(vcf)
    path fasta_fai
	path target_bed

    output:
    tuple val(patient), path("*concatenated.vcf.gz"), path("*concatenated.vcf.gz.tbi"), emit: vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}"
    options = params.intervals ? "-t ${target_bed}" : ""
	"""
	concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}_${args}_concatenated.vcf ${options}
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}"
    options = params.intervals ? "-t ${target_bed}" : ""
    """
    echo -e "concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}_${args}_concatenated.vcf ${options}"
    touch ${prefix}_${args}_concatenated.vcf.gz
    touch ${prefix}_${args}_concatenated.vcf.gz.tbi
    """
}
