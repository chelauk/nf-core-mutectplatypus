// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CONCAT_VCF {
    tag "$patient"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::htslib=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        //TODO: No singularity container at the moment, use docker container for the moment
		https://depot.galaxyproject.org/singularity/htslib%3A1.12--hd3b49d5_0
        container "quay.io/biocontainers/htslib:1.12--h9093b5e_1"
    } else {
        container "quay.io/biocontainers/htslib:1.12--hd3b49d5_0"
    }


    input:
    tuple val(patient), path(vcf)
    path fasta_fai

    output:
    tuple val(patient), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    script:
    def prefix           = options.suffix ? "${options.suffix}_${meta.id}" : "${patient}"
    """
    for i in *vcf.gz
	  do
        bgzip -d \$i
	  done
	concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}.vcf
    """
    stub:
    def prefix           = options.suffix ? "${options.suffix}_${meta.id}" : "${patient}"
    """
    echo -e "concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}.vcf"
    touch ${prefix}_concatenated.vcf.gz
    touch ${prefix}_concatenated.vcf.gz.tbi
    """
}
