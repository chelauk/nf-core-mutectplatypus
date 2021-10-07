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
        container "quay.io/biocontainers/htslib:1.12--h9093b5e_1"
    } else {
        container "quay.io/biocontainers/htslib:1.12--hd3b49d5_0"
    }


    input:
    tuple val(patient), path(vcf)
    path fasta_fai
    path bed

    output:
    tuple val(patient), path("*_*.vcf.gz"), path("*_*.vcf.gz.tbi"), emit: vcf

    script:
    def prefix           = options.suffix ? "${options.suffix}_${meta.id}" : "${patient}"
    def target_options   = params.target_bed ? "-t ${bed}" : ""
    """
    concatenateVCFs.sh -i ${fai} -c ${task.cpus} -o ${prefix}.vcf ${target_options}
    """
    stub:
    def prefix           = options.suffix ? "${options.suffix}_${meta.id}" : "${patient}"
    def target_options   = params.target_bed ? "-t ${bed}" : ""
    """
    echo -e "concatenateVCFs.sh -i ${fai} -c ${task.cpus} -o ${prefix}.vcf ${target_options}"
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}