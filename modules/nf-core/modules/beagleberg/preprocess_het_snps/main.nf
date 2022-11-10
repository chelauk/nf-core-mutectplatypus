process PREPROCESS_HETSNPS {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/r-sequenza:2.1.2--r351h6115d3f_2' :
    'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(patient), val(id), val(gender), path(het_snps)
    val(min_reads)

    output:
    tuple val(patient), val(id), val(gender), path("*BAF.rds"), path("*LRR.rds"), emit: baf_lrr
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    preprocessing_data_seqz.R ${patient} ${id} ${gender} ${min_reads}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preprocess_het_snps: 1.0.0
    END_VERSIONS
    """

	stub:
	def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo preprocessingData_Sequenza.R ${patient} ${id} ${gender} ${min_reads}
    touch ${patient}_${id}.BAF.rds
    touch ${patient}_${id}.LRR.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preprocess_het_snps: 3.0.0
    END_VERSIONS
    """
}
