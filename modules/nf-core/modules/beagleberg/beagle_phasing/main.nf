process BEAGLE_PHASING {
    tag "$id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beagle:4.1_21Jan17.6cc.jar--0' :
        'quay.io/biocontainers/beagle:4.1_21Jan17.6cc.jar--0' }"

    input:
    tuple val(patient), val(sample), val(status), val(id), val(gender), val(chr), path(vcf)
    path genome_ref
    path beagle_plink

    output:
    tuple val(patient), val(sample), val(status), val(id), val(gender), val(chr), path("*phased"), emit: beagle_phase
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    beagle ref=${genome_ref}/ALL.chr${chr}_GRCh38.genotypes.20170504_amended.vcf.gz \
        gt=${vcf} \
        out=${id}.phased.${chr}.vcf.gz \
        nthreads=${task.cpus} \
        map=${beagle_plink}/plink.chr${chr}.GRCh38_amended.map \
        impute=false \
        burnin=5000 \
        iterations=10000 \
        seed=1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: 4.1.0
    END_VERSIONS
    """

	stub:
	def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo -e "java -jar beagle.jar ref=${genome_ref}/ALL.chr${chr}_GRCh38.genotypes.20170504_amended.vcf.gz \\
        gt=${vcf} \\
        out=${id}.phased.chr${chr} \\
        nthreads=${task.cpus} \\
        map=${beagle_plink}/plink.chr${chr}.GRCh38_amended.map \\
        impute=false \\
        burnin=5000 \\
        iterations=10000 \\
        seed=1 "

    touch ${id}.chr${chr}.phased

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: 4.1.0
    END_VERSIONS
    """
}


