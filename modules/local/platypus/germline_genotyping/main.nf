process PLATYPUS_GERMLINE_GENOTYPING {
    tag "$patient"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::platypus-variant:0.8.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platypus-variant:0.8.1.2--py27hb763d49_0' :
        'quay.io/biocontainers/platypus-variant:0.8.1.2--py27hb763d49_0' }"

    input:
    tuple val(patient), val(sample), val(status), val(id), val(gender), path(bam), path(bai)
    path genotype_ref
    path fasta
    path fasta_fai
    each chr

    output:
    tuple val(patient), val(sample), val(status), val(id), val(gender), val(chr), path("*vcf"), emit: gen_platypus
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[Platypus callvariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    platypus callVariants --refFile=${fasta} \
    --bamFiles=${bam} \
    --output=vcftemp \
    --source=${genotype_ref}/ALL.chr${chr}_GRCh38.genotypes.20170504_amended_good_mappability_bgzip.vcf.gz \
    --filterReadPairsWithSmallInserts=0 \
    --maxReads=100000000 \
    --maxVariants=100 \
    --minPosterior=0 \
    --getVariantsFromBAMs=0 \
    --nCPU=${task.cpus} \
    --regions=chr${chr} \
    --logFileName ${prefix}_${chr}.log
    grep \"^#\\|PASS\" vcftemp > ${id}.${chr}.genotype.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo    "grep \"^#\\|PASS\" > vcftemp ${id}.${chr}.genotype.vcf "
    echo    "platypus callVariants --refFile=${fasta} "
    echo    "--bamFiles=${bam} "
    echo    "--output=${id}.${chr}.genotype.vcf "
    touch   ${id}.chr${chr}.genotype.vcf
    echo    "--source=${genotype_ref}/ALL.${chr}_GRCh38.genotypes.20170504_amended_good_mappability_bgzip.vcf.gz "
    echo    "--filterReadPairsWithSmallInserts=0 "
    echo    "--maxReads=100000000 "
    echo    "--maxVariants=100 "
    echo    "--minPosterior=0 "
    echo    "--getVariantsFromBAMs=0 "
    echo    "--nCPU=${task.cpus} "
    echo    "--logFileName {output.out_log}"

    touch versions.yml
    """
}

