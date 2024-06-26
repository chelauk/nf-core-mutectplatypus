process ENSEMBLVEP {
    tag "${patient}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ensembl-vep:108.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:108.1--pl5321h4a94de4_0' :
        'quay.io/biocontainers/ensembl-vep:108.1--pl5321h4a94de4_0' }"

    input:
    tuple val(patient), path(vcf)
    val   genome
    val   species
    val   cache_version
    path  cache
    path  extra_files

    output:
    tuple val(patient), path("*.ann.vcf"), emit: vcf
    path "*.summary.html"                , emit: report
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    """
    if [ ! -d $prefix ]; then
        mkdir $prefix
    fi

    vep \\
        -i $vcf \\
        --fasta ${dir_cache}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \\
        --everything \\
        --shift_hgvs 1 \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --xref_refseq \\
        --failed 1 \\
        --flag_pick_allele \\
        --pick_order canonical,tsl,biotype,rank,ccds,length \\
        -o ${prefix}_${args}.ann.vcf \\
        --force_overwrite \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --vcf \\
        --offline \\
        --verbose \\
        --format vcf \\
        --stats_file ${prefix}_${args}.summary.html

    rm -rf $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${patient}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    """
        mkdir $prefix
        echo -e "
        vep \\
        -i $vcf \\
        -o ${prefix}_${args}.ann.vcf \\
        $args \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --vcf \\
        --stats_file ${prefix}.summary.html"

        rm -rf ${prefix}
        touch ${prefix}_${args}.ann.vcf
        touch ${prefix}_${args}.summary.html
        touch versions.yml
        """
}
