// bam2seqz exits with code 0 even if it fails
// so I include if statement to check the size of the output file

process SEQUENZAUTILS_BAM2SEQZ {
    tag  "${meta.id}_$chromosome"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h67e14b5_5' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39h67e14b5_5' }"

    input:
    //tuple val(patient), val(id), val(chr), path(tumourbam), path(normalbam)
    tuple val(meta), path(tumourbam), path(normalbam)
    path fasta
    path fasta_fai
    val  het
    path wigfile
    each chromosome

    output:
    tuple val(meta), path("*seqz.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-C ${chromosome}"
    def prefix = task.ext.prefix ?: "${meta.id}_${chromosome}"
    """
    chromosome=$chromosome
    sequenza-utils \\
        bam2seqz \\
        $args \\
        -n ${normalbam[0]} \\
        -t ${tumourbam[0]} \\
        --fasta $fasta \\
        --het $het \\
        -gc $wigfile \\
        -o ${prefix}.seqz.gz
    
     
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    
    file_size=\$(head -n 10 < <( zcat ${prefix}.seqz.gz ) | wc -l )

    if [ \$file_size -eq 1 ] && [ \$chromosome != "chrY" ]
      then
        exit 104
      fi
    """
    stub:
    def args = task.ext.args ?: "-C ${chromosome}"
    def prefix = task.ext.prefix ?: "${meta.id}_${chromosome}"
    """
    echo -e "sequenza-utils \\
        bam2seqz \\
        $args \\
        -n ${normalbam[0]} \\
        -t ${tumourbam[0]} \\
        --fasta $fasta \\
        --het $het \\
        -gc $wigfile \\
        -o ${prefix}.seqz.gz"
    touch ${prefix}.seqz.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
