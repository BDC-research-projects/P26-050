process FASTREADCOUNTER {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openjdk:17.0.1' :
        'openjdk:17-jdk-slim' }"

    input:
    tuple val(meta), path(bam)
    path(genomic_ranges)

    output:
    tuple val(meta), path("*.tsv"), emit: counts
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format_flag = genomic_ranges.name.endsWith('.gtf') ? '--gtf' :
                      genomic_ranges.name.endsWith('.bed') ? '--bed' :
                      genomic_ranges.name.endsWith('.vcf') ? '--vcf' : '--gtf'
    """
    java -jar \$FASTREADCOUNTER_JAR \\
        --bam ${bam} \\
        ${format_flag} ${genomic_ranges} \\
        -o . \\
        ${args}

    # Rename output TSV to use the prefix
    for f in *.tsv; do
        mv "\$f" "${prefix}.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastreadcounter: 1.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastreadcounter: 1.1.0
    END_VERSIONS
    """
}
