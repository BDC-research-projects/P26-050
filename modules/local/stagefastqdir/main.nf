process STAGEFASTQDIR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("fastq_dir"), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def read_names = (reads instanceof List ? reads : [reads]).collect {file -> file.name }.join(' ')
    """
    mkdir fastq_dir
    for fastq in ${read_names}; do
        ln -s "\$(readlink -f "\${fastq}")" "fastq_dir/\${fastq}"
    done
    """

    stub:
    """
    mkdir fastq_dir
    """
}
