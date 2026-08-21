process BARCODESWHITELIST {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(barcodes_tsv)

    output:
    tuple val(meta), path("*.whitelist.txt"), emit: whitelist

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // STARsolo's --soloCBwhitelist only accepts a plain, header-less list of
    // barcodes, so the sample_id column has to be dropped from the TSV.
    """
    awk -F '\\t' 'NR==1 { for (i=1; i<=NF; i++) if (\$i == "barcode") col=i; next } { print \$col }' ${barcodes_tsv} > ${prefix}.whitelist.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.whitelist.txt
    """
}
