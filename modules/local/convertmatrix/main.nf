// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process CONVERTMATRIX {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f3c57b3530f2013369ff7ff9732cfe60ed42f63086224c61a469ba978ef3bc1f/data':
        'community.wave.seqera.io/library/r-data.table_r-matrix_r-r.utils:b228e2390ddbd670' }"

    input:
    tuple val(meta), path(bamdir), path(barcodes_tsv)

    output:
    tuple val(meta)             , path("*.tsv")                                                                                 , emit: tsv
    tuple val("${task.process}"), val('R')         , val("4.5.3"), topic: versions, emit: versions_r
    tuple val("${task.process}"), val('data.table'), val("1.17.8")                                                              , topic: versions, emit: versions_data_table
    tuple val("${task.process}"), val('Matrix')    , val("1.7_4")                                                               , topic: versions, emit: versions_matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """#!Rscript
    library(data.table)
    library(Matrix)
    matrix_dir <- "${bamdir}/Gene/raw"
    feature.names = fread(file.path(matrix_dir, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
    barcode.names = fread(file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
    barcode_to_sample = fread("${barcodes_tsv}", header = TRUE, stringsAsFactors = FALSE, data.table = F)
    sample.names = barcode_to_sample\$sample_id[match(barcode.names\$V1, barcode_to_sample\$barcode)]
    f <- gzfile(file.path(matrix_dir, "umiDedup-1MM_Directional.mtx.gz"), "r")
    mat <- as.data.frame(as.matrix(readMM(f)))
    close(f)
    colnames(mat) <- sample.names
    rownames(mat) <- feature.names\$V1
    fwrite(mat, file = "${prefix}.umi_counts.tsv", sep = "\\t", quote = F, row.names = T,
    col.names = T)
    f <- gzfile(file.path(matrix_dir, "umiDedup-NoDedup.mtx.gz"), "r")
    mat <- as.data.frame(as.matrix(readMM(f)))
    close(f)
    colnames(mat) <- sample.names
    rownames(mat) <- feature.names\$V1
    fwrite(mat, file = "${prefix}.read_counts.tsv", sep = "\\t", quote = F, row.names = T,
    col.names = T)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
