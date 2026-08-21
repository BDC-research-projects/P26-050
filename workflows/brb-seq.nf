/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BARCODESWHITELIST      } from '../modules/local/barcodeswhitelist/main'
include { CONVERTMATRIX          } from '../modules/local/convertmatrix/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FQTK                   } from '../modules/nf-core/fqtk/main'
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../modules/nf-core/gunzip/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { STAGEFASTQDIR          } from '../modules/local/stagefastqdir/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STARSOLO               } from '../modules/nf-core/star/starsolo/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_brb-seq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BRB_SEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel: FASTA file path from --fasta
    ch_gtf         // channel: GTF file path from --gtf
    unzip_fasta    // boolean parameter: whether to unzip FASTA file (if gzipped) for STAR genome generation
    unzip_gtf      // boolean parameter: whether to unzip GTF file (if gzipped) for STAR genome generation

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Barcodes TSV (sample_id, barcode) referenced per-sample in the samplesheet.
    // Used to: (1) build a plain STARsolo whitelist, (2) demultiplex with FQTK,
    // and (3) label CONVERTMATRIX output columns with sample names.
    //
    ch_barcodes = ch_samplesheet
        .map { meta, _reads1, _reads2, barcodes_file -> [meta, barcodes_file] }

    // STARsolo's --soloCBwhitelist does not support sample names, so the
    // sample_id column has to be stripped before it can be used as a whitelist.
    BARCODESWHITELIST ( ch_barcodes )

    ch_samplesheet
        .join( BARCODESWHITELIST.out.whitelist )
        .multiMap { meta, reads1, reads2, _barcodes_file, whitelist ->
            star_fq: [meta, "CB_UMI_Simple", [reads1, reads2].transpose().flatten()]
            star_barcodes: whitelist
        }
        .set { ch_input }

    ch_samplesheet
        .map { meta, reads1, reads2, barcodes_file ->
            [meta, reads1 + reads2]
        }
        .transpose()
        .set { ch_fastqc_input }
      

    FASTQC ( ch_fastqc_input )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{_meta, file -> file})

    //
    // Demultiplex the multiplexed FASTQs with FQTK for QC/archival purposes only.
    // STARsolo further below still consumes the original multiplexed FASTQs.
    //
    ch_samplesheet
        .map { meta, reads1, reads2, _barcodes_file ->
            [meta, reads1 + reads2]
        }
        .set { ch_fqtk_reads }

    STAGEFASTQDIR ( ch_fqtk_reads )

    ch_samplesheet
        .map { meta, reads1, reads2, barcodes_file ->
            def read_structure_pairs = [reads1, reads2].transpose().collectMany { reads1_file, reads2_file ->
                [[reads1_file.name, '14B14M'], [reads2_file.name, '90T']]
            }
            [meta, barcodes_file, read_structure_pairs]
        }
        .join( STAGEFASTQDIR.out.dir )
        .map { meta, barcodes_file, read_structure_pairs, fastq_dir ->
            [meta, barcodes_file, fastq_dir, read_structure_pairs]
        }
        .set { ch_fqtk_input }

    FQTK ( ch_fqtk_input )
    ch_multiqc_files = ch_multiqc_files.mix(FQTK.out.metrics.map { _meta, file -> file })

    if (unzip_fasta) {
        GUNZIP_FASTA ( ch_fasta )
        ch_fasta = GUNZIP_FASTA.out.gunzip.collect()
    }

    if (unzip_gtf) {
        GUNZIP_GTF ( ch_gtf )
        ch_gtf = GUNZIP_GTF.out.gunzip.collect()
    }

    STAR_GENOMEGENERATE (
        ch_fasta,
        ch_gtf
    )

    STARSOLO (
        ch_input.star_fq,
        ch_input.star_barcodes,
        STAR_GENOMEGENERATE.out.index.collect(),
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.log_final.map { _meta, file -> file } )

    // All STARsolo outputs, published together under a single "starsolo" directory
    // (replicates the process-wide publishDir behavior previously configured).
    ch_starsolo_outputs = STARSOLO.out.counts
        .mix(STARSOLO.out.log_final)
        .mix(STARSOLO.out.log_out)
        .mix(STARSOLO.out.log_progress)
        .mix(STARSOLO.out.summary)

    CONVERTMATRIX (
        STARSOLO.out.counts.join( ch_barcodes )
    )

    // All FQTK outputs (demuxed FASTQs, metrics, unmatched reads), published together.
    ch_fqtk_outputs = FQTK.out.sample_fastq
        .mix(FQTK.out.metrics)
        .mix(FQTK.out.most_frequent_unmatched)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'brb-seq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    // All MultiQC outputs (report, data, plots), published together.
    ch_multiqc_outputs = MULTIQC.out.report
        .mix(MULTIQC.out.data)
        .mix(MULTIQC.out.plots)

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc        = ch_multiqc_outputs          // channel: multiqc report + data + plots, for publishing
    starsolo       = ch_starsolo_outputs         // channel: STARsolo alignment + count outputs, for publishing
    umi_counts     = CONVERTMATRIX.out.tsv       // channel: per-sample UMI/read count matrices, for publishing
    fqtk           = ch_fqtk_outputs             // channel: FQTK demultiplexed FASTQs + metrics, for publishing
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
