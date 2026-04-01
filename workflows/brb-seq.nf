/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CONVERTMATRIX          } from '../modules/local/convertmatrix/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../modules/nf-core/gunzip/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
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

    ch_samplesheet
        .multiMap { meta, reads1, reads2, barcodes_file ->
            star_fq: [meta, "CB_UMI_Simple", [reads1, reads2].transpose().flatten()]
            star_barcodes: barcodes_file
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
        ch_input.star_fq.dump(),
        ch_input.star_barcodes,
        STAR_GENOMEGENERATE.out.index.collect(),
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.log_final.map { _meta, file -> file } )

    CONVERTMATRIX (
        STARSOLO.out.counts
    )

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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
