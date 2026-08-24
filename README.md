# ClinicalGenomicsGBG/brb-seq

## Introduction

**ClinicalGenomicsGBG/brb-seq** is a bioinformatics pipeline that preprocesses raw sequencing data from BRB-seq and computes count matrices.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

This pipeline processes pooled BRB-seq FASTQ files and demultiplexes them using a sample barcode TSV before alignment and counting. Prepare a samplesheet listing each UDI run and the corresponding FASTQ pair plus barcode file.

`samplesheet.csv`:

```csv
udi,fastq_1,fastq_2,barcodes
MQ-UDI-1,/path/to/fastq/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/AEG588A1_S1_L002_R2_001.fastq.gz,/path/to/barcodes/AEG588A1.barcodes.tsv
MQ-UDI-2,/path/to/fastq/AEG588A2_S2_L002_R1_001.fastq.gz,/path/to/fastq/AEG588A2_S2_L002_R2_001.fastq.gz,/path/to/barcodes/AEG588A2.barcodes.tsv
```

Each row describes one BRB-seq sample set. The first four columns are mandatory and must be provided in this order:

- `udi`: unique dual index name
- `fastq_1`: path to read 1 FASTQ (`.fq.gz` or `.fastq.gz`)
- `fastq_2`: path to read 2 FASTQ (`.fq.gz` or `.fastq.gz`)
- `barcodes`: path to a tab-separated barcode file with columns `sample_id` and `barcode`

The pipeline accepts additional columns if needed, but the first four columns must match the schema above. A full example with barcode layout is included in [`assets/samplesheet.csv`](assets/samplesheet.csv) and [`assets/barcodes.tsv`](assets/barcodes.tsv).

Run the pipeline with:

```bash
nextflow run ClinicalGenomicsGBG/brb-seq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta reference_genome.fa \
   --gtf reference_genome.gtf \
   -output-dir <OUTDIR>
```

If you have a pre-computed STAR index for your genome, supply it using `--star_index` and omit `--fasta` and `--gtf`.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

ClinicalGenomicsGBG/brb-seq was originally written by Daniel Schmitz.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use ClinicalGenomicsGBG/brb-seq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
