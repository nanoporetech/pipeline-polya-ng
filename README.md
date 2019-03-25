![ONT_logo](/ONT_logo.png)
-----------------------------

Pipeline for calling poly(A) tail lengths from nanopore direct RNA data using nanopolish
========================================================================================

This pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/), [minimap2](https://github.com/lh3/minimap2) and [nanopolish](https://github.com/jts/nanopolish) to call poly(A) tails from Oxford Nanopore direct RNA data.

The [pipeline-polya-diff](https://github.com/nanoporetech/pipeline-polya-diff) pipeline takes the output file `tails/filtered_tails.tsv` from multiple controls and treated samples and performs analysis of shifts in poly(A) tail lengths.

Getting Started
===============

## Input

The input files and parameters are specified in `config.yml`:

- `transcriptome` - the input transcriptome.
- `fast5_dir` - directory with pass FAST5 files.
- `fastq_dir` - directory with the fastq files.
- `summary_dir` - directory with the sequencing summary files.
- `spikein_fasta` - (optional) fasta file with spike-inf on known poly(A) tails length. The sequence names must end in *_<tail_length>* (for example "_50").
- `min_mapping_qual` - filter out reads with mapping quality less than this parameter.
- `per_transcript_plots` - plot the distribution of estimated tails lengths for all transcript (true or false).
- `threads` - number of threads to use for the analyses.


## Output

- `alignment/`:
    - `aligned_reads_sorted.bam` - sorted indexed alignment of reads to the transcriptome.

- `input/`:
    - `reads.fastq&ast` - concatenated input reads and nanopolish index files.
    - `reference.fas` - reference fasta (including spike-ins).
    - `summaries.fofn` - list of sequencing summary files.

- `reports/`:
    - `filtering_report.pdf` and `filtering_report.tsv` - nanopolish QC statistics.
    - `spikein_medians.tsv` - expected and estimated medians of spike-ins.
    - `spikein_report.pdf` - plots of distribution of tail lengths in spike-ins.
    - `tails_report.pdf` - global and per-transcript poly(A) tail length distributions.

- `tails/`:
    - `all_tails.tsv` - raw nanopolish output.
    - `filtered_tails.tsv` - nanopolish output - PASS reads only.
    - `spikein_tails.tsv` - results for reads mapping to spike-ins.

## Dependencies

- [miniconda](https://conda.io/miniconda.html) - install it according to the [instructions](https://conda.io/docs/user-guide/install/index.html).
- [snakemake](https://anaconda.org/bioconda/snakemake) install using `conda`.
- The rest of the dependencies are automatically installed using the `conda` feature of `snakemake`.

## Layout

* `README.md`
* `Snakefile`         - master snakefile
* `config.yml`        - YAML configuration file
* `snakelib/`         - snakefiles collection included by the master snakefile
* `lib/`              - python files included by analysis scripts and snakefiles
* `scripts/`          - analysis scripts
* `data/`             - input data needed by pipeline - use with caution to avoid bloated repo
* `results/`          - pipeline results to be commited - use with caution to avoid bloated repo

## Installation

Clone the repository:

```bash
git clone https://github.com/nanoporetech/pipeline-polya-ng
```

## Usage:

Edit `config.yml` to set the input datasets and parameters then issue:

```bash
snakemake --use-conda -j <num_cores> all
```

Help
====

## Licence and Copyright

(c) 2019 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips

## References and Supporting Information

### Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.

