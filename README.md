# recombass

[![PyPI](https://img.shields.io/pypi/v/recombass)](https://pypi.org/project/recombass/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue)]()

A command-line tool for detecting **recombination hotspots and coldspots** from bacterial SNP matrices using wavelet-based denoising and PMR (Pairwise Mismatch Rate) analysis.

Designed for microbial population genomics, `recombass` integrates strain filtering, MAF pruning, distance computation, and signal smoothing to identify genomic regions with elevated or suppressed recombination signals.

## Requirements

- Python 3.9+
- `snp-dists` available on `PATH`


## Installation

### From PyPI

Install the external runtime dependency first, then install `recombass` from PyPI:

```bash
mamba install -c conda-forge -c bioconda snp-dists
pip install recombass
```


## Basic Usage

### Input Format

The input SNP matrix must be a tab-delimited table with:

- SNP positions in the first column
- strain names in the header row
- nucleotide calls in each cell

Example:

```tsv
	strain1	strain2	strain3	strain4
100	A	T	A	G
200	C	C	T	C
300	G	A	G	A
400	T	T	C	T
500	G	G	A	G
```

### Command Line Interface

The main command structure:

```bash
recombass [OPTIONS] -m  -i <input_path> -o <output_prefix>
```

#### Available Options

- `-m`, `--maf` FLOAT: Max major allele frequency (e.g., 0.95) (highly suggested)
- `-n`, `--nonredundant` FLOAT: Non-redundant strain clustering threshold (e.g., 0.01)
- `--cold-criterion` FLOAT: Coldspot detection threshold (default: -0.4)
- `--hot-criterion` FLOAT: Hotspot detection threshold (default: 0.6)
- `-t`, `--threads` INT: Number of threads to use (default: 20)
- `-i`, `--input` PATH: Input file path (required)
- `-o`, `--output` PATH: Output prefix (required)

#### Example Usage

```bash
# Basic usage
recombass -i <snps.tsv> -o <result_prefix>

# With MAF filtering and non-redundant filtering
recombass -i <snps.tsv> -o <result_prefix> -m 0.95 -n 0.01

# With custom thresholds and threads
recombass -i <snps.tsv> -o <result_prefix> -m 0.9 -n 0.02 --cold-criterion -0.5 --hot-criterion 0.7 -t 10
```

## Output Files

The pipeline can generate the following files:

- `<result_prefix>.nr{rate}`: Non-redundant filtered SNP matrix when `-n` is used
- `<result_prefix>.maf{rate}`: MAF-filtered SNP matrix when `-m` is used
- `<result_prefix>.fa`: FASTA converted from the current SNP matrix
- `<result_prefix>.fa.dist`: Pairwise SNP distance matrix from `snp-dists`
- `<result_prefix>.fa.dist.tr`: Triangular distance table used downstream
- `<result_prefix>.fa.dist.png`: Distance distribution plot
- `<result_prefix>.recx.txt`: Recombination index value
- `<result_prefix>.fa.pmr.30.2.wt.l.pdf`: Hotspot plot from the denoised close-strain PMR signal
- `<result_prefix>.fa.pmr.all.2.wt.l.pdf`: Coldspot plot from the denoised all-strain PMR signal
- `<result_prefix>.integrated.tsv`: Integrated per-site output with raw PMR, denoised PMR, and hot/cold labels

When filtering is enabled, the effective output prefix changes as intermediate files are generated. For example:

```bash
recombass -i sample.tsv -o sample_out -m 0.95 -n 0.01
```

will produce final files such as:

- `sample_out.nr01`
- `sample_out.nr01.maf95`
- `sample_out.nr01.maf95.integrated.tsv`
- `sample_out.nr01.maf95.fa.pmr.30.2.wt.l.pdf`
- `sample_out.nr01.maf95.fa.pmr.all.2.wt.l.pdf`

## License

This project is licensed under the MIT License.
