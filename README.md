# recombass

[![PyPI](https://img.shields.io/pypi/v/recombpy)](https://pypi.org/project/recombpy/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)]()

A command-line tool for detecting **recombination hotspots and coldspots** from bacterial SNP matrices using wavelet-based denoising and PMR (Pairwise Mismatch Rate) analysis.

Designed for microbial population genomics, `recombass` integrates strain filtering, MAF pruning, distance computation, and signal smoothing to identify genomic regions with elevated or suppressed recombination signals.

## Requirements

- Python 3.8-3.13


##  Installation

### From bioconda (recommended)
- Use `conda` or `mamba` create a new environment with suitable python version (suggested)
```bash
mamba create -n recombass python=3.10
mamba activate recombass
```
- install `recombass` in the environment or in your terminal

```bash
mamba install -c conda-forge -c bioconda recombass
```


## Basic Usage

### Input Format

The input SNP matrix should have SNP positions as rows and bacterial strains as columns. Here is an example:

```
    strain1  strain2  strain3  strain4
100    A        T        A        G
200    C        C        T        C
300    G        A        G        A
400    T        T        C        T
500    G        G        A        G
```

In this format:
- The first column contains SNP positions
- The first row contains strain names
- Each cell contains the nucleotide at that position for that strain

### Command Line Interface

The main command structure:

```bash
recombass [OPTIONS] -m  -i <input_path> -o <output_prefix>
```

#### Available Options:

- `-m`, `--maf` FLOAT: Max major allele frequency (e.g., 0.95) (highly suggested)
- `-n`, `--nonredundant` FLOAT: Non-redundant strain clustering threshold (e.g., 0.01)
- `--cold-criterion` FLOAT: Coldspot detection threshold (default: -0.4)
- `--hot-criterion` FLOAT: Hotspot detection threshold (default: 0.6)
- `-t`, `--threads` INT: Number of threads to use (default: 20)
- `-i`, `--input` PATH: Input file path (required)
- `-o`, `--output` PATH: Output prefix (required)

#### Example Usage:

```bash
# Basic usage
recombass -i <snps.tsv> -o <result_prefix>

# With MAF filtering and non-redundant filtering
recombass -i <snps.tsv> -o <result_prefix> -m 0.95 -n 0.01

# With custom thresholds and threads
recombass -i <snps.tsv> -o <result_prefix> -m 0.9 -n 0.02 --cold-criterion -0.5 --hot-criterion 0.7 -t 10
```

## Output Files

The pipeline generates several output files:

- `<result_prefix>.maf{rate}`: MAF-filtered SNP matrix (when -m option used)
- `<result_prefix>.nr{rate}`: Non-redundant filtered SNP matrix (when -n option used)
- `<result_prefix>.fa.dist.png`: Distribution of pairwise distances
- `<result_prefix>.integrated.tsv`: PMR values and status for each SNP 
- `<result_prefix>.recx.txt`: Recombination index value
- `<result_prefix>.fa.pmr.30.2.wt.l.pdf`: Plot of hotspots from close strains comparison
- `<result_prefix>.fa.pmr.all.2.wt.l.pdf`: Plot of coldspots from all strains comparison


## License

This project is licensed under the MIT License.

## Citation

If you use this tool in your research, please cite:

> 
