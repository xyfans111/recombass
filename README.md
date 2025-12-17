# recombass

[![PyPI](https://img.shields.io/pypi/v/recombpy)](https://pypi.org/project/recombpy/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)]()

A command-line tool for detecting **recombination hotspots and coldspots** from bacterial SNP matrices using wavelet-based denoising and PMR (Pairwise Mismatch Rate) analysis.

Designed for microbial population genomics, `recombass` integrates strain filtering, MAF pruning, distance computation, and signal smoothing to identify genomic regions with elevated or suppressed recombination signals.

---

##  Features

- **Non-redundant strain filtering**: Cluster highly similar isolates to reduce bias.
- **MAF (Minor Allele Frequency) filtering**: Remove rare variants that may be sequencing errors.
- **SNP distance matrix computation**: Via external `snp-dists` (from Torsten Seemann).
- **PMR profile calculation**: Sliding-window pairwise mismatch rates across the genome.
- **Wavelet denoising**: Smooth noisy PMR signals using discrete wavelet transform (DWT).
- **Hotspot/coldspot calling**: Threshold-based detection with publication-ready plots.

---

##  Installation

### From PyPI (recommended)

```bash
pip install recombpy
