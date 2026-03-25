# Bioconda submission notes

This directory contains a Bioconda-ready recipe template for `recombass`.

Before opening a pull request against `bioconda-recipes`, complete these steps:

1. Build and upload the source distribution to PyPI.
2. Compute the source tarball checksum:

```bash
curl -L https://pypi.io/packages/source/r/recombass/recombass-0.1.6.tar.gz | shasum -a 256
```

3. Replace `REPLACE_WITH_PYPI_SDIST_SHA256` in `meta.yaml`.
4. Replace `REPLACE_WITH_YOUR_GITHUB_ID` in `meta.yaml`.
5. Copy this recipe into the Bioconda repository path `recipes/recombass/`.
6. Validate locally:

```bash
conda build recipes/recombass
```

If `conda build` fails on the integration test, verify that the uploaded sdist includes `tests/data/toy_snps.tsv`.
