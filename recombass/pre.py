#!/usr/bin/env python
import polars as pl
import subprocess as sp

from .inputprep import normalize_snp_matrix


def process_snp_dists(input_path: str, output_path: str, n_jobs: int = 20):
    """Process SNP distances and generate triangular matrix result.
    
    Args:
        input_path: Input file path
        output_path: Output file path
        n_jobs: Number of parallel jobs (default 20)
    """
    # Generate FASTA format data
    df, _, sequence_cols = normalize_snp_matrix(input_path)
    result = df.select(pl.col(sequence_cols).str.join())
    fa = [f'>{sample}\n{result[sample].item()}\n' for sample in sequence_cols]
    
    with open(f'{output_path}.fa', 'w') as f:
        f.write(''.join(fa))
    
    # Execute snp-dists calculation
    try:
        with open(f"{output_path}.fa.dist", "w") as stdout:
            sp.run(
                ["snp-dists", "-j", str(n_jobs), f"{output_path}.fa"],
                stdout=stdout,
                stderr=sp.DEVNULL,
                check=True,
            )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "The 'snp-dists' executable was not found. Install it from Bioconda with "
            "'mamba install -c conda-forge -c bioconda snp-dists'."
        ) from exc
    except sp.CalledProcessError as exc:
        raise RuntimeError(
            "snp-dists failed while processing "
            f"'{input_path}' (FASTA: '{output_path}.fa', exit code {exc.returncode}). "
            "Check that the input is a valid SNP matrix with at least two sequence columns "
            "and that all sequence values are single-base calls or gaps."
        ) from exc
    
    # Process output result
    df = pl.read_csv(f"{output_path}.fa.dist", separator="\t", has_header=True)
    row_names = df.columns[1:]
    output_path = f"{output_path}.fa.dist.tr"
    
    # Build triangular matrix data
    data = []
    for i in range(len(row_names)):
        for j in range(i + 1, len(row_names)):
            data.append({
                "strain1": row_names[i],
                "strain2": row_names[j],
                "distance": df[row_names[i]][j]
            })
    
    # Write output file
    pl.DataFrame(data).write_csv(output_path, separator="\t")
    return output_path
