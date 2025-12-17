# recombass/cli.py
import click
from . import maffilter, nr, pre, pmrcal, wtdenoise
import os

@click.command()
@click.option('-m', '--maf', type=float, help='Max major allele frequency (e.g., 0.95)')
@click.option('-n', '--nonredundant', type=float, help='Non-redundant strain clustering threshold (e.g., 0.01)')
@click.option('-cc', '--cold-criterion', type=float, default=-0.4, show_default=True, help='Cold spot threshold')
@click.option('-hc', '--hot-criterion', type=float, default=0.6, show_default=True, help='Hot spot threshold')
@click.option('-i', '--input', required=True, type=click.Path(exists=True), help='Input SNP matrix (TSV)')
@click.option('-o', '--output', required=True, help='Output prefix')
def main(maf, nonredundant, cold_criterion, hot_criterion, input, output):
    """Recombination detection pipeline from SNP matrix."""
    path = input
    intermediate_path = output

    # Step 1: Non-redundant filtering
    if nonredundant is not None:
        path = nr.process_noredundant(path, intermediate_path, percentage=nonredundant)
        intermediate_path = path  # update for next step

    # Step 2: MAF filtering
    if maf is not None:
        path = maf.filter_maf(path, maf_rate=maf, output_path=intermediate_path)

    # Step 3: SNP distance & PMR calculation
    pre.process_snp_dists(path, intermediate_path)
    pmrcal.process_snp_data(path, intermediate_path)

    # Step 4: Wavelet denoising and plotting
    wtdenoise.plot_and_save(
        wtdenoise.get_crange,
        wtdenoise.get_hrange,
        wtdenoise.wavelet_transform,
        result_path=intermediate_path,
        hc=hot_criterion,
        lc=cold_criterion
    )

if __name__ == '__main__':
    main()