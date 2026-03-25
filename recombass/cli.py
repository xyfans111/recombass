import click

from . import __version__

@click.command()
@click.version_option(__version__, prog_name="recombass")
@click.option('-m', '--maf', type=float, help='Max major allele frequency (e.g., 0.95)')
@click.option('-n', '--nonredundant', type=float, help='Non-redundant strain clustering threshold (e.g., 0.01)')
@click.option('-cc', '--cold-criterion', type=float, default=-0.4, show_default=True)
@click.option('-hc', '--hot-criterion', type=float, default=0.6, show_default=True)
@click.option("-t", "--threads", type=int, default=20, show_default=True)
@click.option('-i', '--input', required=True, type=click.Path(exists=True))
@click.option('-o', '--output', required=True, help='output prefix')
def main(maf, nonredundant, cold_criterion, hot_criterion, input, output, threads):
    """Recombination detection from SNP matrix."""
    path = input
    intermediate_path = output

    # Step 1: Non-redundant filtering
    if nonredundant is not None:
        from . import nr  
        path = nr.process_noredundant(path, intermediate_path, percentage=nonredundant, n_jobs=threads)
        intermediate_path = path

    # Step 2: MAF filtering
    if maf is not None:
        from .maffilter import filter_maf 
        path = filter_maf(path, maf_rate=maf, output_path=intermediate_path)
        intermediate_path = path

    # Step 3: SNP distance & PMR calculation
    from . import pre, pmrcal
    pre.process_snp_dists(path, intermediate_path, n_jobs=threads)
    # Get the PMR data
    pmr_data = pmrcal.process_snp_data(path, intermediate_path)

    # Step 4: Wavelet denoising and plotting (now includes output integration)
    from . import wtdenoise
    wtdenoise.plot_and_save(
        wtdenoise.get_crange,
        wtdenoise.get_hrange,
        wtdenoise.wavelet_transform,
        result_path=intermediate_path,
        pmr_data=pmr_data,  # Pass the pmr_data for integration
        hc=hot_criterion,
        lc=cold_criterion
    )


if __name__ == '__main__':
    main()
