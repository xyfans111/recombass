import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import polars as pl
import subprocess as sp

def process_noredundant(input_path: str, output_path: str, percentage:float=0.01, n_jobs: int = 20):
    # Generate FASTA format data
    df = pl.read_csv(input_path, separator="\t")
    position_col = df.columns[0]
    sequence_cols = df.columns[1:]
    result = df.select(pl.col(sequence_cols).str.join())

    fa=[f'>{sample}\n{result[sample].item()}\n' for sample in sequence_cols]
    
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
    distance=pd.read_table(f'{output_path}.fa.dist',index_col=0,header=0)
    clustering=AgglomerativeClustering(n_clusters=None,metric='precomputed',linkage='average',distance_threshold=percentage*np.median(distance.values))
    cluster=clustering.fit(distance)
    dft=pd.Series(cluster.labels_,index=sequence_cols).drop_duplicates()
    output_path+=f'.nr{int(percentage*100):02d}'
    df.select([position_col] + dft.index.to_list()).write_csv(output_path,include_header=True,separator='\t')
    return output_path
