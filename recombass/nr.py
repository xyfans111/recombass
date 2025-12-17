import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import polars as pl
import subprocess as sp
def process_noredundant(input_path: str, output_path: str,percentage:float=0.01,n_jobs: int = 20):
    # 生成FASTA格式数据
    df = pl.read_csv(input_path, separator="\t", row_index_name='pos')
    result = df.select(pl.all().str.join())

    fa=[f'>{result[series].name}\n{result[series].item()}\n' for series in result.columns]
    
    with open(f'{output_path}.fa', 'w') as f:
        f.write(''.join(fa[2:]))  # 跳过前两列
    
    # 执行snp-dists计算
    sp.run(f"snp-dists -j {n_jobs} {output_path}.fa > {output_path}.fa.dist", shell=True)
    distance=pd.read_table(f'{output_path}.fa.dist',index_col=0,header=0)
    clustering=AgglomerativeClustering(n_clusters=None,metric='precomputed',linkage='average',distance_threshold=percentage*np.median(distance.values))
    cluster=clustering.fit(distance)
    dft=pd.Series(cluster.labels_,index=df.columns[2:]).drop_duplicates()
    output_path+=f'.nr{int(percentage*100):02d}'
    df[['']+dft.index.to_list()].write_csv(output_path,include_header=True,separator='\t')
    return output_path
