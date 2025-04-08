#!/usr/bin/env python
import numpy as np
import pyfastx
import pandas as pd
import sys
from collections import Counter
import polars as pl
def filter_maf(input_path: str, maf_rate: float, output_path: str,) -> str:
    """基于次要等位基因频率(MAF)过滤序列数据
    
    Args:
        fasta_path: 输入FASTA文件路径
        maf_rate: MAF过滤阈值（0-1之间）
        output_path: 输出文件路径

        
    Returns:
        过滤后的DataFrame
    """
    df = pl.read_csv(input_path, separator="\t", row_index_name='pos')
    result = df.select(pl.all().str.join())
    fa=[f'>{result[series].name}\n{result[series].item()}\n' for series in result.columns]
    
    with open(f'{output_path}.fa', 'w') as f:
        f.write(''.join(fa[2:]))  # 跳过前两列å
    # 读取并处理FASTA数据
    
    fa = pyfastx.Fasta(f'{output_path}.fa', build_index=False)
    n = 0
    for name, seq in fa:
        n += 1
    
    index = []
    array = np.ndarray(shape=(n, len(seq)), dtype='U1')
    x = 0
    for name, seq in fa:
        index.append(name)
        array[x, :] = list(seq)
        x += 1
    
    maf_value = n * maf_rate
    array1 = np.where(np.array([Counter(array[:, i]).most_common(1)[0][1] for i in range(array.shape[1])]) <= maf_value)
    
    # 读取并处理原始TSV数据
    df = pd.read_csv(input_path, sep="\t", header=0, index_col=0)
    filtered_df = df.iloc[array1[0].tolist(),:]
    output_path = f'{output_path}.maf{int(maf_rate*100):02d}'
    # 保存结果
    filtered_df.to_csv(path_or_buf=output_path, sep="\t")
    
    return output_path

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit('Usage: python maf.py aln rate(like 0.95) output')
    input_path = sys.argv[1]
    maf_rate = float(sys.argv[2])
    output_path = sys.argv[3]
    filter_maf(input_path, maf_rate, output_path)