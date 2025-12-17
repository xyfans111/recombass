#!/usr/bin/env python
import polars as pl
import subprocess as sp

def process_snp_dists(input_path: str,output_path: str, n_jobs: int = 20):
    """处理SNP距离计算并生成三角矩阵结果
    
    Args:
        input_path: 输入文件路径
        n_jobs: 并行任务数 (默认4)
    """
    # 生成FASTA格式数据
    df = pl.read_csv(input_path, separator="\t", row_index_name='pos')
    result = df.select(pl.all().str.join())
    fa=[f'>{result[series].name}\n{result[series].item()}\n' for series in result.columns]
    
    with open(f'{output_path}.fa', 'w') as f:
        f.write(''.join(fa[2:]))  # 跳过前两列
    
    # 执行snp-dists计算
    sp.run(f"snp-dists -j {n_jobs} {output_path}.fa > {output_path}.fa.dist", shell=True)
    
    # 处理输出结果

    df = pl.read_csv(f"{output_path}.fa.dist", separator="\t", has_header=True)
    row_names = df.columns[1:]
    output_path = f"{output_path}.fa.dist.tr"
    # 构建三角矩阵数据
    data = []
    for i in range(len(row_names)):
        for j in range(i + 1, len(row_names)):
            data.append({
                "strain1": row_names[i],
                "strain2": row_names[j],
                "distance": df[row_names[i]][j]
            })
    
    # 写入输出文件
    pl.DataFrame(data).write_csv(output_path, separator="\t")
    return output_path
