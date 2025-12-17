#!/usr/bin/env python
import numpy as np
import polars as pl
def filter_maf(input_path: str, maf_rate: float, output_path: str,) -> str:
    """基于次要等位基因频率(MAF)过滤序列数据"""
    df = pl.read_csv(input_path, separator="\t")
    
    array = df[:,1:].to_numpy().T
    
    maf_value = array.shape[0] * maf_rate
    
    # 使用 numpy.apply_along_axis 对每列应用计算函数
    def calculate_maf_for_column(column_data):
        # 使用 np.unique 获取唯一值和计数
        unique_vals, counts = np.unique(column_data, return_counts=True)
        
        if len(counts) > 0:
            # 获取最大计数
            max_count = np.max(counts)
            
            # 使用布尔索引查找 '-'（更高效）
            gap_mask = unique_vals == '-'
            gap_count = counts[gap_mask][0] if np.any(gap_mask) else 0
            
            return max_count + gap_count
        else:
            return 0
    
    # 对每列应用函数
    maf_counts = np.apply_along_axis(calculate_maf_for_column, axis=0, arr=array)
    
    array1 = np.where(maf_counts <= maf_value)
    
    filtered_df = df[array1[0].tolist(),:]
    output_path = f'{output_path}.maf{int(maf_rate*100):02d}'
    filtered_df.write_csv(output_path, separator="\t", include_header=True)
    
    return output_path

