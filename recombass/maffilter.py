#!/usr/bin/env python
import numpy as np

from .inputprep import normalize_snp_matrix

def filter_maf(input_path: str, maf_rate: float, output_path: str) -> str:
    """Filter sequence data based on Minor Allele Frequency (MAF)"""
    df, _, _ = normalize_snp_matrix(input_path)
    
    array = df[:,1:].to_numpy().T
    
    maf_value = array.shape[0] * maf_rate
    
    # Use numpy.apply_along_axis to apply calculation function to each column
    def calculate_maf_for_column(column_data):
        # Use np.unique to get unique values and counts
        unique_vals, counts = np.unique(column_data, return_counts=True)
        
        if len(counts) > 0:
            # Get maximum count
            max_count = np.max(counts)
            
            # Use boolean indexing to find '-' (more efficient)
            gap_mask = unique_vals == '-'
            gap_count = counts[gap_mask][0] if np.any(gap_mask) else 0
            
            return max_count + gap_count
        else:
            return 0
    
    # Apply function to each column
    maf_counts = np.apply_along_axis(calculate_maf_for_column, axis=0, arr=array)
    
    array1 = np.where(maf_counts <= maf_value)
    
    filtered_df = df[array1[0].tolist(),:]
    output_path = f'{output_path}.maf{int(maf_rate*100):02d}'
    filtered_df.write_csv(output_path, separator="\t", include_header=True)
    
    return output_path
