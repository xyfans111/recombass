import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
import polars.selectors as cs
from numba import jit



@jit
def pmr_gap(x, y):
    return x * (y - 1) / 2 + x * (y - x) / 2

@jit
def pmr_snp(x, y):
    return (x * x - x) / 2

cal_map = {
    "A": pmr_snp,
    "T": pmr_snp,
    "C": pmr_snp,
    "G": pmr_snp,
    "-": pmr_gap
}

def cal(x, y):
    return cal_map[x[0]](x[1], y)

def pmr_us(df, lenth):
    counts = df.value_counts()
    dist = [cal(x, lenth) for x in counts.items()]
    return np.sum(dist)

def R_dist(x, y):
    if y==x.es:
        return np.nan
    return -1 * (x.pmn - x.es) / (y - x.es)

def remove_top_bottom_10_percent(arr):
    arr = np.array(arr)  # 确保是 numpy 数组
    lower = np.percentile(arr, 10)
    upper = np.percentile(arr, 90)
    return arr[(arr > lower) & (arr < upper)]

def cal_recx(values):
    """
    计算 recx 值
    
    Args:
        values: list 或 array 类型的数据
        
    Returns:
        计算得到的 recx 值
    """
    values = np.array(values)  # 确保是 numpy 数组
    
    # 移除 NaN 值
    values = values[~np.isnan(values)]
    
    # 如果没有有效数据，返回 NaN
    if len(values) == 0:
        return np.nan
    
    # 移除前10%和后10%的数据
    values = remove_top_bottom_10_percent(values)
    
    if len(values) > 0:
        x = np.std(values)
        return np.exp(-10.106*x + 3.388)
    else:
        return np.nan
# 提取核心逻辑为函数
def process_snp_data(oripath,snppath):
    df2 = pl.read_csv(f'{snppath}.fa.dist.tr', separator='\t')
    df2 = df2.with_columns(pl.col('strain1').cast(pl.String))
    df2 = df2.with_columns(pl.col('strain2').cast(pl.String))
    var = df2['distance'].std()
    med = df2['distance'].median()
    cut = 0.9 * med

    sns.displot(df2['distance'].to_numpy())
    plt.axvline(x=cut, c='#FFAC00')
    plt.axvline(x=cut * 0.30, c='#FFAC00')
    plt.axvline(x=cut * 0.60, c='#FFAC00')
    plt.axvline(x=cut * 0.90, c='#FFAC00')
    plt.savefig(f"{snppath}.fa.dist.png", dpi=300)

    df = df2.filter(pl.col('distance') > cut)
    df2 = df2.filter(pl.col('distance') < cut)
    udist = df['distance'].mean()

    data = pl.read_csv(f'{oripath}', separator='\t').select(cs.by_dtype(pl.String).cast(pl.Categorical))
    dft = data.to_pandas().T
    lenth = dft.shape[0]
    pmr_list = np.array([pmr_us(dft[x], lenth) for x in dft.columns])
    list4 = dft.columns.to_list()
    dft = 0

    # 处理不同距离范围的数据
    def process_distance_range(df2, cut, lower, upper):
        df = df2.filter((lower * cut < pl.col('distance')) & (pl.col('distance') <= upper * cut))
        if df.shape[0] == 0:
            return 0, 0, 0
        mdist = df['distance'].mean() 
        list_ur = tuple(map(tuple, df.to_numpy()))
        # Create expressions for each comparison
        expressions = [((data[:, list_ur[t][0]] == data[:, list_ur[t][1]]) |
                   (data[:, list_ur[t][0]] == '-') |
                   (data[:, list_ur[t][1]] == '-')).alias(f'{t}') 
                   for t in range(len(list_ur))]
    
        # Use sum_horizontal for better performance
        p = data.select(expressions).sum_horizontal().to_numpy()
        l = len(list_ur)
        return p, l, mdist

    p3, l3, mdist3 = process_distance_range(df2, cut, 0.6, 0.9)
    p2, l2, mdist2 = process_distance_range(df2, cut, 0.3, 0.6)
    p1, l1, mdist1 = process_distance_range(df2, cut, 0.01, 0.3)
    df = df2.filter((0.01*cut>=pl.col('distance'))|(pl.col('distance')>0.9*cut))
    list_ur=tuple(map(tuple,df.to_numpy()))

    p = data.select([((data[:,list_ur[t][0]]==data[:,list_ur[t][1]])|(data[:,list_ur[t][0]]=='-')|(data[:,list_ur[t][1]]=='-')).alias(f'{t}') for t in range(len(list_ur))]).sum_horizontal().to_numpy()
    l =len(list_ur) 

    # 合并结果
    p=p+p1+p2+p3
    l=l+l1+l2+l3
    lenur=lenth*(lenth-1)/2-l
    pmrur=pmr_list-p
    pmr_urlist=pmrur/lenur
    p4=p1+p2+p3
    l4=l1+l2+l3
    with open(f'{snppath}.fa.pmrus.txt','w') as f:
        f.write('\n'.join([str(x) for x in pmr_urlist.tolist()]))
    def f(x):
        return (1+mdist1/udist*(x-1))*l1
    def g(x,y):
        if y==x.pmn:
            return np.nan
        return (x.es-x.pmn)/(y-x.pmn)
    es1=f(pmr_urlist)
    dfr = pd.DataFrame(zip(p1, es1), columns=['pmn', 'es'])
    list1 = dfr.apply(R_dist, y=l1, axis=1).values.tolist()
    list1r=dfr.apply(g, y=l1, axis=1).values.tolist()       
    def f3(x): return (1+(mdist1*l1+mdist2*l2+mdist3*l3)/l4/udist*(x-1))*l4
    es4=f3(pmr_urlist)
    dfr = pd.DataFrame(zip(p4, es4), columns=['pmn', 'es'])
    with open(f'{snppath}.fa.pmr.all.ox.txt', 'w') as file:
            file.write('\n'.join([str(x) for x in p4]))
    with open(f'{snppath}.fa.pmr.all.ex.txt', 'w') as file:
            file.write('\n'.join([str(x) for x in es4]))
    with open(f'{snppath}.fa.pmr.30.ox.txt', 'w') as file:
            file.write('\n'.join([str(x) for x in p1]))
    with open(f'{snppath}.fa.pmr.30.ex.txt', 'w') as file:
            file.write('\n'.join([str(x) for x in es1]))
    list5 = dfr.apply(R_dist, y=l4, axis=1).values.tolist()
    list5r=dfr.apply(g, y=l4, axis=1).values.tolist()
    # 保存结果到文件
    def save_results(list_r, list_gr, output_file):
        dfa = pd.DataFrame(zip(list4, list_r), columns=['pos', 'r'])
        dfa = dfa[dfa.r <= 0]
        dfb = pd.DataFrame(zip(list4, list_gr), columns=['pos', 'r'])
        dfb = dfb[dfb.r >= 0]
        res = pd.concat([dfa, dfb]).sort_values('pos')
        all_positions = pd.DataFrame({'pos': list4})
        
        # 合并所有位置与筛选后的结果，用0填充缺失值
        res_complete = all_positions.merge(res, on='pos', how='left')
        res_complete['r'] = res_complete['r'].fillna(0)
        res_complete = res_complete.sort_values('pos')
        
        with open(output_file, 'w') as f:
            f.write('\n'.join([str(x) for x in res_complete['r'].values]))
        return res['r'].values.tolist()
    save_results(list1, list1r, f'{snppath}.fa.pmr.30.2.txt')
    res=save_results(list5, list5r, f'{snppath}.fa.pmr.all.2.txt')
    print(f'{snppath} recx = {cal_recx(res)}')

