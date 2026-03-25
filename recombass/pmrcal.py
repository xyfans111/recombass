import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl
import polars.selectors as cs


# Calculate PMR for gaps
def pmr_gap(x, y):
    return x * (y - 1) / 2 + x * (y - x) / 2


# Calculate PMR for SNPs
def pmr_snp(x, _):
    return (x * x - x) / 2


# Mapping of nucleotides to PMR calculation functions
cal_map = {
    "A": pmr_snp,
    "T": pmr_snp,
    "C": pmr_snp,
    "G": pmr_snp,
    "-": pmr_gap,
}


# Calculate PMR based on nucleotide type
def cal(x, y):
    return cal_map[x[0]](x[1], y)


# Calculate PMR for unique sequences
def pmr_us(df, lenth):
    counts = df.value_counts()
    dist = [cal(x, lenth) for x in counts.items()]
    return np.sum(dist)


# Calculate R distance
def R_dist(x, y):
    if y == x.es:
        return np.nan
    return -1 * (x.pmn - x.es) / (y - x.es)


# Remove top and bottom 10% of values
def remove_top_bottom_10_percent(arr):
    arr = np.array(arr)
    lower = np.percentile(arr, 10)
    upper = np.percentile(arr, 90)
    return arr[(arr > lower) & (arr < upper)]


# Calculate recx value
def cal_recx(values):
    values = np.array(values)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return np.nan

    values = remove_top_bottom_10_percent(values)
    if len(values) > 0:
        x = np.std(values)
        return np.exp(-10.106 * x + 3.388)
    return np.nan


# Process SNP data
def process_snp_data(oripath, snppath):
    df2 = pl.read_csv(f"{snppath}.fa.dist.tr", separator="\t")
    df2 = df2.with_columns(pl.col("strain1").cast(pl.String))
    df2 = df2.with_columns(pl.col("strain2").cast(pl.String))
    med = df2["distance"].median()
    cut = 0.9 * med
    plt.ioff()
    sns.displot(df2["distance"].to_numpy())
    plt.axvline(x=cut, c="#FFAC00")
    plt.axvline(x=cut * 0.30, c="#FFAC00")
    plt.axvline(x=cut * 0.90, c="#FFAC00")
    plt.savefig(f"{snppath}.fa.dist.png", dpi=300)
    plt.close()
    df = df2.filter(pl.col("distance") > cut)
    df2 = df2.filter(pl.col("distance") < cut)
    udist = df["distance"].mean()
    if udist is None or np.isnan(udist) or udist == 0:
        udist = med if med is not None and not np.isnan(med) and med != 0 else 1.0

    data = pl.read_csv(f"{oripath}", separator="\t").select(cs.by_dtype(pl.String).cast(pl.Categorical))
    dft = data.to_pandas().T
    lenth = dft.shape[0]
    num_positions = data.height
    pmr_list = np.array([pmr_us(dft[x], lenth) for x in dft.columns])
    list4 = dft.columns.to_list()
    dft = 0

    def process_distance_range(df2, cut, lower, upper):
        df = df2.filter((lower * cut < pl.col("distance")) & (pl.col("distance") <= upper * cut))
        if df.shape[0] == 0:
            return np.zeros(num_positions, dtype=float), 0, np.nan
        mdist = df["distance"].mean()
        list_ur = tuple(map(tuple, df.to_numpy()))
        expressions = [
            (
                (data[:, list_ur[t][0]] == data[:, list_ur[t][1]])
                | (data[:, list_ur[t][0]] == "-")
                | (data[:, list_ur[t][1]] == "-")
            ).alias(f"{t}")
            for t in range(len(list_ur))
        ]

        p = data.select(expressions).sum_horizontal().to_numpy()
        l = len(list_ur)
        return p, l, mdist

    p3, l3, mdist3 = process_distance_range(df2, cut, 0.6, 0.9)
    p2, l2, mdist2 = process_distance_range(df2, cut, 0.3, 0.6)
    p1, l1, mdist1 = process_distance_range(df2, cut, 0.01, 0.3)
    df = df2.filter((0.01 * cut >= pl.col("distance")) | (pl.col("distance") > 0.9 * cut))
    list_ur = tuple(map(tuple, df.to_numpy()))

    p = data.select(
        [
            (
                (data[:, list_ur[t][0]] == data[:, list_ur[t][1]])
                | (data[:, list_ur[t][0]] == "-")
                | (data[:, list_ur[t][1]] == "-")
            ).alias(f"{t}")
            for t in range(len(list_ur))
        ]
    ).sum_horizontal().to_numpy()
    l = len(list_ur)

    p = p + p1 + p2 + p3
    l = l + l1 + l2 + l3
    lenur = lenth * (lenth - 1) / 2 - l
    pmrur = pmr_list - p
    pmr_urlist = pmrur / lenur
    p4 = p1 + p2 + p3
    l4 = l1 + l2 + l3

    def f(x):
        if l1 == 0 or mdist1 is None or np.isnan(mdist1):
            return np.zeros_like(x, dtype=float)
        return (1 + mdist1 / udist * (x - 1)) * l1

    def g(x, y):
        if y == x.pmn:
            return np.nan
        return (x.es - x.pmn) / (y - x.pmn)

    es1 = f(pmr_urlist)
    dfr = pd.DataFrame(zip(p1, es1), columns=["pmn", "es"])
    list1 = dfr.apply(R_dist, y=l1, axis=1).values.tolist()
    list1r = dfr.apply(g, y=l1, axis=1).values.tolist()

    def f3(x):
        if l4 == 0:
            return np.zeros_like(x, dtype=float)
        return (1 + (mdist1 * l1 + mdist2 * l2 + mdist3 * l3) / l4 / udist * (x - 1)) * l4

    es4 = f3(pmr_urlist)
    dfr = pd.DataFrame(zip(p4, es4), columns=["pmn", "es"])
    list5 = dfr.apply(R_dist, y=l4, axis=1).values.tolist()
    list5r = dfr.apply(g, y=l4, axis=1).values.tolist()

    def combine_signed_results(negative_values, positive_values):
        dfa = pd.DataFrame(zip(list4, negative_values), columns=["pos", "r"])
        dfa = dfa[dfa.r <= 0]
        dfb = pd.DataFrame(zip(list4, positive_values), columns=["pos", "r"])
        dfb = dfb[dfb.r >= 0]
        res = pd.concat([dfa, dfb]).sort_values("pos")
        all_positions = pd.DataFrame({"pos": list4})
        res_complete = all_positions.merge(res, on="pos", how="left")
        res_complete["r"] = res_complete["r"].fillna(0)
        return res_complete.sort_values("pos").reset_index(drop=True)

    df_30 = combine_signed_results(list1, list1r).rename(columns={"r": "30_pmr"})
    df_all = combine_signed_results(list5, list5r).rename(columns={"r": "all_pmr"})
    df_us = pd.DataFrame({"pos": list4, "us_pmr": pmr_urlist.tolist()}).sort_values("pos").reset_index(drop=True)
    result_df = df_30.merge(df_all, on="pos").merge(df_us, on="pos")

    recx = cal_recx(result_df["all_pmr"].to_numpy())
    with open(f"{snppath}.recx.txt", "w") as a:
        a.write(str(recx))
    print(f"{snppath} recx = {recx}")

    return result_df
