import polars as pl


VALID_BASES = {"A", "T", "C", "G"}


def normalize_snp_matrix(input_path: str) -> tuple[pl.DataFrame, str, list[str]]:
    """Load a SNP matrix and coerce all non-ATCG symbols to gap."""
    df = pl.read_csv(input_path, separator="\t")
    position_col = df.columns[0]
    sequence_cols = df.columns[1:]

    expressions = [
        pl.when(pl.col(col).is_in(VALID_BASES))
        .then(pl.col(col))
        .otherwise(pl.lit("-"))
        .alias(col)
        for col in sequence_cols
    ]

    df = df.with_columns(expressions)
    return df, position_col, sequence_cols
