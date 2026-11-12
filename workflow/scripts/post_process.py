"""
Post-traitement des fichiers {sample}.{chr}.spatial_variant.parquet

1. Enrichissement somatique : retrait des variants présents dans dbSNP common
   (rsid renseigné) et profondeur strictement supérieure à --min-depth.
2. Score de clustering spatial par variant :

       R(V) = sqrt( 1/N * sum_i || c_i - c_mean ||^2 )
       s(V) = R(V) / R_all

   où les c_i sont les N spots distincts portant le variant V, et R_all le
   rayon de giration de tous les spots portant au moins un variant (avant
   filtrage). s ~ 1 : variant réparti comme la coupe ; s -> 0 : variant
   concentré.
"""

import polars as pl

VARIANT_KEY = ["chr", "pos", "ref", "alt"]
ANNOTATIONS = ["gene", "effect", "impact", "hgvs_c", "hgvs_p", "dp", "ao", "af"]


def gyration_radius() -> pl.Expr:
    """Rayon de giration d'un nuage de points (x, y)."""
    return (pl.col("x").var(ddof=0) + pl.col("y").var(ddof=0)).sqrt()


def filter_somatic(df: pl.LazyFrame, min_depth: int) -> pl.LazyFrame:
    """Variants absents de dbSNP common et de profondeur > min_depth."""
    not_in_dbsnp = pl.col("rsid").is_null() | (pl.col("rsid") == ".")
    return df.filter(not_in_dbsnp & (pl.col("dp") > min_depth))


def gyration_score(somatic: pl.LazyFrame, r_all: float) -> pl.LazyFrame:
    """Une ligne par variant, avec son rayon de giration et son score s."""
    molecules = somatic.group_by(VARIANT_KEY).agg(
        pl.struct("barcode", "umi").n_unique().alias("n_molecules")
    )
    spots = (
        somatic.unique(VARIANT_KEY + ["x", "y"])
        .group_by(VARIANT_KEY)
        .agg(
            *[pl.col(c).first() for c in ANNOTATIONS],
            pl.len().alias("n_spots"),
            pl.col("x").mean().alias("x_mean"),
            pl.col("y").mean().alias("y_mean"),
            gyration_radius().alias("radius"),
        )
    )
    return (
        spots.join(molecules, on=VARIANT_KEY)
        .with_columns((pl.col("radius") / r_all).alias("s"))
        .sort("s", "n_spots", descending=[False, True])
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Filter somatic-like variants and compute the gyration score s."
    )
    parser.add_argument("parquet_files", nargs="+", help="spatial_variant parquet files")
    parser.add_argument("--somatic", required=True, help="Output parquet: filtered molecules")
    parser.add_argument("--score", required=True, help="Output parquet: one row per variant with s")
    parser.add_argument("--min-depth", type=int, default=40, help="Keep variants with DP > min-depth")
    args = parser.parse_args()

    df = pl.scan_parquet(args.parquet_files)

    # Rayon de giration de référence : tous les spots portant un variant
    r_all = df.select("x", "y").unique().select(gyration_radius()).collect().item()

    somatic = filter_somatic(df, args.min_depth).collect()
    somatic.write_parquet(args.somatic)

    score = gyration_score(somatic.lazy(), r_all).collect()
    score.write_parquet(args.score)

    n_all = df.select(VARIANT_KEY).unique().select(pl.len()).collect().item()
    print(f"R_all = {r_all:.2f} (spot units)")
    print(f"{len(score)} / {n_all} variants kept (not in dbSNP, DP > {args.min_depth})")
    print(score.head(10))
