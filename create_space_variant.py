import polars as pl
import gzip

filename = (
    "../../dataset/pancreas/Visium_HD_3prime_Human_Pancreatic_Cancer_fastqs/variants.mutect.vcf.gz"
)


def read_vcf(filename):

    def xopen(filename):
        if filename.endswith(".gz"):
            file = gzip.open(filename, "r")
        else:
            file = open(filename, "rb")
        return file

    # Read header count
    with xopen(filename) as file:
        header_count = 0
        for line in file:
            if line.decode().startswith("#"):
                header_count += 1
                header = line
            else:
                break

    header = header.decode().strip().split("\t")
    sample = header[-1]
    columns = ["CHROM", "POS", "REF", "ALT"]
    header_count -= 1
    df = (
        pl.scan_csv(
            filename,
            skip_rows=header_count,  # Je saute les 40 premieres lignes de commentaires
            separator="\t",  # Séparateur TSV
            dtypes={"#CHROM": pl.Utf8},
            with_column_names=lambda cols: [c.replace("#", "") for c in cols],
        )
        .select(
            columns
            + [  # Je sélectionne les colonnes souhaitées
                pl.col("FORMAT").str.split(":"),
                pl.col("Boby").str.split(":"),
            ]
        )
        .explode(["FORMAT", "Boby"])
        .filter(pl.col("FORMAT").is_in(["GT", "AD", "AF", "DP"]))
        .collect()
        .pivot(on="FORMAT", index=columns, values=[sample])
        .with_columns(pl.lit(sample).alias("SAMPLE"))
    )

    return df


read_vcf(filename)


# print(df)
