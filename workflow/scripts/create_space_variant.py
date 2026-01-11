import polars as pl
import gzip
import pysam
from tqdm import tqdm


def read_vcf(filename):
    """
    Lecture d'un fichier VCF en polars
    """

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
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "INFO"]
    header_count -= 1
    df = (
        pl.scan_csv(
            filename,
            skip_rows=header_count,  # Je saute les 40 premieres lignes de commentaires
            separator="\t",  # Séparateur TSV
            schema_overrides={"#CHROM": pl.Utf8},
            with_column_names=lambda cols: [c.replace("#", "") for c in cols],
        )
        .select(
            columns
            + [  # Je sélectionne les colonnes souhaitées
                pl.col("FORMAT").str.split(":"),
                pl.col(sample).str.split(":"),
            ]
        )
        .explode(["FORMAT", sample])
        .filter(pl.col("FORMAT").is_in(["GT", "AD", "AO", "DP", "RO"]))
        .collect()
        .pivot(on="FORMAT", index=columns, values=[sample])
        .with_columns(pl.lit(sample).alias("SAMPLE"))
    )

    df.columns = [i.lower() for i in df.columns]

    df = df.with_columns(
        pl.col("ao").cast(pl.UInt32()),
        pl.col("dp").cast(pl.UInt32()),
        pl.col("ad").str.split(","),
        pl.col("id").str.replace(".", None),
    ).with_columns((pl.col("ao") / pl.col("dp")).round(2).alias("af"))

    return df


def read_spatial_variant(
    bam_file: str, vcf_file: str, min_dp: int = 0, min_af: float = 0.0
):
    """
    Lecture des spatial variants

    """
    vdf = read_vcf(vcf_file)

    vdf = vdf.filter(pl.col("dp") >= min_dp).filter(pl.col("af") >= min_af)
    vdf = vdf.filter(pl.col("chrom").str.contains(r"^chr[\dXY]+$"))
    vdf = vdf.with_columns(
        (pl.col("pos") - 1).alias("start"), (pl.col("pos")).alias("end")
    )

    # vdf = vdf.filter(pl.col("id").str.starts_with("rs"))

    bam = pysam.AlignmentFile(bam_file)

    for index, row in enumerate(tqdm(vdf.rows(named=True), total=len(vdf))):
        info = row["info"]
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        ref = row["ref"]

        # if index > 50:
        #     break
        for colpil in bam.pileup(chrom, start, end, truncate=True):
            for read in colpil.pileups:
                if (
                    not read.is_del
                    and not read.is_refskip
                    and read.alignment.has_tag("CB")
                    and read.alignment.has_tag("UB")
                ):
                    base = read.alignment.query_sequence[read.query_position]

                    if base != ref:
                        barcode = read.alignment.get_tag("CB")
                        _, _, x, y = barcode.split("_")
                        y = y.replace("-1", "")
                        x = int(x)
                        y = int(y)

                        umi = read.alignment.get_tag("UB")

                        # split annotation
                        ann = []
                        ann_keys = [
                            "allele",
                            "effect",
                            "impact",
                            "gene",
                            "gene_id",
                            "feature_type",
                            "feature_id",
                            "biotype",
                            "rank",
                            "hgvs_c",
                            "hgvs_p",
                            "dna",
                            "cds",
                            "protein",
                            "distance",
                            "_",
                        ]
                        for i in info.split(";"):
                            if i.startswith("ANN="):
                                ann = i.split("=")[1].split("|")
                                ann = dict(zip(ann_keys, ann))
                                break

                        item = {}
                        item["chr"] = chrom
                        item["start"] = start
                        item["end"] = end
                        item["pos"] = row["pos"]
                        item["ref"] = row["ref"]
                        item["alt"] = row["alt"]
                        item["rsid"] = row["id"]
                        item["gt"] = row["gt"]
                        item["ad"] = row["ad"]
                        item["dp"] = row["dp"]
                        item["ao"] = row["ao"]
                        item["ro"] = row["ro"]
                        item["af"] = row["af"]
                        item["gene"] = ann["gene"]
                        item["effect"] = ann["effect"]
                        item["feature_type"] = ann["feature_type"]
                        item["impact"] = ann["impact"]
                        item["biotype"] = ann["biotype"]
                        item["hgvs_c"] = ann["hgvs_c"]
                        item["hgvs_p"] = ann["hgvs_p"]

                        item["base"] = base
                        item["barcode"] = barcode
                        item["x"] = x
                        item["y"] = y
                        item["umi"] = umi

                        yield item


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process BAM and VCF files to produce a Parquet output."
    )
    parser.add_argument("bam_file", type=str, help="Input BAM file")
    parser.add_argument("vcf_file", type=str, help="Input VCF file")
    parser.add_argument("output_file", type=str, help="Output Parquet file")
    parser.add_argument(
        "--min-depth", type=int, default=0, help="Minimum depth for filtering"
    )
    parser.add_argument(
        "--min-af", type=float, default=0.0, help="Minimum allele frequency"
    )

    args = parser.parse_args()

    # bam_file = "../../dataset/pancreas/my_output/outs/possorted_genome_bam.bam"
    # vcf_file = "../../dataset/pancreas/Visium_HD_3prime_Human_Pancreatic_Cancer_fastqs/pancreas.mutect.ann.dbsnp.norm.vcf.gz"

    df = pl.DataFrame(
        read_spatial_variant(
            bam_file=args.bam_file,
            vcf_file=args.vcf_file,
            min_dp=args.min_depth,
            min_af=args.min_af,
        )
    )

    print(df)
    df.write_parquet(args.output_file)

