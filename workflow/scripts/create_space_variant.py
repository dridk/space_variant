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
                pl.col("Boby").str.split(":"),
            ]
        )
        .explode(["FORMAT", "Boby"])
        .filter(pl.col("FORMAT").is_in(["GT", "AD", "AF", "DP"]))
        .collect()
        .pivot(on="FORMAT", index=columns, values=[sample])
        .with_columns(pl.lit(sample).alias("SAMPLE"))
    )

    df.columns = [i.lower() for i in df.columns]

    df = df.with_columns(
        pl.col("af").cast(pl.Float32()),
        pl.col("dp").cast(pl.UInt32()),
        pl.col("id").str.replace(".", None),
    )

    return df


def read_spatial_variant(
    bam_file: str, vcf_file: str, min_dp: int = 20, min_af: float = 0.1
):
    """
    Lecture des spatial variants

    """
    vdf = (
        read_vcf(vcf_file)
        .filter(pl.col("af").is_between(min_af, 1))
        .filter(pl.col("dp") > min_dp)
    )

    vdf = vdf.filter(pl.col("chrom").str.contains(r"^chr[\dXY]$"))

    vdf.write_csv("pancreas.csv", separator="\t")
    bam = pysam.AlignmentFile(bam_file)

    for row in tqdm(vdf.rows(named=True), total=len(vdf)):
        chrom = row["chrom"]
        pos = row["pos"]
        start = row["pos"] - 1
        end = row["pos"]
        ref = row["ref"]
        alt = row["alt"]
        gt = row["gt"]
        ad = row["ad"]
        dp = row["dp"]
        af = row["af"]
        info = row["info"]

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

                        yield {
                            "rsid": id,
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "dp": dp,
                            "ad": ad,
                            "af": af,
                            "gt": gt,
                            "info": info,
                            "base": base,
                            "barcode": barcode,
                            "x": x,
                            "y": y,
                            "umi": umi,
                        }


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
    df.write_parquet(args.output_file)
