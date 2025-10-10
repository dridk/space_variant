import polars as pl
import gzip
import pysam

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

    df.columns = [i.lower() for i in df.columns]
    return df


def read_spatial_variant(bam_file, vcf_file):

    bam_file = "../../dataset/pancreas/my_output/outs/possorted_genome_bam.bam"

    vdf = read_vcf(vcf_file)
    bam = pysam.AlignmentFile(bam_file)

    for row in vdf.rows(named=True):

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
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "dp": dp,
                            "ad": ad,
                            "af": af,
                            "gt": gt,
                            "base": base,
                            "barcode": barcode,
                            "x": x,
                            "y": y,
                            "umi": umi,
                        }


bam_file = "../../dataset/pancreas/my_output/outs/possorted_genome_bam.bam"
vcf_file = (
    "../../dataset/pancreas/Visium_HD_3prime_Human_Pancreatic_Cancer_fastqs/variants.mutect.vcf.gz"
)

for rec in read_spatial_variant(bam_file=bam_file, vcf_file=vcf_file):

    print(rec)
