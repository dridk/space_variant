import polars as pl
import pysam
from tqdm import tqdm
from pathlib import Path


def extract_read(bam_file: str, destination: str = None):
    """
    Lecture des spatial variants

    """

    if destination is None:
        path = Path.cwd() / "output"
        path.mkdir(exist_ok=True)

    else:
        path = Path(destination)
        if not path.is_dir():
            raise Exception(f"{path} is not a directory")

    bam = pysam.AlignmentFile(bam_file)

    chromosomes = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]

    for current_chrom in chromosomes:
        print(current_chrom)
        data = []
        for read in tqdm(bam.fetch(current_chrom)):
            if read.has_tag("CB") and read.has_tag("UB"):
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end
                length = read.reference_length
                barcode = read.get_tag("CB")
                _, _, x, y = barcode.split("_")
                y = y.replace("-1", "")
                x = int(x)
                y = int(y)
                umi = read.get_tag("UB")

                data.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "length": length,
                        "barcode": barcode,
                        "x": x,
                        "y": y,
                        "umi": umi,
                    }
                )

        print("create dataframe")
        df = pl.DataFrame(data)
        print("create parquet file")
        df.write_parquet(path / f"reads_{current_chrom}.parquet")


if __name__ == "__main__":
    extract_read("../../../../dataset/pancreas/my_output/outs/possorted_genome_bam.bam")
