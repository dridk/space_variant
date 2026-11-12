# space_variant

Snakemake pipeline that calls variants from **Visium HD 3' (polyA)** spatial
transcriptomics alignments.   
The output is a set of Parquet tables (one per chromosome) in which each row is
a variant called with its annotation and its position on the slide.
A post-processing step then keeps the somatic-like variants and computes a
spatial-clustering score for each of them.


1. [Download the pancreas dataset](#1-download-the-pancreas-dataset)
2. [Re-run Space Ranger to get the BAM](#2-re-run-space-ranger-to-get-the-bam)
3. [Run the pipeline](#3-run-the-pipeline)
4. [Results: the Parquet files](#4-results-the-parquet-files)
5. [Post-processing: somatic filter and clustering score](#5-post-processing-somatic-filter-and-clustering-score)

## Requirements

| Tool | Version used | Purpose |
|------|--------------|---------|
| [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/downloads) | 4.0.1 | alignment, BAM with spatial barcodes |
| snakemake | 9.13 | workflow engine |
| samtools / htslib (`bgzip`, `tabix`) | 1.22 | BAM handling |
| freebayes | 1.3.10 | variant calling |
| bcftools | 1.22 | splitting multi-allelic sites |
| SnpSift / snpEff | 5.3 | dbSNP and functional annotation |
| python + `polars`, `pysam`, `tqdm` | 3.12 | building the Parquet tables |

Everything except Space Ranger can be installed in one conda environment:

```bash
mamba create -n spacevariant -c conda-forge -c bioconda \
    snakemake-minimal samtools htslib bcftools freebayes snpeff snpsift \
    polars pysam tqdm
mamba activate spacevariant
```

Space Ranger must be downloaded from the 10x Genomics website.

- [Download spaceRanger](https://www.10xgenomics.com/support/software/space-ranger/latest)


## 1. Download the pancreas dataset

Dataset: *Visium HD 3' Gene Expression Library, Human Pancreatic Cancer (Fresh
Frozen)*,
<https://www.10xgenomics.com/datasets/visium-hd-three-prime-human-pancreatic-cancer-fresh-frozen>.

10x Genomics does not distribute the BAM file for this dataset, so the raw
reads and the CytAssist image are needed to re-run Space Ranger.


Then the references:

```bash
mkdir -p hg38 && cd hg38

# 10x GRCh38 reference (11 GB): used by Space Ranger AND as the pipeline FASTA
curl -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz

# dbSNP build 151, common variants, GRCh38 (1.5 GB)
curl -O https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
curl -O https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi
cd ..
```

## 2. Re-run Space Ranger to get the BAM

The only non-default option that matters is `--create-bam=true`. The slide
serial number of this dataset is `H1-CWRFYX3`, capture area `A1`.

```bash
spaceranger count \
    --id=pancreas \
    --fastqs=dataset/pancreas/Visium_HD_3prime_Human_Pancreatic_Cancer_fastqs/ \
    --cytaimage=dataset/pancreas/Visium_HD_3prime_Human_Pancreatic_Cancer_image.tif \
    --slide=H1-CWRFYX3 \
    --area=A1 \
    --transcriptome=hg38/refdata-gex-GRCh38-2024-A \
    --create-bam=true \
    --output-dir=dataset/pancreas/spaceranger
```


Each aligned read in the BAM carries a `CB` tag (the 2 µm spatial barcode,
e.g. `s_002um_00334_02282-1`) and a `UB` tag (the UMI). The pipeline relies on
both; reads without them are ignored.

## 3. Run the pipeline

Clone the repository and link the BAM into the working directory `data/` under
the name `pancreas.bam`

    
```bash
git clone https://github.com/dridk/space_variant.git
cd space_variant
```

Create the configuration file `workflow/config.yml`

```yaml
# FASTA of the Space Ranger reference, so that contig names match the BAM
reference: "/abs/path/to/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
# dbSNP common variants (bgzipped + tabix-indexed)
dbsnp: "/abs/path/to/hg38/00-common_all.vcf.gz"
# snpEff database; downloaded automatically on first use
snpeff: "GRCh38.86"

# freebayes thresholds
min_dp: 10   # --min-coverage
min_ao: 3    # --min-alternate-count

# post-processing: keep variants absent from dbSNP common with DP > somatic_min_dp
somatic_min_dp: 40
```

Run:

```bash
# dry run: show what will be executed
snakemake -n -p -s workflow/Snakefile --configfile workflow/config.yml -d data --cores 20

# real run (equivalent to `make run`)
snakemake -p -s workflow/Snakefile --configfile workflow/config.yml -d data --cores 20 --rerun-incomplete
```

## 4. Results: the Parquet files

The pipeline produces 24 files, `data/pancreas.{chr}.spatial_variant.parquet`


| Column | Type | Description |
|--------|------|-------------|
| `chr` | str | Chromosome (`chr1` … `chr22`, `chrX`, `chrY`) |
| `start` | int | 0-based position of the variant (`pos - 1`) |
| `end` | int | 0-based exclusive end (equal to `pos`) |
| `pos` | int | 1-based position, as in the VCF |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele (one per row after `bcftools norm`) |
| `rsid` | str | dbSNP rsID if the variant is a common dbSNP allele, `.` otherwise |
| `gt` | str | freebayes genotype on the whole section (pseudo-bulk), e.g. `0/1` |
| `ad` | list[str] | Allelic depths `[ref, alt]` (pseudo-bulk) |
| `dp` | int | Total read depth at the site (pseudo-bulk) |
| `ao` | int | Number of reads supporting `alt` (pseudo-bulk) |
| `ro` | str | Number of reads supporting `ref` (pseudo-bulk) |
| `af` | float | Allele fraction `ao / dp`, rounded to 2 decimals |
| `gene` | str | Gene symbol (first snpEff annotation) |
| `effect` | str | snpEff consequence, e.g. `missense_variant` |
| `feature_type` | str | snpEff feature type, e.g. `transcript` |
| `impact` | str | snpEff impact: `HIGH`, `MODERATE`, `LOW`, `MODIFIER` |
| `biotype` | str | Transcript biotype, e.g. `protein_coding` |
| `hgvs_c` | str | HGVS coding notation, e.g. `c.146A>C` |
| `hgvs_p` | str | HGVS protein notation (empty if not coding) |
| `base` | str | Base observed on this read at `pos` |
| `barcode` | str | 2 µm spatial barcode of the read (`CB` tag), e.g. `s_002um_00334_02282-1` |
| `x` | int | First index of the barcode = `array_row` of the 2 µm grid (0–3349) |
| `y` | int | Second index of the barcode = `array_col` of the 2 µm grid (0–3349) |
| `umi` | str | UMI of the read (`UB` tag) |


## 5. Post-processing: somatic filter and clustering score

The last rule of the pipeline (`post_process`, script
`workflow/scripts/post_process.py`) reads the 24 tables above and writes two
files.

**`data/pancreas.somatic.parquet`**: the somatic-like variants. No matched
normal tissue is available, so two filters are used as a proxy: the variant is
absent from dbSNP common (`rsid` is `.`) and its depth is above
`somatic_min_dp` (`dp > 40`). The columns are the same as in section 4.

**`data/pancreas.gyration.parquet`**: one row per somatic-like variant, sorted
from the most to the least spatially clustered. For a variant *V* carried by
*N* distinct spots at positions *c<sub>i</sub>* with centroid *c̄*, the radius
of gyration is

    R(V) = sqrt( 1/N · Σ ‖c_i − c̄‖² )

and the score is `s(V) = R(V) / R_all`, where `R_all` is the radius of gyration
of all spots carrying any variant (before filtering). `s` is near 1 for a
variant spread like the section and tends to 0 for a concentrated one.

| Column | Type | Description |
|--------|------|-------------|
| `chr`, `pos`, `ref`, `alt` | | Variant, as in section 4 |
| `gene`, `effect`, `impact`, `hgvs_c`, `hgvs_p` | str | snpEff annotation, as in section 4 |
| `dp`, `ao`, `af` | | Pseudo-bulk depth, alternate count and allele fraction |
| `n_spots` | int | Number of distinct 2 µm spots carrying the variant (*N*) |
| `x_mean`, `y_mean` | float | Centroid of these spots, in 2 µm grid units |
| `radius` | float | Radius of gyration *R(V)*, in 2 µm grid units |
| `n_molecules` | int | Number of distinct molecules (barcode + UMI) carrying the variant |
| `s` | float | Clustering score `R(V) / R_all` |

A variant seen in a single spot has `s = 0`, so filter on `n_spots` before
ranking (e.g. `n_spots >= 20`).

The script can also be run on its own:

```bash
python workflow/scripts/post_process.py data/pancreas.chr*.spatial_variant.parquet \
    --somatic data/pancreas.somatic.parquet \
    --score data/pancreas.gyration.parquet \
    --min-depth 40
```
