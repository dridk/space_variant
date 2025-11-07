import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import polars as pl 
    import seaborn as sns
    import altair as alt
    from PIL import Image
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.patches as patches
    import hvplot.polars 
    import marimo as mo
    import numpy as np 
    return alt, mo, mpimg, np, patches, pl, plt, sns


@app.cell
def _(mo):
    mo.md(r"""
    # Preprocessing variants
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Récupération de dbsnp
    """)
    return


@app.cell
def _(pl):
    snp_df = (
        pl.scan_csv(
            "/Linux_data/recherche/hg38/dbSnp155Common.bed",
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "start",
                "end",
                "rsid",
                "ref",
                "alt_count",
                "alts",
                "shift",
                "freq_count",
                "minor_allele_freq",
                "major",
                "minor",
                "max_impact",
                "class",
                "note",
                "offset",
                "len"
            ]).select("chrom", pl.col("start").cast(pl.UInt32), "rsid", "ref" ).collect()
    ).insert_column(1, (pl.col("start") + 1).alias("pos"))
    snp_df
    return (snp_df,)


@app.cell
def _(pl, snp_df):
    df = pl.read_parquet("spatial_variants.parquet")

    df  = df.with_columns(pl.col("af").cast(pl.Float32), pl.col("ad").str.split(","))
    # Ajout d'un ID 
    df = df.insert_column(0, pl.concat_str(["chrom", "pos", "ref", "alt"], separator="-").alias("id"))


    df = df.join(snp_df, on=["chrom", "pos"], how="left").select(pl.exclude("start", "ref_right"))
    return (df,)


@app.cell
def _(mo):
    mo.md("""
    ## Frequence allelique
    """)
    return


@app.cell
def _(df, sns):
    fdf = df.select("chrom", "pos", "ref", "alt", "dp", "ad","af", "rsid").unique()
    sns.displot(fdf["af"], bins = 100, kde=True)
    return (fdf,)


@app.cell
def _(fdf, pl, sns):
    sns.displot(fdf.filter(pl.col("rsid").is_null())["af"], bins = 100, kde=True)
    return


@app.cell
def _(pl):
    cosmic_df = (
        pl.read_csv(
            "../../hg38/cosmic/Cosmic_MutantCensus_v102_GRCh38.tsv.gz",
            separator="\t",
            infer_schema_length=100000000,
        )
        .with_columns(
            (
                pl.lit("chr")
                + pl.concat_str(
                    [
                        "CHROMOSOME",
                        "GENOME_START",
                        "GENOMIC_WT_ALLELE",
                        "GENOMIC_MUT_ALLELE",
                    ],
                    separator="-",
                )
            ).alias("id")
        )
        .select(
            "id",
            "GENE_SYMBOL",
            "GENOMIC_MUTATION_ID",
            "MUTATION_ID",
            "MUTATION_DESCRIPTION",
            "HGVSC",
            "MUTATION_SOMATIC_STATUS",
        )
        .unique()
    )

    cosmic_df.write_parquet("Cosmic_MutantCensus_v102_GRCh38.parquet")
    return (cosmic_df,)


@app.cell
def _(cosmic_df):
    cosmic_df
    return


@app.cell
def _(pl):
    pos_df = pl.read_parquet("/Linux_data/recherche/dataset/pancreas/my_output/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet")

    pos_df
    return (pos_df,)


@app.cell
def _(df, pl, pos_df):
    super_df = df.join(pos_df, on="barcode")
    super_df = super_df.insert_column(0, pl.concat_str(["chrom","pos","ref","alt"],separator="-").alias("id"))
    super_df = super_df.rename({"pxl_row_in_fullres": "screen_x", "pxl_col_in_fullres": "screen_y"})
    #Rayon et cercle
    super_df = super_df.with_columns([
        pl.col("screen_x").mean().over("id").alias("cx"),
        pl.col("screen_y").mean().over("id").alias("cy"),

    ]).with_columns([
       ((pl.col("cx")-pl.col("screen_x"))**2 + (pl.col("cy")-pl.col("screen_y"))**2).sqrt().max().over("id").alias("radius")


    ]).with_columns(pl.col("id").count().over("id").alias("count"))


    super_df
    return (super_df,)


@app.cell
def _(cosmic_df, super_df):
    all_df = super_df.join(cosmic_df, on="id").unique()
    all_df
    return (all_df,)


@app.cell
def _(all_df):
    all_df.select(["GENE_SYMBOL", "id", "GENOMIC_MUTATION_ID", "MUTATION_DESCRIPTION", "HGVSC", "af", "dp", "MUTATION_SOMATIC_STATUS"]).unique()
    return


@app.cell
def _(all_df, alt):
    all_df["GENE_SYMBOL"].value_counts().sort("count", descending=True).plot.bar(alt.X("GENE_SYMBOL", sort="-y"), "count")
    return


@app.cell
def _(sns, super_df):
    import datashader as ds
    import pandas as pd 
    import colorcet


    sns.scatterplot(super_df, x="screen_x", y="screen_y", size=1, alpha=0.1)
    return


@app.cell
def _(super_df):
    super_df.hvplot.scatter(x="screen_x", y="screen_y", s=1, rasterize=True, cnorm="eq_hist", cmap="Reds", aspect='equal',dynspread=True, aggregator="count", datashade=True)
    return


app._unparsable_cell(
    r"""
    import holoviews as hv 

    hv_img = hv.Image(data, params




    """,
    name="_"
)


@app.cell
def _():
    # p = "/Linux_data/recherche/dataset/pancreas/my_output/outs/binned_outputs/square_002um/spatial/tissue_hires_image.png"

    # img = mpimg.imread(p)
    # fig, ax = plt.subplots(figsize=(10, 10))
    # print(img.shape)
    # img = np.rot90(img, k=3)
    # img = np.fliplr(img)
    # ax.imshow(img)

    # x = super_df["screen_x"]
    # y = super_df["screen_y"]

    # print(x.max())
    # print(y.max())
    # ax.scatter(x, y,  s=1, marker='x',  color="green", linewidths=1,  alpha=s.value)

    # ax
    return


@app.cell
def _():
    def adn_vers_couleur(sequence):
        """
        Convertit une séquence d'ADN de 9 bases en couleur hexadécimale.

        Encodage: A=00, C=01, G=10, T=11 (2 bits par base)
        9 bases × 2 bits = 18 bits
        On complète avec des zéros pour obtenir 24 bits (format RGB)

        Args:
            sequence: String de 9 caractères contenant A, C, G, ou T

        Returns:
            Couleur au format hexadécimal (#RRGGBB)
        """
        if len(sequence) != 9:
            raise ValueError("La séquence doit contenir exactement 9 bases")

        # Mapping des bases vers leurs valeurs binaires
        base_to_bits = {
            'A': 0b00,  # 0
            'C': 0b01,  # 1
            'G': 0b10,  # 2
            'T': 0b11   # 3
        }

        # Conversion de la séquence en valeur binaire
        valeur_binaire = 0
        for base in sequence.upper():
            if base not in base_to_bits:
                raise ValueError(f"Base invalide: {base}. Utilisez A, C, G ou T")
            valeur_binaire = (valeur_binaire << 2) | base_to_bits[base]

        # On a maintenant 18 bits, on décale de 6 bits pour obtenir 24 bits
        valeur_24bits = valeur_binaire << 6

        # Extraction des composantes RGB (8 bits chacune)
        rouge = (valeur_24bits >> 16) & 0xFF
        vert = (valeur_24bits >> 8) & 0xFF
        bleu = valeur_24bits & 0xFF

        # Conversion en format hexadécimal
        couleur_hex = f"#{rouge:02X}{vert:02X}{bleu:02X}"

        return couleur_hex


    adn_vers_couleur("AAAAATACA")
    return (adn_vers_couleur,)


@app.cell
def _(super_df):
    super_df
    return


@app.cell
def _(mo):
    df_ui = mo.ui.range_slider(
        start=0, stop=3304, step=0.05, value=[10, 100], full_width=True
    )
    df_ui
    return (df_ui,)


@app.cell
def _(df_ui):
    print(df_ui.value)
    return


@app.cell
def _(all_df, mo):
    variant = mo.ui.dropdown(
    options=all_df["id"].unique(), searchable=True

    )

    variant
    return (variant,)


@app.cell
def _(adn_vers_couleur, mpimg, np, patches, pl, plt, super_df, variant):
    path = "/Linux_data/recherche/dataset/pancreas/my_output/outs/binned_outputs/square_002um/spatial/tissue_hires_image.png"

    rdf = super_df.filter(pl.col("id") == variant.value)

    x = rdf["screen_x"].to_list()
    y = rdf["screen_y"].to_list()

    count = len(rdf)
    circle = rdf.select(["cx","cy", "radius"]).row(0)
    colors = rdf.select(pl.col("umi").map_elements(adn_vers_couleur, return_dtype=pl.String))["umi"].to_list()



    img = mpimg.imread(path)
    img = np.rot90(img, k=3)
    img = np.fliplr(img)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img)

    hauteur, largeur = img.shape[:2]

    # Drole de rotation ? 
    py = y
    px = x

    for i in range(len(px)):
        ax.plot(px[i], py[i], 'x', markersize=5, markeredgewidth=2,  color=colors[i], label='Points marqués')

    ax.scatter(px, py, c=[(0,1,0, 0.01)], s=1, marker='x', linewidths=1)


    p = patches.Circle((circle[0], circle[1]), circle[2], fill=False, color='red', linewidth=2)
    ax.add_patch(p)
    plt.tight_layout()
    ax.set_title(f"{variant.value} avec {count} evenements")
    ax
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
