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
    import xarray as xr
    import hvplot.xarray
    import hvplot.polars 
    import marimo as mo
    import numpy as np 
    return Image, alt, mo, mpimg, np, patches, pl, plt, sns, xr


@app.cell
def _(pl):
    # cosmic_df = (
    #     pl.read_csv(
    #         "../../hg38/cosmic/Cosmic_MutantCensus_v102_GRCh38.tsv.gz",
    #         separator="\t",
    #         infer_schema_length=100000000,
    #     )
    #     .with_columns(
    #         (
    #             pl.lit("chr")
    #             + pl.concat_str(
    #                 [
    #                     "CHROMOSOME",
    #                     "GENOME_START",
    #                     "GENOMIC_WT_ALLELE",
    #                     "GENOMIC_MUT_ALLELE",
    #                 ],
    #                 separator="-",
    #             )
    #         ).alias("id")
    #     )
    #     .select(
    #         "id",
    #         "GENE_SYMBOL",
    #         "GENOMIC_MUTATION_ID",
    #         "MUTATION_ID",
    #         "MUTATION_DESCRIPTION",
    #         "HGVSC",
    #         "MUTATION_SOMATIC_STATUS",
    #     )
    #     .unique()
    # )

    #cosmic_df.write_parquet("Cosmic_MutantCensus_v102_GRCh38.parquet")

    cosmic_df = pl.read_parquet("/Linux_data/recherche/codes/space_variant/Cosmic_MutantCensus_v102_GRCh38.parquet")
    cosmic_df
    return (cosmic_df,)


@app.cell
def _(mo):
    mo.md(r"""
    # Preprocessing variants
    """)
    return


@app.cell
def _(cosmic_df, pl):
    df = pl.read_parquet("/Linux_data/recherche/codes/space_variant/data/pancreas.*.parquet")
    df = df.insert_column(0, pl.concat_str(["chr", "pos", "ref", "alt"], separator="-").alias("id"))

    # Récupération des cordonnées image
    pos_df = pl.read_parquet("/Linux_data/recherche/dataset/pancreas/my_output/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet")
    df = df.join(pos_df, on="barcode")
    df = df.rename({"pxl_row_in_fullres": "screen_x", "pxl_col_in_fullres": "screen_y"})

    # Add cosmic annotation 
    df = df.join(
        cosmic_df.filter(
            pl.col("MUTATION_SOMATIC_STATUS") == "Confirmed somatic variant"
        ).select("id", pl.col("GENOMIC_MUTATION_ID").alias("cosmic")),
        on="id",
        how="left"
    )


    df
    return df, pos_df


@app.cell
def _(df, pl):
    somatic_df = df.filter(~pl.col("rsid").str.starts_with("rs")).filter(pl.col("af").is_between(0.1, 0.8))
    somatic_df
    return (somatic_df,)


@app.cell
def _(df, somatic_df):
    all_count = len(df["id"].unique())
    somatic_count = len(somatic_df["id"].unique())
    cosmic_count = len(somatic_df["cosmic"].unique().is_not_null())

    return all_count, cosmic_count, somatic_count


@app.cell
def _(cosmic_df):
    del cosmic_df
    return


@app.cell
def _(all_count, cosmic_count, mo, somatic_count):
    mo.md(rf"""
    - Nombre de variants au total {all_count}
    - Nombre de variants sans rsid {somatic_count}
    - Nombre de variants cosmic {cosmic_count}
    """)
    return


@app.cell
def _(alt, df):
    chart_a = df["effect"].value_counts().sort("count", descending=True).head(20).plot.bar(y=alt.Y("effect", sort="-x"), x="count")
    chart_b = df["biotype"].value_counts().sort("count", descending=True).head(20).plot.bar(y=alt.Y("biotype", sort="-x"), x="count")
    chart_c = df["impact"].value_counts().sort("count", descending=True).head(20).plot.bar(y=alt.Y("impact", sort="-x"), x="count")

    chart_a|chart_b|chart_c
    return


@app.cell
def _(df, pl, q, sns, somatic_df):
    dist_a = (
        df.filter(pl.col("dp") > 40)
        .filter(pl.col("biotype") == "protein_coding")
        .filter(pl.col("af") < 0.99)
        .filter(~pl.col("gene").str.contains("HLA"))
        .select(["id", "af"])
        .unique().with_columns(pl.lit("raw").alias("filtered"))
    )

    dist_b = (
        somatic_df.filter(pl.col("dp") > 40)
        .filter(pl.col("biotype") == "protein_coding")
        .filter(pl.col("af") < 0.99)
        .filter(~pl.col("gene").str.contains("HLA"))
        .select(["id", "af"])
        .unique().with_columns(pl.lit("dbsnp").alias("filtered"))
    )
    q
    dist = pl.concat([dist_a, dist_b])

    g = sns.displot(dist, x="af", col="filtered", kde=True, bins = 90)
    g.fig.suptitle("Distribution of allele frequencies (AF) ", y = 1.05)
    for ax in g.axes.flat:
        for line in ax.lines:
            line.set_color("red")
    g
    return


@app.cell
def _():
    return


@app.cell
def _(mo):
    af_slider = mo.ui.range_slider(0, 1, step=0.01, value=[0,1], show_value=True, label="frequence allelique")
    af_slider
    return (af_slider,)


@app.cell
def _(af_slider, alt, pl, somatic_df):
    somatic_df["id", "af", "gene", "biotype"].unique().filter(pl.col("af").is_between(af_slider.value[0], af_slider.value[1])).filter(pl.col("biotype") == "protein_coding")[
        "gene"
    ].value_counts().sort("count", descending=True).head(50).plot.bar(
        y=alt.Y("gene", sort="-x"), x="count"
    ).interactive(False)
    return


@app.cell
def _(Image, np, xr):
    img = Image.open("/Linux_data/recherche/codes/space_variant/tissue_hires_image.png").convert("RGBA")
    img = np.rot90(img, k=3)
    img = np.fliplr(img)

    arr = np.array(img)
    coords = {
        "screen_y" : range(arr.shape[0]),
        "screen_x": range(arr.shape[1]),
        "channel": ["R","G","B","A"]
    }

    d = xr.DataArray(arr, coords=coords, dims=("screen_y","screen_x","channel"), name="image2")

    #d.hvplot.image(x="x",  y="y")

    return (d,)


@app.cell
def _(d, somatic_df):
    w = 500
    h = 500
    image_plot = d.hvplot.rgb(bands="channel", x="screen_x", y="screen_y", aspect="equal", rasterize=True, width=w, height=h)
    image_plot


    space_df = somatic_df.select(["id", "screen_x", "screen_y","x","y", "af", "dp"])
    variant_plot = space_df.hvplot.scatter(x="screen_x", y="screen_y",rasterize=True, aspect='equal', aggregator="count", datashade=True,width=w, height=h)

    overlay = image_plot + variant_plot
    overlay = overlay.opts( shared_axes=True)
    overlay
    return (space_df,)


@app.cell
def _():
    import umap
    return (umap,)


@app.cell
def _(pl, space_df):
    from sklearn.decomposition import PCA

    factor = 50
    down_df = space_df.with_columns((pl.col("x") // factor).alias("X"),
                          (pl.col("y") // factor).alias("Y") )

    down_df=down_df.with_columns(pl.concat_str([pl.lit("xy"),"X","Y"], separator="-").alias("space_id"))
    down_df = down_df.group_by(["space_id", "id"]).count()
    down_df = down_df.pivot(on="id", index="space_id", values="count").fill_null(0)


    down_df = down_df.to_pandas().set_index("space_id")
    down_df

    return (down_df,)


@app.cell
def _(down_df, pl, umap):

    reducer = umap.UMAP(
        n_neighbors=100,      # nombre de voisins locaux
        n_components=2,      # dimensions finales (2D pour visualiser)
        metric='euclidean',  # distance utilisée
        random_state=42
    )

    embedding = reducer.fit_transform(down_df)

    # 4. Mettre le résultat dans un DataFrame
    embedding_df = pl.DataFrame(embedding)
    embedding_df.columns = ["P1", "P2"]
    embedding_df
    return (embedding_df,)


@app.cell
def _(embedding_df):
    embedding_df.plot.scatter(x="P1", y="P2")
    return


@app.cell
def _(pdf):
    pdf
    return


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
def _(fdf):
    fdf
    return


@app.cell
def _():
    return


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
    #super_df = super_df.insert_column(0, pl.concat_str(["chrom","pos","ref","alt"],separator="-").alias("id"))
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
def _(pl, super_df):
    super_df.filter(pl.col("rsid").is_null()).filter(pl.col("af") > 0.6).hvplot.scatter(x="screen_x", y="screen_y", s=1, rasterize=True, aspect='equal',dynspread=True, aggregator="count", datashade=True, width=800, height=800)
    return


@app.cell
def _(super_df):
    super_df.select("id", "af", "barcode").unique().pivot(on="barcode", index="id", values="af") 
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
