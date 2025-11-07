.PHONY: all

run:
	snakemake -ps workflow/Snakefile --configfile workflow/config.yml -d data output.spatial_variant.parquet --cores 20 --rerun-incomplete 
