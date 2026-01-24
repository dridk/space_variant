.PHONY: all

run:
	snakemake -fps workflow/Snakefile --configfile workflow/config.yml -d data --cores 20 --rerun-incomplete 
