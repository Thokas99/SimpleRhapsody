# Command line

<!-- rhapsodist-help -->

`rhapsodist` runs `snakemake --snakefile workflow/Snakefile --use-conda --benchmark-extended` with the given config and cores, and prints the command to standard error. Other arguments go to Snakemake unchanged:

```
rhapsodist --configfile my_config.yaml --cores 10 --rerun-incomplete
rhapsodist --configfile my_config.yaml --config "aligner=['starsolo']" working_dir=output/starsolo_only
```

The same run without the wrapper:

```
snakemake --snakefile workflow/Snakefile --use-conda --benchmark-extended --configfile my_config.yaml --cores 10
```
