# rhapsodist

Snakemake workflow for BD Rhapsody WTA single-cell RNA-seq data on v1, Enhanced and Enhanced V2 beads.

It standardizes the cell barcode and UMI of read 1, builds a whitelist of the observed cell barcodes, and counts with any of STARsolo, kallisto/bustools, salmon/alevin and the BD Rhapsody pipeline. Each aligner writes an HDF5-backed SingleCellExperiment, and optionally an h5ad file. Sample tags are demultiplexed, and HTML reports describe each aligner and compare them.

- [Install](install.md)
- [Quick start](quickstart.md)
- [Configuration](configuration.md)
- [Outputs](outputs.md)
- [Command line](cli.md)

Source: [github.com/imallona/rhapsodist](https://github.com/imallona/rhapsodist)
