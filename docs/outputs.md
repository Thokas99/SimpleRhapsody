# Outputs

Everything is written under `working_dir`:

```
{aligner}/{sample}/{sample}_{aligner}_sce_hdf5/   count object of each aligner
{aligner}/{sample}/{sample}_{aligner}.h5ad        with output_format: h5ad or both
{aligner}/{sample}/by_sampletag/{tag}/            counts of the singlets of each sample tag
{aligner}/{sample}/descriptive_report.html
{sample}_comparison.html                          with two or more aligners
{sample}_biology.html                             with run_biology_report: true
sampletags/{sample}/sampletag_report.html
linker_qc/{sample}_linker_qc.html
benchmarks_report.html                            time and memory per rule
benchmarks/                                       one Snakemake benchmark file per rule
run_info.tsv                                      CPU, memory and Snakemake settings of the run
logs/
data/                                             standardized reads, references and indices
```

`{aligner}` is one of `starsolo`, `kallisto`, `alevin` and `sbg`.

Load a count object in R:

```r
sce <- HDF5Array::loadHDF5SummarizedExperiment("starsolo/my_sample/my_sample_starsolo_sce_hdf5")
```

and an h5ad file in Python:

```python
import anndata
adata = anndata.read_h5ad("starsolo/my_sample/my_sample_starsolo.h5ad")
```
