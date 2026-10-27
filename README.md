# Rhapsodist + SimpleAF

Prepare BD Rhapsody Enhanced V2/V3 whole-transcriptome reads and quantify spliced, unspliced and ambiguous counts for RNA-velocity workflows. The pipeline ends with validated count matrices; it does not estimate velocity, annotate cells or merge the original BD expression matrix.

[![Velocity checks](https://github.com/Thokas99/SimpleRhapsody/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/Thokas99/SimpleRhapsody/actions/workflows/tests.yml) [![Open issues](https://img.shields.io/github/issues/Thokas99/SimpleRhapsody)](https://github.com/Thokas99/SimpleRhapsody/issues) ![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)

## Workflow

The diagram shows the path from paired FASTQs and a prebuilt spliceu reference to the validated H5AD. Downsampling is optional; lanes otherwise stay separate through preparation.

```mermaid
flowchart TD
    ref["Config + prebuilt spliceu reference"] --> preflight["Reference preflight<br/>piscem index + transcript-to-gene map"]
    reads["Paired R1/R2 FASTQ lanes"] --> downsample{"Downsample below 100%?"}
    downsample -->|yes| sampled["Concatenate lanes, then sample pairs"]
    downsample -->|no| lanes["Keep lanes as separate units"]
    preflight --> chemistry["Enhanced V2/V3 linker check"]
    sampled --> chemistry
    lanes --> chemistry
    chemistry --> prepare["Cutadapt + streamed barcode/UMI crop"]
    prepare --> assemble["Assemble prepared gzip streams in lane order"]
    assemble --> whitelist["Observed whitelist + QC receipt"]
    preflight --> quant["SimpleAF USA quantification<br/>both orientations; cr-like UMI resolution"]
    whitelist --> quant
    quant --> native["Native USA matrices: S, U, A"]
    native --> export["BD barcode conversion + H5AD<br/>X / spliced = S + A; unspliced = U"]
    export --> validate["Compare H5AD labels and counts with native matrices"]
    validate --> output["Validated H5AD + QC tables"]
    preflight -->|failure| stop["Stop before FASTQ preparation"]
```

## Count and barcode contract

SimpleAF uses geometry `1{b[27]u[8]}2{r:}`, both read orientations and `cr-like` UMI resolution. The config value `resolution: usa` selects the USA output contract; it is not the SimpleAF UMI-resolution flag.

| Output | Contents |
|---|---|
| Native SimpleAF | Separate spliced (S), unspliced (U) and ambiguous (A) matrices |
| H5AD `X`, `layers['spliced']` | S + A |
| H5AD `layers['unspliced']` | U |
| H5AD `obs_names` | Decimal BD cell index by default |
| H5AD `obs` | Original 27-nt barcode, integer BD cell index and sample |
| H5AD `var` | Original Ensembl IDs/order; reference gene symbols when available |

Quantification uses the observed barcode whitelist with `--unfiltered-pl --min-reads 1`. The whitelist contains prepared-R1 barcodes valid across the ordered BD 384×3 panels; export applies no further barcode correction or cell/gene filtering. It retains zero-count rows/columns and fractional estimates. BD cell identifiers are sample-local: join captures using both `sample` and the cell identifier. An optional `bd_reference` list/TSV checks names but does not import BD expression counts. See [barcode conversion](workflow/scripts/bd_barcode_conversion.py).

## Configure and run

### 1. Configure a sample

Copy the supplied configuration and replace the WTA1 input, output and reference paths. Use absolute paths visible inside the container:

```bash
cp configs/velocity.yaml configs/my_sample.yaml
```

Set these fields in `configs/my_sample.yaml`:

- `working_dir`: output directory available inside the container.
- `simpleaf.index`: prebuilt spliceu reference; the workflow does not build or modify it.
- `simpleaf.resolution: usa`.
- `bead_chemistry`: `enhanced_v2` or `enhanced_v3`.
- `samples[].uses.cb_umi_fq` and `cdna_fq`: R1 and R2 file lists with lanes in the same order.
- `samples[].uses.allowedlist: 384` and `diversity_insets: 'yes'`.
- `samples[].name`: unique; start with a letter or number and use only letters, numbers, `.`, `_` or `-`.
- `downsample`: percentage greater than 0 and at most 100. The default, `100`, keeps all read pairs.

Reference preflight checks the piscem index files and `t2g_3col.tsv` before FASTQ preparation. An identity manifest is optional; set `simpleaf.reference_manifest` to provide one. Without it, genome and release are recorded as unknown. See [reference validation](workflow/scripts/simpleaf_quant.py).

### 2. Build a new image

On Linux x86-64, with Apptainer, fakeroot and Pixi available, run:

```bash
./container_build.sh "$PWD/rhapsodist-simpleaf.performance.sif"
```

The script uses `pixi.lock`, runs development checks and writes a checksum. It refuses to overwrite an existing image. The installed cluster SIF predates the FASTQ-QC changes, so build under a new filename.

### 3. Dry-run, then execute

Copy the new image and config to the cluster paths below (or bind their actual locations). From a Slurm allocation with at least two CPU cores, review the DAG:

```bash
apptainer exec --cleanenv --bind /globalscratch \
  /globalscratch/tsirchi/cemlab/simpleaf/containers/rhapsodist-simpleaf.performance.sif \
  rhapsodist --configfile /globalscratch/tsirchi/cemlab/simpleaf/configs/my_sample.yaml \
  --cores "${SLURM_CPUS_PER_TASK:-16}" --resources mem_mb=64000 --dry-run
```

Remove `--dry-run` to execute. Use a fresh `working_dir` after changing workflow code or configuration; Snakemake may invalidate downstream results. Rules run inside the single container. The supplied [Slurm script](slurm/velocity.sbatch) still names the older image and must be updated before submission.

## Outputs and quality checks

`working_dir` contains persistent prepared FASTQs, the observed whitelist, logs, native USA output, one `{sample}_velocity.h5ad`, and QC tables. Prepared lane units are temporary and removed after successful ordered assembly; completed prepared reads remain so quantification can be retried. Use Snakemake's `--rerun-incomplete` for interrupted outputs.

Routine QC validates prepared R1 structure and barcode/UMI length, checks Cutadapt pair counts, and records a receipt bound to sample, input paths, sizes and modification times. The default does not reread all prepared R2 records to independently check gzip integrity. Set `fastq_qc_exhaustive: true` for a full mate/structure/gzip scan; it can take hours. Linker classification samples a bounded number of reads, not every read. File identity metadata is not a content checksum.

The H5AD is checked against native USA cell/gene order and S+A/U values before publication. See the [performance audit](PERFORMANCE_AUDIT.md) for QC limits, recovery behavior, resources and recorded fixture results. Fixture and synthetic checks are not production-sample or Slurm performance evidence.

## Reproducibility

The run records workflow, config, lockfile, reference and QC-receipt provenance in native output and H5AD metadata. Small reference files are hashed; the large index and FASTQs are not content-hashed. See [pixi.toml](pixi.toml), [pixi.lock](pixi.lock), and the [container definition](rhapsodist.def) for pinned software and image inputs.

Preparation follows [imallona/rhapsodist at commit 7991f870fdfe38e12648673d330cbfd4839e755d](https://github.com/imallona/rhapsodist/tree/7991f870fdfe38e12648673d330cbfd4839e755d), under GPL-3.0; see [LICENSE.txt](LICENSE.txt) and [AUTHORS.txt](AUTHORS.txt). Documentation guidance: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065).
