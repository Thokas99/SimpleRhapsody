# Rhapsodist + SimpleAF

Prepare BD Rhapsody Enhanced V2/V3 WTA reads and quantify spliced, unspliced and ambiguous gene counts for downstream RNA-velocity analysis. The workflow preserves Rhapsodist's barcode/linker preparation and uses a prebuilt spliceu reference for SimpleAF quantification. It produces native USA counts and one H5AD with BD Rhapsody 4.0-compatible cell identifiers.

The workflow ends at validated count matrices. Velocity estimation, cell annotation, lineage inference and merging with the original BD expression matrix are outside its scope.

## Workflow and scientific contract

```text
Prebuilt reference → reference preflight before FASTQ processing
Paired FASTQs → optional sample-wide concatenation/downsampling → linker/chemistry gate
             → streamed preparation units → ordered prepared gzip assembly
             → observed whitelist + QC receipt → SimpleAF USA → H5AD → native-count validation
```

Preparation retains a 27-nucleotide cell barcode and an 8-nucleotide UMI. The observed whitelist contains only barcodes actually found in prepared R1 whose three segments occur in the ordered BD 384-entry panels. Ambiguous UMIs are retained; no additional barcode correction or expression filtering is introduced during export.

Quantification uses `1{b[27]u[8]}2{r:}`, both orientations, `cr-like` UMI resolution, and the observed whitelist through `--unfiltered-pl --min-reads 1`. Here **USA** means spliced (S), unspliced (U) and ambiguous (A) count blocks; the config's `resolution: usa` specifies this output contract, not the SimpleAF UMI-resolution flag.

| H5AD field | Contents |
|---|---|
| `X` and `layers['spliced']` | S + A |
| `layers['unspliced']` | U |
| `obs_names` | Decimal BD cell indices by default |
| `obs` | Original `barcode_27nt`, integer `bd_cell_index`, and `sample` |
| `var` | Original Ensembl gene identifiers/order and reference gene symbols when supplied |
| `uns['pipeline']` | Quantification/conversion provenance, workflow/config/lock hashes and preparation-receipt hash |

Native USA files retain separate S/U/A counts. Export preserves quantified cells, genes, zero-count rows/columns and fractional estimates without rounding. BD identifiers are sample-local: use both `sample` and the cell identifier when joining captures. An optional explicit `bd_reference` identifier list/TSV can supply full cell names; it does not import BD expression counts. See the [conversion implementation](workflow/scripts/bd_barcode_conversion.py) for its accepted formats.

## Build a new container

The runtime targets Linux x86-64 and requires Apptainer with fakeroot for this build recipe. From the source checkout, with Pixi available:

```bash
./container_build.sh "$PWD/rhapsodist-simpleaf.performance.sif"
```

The [build script](container_build.sh) installs the locked environment, runs development checks, builds and tests the image, then writes its checksum. It refuses to overwrite the destination image. The existing SIF predates the FASTQ-QC changes; a source update does not update that image. Keep the installed cluster image intact while validating a newly named artifact.

The container includes Snakemake, Cutadapt, pigz, seqtk, SimpleAF, piscem, alevin-fry and the Python export dependencies. Exact packages/artifact hashes are recorded in [pixi.lock](pixi.lock); primary pins are in [pixi.toml](pixi.toml). The [definition](rhapsodist.def) pins the Debian base image and Pixi installer. These pins support reproducible environments; build timestamps and bootstrap packages can still prevent byte-identical SIFs.

## Configure a sample

Create `configs/my_sample.yaml` using absolute paths accessible inside the container. This example describes paths to replace with your own inputs:

```yaml
working_dir: /globalscratch/tsirchi/cemlab/simpleaf/results_performance
nthreads: 16
bead_chemistry: enhanced_v2  # enhanced_v3 uses the same GTGA/GACA preparation
cb_umi_max_errors: 0
downsample: 100             # percentage; 100 retains all input pairs

export:
  python: python3
simpleaf:
  index: /globalscratch/tsirchi/cemlab/simpleaf/references/GRCh38_gencode_v49_spliceu/output/index
  resolution: usa

samples:
  - name: sample01
    uses:
      cb_umi_fq: [/globalscratch/tsirchi/cemlab/simpleaf/input/sample01_R1.fastq.gz]
      cdna_fq: [/globalscratch/tsirchi/cemlab/simpleaf/input/sample01_R2.fastq.gz]
      allowedlist: 384
      diversity_insets: 'yes'
```

For multiple lanes, extend both lists in the same lane order. Samples require unique, safe path-component names. The supplied [WTA1 configuration](configs/velocity.yaml) describes seven lanes and the GRCh38/GENCODE v49 reference; inspect it before use. The workflow accepts an index directory or its parent, requires a piscem index and an S/U transcript-to-gene map, and does not build or modify the production reference.

Reference preflight validates the index/map before FASTQ processing. Its ordinary dependencies ensure early failure; a reference update can also invalidate preparation, so inspect a dry run before changing a completed workflow. Genome/release labels come from an optional identity manifest; they are recorded as `unknown` when it is absent. The driver no longer substitutes GRCh38/v49 for an arbitrary reference.

The workflow automatically reads `reference_manifest.json` inside the resolved index directory when present. To use a manifest elsewhere, add `reference_manifest: /absolute/path/reference_manifest.json` under `simpleaf`. Its format is:

```json
{
  "version": 1,
  "index_type": "spliceu",
  "genome": "<assembly documented by the reference build>",
  "GENCODE_release": "<release documented by the reference build>",
  "reference_t2g_sha256": "<SHA256 of t2g_3col.tsv>",
  "reference_index_meta_sha256": "<SHA256 of simpleaf_index.json>"
}
```

Replace the placeholders with build-record identities and actual hashes. A supplied manifest must be complete and match both small reference files; otherwise preflight fails before preparation. Matching hashes bind the declared identity to those files, but do not authenticate the biological annotation or checksum the entire large index. The preflight records index-file paths/sizes/mtimes and rechecks them before quantification; see [reference validation](workflow/scripts/simpleaf_quant.py).

## Inspect the DAG before execution

Place the newly tested image and your config at the paths below. From a Slurm allocation, use at least two CPU cores and review a dry run:

```bash
apptainer exec --cleanenv --bind /globalscratch \
  /globalscratch/tsirchi/cemlab/simpleaf/containers/rhapsodist-simpleaf.performance.sif \
  rhapsodist --configfile /globalscratch/tsirchi/cemlab/simpleaf/configs/my_sample.yaml \
  --cores "${SLURM_CPUS_PER_TASK:-16}" --resources mem_mb=64000 --dry-run
```

Remove `--dry-run` when ready to process the configured inputs. Use a fresh `working_dir` when evaluating changed source or migrating completed runs: Snakemake can invalidate downstream results after code/config changes. The supplied [Slurm script](slurm/velocity.sbatch) uses the older image filename and must be adjusted before launching a new artifact. Snakemake rules invoke local tools inside the single container.

## FASTQ quality control

Routine QC reuses exact paired Cutadapt statistics and validates final R1 during the existing whitelist traversal. Each final R1 record must have valid FASTQ structure, matching sequence/quality lengths and exactly 27 barcode + 8 UMI bases. Its total must match Cutadapt's output-pair count. Each prepared R2 unit is written directly by Cutadapt while cropped R1 streams through the unchanged line-wise `cut` operation into pigz. Prepared units are assembled in configured order without decompression; both operations preserve record order and count. Routine prepared-pair consistency therefore relies on this construction, without an independent R2 decompression pass.

For unsampled inputs, each lane runs through Cutadapt independently and supplies its exact raw pair count. The legacy first-token mate-name comparison (including `/1` and `/2` handling) runs inside the existing Cutadapt/dnaio iterator, before filtering, so it requires no second raw traversal. This adapter is guarded for Cutadapt 4.9 and Linux fork workers, matching the pinned runtime. The classifier scans the same bounded R1 sequence prefix across lanes in their configured order; incompatible or unknown chemistry now fails before Cutadapt. Downsampled inputs retain the original sample-wide concat/seqtk selection and independent original-lane count/pairing scans. Retained percentage always uses **all original input pairs** as its denominator. The receipt binds counts to the configured sample, ordered file paths, sizes and mtimes; a content hash protects the observed whitelist. The downstream QC gate checks that receipt without rereading FASTQs.

| Mode | Additional checks | Limit |
|---|---|---|
| Routine, default | Receipt/file association, whitelist validity and exact counts established during preparation | Does not independently reread prepared R2 or verify its post-generation gzip integrity |
| `fastq_qc_exhaustive: true` | Full original/prepared mate traversal, structure, gzip integrity, counts and transformed-header comparison | Can take hours; opt-in for troubleshooting |

The enhanced-linker requirement gates preparation and is checked again before quantification. Linker classification uses a bounded sample of reads, not exhaustive chemistry verification. File size/mtime checks are not cryptographic FASTQ checksums, and configured sample identity cannot authenticate a biologically mislabeled input. Full guarantees and limitations are documented in the [performance audit](PERFORMANCE_AUDIT.md).

## Outputs, resources and recovery

```text
results_performance/
├── data/fastq/                       # persistent prepared mates; temporary units/ during assembly
├── data/{sample}_observed_whitelist.txt
├── logs/                            # preparation, Cutadapt JSON, QC and quantification logs
└── {sample}/
    ├── simpleaf/af_quant/alevin/     # native USA matrix and original row/column labels
    ├── {sample}_velocity.h5ad
    └── qc/                          # preparation receipt, FASTQ QC, conversion and validation
```

Prepared FASTQs remain persistent so quantification can be retried without repeating preparation. Unsampled lanes need no combined raw files or cropped R1 intermediate. Prepared lane files are Snakemake temporary outputs: completed lanes remain available while another lane fails, then are removed after successful assembly. Their small Cutadapt reports/logs persist. Multi-lane assembly briefly holds both unit and final prepared files and copies compressed bytes; a single unit is hard-linked without copying. Combined/downsampled raw files remain retained only on the downsampling path. Deleting retained intermediates saves disk but can force regeneration. Completed steps can be reused; interrupted/incomplete jobs may require Snakemake's `--rerun-incomplete`. Failed SimpleAF work is cleaned up and restarted at the quantification step, without internal mapping checkpoints. Prepared outputs, quantification and H5AD use same-filesystem renames after their respective checks. Multiple file renames are not a single transaction: downstream dependencies and the receipt determine readiness. Controlled pipeline failures clean staging files; a hard process/node kill can leave hidden staging directories for later review.

| Stage | CPU declaration | Memory reservation |
|---|---|---|
| Cutadapt + crop/cut, per unit | `nthreads`, capped by allocation; at least 2 | `preprocess_mem_mb`, default 8192 MB |
| Whitelist and FASTQ QC | 1 | `fastq_qc_mem_mb`, default 4096 MB |
| SimpleAF | `nthreads`, limited by available cores | `max_mem_mb`, default 64000 MB |
| H5AD export and validation | 1 | `export_mem_mb`, default 16000 MB |

A streaming unit assigns `threads - 1` Cutadapt workers and one pigz compression thread; I/O/helper processes add overhead. Independent lanes can run concurrently when CPU and memory budgets permit. Reservations guide scheduling and must fit the Slurm allocation; they are not measured production peak-memory guarantees. Whitelist RAM grows with distinct observed valid barcodes, while export/validation RAM grows with sparse matrix size. Whitelist/raw-scan progress is logged every ten million records. The [audit](PERFORMANCE_AUDIT.md) distinguishes removed passes from work that remains: wall-clock and peak-memory improvements have not been benchmarked on production samples.

## Development validation

These checks use deterministic fixtures, mocked tool calls, a placeholder-index dry run and a tiny synthetic reference. They do not read the configured biological samples:

```bash
pixi install --locked
export XDG_CACHE_HOME="$PWD/.cache"
export VELOCITY_SNAKEMAKE="$PWD/.pixi/envs/default/bin/snakemake"
export VELOCITY_SIMPLEAF="$PWD/.pixi/envs/default/bin/simpleaf"
export VELOCITY_CUTADAPT_PYTHON="$PWD/.pixi/envs/default/bin/python"
export VELOCITY_PREPROCESS_BIN="$PWD/.pixi/envs/default/bin"
pixi run --locked python -W error::DeprecationWarning -W error::FutureWarning tests/velocity/test_velocity.py
pixi run --locked python -m unittest discover -s tests -p test_cli.py
pixi run --locked python tests/velocity/synthetic_e2e.py
```

The local performance-change validation passed 24 velocity tests, two CLI tests and synthetic single-input, split-lane and downsampled workflows derived from 40 generated pairs, including completed-run and failed-lane resume checks. Tests compare decompressed prepared streams against upstream, identical whitelists/QC tables, malformed/truncated inputs, strict mate-name mismatches with multiple workers, pipeline failures, early chemistry/reference rejection, manifest bindings, stale receipts, BD conversion and serialized S+A/U counts. This establishes fixture behavior, not production reference performance or Slurm resource adequacy. The current graph index remains incomplete for some working-tree files; graph review is not a separate correctness certification.

## Provenance and documentation

Each workflow run records SHA256 of the runtime workflow sources, Pixi manifest/lock, canonical effective configuration and preparation QC receipt in native `pipeline_provenance.json`, then retains them in H5AD `uns['pipeline']`. Source paths and hashes are aligned arrays for H5AD compatibility. Sources, lock and receipt are checked before and after quantification; a change prevents publication. These are small-file hashes, not sequencing-file content checksums. Direct library calls without config/receipt metadata record those fields as `unavailable`.

Preparation originates from [imallona/rhapsodist, commit 7991f870fdfe38e12648673d330cbfd4839e755d](https://github.com/imallona/rhapsodist/tree/7991f870fdfe38e12648673d330cbfd4839e755d), under GPL-3.0; see [LICENSE.txt](LICENSE.txt) and [AUTHORS.txt](AUTHORS.txt). Crop/linker/error parameters, the line-wise header transformation and ordered BD panels are retained. Streaming, per-lane orchestration, early gates, JSON statistics, QC receipts and resource declarations are described in the [audit](PERFORMANCE_AUDIT.md). Tests require the upstream commit to remain available in Git history.

Documentation was revised with the scientific-writing and biomedical-proposal-writer skills. The scientific-writing guidance is part of Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. [doi:10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065). This credits the writing guidance; workflow behavior is supported by the linked source and tests.
