# Configuration

A config file lists the fields that differ from `configs/config.yaml`. The blocks below show each field with its default.

## Run

```yaml
working_dir: output   # relative paths start at the repository root
nthreads: 10          # threads per rule
max_mem_mb: 20000     # memory given to STAR
```

## References

Local files:

```yaml
gtf_origin: gencode                    # gencode or ensembl: header style of the transcriptome fasta
gtf: /path/to/annotation.gtf           # uncompressed
genome: /path/to/genome.fa             # uncompressed
transcriptome: /path/to/transcriptome.fa.gz
sjdbOverhang: 70                       # read 2 length minus 1
genomeSAindexNbases: 14                # lower for small genomes
```

Or URLs, downloaded and decompressed into `reference_dir`:

```yaml
reference_dir: data/reference/GRCh38_v46
genome_url: "https://.../GRCh38.p14.genome.fa.gz"
gtf_url: "https://.../gencode.v46.primary_assembly.annotation.gtf.gz"
transcriptome_url: "https://.../gencode.v46.transcripts.fa.gz"
```

## Samples

```yaml
samples:
  - name: my_sample
    uses:
      cb_umi_fq: /path/to/R1.fastq.gz   # or a list of files, concatenated in order
      cdna_fq: /path/to/R2.fastq.gz     # same length and order as cb_umi_fq
      # sra_run: SRR24978231            # instead of the two fastqs
      allowedlist: 384                  # 96 or 384 cell label sequences per section
      diversity_insets: yes             # yes or no
      downsample: 100                   # percentage of read pairs kept
      use_sampletags: yes               # overrides skip_sampletags
      species: mouse                    # human or mouse; needed for sample tags
```

The bead version is detected from the linkers of read 1; `allowedlist` and `diversity_insets` are checked against it. v1 beads are `96` without insets, Enhanced beads `96` with insets, Enhanced V2 beads `384` with insets.

## Aligners

```yaml
aligner: ['starsolo', 'kallisto', 'alevin']   # any of starsolo, kallisto, alevin, sbg
```

- `starsolo` aligns to the genome, counts reads with one gene, and writes a BAM file.
- `kallisto` and `alevin` map to the transcriptome and are faster. `bustools correct` needs about 20 GB of memory.
- `sbg` runs the BD Rhapsody pipeline v2.2.1 through Apptainer or Singularity and takes the longest.

Two or more aligners add a comparison report.

STARsolo:

```yaml
soloCellFilter: EmptyDrops_CR   # used with cell_filtering: native
soloMultiMappers: Unique        # Unique, Uniform, Rescue, PropUnique or EM
soloCBmatchWLtype: 1MM_multi    # whitelist matching of the cell barcode
soloUMIdedup: 1MM_CR
extraStarSoloArgs: ""
```

alevin:

```yaml
alevin_sketch: false   # true: alevin-fry with integer counts on the scale of STARsolo Unique
alevin_usa: false      # true: spliced, unspliced and ambiguous assays; needs alevin_sketch: true
```

BD Rhapsody pipeline:

```yaml
sbg_cwl: third_party/cwl/v2.2.1/rhapsody_pipeline_2.2.1.cwl
# sbg_reference_archive: /path/to/Rhapsody_reference.tar.gz   # default: built from the STAR index and GTF
# sbg_reference_url: "http://..."
```

## Barcodes and cells

```yaml
cb_umi_max_errors: 0     # mismatches allowed in the read 1 linkers
cell_filtering: native   # native, emptydrops or none
downsample: 100          # percentage of read pairs kept, for every sample
downsample_seed: 42
```

- `native`: STARsolo uses `soloCellFilter`; kallisto and alevin use the knee of DropletUtils `barcodeRanks`.
- `emptydrops`: DropletUtils `emptyDrops` for every aligner, so the cell numbers are comparable.
- `none`: every observed barcode is kept.

The linker report `linker_qc/{sample}_linker_qc.html` gives the fraction of reads kept at each `cb_umi_max_errors`.

## Sample tags

```yaml
skip_sampletags: false   # true: no sample tag calling
```

Sample tags are read from the WTA library. For a separate tag library, add its reads to the sample:

```yaml
      sampletag_cb_umi_fq: /path/to/sampletag_R1.fastq.gz
      sampletag_cdna_fq: /path/to/sampletag_R2.fastq.gz
      # sampletag_sra_run: SRR...
```

To match only the tags of a sample and split the counts per tag, list them, or map them to labels:

```yaml
      sampletags: [1, 2, 3]
      # sampletags: {1: donor_a, 2: donor_b}
```

STARsolo provides the tag counts when it runs; otherwise, and for a separate tag library, an alignment-free search of the reads does. Human and mouse tag sets are included.

## Outputs and reports

```yaml
output_format: sce                 # sce, h5ad or both
sampletag_split_backend: memory    # memory, or delayed for low memory and slower splits
run_biology_report: false
biology_markers_file: data/markers/skin_markers.tsv   # marker and cell_type columns; relative to workflow/
use_case: sendoel                  # sendoel (marker panels), hela (cell cycle phase) or trajectory
trajectory_root_marker: ""         # gene marking the earliest state, for use_case: trajectory
```
