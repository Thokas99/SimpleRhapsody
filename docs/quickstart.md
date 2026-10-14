# Quick start

## Simulated data

From the repository root, a run on 1,000 simulated cells with three aligners:

```
rhapsodist --configfile configs/sim_config.yaml --config "aligner=['starsolo','kallisto','alevin']" --cores 4
```

Results are written to `output/simul`. The first run takes longer, while Snakemake creates the conda environments.

## Your data

Write a config with the fields that differ from the defaults in `configs/config.yaml`:

```yaml
working_dir: output/my_experiment
nthreads: 10
gtf_origin: gencode
reference_dir: data/reference/GRCh38_v46
genome_url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz"
gtf_url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
transcriptome_url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz"
sjdbOverhang: 70
aligner: ['starsolo', 'kallisto', 'alevin']
skip_sampletags: true
samples:
  - name: my_sample
    uses:
      cb_umi_fq: /path/to/R1.fastq.gz
      cdna_fq: /path/to/R2.fastq.gz
```

Check the planned jobs, then run:

```
rhapsodist --configfile my_config.yaml --cores 10 --dry-run
rhapsodist --configfile my_config.yaml --cores 10
```

## Containers

Run from the directory that holds your config; outputs are written there.

```
docker run --rm --user "$(id -u):$(id -g)" -e USER="$USER" -v "$PWD":/work rhapsodist --configfile /opt/rhapsodist/configs/sim_config.yaml --config "aligner=['starsolo','kallisto','alevin']" working_dir=/work/output/simul --cores 4

apptainer run rhapsodist.sif --configfile /opt/rhapsodist/configs/sim_config.yaml --config "aligner=['starsolo','kallisto','alevin']" working_dir="$PWD/output/simul" --cores 4
```
