# Install

Requirements: Linux, conda or mamba, and Snakemake 8.12 or later. The BD Rhapsody pipeline (aligner `sbg`) also needs Apptainer or Singularity.

```
git clone https://github.com/imallona/rhapsodist
cd rhapsodist
pip install -e .
```

Install in editable mode (`-e`): `rhapsodist` runs the workflow files of the clone, which a regular `pip install .` does not copy. Snakemake creates the conda environments of the workflow at the first run.

## Containers

The images include Snakemake and the conda environments, without the BD Rhapsody pipeline.

```
docker build -t rhapsodist .
apptainer build rhapsodist.sif rhapsodist.def
```

Run them as shown in the [quick start](quickstart.md#containers).
