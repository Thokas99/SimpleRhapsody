#!/usr/bin/env bash
#SBATCH --job-name=bd-simpleaf
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/globalscratch/tsirchi/cemlab/simpleaf
#SBATCH --output=/globalscratch/tsirchi/cemlab/simpleaf/logs/velocity-%j.log
set -euo pipefail
WORKSPACE=/globalscratch/tsirchi/cemlab/simpleaf
CONFIG="$(realpath -- "${1:-$WORKSPACE/configs/velocity.yaml}")"
export APPTAINER_TMPDIR="$WORKSPACE/work"
cd "$WORKSPACE/work"
apptainer exec --cleanenv --bind /globalscratch \
  "$WORKSPACE/containers/rhapsodist-simpleaf.sif" \
  rhapsodist --configfile "$CONFIG" --cores "${SLURM_CPUS_PER_TASK:-16}" \
  --resources mem_mb=64000
