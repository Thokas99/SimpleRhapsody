#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")"
export PIXI_CACHE_DIR="${PIXI_CACHE_DIR:-$PWD/.cache/pixi}"
export RATTLER_CACHE_DIR="$PIXI_CACHE_DIR"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$PWD/.cache/apptainer}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$PWD/.cache/build}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
pixi install --locked
export XDG_CACHE_HOME="$PWD/.cache"
export VELOCITY_SNAKEMAKE="$PWD/.pixi/envs/default/bin/snakemake"
export VELOCITY_SIMPLEAF="$PWD/.pixi/envs/default/bin/simpleaf"
export VELOCITY_CUTADAPT_PYTHON="$PWD/.pixi/envs/default/bin/python"
export VELOCITY_PREPROCESS_BIN="$PWD/.pixi/envs/default/bin"
pixi run --locked python3 -W error::DeprecationWarning -W error::FutureWarning tests/velocity/test_velocity.py
pixi run --locked python3 -m unittest discover -s tests -p test_cli.py
PYTHONWARNINGS=error::DeprecationWarning,error::FutureWarning pixi run --locked python3 tests/velocity/synthetic_e2e.py
# Build a new image before replacing any existing deployment artifact.
image="${1:-$PWD/rhapsodist-simpleaf.sif}"
test ! -e "$image" || { echo "Refusing to overwrite $image" >&2; exit 2; }
apptainer build --fakeroot "$image" rhapsodist.def
apptainer test "$image"
apptainer exec --cleanenv "$image" env PYTHONWARNINGS=error::DeprecationWarning,error::FutureWarning python3 /opt/rhapsodist/tests/velocity/synthetic_e2e.py
apptainer exec --cleanenv "$image" python3 -c 'import numpy, scipy, pandas, anndata, h5py, yaml'
(cd "$(dirname -- "$image")" && sha256sum "$(basename -- "$image")") > "$image.sha256"
stat --printf='%n: %s bytes\n' "$image"
