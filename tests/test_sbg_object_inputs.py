"""Dry runs of workflow/Snakefile: generate_sce_sbg reads the STARsolo
features.tsv, which gives the gene ids of the sbg object."""

import os
import re
import shutil
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STARSOLO_FEATURES = os.path.join('starsolo', 'simulated', 'Solo.out', 'Gene', 'filtered',
                                 'features.tsv')

pytestmark = pytest.mark.skipif(shutil.which('snakemake') is None, reason='needs snakemake')


def sbg_object_inputs(aligners, directory):
    result = subprocess.run(
        ['snakemake', '--snakefile', os.path.join(REPO_ROOT, 'workflow', 'Snakefile'),
         '--configfile', os.path.join(REPO_ROOT, 'configs', 'sim_config.yaml'),
         '--config', f'aligner={aligners}', f'working_dir={directory}',
         '--dry-run', '--nolock', '--cores', '1'],
        cwd=directory, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    rule = re.search(r'^rule generate_sce_sbg:\n\s+input: (.*)$', result.stdout, re.MULTILINE)
    assert rule, result.stdout
    return rule.group(1).split(', ')


def test_sbg_object_waits_for_starsolo(tmp_path):
    inputs = sbg_object_inputs(['starsolo', 'sbg'], tmp_path)
    assert any(path.endswith(STARSOLO_FEATURES) for path in inputs)


def test_sbg_object_alone_has_no_starsolo_input(tmp_path):
    inputs = sbg_object_inputs(['sbg'], tmp_path)
    assert not any('starsolo' in path for path in inputs)
