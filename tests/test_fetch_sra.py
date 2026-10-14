"""Runs fetch_sra on stub SRA tools: a truncated compressed read file fails
the rule and keeps the uncompressed reads."""

import os
import shutil
import stat
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLE = 'hela_unmod'
SRA_RUN = 'SRR35038871'

pytestmark = pytest.mark.skipif(
    shutil.which('snakemake') is None or shutil.which('pigz') is None,
    reason='needs snakemake and pigz')

STUBS = {
    'prefetch': 'mkdir -p "$4/{run}"\n',
    'fasterq-dump': ('for mate in 1 2; do\n'
                     '  printf "@r\\nACGT\\n+\\nIIII\\n" > "$5/{run}_$mate.fastq"\n'
                     'done\n'),
    'pigz': ('if [ -n "${{TRUNCATE_MATE:-}}" ] && [[ "$4" == *_${{TRUNCATE_MATE}}.fastq ]]; then\n'
             '  "{pigz}" -c "$4" | head -c 10\n'
             'else\n'
             '  exec "{pigz}" "$@"\n'
             'fi\n'),
}


def write_stubs(bin_dir):
    bin_dir.mkdir()
    for name, body in STUBS.items():
        path = bin_dir / name
        path.write_text('#!/usr/bin/env bash\n' + body.format(run=SRA_RUN, pigz=shutil.which('pigz')))
        path.chmod(path.stat().st_mode | stat.S_IEXEC)


def run_fetch(tmp_path, truncate_mate=''):
    write_stubs(tmp_path / 'bin')
    outdir = tmp_path / 'data' / 'fastq' / 'sra' / SAMPLE
    env = dict(os.environ, PATH=f"{tmp_path / 'bin'}:{os.environ['PATH']}",
               TRUNCATE_MATE=truncate_mate)
    result = subprocess.run(
        ['snakemake', '--snakefile', os.path.join(REPO_ROOT, 'workflow', 'Snakefile'),
         '--configfile', os.path.join(REPO_ROOT, 'configs', 'moro_mallona2025_config.yaml'),
         '--config', f'working_dir={tmp_path}', 'nthreads=1',
         '--nolock', '--cores', '1', str(outdir / f'{SAMPLE}_R1.fastq.gz')],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=False)
    return result, outdir


def test_complete_reads_are_kept_compressed(tmp_path):
    result, outdir = run_fetch(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (outdir / f'{SAMPLE}_R2.fastq.gz').exists()
    assert not (outdir / f'{SRA_RUN}_1.fastq').exists()


def test_truncated_reads_fail_and_keep_fastq(tmp_path):
    result, outdir = run_fetch(tmp_path, truncate_mate='2')
    assert result.returncode != 0
    assert (outdir / f'{SRA_RUN}_1.fastq').exists()
    assert (outdir / f'{SRA_RUN}_2.fastq').exists()
