"""Runs paper/zumis/run_zumis.sh with a stub zUMIs: a count table left by an
earlier run does not count as output of a failed run."""

import gzip
import os
import shutil
import stat
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, 'paper', 'zumis', 'run_zumis.sh')

pytestmark = pytest.mark.skipif(shutil.which('pigz') is None, reason='needs pigz')


def make_inputs(tmp_path, zumis_body):
    zumis_dir = tmp_path / 'zumis'
    (zumis_dir / 'zUMIs-env' / 'bin').mkdir(parents=True)
    zumis = zumis_dir / 'zUMIs.sh'
    zumis.write_text('#!/usr/bin/env bash\n' + zumis_body)
    zumis.chmod(zumis.stat().st_mode | stat.S_IEXEC)
    for mate in ('R1', 'R2'):
        with gzip.open(tmp_path / f'{mate}.fastq.gz', 'wt') as fh:
            fh.write('@r comment\nACGT\n+\nIIII\n')
    (tmp_path / 'genome.fa').write_text('>c\nACGT\n')
    (tmp_path / 'genes.gtf').write_text('')
    workdir = tmp_path / 'work'
    (workdir / 'star_index').mkdir(parents=True)
    (workdir / 'star_index' / 'SA').write_text('')
    return zumis_dir, workdir


def run_wrapper(tmp_path, zumis_body):
    zumis_dir, workdir = make_inputs(tmp_path, zumis_body)
    counts = workdir / 'zUMIs_output' / 'expression' / 's.dgecounts.rds'
    counts.parent.mkdir(parents=True)
    counts.write_text('earlier run')
    result = subprocess.run(
        ['bash', SCRIPT, str(zumis_dir), str(workdir), 's', str(tmp_path / 'R1.fastq.gz'),
         str(tmp_path / 'R2.fastq.gz'), str(tmp_path / 'genome.fa'), str(tmp_path / 'genes.gtf'),
         '1', '1'],
        capture_output=True, text=True, check=False)
    return result, counts


def test_failed_run_with_an_earlier_count_table_fails(tmp_path):
    result, counts = run_wrapper(tmp_path, 'exit 1\n')
    assert result.returncode != 0
    assert not counts.exists()


def test_successful_run_keeps_its_count_table(tmp_path):
    body = 'mkdir -p zUMIs_output/expression\necho new > zUMIs_output/expression/s.dgecounts.rds\n'
    result, counts = run_wrapper(tmp_path, body)
    assert result.returncode == 0, result.stderr
    assert counts.read_text() == 'new\n'


def test_nonzero_exit_with_a_count_table_fails(tmp_path):
    body = 'mkdir -p zUMIs_output/expression\necho partial > zUMIs_output/expression/s.dgecounts.rds\nexit 2\n'
    result, _ = run_wrapper(tmp_path, body)
    assert result.returncode != 0
    assert 'exited with 2' in result.stderr
