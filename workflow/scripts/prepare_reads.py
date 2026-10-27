#!/usr/bin/env python3
"""Stream the upstream enhanced crop/cut transform and assemble ordered lane outputs."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import sys

from simpleaf_quant import file_identity

ADAPTER = 'NNNNNNNNNGTGANNNNNNNNNGACANNNNNNNNNNNNNNNNN;min_overlap=43;noindels'
CUT_COLUMNS = '1-9,14-22,27-'


def strict_cutadapt(arguments):
    """Keep legacy mate-name checks inside the pinned Cutadapt parser's existing pass."""
    import multiprocessing
    import cutadapt
    from dnaio.pairedend import TwoFilePairedEndReader
    if cutadapt.__version__ != '4.9' or multiprocessing.get_start_method() != 'fork':
        raise RuntimeError('Strict streaming requires pinned Cutadapt 4.9 and Linux fork workers')
    from cutadapt.cli import main
    original = TwoFilePairedEndReader.__iter__
    def checked_pairs(reader):
        for left, right in original(reader):
            a, b = left.name.split()[0], right.name.split()[0]
            if a.endswith('/1'):
                a = a[:-2]
            if b.endswith('/2'):
                b = b[:-2]
            if a != b:
                raise ValueError(f'FASTQ mates out of order: {a!r}, {b!r}')
            yield left, right
    # Forked Cutadapt workers inherit this check before adapter filtering.
    TwoFilePairedEndReader.__iter__ = checked_pairs
    try:
        main(arguments)
    finally:
        TwoFilePairedEndReader.__iter__ = original


def standardize(r1, r2, output1, output2, report, log, sample, unit, threads, max_errors=0):
    if threads < 2:
        raise ValueError('Streaming preparation requires at least two allocated cores')
    before = file_identity([r1, r2])
    output1, output2, report = map(Path, (output1, output2, report))
    for path in (output1, output2, report, Path(log)):
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=output1.parent) as directory:
        directory = Path(directory)
        left, right, stats = [directory/name for name in ('r1.fq.gz', 'r2.fq.gz', 'cutadapt.json')]
        command = [sys.executable, str(Path(__file__).resolve()), '--cutadapt', '-g', ADAPTER, '--json', str(stats), '--action=crop',
                   '--discard-untrimmed', '--cores', str(threads-1), '-e', str(max_errors),
                   '--pair-filter=both', '-o', '-', '-p', str(right), str(r1), str(r2)]
        pipeline = (shlex.join(command) + ' | cut -c' + CUT_COLUMNS + ' | pigz -p 1')
        with left.open('wb') as reads, Path(log).open('w') as handle:
            handle.write(f'Sample {sample}, preparation unit {unit}\n{pipeline}\n')
            handle.flush()
            # Cutadapt sends its report to stderr when FASTQ is written to stdout.
            subprocess.run(['bash', '-o', 'pipefail', '-c', pipeline],
                           stdout=reads, stderr=handle, check=True)
        if file_identity([r1, r2]) != before:
            raise ValueError('Raw inputs changed during preparation')
        metadata = json.loads(stats.read_text())
        metadata['rhapsodist_unit'] = dict(version=2, sample=sample, unit=unit,
                                          raw_pair_validation='legacy_exact', source_inputs=before)
        stats.write_text(json.dumps(metadata, indent=2)+'\n')
        for source, destination in ((left, output1), (right, output2), (stats, report)):
            os.replace(source, destination)


def combine_prepared(r1s, r2s, reports, raw1, raw2, output1, output2, report,
                     sample, units, downsampled=False):
    """Concatenate gzip members without decompression; link a single unit."""
    if not units or len({len(x) for x in (r1s, r2s, reports, raw1, raw2, units)}) != 1:
        raise ValueError('Preparation units and mate/report lists must be nonempty and paired')
    metadata = []
    for unit, path, left, right in zip(units, reports, raw1, raw2):
        data = json.loads(Path(path).read_text())
        expected = dict(version=2, sample=sample, unit=unit,
                        raw_pair_validation='legacy_exact', source_inputs=file_identity([left, right]))
        if data.get('rhapsodist_unit') != expected:
            raise ValueError('Stale, reordered or wrong-sample preparation unit')
        if data.get('cutadapt_version') != '4.9' or data.get('input', {}).get('paired') is not True:
            raise ValueError('Expected paired Cutadapt 4.9 unit report')
        counts = data['read_counts']
        if (any(type(counts[k]) is not int or counts[k] < 0 for k in ('input', 'output'))
                or counts['output'] > counts['input']):
            raise ValueError('Invalid Cutadapt unit counts')
        metadata.append(data)
    output1, output2, report = map(Path, (output1, output2, report))
    for path in (output1, output2, report):
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.assemble-', dir=output1.parent) as directory:
        directory = Path(directory)
        for sources, destination, name in ((r1s, output1, 'r1.gz'), (r2s, output2, 'r2.gz')):
            temporary = directory/name
            if len(sources) == 1:
                os.link(sources[0], temporary)
            else:
                with temporary.open('wb') as handle:
                    for source in sources:
                        with open(source, 'rb') as reads:
                            shutil.copyfileobj(reads, handle, length=1024*1024)
            os.replace(temporary, destination)
        result = dict(cutadapt_version='4.9', input=dict(paired=True),
                      read_counts={k: sum(m['read_counts'][k] for m in metadata)
                                   for k in ('input', 'output')},
                      preparation=dict(version=2, sample=sample, downsampled=downsampled,
                                       source_inputs=file_identity([*raw1, *raw2]),
                                       prepared_paths=[str(p.resolve()) for p in (output1, output2)],
                                       units=metadata))
        temporary = directory/'report.json'
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        os.replace(temporary, report)


if __name__ == '__main__':
    if sys.argv[1:2] == ['--cutadapt']:
        strict_cutadapt(sys.argv[2:])
        sys.exit(0)
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('r1', 'r2', 'output1', 'output2', 'report', 'log', 'sample', 'unit'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--max-errors', type=float, default=0)
    args = parser.parse_args()
    standardize(args.r1, args.r2, args.output1, args.output2, args.report, args.log,
                args.sample, args.unit, args.threads, args.max_errors)
