#!/usr/bin/env python3
"""Checked SimpleAF 0.30 interface. Never build or modify the reference."""
import argparse
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

GEOMETRY = '1{b[27]u[8]}2{r:}'
EXPECTED = {'simpleaf': '0.30.0', 'piscem': '0.23.0', 'alevin-fry': '0.18.3'}
WORKFLOW_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_SOURCES = [WORKFLOW_ROOT/name for name in (
    'workflow/Snakefile', 'workflow/src/prepare_fastq.snmk',
    'workflow/src/required_helpers.py', 'rhapsodist/cli.py', 'pixi.toml')]
WORKFLOW_SOURCES += sorted(Path(__file__).resolve().parent.glob('*.py'))


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            name = header.rstrip('\r\n')
            seq = handle.readline().rstrip('\r\n')
            plus = handle.readline().rstrip('\r\n')
            qual = handle.readline().rstrip('\r\n')
            if not name.startswith('@') or not plus.startswith('+') or len(seq) != len(qual) or not seq:
                raise ValueError(f'Malformed FASTQ record in {path}: {name}')
            yield name, seq


def paired_count(r1, r2, standardized=False):
    count = 0
    for left, right in zip_longest(records(r1), records(r2)):
        if left is None or right is None:
            raise ValueError('Different FASTQ mate counts')
        # Upstream cut removes linker columns from every R1 line, including headers.
        right_name = right[0]
        if standardized:
            right_name = right_name[:9] + right_name[13:22] + right_name[26:]
        a, b = left[0].split()[0], right_name.split()[0]
        if a.endswith('/1'):
            a = a[:-2]
        if b.endswith('/2'):
            b = b[:-2]
        if a != b:
            raise ValueError(f'FASTQ mates out of order: {a!r}, {b!r}')
        if standardized:
            check_geometry(left[1])
        count += 1
        if count % 10000000 == 0:
            print(f'Validated {count:,} read pairs', flush=True)
    return count


def check_geometry(sequence):
    if len(sequence) != 35:
        raise ValueError('Prepared R1 must contain exactly 27 barcode + 8 UMI bases')
    if sequence.strip('ACGTNacgtn'):
        raise ValueError('Prepared barcode/UMI contains non-DNA bases')


def file_identity(paths):
    # shortcut: stat identity detects ordinary edits, use exhaustive QC for archival integrity.
    result = []
    for path in paths:
        path = Path(path).resolve(strict=True)
        stat = path.stat()
        result.append([str(path), stat.st_size, stat.st_mtime_ns])
    return result


def prepare_qc(raw1, raw2, r1, r2, panels, whitelist, cutadapt_json,
               receipt, sample, downsampled=False):
    """Validate final R1 while deriving the unchanged observed whitelist."""
    if not raw1 or len(raw1) != len(raw2):
        raise ValueError('Raw FASTQ lists must be nonempty and paired')
    paths = [*raw1, *raw2, r1, r2, *panels, cutadapt_json]
    before = file_identity(paths)
    metadata = json.loads(Path(cutadapt_json).read_text())
    if metadata.get('cutadapt_version') != '4.9' or metadata.get('input', {}).get('paired') is not True:
        raise ValueError('Expected paired Cutadapt 4.9 report')
    preparation = metadata.get('preparation')
    if preparation is not None:
        if (preparation.get('version') != 2 or preparation.get('sample') != sample
                or preparation.get('downsampled') != downsampled
                or preparation.get('prepared_paths') != [str(Path(p).resolve()) for p in (r1, r2)]):
            raise ValueError('Wrong-sample or incompatible preparation report')
        if not downsampled and preparation['source_inputs'] != file_identity([*raw1, *raw2]):
            raise ValueError('Raw lane inputs disagree with preparation report')
    else:
        if not downsampled and len(raw1) == 1:
            reported = [Path(metadata['input'][key]).resolve() for key in ('path1', 'path2')]
            if reported != [Path(p).resolve() for p in (raw1[0], raw2[0])]:
                raise ValueError('Cutadapt report belongs to different raw FASTQs')
        arguments = metadata.get('command_line_arguments', [])
        if '-p' not in arguments or Path(arguments[arguments.index('-p')+1]).resolve() != Path(r2).resolve():
            raise ValueError('Cutadapt report does not describe final prepared R2')
    counts = metadata['read_counts']
    raw, written = counts['input'], counts['output']
    if any(type(n) is not int or n <= 0 for n in (raw, written)) or written > raw:
        raise ValueError('Invalid Cutadapt pair counts')
    scan_raw = downsampled or (len(raw1) > 1 and preparation is None)
    if scan_raw:
        print('Validating original raw lane pairs (downsampling/multiple lanes)', flush=True)
        original = sum(paired_count(a, b) for a, b in zip(raw1, raw2))
        if (not downsampled and original != raw) or original < raw:
            raise ValueError('Raw lane counts disagree with Cutadapt input')
        raw = original
    keys = [set(Path(p).read_text().splitlines()) for p in panels]
    if (len(keys) != 3 or any(len(k) != 384 or any(
            len(b) != 9 or set(b)-set('ACGT') for b in k) for k in keys)):
        raise ValueError('Invalid Enhanced barcode panels')
    observed = set()
    print('Validating prepared R1 and deriving observed whitelist', flush=True)
    # Cutadapt writes paired R2 directly; crop/cut preserves record order and count.
    prepared = 0
    for _, sequence in records(r1):
        check_geometry(sequence)
        bc = sequence[:27]
        if bc[:9] in keys[0] and bc[9:18] in keys[1] and bc[18:27] in keys[2]:
            observed.add(bc)
        prepared += 1
        if prepared % 10000000 == 0:
            print(f'Validated {prepared:,} prepared R1 records', flush=True)
    if prepared != written:
        raise ValueError('Final prepared pairs disagree with Cutadapt output')
    if not observed:
        raise ValueError('Empty observed whitelist')
    if file_identity(paths) != before:
        raise ValueError('FASTQ/QC inputs changed during validation')
    Path(whitelist).parent.mkdir(parents=True, exist_ok=True)
    Path(whitelist).write_text(''.join(b+'\n' for b in sorted(observed)))
    data = dict(version=1, sample=sample, inputs=before,
                whitelist_sha256=hashlib.sha256(Path(whitelist).read_bytes()).hexdigest(), input_read_pairs=raw,
                prepared_read_pairs=prepared, observed_valid_barcodes=len(observed),
                raw_count_source='original lane scan' if scan_raw else 'paired Cutadapt input per preparation unit',
                prepared_count_source='final R1 records matched to paired Cutadapt output',
                prepared_pair_validation='Cutadapt paired output + record-preserving crop/cut; R2 not reread')
    path = Path(receipt)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    os.replace(temporary, path)


def fastq_qc(raw1, raw2, r1, r2, whitelist, output, receipt=None,
             sample=None, panels=(), cutadapt_json=None, exhaustive=False):
    if len(raw1) != len(raw2) or not raw1:
        raise ValueError('Raw FASTQ lists must be nonempty and paired')
    if receipt is not None:
        data = json.loads(Path(receipt).read_text())
        paths = [*raw1, *raw2, r1, r2, *panels, cutadapt_json]
        if (data.get('version') != 1 or data.get('sample') != sample
                or data['inputs'] != file_identity(paths)
                or data['whitelist_sha256'] != hashlib.sha256(Path(whitelist).read_bytes()).hexdigest()):
            raise ValueError('Stale or wrong-sample QC receipt; rebuild whitelist/QC')
        raw, prepared = data['input_read_pairs'], data['prepared_read_pairs']
    if receipt is None or exhaustive:
        scanned_raw = sum(paired_count(a, b) for a, b in zip(raw1, raw2))
        scanned_prepared = paired_count(r1, r2, standardized=True)
        if receipt is not None and (raw, prepared) != (scanned_raw, scanned_prepared):
            raise ValueError('QC receipt counts disagree with exhaustive scan')
        raw, prepared = scanned_raw, scanned_prepared
    whitelist_values = Path(whitelist).read_text().splitlines()
    if not whitelist_values or len(whitelist_values) != len(set(whitelist_values)):
        raise ValueError('Empty or duplicated observed whitelist')
    if any(len(b) != 27 or set(b) - set('ACGT') for b in whitelist_values):
        raise ValueError('Observed whitelist must contain canonical 27nt sequences')
    if not raw or not prepared or prepared > raw:
        raise ValueError('Invalid raw/prepared read counts')
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['input_read_pairs', 'prepared_read_pairs', 'retained_percent',
                         'observed_valid_barcodes', 'corrected_barcodes'])
        writer.writerow([raw, prepared, 100*prepared/raw, len(whitelist_values), 'unavailable'])


def index_directory(root):
    """Resolve the reference root or a self-contained copied index by its metadata."""
    root = Path(root).resolve()
    return root if (root/'simpleaf_index.json').is_file() else root/'index'


def reference_files(root):
    index = index_directory(root)
    meta = json.loads((index/'simpleaf_index.json').read_text())
    if meta.get('index_type') != 'piscem':
        raise ValueError('Reference must be a SimpleAF piscem index')
    t2g = index/'t2g_3col.tsv'
    states = set()
    with t2g.open() as handle:
        for line in handle:
            row = line.rstrip('\n').split('\t')
            if len(row) != 3 or not row[0] or not row[1] or row[2] not in ('S', 'U'):
                raise ValueError('Expected a three-column spliceu transcript/gene/S-or-U map')
            states.add(row[2])
    if states != {'S', 'U'}:
        raise ValueError('Reference must contain both S and U targets')
    for suffix in ('ctab', 'refinfo', 'ssi'):
        if not (index/f'piscem_idx.{suffix}').is_file():
            raise ValueError(f'Missing piscem_idx.{suffix}')
    return index, t2g


def inspect_reference(root, manifest=None):
    """Validate reference content and bind optional declared identity to small-file hashes."""
    index = index_directory(root)
    t2g = index/'t2g_3col.tsv'
    paths = [index/'simpleaf_index.json', t2g,
             *[index/f'piscem_idx.{suffix}' for suffix in ('ctab', 'refinfo', 'ssi')]]
    if manifest is None and (index/'reference_manifest.json').is_file():
        manifest = index/'reference_manifest.json'
    if manifest is not None:
        paths.append(Path(manifest))
    before = file_identity(paths)
    reference_files(root)
    def digest(path):
        with open(path, 'rb') as handle:
            return hashlib.file_digest(handle, 'sha256').hexdigest()
    result = dict(reference_index_path=str(index), reference_t2g_path=str(t2g),
                  reference_t2g_sha256=digest(t2g),
                  reference_index_meta_sha256=digest(index/'simpleaf_index.json'),
                  genome='unknown', GENCODE_release='unknown', index_type='spliceu',
                  reference_identity_status='unverified: no reference manifest')
    if manifest is not None:
        data = json.loads(Path(manifest).read_text())
        if (not isinstance(data, dict) or data.get('version') != 1 or data.get('index_type') != 'spliceu'
                or any(not isinstance(data.get(k), str) or not data[k].strip()
                       for k in ('genome', 'GENCODE_release'))):
            raise ValueError('Invalid reference identity manifest')
        for key in ('reference_t2g_sha256', 'reference_index_meta_sha256'):
            if data.get(key) != result[key]:
                raise ValueError(f'Reference manifest hash mismatch: {key}')
        result.update(genome=data['genome'], GENCODE_release=data['GENCODE_release'],
                      reference_identity_status='manifest identity; t2g and index metadata hashes matched',
                      reference_manifest_path=str(Path(manifest).resolve()),
                      reference_manifest_sha256=digest(manifest))
    if file_identity(paths) != before:
        raise ValueError('Reference inputs changed during validation')
    result.update(version=1, reference_files=before)
    return result


def load_reference_validation(root, path):
    data = json.loads(Path(path).read_text())
    if (data.get('version') != 1 or data.get('reference_index_path') != str(index_directory(root))
            or data['reference_files'] != file_identity([p[0] for p in data['reference_files']])):
        raise ValueError('Stale or wrong-index reference validation; rerun reference preflight')
    return data


def quant_command(index, t2g, r1, r2, whitelist, output, threads):
    if threads < 2:
        raise ValueError('SimpleAF requires at least two threads')
    return ['simpleaf', 'quant', '--chemistry', GEOMETRY,
            '--index', str(index), '--t2g-map', str(t2g),
            '--reads1', str(r1), '--reads2', str(r2),
            '--unfiltered-pl', str(whitelist), '--min-reads', '1',
            '--expected-ori', 'both', '--resolution', 'cr-like',
            '--threads', str(threads), '--output', str(output)]


def execution_provenance(config_sha256=None, preparation_qc=None):
    """Hash small runtime sources/metadata, never sequencing files."""
    if (config_sha256 is None) != (preparation_qc is None):
        raise ValueError('Supply both configuration SHA256 and preparation receipt')
    if config_sha256 is not None and not re.fullmatch(r'[0-9a-f]{64}', config_sha256):
        raise ValueError('Invalid configuration SHA256')
    def digest(path):
        with open(path, 'rb') as handle:
            return hashlib.file_digest(handle, 'sha256').hexdigest()
    return dict(workflow_source_paths=[str(path.relative_to(WORKFLOW_ROOT)) for path in WORKFLOW_SOURCES],
                workflow_source_sha256=[digest(path) for path in WORKFLOW_SOURCES],
                workflow_lock_sha256=digest(WORKFLOW_ROOT/'pixi.lock'),
                workflow_config_sha256=config_sha256 or 'unavailable',
                preparation_qc_sha256=digest(preparation_qc) if preparation_qc else 'unavailable')


def run_quant(reference, r1, r2, whitelist, output, threads, log,
              upstream_provenance, reference_validation=None, config_sha256=None, preparation_qc=None):
    checked_reference = (load_reference_validation(reference, reference_validation)
                         if reference_validation else inspect_reference(reference))
    index, t2g = map(Path, (checked_reference['reference_index_path'], checked_reference['reference_t2g_path']))
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to replace quantification directory: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    r1, r2, whitelist = [Path(p).resolve(strict=True) for p in (r1, r2, whitelist)]
    execution = execution_provenance(config_sha256, preparation_qc)
    temporary = Path(tempfile.mkdtemp(prefix='.simpleaf-', dir=output.parent))
    try:
        home = temporary/'af_home'
        home.mkdir()
        env = dict(os.environ, ALEVIN_FRY_HOME=str(home))
        (home/'chemistries.json').write_text('{}\n')
        Path(log).parent.mkdir(parents=True, exist_ok=True)
        with Path(log).open('w') as handle:
            versions = {}
            for tool, expected in EXPECTED.items():
                actual = subprocess.check_output([tool, '--version'], env=env, text=True, stderr=subprocess.STDOUT).strip()
                if not re.search(r'(?<![\d.])'+re.escape(expected)+r'(?![\d.])', actual):
                    raise ValueError(f'Unsupported {tool} version: {actual}; expected {expected}')
                versions[tool] = actual
                handle.write(actual+'\n')
            help_text = subprocess.check_output(['simpleaf', 'quant', '--help'], env=env, text=True)
            for flag in ('--unfiltered-pl', '--min-reads', '--chemistry', '--t2g-map'):
                if flag not in help_text:
                    raise ValueError(f'Unsupported SimpleAF CLI: missing {flag}')
            subprocess.run(['simpleaf', 'set-paths', '--piscem', shutil.which('piscem') or 'piscem',
                                         '--alevin-fry', shutil.which('alevin-fry') or 'alevin-fry'],
                           env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
            command = quant_command(index, t2g, r1, r2, whitelist, temporary/'quant', threads)
            handle.write(json.dumps(command)+'\n')
            handle.flush()
            subprocess.run(command, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
        quant = temporary/'quant/af_quant'
        for name in ('quant.json', 'alevin/quants_mat.mtx', 'alevin/quants_mat_rows.txt', 'alevin/quants_mat_cols.txt'):
            if not (quant/name).is_file():
                raise ValueError(f'SimpleAF did not produce expected USA output: {quant/name}')
        if json.loads((quant/'quant.json').read_text()).get('usa_mode') is not True:
            raise ValueError('SimpleAF output does not confirm USA mode')
        if execution_provenance(config_sha256, preparation_qc) != execution:
            raise ValueError('Workflow source, lockfile or preparation receipt changed during quantification')
        meta = json.loads(Path(upstream_provenance).read_text())
        provenance = dict(rhapsodist_commit=meta['upstream_commit'], versions=versions, **execution,
                          **{k:v for k,v in checked_reference.items() if k not in ('version', 'reference_files')},
                          quantification_resolution='USA', umi_resolution='cr-like',
                          geometry=GEOMETRY, min_reads=1, expected_orientation='both',
                          command=command,
                          upstream_panel_sha256={Path(k).name:v for k,v in meta['sha256'].items() if 'BD_CLS' in k})
        correction = re.search(r'total number of distinct corrected barcodes\s*:\s*([\d,]+)', Path(log).read_text())
        if correction:
            provenance['corrected_distinct_barcodes'] = int(correction.group(1).replace(',', ''))
        (temporary/'quant/pipeline_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
        os.replace(temporary/'quant', output)
    finally:
        shutil.rmtree(temporary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'r1', 'r2', 'whitelist', 'output', 'log', 'upstream-provenance'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--reference-validation')
    parser.add_argument('--config-sha256')
    parser.add_argument('--preparation-qc')
    parser.add_argument('--threads', type=int, required=True)
    args = parser.parse_args()
    run_quant(args.reference, args.r1, args.r2, args.whitelist,
              args.output, args.threads, args.log, args.upstream_provenance, args.reference_validation,
              args.config_sha256, args.preparation_qc)
