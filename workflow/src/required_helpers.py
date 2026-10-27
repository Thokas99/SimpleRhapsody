#!/usr/bin/env python
## Extracted from imallona/rhapsodist at 7991f870fdfe38e12648673d330cbfd4839e755d.
## GPLv3; retain upstream attribution.
import os.path as op

def get_sample_names():
    return [sample['name'] for sample in config['samples']]


def _sample_uses(name):
    return next((sample['uses'] for sample in config['samples']
                 if sample['name'] == name), None)


def _as_fastq_list(value):
    if value is None:
        return None
    return [value] if isinstance(value, str) else list(value)


def _combined_fastq(name, mate):
    return op.join(config['working_dir'], 'data', 'fastq', 'combined',
                   f'{name}_{mate}.fastq.gz')


def get_cbumi_inputs(name):
    return _as_fastq_list((_sample_uses(name) or {}).get('cb_umi_fq'))


def get_cdna_inputs(name):
    return _as_fastq_list((_sample_uses(name) or {}).get('cdna_fq'))


def _raw_fastq_path(name, field, mate):
    uses = _sample_uses(name) or {}
    files = _as_fastq_list(uses.get(field))
    if len(files) > 1:
        return _combined_fastq(name, mate)
    return files[0]


def _raw_cbumi_path(name):
    return _raw_fastq_path(name, 'cb_umi_fq', 'R1')


def _raw_cdna_path(name):
    return _raw_fastq_path(name, 'cdna_fq', 'R2')


def validate_fastq_lists():
    for sample in config['samples']:
        uses = sample['uses']
        if 'cb_umi_fq' not in uses or 'cdna_fq' not in uses:
            continue
        r1 = _as_fastq_list(uses['cb_umi_fq'])
        r2 = _as_fastq_list(uses['cdna_fq'])
        if len(r1) != len(r2):
            raise ValueError(
                f"Sample '{sample['name']}': cb_umi_fq has {len(r1)} file(s) but "
                f"cdna_fq has {len(r2)}. Both must list the same number of "
                "fastqs in matching order."
            )


def get_downsample_fraction(name):
    uses = _sample_uses(name) or {}
    pct = float(uses.get('downsample', config.get('downsample', 100)))
    if pct <= 0 or pct > 100:
        raise ValueError(f"Sample '{name}': downsample must be in (0, 100], got {pct}")
    return pct / 100.0


def _downsampled_fastq(name, mate):
    return op.join(config['working_dir'], 'data', 'fastq', 'downsampled',
                   f'{name}_{mate}.fastq.gz')


def get_cbumi_by_name(name):
    return (_downsampled_fastq(name, 'R1') if get_downsample_fraction(name) < 1.0
            else _raw_fastq_path(name, 'cb_umi_fq', 'R1'))


def get_cdna_by_name(name):
    return (_downsampled_fastq(name, 'R2') if get_downsample_fraction(name) < 1.0
            else _raw_fastq_path(name, 'cdna_fq', 'R2'))


def get_declared_allowedlist(name):
    value = (_sample_uses(name) or {}).get('allowedlist')
    if value is None:
        return None
    value = str(value)
    if value not in ('96', '384'):
        raise ValueError(f"Sample '{name}': allowedlist must be 96 or 384, got '{value}'")
    return value


def _yesno(name, field, value):
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    value = str(value).lower()
    if value in ('yes', 'true'):
        return 'yes'
    if value in ('no', 'false'):
        return 'no'
    raise ValueError(f"Sample '{name}': {field} must be yes or no, got '{value}'")


def get_declared_diversity_insets(name):
    value = (_sample_uses(name) or {}).get('diversity_insets')
    if value is None:
        return None
    return _yesno(name, 'diversity_insets', value)


def get_bead_type_by_name(name):
    allowedlist = get_declared_allowedlist(name)
    insets = get_declared_diversity_insets(name)
    if allowedlist is None and insets is None:
        return 'enhanced_v2'
    allowedlist = allowedlist or '384'
    insets = insets or 'yes'
    if allowedlist == '96' and insets == 'no':
        return 'v1'
    if allowedlist == '96' and insets == 'yes':
        return 'enhanced'
    if allowedlist == '384' and insets == 'yes':
        return 'enhanced_v2'
    raise ValueError(
        f"Sample '{name}': unsupported combination allowedlist={allowedlist}, "
        f"diversity_insets={insets} (valid: 96/no, 96/yes, 384/yes)"
    )


def get_barcode_whitelist_by_name(name):
    allowedlist = get_declared_allowedlist(name)
    if allowedlist is None:
        return '384x3' if get_bead_type_by_name(name) == 'enhanced_v2' else '96x3'
    return allowedlist + 'x3'


def scan_r1_linkers(cb_umi_path, n_reads=10000):
    """Return upstream linker-class fractions and error histograms for R1."""
    import gzip as _gzip
    v1_l1 = 'ACTGGCCTGCGA'; v1_l1_start = 9
    v1_l2 = 'GGTAGCGGTGACA'; v1_l2_start = 30
    enh_l1 = 'GTGA'; enh_l1_start = 9
    enh_l2 = 'GACA'; enh_l2_start = 22
    v1_l2_end = v1_l2_start + len(v1_l2)
    enh_l2_end = enh_l2_start + 3 + len(enh_l2)

    def _hamming(a, b):
        return sum(x != y for x, y in zip(a, b))

    total = 0
    v1_exact = 0
    enh_exact = 0
    v1_errs = {}
    enh_errs = {}
    paths = cb_umi_path if isinstance(cb_umi_path, (list, tuple)) else [cb_umi_path]
    def lines():
        for path in paths:
            opener = _gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rt') as handle:
                yield from handle
    from contextlib import closing
    with closing(lines()) as handle:
        for i, line in enumerate(handle):
            if i % 4 != 1:
                continue
            if total >= n_reads:
                break
            seq = line.strip()
            total += 1
            if len(seq) >= v1_l2_end:
                errors = (_hamming(seq[v1_l1_start:v1_l1_start + len(v1_l1)], v1_l1)
                          + _hamming(seq[v1_l2_start:v1_l2_start + len(v1_l2)], v1_l2))
                v1_errs[errors] = v1_errs.get(errors, 0) + 1
                if errors == 0:
                    v1_exact += 1
            best_enh = None
            for offset in range(4):
                s1, s2 = enh_l1_start + offset, enh_l2_start + offset
                if len(seq) >= s2 + len(enh_l2):
                    errors = (_hamming(seq[s1:s1 + len(enh_l1)], enh_l1)
                              + _hamming(seq[s2:s2 + len(enh_l2)], enh_l2))
                    if best_enh is None or errors < best_enh:
                        best_enh = errors
            if best_enh is not None:
                enh_errs[best_enh] = enh_errs.get(best_enh, 0) + 1
                if best_enh == 0:
                    enh_exact += 1
    if total == 0:
        return {'class': 'unknown', 'n_reads': 0, 'v1_frac': 0.0, 'enh_frac': 0.0,
                'v1_errors': {}, 'enh_errors': {}}
    v1_frac, enh_frac = v1_exact / total, enh_exact / total
    klass = ('v1' if v1_frac >= enh_frac and v1_frac > 0.1
             else 'enhanced' if enh_frac > 0.1 else 'unknown')
    return {'class': klass, 'n_reads': total, 'v1_frac': v1_frac,
            'enh_frac': enh_frac, 'v1_errors': v1_errs, 'enh_errors': enh_errs}
