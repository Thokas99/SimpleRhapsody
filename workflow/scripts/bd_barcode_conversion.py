#!/usr/bin/env python3
"""Enhanced V2/V3 BD integer IDs <-> 27nt sequences; standard library only."""
import argparse
import csv
from pathlib import Path


def load_keys(directory):
    keys = []
    for position in (1, 2, 3):
        values = (directory / f"BD_CLS{position}.txt").read_text().splitlines()
        if (len(values) != 384 or len(set(values)) != 384
                or any(len(s) != 9 or set(s) - set("ACGT") for s in values)):
            raise ValueError(f"Invalid ordered 384 x 9nt table at position {position}")
        keys.append(values)
    return keys


def encode(value, keys):
    value = str(value)
    if not value.isascii() or not value.isdecimal() or not 1 <= int(value) <= 384**3:
        raise ValueError(f"Invalid BD integer cell index: {value!r}")
    n = int(value) - 1
    return keys[0][n // 384**2] + keys[1][(n // 384) % 384] + keys[2][n % 384]


def decode(sequence, inverse):
    sequence = sequence.strip().upper()
    if len(sequence) != 27 or set(sequence) - set("ACGT"):
        raise ValueError(f"Invalid 27nt barcode: {sequence!r}")
    try:
        i, j, k = (inverse[p][sequence[p*9:(p+1)*9]] for p in range(3))
    except KeyError as exc:
        raise ValueError(f"Barcode segment absent from the supplied tables: {sequence}") from exc
    return str(i * 384**2 + j * 384 + k + 1)


def read_mapping(path):
    records = []
    bd_seen, barcode_seen = set(), set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"bd_cell_index", "barcode_27nt"}.issubset(reader.fieldnames or []):
            raise ValueError("Mapping requires bd_cell_index and barcode_27nt columns")
        for row in reader:
            value = row["bd_cell_index"]
            if not value.isascii() or not value.isdecimal() or not 1 <= int(value) <= 384**3:
                raise ValueError(f"Invalid BD cell index: {value!r}")
            bd, barcode = str(int(value)), row["barcode_27nt"].strip().upper()
            if len(barcode) != 27 or set(barcode) - set("ACGT"):
                raise ValueError(f"Invalid barcode: {barcode!r}")
            if bd in bd_seen or barcode in barcode_seen:
                raise ValueError(f"Duplicate mapping: {bd}, {barcode}")
            bd_seen.add(bd)
            barcode_seen.add(barcode)
            records.append((bd, barcode))
    if not records:
        raise ValueError("Mapping is empty")
    return records


def write_tsv(path, columns, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keys", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True,
                        help="One ID per line, no header; no implicit suffix removal")
    parser.add_argument("--direction", choices=["bd-to-barcode", "barcode-to-bd"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    keys = load_keys(args.keys)
    inverse = [{s: i for i, s in enumerate(key)} for key in keys]
    records = []
    bd_seen, barcode_seen = set(), set()
    for value in args.input.read_text().splitlines():
        if args.direction == "bd-to-barcode":
            barcode = encode(value, keys)
            bd = str(int(value))
        else:
            bd = decode(value, inverse)
            barcode = encode(bd, keys)
        if bd in bd_seen or barcode in barcode_seen:
            raise ValueError(f"Duplicate input identity: {value!r}")
        bd_seen.add(bd)
        barcode_seen.add(barcode)
        records.append((bd, barcode))
    if not records:
        raise ValueError("Input cell list is empty")
    if args.direction == "barcode-to-bd":
        write_tsv(args.output, ["barcode_27nt", "bd_cell_index"],
                  [(barcode, bd) for bd, barcode in records])
    else:
        write_tsv(args.output, ["bd_cell_index", "barcode_27nt"], records)
    print(f"Wrote {len(records):,} identities to {args.output}")




def convert_barcodes(barcodes, keys_dir, report_path, reference=None, differences=None):
    """Convert exact canonical barcodes; retain cells absent from a BD reference.

    A reference is either headerless decimal cell IDs, or a TSV with
    bd_cell_index and bd_cell_id columns supplying full, explicit identifiers.
    """
    if differences:
        Path(differences).parent.mkdir(parents=True, exist_ok=True)
    keys = load_keys(Path(keys_dir))
    inverse = [{s: i for i, s in enumerate(key)} for key in keys]
    report = dict(total_barcodes=len(barcodes), valid_barcodes=0,
                  converted_barcodes=0, unknown_barcodes=0, duplicate_ids=0,
                  roundtrip_failures=0, bd_reference_matches='unverified',
                  bd_reference_unmatched='unverified', reference_only='unverified',
                  naming='decimal Cell_index', reference_compatibility='unverified')
    if reference and not Path(reference).is_file():
        report['reference_compatibility'] = 'unverified: missing reference file'
        reference = None
    indexes, errors = [], []
    for barcode in barcodes:
        try:
            if len(barcode) != 27 or set(barcode) - set('ACGT'):
                raise ValueError(f'Noncanonical barcode: {barcode!r}')
            index = decode(barcode, inverse)
            report['valid_barcodes'] += 1
            if encode(index, keys) != barcode:
                report['roundtrip_failures'] += 1
                raise ValueError(f'Round-trip failed: {barcode}')
            indexes.append(index)
            report['converted_barcodes'] += 1
        except ValueError as exc:
            report['unknown_barcodes'] += 1
            errors.append(str(exc))
    report['duplicate_ids'] = len(indexes) - len(set(indexes))
    if report['duplicate_ids']:
        errors.append('Duplicate converted BD indexes')
    names = indexes.copy()
    if reference:
        lines = Path(reference).read_text().splitlines()
        if not lines:
            raise ValueError('BD reference is empty')
        if lines[0].split('\t') == ['bd_cell_index', 'bd_cell_id']:
            rows = list(csv.DictReader(lines, delimiter='\t'))
            ref_indices = [r['bd_cell_index'] for r in rows]
            ref_names = [r['bd_cell_id'] for r in rows]
            if len(set(ref_names)) != len(ref_names) or any(not n for n in ref_names):
                raise ValueError('Empty or duplicate full BD reference IDs')
            lookup = dict(zip(ref_indices, ref_names))
            if any(index not in lookup for index in indexes):
                errors.append('Explicit full-ID mapping must cover every quantifier cell')
            else:
                names = [lookup[index] for index in indexes]
            report['naming'] = 'explicit bd_cell_id from reference TSV'
        else:
            ref_indices = lines
            ref_names = lines
        for index in ref_indices:
            encode(index, keys)
            if str(int(index)) != index:
                raise ValueError(f'Noncanonical BD reference index: {index!r}')
        if len(set(ref_indices)) != len(ref_indices):
            raise ValueError('Duplicate BD reference indexes')
        matched = set(indexes) & set(ref_indices)
        report.update(bd_reference_matches=len(matched),
                      bd_reference_unmatched=len(set(indexes) - set(ref_indices)),
                      reference_only=len(set(ref_indices) - set(indexes)),
                      reference_compatibility='compared; unmatched retained')
        if differences:
            Path(differences).parent.mkdir(parents=True, exist_ok=True)
            with Path(differences).open('w', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['status', 'bd_cell_index'])
                writer.writerows(('quantifier_only', i) for i in sorted(set(indexes)-set(ref_indices), key=int))
                writer.writerows(('reference_only', i) for i in sorted(set(ref_indices)-set(indexes), key=int))
    elif differences:
        Path(differences).write_text('status\tbd_cell_index\nunverified\t\n')
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(report_path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(report.keys())
        writer.writerow(report.values())
    if errors:
        raise ValueError('; '.join(errors[:10]))
    return indexes, names, report


if __name__ == "__main__":
    main()
