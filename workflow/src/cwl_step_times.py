#!/usr/bin/env python3
"""Per-step wall-clock times from a cwltool log written with --timestamps.

A scattered step logs one start line per job and one completion line, so a step
runs from its first start to its last completion. Only top-level steps are
listed: a nested workflow is timed as its own step, not as its inner steps.
"""

import argparse
import re
from datetime import datetime

ANSI_ESCAPE = re.compile(r'\x1b\[[0-9;]*m')
STEP_EVENT = re.compile(
    r'^\[(?P<time>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\] \w+ '
    r'\[step (?P<step>[^\]]+)\] (?P<event>start|completed (?P<status>\w+))$'
)
NESTED_STEP = re.compile(r'\[workflow [^\]\s][^\]]*\] starting step (?P<step>\S+)$')
TIME_FORMAT = '%Y-%m-%d %H:%M:%S'


def parse_step_times(lines):
    """Return one dict per step, in order of first start.

    A step with no completion line keeps an empty end, seconds and status.
    """
    steps = {}
    nested = set()
    for line in lines:
        line = ANSI_ESCAPE.sub('', line).strip()
        nested_match = NESTED_STEP.search(line)
        if nested_match:
            nested.add(nested_match['step'])
        match = STEP_EVENT.match(line)
        if not match:
            continue
        record = steps.setdefault(
            match['step'],
            {'step': match['step'], 'start': '', 'end': '', 'seconds': '', 'status': ''})
        if match['event'] == 'start':
            record['start'] = record['start'] or match['time']
        else:
            record['end'] = match['time']
            record['status'] = match['status']
    for record in steps.values():
        if record['start'] and record['end']:
            elapsed = (datetime.strptime(record['end'], TIME_FORMAT)
                       - datetime.strptime(record['start'], TIME_FORMAT))
            record['seconds'] = int(elapsed.total_seconds())
    return [record for name, record in steps.items() if name not in nested]


def write_step_times(records, path):
    columns = ['step', 'start', 'end', 'seconds', 'status']
    with open(path, 'w') as fh:
        fh.write('\t'.join(columns) + '\n')
        for record in records:
            fh.write('\t'.join(str(record[c]) for c in columns) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--log', required=True, help='cwltool log written with --timestamps')
    ap.add_argument('--out', required=True, help='output TSV')
    args = ap.parse_args()
    with open(args.log, errors='replace') as fh:
        records = parse_step_times(fh)
    write_step_times(records, args.out)


if __name__ == '__main__':
    main()
