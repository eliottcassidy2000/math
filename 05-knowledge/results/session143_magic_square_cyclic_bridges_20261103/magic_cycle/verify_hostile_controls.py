#!/usr/bin/env python3
"""Portable destructive controls run only on disposable copies of this bundle."""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def run(path, optimized):
    args = [sys.executable] + (['-O'] if optimized else []) + [str(path)]
    return subprocess.run(args, capture_output=True, text=True, cwd='/tmp')


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


receipts = []
for program in ('verify.py', 'full_cubic_audit.py', 'audit/full_cubic_nonodd.py', 'verify_manifest.py'):
    normal, optimized = (run(HERE/program, flag) for flag in (False, True))
    need(normal.returncode == optimized.returncode == 0, 'baseline failed: '+program)
    need(normal.stdout == optimized.stdout, 'normal/optimized mismatch: '+program)

mutations = [
    ('source_snapshot_byte', 'sources/S57_square_totals.md', None, '\nHOSTILE MUTATION\n',
     ('verify.py', 'verify_manifest.py')),
    ('duplicated_cycle_state', 'verify.py', 'cycle = [4, 8, 3,', 'cycle = [4, 8, 2,',
     ('verify.py', 'verify_manifest.py')),
    ('damaged_cubic_witness', 'full_cubic_audit.py', 'cubic=[[4,-22,18]', 'cubic=[[5,-22,18]',
     ('full_cubic_audit.py', 'verify_manifest.py')),
]
with tempfile.TemporaryDirectory(prefix='magic-cycle-hostiles-') as temp:
    for name, relative, old, new, programs in mutations:
        target = Path(temp)/name
        shutil.copytree(HERE, target)
        file = target/relative
        original = file.read_text()
        if old is None:
            altered = original+new
        else:
            need(original.count(old) == 1, 'mutation target is not unique: '+name)
            altered = original.replace(old, new, 1)
        file.write_text(altered)
        for program in programs:
            failures = []
            for optimized in (False, True):
                result = run(target/program, optimized)
                need(result.returncode != 0, 'mutation incorrectly accepted: '+name+' / '+program)
                tail = result.stderr.strip().splitlines()[-1]
                need(tail.startswith('RuntimeError:'), 'unexpected failure mechanism: '+tail)
                failures.append({'optimized': optimized, 'returncode': result.returncode, 'failure': tail})
            need(failures[0]['failure'] == failures[1]['failure'], 'optimized hostile mechanism changed')
            receipts.append({'mutation': name, 'program': program, 'outcomes': failures})

print(json.dumps({'status': 'PASS', 'baseline_normal_optimized_pairs': 4,
                  'mutations': 3, 'expected_rejection_runs': 12,
                  'original_bundle_mutated': False, 'receipts': receipts}, indent=2, sort_keys=True))
