#!/usr/bin/env python3
"""Replay each exact checker, compare or write its receipt, hash its sources."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CASES = [('check_paired_budget.py', 'paired_budget.json'),
         ('check_five_interiors.py', 'five_interiors.json'),
         ('check_lift_gram.py', 'lift_gram.json'),
         ('check_depth_two.py', 'depth_two.json')]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    reports = []
    for script, receipt in CASES:
        command = [sys.executable, '-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/script)]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        data = json.loads(result.stdout)
        if args.record:
            (HERE/receipt).write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
        elif data != json.loads((HERE/receipt).read_text()):
            raise RuntimeError('Receipt changed: '+receipt)
        reports.append(dict(script=script, source_sha256=hashlib.sha256((HERE/script).read_bytes()).hexdigest(),
                            receipt=receipt, receipt_sha256=hashlib.sha256((HERE/receipt).read_bytes()).hexdigest()))
    inherited = [
        'root_ports_closed_carriers_20261028/check_root_rays.py',
        'root_ports_closed_carriers_20261028/check_carriers.py',
        'session87_four_vertices_enumeration_20261028/check_generation.py',
        'session88_symmetry_certificates_20261028/blockfill/check_blockfill.py',
        'session89_consecutive_tournament_pairs_20261028/generation/check_pairs.py',
        'session89_consecutive_tournament_pairs_20261028/generation/check_lift.py']
    return dict(status='PASS', exact_receipts=reports,
        inherited_source_sha256={p: hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest() for p in inherited},
        audit='same-agent alternate-path exact checks; not independent human review or a proof-assistant formalization',
        limits=['one-step ROOT theorem, not universal convergence',
                'five-interior paired-port architecture, not global plane chromatic number',
                'full order-nine self-converse image and scoped algebra, not an all-order canonical-label algorithm'])


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
