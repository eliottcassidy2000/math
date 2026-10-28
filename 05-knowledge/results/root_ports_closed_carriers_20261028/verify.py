#!/usr/bin/env python3
"""Reproduce deterministic receipts in normal and optimized Python modes."""
from pathlib import Path
from hashlib import sha256
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = (
    ('check_root_rays.py', (), 'root_receipt.json'),
    ('check_selectors.py', (), 'selector_receipt.json'),
    ('check_bank_bridge.py', (), 'bank_bridge_receipt.json'),
    ('check_carriers.py', ('--max-order', '9'), 'carrier_receipt.json'),
    ('check_paired_probe.py', (), 'paired_probe_receipt.json'),
)


def main():
    outputs = {}; report = {}
    for script, args, receipt in JOBS:
        values = []
        for options in ((), ('-O',)):
            command = [sys.executable, '-B', *options, str(HERE/script), *args]
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            json.loads(result.stdout)
            values.append(result.stdout)
        if values[0] != values[1]:
            raise RuntimeError(('normal/optimized mismatch', script))
        outputs[receipt] = values[0]
        report[script] = dict(normal_optimized_byte_identical=True,
                              receipt_sha256=sha256(values[0]).hexdigest())
    # Write only after every complete normal/optimized computation passes.
    for receipt, raw in outputs.items():
        (HERE/receipt).write_bytes(raw)
    dependency_paths = [HERE.parent/'session87_four_vertices_enumeration_20261028/check_generation.py',
                        HERE.parent/'relation_tower_frobenius_ports_20261028/tournament_bank.json.gz']
    dependencies = {str(p.relative_to(HERE.parents[2])): sha256(p.read_bytes()).hexdigest()
                    for p in dependency_paths}
    own = {p.name: sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))}
    result = dict(status='PASS', semantic_jobs=report, dependency_sha256=dependencies,
                  own_script_sha256=own, independent_agent_audit=False,
                  audit_scope='same-agent proofs, literal versus regrouped arithmetic, exact reference-bank bijection, independent full-labelled/arc searches')
    (HERE/'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
