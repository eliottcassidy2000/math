#!/usr/bin/env python3
"""Validate frozen session evidence; --rerun freshly executes six checker pairs."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOBS = (
    ('h5/check.py', 'h5/check.json', 'h5/check.optimized.json', False),
    ('targets/check.py', 'targets/receipt.json', 'targets/receipt.optimized.json', False),
    ('golden/check.py', 'golden/CHECK.json', 'golden/CHECK-optimized.json', True),
    ('golden/audit_targets.py', 'golden/TARGET-AUDIT.json', 'golden/TARGET-AUDIT-optimized.json', True),
    ('h5/audit_golden.py', 'h5/audit_golden.json', 'h5/audit_golden.optimized.json', False),
    ('resets/check.py', 'resets/check.json', 'resets/check.optimized.json', False),
)

def require(q, message):
    if not q:
        raise RuntimeError(message)

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def object_bytes(ref):
    return subprocess.run(['git', 'show', ref], cwd=ROOT, capture_output=True, check=True).stdout

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    pins = set()
    historical = []
    for pins_file in (HERE/'SOURCE_PINS.json', HERE/'h5/SOURCE_PINS.json'):
        data = json.loads(pins_file.read_text())
        for name, expected in data.get('source_sha256', {}).items():
            path = ROOT/name
            if path.is_file() and digest(path) == expected:
                pins.add(name)
                continue
            # Mutable navigation/ledger files can change in another session.
            # The declared immutable source commit still pins exactly what was read.
            ref = data['source_HEAD'] + ':' + name
            require(sha256(object_bytes(ref)).hexdigest() == expected, ('historical source pin', ref))
            historical.append(ref)
            pins.add(name)
        for ref, expected in data.get('git_object_sha256', {}).items():
            require(sha256(object_bytes(ref)).hexdigest() == expected, ('Git source pin', ref))
            pins.add(ref)
    for script_name, normal_name, optimized_name, writes_output in JOBS:
        script, normal, optimized = (HERE/s for s in (script_name, normal_name, optimized_name))
        require(normal.read_bytes() == optimized.read_bytes(), ('optimized mismatch', script_name))
        data = json.loads(normal.read_text())
        require(str(data['status']).startswith(('PASS', 'FINITE-EXACT')), ('failed receipt', script_name))
        if 'script_sha256' in data:
            require(data['script_sha256'] == digest(script), ('changed script', script_name))
        for key in ('source_sha256', 'input_sha256'):
            for name, expected in data.get(key, {}).items():
                candidates = {p.resolve() for p in (ROOT/name, HERE/name) if p.is_file()}
                require(len(candidates) == 1, ('ambiguous or absent input', name))
                path = next(iter(candidates))
                require(digest(path) == expected, ('changed input', script_name, name))
                pins.add(str(path.relative_to(ROOT)))
        if 'bank_sha256' in data:
            bank = ROOT/'04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json'
            require(digest(bank) == data['bank_sha256'], ('changed strong bank', script_name))
            pins.add(str(bank.relative_to(ROOT)))
        if 'golden_source_file_sha256' in data:
            path = ROOT/'05-knowledge/results/golden_hamilton_transfer_20261103/verification.json'
            require(digest(path) == data['golden_source_file_sha256'], 'changed golden source control')
            pins.add(str(path.relative_to(ROOT)))
        if args.rerun:
            for optimized_mode in (False, True):
                with tempfile.TemporaryDirectory(prefix='chanpak-replay-') as directory:
                    output = Path(directory)/'receipt.json'
                    command = [sys.executable, '-B'] + (['-O'] if optimized_mode else []) + [str(script)]
                    if writes_output:
                        command += ['--output', str(output)]
                    run = subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
                    fresh = json.loads(output.read_text() if writes_output else run.stdout)
                    require(fresh == data, ('fresh replay differs', script_name, optimized_mode))
    links = 0
    for md in HERE.rglob('*.md'):
        for link in re.findall(r'\]\(([^\s)]+)\)', md.read_text()):
            if '://' in link or link.startswith('#'):
                continue
            target = (md.parent/link.split('#', 1)[0]).resolve()
            require(target.is_file(), ('broken local link', str(md.relative_to(HERE)), link))
            links += 1
    actual = {str(p.relative_to(HERE)): digest(p) for p in sorted(HERE.rglob('*'))
              if p.is_file() and p.name != 'manifest.json' and '__pycache__' not in p.parts}
    if args.write_manifest:
        (HERE/'manifest.json').write_text(json.dumps(actual, indent=2, sort_keys=True)+'\n')
    else:
        require(json.loads((HERE/'manifest.json').read_text()) == actual, 'manifest mismatch')
    print(json.dumps({'status': 'PASS', 'normal_optimized_pairs': len(JOBS),
                      'pinned_inputs': len(pins), 'historical_pins': historical,
                      'local_links': links, 'artifacts': len(actual),
                      'mode': 'fresh reruns' if args.rerun else 'saved evidence validation'},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
