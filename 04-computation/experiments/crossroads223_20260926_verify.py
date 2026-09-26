"""Reproduce the 223-session controls and hash their retained artifacts.

All selected programs use explicit validation gates. Each is run normally
and under -O; their complete UTF-8 outputs must agree before being saved.
This runner verifies reproduction, not the general mathematical proofs.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREFIX = 'crossroads223_20260926_'
CASES = [
    ('automata', []), ('flow', []), ('geometry', ['--rank']),
    ('robin', []), ('bridge', []), ('bridge_audit', []),
    ('sextic_audit', []), ('boundary_audit', []), ('flow_global', []),
    ('geometry_growth', []), ('flow_series', []), ('geometry_growth_audit', []),
]


def run_case(case):
    stem, args = case
    source = Path('04-computation/experiments') / (PREFIX + stem + '.py')
    output = Path('05-knowledge/results') / (PREFIX + stem + '.out')
    target = ROOT / source
    # Match the repository's committed LF hash basis on Windows as well.
    raw = target.read_bytes()
    normalized = raw.replace(b'\r\n', b'\n')
    if raw != normalized:
        target.write_bytes(normalized)
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    outputs = []
    for flags in ([], ['-O']):
        result = subprocess.run([sys.executable, *flags, str(source), *args],
                                cwd=ROOT, env=env, capture_output=True,
                                text=True, encoding='utf-8', timeout=240)
        if result.returncode:
            raise RuntimeError(f'{stem} {flags}: {result.stderr}\n{result.stdout}')
        outputs.append(result.stdout.replace('\r\n', '\n'))
    if outputs[0] != outputs[1]:
        raise RuntimeError(f'{stem}: normal and -O outputs differ')
    (ROOT / output).write_text(outputs[0], encoding='utf-8', newline='\n')
    print(f'PASS normal/-O: {stem}', flush=True)
    return {'source': source.as_posix(), 'args': args, 'output': output.as_posix(),
            'normal_equals_optimized': True}


def main():
    with ThreadPoolExecutor(max_workers=3) as pool:
        completed = [future.result() for future in
                     as_completed([pool.submit(run_case, case) for case in CASES])]
    files = []
    for directory in ('04-computation/experiments', '05-knowledge/results'):
        files.extend(p for p in (ROOT / directory).glob(PREFIX + '*')
                     if p.is_file() and p.suffix in ('.py', '.md', '.out'))
    for number in (4488, 4489, 4490, 4491):
        files.extend((ROOT / '01-canon/theorems').glob(f'THM-{number}-*.md'))
    artifacts = []
    for path in sorted(files):
        raw = path.read_bytes()
        artifacts.append({'path': path.relative_to(ROOT).as_posix(),
                          'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)})
    manifest = {
        'status': 'REPRODUCTION VERIFIED; mathematical scopes remain those of the proof notes',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'python': sys.version.split()[0], 'hash_basis': 'raw UTF-8 LF bytes',
        'cases': sorted(completed, key=lambda row: row['source']), 'artifacts': artifacts,
    }
    path = ROOT / '05-knowledge/results' / (PREFIX + 'manifest.json')
    path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f'PASS: {len(completed)} programs in both modes; {len(artifacts)} artifact hashes', flush=True)


if __name__ == '__main__':
    main()
