#!/usr/bin/env python3
"""Negative controls for the packet integrity gate, normal and optimized.

Every mutation occurs in a private temporary copy. No native file is changed.
"""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    manifest = json.loads((ROOT / 'FILE_MANIFEST.json').read_text())
    script = next(row['path'] for row in manifest['files']
                  if row['path'].endswith('.py') and row['path'] != 'verify_packet.py')
    receipt = manifest['replays'][0]['normal']
    cases = ('source-byte', 'receipt-byte', 'missing-file', 'extra-file',
             'symlink-file', 'unsafe-manifest-path')
    results = []
    with tempfile.TemporaryDirectory(prefix='session69-mutations-') as temp:
        for case in ('unchanged', *cases):
            packet = Path(temp) / case
            shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))
            if case == 'source-byte':
                p = packet / script
                p.write_bytes(p.read_bytes() + b'\n# mutation control\n')
            elif case == 'receipt-byte':
                p = packet / receipt
                data = p.read_bytes()
                p.write_bytes(bytes([data[0] ^ 1]) + data[1:])
            elif case == 'missing-file':
                (packet / receipt).unlink()
            elif case == 'extra-file':
                (packet / 'UNLISTED_CONTROL.txt').write_text('unlisted\n')
            elif case == 'symlink-file':
                p = packet / receipt
                p.unlink()
                p.symlink_to(ROOT / receipt)
            elif case == 'unsafe-manifest-path':
                altered = json.loads((packet / 'FILE_MANIFEST.json').read_text())
                altered['files'][0]['path'] = '../OUTSIDE_CONTROL'
                altered['files'].sort(key=lambda row: row['path'])
                (packet / 'FILE_MANIFEST.json').write_text(json.dumps(altered))
            for optimized in (False, True):
                cmd = [sys.executable, '-I', '-B'] + (['-O'] if optimized else [])
                proc = subprocess.run(cmd + [str(packet / 'verify_packet.py')],
                                      cwd=packet, capture_output=True, check=False)
                require((proc.returncode == 0) == (case == 'unchanged'),
                        'unexpected verification result: ' + case)
                results.append({'case': case, 'optimized': optimized,
                                'outcome': 'accepted' if proc.returncode == 0 else 'rejected'})
    print(json.dumps({'status': 'PASS', 'mutation_cases': len(cases),
                      'normal_and_optimized_checks': len(results), 'results': results},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
