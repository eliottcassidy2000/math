#!/usr/bin/env python3
"""Actually copy only manifest-listed bytes and replay in an unrelated cwd."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def execute(script, flags, cwd, args=()):
    p = subprocess.run([sys.executable, *flags, str(script), *map(str, args)], cwd=cwd, capture_output=True, text=True, check=True)
    return {'stdout': p.stdout, 'stderr': p.stderr, 'exit_code': p.returncode}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='session49-independent-repair-audit-') as raw:
        root = Path(raw); clone = root / 'packet'; clone.mkdir(); cwd = root / 'unrelated'; cwd.mkdir()
        for name in ['MANIFEST.json', *manifest['files']]:
            destination = clone / name; destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(HERE / name, destination)
        result = {'schema': 'session49-manifest-only-relocated-theory-audit-v1', 'copied_manifest_listed_files': len(manifest['files'])}
        result['initial_manifest'] = execute(clone / 'verify_manifest.py', ['-I', '-B'], cwd)
        for label, flags in [('normal', ['-I', '-B']), ('optimized', ['-I', '-O', '-B'])]:
            receipt = root / (label + '.json')
            result[label] = execute(clone / 'independent_check.py', flags, cwd, ['--output', receipt])
            if receipt.read_bytes() != (clone / (label + '.json')).read_bytes():
                raise RuntimeError('Relocated independent receipt mismatch: ' + label)
            result[label]['receipt_sha256'] = sha256(receipt.read_bytes()).hexdigest()
        if result['normal'] != result['optimized']:
            raise RuntimeError('Relocated interpreter modes differ')
        result['producer_replay'] = execute(clone / 'replay_producers.py', ['-I', '-B'], cwd)
        result['final_manifest'] = execute(clone / 'verify_manifest.py', ['-I', '-B'], cwd)
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.write:
        (HERE / 'PORTABILITY.json').write_text(payload)
    elif (HERE / 'PORTABILITY.json').read_text() != payload:
        raise RuntimeError('Relocated audit differs from frozen portability receipt')
    print('PASS actual manifest-only copy: independent and all producer normal/optimized replays')


if __name__ == '__main__':
    main()
