#!/usr/bin/env python3
"""Verify frozen bytes; this is a provenance check, not a replacement proof audit."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    count = 0
    for name, expected in manifest['files_sha256'].items():
        if digest((HERE / name).read_bytes()) != expected:
            raise RuntimeError(f'packet mismatch: {name}')
        count += 1
    for entry in manifest['inherited_proof_inputs']:
        # Always check the pinned object, not a possibly newer workspace version.
        raw = subprocess.check_output(
            ['git', 'show', f"{entry['commit']}:{entry['path']}"], cwd=ROOT)
        if digest(raw) != entry['sha256']:
            raise RuntimeError(f"inherited mismatch: {entry['path']}")
        count += 1
    print(json.dumps({'status': 'passed', 'hashed_files': count,
                      'scope': 'Packet bytes and pinned input Git objects only.'},
                     sort_keys=True))


if __name__ == '__main__':
    main()
