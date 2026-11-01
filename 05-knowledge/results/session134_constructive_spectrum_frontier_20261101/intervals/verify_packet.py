#!/usr/bin/env python3
"""Verify the frozen packet, with no imports from the producer."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'MANIFEST.json').read_text())
for name, expected in manifest['sha256'].items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError('checksum mismatch: ' + name)
if (root / 'normal.json').read_bytes() != (root / 'optimized.json').read_bytes():
    raise RuntimeError('normal/optimized receipt mismatch')
print(json.dumps({'status': 'passed', 'sealed_files': len(manifest['sha256'])}, sort_keys=True))
