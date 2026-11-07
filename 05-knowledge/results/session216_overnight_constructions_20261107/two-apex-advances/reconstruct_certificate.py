"""Verify ordered transport parts and reconstruct the original gzip."""
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent / 'certificate'
manifest = json.loads((root / 'MANIFEST.json').read_text())
if len(sys.argv) != 2:
    raise SystemExit('usage: python3 reconstruct_certificate.py OUTPUT_GZIP_PATH')
parts = manifest['parts']
expected_names = [f"coverage2327673.bin.gz.part{i:03d}" for i in range(len(parts))]
if [part['name'] for part in parts] != expected_names:
    raise RuntimeError('Nonconsecutive or incorrectly ordered part names')
if sorted(path.name for path in root.glob('*.part*')) != expected_names:
    raise RuntimeError('Missing or extra certificate part')
blob = bytearray()
for i, part in enumerate(parts):
    data = (root / part['name']).read_bytes()
    if len(data) != part['bytes'] or hashlib.sha256(data).hexdigest() != part['sha256']:
        raise RuntimeError('Part hash or size mismatch: ' + part['name'])
    if i < len(parts) - 1 and len(data) != manifest['part_size']:
        raise RuntimeError('Nonfinal part has incorrect fixed size')
    blob.extend(data)
if len(blob) != manifest['bytes'] or hashlib.sha256(blob).hexdigest() != manifest['sha256']:
    raise RuntimeError('Reconstructed gzip hash or size mismatch')
output = Path(sys.argv[1])
# Verify everything before writing. Refuse to overwrite an existing file.
with output.open('xb') as handle:
    handle.write(blob)
print(f"Reconstructed {len(blob)} bytes, SHA256 {manifest['sha256']}")
