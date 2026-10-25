"""Verify complete packet file set and exact bytes against MANIFEST.json."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
manifest = json.loads((HERE/'MANIFEST.json').read_text(encoding='utf-8'))
files = {p.relative_to(HERE).as_posix(): p for p in HERE.rglob('*')
         if p.is_file() and p.name != 'MANIFEST.json' and '__pycache__' not in p.parts}
expected = {row['path']: row for row in manifest['files']}
if files.keys() != expected.keys():
    raise RuntimeError(f'file set mismatch: missing={sorted(expected.keys()-files.keys())}, extra={sorted(files.keys()-expected.keys())}')
for name, path in files.items():
    data = path.read_bytes()
    if len(data) != expected[name]['bytes'] or hashlib.sha256(data).hexdigest() != expected[name]['sha256']:
        raise RuntimeError(f'byte mismatch: {name}')
print(f'PASS: {len(files)} packet files, exact bytes and file set')
