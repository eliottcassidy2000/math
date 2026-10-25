"""Verify the complete packet file set and SHA256 bytes, including this tool."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'MANIFEST.json').read_text(encoding='utf-8'))
files={p.relative_to(HERE).as_posix():p for p in HERE.rglob('*')
       if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
expected={r['path']:r for r in manifest['files']}
if len(expected)!=len(manifest['files']) or files.keys()!=expected.keys():
    raise RuntimeError(f'file-set mismatch missing={sorted(expected.keys()-files.keys())}, '
                       f'extra={sorted(files.keys()-expected.keys())}')
for name,p in files.items():
    data=p.read_bytes()
    if len(data)!=expected[name]['bytes'] or hashlib.sha256(data).hexdigest()!=expected[name]['sha256']:
        raise RuntimeError(f'byte mismatch: {name}')
print(f'PASS: {len(files)} packet files; complete file set and exact bytes')
