"""Verify the complete packet with portable normalized-text hashes."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
TARGET = HERE/'MANIFEST.json'
actual = {}
for path in sorted(HERE.rglob('*')):
    if not path.is_file() or path == TARGET or '__pycache__' in path.parts or path.suffix == '.pyc':
        continue
    raw = path.read_bytes()
    if path.suffix in {'.md','.py','.json','.txt'}:
        raw = raw.replace(b'\r\n', b'\n')
    actual[path.relative_to(HERE).as_posix()] = hashlib.sha256(raw).hexdigest()

if '--write' in sys.argv:
    TARGET.write_text(json.dumps(actual, indent=2, sort_keys=True)+'\n', encoding='utf-8')
else:
    expected = json.loads(TARGET.read_text(encoding='utf-8'))
    if actual != expected:
        changed = [p for p in sorted(set(actual)|set(expected)) if actual.get(p) != expected.get(p)]
        raise RuntimeError('Packet hash mismatch: '+', '.join(changed))
print(f'PASS {len(actual)} packet hashes; text CRLF normalized')
