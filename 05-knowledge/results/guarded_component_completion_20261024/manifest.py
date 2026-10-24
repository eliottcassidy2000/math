"""Verify packet files without treating byte order or line endings as proofs."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
TARGET = HERE/'MANIFEST.json'


def hashes():
    result = {}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p == TARGET or '__pycache__' in p.parts or p.suffix == '.pyc':
            continue
        raw = p.read_bytes()
        if p.suffix in {'.md', '.py', '.json', '.txt'}:
            raw = raw.replace(b'\r\n', b'\n')
        result[p.relative_to(HERE).as_posix()] = hashlib.sha256(raw).hexdigest()
    return result


actual = hashes()
if '--write' in sys.argv:
    TARGET.write_text(json.dumps(actual, indent=2, sort_keys=True)+'\n', encoding='utf-8')
elif actual != json.loads(TARGET.read_text(encoding='utf-8')):
    expected = json.loads(TARGET.read_text(encoding='utf-8'))
    changes = [s for s in sorted(set(expected)|set(actual)) if expected.get(s) != actual.get(s)]
    raise RuntimeError('Packet hash mismatch: '+', '.join(changes))
print(f'PASS {len(actual)} packet file hashes (text CRLF normalized)')
