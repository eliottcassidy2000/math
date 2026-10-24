"""Validate a closed audit packet. No external dependencies or disabled assertions."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import sys

def validate(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'MANIFEST.json').read_text())
    if manifest.get('schema') != 's7-independent-audit-manifest-v1':
        raise ValueError('unknown schema')
    seen = set()
    for row in manifest['files']:
        name = row['path']
        path = PurePosixPath(name)
        if not isinstance(name, str) or path.is_absolute() or '..' in path.parts or str(path) != name or name in seen or name == 'MANIFEST.json':
            raise ValueError('unsafe or duplicate path: ' + str(name))
        seen.add(name)
        target = root / name
        if target.is_symlink() or not target.is_file() or root not in target.resolve().parents:
            raise ValueError('missing or unsafe file: ' + name)
        data = target.read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError('content mismatch: ' + name)
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() or p.is_symlink()}
    if actual != seen | {'MANIFEST.json'}:
        raise ValueError('undeclared payload: ' + str(sorted(actual - seen - {'MANIFEST.json'})))
    return len(seen)

if __name__ == '__main__':
    try:
        count = validate(Path(__file__).parent)
        print('PASS: manifest closure, paths, sizes, SHA256; %d payload files' % count)
    except Exception as exc:
        print('REJECT: ' + str(exc), file=sys.stderr)
        sys.exit(1)
