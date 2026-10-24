#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import json

ROOT = Path(__file__).resolve().parent


def main():
    manifest = json.loads((ROOT/'MANIFEST.json').read_text())
    for row in manifest['files']:
        p = ROOT/row['path']
        data = p.read_bytes()
        if len(data) != row['bytes'] or sha256(data).hexdigest() != row['sha256']:
            raise RuntimeError('manifest mismatch: '+row['path'])
    pins = json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for row in pins:
        data = (ROOT/row['archive_path']).read_bytes()
        if len(data) != row['bytes'] or sha256(data).hexdigest() != row['sha256']:
            raise RuntimeError('source mismatch: '+row['archive_path'])
    expected = ''.join(f"{r['sha256']}  {r['bytes']}  {r['archive_path']}\n" for r in pins)
    if (ROOT/'bytepinned.txt').read_text() != expected:
        raise RuntimeError('source pin table mismatch')
    print('PASS: packet manifest and byte-pinned sources')


if __name__ == '__main__':
    main()
