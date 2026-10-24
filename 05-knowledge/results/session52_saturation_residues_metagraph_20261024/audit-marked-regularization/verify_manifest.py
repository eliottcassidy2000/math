#!/usr/bin/env python3
"""Verify only manifest-declared local artifacts, rejecting unsafe paths."""
import hashlib
import json
from pathlib import Path

H=Path(__file__).resolve().parent

def verify(root):
    m=json.loads((root/'MANIFEST.json').read_text())
    names=[]
    for entry in m['files']:
        p=Path(entry['path'])
        if p.is_absolute() or '..' in p.parts or str(p) in names:
            raise ValueError('unsafe/duplicate manifest entry')
        names.append(str(p))
        path=root/p
        if path.is_symlink() or not path.is_file():
            raise ValueError('missing/nonregular file: '+str(p))
        data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=entry['sha256'] or len(data)!=entry['bytes']:
            raise ValueError('manifest mismatch: '+str(p))
    return len(names)

if __name__ == '__main__':
    print(json.dumps({'status':'PASS','files_verified':verify(H)},sort_keys=True))
