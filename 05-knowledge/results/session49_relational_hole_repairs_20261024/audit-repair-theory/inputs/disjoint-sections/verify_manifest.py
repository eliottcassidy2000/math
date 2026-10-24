#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent

def verify(root=HERE):
    m=json.loads((root/'MANIFEST.json').read_text())
    for row in m['files']:
        p=root/row['path']
        if not p.is_file():raise RuntimeError('missing '+row['path'])
        b=p.read_bytes()
        if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise RuntimeError('mismatch '+row['path'])
    return {'status':'PASS','files':len(m['files'])}
if __name__=='__main__':print(json.dumps(verify(),sort_keys=True))
