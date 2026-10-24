#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
manifest=json.loads((p/'MANIFEST.json').read_text())
for name,expected in manifest['files'].items():
 b=(p/name).read_bytes()
 if name.endswith(('.json','.md','.py','.log')):b=b.replace(b'\r\n',b'\n')
 actual=hashlib.sha256(b).hexdigest()
 if actual!=expected:raise RuntimeError('Manifest mismatch: '+name)
print('PASS independent packet manifest:',len(manifest['files']),'files')
