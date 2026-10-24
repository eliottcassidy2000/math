#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
m=json.loads((P/'MANIFEST.json').read_text())
for r in m['files']:
 p=P/r['path']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError('manifest mismatch: '+r['path'])
print(json.dumps({'status':'PASS','verified_files':len(m['files'])}))
