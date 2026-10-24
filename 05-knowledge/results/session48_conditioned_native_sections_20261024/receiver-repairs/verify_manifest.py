#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import argparse,json
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
files={p.relative_to(P).as_posix():sha256(p.read_bytes()).hexdigest() for p in sorted(P.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
target=P/'MANIFEST.json'
if a.write:target.write_text(json.dumps(files,sort_keys=True,indent=2)+'\n')
else:
 expected=json.loads(target.read_text())
 if expected!=files:raise RuntimeError('Packet hash mismatch: '+', '.join(k for k in sorted(set(files)|set(expected)) if files.get(k)!=expected.get(k)))
pins=json.loads((P/'SOURCE_PINS.json').read_text())
for row in pins['files']:
 data=(P/row['local_path']).read_bytes()
 if sha256(data).hexdigest()!=row['sha256'] or len(data)!=row['bytes']:raise RuntimeError('Original input pin mismatch: '+row['local_path'])
print(json.dumps({'status':'PASS','packet_files':len(files),'original_source_files':len(pins['files']),'source_commit':pins['source_commit']},sort_keys=True))
