#!/usr/bin/env python3
"""Verify the entire frozen no-Git audit packet, independent of optimization."""
from pathlib import Path
import hashlib,json
H=Path(__file__).resolve().parent
M=json.loads((H/'MANIFEST.json').read_text())
listed={r['path'] for r in M['files']}
pins=json.loads((H/'SOURCE_PINS.json').read_text())['files']
for r in pins:
 if r['snapshot'] not in listed:raise RuntimeError('Source-pin snapshot missing from delivery manifest '+r['snapshot'])
for r in M['files']:
 p=H/r['path'];b=p.read_bytes()
 if len(b)!=r['bytes'] or hashlib.sha256(b).hexdigest()!=r['sha256']:raise RuntimeError('Manifest mismatch '+r['path'])
for stem in ['','producer-native.','producer-audit.','guard-completeness.','odd-part.','helper-raw.']:
 a=H/(stem+'normal.json');b=H/(stem+'optimized.json')
 if a.read_bytes()!=b.read_bytes():raise RuntimeError('Normal/-O mismatch '+stem)
print(json.dumps({'status':'PASS','manifest_files':len(M['files']),'normal_optimized_pairs':6},sort_keys=True))
