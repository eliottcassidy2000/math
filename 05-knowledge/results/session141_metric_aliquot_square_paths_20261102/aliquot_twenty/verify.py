#!/usr/bin/env python3
"""Closed, relocation-safe verifier; no original repository required."""
import argparse,hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
args=argparse.ArgumentParser();args.add_argument('--replay',action='store_true');a=args.parse_args()
def need(condition,message):
    if not condition:raise RuntimeError(message)
man=json.loads((root/'MANIFEST.json').read_text())
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
need(actual==set(man['files'])|{'MANIFEST.json'},'closed manifest member mismatch')
for name,meta in man['files'].items():
    b=(root/name).read_bytes()
    need(len(b)==meta['bytes'] and hashlib.sha256(b).hexdigest()==meta['sha256'],'manifest hash: '+name)
source=json.loads((root/'sourcepins/index.json').read_text())
for entry in source['sources']:
    b=(root/'sourcepins'/entry['snapshot']).read_bytes()
    need(hashlib.sha256(b).hexdigest()==entry['sha256'],'source sha256: '+entry['original_path'])
    git_blob=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    need(git_blob==entry['git_blob'],'source git blob: '+entry['original_path'])
need((root/'normal.json').read_bytes()==(root/'optimized.json').read_bytes(),'normal/optimized mismatch')
if a.replay:
    for optimized in [False,True]:
        cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(root/'check.py')]
        result=subprocess.run(cmd,capture_output=True,check=True)
        need(result.stdout==(root/'normal.json').read_bytes(),'replay mismatch')
        need(not result.stderr,'replay stderr')
print(json.dumps({'status':'PASS','manifest_members':len(man['files']),
 'source_pins':len(source['sources']),'source_commit':source['commit'],
 'replays':2 if a.replay else 0},sort_keys=True))
