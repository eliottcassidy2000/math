#!/usr/bin/env python3
"""Verify the precommitted sampling recipe against its exact Git tree."""
import sys,json,subprocess,random,hashlib,pathlib
here=pathlib.Path(__file__).resolve().parent
repo=sys.argv[1] if len(sys.argv)>1 else '../session67-repo'
m=json.loads((here/'random_manifest.json').read_text())
paths=subprocess.check_output(['git','-C',repo,'ls-tree','-r','--name-only',m['pin']],text=True).splitlines()
rng=random.Random(m['seed'])
for s in m['strata']:
 u=sorted(p for p in paths if p.startswith(s['prefix']) and p.endswith(s['suffix']))
 if len(u)!=s['count'] or hashlib.sha256(('\n'.join(u)+'\n').encode()).hexdigest()!=s['universe_sha256']:raise RuntimeError('universe mismatch '+s['name'])
 if rng.sample(u,6)!=s['draws']:raise RuntimeError('draw mismatch '+s['name'])
 for p in s['draws']:
  b=subprocess.check_output(['git','-C',repo,'show',m['pin']+':'+p])
  if hashlib.sha256(b).hexdigest()!=s['draw_blob_sha256'][p]:raise RuntimeError('blob mismatch '+p)
print('PASS: 24 original stratified draws and blob hashes at '+m['pin'])
