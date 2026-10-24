#!/usr/bin/env python3
"""Manifest-only standalone relocation and independent replay."""
import hashlib,json,subprocess,sys,tempfile,shutil
from pathlib import Path
P=Path(__file__).resolve().parent;m=json.loads((P/'MANIFEST.json').read_text());results=[]
with tempfile.TemporaryDirectory(prefix='blue-wiring-audit-') as td:
 q=Path(td)
 for r in m['files']:
  source=P/r['path'];target=q/r['path']
  if hashlib.sha256(source.read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError('source hash '+r['path'])
  target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 shutil.copy2(P/'MANIFEST.json',q/'MANIFEST.json')
 for optimized in [False,True]:
  prefix=[sys.executable]+(['-O'] if optimized else [])
  r=subprocess.run(prefix+[str(q/'audit.py')],cwd=q,capture_output=True,text=True,check=True)
  d=json.loads(r.stdout)
  if d!=json.loads((P/'independent-normal.json').read_text()):raise RuntimeError('relocated independent receipt mismatch')
  t=subprocess.run(prefix+[str(q/'verify_topology.py')],cwd=q,capture_output=True,text=True,check=True)
  distance=subprocess.run(prefix+[str(q/'distances.py')],cwd=q,capture_output=True,text=True,check=True)
  if json.loads(distance.stdout)!=json.loads((P/'distances-normal.json').read_text()):raise RuntimeError('distance replay mismatch')
  results.append({'distance_status':json.loads(distance.stdout)['status'],'optimized':optimized,'audit_status':d['status'],'topology_status':json.loads(t.stdout)['status']})
 subprocess.run([sys.executable,str(q/'verify_manifest.py')],cwd=q,capture_output=True,text=True,check=True)
print(json.dumps({'status':'PASS','manifest_only_files':len(m['files']),'outside_dependencies':['Python3','NumPy'],'runs':results},indent=2))
