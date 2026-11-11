"""Replay the three compact Session252 engines and compare their retained outputs."""
from pathlib import Path
import subprocess,sys,hashlib,json
root=Path(__file__).resolve().parent
jobs=[('algebra/probe.py','algebra/results.json',False),('golden/probe.py','golden/results.txt',False),('arithmetic/derive.py','arithmetic/results.json',True)]
results=[]
for script,target,writes in jobs:
 expected=(root/target).read_bytes();outputs=[]
 for flags in [[],['-O']]:
  r=subprocess.run([sys.executable,*flags,str(root/script)],capture_output=True,check=True)
  produced=(root/target).read_bytes() if writes else r.stdout
  if produced!=expected:raise AssertionError(('artifact mismatch',script,flags))
  outputs.append(r.stdout)
 if outputs[0]!=outputs[1]:raise AssertionError(('optimized mismatch',script))
 results.append({'script':script,'ordinary':'PASS','optimized':'PASS','artifact_sha256':hashlib.sha256(expected).hexdigest()})
print(json.dumps({'scope':'new compact identities, not an all-height coverage test','engines':results},indent=2))
