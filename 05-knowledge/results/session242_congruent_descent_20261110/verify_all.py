#!/usr/bin/env python3
"""Replay every canonical Session242 producer and independent audit, twice."""
import hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[('arithmetic/verify.py','arithmetic/verify.json'),('descent/verify.py','descent/CHECK.json'),('exception/verify.py','exception/CHECK.json'),('audit/verify.py','audit/verify.json'),('graft_audit/verify.py','graft_audit/CHECK.json')]
records=[]
for script,receipt in CASES:
 expected=(ROOT/receipt).read_bytes()
 for optimized in (False,True):
  cmd=[sys.executable,'-I','-B']+(['-O'] if optimized else[])+[str(ROOT/script)]
  r=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=False)
  if r.returncode or r.stdout!=expected:
   raise RuntimeError(f'Replay mismatch: {script} optimized={optimized}: {r.stderr.decode()}')
  records.append({'script':script,'optimized':optimized,'sha256':hashlib.sha256(r.stdout).hexdigest(),'status':'PASS'})
negative=[]
for script in ['exception/verify.py','graft_audit/verify.py']:
 for optimized in (False,True):
  cmd=[sys.executable,'-I','-B']+(['-O'] if optimized else[])+[str(ROOT/script),'--fail']
  r=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=False)
  if r.returncode==0 or b'deliberate' not in r.stderr:
   raise RuntimeError(f'Negative gate did not fail intentionally: {script}')
  negative.append({'script':script,'optimized':optimized,'status':'EXPECTED_FAILURE'})
print(json.dumps({'status':'PASS','canonical_replays':records,'negative_controls':negative},indent=2,sort_keys=True))
