#!/usr/bin/env python3
"""Portable replay, including copied scripts executed from an unrelated cwd."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
rows=[]
with tempfile.TemporaryDirectory(prefix='relocation-',dir=ROOT) as temp:
 base=Path(temp);moved=base/'moved';cwd=base/'unrelated-cwd'
 moved.mkdir();cwd.mkdir()
 for script,normal,optimized in [('verify.py','normal.json','optimized.json'),('independent.py','independent.normal.json','independent.optimized.json')]:
  (moved/script).write_bytes((ROOT/script).read_bytes())
  outputs=[]
  for flags,expected in [([],normal),(['-O'],optimized)]:
   run=subprocess.run([sys.executable,'-B',*flags,str(moved/script)],cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
   if run.stderr:raise RuntimeError(run.stderr.decode())
   if run.stdout!=(ROOT/expected).read_bytes():raise RuntimeError('relocation differs from frozen '+expected)
   outputs.append(run.stdout)
  if outputs[0]!=outputs[1]:raise RuntimeError('normal/optimized mismatch')
  result=json.loads(outputs[0]);rows.append({'script':script,'checks':result['checks'],'stdout_sha256':sha(outputs[0]),'row_digest':result['row_digest'],'normal_equals_optimized':True,'relocated_matches_frozen':True})
if rows[0]['row_digest']!=rows[1]['row_digest']:raise RuntimeError('independent row digest mismatch')
print(json.dumps({'status':'PASS','relocation':'Copied each script to a fresh directory; ran from a separate unrelated working directory','replays':rows},sort_keys=True,indent=2))
