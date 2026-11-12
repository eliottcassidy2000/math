#!/usr/bin/env python3
"""Replay every session engine, optimized engine and committed receipt."""
import hashlib,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
TASKS=[('reset/verify.py','reset/RESULTS.json'),('symmetry/verify.py','symmetry/REPLAY.json'),('wildcard/verify.py','wildcard/RESULTS.json'),('audit/independent.py','audit/independent.out')]
def main():
 rows=[]
 for script,receipt in TASKS:
  results=[]
  for options in ([],['-O']):
   p=subprocess.run([sys.executable,*options,str(BASE/script)],cwd=BASE,capture_output=True,text=True,check=True)
   results.append(json.loads(p.stdout))
  stored=json.loads((BASE/receipt).read_text())
  if not results[0]==results[1]==stored:raise RuntimeError('receipt mismatch: '+script)
  rows.append({'script':script,'normal_optimized_saved_agree':True,'sha256':hashlib.sha256((BASE/script).read_bytes()).hexdigest()})
  print('PASS '+script,file=sys.stderr,flush=True)
 print(json.dumps({'status':'PASS','tasks':rows,'scope':'Exact controls for the specified families and proofs, not Chan-Pak completeness'},indent=2,sort_keys=True))
if __name__=='__main__':main()
