#!/usr/bin/env python3
"""Run exact standalone proof controls twice and compare stored receipts."""
import hashlib,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
TASKS=[('automata/verify.py','automata/replay.json'),('additive/verify.py','additive/RESULTS.json'),('additive/mixed.py','additive/MIXED.json'),('additive/padded_bank.py','additive/PADDED_BANK.json'),('reset/verify.py','reset/RESULTS.json'),('wildcard/verify.py','wildcard/replay.json'),('audit/independent.py','audit/independent.json'),('audit/reset_independent.py','audit/reset_independent.json'),('audit/homometry_independent.py','audit/homometry_independent.json'),('audit/decoder_independent.py','audit/decoder_independent.json'),('audit/padded_independent.py','audit/padded_independent.json')]
def main():
 rows=[]
 for script,receipt in TASKS:
  outputs=[]
  for flags in ([],['-O']):
   p=subprocess.run([sys.executable,*flags,str(BASE/script)],cwd=BASE,capture_output=True,text=True,check=True)
   outputs.append(json.loads(p.stdout))
  saved=json.loads((BASE/receipt).read_text())
  if not outputs[0]==outputs[1]==saved:raise RuntimeError('receipt mismatch: '+script)
  rows.append({'script':script,'normal_optimized_and_saved_agree':True,'source_sha256':hashlib.sha256((BASE/script).read_bytes()).hexdigest()})
  print('PASS '+script,file=sys.stderr,flush=True)
 print(json.dumps({'status':'PASS','tasks':rows,'scope':'Exact finite controls of stated proofs; not a Chan–Pak completeness certificate'},indent=2,sort_keys=True))
if __name__=='__main__':main()
