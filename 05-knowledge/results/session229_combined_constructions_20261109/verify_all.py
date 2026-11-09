#!/usr/bin/env python3
"""Portable exact replay; no third-party packages, no assertions."""
import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[('recursive','verify.py','CHECK'),('recursive_audit','audit.py','CHECK'),('recursive_complete','verify.py','CHECK'),('recursive_complete','ternary.py','TERNARY'),('ternary_audit','audit.py','CHECK'),('overlap','verify.py','CHECK'),('overlap_audit','audit.py','CHECK'),('inverse','verify.py','CHECK'),('inverse_audit','check.py','CHECK')]
def main():
 rows=[]
 for directory,script,receipt in CASES:
  for optimized in (False,True):
   source=ROOT/directory/script
   expected=ROOT/directory/(receipt+('.optimized' if optimized else '')+'.json')
   args=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(source)]
   p=subprocess.run(args,cwd=source.parent,capture_output=True)
   if p.returncode:raise RuntimeError(f'{directory}/{script} failed: {p.stderr.decode()}')
   if p.stdout!=expected.read_bytes():raise RuntimeError(f'{directory}/{script} receipt mismatch, optimized={optimized}')
   rows.append({'script':str(source.relative_to(ROOT)),'optimized':optimized,'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'receipt_sha256':hashlib.sha256(p.stdout).hexdigest()})
 print(json.dumps({'status':'PASS','canonical_runs':len(rows),'runs':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
