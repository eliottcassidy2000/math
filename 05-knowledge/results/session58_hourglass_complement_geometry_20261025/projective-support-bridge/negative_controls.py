#!/usr/bin/env python3
"""Prove representative false formulas fail in normal and optimized modes."""
from pathlib import Path
import subprocess, sys, tempfile, json
source=Path(__file__).with_name('verify_support.py').read_text()
mutations={
 'wrong_support_count':('hist[S]==2**(r-1)','hist[S]==2**r'),
 'wrong_gcd':('gcd(e*e,d)==t*t','gcd(e*e,d)==t'),
 'wrong_torus_count':('hist[S]==(q-1)**(r-1)','hist[S]==q**(r-1)'),
}
results=[]
with tempfile.TemporaryDirectory() as td:
 for name,(old,new) in mutations.items():
  if source.count(old)!=1:raise RuntimeError('Mutation target not unique: '+name)
  path=Path(td)/(name+'.py');path.write_text(source.replace(old,new))
  for mode in ([],['-O']):
   r=subprocess.run([sys.executable,*mode,str(path)],capture_output=True,text=True)
   if r.returncode==0 or 'RuntimeError' not in r.stderr:raise RuntimeError('False formula accepted: '+name)
   results.append({'mutation':name,'mode':'optimized' if mode else 'normal','outcome':'rejected','last_error_line':r.stderr.strip().splitlines()[-1]})
print(json.dumps(results,indent=2))
