#!/usr/bin/env python3
"""Reproduce this scoped constructive packet; no inherited census replay."""
import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
jobs=[('compiler.py',[],'results.json'),('explore.py',['--bound','262143'],'explore-results.json'),('repair_holes.py',[],'hole-repairs.json')]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
normal={}
for optimized in (False,True):
 for script,args,receipt in jobs:
  cmd=[sys.executable]+(['-O'] if optimized else [])+['-B',str(root/script),*args]
  subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
  value=digest(root/receipt)
  if not optimized:normal[receipt]=value
  elif normal[receipt]!=value:raise RuntimeError('Normal/optimized mismatch: '+receipt)
report={'normal_optimized_identical':True,'receipt_sha256':normal,'scripts_sha256':{s:digest(root/s) for s,_,_ in jobs}}
(root/'replay.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(json.dumps(report,sort_keys=True,indent=2))
