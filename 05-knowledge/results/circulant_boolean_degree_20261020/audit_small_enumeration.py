#!/usr/bin/env python3
"""Independent ordinary coefficient traversal, not the discovery Gray walk."""
from audit_basis import generate
from pathlib import Path
import json
rows=[]
for m in [2,3,4,6]:
 for d in [1,2]:
  reps,tgs,mons,pts,basis,can=generate(m,d);count=valid=0
  for c in range(1<<len(mons)):
   out=0
   for k,b in enumerate(basis):
    if c>>k&1:out^=b
   words=[out>>(m*j)&((1<<m)-1) for j in range(len(reps))]
   types=[can[y] for y in words]
   if min(types)>=0:
    valid+=1
    if len(set(types))==len(tgs):count+=1
  rows.append({'m':m,'degree_bound':d,'maps':1<<len(mons),'target_valued':valid,'type_bijections':count})
old=json.loads(Path(__file__).with_name('probe.json').read_text())
if len(rows)!=len(old) or any(any(r[k]!=o[k] for k in r) for r,o in zip(rows,old)):raise RuntimeError('Enumeration mismatch')
Path(__file__).with_name('small_enumeration.audit.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
