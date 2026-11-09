#!/usr/bin/env python3
from check import adj,attach,paths,cut,mod,require
import json
from pathlib import Path

def iso(a,b):
 n=len(a); da=[x.bit_count() for x in a];db=[x.bit_count() for x in b]
 if sorted(da)!=sorted(db):return False
 candidates=[[v for v in range(n) if da[u]==db[v]] for u in range(n)]
 order=sorted(range(n),key=lambda u:len(candidates[u]));mapping={};used=set()
 def rec(i):
  if i==n:return True
  u=order[i]
  for v in candidates[u]:
   if v in used:continue
   if any(((a[u]>>z)&1)!=((b[v]>>w)&1) for z,w in mapping.items()):continue
   mapping[u]=v;used.add(v)
   if rec(i+1):return True
   del mapping[u];used.remove(v)
  return False
 return rec(0)

reps=[];pairs=0
for line in Path(__file__).with_name('witnesses_source.txt').read_text().splitlines():
 if not line.startswith('W '):continue
 _,c,x,h,status=line.split();a=attach(adj(10,int(c)),int(x));h=int(h)
 require(sum(paths(a)[-1])==h, f"Hamiltonian count mismatch: {line}")
 require(cut(a)>=3, f"Connectivity mismatch: {line}")
 if status=='inS10':require(mod(a), f"Missing nontrivial module: {line}")
 pairs+=1
 for r in reps:
  if r['h']==h and iso(a,r['a']):r['multiplicity']+=1;break
 else:reps.append({'h':h,'status':status,'a':a,'multiplicity':1,'has_module':mod(a)})
require(len(reps)==11, f"Expected 11 classes, found {len(reps)}")
require(sum(r["status"]=="inS10" for r in reps)==7, "Expected 7 old-valued classes")
print(json.dumps({'pairs':pairs,'classes':reps,'scope':'Direct independent checks of supplied witnesses only, not exhaustive completeness of order-11 input universe or membership in S10.'},indent=2))
