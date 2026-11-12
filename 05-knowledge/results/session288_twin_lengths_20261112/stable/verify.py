from pathlib import Path
from itertools import permutations,combinations,product
from tensor import *
from cones import run as cones_run,expand
from normalforms import run as normal_run,need,canon
from explore import val
import json
HERE=Path(__file__).resolve().parent

def brute(rows,O):
 H=W=0;n=len(rows)
 for p in permutations(range(n)):
  ok=[bool(rows[p[i]]>>p[i+1]&1) for i in range(n-1)]
  H+=all(ok)
  for j in range(n-1):
   if O>>p[j]&1 and O>>p[j+1]&1 and all(ok[:j]+ok[j+1:]):W+=1
 return H-W

def run():
 formal={(2,3):1,(2,1):-3,(1,3):-2,(1,1):7}
 A=anchored(formal);need(A=={(0,0):1,(2,3):1},'anchored formal interaction')
 need(stable_faces(A,2)==[0,1,2],'formal minimal forbidden face')
 r=[14,28,8,16,37,79,31];O=3;ports=[2,6]
 K=kernel(r,O,ports);need(K=={(1,1):1},'two-clock exact tensor')
 controls=0;keys={}
 for a,b in product(range(1,4),repeat=2):
  need(val(r,O,{2:a,6:b})==1,'two-clock literal delta');controls+=1
  rr,oo,typ=expand(r,O,[a if v==2 else b if v==6 else 1 for v in range(7)])
  key=canon(rr,oo);need(key not in keys,'faithful two-clock lengths');keys[key]=(a,b)
 for rows,sel in [([2,4,1],3),(r,O),([6,4,0],1)]:
  need(brute(rows,sel)==val(rows,sel),'independent permutation split count');controls+=1
 hostile=[90,28,9,16,36,79,30];hO=72
 hK=kernel(hostile,hO,[2]);need(hK=={(1,):1,(3,):-4},'screen hostile full tensor')
 need(val(hostile,hO,{2:2})==1 and val(hostile,hO,{2:3})==-23,'hostile literal');controls+=2
 try:kernel([2,4,1],0,[0,1,2])
 except ValueError:pass
 else:raise RuntimeError('cyclic port accepted')
 cc=cones_run();nn=normal_run()
 return {'literal_and_independent_controls':controls,'formal_hostile_not_claimed_tournament':True,'all_cone_tensors_verified':len(cc['full_tensor_certificates']),'normalforms':nn['states'],'cone_sequence':cc['counts_orders1_to30'],'eventual_linear_formula':'19*n-49 for n>=8','scope':'Exact fixed-bank cylinder union; no global stable-simplex or adaptive closure claim.'}
if __name__=='__main__':
 out=run();print(json.dumps(out,sort_keys=True,indent=2));(HERE/('replay.optimized.json' if not __debug__ else 'replay.json')).write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
