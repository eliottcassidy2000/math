from itertools import combinations
from math import isqrt
import json

def scan(N):
 seen={};triples=0
 for v in range(1,N+1):
  neighbors=[s*s-v for s in range(isqrt(v)+1,isqrt(v+N)+1) if s*s-v!=v]
  for triple in combinations(neighbors,3):
   triples+=1
   old=seen.setdefault(triple,[]);old.append(v)
   if len(old)>=3:return {'N':N,'contains_K33':True,'parts':[list(triple),old[:3]],'tested_neighbor_triples_before_witness':triples}
 return {'N':N,'contains_K33':False,'all_neighbor_triples':triples,'distinct_triples':len(seen)}
lo=scan(672);hi=scan(673)
if lo['contains_K33'] or not hi['contains_K33']:raise RuntimeError('threshold changed')
A,B=hi['parts']
for part in [A,B]:
 for u,v in combinations(part,2):
  if isqrt(u+v)**2==u+v:raise RuntimeError('unwanted edge')
for a in A:
 for b in B:
  if isqrt(a+b)**2!=a+b:raise RuntimeError('missing edge')
hi['induced']=True;hi['cross_roots']=[[isqrt(a+b)for b in B]for a in A]
print(json.dumps({'status':'PASS','lower_complete_scan':lo,'upper_witness':hi,'scope':'Finite exact threshold; absence of any K3,3 subgraph at672 and induced witness at673'},indent=2))
