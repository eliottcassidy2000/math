#!/usr/bin/env python3
import random,json,math,hashlib,pathlib
from check import tour,hp,substitute,triangle_row,check

def profile_subset(a):
 n=len(a); N=1<<n
 end=[[0]*n for _ in range(N)]; h=[0]*N
 for v in range(n): end[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   if not s>>v&1: continue
   prev=s^(1<<v)
   if prev: end[s][v]=sum(end[prev][u] for u in range(n) if prev>>u&1 and a[u][v])
  h[s]=sum(end[s])
 pc=[[0]*(n+1) for _ in range(N)];pc[0][0]=1
 for s in range(1,N):
  anchor=s&-s; rest=s^anchor; sub=rest
  while True:
   block=sub|anchor; rem=s^block
   for k in range(1,n+1): pc[s][k]+=h[block]*pc[rem][k-1]
   if not sub:break
   sub=(sub-1)&rest
 return [math.factorial(k)*pc[-1][k] for k in range(n+1)]

def main():
 rng=random.Random(330071); seen={}; n=7; result=None; masks=[]; profiles=set()
 for t in range(20000):
  mask=rng.randrange(1<<21); masks.append(mask); f=profile_subset(tour(n,mask)); profiles.add(tuple(f)); key=tuple(f[1:3])
  if key in seen and seen[key][1]!=f:
   old,oldf=seen[key]; rows=[triangle_row(j) for j in range(1,4)]
   result={'n':n,'masks':[old,mask],'profiles':[oldf,f],'samples_until_hit':t+1,'responses':[[sum(x*y for x,y in zip(row,g)) for row in rows] for g in [oldf,f]]}
   c3=[[0,1,0],[0,0,1],[1,0,0]]; T2=[[0,1],[0,0]]
   direct=[hp(substitute(c3,[tour(n,b),T2,T2])) for b in [old,mask]]
   check(direct==[r[1] for r in result['responses']],'direct eleven-vertex twin response')
   result['direct_11_vertex_responses']=direct;break
  seen[key]=(mask,f)
 if result is None:
  result={'status':'NO_WITNESS_IN_SAMPLED_UNIVERSE','n':n,'draws':len(masks),'distinct_masks':len(set(masks)),'distinct_profiles':len(profiles),'distinct_prefixes':len(seen),'sample_masks_sha256':hashlib.sha256(json.dumps(masks).encode()).hexdigest(),'seed':330071,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'profiles_observed':sorted(profiles),'scope':'Finite deterministic sample, not exhaustive classification or nonexistence theorem.'}
  print(json.dumps(result,indent=2,sort_keys=True));return
 # Independent permutation-and-cut check of just the final witnesses.
 from check import profile
 for b,f in zip(result['masks'],result['profiles']):check(profile(tour(n,b))==f,'independent permutation witness profile')
 result['status']='PASS'; result['seed']=330071;result['script_sha256']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
