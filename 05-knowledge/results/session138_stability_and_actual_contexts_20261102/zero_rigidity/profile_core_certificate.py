"""Finite local certificates for the all-order ordinal profile dichotomy.

Only six six-vertex and nine nine-vertex cores are load-bearing.
The path-cover subset DP is independently checked on all six-vertex cores
by direct permutation/descent enumeration.
"""
import sys,itertools,json,hashlib,math
from pathcover_probe import profile,induce
from reference_paths import c3_blowup,extend

def tt(k):return tuple(sum(1<<j for j in range(i+1,k)) for i in range(k))
def make(k,r,t):
 T=c3_blowup(tt(k));vs=[i for i in range(3*k) if i!=3*r];B=induce(T,vs);S=0
 for i,n in enumerate(t):
  if i==r:
   for j in range(2):
    if n>>j&1:S|=1<<vs.index(3*i+j+1)
  else:
   for j in range(n):S|=1<<vs.index(3*i+j)
 return T,extend(B,S)
def literal_profile(a):
 n=len(a);f=[0]*n
 for p in itertools.permutations(range(n)):
  d=sum(not(a[p[i]]>>p[i+1]&1) for i in range(n-1))
  for r in range(d+1,n+1):f[r-1]+=math.comb(n-1-d,r-1-d)
 return f
rows=[]
for k,r,t in [(2,0,(2,t)) for t in range(3)]+[(2,1,(t,0)) for t in range(1,4)]+[(3,0,(0,s,t)) for s in range(1,4) for t in range(3)]:
 T,U=make(k,r,t);f=profile(T);ff=profile(U)
 if k==2 and (f!=literal_profile(T) or ff!=literal_profile(U)):raise RuntimeError('literal cover audit')
 d=[b-a for a,b in zip(f,ff)]
 if min(d)<0 or not any(d):raise RuntimeError(('core domination fails',k,r,t,d))
 row={'k':k,'root':r,'cut_profile':t,'core_masks':U,'difference':d};rows.append(row);print(json.dumps(row,sort_keys=True))
print('PASS '+json.dumps({'six_vertex_cores':6,'nine_vertex_cores':9,'all_entrywise_nonnegative':True,'sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()},sort_keys=True))
