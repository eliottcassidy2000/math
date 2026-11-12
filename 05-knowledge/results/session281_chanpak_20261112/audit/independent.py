#!/usr/bin/env python3
"""Restored independent root-split audit: no producer module imports."""
import json
from pathlib import Path
from functools import lru_cache
P=Path(__file__).resolve().parents[1]
def check(x,label):
 if not x:raise RuntimeError(label)
def profile(T):
 n=len(T);F=(1<<n)-1
 incoming=[sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n)]
 @lru_cache(None)
 def ending(S,v):
  if S==1<<v:return 1
  R=S^(1<<v);return sum(ending(R,u) for u in range(n) if R>>u&1 and T[u]>>v&1)
 @lru_cache(None)
 def starting(S,v):
  if S==1<<v:return 1
  R=S^(1<<v);return sum(starting(R,u) for u in range(n) if R>>u&1 and T[v]>>u&1)
 def H(S):return sum(ending(S,v) for v in range(n) if S>>v&1) if S else 1
 hs=[H(S) for S in range(F+1)];R=F^1
 Q=sum(hs[S]*hs[R^S] for S in range(F+1) if not S&1)
 U=V=d=e=0
 for S in range(1,F+1,2):
  u=starting(S,0)*hs[F^S];v=ending(S,0)*hs[F^S];sgn=1 if S.bit_count()%2 else -1
  U+=u;V+=v;d+=sgn*u;e+=sgn*v
 return dict(H=hs[F],h=hs[R],Q=Q,U=U,V=V,d=d,e=e,k=int(n==1))
def fusion(A,B):
 n=len(A);out=[0]*(n+len(B)-1)
 for T,ids in ((A,range(n)),(B,[0]+list(range(n,len(out))))):
  for i,row in enumerate(T):
   for j in range(len(T)):
    if row>>j&1:out[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,len(out)):out[i]|=1<<j
 return out
r=json.loads((P/'receiver/RESULTS.json').read_text())
raw=r['raw_rows'];done=[]
for i,T in enumerate(raw):
 p=profile(tuple(T));q=profile(tuple(fusion((2,4,1),T)))
 check((p['H'],p['h'],p['V'],p['e'])==(979+2*i,15,1203,13 if i==3 else 15),'raw receiver')
 check((q['H'],q['h'],q['V'],q['e'])==(3385+2*i,15,4797,15),'prepared receiver')
 done.append(dict(raw=p,prepared=q))
A=(0,41,3,5,15,29);B=(0,33,3,7,15,29)
check(profile(A)==profile(B),'full fusion collision')
for T,H in [((96,65,3,7,79,94,12),65),((102,8,66,69,79,30,34),69)]:
 # Given witnesses use root6; put it first.
 ids=[6,0,1,2,3,4,5]
 T=tuple(sum(1<<j for j in range(7) if T[ids[i]]>>ids[j]&1) for i in range(7))
 p=profile(T);check((p['H'],p['h'],p['Q'],p['U'],p['V'],p['d'],p['e'])==(H,17,224,89,152,11,6),'coprime central pair')
print(json.dumps(dict(status='PASS',method='independent recursive root-split counting; no producer imports',receiver=done,central_pair_h17=True,full_carrier_alias=True),indent=2,sort_keys=True))
