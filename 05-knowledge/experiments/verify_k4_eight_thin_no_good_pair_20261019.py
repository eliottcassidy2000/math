#!/usr/bin/env python3
"""Exact one-poset control; no poset census. Explicit guards survive python -O."""
from functools import lru_cache
from fractions import Fraction
import json

def require(test,message):
 if not test:raise RuntimeError(message)
n=12;full=(1<<n)-1
missing=[{0,2,4},{0,3,5},{1,2,5},{1,3,4},{0,1},{2,3}]
p=[0]*6+[sum(1<<i for i,r in enumerate(missing) if j not in r) for j in range(6)]
U=[sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n)]
@lru_cache(None)
def suffix(mask):
 if mask==full:return 1
 return sum(suffix(mask|1<<i) for i in range(n) if not mask>>i&1 and not p[i]&~mask)
E=suffix(0);forward={0:1};N=[[0]*n for _ in range(n)]
for mask in range(full+1):
 f=forward.get(mask,0)
 if not f:continue
 for b in range(n):
  if mask>>b&1 or p[b]&~mask:continue
  nxt=mask|1<<b;forward[nxt]=forward.get(nxt,0)+f
  w=f*suffix(nxt)
  for a in range(n):
   if mask>>a&1:N[a][b]+=w
require(forward[full]==E==984960,'Extension total')
def chain(mask,D):
 return all(D[a]>>b&1 or D[b]>>a&1 for a in range(n) for b in range(a) if mask>>a&1 and mask>>b&1)
rows=[];very=[]
for dual in (False,True):
 D,V=(U,p) if dual else (p,U)
 for a in range(n):
  for b in range(n):
   if a==b or (D[a]|V[a])>>b&1:continue
   if D[a]==D[b] and chain(V[a]&~V[b],D) and chain(V[b]&~V[a],D):very.append((dual,a,b))
   if D[a]&~D[b] or not chain(V[b]&~V[a],D):continue
   num=N[b][a] if dual else N[a][b]
   require(2*num>E,'A good pair exists')
   rows.append({'dual':dual,'a':a,'b':b,'probability':str(Fraction(num,E))})
require(not very,'A very good pair exists')
deg=[n-1-(p[i]|U[i]).bit_count() for i in range(n)]
require(max(deg)==8,'Thinness')
require(2*N[0][1]==E,'Half comparison')
# An independent ideal DP with the extra comparison 0<1.
q=p.copy();q[1]|=1
@lru_cache(None)
def added(mask):
 if mask==full:return 1
 return sum(added(mask|1<<i) for i in range(n) if not mask>>i&1 and not q[i]&~mask)
require(added(0)==N[0][1],'Added-comparison check')
from collections import Counter
hist=Counter(r['probability'] for r in rows)
require(hist=={'2377/4560':16,'717/760':16,'109/114':16},'Eligible probability census')
perm=[1,0,3,2,4,5,6,7,9,8,11,10]
require(all(bool(p[v]>>u&1)==bool(p[perm[v]]>>perm[u]&1) for u in range(n) for v in range(n)),'Automorphism')
print(json.dumps({'status':'EXACT_K4_HEIGHT_TWO_CONTROLS_PASSED','pred':p,'extensions':E,'incomparability_degrees':deg,'good_pairs':0,'very_good_pairs':0,'structurally_eligible_orientations':rows,'eligible_probability_histogram':dict(hist),'automorphism':perm,'half_pair':[0,1],'half_numerator':N[0][1],'weak_gold_strategy':'Compare 0 and 1, then repeat the determined comparison; each branch has E/2 extensions.'},indent=2))
