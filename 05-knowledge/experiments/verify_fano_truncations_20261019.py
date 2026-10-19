#!/usr/bin/env python3
from functools import lru_cache
from fractions import Fraction
import json
base=[0]*16
for j in range(7):base[8+j]=1|sum(1<<(1+i) for i in range(7) if (j-i)%7 in (0,1,3))
base[15]=254

def check(deleted):
 labels=[i for i in range(16) if i not in deleted];n=len(labels);full=(1<<n)-1
 p=[sum(1<<i for i,a in enumerate(labels) if base[b]>>a&1) for b in labels]
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
   nxt=mask|1<<b;forward[nxt]=forward.get(nxt,0)+f;w=f*suffix(nxt)
   for a in range(n):
    if mask>>a&1:N[a][b]+=w
 def chain(mask,D):return all(D[a]>>b&1 or D[b]>>a&1 for a in range(n) for b in range(a) if mask>>a&1 and mask>>b&1)
 good=[];very=[];eligible=0
 for dual in (False,True):
  D,V=(U,p) if dual else (p,U)
  for a in range(n):
   for b in range(n):
    if a==b or (D[a]|V[a])>>b&1:continue
    if a<b and D[a]==D[b] and chain(V[a]&~V[b],D) and chain(V[b]&~V[a],D):very.append({'dual':dual,'pair_original_labels':[labels[a],labels[b]]})
    if D[a]&~D[b] or not chain(V[b]&~V[a],D):continue
    eligible+=1;num=N[b][a] if dual else N[a][b]
    if 2*num<=E:good.append({'dual':dual,'pair_original_labels':[labels[a],labels[b]],'probability':str(Fraction(num,E))})
 return {'deleted_original_labels':deleted,'pred':p,'labels':labels,'extensions':E,'degrees':[n-1-(p[i]|U[i]).bit_count() for i in range(n)],'eligible_orientations':eligible,'good_pairs':good,'very_good_pairs':very}
print(json.dumps({'status':'EXACT_TWO_CONSTRUCTION_CHECK','incident':check([1,8]),'nonincident':check([1,10]),'scope':'These truncations have root degrees seven but ambient maximum degree ten, not seven-thin.'},indent=2))
