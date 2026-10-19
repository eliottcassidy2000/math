#!/usr/bin/env python3
"""Exact one-poset control; no poset census. Explicit guards survive python -O."""
from functools import lru_cache
from fractions import Fraction
import json

def require(test,message):
 if not test:raise RuntimeError(message)
n=16;full=(1<<n)-1
# x=0, L_i=1+i, R_j=8+j, y=15; cyclic Fano difference set.
p=[0]*n
for j in range(7):p[8+j]=1|sum(1<<(1+i) for i in range(7) if (j-i)%7 in (0,1,3))
p[15]=sum(1<<i for i in range(1,8))
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
require(forward[full]==E==3557240064,'Extension total')
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
require(all(Fraction(N[0][i],E)==Fraction(216,389) for i in range(1,8)),'Root/minimum probability')
require(all(Fraction(N[j][15],E)==Fraction(216,389) for j in range(8,15)),'Dual contact probability')
# Independently count with the single extra comparison x<L0.
q=p.copy();q[1]|=1
@lru_cache(None)
def added(mask):
 if mask==full:return 1
 return sum(added(mask|1<<i) for i in range(n) if not mask>>i&1 and not q[i]&~mask)
require(added(0)==N[0][1]==1975228416,'Added-relation cross-check')
require(not p[0]&~p[15] and not U[15]&~U[0],'Root criticality')
common=[v for v in range(1,15) if not ((p[0]|U[0])>>v&1) and not ((p[15]|U[15])>>v&1)]
require(not common,'Common incomparable neighbor')
deg=[n-1-(p[i]|U[i]).bit_count() for i in range(n)]
perm=[0]+[1+(i+1)%7 for i in range(7)]+[8+(j+1)%7 for j in range(7)]+[15]
require(all(bool(p[b]>>a&1)==bool(p[perm[b]]>>perm[a]&1) for a in range(n) for b in range(n)),'Rotation automorphism')
print(json.dumps({'status':'EXACT_SINGLE_POSET_CONTROLS_PASSED','pred':p,'extensions':E,'root_pair':[0,15],'root_probability':str(Fraction(N[0][15],E)),'common_incomparable_neighbors':common,'incomparability_degrees':deg,'good_pairs':0,'very_good_pairs':0,'structurally_eligible_orientations':rows,'contact_probability':'216/389','contact_numerator':N[0][1],'nontrivial_automorphism':perm,'scope':'Refutes only unrestricted c=0 implies existence of a good pair. Eleven-thin, with Gold already guaranteed by nontrivial automorphism; no conjecture counterexample.'},indent=2))
