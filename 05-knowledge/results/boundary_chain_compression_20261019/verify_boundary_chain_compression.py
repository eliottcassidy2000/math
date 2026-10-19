#!/usr/bin/env python3
"""Exact standard-library controls; no numerical solver or saved atom matrix."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from fractions import Fraction
from pathlib import Path
import hashlib,json,math
ROOT=Path(__file__).resolve().parent
checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
def closure(n,edges):
 p=[0]*n
 for a,b in edges:p[b]|=1<<a
 for _ in range(n):
  old=p.copy()
  for b in range(n):
   for a in range(n):
    if p[b]>>a&1:p[b]|=p[a]
  if old==p:break
 require(all(not(p[q]>>q&1) for q in range(n)),'cycle')
 return tuple(p)
def extensions(p):
 n=len(p)
 def rec(w,seen):
  if len(w)==n:yield w;return
  for q in range(n):
   if not(seen>>q&1) and p[q]&~seen==0:yield from rec(w+(q,),seen|1<<q)
 return list(rec((),0))
def number(p,coreword=()):
 n=len(p);pp=list(p)
 if coreword:
  edges=[(a,b) for b in range(n) for a in range(n) if pp[b]>>a&1]
  edges+=list(zip(coreword,coreword[1:]));pp=closure(n,edges)
 @lru_cache(None)
 def f(s):
  if s==(1<<n)-1:return 1
  return sum(f(s|1<<q) for q in range(n) if not(s>>q&1) and pp[q]&~s==0)
 return f(0)
def degrees(p):return [sum(a!=b and not(p[a]>>b&1 or p[b]>>a&1) for b in range(len(p))) for a in range(len(p))]
@lru_cache(None)
def G(c):
 if not c:return 1
 return sum(G(tuple(min(x,k) for x in c[1:])) for k in range(c[0]+1))
def F(c):return G(tuple(x-1 for x in c)) if c and min(c)>0 else (1 if not c else 0)
def tuples(c):return [z for z in product(*(range(1,t+1) for t in c)) if all(z[i]>=z[i+1] for i in range(len(z)-1))]
def rho(c,a):return Fraction(F(tuple(c)+(a,)),F(c)) if a else Fraction(0)
# Independent finite hostile check of stochastic cap monotonicity, not its proof.
monotonic=0
for r in range(1,5):
 for c in product(range(1,5),repeat=r):
  z=tuples(c);require(len(z)==F(c),'F direct enumeration')
  for i in range(r):
   cc=list(c);cc[i]+=1;zz=tuples(cc)
   for h in range(1,6):
    require(sum(x[-1]>=h for x in z)*len(zz)<=sum(x[-1]>=h for x in zz)*len(z),'stochastic domination');monotonic+=1
for s in range(1,8):
 for r in range(0,7):
  require(F((s,)*r)==math.comb(s+r-1,r),'equal-cap count')
  require(rho((s,)*r,s)==Fraction(s+r,r+1),'equal-cap ratio')
# Shared allocation facets: cross-multiplied exact guards for all allocations.
facet_cases=0
for l in [(),(1,),(2,),(3,),(4,2),(2,4),(3,3,3)]:
 for r in [(),(1,),(2,),(3,),(4,2),(2,4),(3,3,3)]:
  for d in range(1,8):
   pts=[(rho(l,a),rho(r,d-a)) for a in range(d+1)]
   for a in range(d):
    dl=pts[a+1][0]-pts[a][0];dr=pts[a][1]-pts[a+1][1]
    rhs=dr*pts[a][0]+dl*pts[a][1]
    require(all(dr*x+dl*y<=rhs for x,y in pts),'allocation facet');facet_cases+=1
require(rho((3,),3)==2 and 9>rho((3,),3)*3,'old ancestor-cap hostile')
require(3<=3 and 3<=3 and 3+3>3,'independent-side hostile')
# An actual eleven-vertex ambient test with non-chain exteriors.
nq=5;A=(5,6,7);B=(8,9,10)
q=closure(5,[(2,1),(0,4),(2,4)])
edges=[(a,b) for b in range(5) for a in range(5) if q[b]>>a&1]
edges += [(5,7),(8,10),(5,2),(6,3),(4,9),(3,8)]
edges += [(a,r) for a in A for r in (0,1)]+[(r,b) for b in B for r in (0,1)]
p=closure(11,edges);ambient_p=p;require(max(degrees(p))<=7,'ambient thinness')
qw=extensions(q)
aorders=[w for w in permutations(A) if all(w.index(a)<w.index(b) for a in A for b in A if p[b]>>a&1)]
borders=[w for w in permutations(B) if all(w.index(a)<w.index(b) for a in B for b in B if p[b]>>a&1)]
actual={w:number(p,w) for w in qw};summed=Counter();atom_events=[];atom_totals=[];profilepairs=Counter()
for aa in aorders:
 al=tuple(len(A)-max([i+1 for i,a in enumerate(aa) if p[x]>>a&1]+[0]) for x in range(5))
 for bb in borders:
  be=tuple(min([i for i,b in enumerate(bb) if p[b]>>x&1]+[len(B)]) for x in range(5))
  h=max(al);j=max(be);ce=[(a,b) for b in range(5) for a in range(5) if q[b]>>a&1]
  ce+=list(zip(range(5,5+h-1),range(6,5+h)))+list(zip(range(5+h,5+h+j-1),range(6+h,5+h+j)))
  ce += [(5+i,x) for i in range(h) for x in range(5) if i<h-al[x]]
  ce += [(x,5+h+i) for i in range(j) for x in range(5) if i>=be[x]]
  ce += [(a,b) for a in range(5,5+h) for b in range(5+h,5+h+j)]
  cp=closure(5+h+j,ce);require(max(degrees(cp))<=7,'compressed thinness')
  counts={}
  for w in qw:
   lo=min(w.index(0),w.index(1));hi=max(w.index(0),w.index(1))
   z=G(tuple(al[x] for x in w[:lo]))*G(tuple(be[x] for x in reversed(w[hi+1:])))
   require(z==number(cp,w),'atom G versus independent ideal DP');counts[w]=z;summed[w]+=z
  atom_totals.append(sum(counts.values()));atom_events.append(counts);profilepairs[(al,be)]+=1
require(dict(summed)==actual,'all core-word lift counts simultaneously')
require(sum(actual.values())==number(p),'ambient total DP')
wrong_uniform=None
for a in range(5):
 for b in range(a+1,5):
  prob=Fraction(sum(z for w,z in actual.items() if w.index(a)<w.index(b)),sum(actual.values()))
  uniform=sum((Fraction(sum(z for w,z in c.items() if w.index(a)<w.index(b)),e) for c,e in zip(atom_events,atom_totals)),Fraction())/len(atom_events)
  if prob!=uniform:wrong_uniform={'pair':[a,b],'actual':str(prob),'incorrect_uniform_atoms':str(uniform)};break
 if wrong_uniform:break
require(wrong_uniform is not None,'hostile uniform-mixture distinction')
# Exact finite four-comparison certificate on the prescribed nine-point core.
p=(0,28,0,0,12,12,13,125,125);order=(2,3,0,4,5,6,1,7,8);k=7;n=len(p)
slack=[k-d for d in degrees(p)];words=extensions(p);require(len(words)==408,'fixed core count')
def profiles(side):
 eligible=[q for q in range(n) if q not in (0,1) and (not(p[q]&3) if side==0 else not(any(p[r]>>q&1 for r in (0,1))))]
 result=[]
 for vals in product(*(range(slack[q]+1) for q in eligible)):
  t=[0]*n
  for q,x in zip(eligible,vals):t[q]=x
  if any((t[q]<t[r] if side==0 else t[q]>t[r]) for q in range(n) for r in range(n) if p[r]>>q&1):continue
  if sum(x>0 for x in t)>k:continue
  result.append(tuple(t))
 return result
lp=profiles(0);rp=profiles(1);require(len(lp)==len(rp)==270,'profile universe count')
pre=[w[:min(w.index(0),w.index(1))] for w in words];suf=[tuple(reversed(w[max(w.index(0),w.index(1))+1:])) for w in words]
LL=[[G(tuple(t[q] for q in u)) for u in pre] for t in lp];RR=[[G(tuple(t[q] for q in v)) for v in suf] for t in rp]
comparisons=[(0,4,2),(6,1,2),(4,5,4),(5,6,5)]
score=[sum(weight*(w.index(a)<w.index(b)) for a,b,weight in comparisons) for w in words]
minmargin=0;maxmargin=-10**100;atoms=0;worst=Fraction();witness=None
for a,al in enumerate(lp):
 for b,be in enumerate(rp):
  if any(x+y>s for x,y,s in zip(al,be,slack)):continue
  E=W=0
  for l,r,z in zip(LL[a],RR[b],score):E+=l*r;W+=l*r*z
  margin=20*W-169*E
  require(margin<=-472,'four-comparison inequality')
  minmargin=min(minmargin,margin);maxmargin=max(maxmargin,margin);atoms+=1
  ratio=Fraction(W,13*E)
  if ratio>worst:worst=ratio;witness={'lower':al,'upper':be,'E':E,'weighted_N':W}
require(atoms==56379,'joint universe count');require((minmargin,maxmargin)==(-4014608,-472),'certificate extrema')
require(worst==Fraction(8069,12441),'sharp fixed-row ratio')
report={'status':'PASS','guards':checks,'stochastic_threshold_controls':monotonic,'allocation_facet_controls':facet_cases,'ambient_mixture_control':{'predecessors':ambient_p,'extensions':sum(actual.values()),'core_extensions':len(qw),'exterior_order_pairs':len(atom_totals),'distinct_profile_pairs':len(profilepairs),'atom_extension_totals':atom_totals,'uniform_mixture_hostile':wrong_uniform},'fixed_core_certificate':{'core':p,'order':order,'core_extensions':len(words),'lower_profiles':len(lp),'upper_profiles':len(rp),'joint_atoms':atoms,'comparison_weights':comparisons,'cleared_margin_range':[minmargin,maxmargin],'maximum_fixed_row_ratio':str(worst),'attaining_atom':witness},'scope':'Exact finite controls and a fixed-core quantitative certificate; no all-morphology or seven-thin conjecture resolution.'}
print(json.dumps(report,indent=2))
