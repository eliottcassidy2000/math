#!/usr/bin/env python3
"""Exact finite controls for tied suffix weights and passive-attachment bias."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from fractions import Fraction
import json,math
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
  if p==old:break
 require(all(not(p[a]>>a&1) for a in range(n)),'cycle')
 return tuple(p)
def words(p,labels):
 labels=tuple(labels);full=sum(1<<a for a in labels)
 def rec(w,s):
  if s==full:yield w;return
  for a in labels:
   if not(s>>a&1) and not(p[a]&full&~s):yield from rec(w+(a,),s|1<<a)
 return list(rec((),0))
def inc(p,a,b):return a!=b and not(p[a]>>b&1 or p[b]>>a&1)
def e(p,labels):
 full=sum(1<<a for a in labels)
 @lru_cache(None)
 def f(s):
  if s==full:return 1
  return sum(f(s|1<<a) for a in labels if not(s>>a&1) and not(p[a]&full&~s))
 return f(0)
def examine(name,p,Q,A,k,expected=None):
 Q=tuple(Q);A=tuple(A);T=tuple(a for a in A if any(inc(p,a,q) for q in Q));A0=tuple(a for a in A if a not in T)
 require(all(sum(inc(p,a,b) for b in range(len(p)))<=k for a in range(len(p))),'thinness')
 require(all(not(p[z]>>a&1) for a in T for z in A0),'no active below passive')
 delta=[k-sum(inc(p,q,x) for x in Q) for q in Q];r=max(delta,default=0)
 aw=words(p,A);require(len(aw)==e(p,A),'independent exterior DP')
 @lru_cache(None)
 def H(S):return e(p,tuple(a for a in A if a not in S))
 def M(S):return tuple(a for a in T if a not in S and not(any(p[b]>>a&1 for b in T if b not in S)))
 def profile(suffix):
  out=[]
  for q in Q:
   h=0
   for a in reversed(suffix):
    if not inc(p,a,q):break
    h+=1
   out.append(h)
  return tuple(out)
 actual_states=Counter();actual_profiles=Counter()
 for w in aw:
  h=0
  for a in reversed(w):
   if a not in T or h==r:break
   h+=1
  suffix=w[len(w)-h:] if h else ()
  actual_states[suffix]+=1;actual_profiles[profile(w)]+=1
  require(profile(suffix)==profile(w),'truncation preserves actual profile')
 states=Counter();profiles=Counter();filters=set();rows=[]
 def visit(deleted):
  S=tuple(sorted(deleted));filters.add(S);hh=H(S);choices=M(S)
  # A fixed ordered suffix has exactly H(S) lifts, checked independently.
  suffix=tuple(reversed(deleted))
  require(sum(w[-len(suffix):]==suffix if suffix else True for w in aw)==hh,'same-set suffix multiplicity')
  require(hh-sum(H(tuple(sorted(S+(a,)))) for a in choices)>=0,'nonnegative stopping flow')
  for a in choices:
   Sa=tuple(sorted(S+(a,)));tq=sum(inc(p,a,q) for q in Q);removed=sum(inc(p,a,b) for b in S);bound=k-tq-removed+1
   require(bound>=1 and H(Sa)<=hh<=bound*H(Sa),'branch lower capacity')
  for a,b in combinations(choices,2):
   sa=tuple(sorted(S+(a,)));sb=tuple(sorted(S+(b,)));sab=tuple(sorted(S+(a,b)))
   left=Fraction(H(sa),hh)*Fraction(H(sab),H(sa));right=Fraction(H(sb),hh)*Fraction(H(sab),H(sb))
   require(left==right,'diamond equality')
  if len(deleted)==r:
   weight=hh
  else:
   weight=hh-sum(H(tuple(sorted(S+(a,)))) for a in choices)
   for a in choices:visit(deleted+(a,))
  if weight:states[suffix]+=weight;profiles[profile(suffix)]+=weight
  rows.append({'filter':S,'H':hh,'J':hh-sum(H(tuple(sorted(S+(a,)))) for a in choices)})
 visit(())
 require(states==actual_states,'terminal state counts match actual extension words')
 require(profiles==actual_profiles,'profile multiplicity reconstruction')
 require(sum(profiles.values())==len(aw),'profile total')
 for j,q in enumerate(Q):
  I=tuple(sorted(a for a in A if inc(p,a,q)));ss=len(I);mass=sum(v for al,v in profiles.items() if al[j]==ss)
  require(max(al[j] for al in profiles)==ss,'support maximum recovers actual allocation')
  require(mass==H(I)*e(p,I),'maximal-gap exact mass')
  require(mass*k**ss>=e(p,I)*len(aw),'maximal-gap positive capacity')
 require(len(filters)<=math.comb(k+r,k),'filter count bound')
 if expected:
  for S,v in expected.items():require(H(S)==v,'prescribed sharp H value')
 unique={tuple(z['filter']):(z['H'],z['J']) for z in rows}
 return {'name':name,'ambient_predecessors':p,'k':k,'core':Q,'active':T,'passive':A0,'depth':r,'exterior_extensions':len(aw),'filters':len(filters),'terminal_states':[{'suffix':s,'count':v} for s,v in sorted(states.items())],'profiles':[{'profile':s,'count':v} for s,v in sorted(profiles.items())],'H_rows':[{'filter':s,'H':v[0],'J':v[1]} for s,v in sorted(unique.items())]}
reports=[]
# Sharp k=5 example, and k=4 variant exercising truncation before all active vertices.
for k,m in [(5,3),(4,2)]:
 # Q=0,1,2,3; active a=4,b=5; passive chain6,...
 edges=[(4,q) for q in (0,1,3)]+[(5,q) for q in (0,1,2)]
 edges+=list(zip(range(6,6+m-1),range(7,6+m)))
 edges += [(c,q) for c in range(6,6+m) for q in range(4)]+[(c,4) for c in range(6,6+m)]
 p=closure(6+m,edges)
 reports.append(examine(f'sharp_passive_bias_k{k}',p,range(4),range(4,6+m),k,{():m+2,(4,):m+1,(5,):1,(4,5):1}))
# Non-chain active/passive interface from the independent eleven-vertex mixture control.
Q=range(5);A=(5,6,7);B=(8,9,10)
edges=[(2,1),(0,4),(2,4),(5,7),(8,10),(5,2),(6,3),(4,9),(3,8)]
edges += [(a,r) for a in A for r in (0,1)]+[(r,b) for b in B for r in (0,1)]
p=closure(11,edges);reports.append(examine('mixed_interface_lower',p,Q,A,7))
# Exact upper-side dual control.
up=tuple(sum(1<<b for b in range(11) if p[b]>>a&1) for a in range(11))
reports.append(examine('mixed_interface_upper_dual',up,Q,B,7))
# Empty active interface and depth-zero boundary cases.
p=closure(5,[(2,3),(3,4)]+[(a,q) for a in (2,3,4) for q in (0,1)])
reports.append(examine('zero_active',p,(0,1),(2,3,4),1))
# Profile-factorization minors before normalization, not normalized atom-weight minors.
ma={(0,1):4,(1,0):1};mb={(0,):2,(1,):3};E={(a,b):i+2 for i,(a,b) in enumerate((a,b) for a in ma for b in mb)}
total=sum(ma[a]*mb[b]*E[a,b] for a in ma for b in mb);pi={(a,b):Fraction(ma[a]*mb[b]*E[a,b],total) for a in ma for b in mb}
a,aa=ma;b,bb=mb
require(pi[a,b]/E[a,b]*pi[aa,bb]/E[aa,bb]==pi[a,bb]/E[a,bb]*pi[aa,b]/E[aa,b],'correct reweighted rank one')
require(pi[a,b]*pi[aa,bb]!=pi[a,bb]*pi[aa,b],'normalized mixture generally not rank one')
print(json.dumps({'status':'PASS','guards':checks,'cases':reports,'rank_one_control':{'multiplicities_A':[4,1],'multiplicities_B':[2,3],'atom_totals':list(E.values()),'total':total,'normalized_weights':[str(x) for x in pi.values()]},'scope':'Exact selected suffix-state, passive-bias, duality and rank-one controls. No characterization of all realizable H arrays or resolution of balance/Gold.'},indent=2))
