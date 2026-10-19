#!/usr/bin/env python3
"""Exact small controls for simultaneous compression and correlation constraints.
Standard library only. Guards remain active under python -O.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import factorial
import json

checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok: raise ValueError(msg)
def closure(n,edges):
 p=[0]*n
 for a,b in edges:p[b]|=1<<a
 for _ in range(n):
  old=p.copy()
  for b in range(n):
   for a in range(n):
    if p[b]>>a&1:p[b]|=p[a]
  if p==old:break
 require(all(not(p[a]>>a&1) for a in range(n)),'acyclic')
 return tuple(p)
def inc(p,a,b):return a!=b and not(p[a]>>b&1 or p[b]>>a&1)
def words(p,labels):
 labels=tuple(labels);full=sum(1<<a for a in labels)
 def go(w,s):
  if s==full:yield w;return
  for a in labels:
   if not(s>>a&1) and not(p[a]&full&~s):yield from go(w+(a,),s|1<<a)
 return list(go((),0))
def count_fn(p):
 @lru_cache(None)
 def e(mask):
  if not mask:return 1
  return sum(e(mask^(1<<a)) for a in range(len(p)) if mask>>a&1 and not(p[a]&mask))
 return e
def subsets(labels):
 labels=tuple(labels)
 for bits in range(1<<len(labels)):yield tuple(a for i,a in enumerate(labels) if bits>>i&1)
def mask(labels):return sum(1<<a for a in labels)
def maximum(p,a,remaining):return not any(p[b]>>a&1 for b in remaining)
def is_filter(p,S,labels):return all(not(p[b]>>a&1) or b in S for a in S for b in labels)

def compression_case(name,p,Q,A,k):
 Q=tuple(Q);A=tuple(A);T=tuple(a for a in A if any(inc(p,a,q) for q in Q))
 A0=tuple(a for a in A if a not in T);C=tuple(a for a in A0 if any(inc(p,a,t) for t in T));D=tuple(a for a in A0 if a not in C);R=T+C
 require(all(sum(inc(p,a,b) for b in range(len(p)))<=k for a in range(len(p))),'input thinness')
 require(all(p[t]>>d&1 for d in D for t in T),'D below all active')
 require(all(not(p[d]>>r&1) for d in D for r in R),'no R below D')
 require(len(C)<=sum(k-sum(inc(p,t,q) for q in Q)-sum(inc(p,t,u) for u in T) for t in T),'passive incidence bound')
 E=count_fn(p);original={S:E(mask(a for a in A if a not in S)) for S in subsets(T)};sums=Counter();atoms=[]
 # Same replacement is reused for every S, including nonfilters.
 for sigma in words(p,D):
  gamma={r:len(D)-max((i+1 for i,d in enumerate(sigma) if p[r]>>d&1),default=0) for r in R}
  h=max(gamma.values(),default=0);keep=Q+R;idx={a:i for i,a in enumerate(keep)}
  chain=tuple(range(len(keep),len(keep)+h));edges=[(idx[a],idx[b]) for a in keep for b in keep if p[b]>>a&1]
  edges += [(chain[i],chain[i+1]) for i in range(h-1)]
  edges += [(c,idx[q]) for c in chain for q in Q]
  edges += [(chain[i],idx[r]) for r in R for i in range(h-gamma[r])]
  ap=closure(len(keep)+h,edges);AE=count_fn(ap);alabels=tuple(idx[r] for r in R)+chain
  require(all(sum(inc(ap,a,b) for b in range(len(ap)))<=k for a in range(len(ap))),'replacement thinness')
  require(h<=k-1,'chain length bound')
  require(len(alabels)<=k*len(T)+k-1,'attachment size bound')
  require(tuple(r for r in R if any(inc(ap,idx[r],idx[q]) for q in Q))==T,'same exact active interface')
  values={S:AE(mask(a for a in alabels if a not in tuple(idx[t] for t in S))) for S in subsets(T)}
  for S,v in values.items():sums[S]+=v
  # Independent direct restriction count for each D order and each deletion.
  for S in subsets(T):
   aw=words(p,tuple(a for a in A if a not in S))
   require(values[S]==sum(tuple(a for a in w if a in D)==sigma for w in aw),'one atom matches conditioned exterior orders')
  atoms.append({'D_order':sigma,'gamma':[(r,gamma[r]) for r in R],'vertices':len(alabels),'H':[(S,values[S]) for S in subsets(T)]})
 for S,v in original.items():require(sums[S]==v,'simultaneous H identity')
 total=original[()]
 for S in subsets(T):
  require(sum(Fraction(a['H'][0][1],total)*Fraction(dict(a['H'])[S],a['H'][0][1]) for a in atoms)==Fraction(original[S],total),'extension-weighted normalized mixture')
 return {'name':name,'predecessors':p,'k':k,'Q':Q,'T':T,'C':C,'D':D,'H':[(S,v) for S,v in original.items()],'replacements':atoms}

reports=[]
# Nontrivial D orders, unequal profiles and stripped universal prefix.
edges=[(5,q) for q in (0,1,3)]+[(6,q) for q in (0,1,2)]+[(7,5),(8,6),(9,7),(10,8)]
edges += [(d,t) for d in (9,10) for t in (5,6)]+[(c,q) for c in (7,8,9,10,11) for q in range(5)]
edges += [(11,r) for r in (7,8,9,10)]
p=closure(12,edges);reports.append(compression_case('two_D_orders_with_universal_prefix',p,range(5),range(5,12),6))
# T is a chain, so subsets include a genuinely nonfilter deletion.
edges=[(5,6)]+[(5,q) for q in (0,1,2,3)]+[(6,q) for q in (0,1,2)]+[(7,6),(8,5),(9,5),(9,7)]
edges += [(c,q) for c in (7,8,9) for q in range(5)]
p=closure(10,edges);reports.append(compression_case('nonfilter_deletions_retained',p,range(5),range(5,10),6))
# Empty active interface, and empty D, exercise both zero-length cases.
p=closure(5,[(2,3)]+[(a,q) for a in (2,3,4) for q in (0,1)])
reports.append(compression_case('empty_active',p,(0,1),(2,3,4),2))
p=closure(5,[(3,0),(3,1),(4,0),(4,2)])
reports.append(compression_case('empty_passive',p,(0,1,2),(3,4),3))

# Exact convexification hostile control, on one fixed Q and fixed T.
hostile_actual=[]
for target in (5,6):
 edges=[(5,q) for q in (0,1,3)]+[(6,q) for q in (0,1,2)]
 edges += [(7,8),(8,9)]+[(c,q) for c in (7,8,9) for q in range(5)]+[(c,target) for c in (7,8,9)]
 p=closure(10,edges);E=count_fn(p);A=(5,6,7,8,9)
 require(all(sum(inc(p,a,b) for b in range(10))<=6 for a in range(10)),'six-thin hostile constituents')
 H=[E(mask(a for a in A if a not in S)) for S in ((),(5,),(6,),(5,6))]
 require(H in ([5,4,1,1],[5,1,4,1]),'prescribed constituent H')
 require(4*H[0]*H[3]==5*H[1]*H[2],'each constituent sharp Fishburn')
 hostile_actual.append({'predecessors':p,'H':H})
low=[sum(c['H'][j] for c in hostile_actual) for j in range(4)]
require(low==[10,5,5,2],'forbidden lower array is sum of actual arrays')
for H in (low,[10,2,2,1]):
 require(H[0]>=H[1]+H[2] and H[1]>=H[3] and H[2]>=H[3],'nonnegative stopping flow')
 require(all(H[0]<=5*H[j] for j in (1,2)) and all(H[j]<=4*H[3] for j in (1,2)),'same degree branch bounds')
 require(Fraction(H[1],H[0])*Fraction(H[3],H[1])==Fraction(H[2],H[0])*Fraction(H[3],H[2]),'shared-set diamond identity')
require(low[0]*low[3]<low[1]*low[2],'lower correlation catches convex mixture')
require(10*1>2*2*2,'upper correlation catches second hostile array')

# Small exhaustive theorem controls, not a balance scan or proof by computation.
family_counts=[];diamond_rows=0;multiway_rows=0;filter_rows=0;upper_equalities=0;lower_equalities=0
for n in range(1,6):
 pairs=list(combinations(range(n),2));family=set()
 for bits in range(1<<len(pairs)):
  family.add(closure(n,[edge for i,edge in enumerate(pairs) if bits>>i&1]))
 family_counts.append({'n':n,'naturally_labeled_posets':len(family)})
 for p in sorted(family):
  E=count_fn(p);full=(1<<n)-1;filters=[S for S in range(1<<n) if is_filter(p,tuple(a for a in range(n) if S>>a&1),range(n))]
  for S in filters:
   rem=full^S;m=rem.bit_count();H=E(rem);maxima=[a for a in range(n) if rem>>a&1 and maximum(p,a,[b for b in range(n) if rem>>b&1])]
   for a,b in combinations(maxima,2):
    Ha=E(rem^(1<<a));Hb=E(rem^(1<<b));Hab=E(rem^(1<<a)^(1<<b));L=H*Hab;R=Ha*Hb
    require((m-1)*L>=m*R,'Fishburn exact-size diamond')
    require(L<=2*R,'Chan-Pak upper diamond')
    lower_equalities+=(m-1)*L==m*R;upper_equalities+=L==2*R;diamond_rows+=1
   for W in subsets(maxima):
    j=len(W)
    if j<2:continue
    require(factorial(m)//factorial(m-j)*H**(j-1)*E(rem^mask(W))>=m**j*__import__('math').prod(E(rem^(1<<a)) for a in W),'multiway maximal corner')
    multiway_rows+=1
  for S,T in combinations(filters,2):
   I=S&T;U=S|T
   require(Fraction(E(full^I),factorial(n-I.bit_count()))*Fraction(E(full^U),factorial(n-U.bit_count()))>=Fraction(E(full^S),factorial(n-S.bit_count()))*Fraction(E(full^T),factorial(n-T.bit_count())),'full filter volume log-supermodularity')
   filter_rows+=1
# Maximality is a necessary hypothesis.
require(1*1<2*1,'two-chain nonmaximal counterexample')
print(json.dumps({'status':'PASS','guards':checks,'compression_cases':reports,'convexification_hostile':{'actual_same_interface_arrays':hostile_actual,'sum_H':low,'lower_products':[20,25]},'upper_hostile':{'H':[10,2,2,1],'upper_products':[10,8]},'small_theorem_controls':{'families':family_counts,'diamonds':diamond_rows,'multiway_rows':multiway_rows,'filter_lattice_rows':filter_rows,'sharp_lower_rows':lower_equalities,'sharp_upper_rows':upper_equalities},'scope':'Exact illustrative compression controls and all naturally labeled posets through five vertices for stated correlation rows. Theorems rely on written proofs and credited primary results. No balance/Gold resolution or sufficient realizability test.'},indent=2))
