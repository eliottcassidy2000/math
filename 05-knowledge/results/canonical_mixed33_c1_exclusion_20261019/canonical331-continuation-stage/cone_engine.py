#!/usr/bin/env python3
"""Exact finite morphology certificate search; floats suggest, integers decide."""
from itertools import combinations,product
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import math,json,gzip,time,hashlib
ROOT=Path(__file__).resolve().parent

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def closure(p):
 p=list(p)
 for k in range(len(p)):
  for j in range(len(p)):
   if p[j]>>k&1:p[j]|=p[k]
 return tuple(p)
def degrees(p):
 return tuple(len(p)-1-p[v].bit_count()-sum(bool(t>>v&1) for t in p) for v in range(len(p)))
def extensions(p):
 out=[];n=len(p)
 def visit(w,m):
  if len(w)==n:out.append(tuple(w));return
  for v in range(n):
   if not(m>>v&1) and not(p[v]&~m):visit(w+[v],m|1<<v)
 visit([],0);return out
def generators(p,words,inc):
 gens={};rankwords=[]
 for o in words:
  r=[o.index(v) for v in range(len(p))];rankwords.append((o,r));h=min(r[:2]);k=max(r[:2])+1
  counts=[1]+[int(r[a]<r[b]) for a,b in inc]
  for j in range(h+1):
   for t in range(k,len(p)+1):
    key=(o[:j],o[t:]);col=gens.setdefault(key,[0]*len(counts))
    for i,v in enumerate(counts):col[i]+=v
 keys=sorted(gens,key=lambda x:(len(x[0]),x[0],len(x[1]),x[1]))
 return keys,[gens[k] for k in keys],rankwords

def defect_bound(p,u,v):
 n=len(p);up=[sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n)];d=degrees(p)
 seeds=3|sum(1<<i for i in range(n) if d[i]==7)
 qa=0 if ((1<<u)|p[u])&seeds else 7-d[u]
 qb=0 if ((1<<v)|up[v])&seeds else 7-d[v]
 return (p[v]&~p[u]).bit_count()+(up[u]&~up[v]).bit_count()+qa+qb

def forced_cycle(p,bounds):
 n=len(p);reach=[v for v in p]
 arcs=[(a,b) for b in range(n) for a in range(n) if p[b]>>a&1]
 arcs+=[(b,a) for (a,b),h in bounds.items() if h<=1]
 common=[z for z in range(2,n) if not any(p[z]>>r&1 or p[r]>>z&1 for r in (0,1))]
 arcs+=[(0,1)]+[(0,z) for z in common]+[(z,1) for z in common]
 for a,b in arcs:reach[b]|=1<<a
 reach=closure(reach)
 if any(reach[v]>>v&1 for v in range(n)):return arcs
 return None

def exact_dual(rows,candidate):
 support=[i for i,w in enumerate(candidate) if w>1e-9]
 for cap in (1000,1000000,1000000000):
  w={i:Fraction(float(candidate[i])).limit_denominator(cap) for i in support}
  if not w:continue
  den=math.lcm(*(v.denominator for v in w.values()));w={i:int(v*den) for i,v in w.items() if v>0};g=math.gcd(*w.values());w={i:v//g for i,v in w.items()}
  vals=[sum(v*int(rows[i,j]) for i,v in w.items()) for j in range(rows.shape[1])]
  if max(vals)<=0:return w
 return None

def split_certificate(p,keys,signed,ids):
 """A complete binary slack-one split tree; failure preserves an obstacle."""
 slack=[v for v,d in enumerate(degrees(p)) if d==6];failures=[]
 def solve(keep,remaining,fa,fb):
  rows=signed[np.ix_(ids,keep)];w=None
  good=np.flatnonzero(rows.max(axis=1)<=0)
  if len(good):w={int(good[0]):1}
  if w is None:
   # Positive small integer suggestions, checked exactly below.
   for i,j,k in combinations(range(len(ids)),3):
    for a,b,c in ((1,1,1),(2,1,2),(2,2,1),(1,2,2)):
     if max(a*rows[i]+b*rows[j]+c*rows[k])<=0:w={i:a,j:b,k:c};break
    if w is not None:break
  if w is None:
   sol=linprog(np.r_[np.zeros(len(keep)),-1.],A_ub=np.c_[-rows,np.ones(len(ids))],b_ub=np.zeros(len(ids)),A_eq=np.r_[np.ones(len(keep)),0.][None,:],b_eq=[1.],bounds=[(0,None)]*len(keep)+[(None,None)],method='highs')
   require(sol.success,sol.message)
   if sol.x[-1]<=1e-9:w=exact_dual(rows,-sol.ineqlin.marginals)
  if w is not None:
   full={ids[i]:weight for i,weight in w.items()}
   vals=[sum(weight*int(signed[i,j]) for i,weight in full.items()) for j in keep]
   require(max(vals)<=0,'split exact inequality')
   return {'leaf':True,'forced_A':sorted(fa),'forced_B':sorted(fb),'generator_indices':keep,'weights':full,'max_margin':max(vals)}
  if not remaining:
   support=np.flatnonzero(sol.x[:-1]>1e-9);witness=None
   for cap in (1000,1000000,1000000000):
    mix={int(j):Fraction(float(sol.x[j])).limit_denominator(cap) for j in support};mass=sum(mix.values());mix={j:v/mass for j,v in mix.items()}
    vals=[sum(v*int(row[j]) for j,v in mix.items()) for row in rows]
    if min(vals)>0:
     witness={'forced_A':sorted(fa),'forced_B':sorted(fb),'generator_indices':keep,'generator_coefficients':[[keep[j],str(v)] for j,v in mix.items()],'signed_raw_margins':list(map(str,vals))};break
   require(witness is not None,'uncertified terminal split branch')
   failures.append(witness);return None
  v=remaining[0];rest=remaining[1:]
  na=fa|{t for t in range(len(p)) if t==v or p[t]>>v&1}
  nb=fb|{t for t in range(len(p)) if t==v or p[v]>>t&1}
  left=solve([j for j in keep if not na.intersection(keys[j][0])],rest,na,fb)
  if left is None:return None
  right=solve([j for j in keep if not nb.intersection(keys[j][1])],rest,fa,nb)
  if right is None:return None
  return {'split_vertex':v,'A_below':left,'B_above':right}
 # A failed unsplit cone has already been checked; split at the first slack vertex.
 if not slack:return None,None
 tree=solve(list(range(len(keys))),slack,set(),set())
 return tree,(failures[0] if failures else None)

def tree_leaves(t):
 if t.get('leaf'):return [t]
 return tree_leaves(t['A_below'])+tree_leaves(t['B_above'])

def scan(p):
 n=len(p);inc=[(a,b) for a,b in combinations(range(n),2) if not(p[b]>>a&1 or p[a]>>b&1)]
 bounds={(a,b):defect_bound(p,a,b) for a,b in inc for a,b in ((a,b),(b,a))}
 for a,b in inc:
  if bounds[a,b]<=1 and bounds[b,a]<=1:return {'analytic_balanced':[a,b]}
 cyc=forced_cycle(p,bounds)
 if cyc:return {'analytic_cycle':cyc}
 words=extensions(p);valid=[];anal=Counter();short=0;canonical=0
 im=[sum(1<<v for v in range(n) if v!=u and not(p[u]>>v&1 or p[v]>>u&1)) for u in range(n)];root_common=(im[0]&im[1]).bit_count();cc={(a,b):(im[a]&im[b]).bit_count() for a,b in inc}
 for o in words:
  r=[o.index(v) for v in range(n)];span=r[1]-r[0]
  if span<=0 or any(abs(r[a]-r[b])>span for a,b in inc):continue
  if span<=3:short+=1;continue
  wit=next(((a,b) if r[a]<r[b] else (b,a) for a,b in inc if bounds[(a,b) if r[a]<r[b] else (b,a)]<=1),None)
  if wit:anal[str(wit)]+=1
  else:
   if any(abs(r[a]-r[b])==span and cc[a,b]>root_common for a,b in inc):canonical+=1
   else:valid.append((o,r))
 rec={'canonical_maxspan_maxcommon':True,'canonical_excluded_orders':canonical,'canonical_fully_removed':bool(canonical and not valid),'six_thin_theorem':True,'short_span_orders':short,'analytic_directed_orders':sum(anal.values()),'residual_orders':len(valid),'certificates':[],'LP_calls':0}
 if not valid:return rec
 keys,cols,_=generators(p,words,inc);sat={v for v,d in enumerate(degrees(p)) if d==7}
 keep=[j for j,(u,v) in enumerate(keys) if not(sat.intersection(u) or sat.intersection(v))];keys=[keys[j] for j in keep];cols=[cols[j] for j in keep]
 arr=np.array(cols,dtype=np.int64);E=arr[:,0];signed=np.empty((2*len(inc),len(cols)),dtype=np.int64);signed[0::2]=(3*arr[:,1:]-2*E[:,None]).T;signed[1::2]=(E[:,None]-3*arr[:,1:]).T
 cache=[];hist=Counter();stream=hashlib.sha256()
 for o,r in valid:
  ids=[2*i+int(r[a]>r[b]) for i,(a,b) in enumerate(inc)];idset=set(ids);hit=None
  for k,(support,w) in enumerate(cache):
   if support<=idset:hit=k;break
  if hit is None:
   rows=signed[ids];w=None
   good=np.flatnonzero(rows.max(axis=1)<=0)
   if len(good):w={int(good[0]):1};kind='one'
   if w is None:
    for i in range(len(rows)):
     good=np.flatnonzero((rows[i]+rows[i:]).max(axis=1)<=0)
     if len(good):w=dict(Counter((i,i+int(good[0]))));kind='two';break
   if w is None:
    rec['LP_calls']+=1
    sol=linprog(np.r_[np.zeros(len(cols)),-1.],A_ub=np.c_[-rows/E[None,:],np.ones(len(inc))],b_ub=np.zeros(len(inc)),A_eq=np.r_[np.ones(len(cols)),0.][None,:],b_eq=[1.],bounds=[(0,None)]*len(cols)+[(None,None)],method='highs')
    require(sol.success,sol.message)
    if sol.x[-1]<=1e-9:w=exact_dual(rows,-sol.ineqlin.marginals);kind='dual'
    if w is None:
     tree,failed_branch=split_certificate(p,keys,signed,ids)
     if tree is not None:
      leaves=tree_leaves(tree);full=set()
      for leaf in leaves:
       weights=leaf.pop('weights');full.update(weights)
       leaf['terms']=[[*(inc[i//2] if i%2==0 else inc[i//2][::-1]),weight] for i,weight in weights.items()]
      hit=len(cache);cache.append((full,{}));rec['certificates'].append({'kind':'slack_one_tree','tree':tree})
      hist['slack_one_tree']+=1;stream.update(json.dumps([o,hit]).encode());continue
     support=np.flatnonzero(sol.x[:-1]>1e-9)
     for cap in (1000,1000000,1000000000):
      mix={int(j):Fraction(float(sol.x[j])).limit_denominator(cap) for j in support};mass=sum(mix.values());mix={j:w/mass for j,w in mix.items()}
      vals=[sum(w*Fraction(int(row[j]),int(E[j])) for j,w in mix.items()) for row in rows]
      if min(vals)>0:
       return {'obstruction':{'failed_slack_branch':failed_branch,'order':o,'mixture':[{'prefix':keys[j][0],'suffix':keys[j][1],'E':int(E[j]),'weight':str(w)} for j,w in mix.items()],'forward_pairs':[inc[i] if r[inc[i][0]]<r[inc[i][1]] else inc[i][::-1] for i in range(len(inc))],'margins':list(map(str,vals))},**rec}
     raise RuntimeError('uncertified LP')
   full={ids[i]:v for i,v in w.items()};vals=[sum(v*int(signed[i,j]) for i,v in full.items()) for j in range(len(cols))];require(max(vals)<=0,'exact replay')
   hit=len(cache);cache.append((set(full),full));rec['certificates'].append({'terms':[[*(inc[i//2] if i%2==0 else inc[i//2][::-1]),w] for i,w in full.items()],'kind':kind,'max_margin':max(vals)})
  hist[rec['certificates'][hit]['kind']]+=1;stream.update(json.dumps([o,hit]).encode())
 rec.update(hist=dict(hist),case_digest=stream.hexdigest(),generators=len(cols),extensions=len(words));return rec
