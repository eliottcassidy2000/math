#!/usr/bin/env python3
"""Exact independent vertex-level controls for the weighted-chain arc-edit laws."""
from itertools import combinations, product
from collections import Counter, defaultdict
from math import comb
import json
from pathlib import Path

checks=0
def check(p,msg='check failed'):
 global checks
 checks+=1
 if not p: raise AssertionError(msg)

def rows(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else (j,i)
  a[u]|=1<<v
 return tuple(a)
def h(a,u,v):
 return ((a[u]^a[v])&~((1<<u)|(1<<v))).bit_count()
def chains(a):
 n=len(a); adj=[[] for _ in a]
 for u,v in combinations(range(n),2):
  if h(a,u,v)==0: adj[u].append(v);adj[v].append(u)
 seen=set(); ans=[]
 for u in range(n):
  if u in seen: continue
  stack=[u]; component=[];seen.add(u)
  while stack:
   v=stack.pop();component.append(v)
   for w in adj[v]:
    if w not in seen: seen.add(w);stack.append(w)
  # Internal outdegree strictly orders each transitive component.
  S=sum(1<<v for v in component)
  component.sort(key=lambda v:-(a[v]&S).bit_count())
  ans.append(component)
 return ans
def flip(a,u,v):
 out=list(a);out[u]^=1<<v;out[v]^=1<<u;return tuple(out)
def core(a,cs):
 return tuple(sum(1<<j for j,cj in enumerate(cs) if i!=j and a[ci[0]]>>cj[0]&1) for i,ci in enumerate(cs))
def predicted_h(q,ls):
 out=Counter()
 for l in ls:
  for d in range(1,l): out[d-1]+=l-d
 for i,j in combinations(range(len(q)),2):
  D=sum(ls[k] for k in range(len(q)) if k not in (i,j) and ((q[i]^q[j])>>k&1))
  for s in range(ls[i]):
   for t in range(ls[j]): out[D+s+t]+=1
 return +out
def blow(q,ls):
 ids=[];s=0
 for l in ls: ids.append(list(range(s,s+l)));s+=l
 a=[0]*s
 for i,C in enumerate(ids):
  for x,u in enumerate(C):
   a[u]|=sum(1<<v for v in C[x+1:])
   for j,D in enumerate(ids):
    if q[i]>>j&1: a[u]|=sum(1<<v for v in D)
 return tuple(a),ids

def internal_jump(l):
 return +Counter({0:max(l-1,0),2:int(l>=3),3:2*max(l-3,0),4:comb(l-3,2) if l>=5 else 0})
def endpoint_types(l):
 if l==1: return [(0,1)]
 if l==2: return [(0,1),(1,1)]
 return [(0,1),(1,l-2),(2,1)]
def template_jump(q,ls):
 """Rank-transition histogram, no construction of the full weighted tournament."""
 out=Counter();b=len(q)
 for l in ls: out.update(internal_jump(l))
 for i,j in combinations(range(b),2):
  lens=[1]*b;lens[i]=min(ls[i],3);lens[j]=min(ls[j],3)
  small,ids=blow(q,lens)
  for a,ca in endpoint_types(ls[i]):
   for c,cc in endpoint_types(ls[j]):
    delta=len(chains(flip(small,ids[i][a],ids[j][c])))-b
    out[delta]+=ca*cc
 return +out

def run():
 by_nr=defaultdict(dict);by_hist={};hostile=None; histogram_hostile=None; totals={}
 for n in range(1,6):
  for code in range(1<<comb(n,2)):
   a=rows(n,code);cs=chains(a);b=len(cs);ls=list(map(len,cs));q=core(a,cs)
   hist=Counter(h(a,u,v) for u,v in combinations(range(n),2))
   check(hist==predicted_h(q,ls),'separation polynomial')
   check(hist[0]==n-b,'forest rank')
   e1pred=sum(max(l-2,0) for l in ls)
   for i,j in combinations(range(b),2):
    D=sum(ls[k] for k in range(b) if k not in (i,j) and (q[i]^q[j])>>k&1)
    check(D>=1,'core is pair-free')
    e1pred+=D==1
   check(hist[1]==e1pred,'singleton separator formula')
   jump=Counter()
   for u,v in combinations(range(n),2):
    d=len(chains(flip(a,u,v)))-b;jump[d]+=1;check(abs(d)<=4,'rank jump')
   check(sum(d*c for d,c in jump.items())==2*((n-2)*(n-b)-hist[1]),'drift sum')
   check(jump==template_jump(q,ls),'all exact templates')
   key=(n,tuple(sorted(hist.items())))
   if histogram_hostile is None and key in by_hist and by_hist[key]['jump_histogram']!=dict(jump):
    histogram_hostile={'n':n,'pair_histogram':dict(hist),'first':by_hist[key],'second':{'rows':a,'jump_histogram':dict(jump)}}
   by_hist.setdefault(key,{'rows':a,'jump_histogram':dict(jump)})
   old=by_nr[n,b]
   if hostile is None and old and tuple(sorted(jump.items())) not in old:
    other=next(iter(old.values()));hostile={'n':n,'core_rank':b,'first':other,'second':{'rows':a,'e1':hist[1],'jump_histogram':dict(jump)}}
   old.setdefault(tuple(sorted(jump.items())),{'rows':a,'e1':hist[1],'jump_histogram':dict(jump)})
  totals[n]=1<<comb(n,2)
 # Positive controls with potentially large weights: all 4-vertex cores, triangle,
 # and explicit homogeneous-free five-vertex examples.
 Qs=[(2,4,1)]
 for n in (4,5):
  for code in range(1<<comb(n,2)):
   a=rows(n,code)
   if len(chains(a))==n and len(Qs)<12: Qs.append(a)
 weighted=0
 for q in Qs:
  for ls in ([3]*len(q),list(range(1,len(q)+1)),[1+(i%3) for i in range(len(q))]):
   a,cs=blow(q,ls);b=len(q);n=len(a);weighted+=1
   check(len(chains(a))==b,'blowup core')
   check(Counter(h(a,u,v) for u,v in combinations(range(n),2))==predicted_h(q,ls))
   for C in cs:
    internal=Counter()
    for x,y in combinations(range(len(C)),2):
     d=len(chains(flip(a,C[x],C[y])))-b
     expected=0 if y==x+1 else 2+(x>0)+(y<len(C)-1)
     check(d==expected,'internal flip')
     internal[d]+=1
    l=len(C)
    pred=Counter({0:max(l-1,0),2:int(l>=3),3:2*max(l-3,0),4:comb(l-3,2) if l>=5 else 0})
    check(+internal==+pred,'internal polynomial')
   for C,D in combinations(cs,2):
    for u in C[1:-1]:
     for v in D[1:-1]: check(len(chains(flip(a,u,v)))-b==4,'interior cross flip')
   jump=Counter(len(chains(flip(a,u,v)))-b for u,v in combinations(range(n),2))
   check(jump==template_jump(q,ls),'weighted exact templates')
   e1=predicted_h(q,ls)[1]
   check(sum(d*c for d,c in jump.items())==2*((n-2)*(n-b)-e1))
   I=sum(max(l-2,0) for l in ls)
   guaranteed=comb(I,2)-sum(max(l-3,0) for l in ls)
   check(jump[4]>=guaranteed,'interior bulk')
 # Sharpness of absolute bound.
 tt=tuple(sum(1<<j for j in range(i+1,5)) for i in range(5))
 sharp=flip(tt,1,3);check(len(chains(tt))==1 and len(chains(sharp))==5)
 large_weights=[10**6,10**6+1,10**6+2]
 large=template_jump((2,4,1),large_weights)
 check(sum(large.values())==comb(sum(large_weights),2),'large symbolic total')
 out={'large_weight_example':{'core_rows':[2,4,1],'weights':large_weights,'jump_histogram':dict(large)},'status':'FINITE-EXACT controls; proofs in PROOF.md','labelled_universes':totals,'weighted_cases':weighted,'checks':checks,'rank_jump_sharp':{'before':tt,'after':sharp},'non_lumpability':hostile,'histogram_non_lumpability':histogram_hostile}
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':run()
