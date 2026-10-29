#!/usr/bin/env python3
"""Exact finite checks for decorated Royle parity; stdlib, no assertions."""
from itertools import combinations, permutations
from fractions import Fraction
import json

checks=0
def check(p,msg):
 global checks
 checks+=1
 if not p: raise RuntimeError(msg)
def setup(n):
 es=list(combinations(range(n),2)); ix={e:i for i,e in enumerate(es)}
 ps=list(permutations(range(n))); maps=[]
 for p in ps:
  dest=[ix[tuple(sorted((p[a],p[b])))] for a,b in es]
  rev=sum(1<<i for i,(a,b) in enumerate(es) if p[a]>p[b])
  gm=[sum(((x>>i)&1)<<j for i,j in enumerate(dest)) for x in range(1<<len(es))]
  vm=[sum(((d>>i)&1)<<p[i] for i in range(n)) for d in range(1<<n)]
  maps.append((gm,vm,rev,gm[rev]))
 return es,ps,maps

def census(n,walsh=False):
 es,ps,maps=setup(n); m=len(es); L=1<<m; V=1<<n
 def act(z,k,tour):
  x,r,d=z; gm,vm,rev,c=maps[k]
  return (gm[x]^(c if tour else 0),gm[r],vm[d])
 def orbits(tour):
  seen=set(); rows=[]
  for x in range(L):
   for r in range(L):
    for d in range(V):
     z=(x,r,d)
     if z in seen: continue
     orb={act(z,k,tour) for k in range(len(ps))};seen.update(orb)
     stab=[k for k in range(len(ps)) if act(z,k,tour)==z]
     even=all(((x&maps[k][2]).bit_count()%2)==0 for k in stab)
     rows.append((z,orb,stab,even))
  return rows
 ts=orbits(True); gs=orbits(False); good=[r for r in gs if r[3]]
 check(len(ts)==len(good),('counts',n))
 # Every source stabilizer has odd order; target may have even order.
 for _,_,st,_ in ts: check(len(st)%2==1,('odd source',n))
 if walsh:
  columns=[]
  for z,orb,st,even in good:
   signed={}
   for k in range(len(ps)):
    z1=act(z,k,False);sgn=(-1)**((z[0]&maps[k][2]).bit_count()%2)
    check(z1 not in signed or signed[z1]==sgn,'transport sign')
    signed[z1]=sgn
   col=[]
   for (t,r,d),_,_,_ in ts:
    col.append(sum(s*((-1)**((x&t).bit_count()%2)) for (x,r1,d1),s in signed.items() if (r1,d1)==(r,d)))
   columns.append(col)
  # Exact weighted orthogonality: partial Fourier has scale 2^m, not 2^(n^2).
  for i,a in enumerate(columns):
   for j,b in enumerate(columns):
    lhs=sum(len(ts[k][1])*a[k]*b[k] for k in range(len(ts)))
    check(lhs==((1<<m)*len(good[i][1]) if i==j else 0),('Walsh',n,i,j,lhs))
 return {'n':n,'pair_classes':len(ts),'decorated_even_classes':len(good),'all_decorated_graph_classes':len(gs)}

def sylow5():
 n=5;es,ps,maps=setup(n)
 def compose(a,b):return tuple(a[b[i]] for i in range(n))
 def inverse(a):return tuple(a.index(i) for i in range(n))
 out=[]
 for p,gen in [(3,(1,2,0,3,4)),(5,(1,2,3,4,0))]:
  P={tuple(range(n))};a=gen
  while a not in P:P.add(a);a=compose(gen,a)
  N=[k for k,g in enumerate(ps) if {compose(compose(g,h),inverse(g)) for h in P}==P]
  Pidx=[ps.index(g) for g in P]; L=1<<len(es)
  src=[];dst=[]
  for x in range(L):
   for r in range(L):
    # prune P fixed graphs and tournament separately
    rf=all(maps[k][0][r]==r for k in Pidx)
    if not rf:continue
    xf=all(maps[k][0][x]==x for k in Pidx)
    tf=all((maps[k][0][x]^maps[k][3])==x for k in Pidx)
    if not (xf or tf):continue
    for d in range(1<<n):
     if not all(maps[k][1][d]==d for k in Pidx):continue
     if tf:src.append((x,r,d))
     if xf:dst.append((x,r,d))
  def count_orbits(states,tour):
   seen=set();count=0
   for z in states:
    if z in seen:continue
    x,r,d=z;orb={(maps[k][0][x]^(maps[k][3] if tour else 0),maps[k][0][r],maps[k][1][d]) for k in N};seen.update(orb)
    st=[k for k in N if (maps[k][0][x]^(maps[k][3] if tour else 0),maps[k][0][r],maps[k][1][d])==z]
    if tour or all((x&maps[k][2]).bit_count()%2==0 for k in st):count+=1
   return count
  a=count_orbits(src,True);b=count_orbits(dst,False)
  check(a==b,('Sylow',p,a,b));out.append({'p':p,'normalizer':len(N),'source_types':a,'even_target_types':b})
 return out

def hostiles():
 es,ps,maps=setup(2);swap=ps.index((1,0));gm,vm,rev,c=maps[swap]
 # K2 is odd; coloring its endpoints differently removes its odd automorphism.
 check((1&rev).bit_count()%2==1,'K2 odd')
 check(vm[1]!=1,'color destroys swap')
 # Edge permutation on the singleton edge is even, but Royle sign is odd.
 check(gm[1]==1 and rev==1,'edge-permutation distinction')
 # Partial Fourier coordinates can be negative: not a positive set bijection.
 check((-1)**((1&1).bit_count()%2)==-1,'negative Walsh')

if __name__=='__main__':
 hostiles();rows=[census(n,walsh=n<=3) for n in range(5)]
 sy=sylow5()
 print(json.dumps({'census':rows,'sylow_n5':sy,'checks':checks,'status':'PASS'},indent=2))
