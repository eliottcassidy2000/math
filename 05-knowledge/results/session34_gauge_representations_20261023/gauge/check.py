#!/usr/bin/env python3
"""Exact independent small controls. No producer imports, floating point or assert."""
from itertools import combinations,permutations,product
from collections import Counter,defaultdict
from math import comb,factorial
import json
checks=0
def eq(a,b,label):
 global checks
 checks+=1
 if a!=b:raise RuntimeError((label,a,b))
def xor(xs):
 z=0
 for x in xs:z^=x
 return z
def parity(x):return x.bit_count()%2
def hamilton(n,mask,edges):
 idx={e:k for k,e in enumerate(edges)}
 dp={(1<<v,v):1 for v in range(n)}
 for s in range(1,1<<n):
  for v in range(n):
   c=dp.get((s,v),0)
   if not c:continue
   for w in range(n):
    if (s>>w)&1:continue
    bit=(mask>>idx[tuple(sorted((v,w)))])&1
    if bool(v<w)!=bool(bit):dp[s|1<<w,w]=dp.get((s|1<<w,w),0)+c
 return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
small={}
for n in range(2,6):
 E=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(E)}
 F=list(combinations(range(n),3)); face=[sum(1<<idx[e] for e in combinations(f,2)) for f in F]
 hist=Counter(); fluxfib=Counter(); gavg=Counter()
 for u in range(1<<len(E)):
  flux=sum(parity(u&f)<<k for k,f in enumerate(face));fluxfib[flux]+=1
  hist[flux.bit_count()]+=1
  for tet in combinations(range(n),4):
   eq(sum((flux>>F.index(f))&1 for f in combinations(tet,3))%2,0,'Bianchi')
  gauge=sum(((u>>idx[(0,v)])&1)<<v for v in range(1,n))
  fixed=u^sum((((gauge>>a)^(gauge>>b))&1)<<k for k,(a,b) in enumerate(E))
  reconstructed=sum(((flux>>F.index((0,a,b)))&1)<<k for k,(a,b) in enumerate(E) if a>0)
  eq(fixed,reconstructed,'root gauge reconstruction')
  gavg[fixed]+=hamilton(n,u,E)
 eq(set(fluxfib.values()),{1<<(n-1)},'gauge fibers')
 eq(set(gavg.values()),{factorial(n)},'Hamilton gauge average')
 closed=Counter(); fillings=Counter(); C=0
 if n>=4:C=sum(1<<idx[e] for e in [(0,1),(1,2),(2,3),(0,3)])
 for s in range(1<<len(F)):
  bd=xor(face[k] for k in range(len(F)) if (s>>k)&1)
  if bd==0:closed[s.bit_count()]+=1
  if n>=4 and bd==C:fillings[s.bit_count()]+=1
 # independent formal high temperature identity
 lhs=[0]*(len(F)+1)
 for k,h in hist.items():
  for a in range(len(F)-k+1):
   for b in range(k+1):lhs[a+b]+=h*comb(len(F)-k,a)*comb(k,b)*(-1)**b
 eq(lhs,[(1<<len(E))*closed[j] for j in range(len(F)+1)],'partition polynomial')
 small[n]={'configurations':1<<len(E),'gauge_classes':len(fluxfib),'closed_surface_weights':dict(sorted(closed.items())),'C4_filling_weights':dict(sorted(fillings.items()))}
# Six vertices: all closed surfaces from root tetrahedron basis
n=6;F=list(combinations(range(n),3));tet=[sum(1<<F.index(f) for f in combinations((0,)+abc,3)) for abc in combinations(range(1,n),3)]
surf={xor(tet[i] for i in range(len(tet)) if (s>>i)&1) for s in range(1<<len(tet))}
eq(len(surf),1024,'independent tetra boundaries');small[6]={'closed_surface_weights':dict(sorted(Counter(s.bit_count() for s in surf).items()))}
# Fano rainbow Wilson theorem
n=7;E=list(combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
L=[{i,(i+1)%7,(i+3)%7} for i in range(7)];A=[4,1,3,7,6,5,2]
el={e:next(i for i,l in enumerate(L) if set(e)<=l) for e in E}
line=[sum(1<<idx[e] for e in E if el[e]==i) for i in range(7)]
base=sum((0 if (b-a)%7 in (1,2,4) else 1)<<i for i,(a,b) in enumerate(E))
dual={sum(parity(b&a)<<i for i,a in enumerate(A)) for b in range(1,8)}
cycles=[];rainbow=defaultdict(list)
for vs in combinations(range(n),4):
 for tail in permutations(vs[1:]):
  c=(vs[0],)+tail
  if c[1]>c[-1]:continue
  edges=[tuple(sorted((c[i],c[(i+1)%4]))) for i in range(4)]
  em=sum(1<<idx[e] for e in edges); lm=xor(1<<el[e] for e in edges)
  cycles.append((em,lm))
  if len({el[e] for e in edges})==4:rainbow[lm].append((c,em))
eq(set(rainbow),dual,'exact dual support');eq(set(map(len,rainbow.values())),{3},'three cycles per support')
fano_hist=Counter();saturation=set();fullH_disagreement=None
for m in range(128):
 delta=xor(line[i] for i in range(7) if m>>i&1);u=base^delta
 h=hamilton(7,u,E);fano_hist[h]+=1
 w=sum((-1)**parity(delta&em) for rows in rainbow.values() for c,em in rows)
 eq(4*h,693+3*w,'Fano Hamilton Wilson identity')
 for g in range(64):
  cut=sum((((g>>a)^(g>>b))&1)<<k for k,(a,b) in enumerate(E))
  v=u^cut;saturation.add(v)
  eq(sum((-1)**parity((v^base)&em) for rows in rainbow.values() for c,em in rows),w,'Wilson switching invariant')
  if m==0 and fullH_disagreement is None and hamilton(7,v,E)!=h:fullH_disagreement={'base_H':h,'switched_H':hamilton(7,v,E),'gauge_bits':g}
eq(len(saturation),8192,'Fano transverse saturation')
# Nonabelian S3 tetrahedral Bianchi and hostile to unordered abelian version
S=list(permutations(range(3)));one=(0,1,2)
def mul(a,b):return tuple(a[b[i]] for i in range(3))
def inv(a):return tuple(a.index(i) for i in range(3))
def prod(*xs):
 z=one
 for x in xs:z=mul(z,x)
 return z
nonabelian_wrong=0;flat=0
for xs in product(S,repeat=6):
 u=dict(zip(combinations(range(4),2),xs))
 def U(a,b):return u[a,b] if a<b else inv(u[b,a])
 def F(a,b,c):return prod(U(a,b),U(b,c),U(c,a))
 eq(prod(F(0,1,2),F(0,2,3)),prod(U(0,1),F(1,2,3),U(1,0),F(0,1,3)),'nonabelian Bianchi')
 nonabelian_wrong+=prod(F(0,1,2),F(0,2,3),inv(F(0,1,3)),inv(F(1,2,3)))!=one
 flat+=all(F(*f)==one for f in combinations(range(4),3))
eq(flat,216,'S3 flat count')
print(json.dumps({'status':'PASS','checks':checks,'small_simplex':small,'fano':{'H_histogram':dict(fano_hist),'rainbow_cycles':sum(map(len,rainbow.values())),'rainbow_supports':{str(k):[list(c) for c,e in rows] for k,rows in sorted(rainbow.items())},'saturation':len(saturation),'hostile_full_H':fullH_disagreement},'S3':{'configurations':6**6,'flat':flat,'naive_untransported_Bianchi_failures':nonabelian_wrong}},indent=2,sort_keys=True))
