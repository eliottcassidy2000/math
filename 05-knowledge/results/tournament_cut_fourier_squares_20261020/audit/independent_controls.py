#!/usr/bin/env python3
"""Independent direct-signing/forest audit; no producer imports or Walsh routine."""
from itertools import combinations,permutations
from collections import Counter,defaultdict
from math import factorial,comb
from pathlib import Path
import json,hashlib
checks=0
def check(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)
def setup(n):
 E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
 bit=lambda u,v:1<<idx[tuple(sorted((u,v)))]
 paths=[sum(bit(u,v) for u,v in zip(p,p[1:])) for p in permutations(range(n))]
 cuts={sum(1<<k for k,(u,v) in enumerate(E) if bool(s>>u&1)!=bool(s>>v&1)) for s in range(1<<n)}
 return E,bit,paths,cuts

def forest_coefficient(n,E,F):
 adj={v:set() for v in range(n)}
 for k,(u,v) in enumerate(E):
  if F>>k&1:adj[u].add(v);adj[v].add(u)
 if any(len(v)>2 for v in adj.values()):return 0
 unused=set(range(n));blocks=0
 while unused:
  root=unused.pop();visited={root};todo=[root]
  while todo:
   u=todo.pop()
   for v in adj[u]-visited:visited.add(v);unused.discard(v);todo.append(v)
  edges=sum(len(adj[u]) for u in visited)//2
  if edges!=len(visited)-1:return 0
  if edges:blocks+=1
 return (2**blocks)*factorial(n-F.bit_count())

out={'status':'PASS','orders':[],'all_signing_controls':[]}
for n in range(1,7):
 E,bit,paths,cuts=setup(n);N=len(E);full=(1<<N)-1
 canonical=lambda x:min(x^c for c in cuts)
 kernel=Counter(canonical(full^p) for p in paths)
 basis_cuts=[sum(bit(v,w) for w in range(n) if w!=v) for v in range(n)]
 dual=[Z for Z in range(1<<N) if all((Z&c).bit_count()%2==0 for c in basis_cuts)]
 check(len(dual)==2**(N-n+1),'dual dimension')
 spectrum=Counter();values={}
 for Z in dual:
  signed=sum(1 if (Z&p).bit_count()%2==0 else -1 for p in paths)
  eig=(-1)**Z.bit_count()*signed
  direct=sum(w*(1 if (Z&g).bit_count()%2==0 else -1) for g,w in kernel.items())
  check(eig==direct,'direct signing versus cut kernel')
  values[Z]=eig;spectrum[eig]+=1
 check(sum(k*v for k,v in spectrum.items())==len(dual)*kernel[0],'trace')
 check(sum(k*k*v for k,v in spectrum.items())==len(dual)*sum(w*w for w in kernel.values()),'second moment')
 if n>=3:
  check(kernel[0]==0,'no loops')
  nonzero=[g for g in kernel if kernel[g]]
  w=kernel[nonzero[0]]
  check({kernel[g] for g in nonzero}=={w},'uniform incidence multiplicity')
  check(sum(k*k*v for k,v in spectrum.items())==len(dual)*factorial(n)*w,'weighted second moment')
  check(spectrum[factorial(n)]==(2 if n%4==1 else 1),'top multiplicity')
  check(spectrum[-factorial(n)]==(1 if n%4==3 else 0),'bottom multiplicity')
  if n>=4:
   for g in nonzero:
    h=next(h for h in nonzero if h!=g)
    cyc=[0,g,g^h,h]
    check(len(set(cyc))==4 and all(canonical(cyc[i]^cyc[(i+1)%4]) in nonzero for i in range(4)),'square certificate')
  else:check(len(nonzero)==1 and w==6,'order-three support/multiplicity')
  transitive_paths=[p for p in permutations(range(n)) if all(u<v for u,v in zip(p,p[1:]))]
  check(len(transitive_paths)==1,'transitive leaf')
  check(full^sum(bit(u,v) for u,v in zip(transitive_paths[0],transitive_paths[0][1:]))!=0,'leaf edge not loop')
 else:check(kernel=={0:factorial(n)},'small loop boundary')
 if n==5:
  paired={min(Z,full^Z) for Z in dual if values[Z]==24}
  check(spectrum[24]==20 and len(paired)==10,'n5 full/component multiplicity')
 if n<=5:
  coeff={F:forest_coefficient(n,E,F) for F in range(1<<N)}
  for F,count in coeff.items():check(count==sum((p&F)==F for p in paths),'literal containing-word coefficient')
  for Z in range(1<<N):
   total=0;F=Z
   while True:
    total+=(-2)**F.bit_count()*coeff[F]
    if not F:break
    F=(F-1)&Z
   direct=sum(1 if (p&Z).bit_count()%2==0 else -1 for p in paths)
   check(total==direct,'all-signing forest expansion')
   if n>=2:
    # Each unoriented cycle with the positive apex has one representative
    # with the two apex neighbours in increasing order.
    cycle_sum=0
    for p in permutations(range(n)):
     if p[0]>p[-1]:continue
     pm=sum(bit(u,v) for u,v in zip(p,p[1:]))
     cycle_sum+=1 if (pm&Z).bit_count()%2==0 else -1
    check(direct==2*cycle_sum,'positive-apex factor two')
  out['all_signing_controls'].append({'n':n,'signings':1<<N,'containing_word_coefficients':1<<N})
 out['orders'].append({'n':n,'characters':len(dual),'spectrum':dict(sorted(spectrum.items()))})
# Triangle sum by the positions of its three distinguished vertices, independently
# of both the forest formula and the producer subset DP.
triangle=[]
for n in range(3,12):
 signed=6*factorial(n-3)*sum((-1)**sum((b==a+1) for a,b in zip(slots,slots[1:])) for slots in combinations(range(n),3))
 expected=-factorial(n-2)*(n-4)*(n-9)
 check(-signed==expected,'triangle marked-position count')
 triangle.append({'n':n,'triangle_eigenvalue':-signed})
out['triangle_modes']=triangle
# Pair restrictions on the even displacement hyperplane without enumerating Q9.
n=9;E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)};full=(1<<len(E))-1
triangles={sum(1<<idx[e] for e in combinations(vs,2)) for vs in combinations(range(n),3)}
zero_modes=triangles|{full^Z for Z in triangles}
check(len(triangles)==84 and len(zero_modes)==168,'n9 full lower bound')
check(all(all(sum(bool(Z>>idx[tuple(sorted((u,v)))]&1) for v in range(n) if v!=u)%2==0 for u in range(n)) for Z in zero_modes),'n9 Eulerian modes')
even_generators=[1^(1<<k) for k in range(1,len(E))]
signature=lambda Z:tuple((Z&g).bit_count()%2 for g in even_generators)
classes=defaultdict(list)
for Z in zero_modes:classes[signature(Z)].append(Z)
check(len(classes)==84 and all(len(v)==2 and v[0]^v[1]==full for v in classes.values()),'n9 component restrictions')
check(all(((Z&1).bit_count()%2)!=(((Z^full)&1).bit_count()%2) for Z in triangles),'complements differ by sign on odd component')
out['n9_full_nullity_lower_bound']=168;out['n9_each_component_nullity_lower_bound']=84
out['odd_component_note']='Complementary characters restrict to scalar negatives on the odd coset; they represent the same one-dimensional eigenspace.'
# Hostile to transferring signing-switch invariance from cycles to paths.
E,bit,paths,cuts=setup(3);star=bit(0,1)|bit(0,2)
switched=sum((-1)**(star&p).bit_count() for p in paths)
check(switched==-2 and len(paths)==6,'signing-switch hostile')
out['signing_switch_hostile']={'n':3,'empty_signed_path_sum':6,'vertex_cut_signed_path_sum':switched}
out['checks']=checks
out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print(json.dumps(out,indent=2,sort_keys=True))
