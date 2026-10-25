#!/usr/bin/env python3
"""Exact torus/Fano/hourglass comparison; no external files or float predicates."""
from itertools import combinations, permutations
from collections import Counter
import json
import sympy as sp
checks=0
def need(c,msg):
 global checks
 checks+=1
 if not c: raise RuntimeError(msg)
def design(ts):return frozenset(frozenset(t) for t in ts)
def image(f,p):return design({p[x] for x in t} for t in f)
def rows(f):return sorted(map(sorted,f))
U=design({t,(t+1)%7,(t+3)%7} for t in range(7))
D=design({t,(t+2)%7,(t+3)%7} for t in range(7))
faces=U|D; edges=list(combinations(range(7),2))
for name,f in [('U',U),('D',D)]:
 need(len(f)==7,name+' seven triples')
 count=Counter(p for t in f for p in combinations(sorted(t),2))
 need(set(count)==set(edges) and set(count.values())=={1},name+' STS')
need(not U&D,'two disjoint systems')
counts=Counter(p for t in faces for p in combinations(sorted(t),2))
need(len(faces)==14 and set(counts.values())=={2},'closed complex')
for v in range(7):
 link=[tuple(sorted(t-{v})) for t in faces if v in t]
 adj={w:set() for w in range(7) if w!=v}
 for a,b in link: adj[a].add(b);adj[b].add(a)
 need(all(len(ns)==2 for ns in adj.values()),'link degrees')
 seen=set();todo=[next(iter(adj))]
 while todo:
  x=todo.pop()
  if x not in seen:seen.add(x);todo.extend(adj[x]-seen)
 need(len(seen)==6,'link is one cycle')
B1=sp.zeros(7,21);B2=sp.zeros(21,14)
for j,(a,b) in enumerate(edges):B1[a,j]=-1;B1[b,j]=1
for j,(a,b,c) in enumerate(sorted(tuple(sorted(t)) for t in faces)):
 for e,s in [((a,b),1),((a,c),-1),((b,c),1)]:B2[edges.index(e),j]=s
need(B1*B2==sp.zeros(7,14),'chain complex')
need((B1.rank(),B2.rank())==(6,13),'rational Betti numbers 1,2,1')
orientation=B2.nullspace()[0]
need(set(abs(x) for x in orientation)=={1},'integral fundamental orientation')
# All Fano markings and all designs, independent of cyclic difference set test.
base=design(t for t in combinations(range(7),3) if (t[0]+1)^(t[1]+1)^(t[2]+1)==0)
all_fano={image(base,p) for p in permutations(range(7))}
need(len(all_fano)==30,'thirty Fano structures')
for f,vs in [(U,[1,2,4,3,6,7,5]),(D,[1,2,4,5,7,3,6])]:
 need(f in all_fano,'cyclic design is Fano')
 need(set(vs)==set(range(1,8)),'nonzero binary bijection')
 need(design(t for t in combinations(range(7),3) if vs[t[0]]^vs[t[1]]^vs[t[2]]==0)==f,'explicit marking')
 # even-weight representative for projective signs, n odd
 ker=[m for m in range(128) if m.bit_count()%2==0 and not __import__('functools').reduce(int.__xor__,(vs[i] for i in range(7) if m>>i&1),0)]
 need(len(ker)==8 and all(m==0 or m.bit_count()==4 for m in ker),'simplex sign kernel')
R=[(0,0),(1,0),(0,1),(1,1),(-1,0),(0,-1),(-1,-1)]
c=lambda v,k=2:(v[0]+k*v[1])%7
rot=lambda v:(v[0]-v[1],v[0])
ref=lambda v:(v[0]-v[1],-v[1])
group=[]
for flip in (False,True):
 for k in range(6):
  p={}
  for v in R:
   w=ref(v) if flip else v
   for _ in range(k):w=rot(w)
   p[c(v)]=c(w)
  group.append(tuple(p[i] for i in range(7)))
need(len(set(group))==12,'D12 size')
need(all(tuple(p[q[i]] for i in range(7)) in group for p in group for q in group),'group closure')
star=design([{0,1,6},{0,2,5},{0,3,4}])
need(all(image(star,p)==star for p in group),'actual crossing invariance')
orbit={image(U,p) for p in group}
need(len(orbit)==4 and D in orbit,'four orbit containing both faces')
need(all(not f&star for f in orbit),'all four codes free at crossings')
# Independently produce both chiral quotient face systems on fixed c2 labels.
constructions=[]
for k in (2,4):
 # inverse c_k, then re-express on c_2 of the same branch seven-set
 p={c(v,k):c(v,2) for v in R}
 up=design({t,(t+1)%7,(t+1+k)%7} for t in range(7))
 down=design({t,(t+k)%7,(t+1+k)%7} for t in range(7))
 constructions.extend([image(up,p),image(down,p)])
need(len(set(constructions))==4 and set(constructions)==orbit,'two ideals x two faces exact four')
remaining=set(all_fano); census=[]
while remaining:
 f=min(remaining,key=lambda f:rows(f));o={image(f,p) for p in group}
 need(o<=remaining,'orbits partition')
 census.append((len(f&star),len(o)));remaining-=o
need(sorted(census)==[(0,4),(0,12),(1,6),(1,6),(3,2)],'complete orbit census')
need(Counter(len(f&star) for f in all_fano)=={0:16,1:12,3:2},'crossing type census')
stab=[p for p in group if image(U,p)==U]
need(len(stab)==3,'stabilizer order3')
inv=tuple((-i)%7 for i in range(7))
need(image(U,inv)==D and not any(image(f,inv)==f for f in all_fano),'inversion hostile')
# Source lattice is 3-colored, its K7 quotient is 7-colored.
units=R[1:]
need(all((a+b)%3!=0 for a,b in units),'source proper 3-color')
need(set(c(v) for v in units)==set(range(1,7)),'quotient all differences K7')
# Pure Euclidean concurrency is not this Fano structure.
need(star.isdisjoint(U) and star.isdisjoint(D),'incidence changes rather than preserves')
print(json.dumps({'status':'PASS','checks':checks,'vertices':7,'edges':21,'faces':14,'betti':[1,2,1],'up':rows(U),'down':rows(D),'orbit_size':len(orbit),'orbit':sorted((rows(f) for f in orbit)),'all_code_orbits':sorted(census),'automorphism_consequence':'C2 x A4 via inherited Session60 classification','scope':'no Euclidean plane lower bound'},indent=2))
