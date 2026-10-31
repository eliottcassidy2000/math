#!/usr/bin/env python3
"""Exact finite controls for AffRel(F2), including empty relations."""
import itertools,json
from functools import lru_cache
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise ValueError(msg)
@lru_cache(None)
def affine_sets(d):
 # Enumerate subspaces by adjoining one vector; then all their distinct cosets.
 zero=frozenset([0]);spaces={zero};todo=[zero]
 while todo:
  S=todo.pop()
  for v in range(1<<d):
   if v not in S:
    T=S|frozenset(x^v for x in S)
    if T not in spaces:spaces.add(T);todo.append(T)
 result={frozenset()}
 for S in spaces:
  unused=set(range(1<<d))
  while unused:
   v=min(unused);T=frozenset(x^v for x in S);result.add(T);unused-=T
 return tuple(sorted(result,key=lambda S:(len(S),tuple(sorted(S)))))
def rows(n,m,R):
 return tuple(sum(1<<y for y in range(1<<m) if x|(y<<n) in R) for x in range(1<<n))
def points(n,m,R):return frozenset(x|(y<<n) for x,ys in enumerate(R) for y in range(1<<m) if ys>>y&1)
def compose(R,S):
 out=[]
 for choices in R:
  row=0
  while choices:
   bit=choices&-choices;choices-=bit;row|=S[bit.bit_length()-1]
  out.append(row)
 return tuple(out)
def identity(n):return tuple(1<<x for x in range(1<<n))
def dagger(n,m,R):return tuple(sum(1<<x for x in range(1<<n) if R[x]>>y&1) for y in range(1<<m))
def tensor(n,m,R,p,q,S):
 return tuple(sum(1<<(y|(z<<m)) for y in range(1<<m) for z in range(1<<q) if R[x]>>y&1 and S[t]>>z&1) for t in range(1<<p) for x in range(1<<n))
relations={(n,m):[rows(n,m,R) for R in affine_sets(n+m)] for n in range(3) for m in range(3)}
composition_cases=0
for n,m,k in itertools.product(range(3),repeat=3):
 allowed=set(relations[n,k])
 for R in relations[n,m]:
  need(compose(identity(n),R)==R and compose(R,identity(m))==R,'identity laws')
  for S in relations[m,k]:
   T=compose(R,S);need(T in allowed,'affine composition closure')
   need(dagger(n,k,T)==compose(dagger(m,k,S),dagger(n,m,R)),'dagger reverses composition');composition_cases+=1
snake_cases=0
for n in range(3):
 cup=(sum(1<<(x|(x<<n)) for x in range(1<<n)),)
 cap=dagger(0,2*n,cup)
 left=compose(tensor(0,2*n,cup,n,n,identity(n)),tensor(n,n,identity(n),2*n,0,cap))
 right=compose(tensor(n,n,identity(n),0,2*n,cup),tensor(2*n,0,cap,n,n,identity(n)))
 need(left==identity(n) and right==identity(n),'two snake identities');snake_cases+=2
currying_cases=0
for u,v,w in itertools.product(range(3),repeat=3):
 # With low-bit input-first coordinate order, currying is literal reassociation.
 evaluation=tuple(1<<((x>>v)&((1<<w)-1)) if (x&((1<<v)-1))==(x>>(v+w)) else 0 for x in range(1<<(2*v+w)))
 for R in affine_sets(u+v+w):
  original=rows(u+v,w,R);curried=rows(u,v+w,R)
  rebuilt=compose(tensor(u,v+w,curried,v,v,identity(v)),evaluation)
  need(rebuilt==original,'evaluation after currying')
  need(points(u,v+w,curried)==points(u+v,w,original),'currying bijection inverse');currying_cases+=1
associativity_cases=0
for n,m,k,l in itertools.product(range(2),repeat=4):
 for R in relations[n,m]:
  for S in relations[m,k]:
   for T in relations[k,l]:
    need(compose(compose(R,S),T)==compose(R,compose(S,T)),'associativity');associativity_cases+=1
need(compose((1,),(0,1))==(0,),'inconsistent affine constraints compose to empty')
need(frozenset([1]) in affine_sets(1) and 0 not in frozenset([1]),'nonzero affine point is not linear')
# General Hamiltonian feasibility is not an affine condition: K4 has3 cycle edge-sets.
vertices=range(4);edges=list(itertools.combinations(vertices,2));cycles=set()
for tail in itertools.permutations([1,2,3]):
 P=(0,)+tail;chosen={tuple(sorted((P[i],P[(i+1)%4]))) for i in range(4)};cycles.add(sum(1<<edges.index(e) for e in chosen))
need(len(cycles)==3 and len(cycles)&(len(cycles)-1)!=0,'Hamiltonian cycle family not affine in general')
print(json.dumps({'status':'passed','checks':checks,'objects_tested_dimensions':[0,1,2],'composition_pairs':composition_cases,'snake_equations':snake_cases,'currying_relations':currying_cases,'associativity_triples_dimensions_at_most1':associativity_cases,'affine_set_counts_by_ambient_dimension':[len(affine_sets(d)) for d in range(7)],'K4_Hamiltonian_cycle_edge_sets':len(cycles)},indent=2))
