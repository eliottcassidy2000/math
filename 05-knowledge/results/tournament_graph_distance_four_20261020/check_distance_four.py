#!/usr/bin/env python3
"""Exact n=5 cut-aggregate model. No producer imports; checks survive -O."""
from itertools import combinations, permutations
from collections import Counter
from functools import reduce
from pathlib import Path
import hashlib, json
checks=0
def require(condition,message):
 global checks
 checks+=1
 if not condition: raise RuntimeError(message)
E=tuple(combinations(range(5),2)); index={e:i for i,e in enumerate(E)}
F=tuple(combinations(range(1,5),2)); full=(1<<10)-1
rows=(13,14,19,22,26,28)
def edge(u,v): return 1<<index[tuple(sorted((u,v)))]
def project(x):
 return sum(1<<k for k,(a,b) in enumerate(F) if (bool(x&edge(a,b))^bool(x&edge(0,a))^bool(x&edge(0,b))))
def lift(q): return sum(edge(a,b) for k,(a,b) in enumerate(F) if q>>k&1)
def encode(q): return sum(((q&r).bit_count()%2)<<i for i,r in enumerate(rows))
def inverse(y):
 z=[(y>>i)&1 for i in range(6)]
 a,b,c,d,e,f=z
 return ((a^d^e)<<0)|((b^d^e)<<1)|((a^c^e)<<2)|((a^c^d)<<3)|((a^b^c)<<4)|((c^d^e)<<5)
def permute(x,p): return sum(edge(p[a],p[b]) for i,(a,b) in enumerate(E) if x>>i&1)
def cycle_from_row(row):
 x=lift(row)
 for v in range(1,5):
  if sum(bool(x&edge(v,w)) for w in range(1,5) if w!=v)%2: x ^=edge(0,v)
 return x
cycles=tuple(map(cycle_from_row,rows))
pairs=tuple(frozenset((c,c^full)) for c in cycles)
require(len(set(pairs))==6,'cycle pairs not distinct')
for C in cycles:
 require(C.bit_count()==5,'cycle edge count')
 require(all(sum(bool(C&edge(v,w)) for w in range(5) if w!=v)==2 for v in range(5)),'cycle degree')
# The actual path-complement displacement generators.
generators=Counter()
for p in permutations(range(5)):
 support=full
 for a,b in zip(p,p[1:]): support ^=edge(a,b)
 generators[project(support)]+=1
require(len(generators)==15 and set(generators.values())=={8},'generator count/multiplicity')
V=tuple(q for q in range(64) if q.bit_count()%2==0)
W2={q for q in V if q.bit_count()==2}; W4={q for q in V if q.bit_count()==4}
require({encode(g) for g in generators}==W4,'generator image is not exactly weight four')
require({encode(q) for q in V}==set(V),'linear map not onto')
for q in V:
 require(inverse(encode(q))==q and encode(inverse(q))==q,'explicit inverse')
 require(all(((lift(q)&C).bit_count()%2)==((q&r).bit_count()%2) for C,r in zip(cycles,rows)),'cycle character')
 for r in V:
  require(((q^r) in generators)==((encode(q)^encode(r)).bit_count()==4),'adjacency equivalence')
  require(encode(q^r)==encode(q)^encode(r),'linearity')
# Fourier characters of V: dual masks modulo the all-one vector.
characters=tuple(q for q in range(64) if q<(q^63))
spectrum=Counter(sum(1 if (h&g).bit_count()%2==0 else -1 for g in generators) for h in characters)
halfcube_spectrum=Counter(sum(1 if (h&g).bit_count()%2==0 else -1 for g in W2) for h in characters)
require(spectrum=={15:1,3:10,-1:15,-5:6},'actual spectrum')
require(halfcube_spectrum=={15:1,5:6,-1:15,-3:10},'halfcube spectrum')
require(spectrum!=halfcube_spectrum,'spectral hostile failed')
local_degrees=Counter(sum(a^b in generators for b in generators) for a in generators)
require(local_degrees=={6:15},'local KG degree')
# S5 acts on complementary Hamilton-five-cycle pairs, not as a fixed point in S6.
actions=set()
for p in permutations(range(5)):
 action=tuple(pairs.index(frozenset(permute(C,p) for C in pair)) for pair in pairs)
 actions.add(action)
 for q in V:
  target=sum(((encode(q)>>i)&1)<<action[i] for i in range(6))
  require(encode(project(permute(lift(q),p)))==target,'six-point equivariance on displacement space')
require(len(actions)==120,'six-point action not faithful')
require({a[0] for a in actions}==set(range(6)),'six-point action not transitive')
require({(a[0],a[1]) for a in actions}=={(i,j) for i in range(6) for j in range(6) if i!=j},'six-point action not doubly transitive')
# Clique obstruction independently verifies no local triangle larger than size3.
maximum_local_clique=0
for k in range(1,6):
 found=any(all(a^b in generators for a,b in combinations(subset,2)) for subset in combinations(generators,k))
 if found: maximum_local_clique=k
require(maximum_local_clique==3,'maximum clique is not four')
seed=(0,1,2,3,4)
path=sum(edge(a,b) for a,b in zip(seed,seed[1:]))
intersections=[(path&C).bit_count() for C in cycles]
require(intersections==[2,1,1,0,2,2],'seed cycle incidence table')
receipt={'status':'PASS','checks':checks,'component_vertices':32,'component_degree':15,
 'aggregate_edge_weight':8,'model':'Hamming-distance-four graph on even-weight six-bit vectors',
 'not_model':'halved six-cube, which uses Hamming distance two','parity_rows':rows,
 'generator_masks':dict(sorted(generators.items())),'generator_image_masks':sorted(encode(g) for g in generators),
 'spectrum':dict(sorted(spectrum.items())),'halved_cube_spectrum':dict(sorted(halfcube_spectrum.items())),
 'local_graph':'Kneser KG(6,2)','maximum_clique':4,'cycle_masks':cycles,'seed_cycle_intersections':intersections,
 'six_point_action_order':len(actions),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(receipt,indent=2,sort_keys=True))
