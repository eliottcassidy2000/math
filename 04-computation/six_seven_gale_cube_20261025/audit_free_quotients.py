#!/usr/bin/env python3
"""Independent incidence/code audit; does not import producer checker."""
from collections import Counter
from itertools import combinations,permutations
import json
N=0

def check(ok,name):
 global N
 N+=1
 if not ok: raise RuntimeError(name)

T=list(combinations(range(7),3))
F=set()
for p in permutations(range(1,8)):
 F.add(frozenset(t for t in T if p[t[0]]^p[t[1]]^p[t[2]]==0))
check(len(F)==30,'thirty Fano systems')
actual=frozenset([(0,1,2),(0,3,4),(0,5,6)])
bad=Counter(len(f&actual) for f in F)
check(bad=={0:16,1:12,3:2},'hourglass bad-triple census')
pairs=[(f,g) for f,g in combinations(F,2) if not f&g]
check(len(pairs)==120,'disjoint Fano pairs')
pairc=Counter(tuple(sorted((len(f&actual),len(g&actual)))) for f,g in pairs)
check(pairc=={(0,0):32,(0,1):48,(0,3):16,(1,1):24},'complementary kernel count')
# Use 7-bit even representative of each projective sign class.
def canonical(mask):
 return mask if mask.bit_count()%2==0 else mask^127

def kernel(f):
 return frozenset([0]+[canonical(sum(1<<i for i in t)) for t in f])
for f in F:
 h=kernel(f)
 check(len(h)==8 and all(a^b in h for a in h for b in h),'kernel linear subgroup')
 inertia=[]
 for t in actual:
  generators=[canonical(1<<i) for i in t]
  elements={0}
  for g in generators: elements|={e^g for e in list(elements)}
  inertia.append(elements)
 check(all(h.intersection(I)=={0} for I in inertia)==(not f&actual),'direct inertia equivalence')
for f,g in pairs:
 H,H1=kernel(f),kernel(g)
 check(H&H1=={0},'complementary kernel intersection')
 check(len({a^b for a in H for b in H1})==64,'complementary span')
# Six-root standard: enumerate 10 unordered 3+3 partitions without code names.
Tsix=set([(0,1,2),(0,3,4)])
parts=[frozenset((frozenset(t),frozenset(set(range(6))-set(t)))) for t in combinations(range(6),3)]
parts=set(parts)
free=[P for P in parts if not any(frozenset(t) in P for t in Tsix)]
check(len(parts)==10 and len(free)==8,'six-root free involutions')
# Every sum of two distinct free size3 representatives has size2 modulo complement.
reps=[sum(1<<j for j in min((tuple(sorted(s)) for s in P))) for P in free]
for a,b in combinations(reps,2):
 c=a^b
 check(min(c.bit_count(),6-c.bit_count())==2,'no larger free six-root subgroup')
print(json.dumps({'status':'PASS','checks':N,'Fano_systems':30,'hourglass_free_order8':16,'complementary_pairs':120,'both_free_pairs':32,'six_root_free_involutions':8},sort_keys=True,indent=2))
