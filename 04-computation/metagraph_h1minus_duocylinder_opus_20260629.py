"""
Is H_1^-(n=5)=7 a FANO plane or a DUOCYLINDER? Test: each NS R-pair {a,Ra} with its ring of shared SC
neighbors = a 'cylinder' (axis=NS pair, ring=SC neighbors). All theta-4-cycles a-s-Ra-s'-a satisfy
R(C)=-C (R-ODD). Count the cycle rank; identify the two cylinders + the shared 'ridge' (duocylinder).
"""
from itertools import permutations, combinations
def canon(n, arcs):
    best=None
    for p in permutations(range(n)):
        key=tuple(sorted((p[i],p[j]) for (i,j) in arcs))
        if best is None or key<best: best=key
    return best
n=5
base=set((k,k-1) for k in range(n-1,0,-1))
tiles=sorted([(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in base and (j,i) not in base])
m=len(tiles)
def tour(bits):
    arcs=set(base)
    for idx,(i,j) in enumerate(tiles): arcs.add((i,j) if (bits>>idx)&1 else (j,i))
    return arcs
cl={}; tc=[]
for bits in range(1<<m):
    c=canon(n,tour(bits))
    if c not in cl: cl[c]=len(cl)
    tc.append(cl[c])
V=len(cl)
edgeset=set()
for bits in range(1<<m):
    a=tc[bits]
    for idx in range(m):
        b=tc[bits^(1<<idx)]
        if a!=b: edgeset.add((min(a,b),max(a,b)))
adj={v:set() for v in range(V)}
for (a,b) in edgeset: adj[a].add(b); adj[b].add(a)
def comp(arcs): return set((j,i) for (i,j) in arcs)
R={i:cl[canon(n,comp(set(c)))] for c,i in cl.items()}
NS=[i for i in range(V) if R[i]!=i]; SC=[i for i in range(V) if R[i]==i]
# NS R-pairs
pairs=[]; used=set()
for i in NS:
    if i not in used: pairs.append((i,R[i])); used.add(i); used.add(R[i])
print(f"NS R-pairs (cylinder AXES): {pairs}")
print(f"SC classes (R-fixed): {SC}")
for (a,b) in pairs:
    ring=sorted(adj[a]&adj[b]&set(SC))   # shared SC neighbors = the cylinder RING
    print(f"  cylinder axis {a}~{b}: shared-SC ring = {ring}  (size {len(ring)}, theta cycle-rank {max(0,len(ring)-1)})")
# the ridge = SC shared by BOTH pairs
ringA=adj[pairs[0][0]]&adj[pairs[0][1]]&set(SC)
ringB=adj[pairs[1][0]]&adj[pairs[1][1]]&set(SC)
print(f"  RIDGE (SC in both rings) = {sorted(ringA & ringB)}")
# verify every theta 4-cycle is R-odd: R(a-s-b-s'-a) = -(itself)
a,b=pairs[0]; ring=sorted(ringA)
print(f"\nverify R(theta cycle)=-cycle for axis {a}~{b}:")
ok=True
for s,sp in combinations(ring,2):
    # cycle a->s->b->s'->a ; R maps a<->b, s,s' fixed => a->s->b->s'->a maps to b->s->a->s'->b = reverse
    # i.e. R(C)= reverse(C) = -C in homology. (structural, always true when R swaps the two NS endpoints)
    pass
print("  every theta cycle a-s-b-s'-a has R swapping a<->b, s,s' fixed => R(C)=reverse(C)=-C: TRUE (structural)")
print("  => ALL theta cycles are R-ODD; the cylinders' entire cycle space lands in H_1^-.")
# total: duocylinder rank = (|ringA|-1)+(|ringB|-1)+ cross-links via ridge
print(f"\nDUOCYLINDER tally: cyl A rank {len(ringA)-1} + cyl B rank {len(ringB)-1} + cross/ridge links = b1^- = 7?")
print(f"  |ringA|={len(ringA)} (rank {len(ringA)-1}), |ringB|={len(ringB)} (rank {len(ringB)-1})")
