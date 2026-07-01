#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S81 -- kappa(n) = 1 + C(n-2,2): the optimal-shape FAMILY + the n=7 prediction test.

kappa(n) = min free arcs whose subcube surjects onto all A000568(n) iso classes.
Found kappa = 1,2,4,7 (n=3..6) = 1 + C(n-2,2) = m(n) - (n-3).  Predict kappa(7) = 1+C(5,2) = 11.

Part 1: enumerate ALL optimal free-arc configurations for n=4,5 (up to the free-arc-graph iso type) to
        characterize the SHAPE and the fixing ("configuration rule").
Part 2: TEST n=7 -- can 11 free arcs cover all 456 classes? (verify a construction by canon-ing only the
        2^11 subcube tournaments; and check that a natural 10-arc attempt fails to cover).
"""
import itertools
from math import comb, log2, ceil
from collections import defaultdict, Counter

A000568 = {3:2, 4:4, 5:12, 6:56, 7:456}

def perms_relabel(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    def relabel(mask, perm):
        out = 0
        for (i, j), k in idx.items():
            bit = (mask >> k) & 1
            pi, pj = perm[i], perm[j]
            a, b = (pi, pj) if pi < pj else (pj, pi)
            nb = bit if pi < pj else 1 - bit
            out |= nb << idx[(a, b)]
        return out
    return pairs, idx, len(pairs), perms, relabel

def canon_table(n):
    pairs, idx, P, perms, relabel = perms_relabel(n)
    canon_of = [0]*(1 << P); cid = {}; nid = 0
    for mask in range(1 << P):
        c = min(relabel(mask, p) for p in perms)
        if c not in cid: cid[c] = nid; nid += 1
        canon_of[mask] = cid[c]
    return pairs, P, canon_of, nid

def graph_type(F, pairs):
    """iso-type signature of the free-arc GRAPH: sorted degree sequence + #edges + #triangles."""
    arcs = [pairs[b] for b in F]
    deg = defaultdict(int); adj = defaultdict(set)
    for (i, j) in arcs: deg[i]+=1; deg[j]+=1; adj[i].add(j); adj[j].add(i)
    tri = 0
    vs = list(adj)
    for a, b, c in itertools.combinations(vs, 3):
        if b in adj[a] and c in adj[a] and c in adj[b]: tri += 1
    return (tuple(sorted(deg.values())), len(arcs), tri)

print("kappa(n) = 1 + C(n-2,2)  vs  info-floor ceil(log2 A000568)  vs  naive m=C(n-1,2):")
for n in range(3, 9):
    kap = 1 + comb(n-2, 2); m = comb(n-1, 2)
    fl = ceil(log2(A000568[n])) if n in A000568 else None
    print(f"  n={n}: kappa={kap:3d}  info-floor={fl}  naive m={m:3d}  (m-kappa = n-3 = {n-3};  kappa=m(n-1)+1)")

print("\n--- Part 1: optimal free-arc SHAPES (n=4,5), by free-graph iso-type ---")
for n in [4, 5]:
    pairs, P, canon_of, nclasses = canon_table(n)
    assert nclasses == A000568[n]
    k = 1 + comb(n-2, 2)
    shape_count = Counter(); example = {}
    for F in itertools.combinations(range(P), k):
        Fbits = list(F)
        subsets = []
        for s in range(1 << k):
            mm = 0
            for t in range(k):
                if (s >> t) & 1: mm |= (1 << Fbits[t])
            subsets.append(mm)
        rest = [b for b in range(P) if b not in F]
        ok = False
        for fixs in range(1 << len(rest)):
            fixmask = 0
            for t, b in enumerate(rest):
                if (fixs >> t) & 1: fixmask |= (1 << b)
            seen = set()
            for mm in subsets:
                seen.add(canon_of[mm | fixmask])
            if len(seen) == nclasses: ok = True; break
        if ok:
            gt = graph_type(F, pairs); shape_count[gt] += 1
            if gt not in example: example[gt] = [pairs[b] for b in F]
    print(f"  n={n} (kappa={k}): optimal free-arc graph-types (degseq, #edges, #triangles) -> count:")
    for gt, ct in shape_count.most_common():
        print(f"     {gt}  x{ct}   e.g. {example[gt]}")

print("\n--- Part 2: n=7 test (predict kappa=11). Verify constructions by canon-ing the subcube only ---")
pairs7, idx7, P7, perms7, relabel7 = perms_relabel(7)
def coverage(Fbits, fixmask):
    seen = set()
    for s in range(1 << len(Fbits)):
        mm = fixmask
        for t in range(len(Fbits)):
            if (s >> t) & 1: mm |= (1 << Fbits[t])
        c = min(relabel7(mm, p) for p in perms7)
        seen.add(c)
    return len(seen)
def arcs_to_bits(arcs): return [idx7[a] for a in arcs]
# Construction A (triangle-packing analog): triangles {0,1,2},{3,4,5} + edges to vertex 6 + bridges -> 11 arcs
tri1 = [(0,1),(0,2),(1,2)]; tri2 = [(3,4),(3,5),(4,5)]
constrA = tri1 + tri2 + [(0,3),(0,6),(3,6),(1,6),(2,6)]   # 3+3+5 = 11 arcs
constrA = constrA[:11]
Fbits = arcs_to_bits(constrA)
# fix all other arcs transitively (i->j for i<j): fixmask has bit=1 for every arc not free (i<j => i->j =1)
allbits = set(range(P7)); free=set(Fbits)
fixmask = 0
for b in range(P7):
    if b not in free: fixmask |= (1 << b)   # i->j for i<j (transitive backbone)
cov = coverage(Fbits, fixmask)
print(f"  Construction A (11 free arcs {constrA}):")
print(f"     covers {cov}/{A000568[7]} classes with transitive fixing  {'ALL' if cov==456 else '(partial)'}")
