#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S81 -- THE MINIMUM-FREE-ARCS (transversal subcube) invariant kappa(n).

Reinterpretation of the owner's n=4 observation ("all 4 iso classes via 3 tile-flips naively, but 2 arcs
suffice given a configuration rule on the fixed arcs"): NOT the Hamiltonian-path tiling (fix n-1 path arcs,
flip m=C(n-1,2) tiles), but the MINIMUM number k of FREE arcs -- choose k of the C(n,2) arcs to vary and FIX
the other C(n,2)-k to specific orientations -- such that the 2^k resulting tournaments realize EVERY one of
the A000568(n) isomorphism classes.

    kappa(n) := min { k : exists a set F of k arcs and a fixing of the other arcs
                          such that {tournaments in the subcube} hits all iso classes }.

Info-theoretic floor: kappa(n) >= ceil(log2 A000568(n)).  Naive tiling uses k = C(n-1,2).
We find kappa(n), the optimal free-arc set F (its SHAPE), and the fixing (the "configuration rule").
"""
import itertools
from math import comb, log2, ceil
from collections import defaultdict

A000568 = {3:2, 4:4, 5:12, 6:56, 7:456}

def build(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    P = len(pairs)
    # canon-id for every tournament bitmask
    perms = list(itertools.permutations(range(n)))
    def relabel(mask, perm):
        out = 0
        for (i, j), k in idx.items():
            bit = (mask >> k) & 1                # 1 => i->j
            pi, pj = perm[i], perm[j]
            a, b = (pi, pj) if pi < pj else (pj, pi)
            nb = bit if pi < pj else 1 - bit     # orientation after relabel
            out |= nb << idx[(a, b)]
        return out
    canon_of = [0]*(1 << P); cid = {}; nid = 0
    for mask in range(1 << P):
        c = min(relabel(mask, p) for p in perms)
        if c not in cid: cid[c] = nid; nid += 1
        canon_of[mask] = cid[c]
    return pairs, idx, P, canon_of, nid

def min_free_arcs(n, kmax=None):
    pairs, idx, P, canon_of, nclasses = build(n)
    assert nclasses == A000568[n], (n, nclasses)
    lo = ceil(log2(nclasses))
    kmax = kmax or P
    for k in range(lo, kmax + 1):
        for F in itertools.combinations(range(P), k):
            Fmask = 0
            for b in F: Fmask |= (1 << b)
            restbits = [b for b in range(P) if b not in F]
            # subsets of F (the 2^k varying assignments)
            Fbits = list(F)
            subsets = []
            for s in range(1 << k):
                mm = 0
                for t in range(k):
                    if (s >> t) & 1: mm |= (1 << Fbits[t])
                subsets.append(mm)
            # try each fixing of the rest
            for fixs in range(1 << len(restbits)):
                fixmask = 0
                for t, b in enumerate(restbits):
                    if (fixs >> t) & 1: fixmask |= (1 << b)
                seen = set()
                for mm in subsets:
                    seen.add(canon_of[mm | fixmask])
                    if len(seen) == nclasses: break
                if len(seen) == nclasses:
                    return k, F, fixmask, pairs
    return None

def shape_desc(F, pairs, n):
    arcs = [pairs[b] for b in F]
    # describe the free-arc graph: degrees, is-matching, spans
    deg = defaultdict(int)
    for (i, j) in arcs: deg[i] += 1; deg[j] += 1
    verts = sorted(deg)
    matching = all(v <= 1 for v in deg.values())
    return arcs, dict(deg), matching, verts

print(f"{'n':>2} {'A000568':>8} {'log2':>6} {'ceil':>5} {'kappa':>6} {'naive m':>8}  free-arc shape")
for n in [3, 4, 5, 6]:
    res = min_free_arcs(n)
    if res is None:
        print(f"{n:>2}  (search exceeded)"); continue
    k, F, fixmask, pairs = res
    arcs, deg, matching, verts = shape_desc(F, pairs, n)
    m = comb(n-1, 2)
    print(f"{n:>2} {A000568[n]:>8} {log2(A000568[n]):>6.2f} {ceil(log2(A000568[n])):>5} {k:>6} {m:>8}  "
          f"free arcs={arcs} deg={deg} matching={matching}")
