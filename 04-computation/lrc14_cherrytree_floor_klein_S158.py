#!/usr/bin/env python3
"""
klein-2026-07-07-S158 (part 2) -- the CHERRY-TREE floor at the k=8 criticality.

Bukszar-Prekopa t-cherry-tree bound (classical, provable by leaf induction like Hunter):
for events A_1..A_n and a t-cherry tree structure (start with an edge; each new vertex
attaches to an existing EDGE, contributing that edge's two pairs... precisely: the cherry
tree has 2n-3 pairs and n-2 triples):
   P(U A_i) <= S1 - Sum_{pairs in cherry} P(A_i n A_j) + Sum_{cherry triples} P(A_i n A_j n A_k).
At the k=8 criticality (S1 = 1 exactly, 7 events = the endpoint differences):
   W_end >= Sum_{11 cherry pairs} m2 - Sum_{5 cherry triples} m3.
By THM-638 each m2 >= theta^2 = 1/49, so the floor is >= 11/49 - Sum(5 chosen m3) --
and WE CHOOSE the cherry structure: pick triples with provably small mass.

THIS SCRIPT: per shape, greedily optimize the cherry tree (maximize Sum m2 - Sum m3,
numeric masses); report the cherry floor vs Hunter (6/49), vs the R-route bar 0.197,
vs true W; adversarial minimization over 8-sets (jump moves); also the AP.
"""
import numpy as np
from itertools import combinations

TH = 1/7

def masses(D, x):
    hits = [(((d*x) % 1.0) > 0) & (((d*x) % 1.0) <= TH) for d in D]
    H = np.stack(hits)
    n = len(D)
    m2 = {}; m3 = {}
    for a in range(n):
        for b in range(a+1, n):
            m2[(a,b)] = (H[a] & H[b]).mean()
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                m3[(a,b,c)] = (H[a] & H[b] & H[c]).mean()
    W = (~H.any(axis=0)).mean()
    return m2, m3, W

def cherry_floor(D, x):
    """greedy cherry-tree maximizing sum(m2 pairs) - sum(m3 triples)."""
    n = len(D)
    m2, m3, W = masses(D, x)
    best = None
    # try each starting edge (greedy attach afterwards)
    for e0 in combinations(range(n), 2):
        used = set(e0)
        edges = {tuple(sorted(e0))}
        val = m2[tuple(sorted(e0))]
        triples = []
        while len(used) < n:
            bestadd = None
            for v in range(n):
                if v in used: continue
                for (a, b) in list(edges):
                    tri = tuple(sorted((a, b, v)))
                    add = m2[tuple(sorted((a,v)))] + m2[tuple(sorted((b,v)))] - m3[tri]
                    if bestadd is None or add > bestadd[0]:
                        bestadd = (add, v, (a,b), tri)
            add, v, (a,b), tri = bestadd
            val += add
            used.add(v)
            edges.add(tuple(sorted((a,v)))); edges.add(tuple(sorted((b,v))))
            triples.append(tri)
        if best is None or val > best[0]:
            best = (val, triples)
    return best[0], W

def grid(NG): return (np.arange(NG)+0.5)/NG

if __name__ == "__main__":
    rng = np.random.default_rng(21580)
    x = grid(40009)
    bank = {
        "AP {1..8}": [1,2,3,4,5,6,7,8],
        "spread": [0,5,11,17,26,33,41,50],
        "geometric": [0,3,8,17,31,52,80,118],
        "two-cluster": [0,1,2,3,100,101,102,103],
        "mild": [0,2,5,9,14,20,27,35],
    }
    print("=== cherry-tree floor at k=8 (W_end >= sum 11 pairs - sum 5 chosen triples) ===")
    print(f"{'shape':>14} {'cherry floor':>13} {'Hunter 6/49':>12} {'true W':>8}   (R-route bar 0.197)")
    for nm, E in bank.items():
        et = max(E); D = sorted(et - e for e in E if e != et)
        cf, W = cherry_floor(D, x)
        print(f"{nm:>14} {cf:>13.4f} {6/49:>12.4f} {W:>8.4f}   [{'>=0.197 OK' if cf >= 0.197 else 'below'}]")
    # adversarial min of the cherry floor
    xs = grid(8009)
    gmin = (2.0, None)
    for trial in range(24):
        H = int(rng.choice([9, 14, 22, 40]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        def cf_of(EE, xg):
            et = max(EE); D = sorted(et - e for e in EE if e != et)
            return cherry_floor(D, xg)[0]
        cur = cf_of(E, xs)
        for step in range(30):
            i = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 26, 44]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 8: continue
            c = cf_of(cand, xs)
            if c < cur - 1e-4: E, cur = cand, c
        v = cf_of(E, x)
        if v < gmin[0]: gmin = (v, tuple(E))
    print(f"\n  ADVERSARIAL MIN cherry floor = {gmin[0]:.4f} at {gmin[1]}")
    print(f"  (if >= 0.197 uniformly: the k=8 R-route needs only R >= 0.75 + the 5-triple bounds;")
    print(f"   if >= 0 at the AP: cherry replaces Hunter everywhere with a better constant)")
