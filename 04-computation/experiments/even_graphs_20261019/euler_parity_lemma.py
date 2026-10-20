#!/usr/bin/env python3
"""Euler parity lemma test: for odd-order P <= S_n and even-order h in N(P), G = <P,h>, is there a G-invariant EULER graph
X (all degrees even) with r(h,X) = 1 (h reverses an odd number of its edges)?  Also: for odd-order G, is every G-invariant
even-weight degree-parity vector realised by a G-invariant graph?  n <= NMAX, P from subgroup_reps."""
import sys, itertools
sys.path.insert(0, '.')
from subgroups_fix import cyc, mul, closure, pair_orbits, subgroup_reps
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
def order(g):
    idt = tuple(range(len(g))); x = g; k = 1
    while x != idt: x = mul(g, x); k += 1
    return k
def r(h, X):
    s = 0
    for (u, v) in X:
        if h[u] > h[v]: s ^= 1
    return s
def rank_f2(rows):
    rows = [row for row in rows]; rk = 0; ncol = max((len(row) for row in rows), default=0)
    pivots = []
    M = [int(''.join(map(str, row)), 2) if row else 0 for row in rows]
    basis = []
    for v in M:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); rk += 1
    return rk
for n in range(5, NMAX + 1):
    perms = list(itertools.permutations(range(n)))
    bad = 0; tested = 0; real_bad = 0; real_tested = 0
    for name, P in subgroup_reps(n).items():
        N = [s for s in perms if frozenset(mul(mul(s, g), inv(s)) for g in P) == P]
        for h in N:
            G = closure(list(P) + [h], n)
            orbs = pair_orbits(G, n)
            if order(h) % 2 == 0:
                # G-invariant Euler graphs: subsets of orbits with all degrees even -> solve over F_2: variables = orbits
                # degree parity of vertex v in orbit O = number of pairs in O containing v (mod 2)
                nO = len(orbs)
                deg = [[sum(1 for (a, b) in O if v in (a, b)) % 2 for O in orbs] for v in range(n)]
                rv = [r(h, O) for O in orbs]
                # Need x in F_2^nO with deg*x = 0 and rv.x = 1: exists iff rv not in rowspace(deg)
                rk1 = rank_f2(deg); rk2 = rank_f2(deg + [rv])
                tested += 1
                if rk2 == rk1: bad += 1; print(f"  n={n} FAIL: P={name}, h={h}, |G|={len(G)}: no h-odd G-invariant Euler graph")
            else:
                if h in P or len(G) % 2 == 0: continue
                # realisability: G-invariant even-weight parity vectors (constant on G-orbits of vertices)
                vorbs = []; seen = set()
                for v in range(n):
                    if v in seen: continue
                    o = sorted({g[v] for g in G}); seen |= set(o); vorbs.append(o)
                deg = [[sum(1 for (a, b) in O if v in (a, b)) % 2 for O in orbs] for v in range(n)]
                # image of degree map restricted to G-invariant graphs: span of columns (one per orbit)
                cols = [[deg[v][j] for v in range(n)] for j in range(len(orbs))]
                rk = rank_f2(cols)
                # target space: vectors constant on vertex orbits with even total weight
                cand = 0
                for bits in itertools.product([0, 1], repeat=len(vorbs)):
                    vec = [0] * n
                    for b, o in zip(bits, vorbs):
                        for v in o: vec[v] = b
                    if sum(vec) % 2: continue
                    real_tested += 1
                    if rank_f2(cols + [vec]) > rk: real_bad += 1
    print(f"n = {n}: Euler parity lemma: {tested - bad}/{tested} pairs (P,h) OK; realisability: {real_tested - real_bad}/{real_tested} OK")
