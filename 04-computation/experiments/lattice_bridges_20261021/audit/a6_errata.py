#!/usr/bin/env python3
"""AUDIT of MISTAKE-659 and the errata banner of 07-reflections/sphere-packing-eight-twentyfour.md.
 * P_7: directed 3-, 5-, 7-cycle counts, alpha_2 (vertex-disjoint pairs of odd cycles), OCF H = 1 + 2 a1 + 4 a2;
 * the 14 cyclic triples of P_7 = dev{0,1,3} u dev{0,1,5}, two disjoint Fano planes, one orientation each;
   automorphism group of the oriented Fano plane (extra check: the reflection calls it G_2);
 * regular 5-tournaments: all 24 labelled ones have H = 15 (HYP-1363's 13 is wrong);
 * Coxeter numbers h = |roots|/rank for E_6, E_8, G_2, A_6; dodecahedron/icosahedron counts;
 * Epstein zeta of E_8 with roots of norm 2: sum_v |v|^(-2s) = 2^(-s) 240 zeta(s) zeta(s-3) (numerical check at s = 6),
   and its continued value at s = 2.
Output: a6_errata.out
"""
import os, itertools, math
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a6_errata.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

n = 7; QR = {1, 2, 4}
adj = [[1 if (j - i) % 7 in QR else 0 for j in range(7)] for i in range(7)]
def directed_cycles(adj, k):
    """set of directed k-cycles as (frozenset of vertices, canonical rotation tuple)."""
    n = len(adj); cyc = set()
    def dfs(path):
        if len(path) == k:
            if adj[path[-1]][path[0]]:
                r = min(range(k), key=lambda i: path[i]); cyc.add(tuple(path[r:] + path[:r]))
            return
        for w in range(n):
            if w not in path and adj[path[-1]][w] and w > path[0]: dfs(path + [w])
    for v in range(n): dfs([v])
    return cyc
C = {k: directed_cycles(adj, k) for k in (3, 5, 7)}
odd = [c for k in (3, 5, 7) for c in C[k]]
a1 = len(odd)
a2 = sum(1 for x, y in itertools.combinations(odd, 2) if not set(x) & set(y))
def ham(adj):
    n = len(adj); f = {}
    for v in range(n): f[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = f.get((mask, v))
            if not c: continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]: f[(mask | 1 << w, w)] = f.get((mask | 1 << w, w), 0) + c
    return sum(f.get(((1 << n) - 1, v), 0) for v in range(n))
say(f"P_7: c3 = {len(C[3])}, c5 = {len(C[5])}, c7 = {len(C[7])}, alpha_1 = {a1}, alpha_2 = {a2}; 1 + 2a1 + 4a2 = {1 + 2*a1 + 4*a2}; H (DP) = {ham(adj)}")
tri = sorted(tuple(sorted(c)) for c in C[3])
dev013 = sorted(tuple(sorted(((x + t) % 7) for x in (0, 1, 3))) for t in range(7))
dev015 = sorted(tuple(sorted(((x + t) % 7) for x in (0, 1, 5))) for t in range(7))
def is_fano(blocks):
    return len(blocks) == 7 and all(sum(1 for b in blocks if a in b and c in b) == 1 for a, c in itertools.combinations(range(7), 2))
say(f"cyclic triples = dev{{0,1,3}} u dev{{0,1,5}}: {sorted(dev013 + dev015) == tri}; disjoint: {not set(dev013) & set(dev015)}; "
    f"both Fano planes: {is_fano(dev013)}, {is_fano(dev015)}; distinct underlying triples: {len(set(tri))} (so no triple carries both orientations)")
say(f"out-neighbourhood lines {{1,2,4}}+t form dev{{0,1,3}}: {sorted(tuple(sorted(((x + t) % 7) for x in (1, 2, 4))) for t in range(7)) == dev013}")
# automorphisms of the oriented Fano plane dev{0,1,3} with the P_7 orientations
orient = {c for c in C[3] if tuple(sorted(c)) in set(dev013)}
def rot_canon(t):
    r = min(range(3), key=lambda i: t[i]); return tuple(t[r:] + t[:r])
aut = [p for p in itertools.permutations(range(7)) if {rot_canon(tuple(p[x] for x in c)) for c in orient} == orient]
autu = [p for p in itertools.permutations(range(7)) if {tuple(sorted(p[x] for x in c)) for c in orient} == set(dev013)]
say(f"oriented Fano plane (P_7 orientation): |Aut| = {len(aut)}; unoriented Fano plane: |Aut| = {len(autu)} (= |GL(3,2)| = 168)")
# regular 5-tournaments
regs = []
pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for bits in range(1 << 10):
    A = [[0] * 5 for _ in range(5)]
    for t, (i, j) in enumerate(pairs):
        if (bits >> t) & 1: A[i][j] = 1
        else: A[j][i] = 1
    if all(sum(r) == 2 for r in A): regs.append(A)
say(f"labelled regular 5-tournaments: {len(regs)}; their H values: {sorted({ham(A) for A in regs})}  (HYP-1363 says 13)")
# Coxeter numbers and polyhedra
roots = {'E6': (72, 6), 'E8': (240, 8), 'G2': (12, 2), 'A6': (42, 6)}
say("Coxeter numbers h = |roots|/rank:", {k: v[0] // v[1] for k, v in roots.items()})
say("dodecahedron (V,E,F) = (20,30,12), Euler:", 20 - 30 + 12, "; icosahedron (12,30,20), Euler:", 12 - 30 + 20, "; 240 = 20*12 =",
    "(vertices of dodecahedron) x (vertices of icosahedron)")
# Epstein zeta of E_8 (roots of norm 2): the norm-2m shell has 240 sigma_3(m) vectors
mp.mp.dps = 30
s = 6
def sig3(m): return sum(d ** 3 for d in range(1, m + 1) if m % d == 0)
M = 20000
# sigma_3 via sieve
sig = [0] * (M + 1)
for d in range(1, M + 1):
    for k in range(d, M + 1, d): sig[k] += d ** 3
partial = mp.fsum(240 * sig[m] * mp.mpf(2 * m) ** (-s) for m in range(1, M + 1))
f1 = 240 * mp.mpf(2) ** (-s) * mp.zeta(s) * mp.zeta(s - 3)
f2 = 240 * mp.zeta(s) * mp.zeta(s - 3)
say(f"E_8 Epstein zeta at s = 6: partial sum (m <= {M}) = {mp.nstr(partial, 15)}; 240*2^-s*zeta(s)zeta(s-3) = {mp.nstr(f1, 15)}; "
    f"240*zeta(s)zeta(s-3) = {mp.nstr(f2, 15)}")
say(f"   continued value at s = 2 (roots of norm 2): 240*2^-2*zeta(2)*zeta(-1) = {mp.nstr(240 * mp.mpf(2) ** -2 * mp.zeta(2) * mp.zeta(-1), 12)} = -5 pi^2/6 = {mp.nstr(-5 * mp.pi ** 2 / 6, 12)}; "
    f"residue at the only pole s = 4: {mp.nstr(240 * mp.mpf(2) ** -4 * mp.zeta(4), 12)}; zeta(s)zeta(s-3) at s = 1 is finite (zeta(-2) = 0): "
    f"limit = {mp.nstr(mp.limit(lambda x: mp.zeta(x) * mp.zeta(x - 3), 1), 10)}")
out.close()
