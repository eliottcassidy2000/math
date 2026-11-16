#!/usr/bin/env python3
"""fifteen.py: the typed 'fifteen / twenty-five' checks.  (1) RT_5: its 15 Hamiltonian paths, |Aut|, orbits and step
patterns; the 15 double transpositions of S_5 and their orbits under the same cyclic group.  (2) The square-sum graph on
1..25 as six antidiagonal matchings (i + j = 4, 9, 16, 25, 36, 49), by residue class mod 5, and the unique Hamiltonian path
of 1..15 cut into chunks of five.  (3) The prime square-sum graph on the 25 primes below 100: degrees, components, isolated
vertices, rows of five consecutive primes, the supersingular/non-supersingular split, and its Hamiltonian paths."""
import itertools, math
from collections import Counter
# (1)
S = [[1 if (j - i) % 5 in (1, 2) else 0 for j in range(5)] for i in range(5)]
paths = [p for p in itertools.permutations(range(5)) if all(S[p[i]][p[i+1]] for i in range(4))]
auts = [g for g in itertools.permutations(range(5)) if all(S[g[i]][g[j]] == S[i][j] for i in range(5) for j in range(5))]
orb = Counter(min(tuple(g[v] for v in p) for g in auts) for p in paths)
pat = Counter(tuple((p[i+1] - p[i]) % 5 for i in range(4)) for p in paths)
print(f"(1) RT_5 = circulant {{1,2}}: Hamiltonian paths {len(paths)}, |Aut| = {len(auts)}, orbit sizes {sorted(orb.values())}, step patterns {dict(pat)}")
def cyc_type(g):
    seen = set(); t = []
    for i in range(5):
        if i in seen: continue
        l = 0; j = i
        while j not in seen: seen.add(j); j = g[j]; l += 1
        t.append(l)
    return sorted(t)
dts = [g for g in itertools.permutations(range(5)) if cyc_type(g) == [1, 2, 2]]
rot = [tuple((i + a) % 5 for i in range(5)) for a in range(5)]
def conj(r, g):
    rinv = [0]*5
    for i in range(5): rinv[r[i]] = i
    return tuple(r[g[rinv[i]]] for i in range(5))
orb2 = Counter(min(conj(r, g) for r in rot) for g in dts)
print(f"    double transpositions of S_5: {len(dts)}, orbits under the rotation group Z_5: sizes {sorted(orb2.values())}; involutions of A_5: {sum(1 for g in itertools.permutations(range(5)) if cyc_type(g) == [1,2,2])}")
# (2)
n = 25; sq = {k*k for k in range(2, 8)}
adj = {i: [j for j in range(1, n+1) if j != i and i + j in sq] for i in range(1, n+1)}
print("(2) square-sum graph on 1..25 as antidiagonal matchings i + j = s:", {s: [(i, s - i) for i in range(1, n+1) if i < s - i <= n] for s in sorted(sq)})
print("    matching sizes:", {s: len([1 for i in range(1, n+1) if i < s - i <= n]) for s in sorted(sq)}, "total edges", sum(len(v) for v in adj.values()) // 2)
cls = {r: [i for i in range(1, n+1) if i % 5 == r] for r in range(5)}
inner = {r: [(i, j) for i in cls[r] for j in cls[r] if i < j and j in adj[i]] for r in range(5)}
between = [[sum(1 for i in cls[a] for j in cls[b] if j in adj[i]) for b in range(5)] for a in range(5)]
print("    residue classes mod 5: internal edges", {r: inner[r] for r in range(5)}, "; class-pair edge matrix", between)
def hampaths(m):
    A = {i: [j for j in range(1, m+1) if j != i and i + j in sq] for i in range(1, m+1)}
    start = min(A, key=lambda i: (len(A[i]), i)); res = []
    def dfs(v, vis, path):
        if len(path) == m: res.append(list(path)); return
        for w in A[v]:
            if w not in vis: vis.add(w); path.append(w); dfs(w, vis, path); path.pop(); vis.remove(w)
    dfs(start, {start}, [start]); return res
p15 = hampaths(15)[0]
print("    the unique Hamiltonian path of 1..15:", p15, "; chunks of five:", [p15[i:i+5] for i in range(0, 15, 5)], "chunk sums", [sum(p15[i:i+5]) for i in range(0, 15, 5)])
# (3)
primes = [p for p in range(2, 100) if all(p % q for q in range(2, int(p**0.5) + 1))]
P = {p: [q for q in primes if q != p and int(math.isqrt(p + q))**2 == p + q] for p in primes}
print(f"(3) the 25 primes below 100: {primes}")
print("    prime square-sum degrees:", [len(P[p]) for p in primes], "; isolated:", [p for p in primes if not P[p]])
comps = []; seen = set()
for p in primes:
    if p in seen: continue
    c = {p}; f = [p]
    while f:
        a = f.pop()
        for b in P[a]:
            if b not in c: c.add(b); f.append(b)
    seen |= c; comps.append(sorted(c))
print("    components:", comps)
rows = [primes[5*r:5*r+5] for r in range(5)]
print("    rows of five consecutive primes with an internal edge:", [r+1 for r, R in enumerate(rows) if any(q in P[p] for p in R for q in R)], "; row-pair edge matrix:", [[sum(1 for p in rows[a] for q in rows[b] if q in P[p]) for b in range(5)] for a in range(5)])
ss = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]; non = [p for p in primes if p not in ss]
print(f"    supersingular (Ogg) primes among them: {len(ss)}; non-supersingular: {non}; edges inside the supersingular 15: {sum(1 for p in ss for q in ss if p < q and q in P[p])}, inside the other 10: {sum(1 for p in non for q in non if p < q and q in P[p])}, between: {sum(1 for p in ss for q in non if q in P[p])}")
print("    mod-4 bipartition (odd primes): edges with p ≡ q mod 4:", sum(1 for p in primes for q in P[p] if p < q and p % 2 and q % 2 and p % 4 == q % 4))
res = []
def dfs(v, vis, path):
    if len(path) == 25: res.append(list(path)); return
    for w in P[v]:
        if w not in vis: vis.add(w); path.append(w); dfs(w, vis, path); path.pop(); vis.remove(w)
for s in primes: dfs(s, {s}, [s])
print("    Hamiltonian paths of the prime square-sum graph on the 25 primes (directed count):", len(res))
print("    Fermat primes 3, 5, 17, 257, 65537: supersingular?", {p: p in ss for p in (3, 5, 17, 257, 65537)})
