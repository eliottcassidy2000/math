#!/usr/bin/env python3
"""audit_anf5.py -- independent check of the 'for the record' statement of section 6.3: algebraic normal form over F_2 of
c_odd mod 2 = (c_3 + c_5) mod 2 on a labelled 5-set, in the arc variables x_ij = [i -> j], i < j (Moebius transform of the
truth table; cycles counted directly)."""
import itertools
n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
E = len(pairs)
def codd(t):
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        if (t >> b) & 1: adj[i][j] = 1
        else: adj[j][i] = 1
    c = 0
    for k in (3, 5):
        for S in itertools.combinations(range(n), k):
            for p in itertools.permutations(S[1:]):
                seq = (S[0],) + p
                if all(adj[seq[i]][seq[(i + 1) % k]] for i in range(k)): c += 1
    return c & 1
f = [codd(t) for t in range(1 << E)]
for b in range(E):
    for t in range(1 << E):
        if (t >> b) & 1: f[t] ^= f[t ^ (1 << b)]
def kind(t):
    es = [pairs[b] for b in range(E) if (t >> b) & 1]
    deg = {}
    for i, j in es:
        deg[i] = deg.get(i, 0) + 1; deg[j] = deg.get(j, 0) + 1
    return (len(es), tuple(sorted(deg.values(), reverse=True)))
hist = {}
for t in range(1 << E):
    if f[t]:
        k = kind(t); hist[k] = hist.get(k, 0) + 1
names = {(0, ()): "constant", (1, (1, 1)): "single arc", (2, (2, 1, 1)): "cherry", (2, (1, 1, 1, 1)): "two disjoint arcs",
         (3, (2, 2, 1, 1)): "path on four vertices", (3, (2, 1, 1, 1, 1)): "cherry with a disjoint arc", (4, (2, 2, 2, 1, 1)): "Hamiltonian path of K_5 (or other 4-edge graph with this degree sequence)"}
print("monomials of the normal form of (c_3 + c_5) mod 2 on 5 labelled vertices:", sum(hist.values()))
for k, v in sorted(hist.items()):
    print(f"  degree {k[0]}, degree sequence {k[1]}: {v}   [{names.get(k, '?')}]")
