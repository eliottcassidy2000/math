#!/usr/bin/env python3
"""Audit E (2026-10-14): rigidity claims of the results note (sections 2 and 5).
Generic rank of the rigidity matrix with random integer coordinates, computed mod a large prime by own Gaussian
elimination (rank mod p <= rank over Q = generic rank; equality with high probability, repeated with 3 seeds).
A graph is a (generic) rigidity circuit in R^d iff rank(E) = |E| - 1 and rank(E - e) = |E| - 1 for every edge e.
Checks: K_11 in R^10 (minimally rigid?), K_12 in R^10 (circuit?), K_6 - P_6 in R^2 (circuit?), wheel W_4 in R^2."""
import random
from itertools import combinations

PR = (1 << 61) - 1

def rank_mod(rows, p=PR):
    M = [r[:] for r in rows]
    rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], p - 2, p)
        M[rk] = [x * inv % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk

def rig_rows(n, d, edges, X):
    rows = []
    for (u, v) in edges:
        r = [0] * (n * d)
        for k in range(d):
            r[u * d + k] = X[u][k] - X[v][k]
            r[v * d + k] = X[v][k] - X[u][k]
        rows.append(r)
    return rows

def analyse(name, n, d, edges, seeds=(1, 2, 3)):
    res = []
    for s in seeds:
        rnd = random.Random(s)
        X = [[rnd.randrange(1, 10 ** 9) for _ in range(d)] for _ in range(n)]
        rk = rank_mod(rig_rows(n, d, edges, X))
        minus = [rank_mod(rig_rows(n, d, [f for f in edges if f != e], X)) for e in edges]
        res.append((rk, min(minus), max(minus)))
    E = len(edges); target = d * n - d * (d + 1) // 2
    rk, mn, mx = res[0]
    assert all(r == res[0] for r in res)
    circuit = (rk == E - 1 and mn == mx == E - 1)
    print(f"{name}: n={n}, d={d}, |E|={E}, dn - C(d+1,2) = {target}; generic rank {rk}; rank after deleting any one edge "
          f"in [{mn},{mx}]; rigid={rk == target}; independent={rk == E}; minimally rigid={rk == E == target}; circuit={circuit}")

K = lambda n: list(combinations(range(n), 2))
analyse("K_11 in R^10", 11, 10, K(11))
analyse("K_12 in R^10", 12, 10, K(12))
P6 = [(i, i + 1) for i in range(5)]
K6P6 = [e for e in K(6) if e not in P6]
analyse("K_6 - P_6 in R^2", 6, 2, K6P6)
W4 = [(0, i) for i in range(1, 5)] + [(1, 2), (2, 3), (3, 4), (4, 1)]
analyse("W_4 in R^2", 5, 2, W4)
analyse("K_4 in R^2", 4, 2, K(4))

def alpha(n, edges):
    E = set(edges) | {(v, u) for u, v in edges}
    return max(len(S) for k in range(n + 1) for S in combinations(range(n), k) if all((u, v) not in E for u, v in combinations(S, 2)))
print("independence number of P_6:", alpha(6, P6), "; independence number of the complement K_6 - P_6:", alpha(6, K6P6),
      "; clique number of K_6 - P_6 = alpha(P_6) =", alpha(6, P6))
