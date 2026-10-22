# Tree level (excess -1) of the cluster expansion for NON-regular tournaments, by brute force over (run tuples,
# partitions) without the regular-only pruning, against X_(-1)(z) = -z^2 delta^T (I + z^2 Shat^2)^-1 delta, delta = Shat 1.
import sys, itertools
import numpy as np
from math import factorial
sys.path.insert(0, '.')
from cluster_check import compositions, set_partitions, mobius
rng = np.random.default_rng(3)
for n in (6, 7):
    A = np.triu((rng.random((n, n)) < 0.5).astype(float), 1); A = A + np.triu(1 - A, 1).T * (np.tri(n, n, -1) > 0)
    A = np.triu(A, 1) + (1 - np.triu(A, 1)).T * np.tri(n, n, -1)
    S = A - A.T; Sh = S / n; delta = Sh @ np.ones(n)
    letters = 'abcdefghijklmnop'
    for L in (2, 4, 6):
        tot = {}
        for comp in compositions(L):
            p = len(comp)
            pos = [(j, i) for j in range(p) for i in range(comp[j] + 1)]
            for blocks in set_partitions(pos):
                bid = {q: t for t, b in enumerate(blocks) for q in b}
                edges = [(bid[(j, i)], bid[(j, i + 1)]) for j in range(p) for i in range(comp[j])]
                if any(a == b for a, b in edges): continue
                nb = len(blocks)
                par = list(range(nb))
                def f(x):
                    while par[x] != x: par[x] = par[par[x]]; x = par[x]
                    return x
                for a, b in edges: par[f(a)] = f(b)
                if len({f(x) for x in range(nb)}) != 1: continue
                exc = len(edges) - nb
                if exc != -1: continue
                spec = ','.join(letters[a] + letters[b] for a, b in edges) + '->'
                Z = np.einsum(spec, *([Sh] * len(edges)), optimize='greedy')
                tot[exc] = tot.get(exc, 0.0) + mobius(blocks) * Z / factorial(p)
        m = L // 2
        pred = -delta @ np.linalg.matrix_power(-Sh @ Sh, m - 1) @ delta     # [z^L] of -z^2 d^T (I + z^2 S^2)^-1 d
        alt = -delta @ np.linalg.matrix_power(Sh @ Sh, m - 1) @ delta       # [z^L] of -z^2 d^T (I - z^2 S^2)^-1 d
        print(f"n={n} L={L}: brute excess -1 = {tot.get(-1, 0):+.10f}   -d(I+S^2)^-1 d coeff = {pred:+.10f}   -d(I-S^2)^-1 d coeff = {alt:+.10f}", flush=True)
