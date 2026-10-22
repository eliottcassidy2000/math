"""Brute-force cluster sums of the linear-forest loop expansion (THM-4757), classified by excess AND number of leaves,
for an arbitrary (non-regular) tournament.  Session opus-2026-10-22-S54.

log Phi_T(z) = sum over connected (run tuple, partition of run positions) of mu(pi) Z_pi z^L / p!   (THM-4757 (2)).
A leaf is a block of quotient-degree 1, i.e. a singleton run end (Moebius weight 1); for regular T every cluster with
a leaf vanishes because Shat 1 = 0.  Here we keep them and report [z^L] of the class sums X_(e,k)(z) for excess e and k
leaves.  Clusters containing a run of odd length cancel in pairs within each class (run reversal), so only
compositions of L into even parts are enumerated.
Usage: python cluster_brute_classes.py n LMAX [seed]   (a random almost-regular tournament on even n)
       or import classes(Sh, LMAX)."""
import sys, os, itertools
from math import factorial
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'sierpinski_ocf_20261022'))
from cluster_check import set_partitions, mobius

LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def even_compositions(L):
    if L == 0:
        yield (); return
    for first in range(2, L + 1, 2):
        for rest in even_compositions(L - first):
            yield (first,) + rest


def classes(Sh, LMAX):
    """{L: {(excess, leaves): [z^L] class sum}} for even L <= LMAX."""
    out = {}
    for L in range(2, LMAX + 1, 2):
        tot = defaultdict(float)
        cache = {}
        for comp in even_compositions(L):
            p = len(comp)
            pos = [(j, i) for j in range(p) for i in range(comp[j] + 1)]
            for blocks in set_partitions(pos):
                bid = {}
                for t, b in enumerate(blocks):
                    for q in b:
                        bid[q] = t
                edges = [(bid[(j, i)], bid[(j, i + 1)]) for j in range(p) for i in range(comp[j])]
                if any(a == b for a, b in edges):
                    continue                                   # a loop arc has Shat_vv = 0
                nb = len(blocks)
                par = list(range(nb))

                def f(x):
                    while par[x] != x:
                        par[x] = par[par[x]]; x = par[x]
                    return x
                for a, b in edges:
                    par[f(a)] = f(b)
                if len({f(x) for x in range(nb)}) != 1:
                    continue
                deg = [0] * nb
                for a, b in edges:
                    deg[a] += 1; deg[b] += 1
                k = deg.count(1)
                e = len(edges) - nb
                key = tuple(edges)
                if key not in cache:
                    spec = ','.join(LETTERS[a] + LETTERS[b] for a, b in edges) + '->'
                    cache[key] = float(np.einsum(spec, *([Sh] * len(edges)), optimize='greedy'))
                tot[(e, k)] += mobius(blocks) * cache[key] / factorial(p)
        out[L] = dict(tot)
    return out


def almost_regular(n, rng, steps=None):
    """a random almost-regular tournament on even n (scores n/2, n/2 - 1): delete a vertex of the rotational
    (n+1)-tournament, then random 3-cycle reversals (score preserving)."""
    m = n + 1; d = n // 2
    A = np.array([[1 if 0 < (j - i) % m <= d else 0 for j in range(m)] for i in range(m)])[:n, :n]
    for _ in range(steps or 50 * n * n):
        u, v, w = rng.choice(n, 3, replace=False)
        if A[u, v] and A[v, w] and A[w, u]:
            A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
    return A


if __name__ == '__main__':
    n = int(sys.argv[1]); LMAX = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = np.random.default_rng(seed)
    A = almost_regular(n, rng)
    assert sorted(A.sum(1)) == sorted([n // 2] * (n // 2) + [n // 2 - 1] * (n // 2))
    S = (A - A.T).astype(float); Sh = S / n
    res = classes(Sh, LMAX)
    for L, d in res.items():
        print(f'L={L}: ' + ', '.join(f'(e={e},k={k}) {v:+.12e}' for (e, k), v in sorted(d.items())), flush=True)
