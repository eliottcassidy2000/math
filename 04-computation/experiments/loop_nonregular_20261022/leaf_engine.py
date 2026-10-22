"""Loop engine with leaves: the linear-forest cluster expansion (THM-4757) for NON-regular tournaments.
Session opus-2026-10-22-S54.

For regular T every cluster with a singleton run end vanishes (Shat 1 = 0). For general T such a block is a *leaf*:
a quotient vertex of degree 1, Moebius weight mu(1) = 1, whose free vertex sum turns the propagator of its path into a
vector ending in delta-type data. The core decomposition of loop_engine.py (pass-through pairings at core vertices,
m1/m2 path types, free-run orientations, closed runs excluded) applies verbatim when leaves are allowed as core vertices.
A core of excess e with k leaves is a connected multigraph with exactly k vertices of degree 1 (no loops there) and
all other vertices of degree >= 3, E - V = e. Counting degrees gives at most k + 2e vertices of degree >= 3. The pure
cycle (e = 0, k = 0) is the one-loop sum X_0 of THM-4757 (3); k = 0 and e >= 1 are the regular cores.

X(e, k; z) = sum over such cores of (1/|G|) sum over decompositions of mu * Hom(propagators at z Shat).

Usage: python leaf_engine.py   (self-test against cluster_brute_classes on random almost-regular tournaments)"""
import sys, os, itertools
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'sierpinski_ocf_20261022'))
import loop_engine_fast as LEF


def cores_leaves(excess, k):
    """connected multigraphs, loops allowed, vertices 0..k-1 of degree 1 (no loops), the others of degree >= 3,
    E - V = excess; up to isomorphism preserving the leaf set."""
    found = []
    for V3 in range(0, k + 2 * excess + 1):
        V = V3 + k
        E = V + excess
        if V == 0 or E < 1:
            continue
        if excess == 0 and k == 0:
            continue                                      # the pure cycle is X_0
        pairs = [(a, b) for a in range(V) for b in range(a, V) if not (a == b and a < k)]
        seen = set()
        perms = [lp + tuple(k + x for x in ip) for lp in itertools.permutations(range(k))
                 for ip in itertools.permutations(range(V3))]
        for ms in itertools.combinations_with_replacement(range(len(pairs)), E):
            edges = [pairs[i] for i in ms]
            deg = [0] * V
            for a, b in edges:
                deg[a] += 1; deg[b] += 1
            if any(deg[v] != 1 for v in range(k)) or any(deg[v] < 3 for v in range(k, V)):
                continue
            par = list(range(V))

            def f(x):
                while par[x] != x:
                    par[x] = par[par[x]]; x = par[x]
                return x
            for a, b in edges:
                par[f(a)] = f(b)
            if len({f(x) for x in range(V)}) != 1:
                continue
            canon = min(tuple(sorted((min(p[a], p[b]), max(p[a], p[b])) for a, b in edges)) for p in perms)
            if canon in seen:
                continue
            seen.add(canon); found.append((V, list(canon)))
    return found


_CACHE = {}


def get(excess, k):
    if (excess, k) not in _CACHE:
        _CACHE[(excess, k)] = [LEF.Core(V, e) for V, e in cores_leaves(excess, k)]
    return _CACHE[(excess, k)]


def X(Sh, excess, k, z=1.0, P=None):
    """class sum X_(e,k) at z (complex z allowed)."""
    if P is None:
        P = LEF.stack(Sh.astype(complex), z)
    return complex(sum(c.value(P) for c in get(excess, k)))


def series(Sh, classes, Lmax, K=48, r=0.5):
    """[z^L] of each class, L <= Lmax, by Cauchy sums on |z| = r."""
    vals = {c: np.zeros(K, dtype=complex) for c in classes}
    for kk in range(K):
        z = r * np.exp(2j * np.pi * kk / K)
        P = LEF.stack(Sh.astype(complex), z)
        for c in classes:
            vals[c][kk] = X(Sh, c[0], c[1], z, P)
    return {c: [(np.sum(v * np.exp(-2j * np.pi * np.arange(K) * L / K)) / K / r ** L).real for L in range(Lmax + 1)]
            for c, v in vals.items()}


if __name__ == '__main__':
    from cluster_brute_classes import classes as brute, almost_regular
    for (e, k) in [(-1, 2), (-1, 3), (-1, 4), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (2, 0)]:
        print(f'cores with excess {e}, {k} leaves: {len(cores_leaves(e, k))}', flush=True)
    for n, seed in ((8, 1), (10, 2)):
        rng = np.random.default_rng(seed)
        A = almost_regular(n, rng)
        S = (A - A.T).astype(float); Sh = S / n
        B = brute(Sh, 6)
        cls = [(-1, 2), (-1, 3), (-1, 4), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (2, 0)]
        ser = series(Sh, cls, 6)
        worst = 0.0
        for L in (2, 4, 6):
            for c in cls:
                b = B[L].get(c, 0.0); s = ser[c][L]
                worst = max(worst, abs(b - s))
                print(f'n={n} L={L} class {c}: brute {b:+.12e} engine {s:+.12e} diff {b - s:+.1e}', flush=True)
        print(f'n={n}: max |brute - engine| = {worst:.2e}', flush=True)
