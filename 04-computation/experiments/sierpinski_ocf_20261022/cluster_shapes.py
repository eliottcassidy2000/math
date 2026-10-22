# Excess-1 connected clusters at L = 4, split by core shape, for Paley P_q (exact coefficient extraction).
# figure-eight (2,2): Z = +-D, D = sum_v (Shat^2)_vv^2; theta (2,1,1): Z = +-Th, Th = sum_(u,v) (Shat^2)_uv Shat_uv^2.
# Reports the Moebius-weighted coefficient of each, as a rational (denominators divide 4!).
import sys
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from cluster_check import compositions, set_partitions, mobius, tournament
from math import factorial

def shape_sums(n, kind='paley', L=4):
    A = tournament(kind, n); S = A - A.T; Sh = S / n
    letters = 'abcdefghijklmnop'
    res = {}
    for comp in compositions(L):
        p = len(comp)
        pos = [(j, i) for j in range(p) for i in range(comp[j] + 1)]
        ends = {(j, 0) for j in range(p)} | {(j, comp[j]) for j in range(p)}
        for blocks in set_partitions(pos):
            bid = {q: t for t, b in enumerate(blocks) for q in b}
            if any(len(blocks[bid[e]]) == 1 for e in ends): continue
            edges = [(bid[(j, i)], bid[(j, i + 1)]) for j in range(p) for i in range(comp[j])]
            if any(a == b for a, b in edges): continue
            nb = len(blocks)
            if len(edges) - nb != 1: continue
            par = list(range(nb))
            def f(a):
                while par[a] != a: par[a] = par[par[a]]; a = par[a]
                return a
            for a, b in edges: par[f(a)] = f(b)
            if len({f(a) for a in range(nb)}) != 1: continue
            deg = [0] * nb
            for a, b in edges: deg[a] += 1; deg[b] += 1
            shape = tuple(sorted(deg, reverse=True))
            spec = ','.join(letters[a] + letters[b] for a, b in edges) + '->'
            Z = np.einsum(spec, *([Sh] * len(edges)), optimize='greedy')
            w = Fraction(mobius(blocks), factorial(p))
            res.setdefault(shape, []).append((w, Z, comp, blocks))
    return res, Sh

for q in (43, 47):
    res, Sh = shape_sums(q)
    S2 = Sh @ Sh
    D = np.sum(np.diag(S2) ** 2)
    Th = np.sum(S2 * Sh * Sh)
    for shape, lst in res.items():
        ref = D if shape == (4, 2, 2) else Th
        coef = Fraction(0)
        for w, Z, comp, blocks in lst:
            s = Z / ref
            si = int(round(s))
            assert abs(s - si) < 1e-9, (shape, s)
            coef += w * si
        print(f"q={q} shape={shape}: {len(lst)} configurations, coefficient = {coef} = {float(coef):.6f}")
    print(f"   D = {D:.6e} (q D = {q * D:.6f}),  Th = {Th:.6e} (q^2 Th = {q * q * Th:.6f})")
