# The octahedron identity of HYP-9275 from PGL(2) invariance.  Injective sums T (distinct points) versus HYP-9275's
# all-maps sums U, at primes of both classes mod 4:
#   (1) T_oct = (q - 5) T_W4                                   (PGL(2): p + 1 points of P^1 minus the six vertices),
#   (2) U = sum over kernel partitions into independent blocks: U(W4) = T_W4 - q(q-1)(q-4),
#       U(oct) = T_oct + 3(q-4) T_C4 + 3(q)_4 + (q)_3 with T_C4 = -q(q-1)(q-3),
#   (3) hence HYP-9275's U(K_222) = q(q-1)((5-q)a_q + 3q) follows from HYP-9274(B)'s U(W4) = -q(q-1)(q + a_q).
# Octahedron labelling: antipodal pairs {1,2},{3,4},{5,6}, arcs i -> j for i < j (all in-degrees even: 0,0,2,2,4,4).
# W4 labelling as in paley_graph_sums.py: hub 0, non-adjacent pairs {1,2},{3,4}.  mac-mini-2026-10-25-pyth.
import numpy as np
from sympy import primerange
import sys
sys.argv = [sys.argv[0]]
from graph_sums6 import S as S6, GRAPHS
OCT = GRAPHS['K6-3K2']
W4 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)]

def chi_table(p):
    c = np.full(p, -1.0); c[0] = 0
    c[[(x * x) % p for x in range(1, p)]] = 1
    return c

def U_all_maps(p, E, nv):
    """sum over ALL maps x: [nv] -> F_p (coincidences allowed) of prod chi(x_j - x_i), by brute force over x_1 = 0
    (translation) -- O(p^(nv-1)); only for small p"""
    c = chi_table(p)
    import itertools
    tot = 0.0
    idx = np.arange(p)
    # vectorise the last two coordinates
    for rest in itertools.product(range(p), repeat=nv - 3):
        x = [0] + list(rest)
        val = np.ones((p, p))
        X1, X2 = np.meshgrid(idx, idx, indexing='ij')
        full = x + [X1, X2]
        for i, j in E:
            a, b = full[i], full[j]
            val = val * c[(np.asarray(b) - np.asarray(a)) % p]
        tot += val.sum()
    return int(round(tot)) * p

def T_inj(p, E, nv):
    c = chi_table(p)
    import itertools
    tot = 0
    for x in itertools.permutations(range(p), nv):
        if x[0] != 0: continue
        v = 1
        for i, j in E:
            v *= c[(x[j] - x[i]) % p]
            if v == 0: break
        tot += v
    return int(tot) * p

f8 = {5: -2, 7: 24, 11: -44, 13: 22}
for p in (5, 7, 11, 13):
    a = f8[p]
    # injective octahedron via the O(p^4) matrix routine (vertices 1..6) and injective W4 by brute force
    Soct = S6(p, OCT); To = p * (p - 1) * Soct
    Tw = T_inj(p, W4, 5)
    Uw = U_all_maps(p, W4, 5)
    Uo = U_all_maps(p, [(i - 1, j - 1) for i, j in OCT], 6) if p <= 11 else None
    print(f"p = {p} (mod 4 = {p % 4}): T_oct = {To}, (p-5) T_W4 = {(p - 5) * Tw}  [{'OK' if To == (p - 5) * Tw else 'FAIL'}];"
          f" T_W4 = q(q-1)(-4-a): {Tw == p * (p - 1) * (-4 - a)}; U(W4) = T_W4 - q(q-1)(q-4): {Uw == Tw - p * (p - 1) * (p - 4)};"
          f" U(W4) = -q(q-1)(q+a): {Uw == -p * (p - 1) * (p + a)}"
          + (f"; U(oct) = q(q-1)((5-q)a+3q): {Uo == p * (p - 1) * ((5 - p) * a + 3 * p)}; "
             f"U(oct) = T_oct + 3(q-4)T_C4 + 3(q)_4 + (q)_3: {Uo == To + 3 * (p - 4) * (-p * (p - 1) * (p - 3)) + 3 * p * (p - 1) * (p - 2) * (p - 3) + p * (p - 1) * (p - 2)}" if Uo is not None else ''))

# (1) at more primes of both classes, with the O(p^4) injective routine and the W4 table of paley_graph_sums.out
import re
SW4 = {}
for l in open('paley_graph_sums.out'):
    m = re.match(r'p=(\d+) mod4=\d K4=-?\d+ K5=-?\d+ K5-e=-?\d+ K5-2K2=(-?\d+)', l)
    if m: SW4[int(m.group(1))] = int(m.group(2))
ok = {1: 0, 3: 0}; tot = {1: 0, 3: 0}
for p in primerange(5, 110):
    tot[p % 4] += 1
    if S6(p, OCT) == (p - 5) * SW4[p]: ok[p % 4] += 1
print(f"S_oct = (q - 5) S_W4: q = 1 mod 4: {ok[1]}/{tot[1]} primes, q = 3 mod 4: {ok[3]}/{tot[3]} primes (5 <= q < 110)")
