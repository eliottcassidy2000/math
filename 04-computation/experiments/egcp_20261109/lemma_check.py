"""Random brute-force tests of the new lemmas of THM-5040 (formulas and evenness criteria).

  Lemma 5  (E_2-multi-join)   d(E_2^t v G) = 2^n - 2^g - 2t + d(G); even iff G even.
  Lemma 7  (stable recursion) G even, leaf q, no universal vertex: d(SR(G)) = 3*2^(n-1) + d(G); SR(G) even, leaf,
                              no universal vertex.
  Lemma 8  (jump)             d(J(H,z)) = 2^(m-1) - D(H); if the hub is the unique vertex of maximum degree then
                              J(H,z) is even iff sgn_H(g) = sgn(g) for every automorphism g of H fixing z.
  Complement formula          d(G) = 2^n - D(comp G).
"""
import random
from egcp_lib import (mk, comp, degs, dcount, Dcount, automorphisms, edge_sign, perm_sign, is_even, has_leaf,
                      no_universal, E2join, SR, jump)

def rand_graph(n, p, rng):
    return mk(n, [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < p])

rng = random.Random(20261109)
res = {}
# complement formula
bad = 0
for _ in range(300):
    g = rand_graph(rng.randint(1, 9), rng.random(), rng)
    if dcount(g) != 2**g[0] - Dcount(comp(g)): bad += 1
res['complement formula'] = (300, bad)
# Lemma 5
bad = 0; tot = 0
for _ in range(300):
    gg = rand_graph(rng.randint(1, 5), rng.uniform(0.1, 0.9), rng); t = rng.randint(1, 2)
    J = E2join(gg, t); tot += 1
    if dcount(J) != 2**J[0] - 2**gg[0] - 2 * t + dcount(gg) or is_even(J) != is_even(gg): bad += 1
res['Lemma 5 (E2-multi-join)'] = (tot, bad)
# Lemma 7
bad = 0; tot = 0
for _ in range(3000):
    n = rng.randint(3, 8); g = rand_graph(n, rng.uniform(0.15, 0.7), rng)
    if not (has_leaf(g) and no_universal(g)) or not is_even(g): continue
    q = rng.choice([v for v in range(n) if degs(g)[v] == 1])
    h = SR(g, q); tot += 1
    if not (dcount(h) == 3 * 2**(n - 1) + dcount(g) and is_even(h) and has_leaf(h) and no_universal(h)): bad += 1
    if tot >= 300: break
res['Lemma 7 (stable recursion)'] = (tot, bad)
# Lemma 8
bad = 0; tot = 0; nev = 0
for _ in range(2000):
    rho = rng.randint(3, 6); H = rand_graph(rho, rng.uniform(0.1, 0.6), rng); z = rng.randrange(rho)
    if degs(H)[z] >= rho: continue
    if rho - degs(H)[z] < 1: continue
    G = jump(H, z); dg = sorted(degs(G))
    if dg[-1] == dg[-2]: continue
    tot += 1
    stab = [p for p in automorphisms(H) if p[z] == z]
    coeven = all(edge_sign(H, p) == perm_sign(p) for p in stab)
    ev = is_even(G); nev += ev
    if dcount(G) != 2**(G[0] - 1) - Dcount(H) or ev != coeven or not (has_leaf(G) and no_universal(G)): bad += 1
    if tot >= 400: break
res['Lemma 8 (jump), even cases %d' % nev] = (tot, bad)
for k, (t, b) in res.items():
    print("%-40s tested %4d  failures %d" % (k, t, b))
