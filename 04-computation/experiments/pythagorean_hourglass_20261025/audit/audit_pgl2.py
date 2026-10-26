# Audit check E: the PGL(2) lemma of THM-4780 (4).
#  (1) Direct test of the invariance F(g x) = F(x) on P^1(F_q) for EVERY injective configuration (so every placement of
#      vertices at 0 and at infinity is covered), for the generators x+1, a x (a a non-square), -1/x and a random
#      Moebius map, for the octahedron (order orientation: in-degrees 0,0,2,2,4,4) and for K7 - C7 with an even
#      in-degree orientation found by search; and a control: an orientation with odd in-degrees is NOT invariant at
#      q = 3 (mod 4) but is at q = 1 (mod 4).
#  (2) Numerical tests of T_G = eps (q + 1 - k) T_(G-u) on further all-even-degree graphs, both classes mod 4:
#      K5 vs K4; C4 vs P3; K7 - C7 vs K6 - P6; K7 - K3 vs K6 - e and K6 - K3; K7 - (C3 + C4) vs K6 - (C4 + K2) and
#      K6 - (C3 + P3); K7 - C5 vs K6 - P4 and K6 - (C5 + K1).  The sign eps is predicted from the orientation flips.
# Audit 2026-10-25.
import itertools, random, sys
import numpy as np
sys.path.insert(0, '.')
from audit_graph_sums import S_norm, T_full

INF = 'inf'
def chi_tab(p):
    c = [0] * p
    for x in range(1, p): c[(x * x) % p] = 1
    return [0] + [1 if c[x] else -1 for x in range(1, p)]

def moebius(p, M, x):
    a, b, c, d = M
    if x == INF: return INF if c % p == 0 else (a * pow(c, p - 2, p)) % p
    den = (c * x + d) % p
    if den == 0: return INF
    return ((a * x + b) * pow(den, p - 2, p)) % p

def F(p, ch, arcs, x):
    v = 1
    for i, j in arcs:
        if x[i] == INF or x[j] == INF: continue
        v *= ch[(x[j] - x[i]) % p]
        if v == 0: return 0
    return v

def invariance(p, arcs, nv, gens, exhaustive=True, samples=20000, seed=1):
    ch = chi_tab(p); pts = list(range(p)) + [INF]
    rng = random.Random(seed)
    if exhaustive: configs = itertools.permutations(pts, nv)
    else: configs = (tuple(rng.sample(pts, nv)) for _ in range(samples))
    bad = {name: 0 for name in gens}; n = 0
    for x in configs:
        n += 1
        f0 = F(p, ch, arcs, x)
        for name, M in gens.items():
            if F(p, ch, arcs, tuple(moebius(p, M, xi) for xi in x)) != f0: bad[name] += 1
    return n, bad

def nonsq(p):
    ch = chi_tab(p); return next(a for a in range(2, p) if ch[a] == -1)

def even_indegree_orientation(nv, edges):
    """search orientations (as flips of the order orientation) with all in-degrees even; return the flip set"""
    for r in range(len(edges) + 1):
        for F_ in itertools.combinations(range(len(edges)), r):
            arcs = [(j, i) if n in F_ else (i, j) for n, (i, j) in enumerate(edges)]
            indeg = [0] * nv
            for i, j in arcs: indeg[j] += 1
            if all(d % 2 == 0 for d in indeg): return set(F_), arcs
    return None, None

out = []
OCT = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) not in [(0, 1), (2, 3), (4, 5)]]
C7 = [(i, (i + 1) % 7) for i in range(7)]
C7s = {tuple(sorted(e)) for e in C7}
K7mC7 = [(i, j) for i in range(7) for j in range(i + 1, 7) if (i, j) not in C7s]
flipsC7, arcsC7 = even_indegree_orientation(7, K7mC7)
out.append(f'K7 - C7 (complement cycle 0-1-...-6-0): an even in-degree orientation flips the order arcs {[K7mC7[n] for n in sorted(flipsC7)]}')
odd_or = [(i, j) for (i, j) in OCT]; odd_or[0] = (odd_or[0][1], odd_or[0][0])     # flip one arc: in-degrees become odd
for p in (7, 5, 11, 13):
    a = nonsq(p)
    gens = {'x+1': (1, 1, 0, 1), f'{a}x': (a, 0, 0, 1), '-1/x': (0, -1 % p, 1, 0), 'random': (2, 3, 1, 5 if p != 7 else 4)}
    det = (gens['random'][0] * gens['random'][3] - gens['random'][1] * gens['random'][2]) % p
    assert det != 0
    ex = p <= 7
    n, bad = invariance(p, OCT, 6, gens, exhaustive=ex)
    out.append(f'(1) q = {p} (mod 4 = {p % 4}) octahedron, order orientation: {n} configurations ({"all" if ex else "random"}), non-invariant: {bad}')
    n, bad = invariance(p, odd_or, 6, gens, exhaustive=ex)
    out.append(f'    control (one arc flipped, odd in-degrees): non-invariant: {bad}')
    n, bad = invariance(p, arcsC7, 7, gens, exhaustive=ex)
    out.append(f'    K7 - C7, even in-degree orientation: {n} configurations, non-invariant: {bad}')

# (2) identities
def sgn_pred(flips, edges, u):
    """eps for T_(G, order) = eps (q+1-k) T_(G-u, order) when the even orientation flips the arcs in `flips`"""
    return (-1) ** sum(1 for n in flips if u in edges[n])

def minus7(comp):
    cs = {tuple(sorted(e)) for e in comp}
    return [(i, j) for i in range(7) for j in range(i + 1, 7) if (i, j) not in cs]
def relabel_delete(edges, u):
    return [(i - (i > u), j - (j > u)) for (i, j) in edges if u not in (i, j)]
tests = {
    'K7-C7': (minus7(C7), [6]),
    'K7-K3': (minus7([(0, 1), (1, 2), (0, 2)]), [0, 6]),
    'K7-(C3+C4)': (minus7([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (5, 6), (3, 6)]), [0, 6]),
    'K7-C5': (minus7([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]), [0, 6]),
}
K5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for p in (5, 13, 17, 29, 37, 41):
    s5 = S_norm(p, 5, K5); s4 = S_norm(p, 4, [(i, j) for i in range(4) for j in range(i + 1, 4)])
    out.append(f'(2) q = {p}: S_K5 = {s5}, (q-4) S_K4 = {(p - 4) * s4}: {"OK" if s5 == (p - 4) * s4 else "FAIL"}')
C4 = [(0, 2), (0, 3), (1, 2), (1, 3)]
for p in (7, 11, 13, 17):
    t4 = T_full(p, 4, C4); t3 = T_full(p, 3, [(0, 2), (1, 2)])
    out.append(f'(2) q = {p}: T_C4 (two sources, two sinks) = {t4} = -q(q-1)(q-3)? {t4 == -p * (p - 1) * (p - 3)}; (q-3) T_P3 = {(p - 3) * t3}: '
               f'{"OK" if t4 == (p - 3) * t3 else "FAIL"}')
for name, (E, us) in tests.items():
    flips, arcs = even_indegree_orientation(7, E)
    for p in (7, 11, 19, 23, 13, 17):
        sG = S_norm(p, 7, E)
        parts = []
        for u in us:
            Eu = relabel_delete(E, u)
            su = S_norm(p, 6, Eu)
            eps = sgn_pred(flips, E, u) if p % 4 == 3 else 1
            parts.append(f'u={u}: eps(q-6)S_(G-u) = {eps * (p - 6) * su} [{"OK" if sG == eps * (p - 6) * su else "FAIL"}]')
        out.append(f'(2) {name} q = {p} (mod 4 = {p % 4}): S_G = {sG}; ' + '; '.join(parts))
print('\n'.join(out))
with open('audit_pgl2.out', 'w') as fh: fh.write('\n'.join(out) + '\n')
