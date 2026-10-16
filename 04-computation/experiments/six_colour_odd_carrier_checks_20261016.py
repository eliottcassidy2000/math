#!/usr/bin/env python3
"""Companion checks for THM-4661 (odd-degree carriers and the six-colour boundary).

1. The counting arithmetic of statement 2(d) and its even-palette analogue: for even q the extremal order
   n = 1 + q + q(q-2)/2 of a minimal failed exterior is excluded by the Gallai count (1) and the parity count of
   step (v).  The script recomputes, for each even q, the least number h of high vertices forced by (1), the range
   of |W_q| allowed by step (v), and reports whether they are incompatible.  Exact rational arithmetic.
2. Certificates that sqrt(-2) and sqrt(-35) are not in the cubic carrier L' = Q(zeta55, theta),
   2 theta^3 + 2 theta^2 - 2 theta - 1 = 0: a prime p = 1 (mod 55) at which the cubic has a root modulo p gives a
   prime of L' with residue field F_p; if -2 (resp. -35) is a non-residue modulo p, it is not a square in L'.
   (Statement 4(a) needs only the quadratic-subfield argument; this is an independent exact check for L'.)
3. The cubic: irreducibility over Q, real roots, discriminant 148 = 4 * 37 (so Q(theta) is a totally real
   non-Galois cubic field, unramified at 5 and 11).
4. EI21 (HoG 51378) contains the hub-15 subframework that carries the fish's genus-one curve C3 (THM-4661 4(c)):
   the four unit 4-cycles at vertex 15 and the edges 10-18, 16-17 are present, and every faithful numerical
   realization (S24 census machinery) satisfies |P|^2 + |Q|^2 = 1, Re(P conj Q) = 0 with P = v13 - v20, Q = v14 - v19.
Session opus-2026-10-16-S38.
"""
import json, math, itertools
from fractions import Fraction as Fr
import sympy as sp

out = {}

# ---------- 1. counting arithmetic ----------
rows = []
for q in range(4, 41, 2):
    n = 1 + q + q * (q - 2) // 2              # extremal order from the counting lemma (2)
    N = q * (q - 2) // 2                      # |W| = |N_2(z)| at a tight vertex
    # Gallai count: (q-1) l <= (5/2)(l-1) + q h, l = n - h  ->  h >= ((q-7/2) n + 5/2) / (2q - 7/2)
    hmin = math.ceil(Fr((2 * q - 7) * n + 5, 4 * q - 7))
    # step (v): (q-3)(N - w) = (q-2) w + 2 e,  e >= 0 integer, w = |W_q| >= hmin - 1
    feasible = [w for w in range(max(0, hmin - 1), N + 1)
                if (q - 3) * N - (2 * q - 5) * w >= 0 and ((q - 3) * N - (2 * q - 5) * w) % 2 == 0]
    rows.append(dict(q=q, n_extremal=n, W=N, h_min=hmin, w_min=hmin - 1, feasible_w=feasible,
                     excluded=not feasible))
    print(f'q={q:2d}: extremal n={n:4d}, |W|={N:3d}, h>={hmin:3d}, |W_q|>={hmin-1:3d}, feasible |W_q|: {feasible}')
assert rows[1]['q'] == 6 and rows[1]['n_extremal'] == 19 and rows[1]['h_min'] == 6 and rows[1]['excluded']
out['counting'] = rows

# ---------- 2. sqrt(-2), sqrt(-35) not in L' ----------
th = sp.symbols('th')
f = 2 * th**3 + 2 * th**2 - 2 * th - 1
certs = {}
for a in (-2, -35):
    for p in sp.primerange(2, 10**6):
        if p % 55 != 1:
            continue
        if not any((2 * r**3 + 2 * r * r - 2 * r - 1) % p == 0 for r in range(p)):
            continue
        if sp.legendre_symbol(a % p, p) == -1:
            root = next(r for r in range(p) if (2 * r**3 + 2 * r * r - 2 * r - 1) % p == 0)
            certs[str(a)] = dict(p=p, theta_root_mod_p=root, legendre=-1)
            print(f'{a} is a non-residue mod p={p} (p = 1 mod 55, cubic root {root} mod p): sqrt({a}) not in L\'')
            break
assert set(certs) == {'-2', '-35'}
out['quadratic_certificates'] = certs

# ---------- 3. the cubic ----------
P = sp.Poly(f, th)
disc = sp.discriminant(f, th)
roots = [float(r) for r in sp.real_roots(f)]
out['cubic'] = dict(irreducible=P.is_irreducible, discriminant=int(disc), factor_disc=str(sp.factorint(int(disc))),
                    real_roots=roots)
print(f'cubic irreducible: {P.is_irreducible}; discriminant {disc} = {sp.factorint(int(disc))}; real roots {roots}')
assert P.is_irreducible and disc == 148 and len(roots) == 3

# ---------- 4. EI21 carries the fish's genus-one curve C3 ----------
import importlib.util
from pathlib import Path
import numpy as np
spec = importlib.util.spec_from_file_location('base', Path('04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py'))
base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
EI21 = [[4, 8, 19], [5, 9, 20], [7, 9, 20], [6, 8, 19], [0, 6, 14], [1, 7, 13], [3, 4, 11], [2, 5, 12], [0, 3, 11], [1, 2, 12],
        [13, 14, 18], [6, 8, 17, 18], [7, 9, 16, 18], [5, 10, 15, 17], [4, 10, 15, 16], [13, 14, 19, 20], [12, 14, 17, 20],
        [11, 13, 16, 19], [10, 11, 12, 19, 20], [0, 3, 15, 17, 18], [1, 2, 15, 16, 18]]      # HoG 51378 (S24 script)
S = base.structure(EI21)
cyc = set(frozenset(c) for c in S['cycles'])
need_cycles = [(10, 13, 15, 14), (13, 15, 19, 17), (14, 15, 20, 16), (15, 19, 18, 20)]
assert all(frozenset(c) in cyc for c in need_cycles), 'hub 15: the four unit 4-cycles'
E = set(map(frozenset, S['E']))
assert frozenset((10, 18)) in E and frozenset((16, 17)) in E, 'the two closing edges'
found, nsol = base.census(S, 120, 1)
faith = [v for v in found.values() if v['mind'] > 1e-6]
rows4 = []
for v in faith:
    z = {int(a): complex(b[0], b[1]) for a, b in v['z'].items()}
    a, b, c, d = (z[k] - z[15] for k in (13, 14, 19, 20))
    P, Q = a - d, b - c
    X = abs(P) ** 2
    res = max(abs(abs(P) ** 2 + abs(Q) ** 2 - 1), abs((P * Q.conjugate()).real),
              abs((a * d.conjugate()).real - (1 - X / 2)), abs((b * c.conjugate()).real - (1 - (1 - X) / 2)))
    rows4.append(dict(X=X, residual=res))
assert len(rows4) >= 5 and max(r['residual'] for r in rows4) < 1e-10 and 0 < min(r['X'] for r in rows4) and max(r['X'] for r in rows4) < 1
out['ei21_hub'] = dict(faithful_samples=len(rows4), X_range=[min(r['X'] for r in rows4), max(r['X'] for r in rows4)],
                       max_residual=max(r['residual'] for r in rows4))
print(f"EI21 hub 15: {len(rows4)} faithful realizations, X = |v13 - v20|^2 in [{out['ei21_hub']['X_range'][0]:.4f}, "
      f"{out['ei21_hub']['X_range'][1]:.4f}], |P|^2 + |Q|^2 = 1 and P perpendicular to Q to {out['ei21_hub']['max_residual']:.1e}")

# ---------- 5. which small four-chromatic triangle-free unit graphs contain the C3 hub pattern ----------
# hub h with neighbours a, b, c, d; rhombus completions w_ab, w_cd, w_ac, w_bd; edges w_ab-w_cd and w_ac-w_bd.
src = open('04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py', encoding='utf-8').read()
GR = {}
for name in ('EI17', 'EI19', 'FISH'):
    i = src.index("'" + name + "': [[")
    j = src.index(']]', i)
    GR[name] = eval(src[i + len(name) + 4:j + 2])
GR['EI21'] = EI21
hubs = {}
for name, ADJ in GR.items():
    adj = {v: set(x) for v, x in enumerate(ADJ)}
    comp = lambda h, x, y: [w for w in adj[x] & adj[y] if w != h]
    hits = []
    for h in range(len(ADJ)):
        for a, b, c, d in itertools.permutations(sorted(adj[h]), 4):
            if a < d:
                for wab in comp(h, a, b):
                    for wcd in comp(h, c, d):
                        if wcd in adj[wab]:
                            for wac in comp(h, a, c):
                                for wbd in comp(h, b, d):
                                    if wbd in adj[wac]:
                                        hits.append([h, a, b, c, d])
    hubs[name] = dict(count=len(hits), hub_vertices=sorted(set(x[0] for x in hits)))
    print(f"{name}: {len(hits)} oriented C3 hub patterns at vertices {hubs[name]['hub_vertices']}")
out['c3_hub_patterns'] = hubs

json.dump(out, open('05-knowledge/results/six_colour_odd_carrier_checks_20261016.json', 'w'), indent=1)
print('wrote 05-knowledge/results/six_colour_odd_carrier_checks_20261016.json')
