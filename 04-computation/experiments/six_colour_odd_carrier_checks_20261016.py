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
Session opus-2026-10-16-S38.
"""
import json, math
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

json.dump(out, open('05-knowledge/results/six_colour_odd_carrier_checks_20261016.json', 'w'), indent=1)
print('wrote 05-knowledge/results/six_colour_odd_carrier_checks_20261016.json')
