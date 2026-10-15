#!/usr/bin/env python3
"""A_12 for the Paley tournaments, assembled from the partition-lattice compiler and the core identities.

Input: a12_full.json (written by paley_run_partition_compiler_20261014.py at L = 12: reducible values and core coefficient
vectors at the primes QS).  Core values used:
  K_4, K_5 and the three six-vertex cores with odd-reversal automorphisms: 0 (THM-4652; paley_core_symmetry_20261014.out);
  wheel W_4 (canonical orientation): U = -q(q-1)(q + a_q)  (Ahlgren-Ono corollary; S34);
  octahedron K_{2,2,2} (canonical orientation): U = q(q-1)((5 - q) a_q + 3q)  (FINITE-EXACT at all odd primes 5..107;
  paley_six_vertex_cores_20261014.out),
with a_q the coefficients of eta(2z)^4 eta(4z)^4.  Checks: A_12(7) = A_12(11) = 0 (13 distinct points are needed);
the polynomial part P_12 is interpolated on 8 primes and must predict the remaining primes exactly.
Session opus-2026-10-14-S35.
"""
import json, sys
from fractions import Fraction
import sympy as sp

path = sys.argv[1] if len(sys.argv) > 1 else 'a12_full.json'
D = json.load(open(path))
QS = D['QS']
red = D['reducible']
N = 220
P = [0] * (N + 1); P[0] = 1
for n in range(1, N + 1):
    for step in (2 * n, 4 * n):
        if step > N: continue
        for _ in range(4):
            for i in range(N, step - 1, -1):
                P[i] -= P[i - step]
a = [0] * (N + 2)
for i in range(N + 1): a[i + 1] = P[i]

W4 = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [2, 4], [3, 4]]
OCT = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]]
coefW = [0] * len(QS); coefO = [0] * len(QS); other = []
for edges, n, vec in D['cores']:
    if edges == W4:
        coefW = vec
    elif edges == OCT:
        coefO = vec
    else:
        other.append((n, len(edges)))
print('cores with zero value (K_4, K_5, odd-reversal six-vertex cores):', other)
print('wheel coefficient at QS:', coefW)
print('octahedron coefficient at QS:', coefO)
A12, poly, modc = {}, {}, {}
for j, q in enumerate(QS):
    UW = -q * (q - 1) * (q + a[q])
    UO = q * (q - 1) * ((5 - q) * a[q] + 3 * q)
    A12[q] = red[j] + coefW[j] * UW + coefO[j] * UO
    # split into polynomial part and the coefficient of a_q
    modc[q] = (-coefW[j] * q * (q - 1) + coefO[j] * q * (q - 1) * (5 - q))
    poly[q] = A12[q] - modc[q] * a[q]
print('A_12 at QS:', A12)
print('too-few-points checks: A_12(7) =', A12.get(7), ', A_12(11) =', A12.get(11))
x = sp.symbols('q')
base = QS[:8]
Pp = sp.interpolate([(b, poly[b]) for b in base], x)
ok = all(Pp.subs(x, q) == poly[q] for q in QS[8:])
print('P_12 (degree <= 7 through 8 primes) predicts the remaining primes exactly:', ok)
print('P_12(q) =', sp.factor(sp.expand(Pp)))
Mc = sp.interpolate([(b, modc[b]) for b in base], x)
print('coefficient of a_q:', sp.factor(sp.expand(Mc)), ' (exact at the remaining primes:', all(Mc.subs(x, q) == modc[q] for q in QS[8:]), ')')
print('A_12 at q = 3 from the formula:', sp.expand(Pp).subs(x, 3) + Mc.subs(x, 3) * a[3])
