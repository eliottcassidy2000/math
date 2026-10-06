#!/usr/bin/env python3
"""Part 3: the trace sequence a_p = #{roots of P mod p} - 1 = chi_std(Frob_p) of the 7-dimensional standard Artin
representation of S8; its moments against the S8 predictions (mean 0, second moment 1, third moment 1, fourth 4,
i.e. the moments of (fixed points - 1) for a uniform permutation of 8 letters), the tilt angle, and small primes.

Reproduce: python3 04-computation/experiments/square11_octic_field_part3_20261006.py [prime_bound]
"""
import sys, math
from itertools import permutations
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
import flint

PB = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
Pc = [1, -20, 178, -842, 1923, -496, -6754, 12420, -6865]
disc = -209869400103889499760776960

# S8 predictions for moments of a = fix - 1
fix = [sum(1 for i, x in enumerate(p) if i == x) for p in permutations(range(8))]
pred = {k: Fr(sum((f - 1) ** k for f in fix), len(fix)) for k in range(1, 5)}
print('S8 moments of a = fix - 1:', {k: str(v) for k, v in pred.items()})

S = {k: 0 for k in range(1, 5)}
n = 0
small = []
for p in sp.primerange(3, PB):
    if disc % p == 0:
        continue
    g = flint.nmod_poly([c % p for c in reversed(Pc)], p)
    roots = sum(e for f, e in g.factor()[1] if f.degree() == 1)
    a = roots - 1
    for k in S:
        S[k] += a ** k
    n += 1
    if p < 60:
        small.append((p, a, sorted(f.degree() for f, e in g.factor()[1] for _ in range(e))))
print(f'primes 3 <= p < {PB}, p not dividing disc: {n}')
print('observed moments:', {k: round(S[k] / n, 5) for k in S})
print('small primes (p, a_p, factor degrees):', small)

mp.mp.dps = 40
u = mp.findroot(lambda x: 5*x**8 - 10*x**7 - 2*x**6 + 14*x**5 + 12*x**4 - 6*x**3 + 2*x**2 + 2*x - 1, 0.3657)
theta = 2 * mp.atan(u)
T = (6*u + 4) / (1 + 2*u - u**2)
T2 = (2 + 2*mp.cos(theta) + 3*mp.sin(theta)) / (mp.sin(theta) + mp.cos(theta))
print('u =', u)
print('tilt angle 2 atan(u) =', mp.degrees(theta), 'degrees')
print('T =', T, '   (2 + 2 cos t + 3 sin t)/(sin t + cos t) =', T2)
