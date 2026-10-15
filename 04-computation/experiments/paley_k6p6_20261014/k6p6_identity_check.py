#!/usr/bin/env python3
"""HYP-9276 check at every prime 5 <= q <= 30000 (k6p6_fft.json from k6p6_fft.py; g_trace_a2_half.txt from identify_g.gp).

Identity tested (exact integers):  M_q = -a_(q^2)(g) - 2q - 1,   a_(q^2)(g) = a_q(g)^2 - chi_(-4)(q) q^2 = chi_(-4)(q)(|a_q|^2 - q^2),
equivalently  W_q = |a_q(g)|^2 - q^2 + 7 + 14 chi_(-4)(q) = chi_(-4)(q)(a_(q^2)(g) + 14) + 7  and  V_q = q(q-1) W_q,
with g = q + 4i q^3 + 2q^5 - 8i q^7 - ... the newform of S_3^new(Gamma_0(32), (-4/.)) (McCarthy-Papanikolas f_2).
g_trace_a2_half.txt holds Tr_(Q(i)/Q)(a_p^2)/2 = Re(a_p^2) = chi_(-4)(p)|a_p|^2 for the primes 5..30000 in order.
Also: Riordan-number moments of the normalised adjoint trace x_q = (|a_q|^2 - q^2)/q^2 (SO(3) Sato-Tate), the zero set of M_q,
and the 2-adic valuations of M_q by residue class mod 8.
"""
import json
import sympy
from fractions import Fraction

fft = {int(k): v for k, v in json.load(open('k6p6_fft.json')).items()}
txt = open('g_trace_a2_half.txt').read().strip().strip('[]')
re_a2 = [int(s) for s in txt.split(',')]
P = list(sympy.primerange(5, 30001))
assert len(P) == len(re_a2) == len(fft), (len(P), len(re_a2), len(fft))
bad = []
xs = {1: [], 3: []}
zeros = []
v2min = {}
for q, r in zip(P, re_a2):
    eps = 1 if q % 4 == 1 else -1
    abs2 = eps * r                       # |a_q|^2
    a_q2 = eps * (abs2 - q * q)          # a_(q^2)(g)
    M, W = fft[q]
    if M != -a_q2 - 2 * q - 1 or W != abs2 - q * q + 7 + 14 * eps or W != eps * (a_q2 + 14) + 7:
        bad.append(q)
    xs[q % 4].append(Fraction(abs2 - q * q, q * q))
    if M == 0:
        zeros.append(q)
    if M:
        v = (M & -M).bit_length() - 1
        v2min[q % 8] = min(v2min.get(q % 8, 99), v)
print(f"HYP-9276: M_q = -a_(q^2)(g) - 2q - 1 and W_q = |a_q|^2 - q^2 + 7 + 14 chi(q) at {len(P)} primes 5..30000; failures: {bad}")
riordan = [1, 0, 1, 1, 3, 6, 15, 36, 91]
for cls in (3, 1):
    v = xs[cls]
    mom = [float(sum(x ** k for x in v) / len(v)) for k in range(9)]
    print(f"  q = {cls} mod 4 ({len(v)} primes): moments of x_q = (|a_q|^2 - q^2)/q^2, k = 0..8: "
          + ", ".join(f"{m:.3f}" for m in mom) + f"   (SO(3) Sato-Tate / Riordan: {riordan})")
print(f"  M_q = 0 exactly at {zeros}: all q = 7 mod 8, and there |a_q(g)| = q + 1: "
      f"{all(q % 8 == 7 for q in zeros)}")
print(f"  min v_2(M_q) by q mod 8: {dict(sorted(v2min.items()))}")
