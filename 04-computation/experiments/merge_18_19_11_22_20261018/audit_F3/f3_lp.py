#!/usr/bin/env python3
"""Audit F3, item (b): the two-point (Delsarte / Schrijver theta') LP for alpha(G_q), written independently.

f is constant on the norm circles N(z) = c (orbits of z -> lambda z, lambda in mu_(q+1)); f(0) = 1, f = 0 on N = 1,
f >= 0, and the Fourier transform at a character of norm class m is 1 - sum_c f(c) Kl(c m) >= 0, because
sum_(N(z) = c) psi(Tr(a z)) = -Kl(c N(a)) (THM-4628(b); checked numerically in f3_kappa_check.py).
alpha <= 1 + (q+1) sum_c f(c).  Compared with the Hoffman ratio bound n (-lmin)/(d - lmin).
Usage: python3 f3_lp.py q [q ...]"""
import sys, math, os
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
for q in [int(a) for a in sys.argv[1:]]:
    inv = lambda y: pow(y, q - 2, q)
    Kl = {m: sum(math.cos(2 * math.pi * ((y + m * inv(y)) % q) / q) for y in range(1, q)) for m in range(1, q)}
    cs = list(range(2, q))                     # norm classes c != 0, 1
    A = [[-Kl[(c * m) % q] for c in cs] for m in range(1, q)]   # 1 - sum f(c) Kl(cm) >= 0  <=>  -sum ... >= -1
    A = [[-x for x in row] for row in A]       # sum f(c) Kl(cm) <= 1
    res = linprog(c=[-1.0] * len(cs), A_ub=A, b_ub=[1.0] * (q - 1), bounds=[(0, None)] * len(cs), method='highs')
    lp = 1 + (q + 1) * (-res.fun)
    lmin = min(-Kl[m] for m in range(1, q))
    hoff = q * q * (-lmin) / ((q + 1) - lmin)
    out.append('q = %d: two-point LP alpha <= %.4f ; Hoffman ratio alpha <= %.4f ; chi >= q^2/LP = %.4f'
               % (q, lp, hoff, q * q / lp))
    print(out[-1])
open(os.path.join(HERE, 'f3_lp.out'), 'w').write('\n'.join(out) + '\n')
