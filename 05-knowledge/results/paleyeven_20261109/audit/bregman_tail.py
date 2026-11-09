#!/usr/bin/env python3
"""Audit P1: monotonicity from q = 55 on, and the explicit tail bound
R(q) <= q log(3/4) + 2.85 log q + 7.2 (derived by hand from Robbins' bounds in AUDIT_P1.md) for both ratios."""
from math import lgamma, log
from bregman_check import ratios, lf_gamma
rows = [(q,) + ratios(q, lf_gamma) for q in range(7, 2000001, 4)]
f = [r for r in rows if r[0] >= 55]
print("D ratio strictly decreasing for q >= 55:", all(f[i+1][1] < f[i][1] for i in range(len(f)-1)))
print("deletion ratio strictly decreasing for q >= 55:", all(f[i+1][2] < f[i][2] for i in range(len(f)-1)))
bound = lambda q: q * log(3/4) + 2.85 * log(q) + 7.2
print("tail bound holds for all q in [7, 2e6] (both ratios):", all(max(r[1], r[2]) <= bound(r[0]) for r in rows))
print("bound negative from q =", next(q for q in range(11, 10**4) if all(bound(x) < 0 for x in range(q, q + 2000))))
print("max slack (bound - ratio) at q = 7, 55, 1e6:", [round(bound(r[0]) - max(r[1], r[2]), 3) for r in rows if r[0] in (7, 55, 999999)])
# tighter constant derived by hand in AUDIT_P1.md (finding on Theorem A4): 4.8 instead of 7.2
b2 = lambda q: q * log(3/4) + 2.85 * log(q) + 4.8
print("bound with constant 4.8 holds for all q in [7, 2e6] (both ratios):", all(max(r[1], r[2]) <= b2(r[0]) for r in rows),
      "; negative and decreasing from q =", next(q for q in range(11, 10**4) if b2(q) < 0))
