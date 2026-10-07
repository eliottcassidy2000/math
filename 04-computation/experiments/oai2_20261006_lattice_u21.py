#!/usr/bin/env python3
"""The triangular lattice attains u(21) = 57 (correction to THM-431 claim (2); MISTAKE-575).  mac-mini, 2026-10-06.

Eisenstein integers a + b w, w = e^(i pi/3), norm N(a + b w) = a^2 + ab + b^2.  alpha = 2 + w and alphabar = 3 - w have
norm 7, and their quotient has argument arccos(11/14).  For unit-distance graphs G, H embedded in Z[w] (distance sqrt7)
along alpha and alphabar, the Minkowski sum alpha*G + alphabar*H has |G||H| points and at least
|G| e(H) + |H| e(G) pairs at squared distance 7 (the Erdos product).
  1. G = triangle {0, 1, w}, H = centred hexagon W_6 = {0} u units: 21 points, exactly 57 = u(21) (Alexeev-Mixon-Parshall
     2024) pairs at distance sqrt7.  So the triangular lattice is optimal at N = 21 (THM-431 (2) said 47, gap 10).
  2. G = H = W_6: 49 points with exactly 168 > 3*49 = 147 pairs (a lattice set beating 3N).
Run: python3 oai2_20261006_lattice_u21.py
"""
from itertools import combinations
OK = True
def check(c, m):
    global OK
    print(("  ok   " if c else "  FAIL ") + m); OK &= bool(c)
def mul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c + b * d)       # w^2 = w - 1
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def nrm(x): a, b = x; return a * a + a * b + b * b
alpha, alphab = (2, 1), (3, -1)
units = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
T = [(0, 0), (1, 0), (0, 1)]
W6 = [(0, 0)] + units
def edges(S, D=7):
    return sum(1 for x, y in combinations(S, 2) if nrm((x[0] - y[0], x[1] - y[1])) == D)
check(nrm(alpha) == 7 == nrm(alphab), "N(2 + w) = N(3 - w) = 7")
# cos of the angle between alpha and alphabar: Re(alpha * conj(alphabar)) / 7
import math
def to_c(x): return complex(x[0] + x[1] * 0.5, x[1] * math.sqrt(3) / 2)
c = (to_c(alpha) * to_c(alphab).conjugate()).real / 7
check(abs(c - 11 / 14) < 1e-12, f"cos(angle(alpha, alphabar)) = {c:.12f} = 11/14")
S21 = {add(mul(alpha, p), mul(alphab, q)) for p in T for q in W6}
check(len(S21) == 21 and edges(S21) == 57, f"triangle x W_6: {len(S21)} points, {edges(S21)} pairs at distance sqrt7 = u(21) = 57")
check(edges([mul(alpha, p) for p in T]) == 3 and edges([mul(alphab, q) for q in W6]) == 12, "factors: triangle 3 edges, W_6 12 edges; 3*12 + 7*3 = 57")
S49 = {add(mul(alpha, p), mul(alphab, q)) for p in W6 for q in W6}
check(len(S49) == 49 and edges(S49) == 168, f"W_6 x W_6: {len(S49)} points, {edges(S49)} pairs (> 3N = 147)")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
