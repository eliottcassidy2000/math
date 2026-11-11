#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A8: is the zero set of D_n modulo 12 (THM-5128 (E), quoted "exact for n <= 48")
provable for every n?  The signed 4 x 4 transfer matrix M of q mod 2 on the states (h_0, h_1):
    M[(a, b) -> (d, a)] = (-1)^(d (a + b)),      S_n = e_(0,0)^T M^n 1,   D_n = e_(0,0)^T M^n u,  u(a, b) = (-1)^a.
If M^12 is a scalar matrix the whole pattern is periodic.
"""
import sympy as sp

states = [(a, b) for a in (0, 1) for b in (0, 1)]
idx = {s: i for i, s in enumerate(states)}
M = sp.zeros(4, 4)
for (a, b) in states:
    for d in (0, 1):
        M[idx[(a, b)], idx[(d, a)]] += (-1) ** (d * (a + b))
print("M =", M.tolist())
lam = sp.symbols("lam")
cp = M.charpoly(lam).as_expr()
print("characteristic polynomial:", sp.factor(cp))
print("M M^T =", (M * M.T).tolist())
for k in (3, 4, 6, 12, 24):
    Mk = M ** k
    scal = Mk[0, 0]
    print("M^%d scalar: %s  (%s)" % (k, Mk == scal * sp.eye(4), scal if Mk == scal * sp.eye(4) else "-"))
one = sp.ones(4, 1)
u = sp.Matrix([(-1) ** a for (a, b) in states])
e = sp.zeros(1, 4)
e[0, idx[(0, 0)]] = 1
v, w = one, u
S, D = [], []
for n in range(1, 25):
    v = M * v
    w = M * w
    S.append((e * v)[0])
    D.append((e * w)[0])
print("S_1..S_24:", S)
print("D_1..D_24:", D)
print("D_n = 0 for n in 1..12 at:", [n for n in range(1, 13) if D[n - 1] == 0])
print("S_n = 0 for n in 1..12 at:", [n for n in range(1, 13) if S[n - 1] == 0])
print("D_(n+12) = -64 D_n and S_(n+12) = -64 S_n for n = 1..12:",
      all(D[n + 12 - 1] == -64 * D[n - 1] and S[n + 12 - 1] == -64 * S[n - 1] for n in range(1, 13)))
