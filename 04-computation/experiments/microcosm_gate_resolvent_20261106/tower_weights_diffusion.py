#!/usr/bin/env python3
# HYP-9490 lemmas (mac-mini-2026-11-06-microcosm): tower weights U_(e,r) = (L^r q^e)(1) of the gate operator
# L f(q) = q[f(q+1) + q f(q) - f(q-1)] (Session210), checked exactly against
#   (i)   the first-step recursion U_(e,r) = U_(e+2,r-1) + sum_k 2 C(e,2k+1) U_(e-2k,r-1), U_(e,0) = 1;
#   (ii)  the diffusion form U_(e,r) = d^e/dy^e [M^r e^y] at y = 0, M = d^2/dy^2 + 2 sinh(y) d/dy
#         (M e^(ky) = k e^((k+1)y) + k^2 e^(ky) - k e^((k-1)y), so M^r e^y = sum_k W_r(k) e^(ky));
#   (iii) monotonicity: U_(e,r) <= U_(e+2,r) (fixed r) and U_(e,r) >= U_(e+2,r-1) (fixed N = e + 2r);
#   (iv)  sandwich 1 <= U_(d,r)/U_(d+2,r-1) <= 1 + 2^d rho_(d+2,r-1), rho_(d,r) = U_(d-2,r)/U_(d,r);
#   (v)   tower OCF in diffusion form: H(Gamma^r T) = [P_T(d/dy) M^r e^y]_(y=0) for the cores C3 and TT3,
#         against the direct value recursion P_Gamma(j) = j[P(j+1) + j P(j) - P(j-1)].
# Exact integer arithmetic throughout. Output: tower_weights_diffusion.out
from math import comb
from fractions import Fraction

RMAX, DMAX = 120, 12
EMAX = DMAX + 2 * RMAX + 2

# (i) table by the first-step recursion
U = [[1] * (EMAX + 1)]
for r in range(1, RMAX + 1):
    prev = U[-1]
    row = []
    for e in range(EMAX + 1 - 2 * r):
        v = prev[e + 2]
        k = 0
        while 2 * k + 1 <= e:
            v += 2 * comb(e, 2 * k + 1) * prev[e - 2 * k]
            k += 1
        row.append(v)
    U.append(row)

# (ii) diffusion form via exponential coefficients W_r(k): M acts on sum_k W(k) e^(ky)
def diffusion_weights(r):
    W = {1: 1}
    for _ in range(r):
        V = {}
        for k, w in W.items():
            if k == 0:
                continue  # M e^0 = 0
            V[k + 1] = V.get(k + 1, 0) + k * w
            V[k] = V.get(k, 0) + k * k * w
            V[k - 1] = V.get(k - 1, 0) - k * w
        W = {k: w for k, w in V.items() if w}
    return W

bad = 0
for r in range(0, 31):
    W = diffusion_weights(r)
    for e in range(0, DMAX + 1):
        val = sum(w * k ** e for k, w in W.items())  # e-th derivative at 0
        if val != U[r][e]:
            bad += 1
print(f"(ii) diffusion form d^e[M^r e^y](0) = U_(e,r) for r <= 30, e <= {DMAX}: mismatches = {bad}")

# (iii) monotonicity
m1 = sum(1 for r in range(RMAX + 1) for e in range(len(U[r]) - 2) if U[r][e] > U[r][e + 2])
m2 = sum(1 for r in range(1, RMAX + 1) for e in range(len(U[r])) if U[r][e] < U[r - 1][e + 2])
print(f"(iii) violations of U_(e,r) <= U_(e+2,r): {m1};  of U_(e,r) >= U_(e+2,r-1): {m2}  (r <= {RMAX})")

# (iv) sandwich
viol = 0
for r in range(2, RMAX + 1):
    for d in range(2, DMAX + 1):
        lo = Fraction(U[r][d], U[r - 1][d + 2])
        hi = 1 + 2 ** d * Fraction(U[r - 1][d], U[r - 1][d + 2])
        if not (1 <= lo <= hi):
            viol += 1
print(f"(iv) violations of 1 <= U_(d,r)/U_(d+2,r-1) <= 1 + 2^d rho_(d+2,r-1): {viol}")

# (v) tower OCF in diffusion form for C3 (P = q^3 + 2q) and TT3 (P = q^3)
def value_recursion(P, r, J=None):
    J = J or (r + 3)
    vals = [P(j) for j in range(J + r + 2)]
    for _ in range(r):
        vals = [j * (vals[j + 1] + j * vals[j] - vals[j - 1]) if j >= 1 else 0 for j in range(len(vals) - 1)]
    return vals[1]

okv = 0
for r in range(0, 25):
    W = diffusion_weights(r)
    for name, P in (("C3", lambda q: q ** 3 + 2 * q), ("TT3", lambda q: q ** 3)):
        lhs = sum(w * P(k) for k, w in W.items())
        okv += lhs == value_recursion(P, r)
print(f"(v) [P_T(d/dy) M^r e^y](0) = H(Gamma^r T) for T in (C3, TT3), r < 25: {okv}/50 agree")

# trend of rho_(d,r) at fixed d, and N^2 rho / log(N)^2
import math
print("\nrho_(d,r) for d = 3, 4 (N = d + 2r):")
for r in (10, 20, 40, 60, 80, 100, 120):
    for d in (3, 4):
        rho = U[r][d - 2] / U[r][d]
        N = d + 2 * r
        print(f"  d={d} r={r:3d} N={N:3d} rho={rho:.4e}  N^2 rho={N*N*rho:8.3f}  N^2 rho/log(N)^2={N*N*rho/math.log(N)**2:6.3f}")
print("\nall done")
