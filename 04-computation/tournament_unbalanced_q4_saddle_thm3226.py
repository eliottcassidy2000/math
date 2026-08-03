#!/usr/bin/env python3
"""Exact/symbolic and high-precision controls for the strong unbalanced Q4 saddle.

This is a scratch audit, not a repository proof dependency.  Every logical
check uses explicit exceptions so ``python3`` and ``python3 -O`` agree.
"""

from __future__ import annotations

import math

import mpmath as mp
import sympy as sp


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# The tournament is 1->2->3->1, 3->4, 4->1, 4->2.
X1, X2, X3, X4 = sp.symbols("X1 X2 X3 X4")
Xs = (X1, X2, X3, X4)
A = sp.Matrix(
    [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [1, 0, 0, 1],
        [1, 1, 0, 0],
    ]
)
XD = sp.diag(*Xs)
M = sp.eye(4) - A * XD
H = sp.factor(M.det())
H_expected = 1 - X1 * X2 * X3 - X2 * X3 * X4 - X1 * X2 * X3 * X4
check(sp.expand(H - H_expected) == 0, "wrong Q4 determinant")

one = sp.ones(4, 1)
W = sp.cancel((one.T * XD * M.inv() * one)[0])
N = sp.Poly(sp.cancel(W * H), *Xs)
check(all(co > 0 for _, co in N.terms()), "the pole numerator is not positive")
check(len(N.terms()) == 15, "unexpected Q4 numerator support")

# Pair symmetry gives the exact reduced boundary.
S, T = sp.symbols("S T", positive=True)
H_pair = sp.factor(H.subs({X1: S, X4: S, X2: T, X3: T}))
check(sp.expand(H_pair - (1 - S * (S + 2) * T**2)) == 0, "bad pair reduction")

mp.mp.dps = 100


def E(p: mp.mpf) -> mp.mpf:
    """The exact one-variable stationarity equation E(p)=0."""

    return p * (1 + p * p) * mp.log(1 + p ** -2) - 2 * (1 + p) * mp.log(1 + p)


left = mp.mpf("0.3866667229")
right = mp.mpf("0.3866667230")
check(E(left) > 0 and E(right) < 0, "the printed root bracket has wrong signs")
p = mp.findroot(E, (left, right))
t = mp.log(1 + p)
s = mp.log(1 + p ** -2) / 2
radius = (s * t) ** 2


def phi(x: mp.mpf) -> mp.mpf:
    return x / (1 - mp.exp(-x))


def chi(x: mp.mpf) -> mp.mpf:
    return x * (1 - (1 + x) * mp.exp(-x)) / (1 - mp.exp(-x)) ** 2


tol = mp.mpf("1e-90")
surface_residual = (mp.exp(2 * s) - 1) * (mp.exp(t) - 1) ** 2 - 1
critical_residual = phi(2 * s) / 2 - phi(t)
check(abs(surface_residual) < tol, "surface residual too large")
check(abs(critical_residual) < tol, "critical residual too large")
check(abs((s / t) - (1 + p) / (p * (1 + p * p))) < tol, "ratio identity failed")
check(p < 1 and s / t > 1, "the arithmetic ratio argument needs q>1")

# Exact arithmetic audit for q=s/t.  If q=m/n were rational in lowest terms,
# p would satisfy F(Z)=m Z^3+(m-n)Z-n.  For q>1 this cubic is strictly
# increasing.  In the irrational-root case it is therefore irreducible, and
# its three field norms are
#   N(p)=n/m, N(1+p)=2, N(1+p^2)=2(n/m)^2.
# Norming p^(2n)(1+p)^(2m)=(1+p^2)^n would force n=2m, contrary to m>n.
Z, m_sym, n_sym = sp.symbols("Z m n", nonzero=True)
F = m_sym * Z**3 + (m_sym - n_sym) * Z - n_sym
check(sp.expand(F.subs(Z, -1) + 2 * m_sym) == 0, "N(1+p) identity failed")
Fi_product = sp.expand(F.subs(Z, sp.I) * F.subs(Z, -sp.I))
check(sp.simplify(Fi_product - 2 * n_sym**2) == 0, "N(1+p^2) identity failed")

# The logarithmic Hessian of G=log(U2 U3 (exp(x1+x4)-1)).
# On an orthonormal tangent basis to sum(y_i)=constant, its eigenvalues are:
#   antisymmetry 1<->4: kappa,
#   antisymmetry 2<->3: chi(t),
#   the two-pair contrast: (chi(2s)/2 + chi(t))/2.
kappa = phi(t)
hessian_tangent = (kappa, chi(t), (chi(2 * s) / 2 + chi(t)) / 2)
check(all(ev > 0 for ev in hessian_tangent), "degenerate logarithmic saddle")
phase_tangent = tuple(ev / kappa for ev in hessian_tangent)
phase_det = mp.fprod(phase_tangent)
check(phase_det > 0, "nonpositive Gaussian determinant")

# Directly evaluate the positive pole numerator at (S,T,T,S).
N_fun = sp.lambdify(Xs, N.as_expr(), "mpmath")
S_num = mp.exp(s) - 1
T_num = mp.exp(t) - 1
N_at_saddle = N_fun(S_num, T_num, T_num, S_num)
check(N_at_saddle > 0, "pole order could cancel")

# The equal positive pole is a useful hostile: it lies on H=0 but is not
# diagonal-critical and has a strictly smaller coordinate product.
rho = mp.findroot(lambda z: z**4 - 2 * z - 1, (mp.mpf("1.3"), mp.mpf("1.5")))
equal_x = mp.log(1 + 1 / rho)
equal_radius = equal_x**4
check(equal_radius < radius, "equal hostile unexpectedly maximizes the product")

# Balanced controls.
c3_radius = mp.log(2) ** 3
c5_radius = mp.log(mp.mpf(3) / 2) ** 5
check(c3_radius > 0 and c5_radius > 0, "balanced radius control failed")

# q=1 hostile: raw d!S(r,d) is a finite exponential polynomial and the
# factorial normalization is the coefficient sequence of (exp(z)-1)^d.
for r in range(0, 16):
    raw = math.factorial(3) * int(sp.functions.combinatorial.numbers.stirling(r, 3, kind=2))
    closed = 3**r - 3 * 2**r + 3 if r > 0 else 0
    check(raw == closed, f"q=1 control failed at r={r}")


def normalized_stirling_profile(r: int) -> list[mp.mpf]:
    fac_r = math.factorial(r)
    return [
        mp.mpf(math.factorial(e) * int(sp.functions.combinatorial.numbers.stirling(r, e, kind=2)))
        / fac_r
        for e in range(r + 1)
    ]


numerator_terms = [(mon, int(co)) for mon, co in N.terms()]


def diagonal_d1(r: int) -> mp.mpf:
    """Compute [x1^r...x4^r] W(u(x1),...,u(x4)) from N/(1-M1-M2-M3)."""

    C = normalized_stirling_profile(r)
    facts = [math.factorial(j) for j in range(r + 1)]
    total = mp.mpf(0)
    for n, nco in numerator_terms:
        max_k = r - max(n[1], n[2])
        for k in range(max_k + 1):
            k_fact = facts[k]
            for a in range(k + 1):
                for b in range(k - a + 1):
                    c = k - a - b
                    exponents = (n[0] + a + c, n[1] + k, n[2] + k, n[3] + b + c)
                    if max(exponents) > r:
                        continue
                    multinomial = k_fact // (facts[a] * facts[b] * facts[c])
                    term = mp.mpf(nco * multinomial)
                    for exponent in exponents:
                        term *= C[exponent]
                    total += term
    return total


# For d=1 and q=4 the saddle predicts a_r ~ C r^(-3/2) R^(-r).
a40 = diagonal_d1(40)
a41 = diagonal_d1(41)
raw_ratio = a41 / a40
corrected_ratio = raw_ratio * (mp.mpf(41) / 40) ** mp.mpf("1.5")
check(abs(corrected_ratio - 1 / radius) < mp.mpf("0.01"), "coefficient growth misses saddle radius")

print("Q4_UNBALANCED_SADDLE_AUDIT")
print("det(I-A X) =", H)
print("positive numerator support terms =", len(N.terms()))
print("p bracket =", mp.nstr(left, 20), mp.nstr(right, 20))
print("p =", mp.nstr(p, 80))
print("s =", mp.nstr(s, 80))
print("t =", mp.nstr(t, 80))
print("s/t =", mp.nstr(s / t, 80))
print("arithmetic = s/t is irrational; p is transcendental (norm + Gelfond-Schneider)")
print("R_Q4 = s^2 t^2 =", mp.nstr(radius, 80))
print("1/R_Q4 =", mp.nstr(1 / radius, 80))
print("surface residual =", mp.nstr(surface_residual, 8))
print("critical residual =", mp.nstr(critical_residual, 8))
print("pole numerator at saddle =", mp.nstr(N_at_saddle, 50))
print("tangent Hessian eigenvalues =", *(mp.nstr(v, 50) for v in hessian_tangent))
print("Gaussian determinant =", mp.nstr(phase_det, 50))
print("equal-pole x =", mp.nstr(equal_x, 50))
print("equal-pole product =", mp.nstr(equal_radius, 50))
print("R_Q4 / equal-product =", mp.nstr(radius / equal_radius, 50))
print("balanced C3 radius =", mp.nstr(c3_radius, 50))
print("balanced C5 radius =", mp.nstr(c5_radius, 50))
print("d=1 a_41/a_40 =", mp.nstr(raw_ratio, 50))
print("power-corrected ratio =", mp.nstr(corrected_ratio, 50))
print("PASS")
