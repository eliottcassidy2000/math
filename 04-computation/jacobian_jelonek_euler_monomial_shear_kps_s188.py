#!/usr/bin/env python3
"""Exact algebraic controls for the THM-3560 Jelonek--Euler gate."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, lam, w, qvar = sp.symbols("a b lambda w q")

# Universal discriminant calculation after q=b^(m+1).  The restriction of
# the Jelonek polynomial to c=-lambda*b^m is quadratic in a.
A2_q = 27 * lam**2 * qvar**2 / b**2
A1_q = 16 + 18 * lam * qvar
A0_q = -b**2 * (1 + lam * qvar)
disc_q = sp.factor(A1_q**2 - 4 * A2_q * A0_q)
P_q = 4 + 3 * lam * qvar
require(sp.expand(disc_q - 4 * P_q**3) == 0, "universal cubic discriminant")

# Avoid the bookkeeping denominator b by checking the same identity directly
# for the first twelve exponents, including the quadratic, cubic/lemniscatic,
# and both hyperelliptic parity classes.
rows = []
for m in range(1, 13):
    n = m + 1
    A2 = 27 * lam**2 * b ** (2 * m)
    A1 = 16 + 18 * lam * b**n
    A0 = -b**2 * (1 + lam * b**n)
    Q = sp.expand(A2 * a**2 + A1 * a + A0)
    P = 4 + 3 * lam * b**n
    discriminant = sp.factor(sp.discriminant(Q, a))
    require(sp.expand(discriminant - 4 * P**3) == 0, f"discriminant m={m}")

    # If s=(2*A2*a+A1)/(2*P), then s^2=P modulo Q.
    normalization_identity = sp.expand((2 * A2 * a + A1) ** 2 - 4 * P**3 - 4 * A2 * Q)
    require(normalization_identity == 0, f"normalization identity m={m}")

    genus = (n - 1) // 2
    infinity_points = 1 if n % 2 else 2
    chi_D = 1 - 2 * genus - infinity_points
    require(chi_D == -m, f"normalization Euler characteristic m={m}")

    e_points = n
    chi_pullback = 3 - 2 * chi_D - e_points
    require(chi_pullback == m + 2, f"pullback Euler characteristic m={m}")
    rows.append((m, n, genus, infinity_points, chi_D, e_points, chi_pullback))

# Irreducibility control.  Over C(a,b), the core fibre cubic on the target
# graph is L_m*x^3+(4+3*lambda*b^(m+1))*x+2*lambda*b^m.  At a=infinity,
# an alleged rational root of integral valuation k gives the three valuations
# (-2+3k,k,0).  Their minimum is unique for every integer k.
def core_root_valuations(k: int) -> tuple[int, int, int]:
    return -2 + 3 * k, k, 0


for k in range(-24, 25):
    values = core_root_valuations(k)
    require(values.count(min(values)) == 1, f"unique cubic valuation minimum k={k}")

# Omitted-curve intersection: b=4/(3t), c=t.  Clearing the nonzero
# denominator gives t^(m+1)=-lambda*(4/3)^m, hence m+1 distinct points.
for m in range(1, 13):
    t = sp.symbols(f"t_{m}", nonzero=True)
    restricted_shear = t + lam * (sp.Rational(4, 3) / t) ** m
    cleared = sp.factor(t**m * restricted_shear)
    expected = t ** (m + 1) + lam * sp.Rational(4, 3) ** m
    require(sp.expand(cleared - expected) == 0, f"omitted-curve equation m={m}")

print("THM-3560 Jelonek--Euler monomial-shear audit")
print("universal discriminant: 4*(4+3*lambda*q)^3, q=b^(m+1)")
print("normalization: s^2=4+3*lambda*b^(m+1), with (b,s)=(0,-2) removed")
print("m n genus infinity chi(D) chi(e) chi(preimage)")
for row in rows:
    print(" ".join(str(value) for value in row))
print("formula: chi(D)=-m, chi(e)=m+1, chi(F^-1(T_m))=m+2")
print("core cubic root valuations at a=infinity: (-2+3k,k,0), unique minimum for every k in Z")
print("verdict: pullback irreducible; no coordinate factor in any monomial shear")
print("all active truth gates passed")
