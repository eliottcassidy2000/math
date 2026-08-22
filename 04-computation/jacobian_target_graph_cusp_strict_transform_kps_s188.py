#!/usr/bin/env python3
"""Exact cusp strict-transform compiler for nonlinear target graphs."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, c, r, s, t = sp.symbols("a b c r s t")
coefficients = sp.symbols("p00 p10 p01 p20 p11 p02 p30 p21 p12 p03")
p00, p10, p01, p20, p11, p02, p30, p21, p12, p03 = coefficients

phi = (
    p00
    + p10 * a
    + p01 * b
    + p20 * a**2
    + p11 * a * b
    + p02 * b**2
    + p30 * a**3
    + p21 * a**2 * b
    + p12 * a * b**2
    + p03 * b**3
)

L = 27 * a**2 * c**2 - 18 * a * b * c + b**3 * c + 16 * a - b**2
P = 12 * a - b**2
E = 54 * a**2 * c - 18 * a * b + b**3
require(sp.expand(108 * a**2 * L - P**3 - E**2) == 0, "global cube-plus-square identity")

# Pull back the target graph c=-phi through the cusp normalization
# P=-r^2, E=r^3.  The universal a=0 base component is (b-r)^2; the
# residual bracket is the strict-transform equation.
a_cusp = (b**2 - r**2) / 12
phi_cusp = sp.expand(phi.subs(a, a_cusp))
E_graph = sp.expand(E.subs({a: a_cusp, c: -phi_cusp}))
strict_r = sp.expand(b + 2 * r + sp.Rational(3, 4) * (b + r) ** 2 * phi_cusp)
factor_identity = sp.expand(E_graph - r**3 + (b - r) ** 2 * strict_r / 2)
require(factor_identity == 0, "cusp base-square factorization")

# In s=b+r coordinates, a=s(2b-s)/12 and the strict transform becomes the
# compact compiler G_phi(b,s)=0.
a_s = s * (2 * b - s) / 12
phi_s = sp.expand(phi.subs(a, a_s))
strict_s = sp.expand(2 * s - b + sp.Rational(3, 4) * s**2 * phi_s)
require(sp.expand(strict_r.subs(r, s - b) - strict_s) == 0, "s-coordinate compiler")

# At r=0 (equivalently P=0), the strict transform records exactly the
# intersection with the omitted curve E, parametrized by b!=0 and c=4/(3b).
strict_at_r0 = sp.factor(strict_r.subs(r, 0))
target_on_E = sp.expand(sp.Rational(4, 3) / b + phi.subs(a, b**2 / 12))
require(
    sp.expand(strict_at_r0 - sp.Rational(3, 4) * b**2 * target_on_E) == 0,
    "omitted-curve intersection identity",
)

# The ideal rational section E=0 makes the Jelonek section P=0, but solving
# it as a target graph exhibits the unavoidable a^(-2) denominator.
phi_ideal = sp.factor((b**3 - 18 * a * b) / (54 * a**2))
E_ideal = sp.factor(E.subs(c, -phi_ideal))
L_ideal = sp.factor(L.subs(c, -phi_ideal))
require(E_ideal == 0, "ideal rational section kills E")
require(sp.factor(L_ideal - P**3 / (108 * a**2)) == 0, "ideal section leaves triple P divisor")

# A direct finite-field detector: among every quadratic graph shear over F_7,
# uniform p^2 fibres force all quadratic coefficients to vanish.  This is a
# finite diagnostic only, not a characteristic-zero proof.
def quadratic_balance_scan(p: int) -> tuple[int, int, int]:
    inverse_four = pow(4, -1, p)
    rows: list[tuple[int, int, int, int, int, int]] = []
    for xv in range(p):
        for yv in range(p):
            uv = (1 + xv * yv) % p
            for zv in range(p):
                f1 = (uv**3 * zv + yv * yv * uv * (4 + 3 * xv * yv)) % p
                f2 = (yv + 3 * xv * uv * uv * zv + 3 * xv * yv * yv * (4 + 3 * xv * yv)) % p
                f3 = (2 * xv - 3 * xv * xv * yv - xv**3 * zv) % p
                av = (f1 + inverse_four) % p
                bv = f2
                rows.append((av, bv, av * av % p, av * bv % p, bv * bv % p, f3))

    balanced = 0
    genuinely_quadratic = 0
    for linear_a in range(p):
        for linear_b in range(p):
            for quad_a in range(p):
                for quad_ab in range(p):
                    for quad_b in range(p):
                        histogram = [0] * p
                        for av, bv, a2v, abv, b2v, f3v in rows:
                            value = (
                                f3v
                                + linear_a * av
                                + linear_b * bv
                                + quad_a * a2v
                                + quad_ab * abv
                                + quad_b * b2v
                            ) % p
                            histogram[value] += 1
                        if all(count == p * p for count in histogram):
                            balanced += 1
                            if (quad_a, quad_ab, quad_b) != (0, 0, 0):
                                genuinely_quadratic += 1
    return p, balanced, genuinely_quadratic


scan = quadratic_balance_scan(7)
require(scan == (7, 19, 0), "F7 quadratic-shear balance detector")

print("target-graph Jelonek cusp strict-transform compiler")
print("108*a^2*L=(12*a-b^2)^3+E^2")
print("P=-r^2,E=r^3,a=(b^2-r^2)/12")
print("E_graph-r^3=-(b-r)^2/2*[b+2r+3/4*(b+r)^2*phi]")
print("s=b+r: G_phi=2s-b+3/4*s^2*phi(s(2b-s)/12,b)")
print("G_phi|r=0=(3/4)*b^2*(4/(3b)+phi(b^2/12,b))")
print("ideal rational graph: phi=(b^3-18ab)/(54a^2), L=P^3/(108a^2)")
print("F7 quadratic balance: total=19, genuinely_quadratic=0")
print("all active truth gates passed")
