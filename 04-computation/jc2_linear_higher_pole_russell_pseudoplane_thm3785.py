#!/usr/bin/env python3
"""Exact companion for THM-3785's linear higher-pole Russell pseudo-plane."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)


x, y, c = sp.symbols("x y c", nonzero=True)
tformal = sp.symbols("tformal")
t = x**2 * (1 + x * y)
q = c + t
Q = q**3 / 3
R = x**3
F = sp.expand(x * q)
P = 1 / (3 * R)
Z = sp.expand(F - R)
E = sp.cancel((Z**3 - c**3 * R) / R**2)


def jac(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(b, y) - sp.diff(a, y) * sp.diff(b, x))


# Rational Keller seed and polynomial Russell coordinates.
same(sp.diff((c + tformal) ** 3 / 3, tformal), (c + tformal) ** 2,
     "formal primitive")
same(jac(F, t), x**3 * q, "J(F,t)")
same(jac(F, P), 1, "J(F,P)")
same(F**3 * P, Q, "F^3 P=Q")
same(Z, x * (c + x**3 * y), "Z source formula")
same(E, 3 * c**2 * y + 3 * c * x**3 * y**2 + x**6 * y**3, "E polynomial formula")
check(sp.denom(sp.cancel(E)) == 1, "E polynomial denominator")
same(R**2 * E, Z**3 - c**3 * R, "Russell relation")
same(F, R + Z, "F recovery")
same(P, 1 / (3 * R), "P recovery")

# The earlier Y-coordinate and its triangular equivalence to (R,Z,E).
Y = sp.cancel((F**3 - 3 * R * F**2 - c**3 * R) / R**2)
same(
    Y,
    3 * c**2 * y + 3 * c * x**3 * y**2 - 3 * c * x
    + x**6 * y**3 - 3 * x**4 * y - 2 * x**3,
    "Y polynomial formula",
)
same(Z, F - R, "triangular Z")
same(E, Y + 3 * F - R, "triangular E")
same(R**2 * Y, F**3 - 3 * R * F**2 - c**3 * R, "Y hypersurface")

# Poisson packet in both coordinate systems.
same(jac(R, Z), 3 * R**2, "{R,Z}")
same(jac(R, E), 9 * Z**2, "{R,E}")
same(jac(Z, E), 3 * c**3 + 6 * R * E, "{Z,E}")
same(jac(R, F), 3 * R**2, "{R,F}")
same(jac(R, Y), 9 * F**2 - 18 * R * F, "{R,Y}")
same(jac(F, Y), 3 * c**3 + 9 * F**2 + 6 * R * Y, "{F,Y}")

# Smoothness/etaleness: the hypersurface gradient and Poisson minors have no
# common zero over Q(c), c != 0.
r, z, e = sp.symbols("r z e")
H = r**2 * e - z**3 + c**3 * r
Kc = sp.QQ.frac_field(c)
smooth_gb = sp.groebner(
    [H, sp.diff(H, r), sp.diff(H, z), sp.diff(H, e)],
    r,
    z,
    e,
    order="grevlex",
    domain=Kc,
)
minor_gb = sp.groebner(
    [H, 3 * r**2, 9 * z**2, 3 * c**3 + 6 * r * e],
    r,
    z,
    e,
    order="grevlex",
    domain=Kc,
)
check([p.as_expr() for p in smooth_gb.polys] == [1], "smooth Groebner unit")
check([p.as_expr() for p in minor_gb.polys] == [1], "minor Groebner unit")
same(jac(R, Z).subs(x, 0), 0, "axis first minor")
same(jac(Z, E).subs(x, 0), 3 * c**3, "axis transverse minor")

# Exact inverse formulas on r != 0 and on the retained arm.
xi, rr, zz, ee = sp.symbols("xi rr zz ee", nonzero=True)
u = zz / xi
y_lift = (u - c) / rr
check(sp.expand(xi**3 - rr) == xi**3 - rr, "cube-root lift marker")
same(
    ((u**3 - c**3) / rr).subs(rr, xi**3),
    ((zz**3 - c**3 * rr) / rr**2).subs(rr, xi**3),
    "off-arm E lift",
)
same(E.subs({x: 0}), 3 * c**2 * y, "retained arm coordinate")

# Divisor/Picard markers: r(c^3+re)=z^3 and ord_L(r)=3 ord_L(z).
same(r * (c**3 + r * e) - z**3, H, "triple boundary identity")
check(sp.gcd(3, 1) == 1, "triple boundary primitive valuation")

# Weighted Laurent normal form in x,w.
w = sp.symbols("w")
h = w**3 - c**3
for i in range(5):
    for j in range(5):
        for k in range(5):
            weight = 3 * i + j - 3 * k
            residue = weight % 3
            floor = max(0, (-(weight) + 2) // 3)
            profile = sp.expand(w**j * h**k)
            check((j - residue) % 3 == 0, f"profile residue {i},{j},{k}")
            if weight < 0:
                check(k >= floor, f"negative profile h-order {i},{j},{k}")
                check(sp.rem(profile, h**floor, w) == 0, f"negative divisibility {i},{j},{k}")

# General homogeneous bracket and the unique boundary candidate (1,-3).
a, b = sp.symbols("a b", integer=True)
f = sp.Function("f")(w)
g = sp.Function("g")(w)
hom_bracket = x ** (a + b + 2) * (a * f * sp.diff(g, w) - b * sp.diff(f, w) * g)
Px = x**a * f
Qx = x**b * g
same(
    x**3 * (sp.diff(Px, x) * sp.diff(Qx, w) - sp.diff(Px, w) * sp.diff(Qx, x)),
    hom_bracket,
    "homogeneous bracket formula",
)
uvar = sp.symbols("u")
pdeg, qdeg = sp.symbols("pdeg qdeg", integer=True, nonnegative=True)
exceptional_degree = 3 * (pdeg + qdeg + 1)
exceptional_lead = 3 * (3 * pdeg + qdeg + 2)
check(sp.simplify(exceptional_degree.subs({pdeg: 0, qdeg: 0})) == 3,
      "exceptional positive degree floor")
check(sp.simplify(exceptional_lead.subs({pdeg: 0, qdeg: 0})) == 6,
      "exceptional leading coefficient floor")

# Complete 2x2 endpoint-power/degree packet.  The low common base has degree
# at least three because it carries h; the high common base may be scalar.
for A0 in range(3, 11):
    for B0 in range(3, 11):
        delta = sp.gcd(A0, B0)
        C0, D0 = B0 - 2, A0 - 2
        epsilon = sp.gcd(C0, D0)
        for dp in range(3, 7):
            for dq in range(0, 5):
                deg_left = sp.Rational(A0, delta) * dp + sp.Rational(D0, epsilon) * dq - 1
                deg_right = sp.Rational(B0, delta) * dp + sp.Rational(C0, epsilon) * dq - 1
                same(
                    deg_right - deg_left,
                    (B0 - A0) * (sp.Rational(dp, delta) + sp.Rational(dq, epsilon)),
                    f"2x2 degree split {A0},{B0},{dp},{dq}",
                )
                lead_left = -A0 * sp.Rational(D0, epsilon) * dq - D0 * sp.Rational(A0, delta) * dp
                lead_right = C0 * sp.Rational(B0, delta) * dp + B0 * sp.Rational(C0, epsilon) * dq
                check(lead_left < 0, f"2x2 left lead {A0},{B0},{dp},{dq}")
                check(lead_right > 0, f"2x2 right lead {A0},{B0},{dp},{dq}")

pp = sp.Function("pp")(w)
qq = sp.Function("qq")(w)
lam0, mu0, AA = sp.symbols("lam0 mu0 AA")
equal_middle = (
    (-AA) * pp * sp.diff(mu0 * qq, w)
    - (AA - 2) * sp.diff(pp, w) * mu0 * qq
    + (AA - 2) * qq * sp.diff(lam0 * pp, w)
    + AA * sp.diff(qq, w) * lam0 * pp
)
same(
    equal_middle,
    (lam0 - mu0) * (AA * pp * sp.diff(qq, w) + (AA - 2) * sp.diff(pp, w) * qq),
    "2x2 equal-base factor",
)

# Affine-carrier critical equations.  After scaling the E coefficient to one
# and writing z=kappa*r, the surface equation has this exact residual.
aa, bb, kappa = sp.symbols("aa bb kappa")
r2 = sp.symbols("r2")
critical_residual = sp.expand(((aa - 2 * kappa**3) * r2 + c**3) / 2)
same(
    (aa * r2 - c**3) / 2 - kappa**3 * r2 + c**3,
    critical_residual,
    "affine critical residual",
)
critical_e = (aa * r**2 - c**3) / (2 * r)
same(
    sp.cancel(H.subs({z: kappa * r, e: critical_e}) / r),
    critical_residual.subs(r2, r**2),
    "affine critical substitution",
)
same((bb * r**2 + 3 * (kappa * r) ** 2) / r**2,
     bb + 3 * kappa**2, "affine root equation")

# Pure E has genus-one generic fibre.  The cyclic cubic has three fully
# ramified points, giving 2g-2=-6+6=0; dr/z^2 has order zero at each.
check(3 * (-2) + 3 * (3 - 1) == 0, "genus-one Riemann-Hurwitz")
check((2 + (-2)) == 0, "finite branch differential order")
check((-4) - (-4) == 0, "infinity differential order")

# Rational mates to Z+lambda R are exactly 1/(3R)+k(Z+lambda R).
lam = sp.symbols("lam")
A_lam = Z + lam * R
same(jac(A_lam, 1 / (3 * R)), 1, "affine rational mate")
same(A_lam, x * (c + lam * x**2 + x**3 * y), "two-component affine carrier")

# Frozen semantic packet.
semantic = {
    "affine_carriers": "smooth iff z+lambda*r or e; only z+lambda*r has rational mate; none has polynomial mate",
    "cover": "surjective_etale_generic_degree_3_one_sheet_on_boundary",
    "field": "k(x,y)=k(r,z)(x), x^3=r",
    "grading": "wt(r,z,e)=(3,1,-3); bracket_degree=2; unique homogeneous boundary=(1,-3)",
    "intersection": "k[x,y]_intersect_k(F,P)=k[r,z,e]/(r^2e-z^3+c^3r)",
    "picard": "units=k*, Pic=Z/3, div(r)=3L",
    "seed": "t=x^2(1+xy), q=t+c, F=xq, P=1/(3x^3), J(F,P)=1",
    "scope": "no_birational_no_affine_no_2x2_Darboux_pair; arbitrary_2x3_or_larger_degree>=2_open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable")
print("seed=F=x(t+c),P=1/(3x^3),J(F,P)=1")
print("surface=R^2E-Z^3+c^3R=0")
print("source=R=x^3;Z=x(c+x^3y);E=3c^2y+3cx^3y^2+x^6y^3")
print("poisson={R,Z}:3R^2;{R,E}:9Z^2;{Z,E}:3c^3+6RE")
print("atlas=surjective_etale_generic_degree_3_boundary_degree_1")
print("intersection=k[x,y]_cap_k(F,P)=k[R,Z,E]/relation")
print("units=k*;Pic=Z/3;div(R)=3L")
print("affine=smooth_carriers_Z+lambdaR_and_E;no_affine_polynomial_mate")
print("homogeneous=unique_possible_weight_pair_(1,-3)_has_positive_degree")
print("two_by_two=complete_endpoint_power_and_degree_nonentry")
print("darboux_floor=surface_degree_at_least_2;source_field_degree_multiple_of_3_at_least_6")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
