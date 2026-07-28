#!/usr/bin/env python3
"""Exact referee for THM-2699.

Classify affine source planes and linear target projections of the
THM-2696 S3<S4 quotient whose induced planar Jacobian is a nonzero constant.
The symbolic calculation is supplemented by exhaustive raw finite-field
censuses in characteristics 5 and 7.
"""

from __future__ import annotations

import itertools

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, y, z = sp.symbols("x y z")
n1, n2, n3 = sp.symbols("n1 n2 n3")
l1, l2, l3 = sp.symbols("l1 l2 l3")
delta = sp.symbols("delta")

A = x**2 - 2 * y
B = y**2 - 2 * x * z
d = z
DF = sp.Matrix([A, B, d]).jacobian([x, y, z])
adj = DF.adjugate()
n = sp.Matrix([n1, n2, n3])
ell = sp.Matrix([l1, l2, l3])
K = sp.expand((n.T * adj * ell)[0])

expected_adj = sp.Matrix(
    [
        [2 * y, 2, 4 * x],
        [2 * z, 2 * x, 4 * x**2],
        [0, 0, 4 * x * y - 4 * z],
    ]
)
expected_K = sp.expand(
    4 * l3 * n2 * x**2
    + 4 * l3 * n3 * x * y
    + (2 * l2 * n2 + 4 * l3 * n1) * x
    + 2 * l1 * n1 * y
    + (2 * l1 * n2 - 4 * l3 * n3) * z
    + 2 * l2 * n1
)
require(adj == expected_adj, "adjugate mismatch")
require(sp.expand(K - expected_K) == 0, "universal Jacobian mismatch")

# n3 != 0: eliminate z and expose the four nonconstant coefficients.
K_n3 = sp.Poly(
    sp.together(K.subs(z, (delta - n1 * x - n2 * y) / n3)), x, y
)
coeff_n3 = {
    "x2": sp.factor(K_n3.coeff_monomial(x**2)),
    "xy": sp.factor(K_n3.coeff_monomial(x * y)),
    "x": sp.factor(K_n3.coeff_monomial(x)),
    "y": sp.factor(K_n3.coeff_monomial(y)),
    "one": sp.factor(K_n3.coeff_monomial(1)),
}
require(coeff_n3["x2"] == 4 * l3 * n2, "n3 x2 coefficient")
require(coeff_n3["xy"] == 4 * l3 * n3, "n3 xy coefficient")
require(
    sp.factor(coeff_n3["x"].subs(l3, 0) - 2 * n2 * (l2 * n3 - l1 * n1) / n3)
    == 0,
    "n3 x coefficient",
)
require(
    sp.factor(coeff_n3["y"].subs(l3, 0) - 2 * l1 * (n1 * n3 - n2**2) / n3)
    == 0,
    "n3 y coefficient",
)
require(
    sp.factor(
        coeff_n3["one"].subs(l3, 0)
        - 2 * (l2 * n1 * n3 + l1 * n2 * delta) / n3
    )
    == 0,
    "n3 constant coefficient",
)

# n3=0,n2!=0: eliminate y.  A nonzero constant would force ell=0.
K_n2 = sp.Poly(
    sp.together(K.subs({n3: 0, y: (delta - n1 * x) / n2})), x, z
)
require(K_n2.coeff_monomial(z) == 2 * l1 * n2, "n2 z coefficient")
require(K_n2.coeff_monomial(x**2) == 4 * l3 * n2, "n2 x2 coefficient")
require(
    sp.factor(K_n2.coeff_monomial(x).subs({l1: 0, l3: 0}) - 2 * l2 * n2)
    == 0,
    "n2 x coefficient",
)

# Explicit triangular inverses for all three surviving families.
a, b = sp.symbols("a b", nonzero=True)
L2, L3 = sp.symbols("L2 L3")

# I: n=(n1,0,n3), n1*n3!=0, ell=e2; use target (A,d).
x_I = (delta - n3 * d) / n1
y_I = (x_I**2 - A) / 2
require(sp.factor(n1 * x_I + n3 * d - delta) == 0, "family I plane inverse")
require(sp.factor(x_I**2 - 2 * y_I - A) == 0, "family I A inverse")

# II: x=x0, ell=(0,L2,L3); use (A,C=L3*B-L2*d).
x0, AA, CC, dd = sp.symbols("x0 AA CC dd")
y_II = (x0**2 - AA) / 2
z_II = (L3 * y_II**2 - CC) / (2 * L3 * x0 + L2)
require(
    sp.factor(L3 * (y_II**2 - 2 * x0 * z_II) - L2 * z_II - CC) == 0,
    "family II inverse",
)

# III: n=(a^2,ab,b^2), ell=(b^2,a^2,0); use (d,C=a^2*A-b^2*B).
CC3, dd3 = sp.symbols("CC3 dd3")
y_III = (delta**2 - b**4 * dd3**2 - a**2 * CC3) / (2 * a * (a**3 + b * delta))
x_III = (delta - a * b * y_III - b**2 * dd3) / a**2
require(
    sp.factor(a**2 * x_III + a * b * y_III + b**2 * dd3 - delta) == 0,
    "family III plane inverse",
)
require(
    sp.factor(
        a**2 * (x_III**2 - 2 * y_III)
        - b**2 * (y_III**2 - 2 * x_III * dd3)
        - CC3
    )
    == 0,
    "family III target inverse",
)

# Direct planar Jacobians in the displayed triangular coordinates.
zI, yI = sp.symbols("zI yI")
xI = (delta - n3 * zI) / n1
JI = sp.factor(sp.det(sp.Matrix([xI**2 - 2 * yI, zI]).jacobian([zI, yI])))
require(JI == 2, "family I planar Jacobian")

yII, zII = sp.symbols("yII zII")
AII = x0**2 - 2 * yII
CII = L3 * (yII**2 - 2 * x0 * zII) - L2 * zII
JII = sp.factor(sp.det(sp.Matrix([AII, CII]).jacobian([yII, zII])))
require(sp.expand(JII - 2 * (L2 + 2 * L3 * x0)) == 0, "family II planar Jacobian")

yIII, zIII = sp.symbols("yIII zIII")
xIII = (delta - a * b * yIII - b**2 * zIII) / a**2
CIII = sp.expand(a**2 * (xIII**2 - 2 * yIII) - b**2 * (yIII**2 - 2 * xIII * zIII))
JIII = sp.factor(sp.det(sp.Matrix([zIII, CIII]).jacobian([yIII, zIII])))
require(
    sp.factor(JIII - 2 * (a**3 + b * delta) / a) == 0,
    "family III planar Jacobian",
)


def actual_and_family(p: int, nn: tuple[int, int, int], ll: tuple[int, int, int], de: int):
    """Return (actual_nonzero_constant, family_tag) over F_p."""
    a1, a2, a3 = nn
    q1, q2, q3 = ll
    if a3 % p:
        inv = pow(a3, -1, p)
        conditions = (
            (q3 * a2) % p == 0,
            (q3 * a3) % p == 0,
            (a2 * (q2 * a3 - q1 * a1)) % p == 0,
            (q1 * (a1 * a3 - a2 * a2)) % p == 0,
        )
        value = (2 * (q2 * a1 * a3 + q1 * a2 * de) * inv) % p
        actual = all(conditions) and value != 0
    elif a2 % p:
        # z, x^2, then x force q1=q3=q2=0.
        actual = False
    else:
        value = (2 * q2 * a1 + 4 * q3 * de) % p
        actual = q1 % p == 0 and value != 0

    tag = None
    if a3 % p and a2 % p == 0 and a1 % p and q1 % p == 0 and q3 % p == 0 and q2 % p:
        tag = "I"
    elif a3 % p == 0 and a2 % p == 0 and a1 % p and q1 % p == 0:
        if (2 * q2 * a1 + 4 * q3 * de) % p:
            tag = "II"
    elif a3 % p and a2 % p:
        if (
            (a1 * a3 - a2 * a2) % p == 0
            and q3 % p == 0
            and q1 % p
            and (q2 * a3 - q1 * a1) % p == 0
        ):
            inv = pow(a3, -1, p)
            value = (2 * (q2 * a1 * a3 + q1 * a2 * de) * inv) % p
            if value:
                tag = "III"
    return actual, tag


def finite_census(p: int) -> tuple[int, dict[str, int]]:
    vectors = [v for v in itertools.product(range(p), repeat=3) if v != (0, 0, 0)]
    counts = {"I": 0, "II": 0, "III": 0}
    actual_count = 0
    for nn in vectors:
        for ll in vectors:
            for de in range(p):
                actual, tag = actual_and_family(p, nn, ll, de)
                require(actual == (tag is not None), f"F_{p} classification mismatch: {nn},{ll},{de}")
                if actual:
                    actual_count += 1
                    counts[tag] += 1
    return actual_count, counts


censuses = {p: finite_census(p) for p in (5, 7)}

print("THM-2699 affine-plane/linear-projection Keller slice classification")
print(f"adjugate={adj.tolist()}")
print(f"universal_K={K}")
print("families=I:n2=0,n1*n3!=0,ell~e2;II:n~e1,ell1=0;III:n~(a^2,ab,b^2),ell~(b^2,a^2,0)")
print(f"triangular_jacobians=I:{JI};II:{JII};III:{JIII}")
for p, (total, counts) in censuses.items():
    print(f"F_{p}:raw_nonzero_constant={total}:I={counts['I']}:II={counts['II']}:III={counts['III']}")
print("zero_boundaries=II:L2+2*L3*x0;III:a^3+b*delta")
print("nonlinear_constant_different_slice=outside_affine_plane_scope")
print("ALL CHECKS PASSED")
