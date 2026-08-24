#!/usr/bin/env python3
"""Exact certificate for the two THM-3972 double-resultant seams."""

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


r, t, v, c = sp.symbols("r t v c", nonzero=True)


def z(label, expr):
    expr = sp.factor(sp.cancel(sp.expand(expr)))
    if expr != 0:
        raise AssertionError(f"{label}: {expr}")
    print(f"PASS  {label}")


def rows(c_value):
    return 1 - t * v**3, 3 * r + 3 * v + c_value * v**2


# The complete scalar discriminant identifies the two and only two seams.
D, N = rows(c)
Xi = sp.factor(sp.resultant(D, N, v))
z("general resultant", Xi - (c**3 + 27 * r**3 * t**2 + 27 * (1 - c * r) * t))
z(
    "resultant discriminant",
    sp.discriminant(Xi, t) + 27 * (c * r - 3) ** 2 * (4 * c * r - 3),
)

# Seam 1: cr=3.  The double resultant root comes from two distinct points.
c1 = 3 / r
D1, N1 = rows(c1)
t1 = 1 / r**3
z("cr=3 resultant square", sp.resultant(D1, N1, v) - 27 * (r**3 * t - 1) ** 2 / r**3)
z("cr=3 N factor", N1 - 3 * (v**2 + r * v + r**2) / r)
z(
    "cr=3 D factor on collision fibre",
    D1.subs(t, t1) + (v - r) * (v**2 + r * v + r**2) / r**3,
)
disc_two_points = sp.discriminant(v**2 + r * v + r**2, v)
z("cr=3 two-point discriminant", disc_two_points + 3 * r**2)
J1 = sp.det(
    sp.Matrix(
        [[sp.diff(D1, t), sp.diff(D1, v)], [sp.diff(N1, t), sp.diff(N1, v)]]
    )
)
Q1 = v**2 + r * v + r**2
res_J1 = sp.factor(sp.resultant(Q1, J1.subs(t, t1), v))
z("cr=3 transverse Jacobians", res_J1 - 27 * r**6)

# Seam 2: 4cr=3.  One length-two curvilinear base scheme appears.
c2 = sp.Rational(3, 4) / r
t2 = -sp.Rational(1, 8) / r**3
v2 = -2 * r
T, V = sp.symbols("T V")
D2, N2 = rows(c2)
z(
    "4cr=3 resultant square",
    sp.resultant(D2, N2, v) - 27 * (8 * r**3 * t + 1) ** 2 / (64 * r**3),
)
Dloc = sp.expand(D2.subs({t: t2 + T, v: v2 + V}))
Nloc = sp.expand(N2.subs({t: t2 + T, v: v2 + V}))
z("4cr=3 length-two N", Nloc - 3 * V**2 / (4 * r))
Dlin = (
    sp.diff(Dloc, T).subs({T: 0, V: 0}) * T
    + sp.diff(Dloc, V).subs({T: 0, V: 0}) * V
)
z("4cr=3 linear D", Dlin - (8 * r**3 * T + sp.Rational(3, 2) * V / r))
z(
    "4cr=3 etale (D,V) coordinate determinant",
    sp.diff(Dloc, T).subs({T: 0, V: 0}) - 8 * r**3,
)
res_local = sp.factor(sp.resultant(Dloc, Nloc, V))
z("4cr=3 resultant has order two", res_local - 27 * T**2 * r**3)

# The Rees chart of (L,V^2) is V^2=LZ, the A1 rational double point.
L, Z = sp.symbols("L Z")
rees = V**2 - L * Z
singular_gradient = [
    sp.diff(rees, w).subs({L: 0, V: 0, Z: 0}) for w in (L, V, Z)
]
for idx, val in enumerate(singular_gradient):
    z(f"A1 gradient component {idx}", val)
hessian = sp.hessian(rees, (L, V, Z))
z("A1 Hessian determinant nonzero", sp.det(hessian) + 2)
snf = smith_normal_form(sp.Matrix([[-2]]))
z("A1 resolution Smith determinant", abs(int(snf.det())) - 2)

print("SMITH  minimal-resolution intersection matrix [-2] -> cokernel Z/2")
print("RESULT cr=3: two reduced transverse basepoints in one fibre")
print("RESULT 4cr=3: one curvilinear (L,V^2) point and an A1 graph singularity")
print("ALL EXACT CHECKS PASSED")
