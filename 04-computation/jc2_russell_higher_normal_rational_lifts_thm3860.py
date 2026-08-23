#!/usr/bin/env python3
"""Exact companion for THM-3860.

The script checks polynomial/formal identities only.  The all-degree pole
lemma is proved in the theorem file; here we verify its algebraic exceptional
family and all displayed constructive formulas without numeric sampling.
"""

from __future__ import annotations

import hashlib

import sympy as sp


CHECKS = 0


def check(label: str, condition: object) -> None:
    """Assertion that remains active under python -O."""

    global CHECKS
    CHECKS += 1
    if condition is True or condition == sp.S.true:
        return
    if isinstance(condition, sp.Basic) and sp.simplify(condition) == 0:
        return
    raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def zero(label: str, expression: sp.Expr) -> None:
    check(label, sp.cancel(sp.factor(expression)) == 0)


s, z = sp.symbols("s z")

# ---------------------------------------------------------------------------
# 1. Full coefficient convolution and affine recursion.
# ---------------------------------------------------------------------------

a = [sp.Function(f"a{i}")(s) for i in range(6)]
b = [sp.Function(f"b{i}")(s) for i in range(6)]
A_series = sum(a[i] * z**i for i in range(6))
C_series = sum(b[i] * z**i for i in range(6))
J_series = sp.expand(sp.diff(A_series, z) * sp.diff(C_series, s)
                     - sp.diff(A_series, s) * sp.diff(C_series, z))

for m in range(5):
    convolution = sum(
        i * a[i] * sp.diff(b[j], s) - j * sp.diff(a[i], s) * b[j]
        for i in range(6)
        for j in range(6)
        if i + j == m + 1
    )
    zero(f"coefficient convolution m={m}",
         sp.expand(J_series).coeff(z, m) - convolution)

p, q, x, y, rhs, denom, tau = sp.symbols(
    "p q x y rhs denom tau", nonzero=True
)
u = -rhs * x / denom + p * tau
v = -rhs * y / denom + q * tau
det_uv = sp.expand(u * q - p * v)
zero("affine recursion determinant factor",
     det_uv + rhs / denom * (x * q - p * y))
zero("tangent kernel", (p * tau) * q - p * (q * tau))

# A nontrivial nodal packet exercises the first two recursive rows.
c = sp.symbols("c", nonzero=True)
a0 = 9 * c**6 * s**2
b0 = 27 * c**9 * s**3 - 3 * c**3 * s
a1 = -sp.Rational(1, 3) / c**3
b1 = -sp.Rational(3, 2) * s
w_nodal = sp.Rational(1, 2) / c**3
zero("nodal Bezout row",
     a1 * sp.diff(b0, s) - sp.diff(a0, s) * b1 - 1)
zero("nodal constant Wronskian",
     a1 * sp.diff(b1, s) - sp.diff(a1, s) * b1 - w_nodal)

tau1 = s**2 + 2 * s + 3
a2 = -w_nodal * a1 / 2 + sp.diff(a0, s) * tau1
b2 = -w_nodal * b1 / 2 + sp.diff(b0, s) * tau1
row1 = (a1 * sp.diff(b1, s) - sp.diff(a1, s) * b1
        + 2 * (a2 * sp.diff(b0, s) - sp.diff(a0, s) * b2))
zero("first higher-normal row", row1)

R2 = (a1 * sp.diff(b2, s) - 2 * sp.diff(a1, s) * b2
      + 2 * a2 * sp.diff(b1, s) - sp.diff(a2, s) * b1)
tau2 = s**3 - s + 1
a3 = -R2 * a1 / 3 + sp.diff(a0, s) * tau2
b3 = -R2 * b1 / 3 + sp.diff(b0, s) * tau2
row2 = R2 + 3 * (a3 * sp.diff(b0, s) - sp.diff(a0, s) * b3)
zero("second higher-normal row", row2)

A3 = a0 + a1 * z + a2 * z**2 + a3 * z**3
C3 = b0 + b1 * z + b2 * z**2 + b3 * z**3
J3 = sp.expand(sp.diff(A3, z) * sp.diff(C3, s)
               - sp.diff(A3, s) * sp.diff(C3, z))
zero("truncated bracket constant", J3.coeff(z, 0) - 1)
zero("truncated bracket z", J3.coeff(z, 1))
zero("truncated bracket z2", J3.coeff(z, 2))

# ---------------------------------------------------------------------------
# 2. Complete vertical rational class and its first jets.
# ---------------------------------------------------------------------------

w = sp.symbols("w", nonzero=True)
phi_fun = sp.Function("phi")(z)
f_fun = sp.Function("f")(z)
L_fun = 1 / ((1 + w * phi_fun) * sp.diff(phi_fun, z))
S_fun = L_fun * s + f_fun
J_transform = (sp.diff(phi_fun, z) * sp.diff(S_fun, s)
               - sp.diff(phi_fun, s) * sp.diff(S_fun, z))
zero("arbitrary vertical density identity",
     (1 + w * phi_fun) * J_transform - 1)

phi2, f2, phi3 = sp.symbols("phi2 f2 phi3")
phi_jet = z - w * z**2 / 2 + phi3 * z**3
f_jet = f2 * z**2
L_jet = 1 / ((1 + w * phi_jet) * sp.diff(phi_jet, z))
check("L constant jet", sp.expand(sp.series(L_jet, z, 0, 2).removeO()).coeff(z, 0) == 1)
check("L linear jet vanishes", sp.expand(sp.series(L_jet, z, 0, 2).removeO()).coeff(z, 1) == 0)
check("S preserves first jet",
      sp.expand(sp.series(L_jet * s + f_jet, z, 0, 2).removeO()) == s)

phi_mobius = z / (1 + w * z / 2)
L_mobius = sp.factor(1 / ((1 + w * phi_mobius) * sp.diff(phi_mobius, z)))
L_expected = (1 + w * z / 2) ** 3 / (1 + 3 * w * z / 2)
zero("Mobius L closed form", L_mobius - L_expected)
mobius_Z_series = sp.series(phi_mobius, z, 0, 4).removeO().expand()
mobius_L_series = sp.series(L_mobius, z, 0, 4).removeO().expand()
check("Mobius forced z2", mobius_Z_series.coeff(z, 2) == -w / 2)
check("Mobius tangent z2", mobius_L_series.coeff(z, 2) == 3 * w**2 / 4)
check("Mobius tangent z3", mobius_L_series.coeff(z, 3) == -w**3)

# The first s-dependent normal term cannot occur in Z before order three.
u_fun = sp.Function("u")(s)
h_fun = sp.Function("h")(s)
Z_mix = z - w * z**2 / 2 + u_fun * z**3
S_mix = s + h_fun * z**2
density_mix = sp.series(
    (1 + w * Z_mix)
    * (sp.diff(Z_mix, z) * sp.diff(S_mix, s)
       - sp.diff(Z_mix, s) * sp.diff(S_mix, z)),
    z, 0, 3,
).removeO().expand()
zero("mixed order z coefficient", density_mix.coeff(z, 1))
zero("mixed order z2 law",
     density_mix.coeff(z, 2)
     - (3 * u_fun + sp.diff(h_fun, s) - 3 * w**2 / 2))

# ---------------------------------------------------------------------------
# 3. Explicit nodal rational lift and named prime pole.
# ---------------------------------------------------------------------------

h = 1 + z / (4 * c**3)
g = 1 + 3 * z / (4 * c**3)
Z_nodal = z / h
S_nodal = s * h**3 / g
A_nodal = 9 * c**6 * S_nodal**2 - Z_nodal / (3 * c**3)
C_nodal_raw = (27 * c**9 * S_nodal**3 - 3 * c**3 * S_nodal
               - sp.Rational(3, 2) * S_nodal * Z_nodal)
C_nodal_closed = 27 * c**9 * s**3 * h**9 / g**3 - 3 * c**3 * s * h**2
zero("nodal C simplification", C_nodal_raw - C_nodal_closed)
J_nodal = (sp.diff(A_nodal, z) * sp.diff(C_nodal_raw, s)
           - sp.diff(A_nodal, s) * sp.diff(C_nodal_raw, z))
zero("nodal exact rational Darboux identity", J_nodal - 1)
zero("nodal arm A", A_nodal.subs(z, 0) - a0)
zero("nodal arm C", C_nodal_raw.subs(z, 0) - b0)
zero("nodal first normal A", sp.diff(A_nodal, z).subs(z, 0) - a1)
zero("nodal first normal C", sp.diff(C_nodal_raw, z).subs(z, 0) - b1)

pole_residue = sp.simplify(sp.limit(h * A_nodal, z, -4 * c**3))
check("named h-prime simple pole residue", pole_residue == sp.Rational(4, 3))

r, z0 = sp.symbols("r z0", nonzero=True)
e_on_divisor = (z0**3 - c**3 * r) / r**2
zero("constant-z quotient relation",
     r**2 * e_on_divisor - z0**3 + c**3 * r)
D_on_divisor = c**3 + e_on_divisor * r
s_on_divisor = sp.factor(e_on_divisor / (3 * D_on_divisor))
zero("constant-z divisor s formula",
     s_on_divisor - (1 / (3 * r) - c**3 / (3 * z0**3)))
check("constant-z divisor s nonconstant", sp.diff(s_on_divisor, r) != 0)

# ---------------------------------------------------------------------------
# 4. Exceptional no-bad-finite-point family in the pole lemma.
# ---------------------------------------------------------------------------

d = sp.symbols("d", integer=True, positive=True)
aa = -w / d
phi_exceptional = ((1 + aa * z) ** (-d) - 1) / w
check("exceptional phi value", phi_exceptional.subs(z, 0) == 0)
zero("exceptional phi first derivative",
     sp.diff(phi_exceptional, z).subs(z, 0) - 1)
exceptional_second = sp.simplify(sp.diff(phi_exceptional, z, 2).subs(z, 0))
zero("exceptional phi second derivative formula",
     exceptional_second - w * (d + 1) / d)
zero("exceptional jet contradiction factor",
     exceptional_second + w - w * (2 * d + 1) / d)

SEMANTIC_FACTS = "\n".join(
    [
        "formal-recursion=affine-tangent-torsor",
        "first-higher-normal=-W/2-normal-plus-free-tangent",
        "constant-W-vertical-rational-family=exact",
        "nodal-canonical-square-gate=not-higher-normal-invariant",
        "explicit-nodal-rational-lift=has-z-plus-4c3-prime-pole",
        "all-Z-of-z-rational-lifts=vertical-pole-barrier",
        "next-live-coordinate=Z_s-nonzero-from-order-three",
        "scope=no-global-Darboux-pair-no-JC2-conclusion",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("THM-3860 RUSSELL HIGHER-NORMAL RATIONAL LIFTS")
print("FORMAL_RECURSION=PASS;COEFFICIENT_ORDERS=0..4;FREE_DIRECTION=ARM_TANGENT")
print("VERTICAL_RATIONAL_CLASS=PASS;JACOBIAN=1;ARM_AND_FIRST_JET=PRESERVED")
print("NODAL_CONTROL=PASS;FIELD=Frac(B);NAMED_POLE=z+4c^3;ORDER_A=-1")
print("VERTICAL_POLE_LEMMA=PASS;EXCEPTIONAL_ONE_ROOT_FAMILY=JET_CONTRADICTION")
print("BOUNDARY=Z_s_MUST_BE_NONZERO;FIRST_POSSIBLE_ORDER=3;GLOBAL_ENTRY=OPEN")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
