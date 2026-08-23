#!/usr/bin/env python3
"""Exact scout for the first composite sextic normal-strip frontier.

This is a scoped algebraic packet audit, not a proof of the sextic strip and
not a result about arbitrary planar Keller maps.  It verifies:

* all twelve universal Jacobian buckets for z-depth at most six;
* the gcd/top-factor split and the target-shear rows j=1,2,3;
* the complete depressed rational packets in the residual rows (6,4),(6,5);
* the nonreduced length-six cusp obstruction in the balanced (6,4) face; and
* the reduced fourteen-point apparent scheme in the balanced (6,5) face,
  killed exactly by the next x^1 Jacobian bucket.
"""

from __future__ import annotations

import hashlib
import math
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def check(label: str, condition: object) -> None:
    """Optimization-safe exact gate."""

    global CHECKS
    CHECKS += 1
    if condition is True or condition == sp.S.true:
        return
    if isinstance(condition, sp.Basic) and sp.cancel(sp.factor(condition)) == 0:
        return
    raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def zero(label: str, expression: sp.Expr) -> None:
    check(label, sp.cancel(sp.factor(expression)) == 0)


def jacobian(first: sp.Expr, second: sp.Expr, x: sp.Symbol, base: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, base)
        - sp.diff(first, base) * sp.diff(second, x)
    )


s, z, y, x = sp.symbols("s z y x")

# ---------------------------------------------------------------------------
# 1. Twelve universal buckets, target direction, and the gcd row split.
# ---------------------------------------------------------------------------

a_coeff = [sp.Function(f"a{i}")(s) for i in range(7)]
c_coeff = [sp.Function(f"c{i}")(s) for i in range(7)]
A_generic = sum(a_coeff[i] * z**i for i in range(7))
C_generic = sum(c_coeff[i] * z**i for i in range(7))
J_generic = jacobian(A_generic, C_generic, z, s)

buckets = []
for degree in range(12):
    row = sp.S.Zero
    for i in range(7):
        j = degree + 1 - i
        if 0 <= j <= 6:
            row += i * a_coeff[i] * sp.diff(c_coeff[j], s)
            row -= j * sp.diff(a_coeff[i], s) * c_coeff[j]
    row = sp.expand(row)
    buckets.append(row)
    zero(
        f"universal Jacobian bucket z^{degree}",
        J_generic.coeff(z, degree) - row,
    )

check("exactly twelve universal buckets", len(buckets) == 12)
check("generic Jacobian has z-degree at most eleven", sp.Poly(J_generic, z).degree() <= 11)
zero(
    "top Wronskian bucket",
    buckets[11]
    - 6
    * (
        a_coeff[6] * sp.diff(c_coeff[6], s)
        - sp.diff(a_coeff[6], s) * c_coeff[6]
    ),
)

aa, bb, cc, dd = sp.symbols("aa bb cc dd")
wfun, qfun = a_coeff[6], c_coeff[6]
w_new = aa * wfun + bb * qfun
q_new = cc * wfun + dd * qfun
top_new = 6 * (w_new * sp.diff(q_new, s) - sp.diff(w_new, s) * q_new)
zero(
    "constant target covariance of the top bucket",
    top_new
    - (aa * dd - bb * cc)
    * 6
    * (wfun * sp.diff(qfun, s) - sp.diff(wfun, s) * qfun),
)

gcd_rows = {
    j: (math.gcd(6, j), 6 // math.gcd(6, j), j // math.gcd(6, j))
    for j in range(1, 6)
}
check(
    "gcd-reduced top packets",
    gcd_rows
    == {
        1: (1, 6, 1),
        2: (2, 3, 1),
        3: (3, 2, 1),
        4: (2, 3, 2),
        5: (1, 6, 5),
    },
)
shear_degrees = {j: 6 // j if 6 % j == 0 else None for j in range(1, 6)}
check(
    "exact target-shear rows",
    shear_degrees == {1: 6, 2: 3, 3: 2, 4: None, 5: None},
)

hfun = sp.Function("hfun")(s)
f_profile = sum(sp.Function(f"f{i}")(s) * y**i for i in range(7))
g_profile = sum(sp.Function(f"g{i}")(s) * y**i for i in range(7))
zero(
    "moving common-scale bracket covariance",
    jacobian(f_profile.subs(y, hfun * z), g_profile.subs(y, hfun * z), z, s)
    - hfun * jacobian(f_profile, g_profile, y, s).subs(y, hfun * z),
)
shift = sp.Function("shift")(s)
zero(
    "moving translation bracket invariance",
    jacobian(f_profile.subs(y, x + shift), g_profile.subs(y, x + shift), x, s)
    - jacobian(f_profile, g_profile, y, s).subs(y, x + shift),
)

# ---------------------------------------------------------------------------
# 2. Complete depressed rational (6,4) packet.
# ---------------------------------------------------------------------------

Q, R = sp.symbols("Q R", nonzero=True)
D, E, F, G, I, A0 = sp.symbols("D E F G I A0")
B4 = sp.Function("B4")(s)
N4 = sp.Function("N4")(s)
b4 = sp.Function("b4")(s)

K4 = 3 * R / (2 * Q)
H4 = 5 * D / (4 * Q)
U4 = K4 * B4 + E
M4 = K4 * N4 + H4 * B4 + F
L4 = K4 * b4 + H4 * N4 + K4 * B4**2 / (4 * Q) + E * B4 / Q + G
P4 = (
    5 * D * b4 / (4 * Q)
    + K4 * B4 * N4 / (2 * Q)
    + E * N4 / Q
    + H4 * B4**2 / (8 * Q)
    + 3 * F * B4 / (4 * Q)
    + I
)
a4 = (
    -R * B4**3
    + 5 * D * Q * B4 * N4
    + 8 * G * Q**2 * B4
    + 12 * Q * R * B4 * b4
    + 16 * E * Q**2 * b4
    + 12 * F * Q**2 * N4
    + 6 * Q * R * N4**2
) / (16 * Q**3) + A0

A64 = R * x**6 + D * x**5 + U4 * x**4 + M4 * x**3 + L4 * x**2 + P4 * x + a4
C64 = Q * x**4 + B4 * x**2 + N4 * x + b4
J64 = jacobian(A64, C64, x, s)

for degree in range(3, 9):
    zero(f"integrated (6,4) bucket x^{degree}", J64.coeff(x, degree))

K64_2 = (
    -5 * B4**3 * D
    - 12 * B4**2 * F * Q
    - 24 * B4**2 * N4 * R
    + 40 * B4 * D * Q * b4
    + 32 * B4 * I * Q**2
    + 20 * D * N4**2 * Q
    + 96 * F * Q**2 * b4
    + 64 * G * N4 * Q**2
    + 96 * N4 * Q * R * b4
) / (32 * Q**2)
K64_1 = (
    3 * B4**4 * R
    - 15 * B4**2 * D * N4 * Q
    - 16 * B4**2 * G * Q**2
    - 24 * B4**2 * Q * R * b4
    - 24 * B4 * F * N4 * Q**2
    - 24 * B4 * N4**2 * Q * R
    + 40 * D * N4 * Q**2 * b4
    + 64 * G * Q**3 * b4
    + 32 * I * N4 * Q**3
    + 48 * Q**2 * R * b4**2
) / (32 * Q**3)

zero("(6,4) first conserved bucket", J64.coeff(x, 2) - sp.diff(K64_2, s))
zero("(6,4) second conserved bucket", J64.coeff(x, 1) - sp.diff(K64_1, s))
zero("(6,4) constant one-form", J64.coeff(x, 0) - (P4 * sp.diff(b4, s) - sp.diff(a4, s) * N4))

# Principal balanced pole coordinates.  The C arm eliminates Z.
X4, Y4, Z4 = sp.symbols("X4 Y4 Z4")
C_arm64 = 1 + X4 + Y4 + Z4
conserved64_a = Y4 * (4 * Z4 - X4**2)
conserved64_b = (X4**2 - 4 * Z4) ** 2 - 8 * X4 * Y4**2
A_arm64 = 6 * Y4**2 - (X4 + 2) ** 3
bracket64 = Y4 * (X4**3 - 4 * X4 * Z4 - 6 * Y4**2)

u4, v4 = sp.symbols("u4 v4")
sub64 = {X4: u4 - 2, Y4: v4, Z4: 1 - u4 - v4}
local64 = sp.groebner(
    [
        sp.expand(A_arm64.subs(sub64)),
        sp.expand(conserved64_a.subs(sub64)),
        sp.expand(conserved64_b.subs(sub64)),
    ],
    v4,
    u4,
    order="lex",
)
check("(6,4) local Groebner basis has three rows", len(local64.polys) == 3)
for expected in (
    6 * v4**2 - u4**3,
    u4**2 * (2 * u4 + 3 * v4),
    u4**4,
):
    zero("(6,4) expected local relation", local64.reduce(expected)[1])

leading_exponents64 = {
    poly.LM(order=local64.order).exponents for poly in local64.polys
}
check(
    "(6,4) local initial ideal",
    leading_exponents64 == {(2, 0), (1, 2), (0, 4)},
)
standard64 = [
    (v_exp, u_exp)
    for v_exp in range(2)
    for u_exp in range(4)
    if not (v_exp >= 1 and u_exp >= 2)
]
check("(6,4) nonreduced local length is six", len(standard64) == 6)
zero(
    "(6,4) leading constant bracket vanishes on the local scheme",
    local64.reduce(sp.expand(bracket64.subs(sub64)))[1],
)
bracket64_local = sp.expand(bracket64.subs(sub64))
relation64_1 = 6 * v4**2 - u4**3
relation64_2 = u4**2 * (2 * u4 + 3 * v4)
relation64_3 = u4**4
zero(
    "(6,4) explicit leading-bracket ideal certificate",
    bracket64_local
    - (
        (-v4 + 2 * u4 / 3 - sp.Rational(4, 3)) * relation64_1
        - sp.Rational(2, 3) * relation64_2
        + sp.Rational(2, 3) * relation64_3
    ),
)

# Arms alone have a large false-positive cusp; this rational point is killed
# by both conserved rows.
arms64 = {X4: sp.Rational(4), Y4: sp.Rational(6), Z4: sp.Rational(-11)}
zero("(6,4) arms-only hostile C arm", C_arm64.subs(arms64))
zero("(6,4) arms-only hostile A arm", A_arm64.subs(arms64))
check("(6,4) arms-only first conserved residual", conserved64_a.subs(arms64) == -360)
check("(6,4) arms-only second conserved residual", conserved64_b.subs(arms64) == 2448)

# The unique reduced support is the exact cusp A=T^3,C=T^2.  It pays every
# zero bucket but cannot pay a nonzero Keller constant.
T64 = s**2 * z**2 + 2 * z
A64_cusp = T64**3
C64_cusp = T64**2
zero("(6,4) exact cusp hostile has zero bracket", jacobian(A64_cusp, C64_cusp, z, s))
check("(6,4) exact cusp hostile z-degrees", (sp.degree(A64_cusp, z), sp.degree(C64_cusp, z)) == (6, 4))

# Primitive simple-pole strict transform.  Let t=g^{-1}; a regular C arm may
# begin at t^4.  These first forced coefficients are independent of that arm
# and of the even constants E,G,A0.  They identify the honest next chart;
# they do not exhaust higher pole order or the later strict-transform rows.
tjet = sp.symbols("tjet")
u1, u2, u3 = sp.symbols("u1 u2 u3")
v1, v2, v3 = sp.symbols("v1 v2 v3")
c4j, c5j, c6j = sp.symbols("c4j c5j c6j")
Xj = -2 + u1 * tjet + u2 * tjet**2 + u3 * tjet**3
Yj = v1 * tjet + v2 * tjet**2 + v3 * tjet**3
Zj = (
    1
    - (Xj + 2)
    - Yj
    + c4j * tjet**4
    + c5j * tjet**5
    + c6j * tjet**6
)

Aface64 = sp.expand(
    (
        -2 * Xj**3
        + 5 * Xj**2 * D * tjet
        + 12 * Xj**2
        + 10 * Xj * Yj * D * tjet
        + 24 * Xj * Yj
        + 24 * Xj * Zj
        + 40 * Xj * D * tjet
        + 32 * Xj * E * tjet**2
        + 24 * Xj * F * tjet**3
        + 16 * Xj * G * tjet**4
        + 48 * Xj
        + 12 * Yj**2
        + 40 * Yj * D * tjet
        + 32 * Yj * E * tjet**2
        + 24 * Yj * F * tjet**3
        + 48 * Yj
        + 40 * Zj * D * tjet
        + 32 * Zj * E * tjet**2
        + 48 * Zj
        + 32 * A0 * tjet**6
        + 32 * D * tjet
        + 32 * E * tjet**2
        + 32 * F * tjet**3
        + 32 * G * tjet**4
        + 32 * I * tjet**5
        + 32
    )
    / 32
)
K2face64 = sp.expand(
    (
        -5 * Xj**3 * D * tjet
        - 24 * Xj**2 * Yj
        - 12 * Xj**2 * F * tjet**3
        + 40 * Xj * Zj * D * tjet
        + 32 * Xj * I * tjet**5
        + 20 * Yj**2 * D * tjet
        + 96 * Yj * Zj
        + 64 * Yj * G * tjet**4
        + 96 * Zj * F * tjet**3
    )
    / 32
)

zero("(6,4) primitive jet t1 A coefficient", Aface64.coeff(tjet, 1) - 3 * D / 8)
zero("(6,4) primitive jet t1 K2 coefficient", K2face64.coeff(tjet, 1) + 5 * D / 4)

jet_after_D = {D: 0}
zero(
    "(6,4) primitive jet t2 forces v1",
    Aface64.coeff(tjet, 2).subs(jet_after_D) - 3 * v1**2 / 8,
)

jet_after_v1 = {D: 0, v1: 0}
zero(
    "(6,4) primitive jet t3 K2 forces F",
    K2face64.coeff(tjet, 3).subs(jet_after_v1) - 3 * F / 2,
)
zero(
    "(6,4) primitive jet t3 A then forces u1",
    Aface64.coeff(tjet, 3).subs({D: 0, v1: 0, F: 0}) + u1**3 / 16,
)

jet_after_u1 = {D: 0, v1: 0, F: 0, u1: 0}
zero(
    "(6,4) primitive jet t4 K2 forces v2",
    K2face64.coeff(tjet, 4).subs(jet_after_u1) + 3 * v2**2,
)

jet_after_v2 = {D: 0, v1: 0, F: 0, u1: 0, v2: 0}
zero(
    "(6,4) primitive jet t5 A forces I",
    Aface64.coeff(tjet, 5).subs(jet_after_v2) - I,
)

jet_after_I = {D: 0, v1: 0, F: 0, u1: 0, v2: 0, I: 0}
zero(
    "(6,4) primitive jet t6 K2 forces v3",
    K2face64.coeff(tjet, 6).subs(jet_after_I) + 3 * v3**2,
)

# ---------------------------------------------------------------------------
# 3. Complete depressed rational (6,5) packet.
# ---------------------------------------------------------------------------

H = sp.symbols("H")
B5 = sp.Function("B5")(s)
N5 = sp.Function("N5")(s)
V5 = sp.Function("V5")(s)
b5 = sp.Function("b5")(s)

K5 = 6 * R / (5 * Q)
U5 = K5 * B5 + E
M5 = K5 * N5 + D * B5 / Q + F
L5 = K5 * V5 + D * N5 / Q + K5 * B5**2 / (10 * Q) + 4 * E * B5 / (5 * Q) + G
P5 = (
    K5 * b5
    + D * V5 / Q
    + K5 * B5 * N5 / (5 * Q)
    + 3 * F * B5 / (5 * Q)
    + 4 * E * N5 / (5 * Q)
    + H
)
a5 = (
    -4 * B5**3 * R
    - 10 * B5**2 * E * Q
    + 50 * B5 * G * Q**2
    + 30 * B5 * Q * R * V5
    + 125 * D * Q**2 * b5
    + 100 * E * Q**2 * V5
    + 75 * F * N5 * Q**2
    + 15 * N5**2 * Q * R
) / (125 * Q**3) + A0

A65 = R * x**6 + D * x**5 + U5 * x**4 + M5 * x**3 + L5 * x**2 + P5 * x + a5
C65 = Q * x**5 + B5 * x**3 + N5 * x**2 + V5 * x + b5
J65 = jacobian(A65, C65, x, s)

for degree in range(4, 10):
    zero(f"integrated (6,5) bucket x^{degree}", J65.coeff(x, degree))

K65_3 = (
    -15 * B5**2 * F * Q
    - 12 * B5**2 * N5 * R
    - 20 * B5 * E * N5 * Q
    + 25 * B5 * H * Q**2
    + 30 * B5 * Q * R * b5
    + 100 * E * Q**2 * b5
    + 75 * F * Q**2 * V5
    + 50 * G * N5 * Q**2
    + 30 * N5 * Q * R * V5
) / (25 * Q**2)
K65_2 = (
    9 * B5**4 * R
    + 20 * B5**3 * E * Q
    - 75 * B5**2 * G * Q**2
    - 60 * B5**2 * Q * R * V5
    - 100 * B5 * E * Q**2 * V5
    - 150 * B5 * F * N5 * Q**2
    - 60 * B5 * N5**2 * Q * R
    - 50 * E * N5**2 * Q**2
    + 375 * F * Q**3 * b5
    + 250 * G * Q**3 * V5
    + 125 * H * N5 * Q**3
    + 150 * N5 * Q**2 * R * b5
    + 75 * Q**2 * R * V5**2
) / (125 * Q**3)

zero("(6,5) first conserved bucket", J65.coeff(x, 3) - sp.diff(K65_3, s))
zero("(6,5) second conserved bucket", J65.coeff(x, 2) - sp.diff(K65_2, s))
zero(
    "(6,5) nonexact penultimate one-form",
    J65.coeff(x, 1)
    - (
        2 * L5 * sp.diff(b5, s)
        - 2 * N5 * sp.diff(a5, s)
        + P5 * sp.diff(V5, s)
        - V5 * sp.diff(P5, s)
    ),
)
zero(
    "(6,5) constant one-form",
    J65.coeff(x, 0) - (P5 * sp.diff(b5, s) - V5 * sp.diff(a5, s)),
)

# The x^1 row is genuinely nonexact before imposing the conserved equations.
Bv, Nv, Vv, bv = sp.symbols("Bv Nv Vv bv")
subs_vars65 = {B5: Bv, N5: Nv, V5: Vv, b5: bv}
coeff_B = sp.expand(J65.coeff(x, 1)).coeff(sp.diff(B5, s)).subs(subs_vars65)
coeff_N = sp.expand(J65.coeff(x, 1)).coeff(sp.diff(N5, s)).subs(subs_vars65)
curl_BN = sp.factor(sp.diff(coeff_B, Nv) - sp.diff(coeff_N, Bv))
check("(6,5) x1 one-form has nonzero exterior derivative", curl_BN != 0)

# Principal balanced pole coordinates and exact apparent scheme.
X5, Y5, Z5, W5 = sp.symbols("X5 Y5 Z5 W5")
C_arm65 = 1 + X5 + Y5 + Z5 + W5
conserved65_3 = -2 * X5**2 * Y5 + 5 * X5 * W5 + 5 * Y5 * Z5
conserved65_2 = 3 * X5**4 - 20 * X5**2 * Z5 - 20 * X5 * Y5**2 + 50 * Y5 * W5 + 25 * Z5**2
row65_1 = 25 * W5 * X5**2 + 225 * W5 * Z5 + 8 * X5**3 * Y5 - 65 * X5 * Y5 * Z5 - 30 * Y5**3
A_arm65 = -4 * X5**3 + 15 * X5**2 + 30 * X5 * Y5 + 30 * X5 * Z5 + 15 * Y5**2 - 25
bracket65 = 125 * W5**2 + 25 * W5 * X5 * Y5 + 4 * X5**3 * Z5 - 30 * X5 * Z5**2 - 15 * Y5**2 * Z5

apparent65 = [C_arm65, conserved65_3, conserved65_2, A_arm65]
scheme65 = sp.groebner(apparent65, W5, Z5, Y5, X5, order="lex")
check("(6,5) apparent scheme has triangular four-row basis", len(scheme65.polys) == 4)
terminal65 = sp.Poly(scheme65.polys[-1].as_expr(), X5)
check("(6,5) apparent scheme has length fourteen", terminal65.degree() == 14)
check("(6,5) terminal eliminant is squarefree", sp.gcd(terminal65, terminal65.diff()).degree() == 0)
check("(6,5) apparent scheme is nonempty", terminal65.LC() != 0)
terminal65_coefficients = [int(value) for value in terminal65.all_coeffs()]
check("(6,5) terminal eliminant leading coefficient", terminal65_coefficients[0] == 23548)
check("(6,5) terminal eliminant constant coefficient", terminal65_coefficients[-1] == 15625000)
terminal65_sha256 = hashlib.sha256(
    ",".join(str(value) for value in terminal65_coefficients).encode("ascii")
).hexdigest()

closed65 = sp.groebner(apparent65 + [row65_1], W5, Z5, Y5, X5, order="grevlex")
check("(6,5) x1 row cuts the apparent scheme to empty", len(closed65.polys) == 1 and closed65.polys[0].as_expr() == 1)
nonzero_bracket65 = sp.groebner(apparent65 + [bracket65], W5, Z5, Y5, X5, order="grevlex")
check(
    "(6,5) constant leading bracket is nonzero at all fourteen apparent points",
    len(nonzero_bracket65.polys) == 1 and nonzero_bracket65.polys[0].as_expr() == 1,
)

# Arms alone admit a rational false positive; the conserved and x1 rows see it.
arms65 = {
    X5: sp.Rational(1),
    Y5: sp.Rational(0),
    Z5: sp.Rational(7, 15),
    W5: -sp.Rational(37, 15),
}
zero("(6,5) arms-only hostile C arm", C_arm65.subs(arms65))
zero("(6,5) arms-only hostile A arm", A_arm65.subs(arms65))
check("(6,5) arms-only first conserved residual", conserved65_3.subs(arms65) == -sp.Rational(37, 3))
check("(6,5) arms-only second conserved residual", conserved65_2.subs(arms65) == -sp.Rational(8, 9))
check("(6,5) arms-only x1 residual", row65_1.subs(arms65) == -sp.Rational(962, 3))
check("(6,5) arms-only constant residual", bracket65.subs(arms65) == sp.Rational(6803, 9))

# ---------------------------------------------------------------------------
# 4. Frozen semantic transcript.
# ---------------------------------------------------------------------------

SEMANTIC_FACTS = "\n".join(
    [
        "sextic-strip-twelve-buckets=exact",
        "top-direction=constant-target-SL2",
        "gcd-packets=j1-6:1,j2-3:1,j3-2:1,j4-3:2,j5-6:5",
        "target-shear-rows=j1+j2+j3-to-THM3871",
        "residual-rows=(6,4)+(6,5)",
        "degree-64=two-conserved-polynomials+constant-one-form",
        "degree-64-balanced-support=(-2,0,1)",
        "degree-64-local-scheme=nonreduced-length-6",
        "degree-64-leading-bracket=zero-mod-local-ideal",
        "degree-64-primitive-jet=D-F-I-zero+u-order2+Y-order4",
        "degree-64-next-task=strict-transform-jet+quadratic-descent-sidecar",
        "degree-65=two-conserved-polynomials+nonexact-x1-row+constant-one-form",
        "degree-65-arms-plus-conserved=14-reduced-points",
        "degree-65-x1-row=unit-ideal-closure",
        "degree-65-leading-bracket=nonzero-on-all-14-apparent-points",
        "scope=principal-balanced-faces+exact-rational-packets-not-full-sextic-theorem",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("SEXTIC NORMAL-STRIP COMPOSITE FRONTIER SCOUT")
print("TWELVE_BUCKETS=PASS;TOP_DIRECTION=SL2_CONSTANT")
print("GCD_PACKETS=j1:6/1,j2:3/1,j3:2/1,j4:3/2,j5:6/5")
print("TARGET_SHEARS=j1:C6,j2:C3,j3:C2;REDUCE_TO_THM-3871")
print("RESIDUAL_ROWS=(6,4)+(6,5)")
print("ROW_6_4=PACKET_EXACT;TWO_CONSERVED;LOCAL_SUPPORT=(-2,0,1)")
print("ROW_6_4_LOCAL_SCHEME=NONREDUCED_LENGTH_6;LEADING_BRACKET=ZERO_MOD_IDEAL")
print("ROW_6_4_PRIMITIVE_JET=D=F=I=0;u=O(t^2);Y=O(t^4)")
print("ROW_6_4_ARMS_ONLY_HOSTILE=4,6,-11;RESIDUALS=-360,2448")
print("ROW_6_5=PACKET_EXACT;TWO_CONSERVED+NONEXACT_X1+CONSTANT")
print("ROW_6_5_APPARENT_SCHEME=14_REDUCED_POINTS;X1_BUCKET=UNIT_IDEAL")
print("ROW_6_5_LEADING_BRACKET=NONZERO_ON_ALL_14_APPARENT_POINTS")
print(f"ROW_6_5_TERMINAL_ELIMINANT_SHA256={terminal65_sha256}")
print("ROW_6_5_ARMS_ONLY_HOSTILE=1,0,7/15,-37/15")
print("STATUS=FINITE-EXACT+SCOUT;NOT_A_SEXTIC_STRIP_THEOREM;JC2_OPEN")
print("NEXT=STRICT_TRANSFORM_(6,4)_CUSP_JETS+QUADRATIC_DESCENT_PARITY")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
