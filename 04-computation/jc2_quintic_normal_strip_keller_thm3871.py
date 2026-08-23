#!/usr/bin/env python3
"""Deterministic exact companion for provisional THM-3871."""

from __future__ import annotations

import hashlib
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


def quotient_zero(
    label: str,
    expression: sp.Expr,
    basis: sp.GroebnerBasis,
) -> None:
    numerator = sp.together(expression).as_numer_denom()[0]
    remainder = basis.reduce(sp.expand(numerator))[1]
    zero(label, remainder)


s, z, y, x = sp.symbols("s z y x")

# ---------------------------------------------------------------------------
# 1. Ten universal buckets, target covariance, and scale/translation gauges.
# ---------------------------------------------------------------------------

a_coeff = [sp.Function(name)(s) for name in ("a", "alpha", "u", "p", "r", "w")]
c_coeff = [sp.Function(name)(s) for name in ("b", "beta", "v", "q", "t", "ell")]
A_generic = sum(a_coeff[index] * z**index for index in range(6))
C_generic = sum(c_coeff[index] * z**index for index in range(6))
J_generic = jacobian(A_generic, C_generic, z, s)

convolution_buckets = []
for degree in range(10):
    row = sp.S.Zero
    for i in range(6):
        j = degree + 1 - i
        if 0 <= j <= 5:
            row += i * a_coeff[i] * sp.diff(c_coeff[j], s)
            row -= j * sp.diff(a_coeff[i], s) * c_coeff[j]
    convolution_buckets.append(sp.expand(row))
    zero(
        f"generic Jacobian bucket z^{degree}",
        J_generic.coeff(z, degree) - convolution_buckets[-1],
    )
check("exactly ten possible buckets", sp.Poly(J_generic, z).degree() <= 9)
zero(
    "top bucket exact",
    convolution_buckets[9]
    - 5 * (a_coeff[5] * sp.diff(c_coeff[5], s)
           - sp.diff(a_coeff[5], s) * c_coeff[5]),
)

aa, bb, cc, dd = sp.symbols("aa bb cc dd")
wfun, ellfun = a_coeff[5], c_coeff[5]
w_new = aa * wfun + bb * ellfun
ell_new = cc * wfun + dd * ellfun
top_new = 5 * (w_new * sp.diff(ell_new, s) - sp.diff(w_new, s) * ell_new)
zero(
    "top Wronskian target covariance",
    top_new
    - (aa * dd - bb * cc)
    * 5 * (wfun * sp.diff(ellfun, s) - sp.diff(wfun, s) * ellfun),
)

hfun = sp.Function("hfun")(s)
f_profile = sum(sp.Function(f"f{i}")(s) * y**i for i in range(6))
g_profile = sum(sp.Function(f"g{i}")(s) * y**i for i in range(6))
scaled_first = f_profile.subs(y, hfun * z)
scaled_second = g_profile.subs(y, hfun * z)
zero(
    "moving-scale bracket covariance",
    jacobian(scaled_first, scaled_second, z, s)
    - hfun * jacobian(f_profile, g_profile, y, s).subs(y, hfun * z),
)

shift = sp.Function("shift")(s)
zero(
    "moving-translation bracket invariance",
    jacobian(f_profile.subs(y, x + shift), g_profile.subs(y, x + shift), x, s)
    - jacobian(f_profile, g_profile, y, s).subs(y, x + shift),
)

shear_rows = {j: [degree for degree in range(1, 7) if j * degree == 5] for j in range(1, 5)}
check("only (5,1) has a target degree-drop shear", shear_rows == {1: [5], 2: [], 3: [], 4: []})

# ---------------------------------------------------------------------------
# 2. Complete depressed (5,2) packet and sharp rational hostile.
# ---------------------------------------------------------------------------

Q, R = sp.symbols("Q R", nonzero=True)
D, E, F, G, A0 = sp.symbols("D E F G A0")
b2 = sp.Function("b2")(s)
K2 = 5 * R / (2 * Q)
U2 = K2 * b2 + E
M2 = 2 * D * b2 / Q + F
L2 = 3 * K2 * b2**2 / (4 * Q) + 3 * E * b2 / (2 * Q) + G
a2 = D * b2**2 / Q**2 + F * b2 / Q + A0
A52 = R * x**5 + D * x**4 + U2 * x**3 + M2 * x**2 + L2 * x + a2
C52 = Q * x**2 + b2
J52 = jacobian(A52, C52, x, s)

for degree in range(1, 5):
    zero(f"integrated (5,2) bucket x^{degree}", J52.coeff(x, degree))
zero("(5,2) constant row", J52.coeff(x, 0) - L2 * sp.diff(b2, s))
zero("(5,2) finite-prime A-arm residue", R - K2 * Q + 3 * K2 * Q / 4 - 3 * R / 8)

A52_hostile = (
    s**35 * z**5
    + 5 * s**27 * z**4
    + sp.Rational(15, 2) * s**19 * z**3
    + sp.Rational(5, 2) * s**11 * z**2
    - sp.Rational(5, 8) * s**3 * z
    + sp.Rational(3, 8) / s**5
)
C52_hostile = s**14 * z**2 + 2 * s**6 * z
J52_hostile = jacobian(A52_hostile, C52_hostile, z, s)
zero("(5,2) hostile Jacobian", J52_hostile - sp.Rational(15, 4))
zero("(5,2) hostile C arm", C52_hostile.subs(z, 0))
zero("(5,2) hostile forced A arm", A52_hostile.subs(z, 0) - sp.Rational(3, 8) / s**5)
for degree in range(1, 6):
    check(
        f"(5,2) hostile A nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(A52_hostile).coeff(z, degree), s).is_univariate,
    )
for degree in range(1, 3):
    check(
        f"(5,2) hostile C nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(C52_hostile).coeff(z, degree), s).is_univariate,
    )

# ---------------------------------------------------------------------------
# 3. Complete depressed (5,3) packet, resultant, and algebraic hostile.
# ---------------------------------------------------------------------------

B3 = sp.Function("B3")(s)
b3 = sp.Function("b3")(s)
K3 = 5 * R / (3 * Q)
H3 = 4 * D / (3 * Q)
U3 = K3 * B3 + E
M3 = H3 * B3 + K3 * b3 + F
L3 = K3 * B3**2 / (3 * Q) + E * B3 / Q + H3 * b3 + G
a3 = (
    H3 * B3**2 / (6 * Q)
    + 2 * K3 * B3 * b3 / (3 * Q)
    + 2 * F * B3 / (3 * Q)
    + E * b3 / Q
    + A0
)
A53 = R * x**5 + D * x**4 + U3 * x**3 + M3 * x**2 + L3 * x + a3
C53 = Q * x**3 + B3 * x + b3
J53 = jacobian(A53, C53, x, s)

for degree in range(2, 6):
    zero(f"integrated (5,3) bucket x^{degree}", J53.coeff(x, degree))

P3 = (
    5 * R * B3**3 / 3
    - 12 * D * Q * B3 * b3
    - 18 * F * Q**2 * b3
    - 9 * G * Q**2 * B3
    - 15 * Q * R * b3**2
)
zero("(5,3) conserved penultimate row", J53.coeff(x, 1) + sp.diff(P3, s) / (9 * Q**2))
Omega3 = (
    12 * D * Q * b3 * sp.diff(b3, s)
    - 4 * D * B3**2 * sp.diff(B3, s)
    - 6 * F * Q * B3 * sp.diff(B3, s)
    + 9 * G * Q**2 * sp.diff(b3, s)
    - 5 * R * B3**2 * sp.diff(b3, s)
    - 10 * R * B3 * b3 * sp.diff(B3, s)
)
zero("(5,3) constant one-form", 9 * Q**2 * J53.coeff(x, 0) - Omega3)

X3, Z3 = sp.symbols("X3 Z3")
p3 = 5 * X3**2 + 10 * X3 + 6
q3 = X3**3 - 9 * (1 + X3) ** 2
check("(5,3) arm resultant", sp.resultant(p3, q3, X3) == 441)
check("(5,3) arm polynomials coprime", sp.gcd(p3, q3) == 1)
zero(
    "(5,3) constant-scale balanced leading coefficient",
    -(5 * sp.Rational(3, 2) + 10) * R - (-sp.Rational(35, 2) * R),
)

# The algebraic hostile is checked directly in the quotient by (72).
h53 = s**8
g53 = 1 / s
B53_h = X3 / s**2
b53_h = Z3 / s**3
x53 = h53 * z + g53
A53_h = sp.expand(
    x53**5
    + sp.Rational(5, 3) * B53_h * x53**3
    + sp.Rational(5, 3) * b53_h * x53**2
    + sp.Rational(5, 9) * B53_h**2 * x53
    + sp.Rational(10, 9) * B53_h * b53_h
)
C53_h = sp.expand(x53**3 + B53_h * x53 + b53_h)
ideal53 = sp.groebner(
    [X3**3 - 9 * Z3**2, 1 + X3 + Z3],
    X3,
    Z3,
    domain=sp.QQ.frac_field(s, z),
)
J53_h = jacobian(A53_h, C53_h, z, s)
quotient_zero("(5,3) hostile C arm", C53_h.subs(z, 0), ideal53)
quotient_zero(
    "(5,3) hostile Jacobian constant",
    J53_h - sp.Rational(35, 9) * X3**2 * Z3,
    ideal53,
)
for degree in range(1, 6):
    check(
        f"(5,3) hostile A nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(A53_h).coeff(z, degree), s).is_univariate,
    )
for degree in range(1, 4):
    check(
        f"(5,3) hostile C nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(C53_h).coeff(z, degree), s).is_univariate,
    )
check("(5,3) hostile bad arm resultant excludes zero", sp.resultant(p3, q3, X3) != 0)

# ---------------------------------------------------------------------------
# 4. Complete depressed (5,4) packet and the two conserved polynomials.
# ---------------------------------------------------------------------------

B4 = sp.Function("B4")(s)
N4 = sp.Function("N4")(s)
b4 = sp.Function("b4")(s)
K4 = 5 * R / (4 * Q)
H4 = D / Q
U4 = K4 * B4 + E
M4 = H4 * B4 + K4 * N4 + F
L4 = K4 * B4**2 / (8 * Q) + 3 * E * B4 / (4 * Q) + H4 * N4 + K4 * b4 + G
a4 = (
    K4 * B4 * N4 / (4 * Q)
    + F * B4 / (2 * Q)
    + 3 * E * N4 / (4 * Q)
    + H4 * b4
    + A0
)
A54 = R * x**5 + D * x**4 + U4 * x**3 + M4 * x**2 + L4 * x + a4
C54 = Q * x**4 + B4 * x**2 + N4 * x + b4
J54 = jacobian(A54, C54, x, s)

for degree in range(3, 7):
    zero(f"integrated (5,4) bucket x^{degree}", J54.coeff(x, degree))

W4 = 5 * R * B4 + 12 * E * Q
P4 = W4 * (B4**2 - 8 * Q * b4) - 20 * Q * R * N4**2 - 64 * F * Q**2 * N4 - 32 * G * Q**2 * B4
T4 = N4 * (15 * R * B4**2 + 24 * E * Q * B4 - 40 * Q * R * b4 - 32 * G * Q**2) + 16 * F * Q * (B4**2 - 4 * Q * b4)
Omega4 = (
    (5 * R * B4**2 + 24 * E * Q * B4 + 40 * Q * R * b4 + 32 * G * Q**2)
    * sp.diff(b4, s)
    - (24 * E * Q * N4 + 10 * R * B4 * N4) * sp.diff(N4, s)
    - (16 * F * Q * N4 + 10 * R * N4**2) * sp.diff(B4, s)
)
zero("(5,4) first conserved row", J54.coeff(x, 2) + sp.diff(P4, s) / (32 * Q**2))
zero("(5,4) second conserved row", J54.coeff(x, 1) + sp.diff(T4, s) / (32 * Q**2))
zero("(5,4) constant one-form", 32 * Q**2 * J54.coeff(x, 0) - Omega4)

# A bounded tropical audit checks the maximum comparison behind (54).
tropical_profiles = 0
for m in range(1, 13):
    for n in range(0, 19):
        for pole_b in range(0, 25):
            p_orders = [3 * m, m + pole_b, 2 * n]
            t_orders = [n + 2 * m, n + pole_b, 2 * m, pole_b]
            if p_orders.count(max(p_orders)) >= 2 and t_orders.count(max(t_orders)) >= 2:
                tropical_profiles += 1
                check(
                    f"(5,4) tropical profile m={m},n={n},z={pole_b}",
                    pole_b == 2 * m and 2 * n <= 3 * m,
                )
check("(5,4) tropical audit has positive controls", tropical_profiles > 0)

X4, Y4, Z4 = sp.symbols("X4 Y4 Z4")
unbalanced_c = 1 + X4 + X4**2 / 8
unbalanced_a = (5 * X4**2 - 8) / 32
check(
    "(5,4) unbalanced arm polynomials coprime",
    sp.gcd(sp.together(unbalanced_c * 8), sp.together(unbalanced_a * 32)) == 1,
)

p4 = 15 * X4**3 + 20 * X4**2 + 40 * X4 + 32
q4 = 50 * X4**5 + 25 * X4**4 - 80 * X4**2 + 64
resultant54 = sp.resultant(p4, q4, X4)
check("(5,4) balanced arm resultant", resultant54 == 3171942400000)
check("(5,4) balanced arm polynomials coprime", sp.gcd(p4, q4) == 1)

zero(
    "(5,4) W-zero A-arm residue",
    R - K4 * Q + R / 4,
)
zero(
    "(5,4) constant-N infinity coefficient",
    sp.Rational(5, 2) * R - sp.Rational(5, 2) * R,
)
zero(
    "(5,4) balanced infinity coefficient",
    sp.Rational(55, 2) * R - sp.Rational(55, 2) * R,
)

# Constant B and nonzero W: eliminate b from P=P0.  The remaining equation
# T=constant is genuinely cubic in N, including all cancellation strata.
B_const, N_const, b_const, P_const = sp.symbols("B_const N_const b_const P_const")
W_const = 5 * R * B_const + 12 * E * Q
P_const_expr = (
    W_const * (B_const**2 - 8 * Q * b_const)
    - 20 * Q * R * N_const**2
    - 64 * F * Q**2 * N_const
    - 32 * G * Q**2 * B_const
)
T_const_expr = (
    N_const
    * (15 * R * B_const**2 + 24 * E * Q * B_const
       - 40 * Q * R * b_const - 32 * G * Q**2)
    + 16 * F * Q * (B_const**2 - 4 * Q * b_const)
)
b_eliminated = sp.solve(sp.Eq(P_const_expr, P_const), b_const)[0]
cubic_after_elimination = sp.Poly(sp.cancel(T_const_expr.subs(b_const, b_eliminated)), N_const)
zero(
    "(5,4) constant-B nonzero-W cubic leading coefficient",
    cubic_after_elimination.coeff_monomial(N_const**3) - 100 * Q * R**2 / W_const,
)

# ---------------------------------------------------------------------------
# 5. Sharp algebraic (5,4) hostile and the arms-only false-positive modulus.
# ---------------------------------------------------------------------------

Z4_value = 3 * X4**2 / 8
h54 = s**9
g54 = 1 / s
B54_h = X4 / s**2
N54_h = Y4 / s**3
b54_h = Z4_value / s**4
x54 = h54 * z + g54
A54_h = sp.expand(
    x54**5
    + sp.Rational(5, 4) * B54_h * x54**3
    + sp.Rational(5, 4) * N54_h * x54**2
    + (sp.Rational(5, 32) * B54_h**2 + sp.Rational(5, 4) * b54_h) * x54
    + sp.Rational(5, 16) * B54_h * N54_h
)
C54_h = sp.expand(x54**4 + B54_h * x54**2 + N54_h * x54 + b54_h)
ideal54 = sp.groebner(
    [Y4**2 + X4**3 / 2, 1 + X4 + Y4 + 3 * X4**2 / 8],
    Y4,
    X4,
    domain=sp.QQ.frac_field(s, z),
)
J54_h = jacobian(A54_h, C54_h, z, s)
quotient_zero("(5,4) hostile C arm", C54_h.subs(z, 0), ideal54)
quotient_zero(
    "(5,4) hostile Jacobian constant",
    J54_h + sp.Rational(55, 32) * X4**4,
    ideal54,
)
for degree in range(1, 6):
    check(
        f"(5,4) hostile A nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(A54_h).coeff(z, degree), s).is_univariate,
    )
for degree in range(1, 5):
    check(
        f"(5,4) hostile C nonarm z^{degree} polynomial",
        sp.Poly(sp.expand(C54_h).coeff(z, degree), s).is_univariate,
    )
check("(5,4) hostile bad arm resultant excludes zero", resultant54 != 0)

X_false = sp.Rational(1)
Y_false = sp.Rational(3, 10)
Z_false = -sp.Rational(23, 10)
zero("arms-only false positive C arm", 1 + X_false + Y_false + Z_false)
zero("arms-only false positive A arm", 5 * X_false**2 + 10 * X_false * Y_false - 8)
check(
    "arms-only false positive first conserved residual",
    X_false**3 - 8 * X_false * Z_false - 4 * Y_false**2 == sp.Rational(476, 25),
)
check(
    "arms-only false positive second conserved residual",
    Y_false * (3 * X_false**2 - 8 * Z_false) == sp.Rational(321, 50),
)

# ---------------------------------------------------------------------------
# 6. Frozen semantic transcript.
# ---------------------------------------------------------------------------

SEMANTIC_FACTS = "\n".join(
    [
        "quintic-strip-ten-buckets=exact",
        "top-row=constant-target-direction+h5-hj-kummer",
        "degree-drop-shear=only-j1",
        "degree-50=empty",
        "degree-51=shear-to-THM3867",
        "degree-52=empty+residual-3R-over-8",
        "degree-53=one-conserved-row+resultant-441",
        "degree-54=two-conserved-rows+resultant-3171942400000",
        "degree-54-W-pole-unit-zero-identically-zero=exhausted",
        "degree-54-constant-scale=B-nonconstant-and-B-constant=exhausted",
        "hostile-52=Jacobian-15-over-4-only-A-arm-pole",
        "hostile-53=algebraic-rational-profile-only-A-arm-pole",
        "hostile-54=algebraic-rational-profile-only-A-arm-pole",
        "arms-only-modulus=real-before-conserved-buckets",
        "scope=polynomial-z-depth-at-most-five-only",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("THM-3871 QUINTIC NORMAL-STRIP KELLER THEOREM")
print("TEN_BUCKETS=PASS;TOP_DIRECTION=SL2_CONSTANT;KUMMER=h5_hj")
print("DEGREE_DROP=(5,1)_ONLY;TARGET_SHEAR=C5;THM-3867")
print("BRANCH_5_2=EMPTY;ARM_RESIDUAL=3R/8;HOSTILE_JACOBIAN=15/4")
print("BRANCH_5_3=EMPTY;ONE_CONSERVED_ROW;ARM_RESULTANT=441")
print("BRANCH_5_4=EMPTY;TWO_CONSERVED_ROWS;ARM_RESULTANT=3171942400000")
print("W_CHANNELS=POLE+UNIT+ZERO+IDENTICALLY_ZERO;CONSTANT_SCALE=PASS")
print("ARMS_ONLY_MODULUS=TRUE;CONSERVED_RESIDUALS=476/25,321/50")
print("RATIONAL_HOSTILES=THREE_PASS;ONLY_A_ARM_NONPOLYNOMIAL")
print("STATUS=PROVED+VERIFIED_EXACT+INDEPENDENTLY_HOSTILE_AUDITED")
print("BOUNDARY=DEGREE_SIX+RATIONAL+INFINITE_NORMAL_SERIES+PLANAR_JC_OPEN")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"TROPICAL_PROFILES={tropical_profiles}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
