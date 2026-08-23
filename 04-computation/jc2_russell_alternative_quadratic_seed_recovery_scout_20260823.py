#!/usr/bin/env python3
"""Finite-exact scout for alternative quadratic-normal Russell seeds.

The script normalizes the nodal controls by x=Z/(3c^3), u=3c^3 S and
studies every seed

    A=u^2-x+p(u)x^2,
    C=u^3-u-(3/2)ux+q(u)x^2

with affine p,q.  It classifies monic recovery of the hidden u-coordinate,
freezes the unique affine-coefficient failure, records an exact valuative
hostile, and separates that failure from an earlier x-only recovery escape.
It also checks the recovery-derivative/seed-Jacobian identity and removes
quadratic jets which are merely rational source-coordinate gauges.

This is a FINITE-EXACT/SCOUT artifact.  It constructs no Keller composite
and claims no consequence for JC(2).
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    """Optimization-safe exact gate."""

    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"GATE FAILED: {label}: {condition}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def jacobian(
    left: sp.Expr,
    right: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
) -> sp.Expr:
    return sp.expand(
        sp.diff(left, first) * sp.diff(right, second)
        - sp.diff(left, second) * sp.diff(right, first)
    )


def laurent_order(expression: sp.Expr, variable: sp.Symbol) -> int:
    """Exact order at variable=0 for a nonzero rational expression."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, variable)
    denominator_poly = sp.Poly(denominator, variable)
    gate(not numerator_poly.is_zero, "Laurent numerator nonzero")
    gate(not denominator_poly.is_zero, "Laurent denominator nonzero")
    numerator_order = min(term[0][0] for term in numerator_poly.terms())
    denominator_order = min(term[0][0] for term in denominator_poly.terms())
    return numerator_order - denominator_order


# ---------------------------------------------------------------------------
# 1. Universal quadratic-normal seed and recovery/different identities.
# ---------------------------------------------------------------------------
x, u, A, C = sp.symbols("x u A C")
p0, p1, q0, q1 = sp.symbols("p0 p1 q0 q1")
p = p1 * u + p0
q = q1 * u + q0

A_seed = u**2 - x + p * x**2
C_seed = u**3 - u - sp.Rational(3, 2) * u * x + q * x**2
D_seed = sp.factor(jacobian(A_seed, C_seed, x, u))

expected_D = (
    1
    + (
        sp.Rational(3, 2)
        + (6 * u**2 - 2) * p
        - 4 * u * q
    ) * x
    + (-sp.diff(q, u) - 3 * p + sp.Rational(3, 2) * u * sp.diff(p, u))
    * x**2
    + 2 * (p * sp.diff(q, u) - sp.diff(p, u) * q) * x**3
)
zero(D_seed - expected_D, "universal affine-pq seed Jacobian")

F = A_seed - A
G = C_seed - C
R_affine = sp.factor(sp.resultant(F, G, x))
L = sp.factor(sp.Rational(3, 2) * p * u - q)
M = sp.factor(q * (u**2 - A) - p * (u**3 - u - C))

# This compact formula is valid before any specialization which drops the
# quadratic x-degree.  It is a polynomial identity over the generic ring.
zero(
    p * R_affine - (p * M**2 + L * M + (u**2 - A) * L**2),
    "compact quadratic resultant",
)
zero(
    (
        sp.diff(R_affine, u) - L * D_seed
    ).subs({A: A_seed, C: C_seed}, simultaneous=True),
    "resultant derivative equals sheet factor times seed Jacobian",
)

# When p=0, x recovers polynomially and the direct hidden-u relation has
# derivative exactly minus the seed Jacobian.
q_linear = q1 * u + q0
X_linear = u**2 - A
H_linear = (
    u**3 - u - sp.Rational(3, 2) * u * X_linear
    + q_linear * X_linear**2 - C
)
A_pzero = u**2 - x
C_pzero = u**3 - u - sp.Rational(3, 2) * u * x + q_linear * x**2
D_pzero = jacobian(A_pzero, C_pzero, x, u)
zero(H_linear.subs({A: A_pzero, C: C_pzero}), "p=0 recovery relation")
zero(
    (
        sp.diff(H_linear, u) + D_pzero
    ).subs({A: A_pzero, C: C_pzero}, simultaneous=True),
    "p=0 recovery derivative is minus seed Jacobian",
)


# ---------------------------------------------------------------------------
# 2. Full affine p,q monicity classification.
# ---------------------------------------------------------------------------
# A nonzero linear coefficient of p dominates every other infinity balance.
gate(sp.Poly(R_affine, u).degree() == 8, "affine p nonconstant degree eight")
zero(sp.Poly(R_affine, u).LC() - p1**2,
     "affine p nonconstant scalar leading coefficient")

a, d, e = sp.symbols("a d e", nonzero=True)
p_constant = a
q_affine = d * u + e
A_constant = u**2 - x + p_constant * x**2
C_constant = (
    u**3 - u - sp.Rational(3, 2) * u * x + q_affine * x**2
)
R_constant = sp.factor(
    sp.resultant(A_constant - A, C_constant - C, x)
)

# Generic slope d differs from the A-correction a.
gate(sp.Poly(R_constant, u).degree() == 6,
     "constant p generic affine q degree six")
zero(sp.Poly(R_constant, u).LC() - (a - d) ** 2,
     "constant p generic affine q scalar leading coefficient")

# The first infinity resonance d=a drops the recovery degree to four.
R_resonant = sp.factor(R_constant.subs(d, a))
gate(sp.Poly(R_resonant, u).degree() == 4,
     "first resonance recovery degree four")
zero(
    sp.Poly(R_resonant, u).LC() - (a + 4 * e**2) / 4,
    "first resonance scalar leading coefficient",
)

# The unique second resonance in this cell destroys monicity over k[A,C].
R_exceptional = sp.factor(R_resonant.subs(a, -4 * e**2))
gate(sp.Poly(R_exceptional, u).degree() == 3,
     "exceptional hidden-u recovery degree three")
exceptional_lc = e * (16 * e**2 * (A - 1) - 1) / 2
zero(sp.Poly(R_exceptional, u).LC() - exceptional_lc,
     "exceptional target-dependent leading coefficient")

# The p=0 degree-drop chart must be computed directly, rather than by
# specializing a fixed-degree Sylvester determinant.
H_pzero_qlinear = sp.expand(H_linear)
gate(sp.Poly(H_pzero_qlinear, u).degree() == 5,
     "p=0 q-linear recovery degree five")
zero(sp.Poly(H_pzero_qlinear, u).LC() - q1,
     "p=0 q-linear scalar leading coefficient")
H_pzero_qconstant = sp.expand(H_pzero_qlinear.subs(q1, 0))
gate(sp.Poly(H_pzero_qconstant, u).degree() == 4,
     "p=0 q-constant recovery degree four")
zero(sp.Poly(H_pzero_qconstant, u).LC() - q0,
     "p=0 q-constant scalar leading coefficient")
H_original = sp.expand(H_pzero_qconstant.subs(q0, 0))
gate(sp.Poly(H_original, u).degree() == 3,
     "fixed seed recovery cubic")
zero(sp.Poly(H_original, u).LC() + sp.Rational(1, 2),
     "fixed seed recovery cubic scalar leading coefficient")


# ---------------------------------------------------------------------------
# 3. Unique affine-coefficient hidden-u failure and exact omitted divisor.
# ---------------------------------------------------------------------------
p_exceptional = -4 * e**2
q_exceptional = -4 * e**2 * u + e
A_exceptional = u**2 - x + p_exceptional * x**2
C_exceptional = (
    u**3 - u - sp.Rational(3, 2) * u * x + q_exceptional * x**2
)
D_exceptional = sp.factor(jacobian(A_exceptional, C_exceptional, x, u))
L_exceptional = -e * (1 + 2 * e * u)
zero(
    (
        sp.diff(R_exceptional, u) - L_exceptional * D_exceptional
    ).subs({A: A_exceptional, C: C_exceptional}, simultaneous=True),
    "exceptional cubic derivative mismatch factor",
)

t, v = sp.symbols("t v")
u_tv = 1 / t
x_tv = 1 / (2 * e * t) - 1 / (8 * e**2) - t / (4 * e) + v * t**2
A_star = 1 + 1 / (16 * e**2)
C_star = (1 - 8 * e**2) / (64 * e**3)
A_tv_expected = (
    A_star - 4 * e * t * v - t**2 / 4 + 2 * e * t**3 * v
    - 4 * e**2 * t**4 * v**2
)
C_tv_expected = (
    C_star - 4 * e * v
    + t * (v / 2 - sp.Rational(1, 4) + 1 / (16 * e**2))
    + t**2 * (2 * e * v - v / (4 * e) + 1 / (16 * e))
    + t**3 * (-4 * e**2 * v**2 - v / 2)
    + e * t**4 * v**2
)
A_tv = sp.factor(A_exceptional.subs({u: u_tv, x: x_tv}))
C_tv = sp.factor(C_exceptional.subs({u: u_tv, x: x_tv}))
zero(A_tv - A_tv_expected, "exceptional exact regular A hostile")
zero(C_tv - C_tv_expected, "exceptional exact regular C hostile")
gate(sp.cancel(A_tv).is_polynomial(t, v), "exceptional hostile A polynomial")
gate(sp.cancel(C_tv).is_polynomial(t, v), "exceptional hostile C polynomial")
gate(laurent_order(u_tv, t) == -1, "exceptional hidden u has simple pole")
gate(laurent_order(x_tv, t) == -1, "exceptional hidden x has simple pole")
zero(jacobian(x_tv, u_tv, t, v) - 1,
     "exceptional pole chart is symplectic")
D_tv = sp.factor(D_exceptional.subs({u: u_tv, x: x_tv}))
zero(D_tv - jacobian(A_tv_expected, C_tv_expected, t, v),
     "exceptional chart seed/target Jacobians agree")
zero(D_tv.subs(t, 0) - 16 * e**2 * v,
     "exceptional omitted-divisor Jacobian is generically nonzero")
zero(A_tv_expected.subs(t, 0) - A_star,
     "exceptional nonproper set is a vertical line")
zero(C_tv_expected.subs(t, 0) - (C_star - 4 * e * v),
     "exceptional omitted divisor maps affinely onto line")
zero(
    (16 * e**2 * (A_tv_expected - 1) - 1).subs(t, 0),
    "exceptional nonmonic leading divisor matches nonproper line",
)


# ---------------------------------------------------------------------------
# 4. Earlier x-only recovery failure with u still monic/integral.
# ---------------------------------------------------------------------------
alpha, u0 = sp.symbols("alpha u0", nonzero=True)
g = u - u0
p_xescape = alpha * g
q_xescape = sp.Rational(3, 2) * alpha * u0 * g
A_xescape = u**2 - x + p_xescape * x**2
C_xescape = (
    u**3 - u - sp.Rational(3, 2) * u * x + q_xescape * x**2
)
D_xescape = sp.factor(jacobian(A_xescape, C_xescape, x, u))
w = sp.factor(alpha * g * x)
X_sidecar = sp.factor(u**2 - A_xescape)
f_recovery = sp.expand(
    u**3 + (2 - 3 * A) * u + 2 * C
)
K_recovery = -alpha * f_recovery / 3
H_xescape = sp.expand(
    (K_recovery + alpha * (u - u0) * (u**2 - A)) ** 2
    - K_recovery
)
zero(alpha * g * X_sidecar - w * (1 - w),
     "x-only escape Danielewski sidecar relation")
zero(
    f_recovery.subs({A: A_xescape, C: C_xescape}, simultaneous=True)
    + 3 * w**2 / alpha,
    "x-only escape square recovery packet",
)
zero(
    (K_recovery + alpha * (u - u0) * (u**2 - A)).subs(
        {A: A_xescape, C: C_xescape}, simultaneous=True
    ) - w,
    "x-only escape recovers regular sidecar w",
)
zero(H_xescape.subs({A: A_xescape, C: C_xescape}, simultaneous=True),
     "x-only escape monic hidden-u relation")
gate(sp.Poly(H_xescape, u).degree() == 6,
     "x-only escape primitive recovery degree six")
zero(sp.Poly(H_xescape, u).LC() - 4 * alpha**2 / 9,
     "x-only escape scalar leading coefficient")
zero(
    (
        sp.diff(H_xescape, u) - sp.Rational(2, 3) * alpha * D_xescape
    ).subs({A: A_xescape, C: C_xescape}, simultaneous=True),
    "x-only escape power-basis derivative equals seed Jacobian",
)

R_xescape = sp.factor(
    sp.resultant(A_xescape - A, C_xescape - C, x)
)
zero(
    R_xescape - sp.Rational(9, 4) * (u - u0) ** 2 * H_xescape,
    "x-only escape strips vertical Sylvester factor",
)

y = sp.symbols("y")
u_gy = u0 + t
x_gy = 1 / (alpha * t) + y
A_gy_expected = u_gy**2 + y + alpha * t * y**2
C_gy_expected = (
    u_gy**3 - u_gy - 3 / (2 * alpha)
    + sp.Rational(3, 2) * (2 * u0 - u_gy) * y
    + sp.Rational(3, 2) * alpha * u0 * t * y**2
)
A_gy = sp.factor(A_xescape.subs({u: u_gy, x: x_gy}))
C_gy = sp.factor(C_xescape.subs({u: u_gy, x: x_gy}))
zero(A_gy - A_gy_expected, "x-only exact regular A hostile")
zero(C_gy - C_gy_expected, "x-only exact regular C hostile")
gate(sp.cancel(A_gy).is_polynomial(t, y), "x-only hostile A polynomial")
gate(sp.cancel(C_gy).is_polynomial(t, y), "x-only hostile C polynomial")
gate(laurent_order(x_gy, t) == -1, "x-only hidden x has simple pole")
zero(jacobian(x_gy, u_gy, t, y) + 1,
     "x-only pole chart has constant Jacobian minus one")
D_gy = sp.factor(D_xescape.subs({u: u_gy, x: x_gy}))
zero(D_gy.subs(t, 0) + (3 * y + 2) / 2,
     "x-only pole branch Jacobian generically nonzero")

# On g=0, the sidecar relation has two branches.  X is finite x on w=0
# and equals -y on the pole branch w=1.
X_regular = sp.symbols("X_regular")
zero(
    D_xescape.subs({u: u0, x: X_regular}) - (3 * X_regular + 2) / 2,
    "x-only regular branch different residue",
)
zero(
    D_gy.subs(t, 0).subs(y, -X_regular) - (3 * X_regular - 2) / 2,
    "x-only pole branch different residue",
)

# If a Keller composite existed, h={u,w} would satisfy D*h=-alpha*g.
# Bracketing alpha*g*X=w(1-w) with u then gives
# D*{u,X}=2w-1.  The two boundary residues must therefore be units.
h_bracket, uX_bracket = sp.symbols("h_bracket uX_bracket")
zero(
    (D_xescape * h_bracket + alpha * g).subs(
        h_bracket, alpha * g * uX_bracket / (1 - 2 * w)
    ) * (1 - 2 * w) / (alpha * g)
    - (D_xescape * uX_bracket - (2 * w - 1)),
    "x-only Keller sidecar divisor identity",
)


# ---------------------------------------------------------------------------
# 5. Quadratic gauge warning and an exact rational cosmetic family.
# ---------------------------------------------------------------------------
eta_gauge = sp.factor(sp.Rational(3, 2) * p * u - q)
xi_gauge = sp.factor(2 * u * eta_gauge - p)
x_control = x + xi_gauge * x**2
u_control = u + eta_gauge * x**2
A_fixed_control = sp.expand(u_control**2 - x_control)
C_fixed_control = sp.expand(
    u_control**3 - u_control
    - sp.Rational(3, 2) * u_control * x_control
)
zero(
    A_fixed_control - A_seed - eta_gauge**2 * x**4,
    "unique quadratic gauge A residual begins in order four",
)
zero(
    C_fixed_control - C_seed
    - (
        -sp.Rational(3, 2) * eta_gauge * x**3
        + sp.Rational(3, 2) * eta_gauge * p * x**4
        + eta_gauge**3 * x**6
    ),
    "unique quadratic gauge first genuine residual is cubic",
)

# L=0 is exactly the cosmetic polynomial precomposition x -> x-px^2.
p_cosmetic = p0 + p1 * u
q_cosmetic = sp.Rational(3, 2) * u * p_cosmetic
x_cosmetic = x - p_cosmetic * x**2
A_cosmetic = u**2 - x + p_cosmetic * x**2
C_cosmetic = (
    u**3 - u - sp.Rational(3, 2) * u * x + q_cosmetic * x**2
)
zero(A_cosmetic - (u**2 - x_cosmetic),
     "L=0 cosmetic A is fixed-seed precomposition")
zero(
    C_cosmetic
    - (u**3 - u - sp.Rational(3, 2) * u * x_cosmetic),
    "L=0 cosmetic C is fixed-seed precomposition",
)

# A rational Mobius normal change has apparent nonintegral x-recovery but
# recovers chi=x/(1+lambda*x) by the original monic cubic, so it is not a
# new alternative seed.
lam = sp.symbols("lam", nonzero=True)
chi = x / (1 + lam * x)
A_mobius = u**2 - chi
C_mobius = u**3 - u - sp.Rational(3, 2) * u * chi
D_mobius = sp.factor(jacobian(A_mobius, C_mobius, x, u))
zero(A_mobius.subs(x, 0) - u**2, "Mobius gauge preserves arm A")
zero(C_mobius.subs(x, 0) - (u**3 - u),
     "Mobius gauge preserves arm C")
zero(sp.diff(A_mobius, x).subs(x, 0) + 1,
     "Mobius gauge preserves first normal A")
zero(sp.diff(C_mobius, x).subs(x, 0) + sp.Rational(3, 2) * u,
     "Mobius gauge preserves first normal C")
zero(
    D_mobius - (1 + sp.Rational(3, 2) * chi) / (1 + lam * x) ** 2,
    "Mobius seed Jacobian is fixed different times gauge Jacobian",
)
zero(
    (
        u**3 + (2 - 3 * A) * u + 2 * C
    ).subs({A: A_mobius, C: C_mobius}, simultaneous=True),
    "Mobius cosmetic family retains fixed monic recovery cubic",
)


# ---------------------------------------------------------------------------
# 6. Frozen semantic packet and optimization-safe source audit.
# ---------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "normalization": "x=Z/(3c^3),u=3c^3S;Jacobian-one scaling",
    "seed": "A=u2-x+p(u)x2;C=u3-u-3ux/2+q(u)x2",
    "affine_classification": (
        "p_linear=>u_monic;constant_p_affine_q has unique hidden-u "
        "exception p=q_slope=-4e2,q_intercept=e"
    ),
    "resultant_different": "R_u=L*D after vertical-factor stripping",
    "hidden_u_hostile": (
        "u=1/t;x=1/(2et)-1/(8e2)-t/(4e)+vt2;A,C polynomial"
    ),
    "hidden_u_nonproper_line": "16e2(A-1)-1=0",
    "x_only_escape": (
        "p=alpha(u-u0);q=3alpha*u0*(u-u0)/2;u,w integral;x=w/(alpha*g)"
    ),
    "x_only_boundary": "w=0,1 require (3X+2)/2,(3X-2)/2 units",
    "gauge": "all quadratic jets are gauge;first exact residual is cubic",
    "rational_warning": "Mobius normal changes factor through fixed seed",
    "status": "FINITE-EXACT/SCOUT;no Keller composite;JC2 open",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_RUSSELL_ALTERNATIVE_QUADRATIC_SEED_RECOVERY_SCOUT")
print("status=FINITE-EXACT+SCOUT;NO_KELLER_COMPOSITE;JC2_OPEN")
print("normalized_seed=A:u2-x+p*x2;C:u3-u-3ux/2+q*x2")
print("universal_different=R_u_on_seed=L*D;L=3pu/2-q")
print("affine_pq=p1_nonzero=>degree8_LC=p1^2")
print("constant_p=d!=a=>degree6_LC=(a-d)^2")
print("first_resonance=d=a=>degree4_LC=(a+4e^2)/4")
print("hidden_u_exception=a=d=-4e^2;degree3_LC=e*(16e^2(A-1)-1)/2")
print("hidden_u_hostile=u,x_simple_poles;A,C_polynomial;J_control=1")
print("hidden_u_nonproper_line=16e^2(A-1)-1=0;D_boundary=16e^2*v")
print("x_only_escape=u_integral;w=alpha(u-u0)x_integral;x_may_pole")
print("x_only_recovery=degree6_LC=4alpha^2/9;H_u=2alpha*D/3")
print("x_only_divisor_test=w0:(3X+2)/2_unit;w1:(3X-2)/2_unit")
print("gauge_warning=quadratic_jet_cosmetic;first_exact_residual=C_order3")
print("rational_warning=Mobius_normal_change_factors_through_fixed_seed")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
