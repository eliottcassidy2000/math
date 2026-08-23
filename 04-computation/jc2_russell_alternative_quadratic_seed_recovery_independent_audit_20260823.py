#!/usr/bin/env python3
"""Independent finite-exact audit of the alternative quadratic seed scout.

This audit deliberately reconstructs the elimination polynomial as a 4-by-4
Sylvester determinant instead of calling ``sympy.resultant``.  It checks every
affine leading-coefficient stratum, both pole charts, the conditional boundary
bracket obstruction, and the unique quadratic normal-coordinate jet.

Scope is FINITE-EXACT only.  In particular, no polynomial Keller composite is
constructed, and the conditional boundary calculation is not a JC(2) proof.
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


CHECKS = 0


def check(condition: object, label: str) -> None:
    """Exact check which remains active under python -O."""

    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.factor(left - right)) == 0, label)


def jac(left: sp.Expr, right: sp.Expr, first: sp.Symbol, second: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, first) * sp.diff(right, second)
        - sp.diff(left, second) * sp.diff(right, first)
    )


def sylvester_quadratic(
    f2: sp.Expr,
    f1: sp.Expr,
    f0: sp.Expr,
    g2: sp.Expr,
    g1: sp.Expr,
    g0: sp.Expr,
) -> sp.Expr:
    """Resultant of two quadratics, built independently as a determinant."""

    matrix = sp.Matrix(
        [
            [f2, f1, f0, 0],
            [0, f2, f1, f0],
            [g2, g1, g0, 0],
            [0, g2, g1, g0],
        ]
    )
    return sp.factor(matrix.det(method="domain-ge"))


def order_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_terms = sp.Poly(numerator, variable).terms()
    denominator_terms = sp.Poly(denominator, variable).terms()
    check(bool(numerator_terms), "valuation numerator nonzero")
    check(bool(denominator_terms), "valuation denominator nonzero")
    return min(monomial[0] for monomial, _ in numerator_terms) - min(
        monomial[0] for monomial, _ in denominator_terms
    )


# -------------------------------------------------------------------------
# 1. Universal seed, direct Jacobian, and determinant elimination.
# -------------------------------------------------------------------------
x, u, A, C = sp.symbols("x u A C")
p0, p1, q0, q1 = sp.symbols("p0 p1 q0 q1")
p = p1 * u + p0
q = q1 * u + q0

seed_A = u**2 - x + p * x**2
seed_C = u**3 - u - sp.Rational(3, 2) * u * x + q * x**2
D = sp.factor(jac(seed_A, seed_C, x, u))

D_by_coefficients = (
    1
    + (
        sp.Rational(3, 2)
        + (6 * u**2 - 2) * p
        - 4 * u * q
    ) * x
    + (-sp.diff(q, u) - 3 * p + sp.Rational(3, 2) * u * sp.diff(p, u)) * x**2
    + 2 * (p * sp.diff(q, u) - sp.diff(p, u) * q) * x**3
)
same(D, D_by_coefficients, "direct seed Jacobian coefficient expansion")
same(D.subs(x, 0), 1, "fixed-arm different has constant term one")

U = u**2 - A
V = u**3 - u - C
F2, F1, F0 = p, -1, U
G2, G1, G0 = q, -sp.Rational(3, 2) * u, V
R = sylvester_quadratic(F2, F1, F0, G2, G1, G0)
L = sp.factor(sp.Rational(3, 2) * p * u - q)
M = sp.factor(q * U - p * V)

# qF-pG=Lx+M, so substituting x=-M/L into F predicts this determinant.
same(p * R, p * M**2 + L * M + U * L**2,
     "Sylvester determinant equals compact elimination identity")
same(
    (sp.diff(R, u) - L * D).subs({A: seed_A, C: seed_C}, simultaneous=True),
    0,
    "recovery derivative on the seed equals L times the seed Jacobian",
)

# Two unrelated rational specializations check the determinant convention.
for values in (
    {p0: 2, p1: 3, q0: -1, q1: 5, u: 7, A: 11, C: -2},
    {p0: -3, p1: 1, q0: 4, q1: -2, u: -5, A: 6, C: 9},
):
    f_numeric = sp.Poly((p * x**2 - x + U).subs(values), x)
    g_numeric = sp.Poly(
        (q * x**2 - sp.Rational(3, 2) * u * x + V).subs(values), x
    )
    same(
        R.subs(values),
        sp.resultant(f_numeric.as_expr(), g_numeric.as_expr(), x),
        "numeric determinant/resultant convention",
    )


# -------------------------------------------------------------------------
# 2. Exhaustive affine-p,q leading-coefficient strata.
# -------------------------------------------------------------------------
R_affine_poly = sp.Poly(sp.expand(R), u)
check(R_affine_poly.degree() == 8, "p1 chart has generic degree eight")
same(R_affine_poly.coeff_monomial(u**8), p1**2,
     "p1 chart leading coefficient")

a, d, e = sp.symbols("a d e")
R_constant_p = sp.factor(R.subs({p1: 0, p0: a, q1: d, q0: e}))
constant_poly = sp.Poly(sp.expand(R_constant_p), u)
check(constant_poly.degree() == 6, "constant-p chart generic degree six")
same(constant_poly.coeff_monomial(u**6), (a - d) ** 2,
     "constant-p chart leading coefficient")

R_first_resonance = sp.factor(R_constant_p.subs(d, a))
first_poly = sp.Poly(sp.expand(R_first_resonance), u)
check(first_poly.degree() == 4, "d=a resonance has degree four")
same(first_poly.coeff_monomial(u**4), (a + 4 * e**2) / 4,
     "d=a resonance leading coefficient")
same(first_poly.coeff_monomial(u**6), 0, "d=a cancels degree six")
same(first_poly.coeff_monomial(u**5), 0, "d=a cancels degree five")

R_exceptional = sp.factor(R_first_resonance.subs(a, -4 * e**2))
exceptional_poly = sp.Poly(sp.expand(R_exceptional), u)
check(exceptional_poly.degree() == 3, "second resonance has degree three")
exceptional_lc = e * (16 * e**2 * (A - 1) - 1) / 2
same(exceptional_poly.coeff_monomial(u**3), exceptional_lc,
     "second resonance target-dependent leading coefficient")
same(exceptional_poly.coeff_monomial(u**4), 0,
     "second resonance cancels degree four")

# p=0 is a different degree chart: x=u^2-A, so eliminate directly.
H_pzero = sp.expand(
    u**3 - u - sp.Rational(3, 2) * u * (u**2 - A)
    + (q1 * u + q0) * (u**2 - A) ** 2 - C
)
pzero_poly = sp.Poly(H_pzero, u)
check(pzero_poly.degree() == 5, "p=0,q1 nonzero chart degree five")
same(pzero_poly.coeff_monomial(u**5), q1,
     "p=0,q1 nonzero scalar leading coefficient")
H_qconstant = sp.expand(H_pzero.subs(q1, 0))
check(sp.Poly(H_qconstant, u).degree() == 4,
      "p=0,q constant nonzero chart degree four")
same(sp.Poly(H_qconstant, u).coeff_monomial(u**4), q0,
     "p=0,q constant scalar leading coefficient")
H_fixed = sp.expand(H_qconstant.subs(q0, 0))
check(sp.Poly(H_fixed, u).degree() == 3, "fixed seed chart degree three")
same(sp.Poly(H_fixed, u).coeff_monomial(u**3), -sp.Rational(1, 2),
     "fixed seed scalar leading coefficient")

# These mutually exclusive cases exhaust affine p,q over a characteristic-zero
# field.  Only the second resonance has a leading coefficient depending on A.
same(
    exceptional_lc.subs(A, 1 + 1 / (16 * e**2)),
    0,
    "exceptional leading divisor",
)


# -------------------------------------------------------------------------
# 3. Independent valuation witness for nonintegrality of hidden u.
# -------------------------------------------------------------------------
p_exc = -4 * e**2
q_exc = -4 * e**2 * u + e
A_exc = u**2 - x + p_exc * x**2
C_exc = u**3 - u - sp.Rational(3, 2) * u * x + q_exc * x**2
D_exc = sp.factor(jac(A_exc, C_exc, x, u))

t, v = sp.symbols("t v")
u_pole = 1 / t
x_pole = 1 / (2 * e * t) - 1 / (8 * e**2) - t / (4 * e) + v * t**2
A_pole = sp.cancel(A_exc.subs({u: u_pole, x: x_pole}))
C_pole = sp.cancel(C_exc.subs({u: u_pole, x: x_pole}))
A_boundary = 1 + 1 / (16 * e**2)
C_boundary = (1 - 8 * e**2) / (64 * e**3)

A_pole_formula = (
    A_boundary - 4 * e * t * v - t**2 / 4 + 2 * e * t**3 * v
    - 4 * e**2 * t**4 * v**2
)
C_pole_formula = (
    C_boundary - 4 * e * v
    + t * (v / 2 - sp.Rational(1, 4) + 1 / (16 * e**2))
    + t**2 * (2 * e * v - v / (4 * e) + 1 / (16 * e))
    + t**3 * (-4 * e**2 * v**2 - v / 2)
    + e * t**4 * v**2
)
same(A_pole, A_pole_formula, "hidden-u pole chart A expansion")
same(C_pole, C_pole_formula, "hidden-u pole chart C expansion")
check(A_pole.is_polynomial(t, v), "hidden-u pole chart A is regular")
check(C_pole.is_polynomial(t, v), "hidden-u pole chart C is regular")
check(order_at_zero(u_pole, t) == -1, "hidden u has negative valuation")
check(order_at_zero(x_pole, t) == -1, "hidden x has negative valuation")
same(jac(x_pole, u_pole, t, v), 1, "hidden-u pole chart is symplectic")
D_pole = sp.cancel(D_exc.subs({u: u_pole, x: x_pole}))
same(D_pole, jac(A_pole_formula, C_pole_formula, t, v),
     "pole-chart and seed Jacobians agree")
same(D_pole.subs(t, 0), 16 * e**2 * v,
     "hidden-u boundary different")
same(A_pole.subs(t, 0), A_boundary,
     "hidden-u boundary maps to constant A")
same(C_pole.subs(t, 0), C_boundary - 4 * e * v,
     "hidden-u boundary maps affinely in C")

# A,C lie in the t-adic valuation ring Q(e,v)[[t]], while u=t^-1 does not.
# Hence this witness proves nonintegrality of u over k[A,C] in the exceptional
# seed, not merely failure of the displayed resultant to be monic.


# -------------------------------------------------------------------------
# 4. The distinct x-only failure: u and w remain integral.
# -------------------------------------------------------------------------
alpha, u0 = sp.symbols("alpha u0")
g = u - u0
p_escape = alpha * g
q_escape = sp.Rational(3, 2) * alpha * u0 * g
A_escape = u**2 - x + p_escape * x**2
C_escape = u**3 - u - sp.Rational(3, 2) * u * x + q_escape * x**2
D_escape = sp.factor(jac(A_escape, C_escape, x, u))
w = alpha * g * x
X = u**2 - A_escape
same(alpha * g * X, w * (1 - w), "x-only Danielewski relation")

f = u**3 + (2 - 3 * A) * u + 2 * C
K = -alpha * f / 3
W_recovered = K + alpha * g * (u**2 - A)
H = sp.expand(W_recovered**2 - K)
same(f.subs({A: A_escape, C: C_escape}, simultaneous=True),
     -3 * w**2 / alpha, "x-only square packet")
same(W_recovered.subs({A: A_escape, C: C_escape}, simultaneous=True),
     w, "regular sidecar w recovery")
same(H.subs({A: A_escape, C: C_escape}, simultaneous=True), 0,
     "hidden-u monic relation vanishes")
check(sp.Poly(H, u).degree() == 6, "x-only hidden-u relation degree six")
same(sp.Poly(H, u).LC(), 4 * alpha**2 / 9,
     "x-only hidden-u relation scalar leading coefficient")
same(
    (sp.diff(H, u) - sp.Rational(2, 3) * alpha * D_escape).subs(
        {A: A_escape, C: C_escape}, simultaneous=True
    ),
    0,
    "x-only recovery derivative/different identity",
)

# Rebuild its resultant from the same determinant and expose the spurious
# vertical factor introduced when the quadratic leading coefficient vanishes.
R_escape = sylvester_quadratic(
    p_escape, -1, u**2 - A,
    q_escape, -sp.Rational(3, 2) * u, u**3 - u - C,
)
same(R_escape, sp.Rational(9, 4) * g**2 * H,
     "x-only Sylvester factor stripping")

y = sp.symbols("y")
u_escape_pole = u0 + t
x_escape_pole = 1 / (alpha * t) + y
A_escape_pole = sp.cancel(A_escape.subs({u: u_escape_pole, x: x_escape_pole}))
C_escape_pole = sp.cancel(C_escape.subs({u: u_escape_pole, x: x_escape_pole}))
A_escape_formula = u_escape_pole**2 + y + alpha * t * y**2
C_escape_formula = (
    u_escape_pole**3 - u_escape_pole - 3 / (2 * alpha)
    + sp.Rational(3, 2) * (2 * u0 - u_escape_pole) * y
    + sp.Rational(3, 2) * alpha * u0 * t * y**2
)
same(A_escape_pole, A_escape_formula, "x-only pole chart A expansion")
same(C_escape_pole, C_escape_formula, "x-only pole chart C expansion")
check(A_escape_pole.is_polynomial(t, y), "x-only pole chart A regular")
check(C_escape_pole.is_polynomial(t, y), "x-only pole chart C regular")
check(order_at_zero(x_escape_pole, t) == -1, "x-only hidden x pole")
same(jac(x_escape_pole, u_escape_pole, t, y), -1,
     "x-only pole chart constant Jacobian")
D_escape_pole = sp.cancel(D_escape.subs({u: u_escape_pole, x: x_escape_pole}))
same(D_escape_pole.subs(t, 0), -(3 * y + 2) / 2,
     "x-only pole-branch different")

Xline = sp.symbols("Xline")
same(D_escape.subs({u: u0, x: Xline}), (3 * Xline + 2) / 2,
     "x-only regular-branch different")
same(D_escape_pole.subs({t: 0, y: -Xline}), (3 * Xline - 2) / 2,
     "x-only pole-branch different in X coordinate")
check(sp.Poly((3 * Xline + 2) / 2, Xline).degree() == 1,
      "regular boundary residue is nonconstant")
check(sp.Poly((3 * Xline - 2) / 2, Xline).degree() == 1,
      "pole boundary residue is nonconstant")

# Conditional bracket algebra only: for a Keller realization, chain rule gives
# D{u,w}=-alpha*g.  Bracketing alpha*g*X=w(1-w) then gives
# D{u,X}=2w-1.  The script checks the algebra, but does not establish that a
# source divisor realizes either full X-line, so this is not a JC conclusion.
h, bracket_uX = sp.symbols("h bracket_uX")
same(
    (D_escape * h + alpha * g).subs(
        h, alpha * g * bracket_uX / (1 - 2 * w)
    ) * (1 - 2 * w) / (alpha * g),
    D_escape * bracket_uX - (2 * w - 1),
    "conditional Keller boundary bracket identity",
)


# -------------------------------------------------------------------------
# 5. Unique quadratic jet gauge and exact-factor warnings.
# -------------------------------------------------------------------------
eta, xi = sp.symbols("eta xi")
jet_matrix = sp.Matrix(
    [[-1, 2 * u], [-sp.Rational(3, 2) * u, 3 * u**2 - 1]]
)
same(jet_matrix.det(), 1, "quadratic normal-jet system determinant")
jet_solution = jet_matrix.inv() * sp.Matrix([p, q])
xi_unique = sp.factor(jet_solution[0])
eta_unique = sp.factor(jet_solution[1])
same(eta_unique, sp.Rational(3, 2) * p * u - q,
     "unique quadratic u-gauge coefficient")
same(xi_unique, 2 * u * eta_unique - p,
     "unique quadratic x-gauge coefficient")

x_tilde = x + xi_unique * x**2
u_tilde = u + eta_unique * x**2
fixed_A_after_gauge = sp.expand(u_tilde**2 - x_tilde)
fixed_C_after_gauge = sp.expand(
    u_tilde**3 - u_tilde - sp.Rational(3, 2) * u_tilde * x_tilde
)
same(fixed_A_after_gauge - seed_A, eta_unique**2 * x**4,
     "quadratic gauge A residual starts at order four")
same(
    fixed_C_after_gauge - seed_C,
    -sp.Rational(3, 2) * eta_unique * x**3
    + sp.Rational(3, 2) * eta_unique * p * x**4
    + eta_unique**3 * x**6,
    "quadratic gauge C residual starts at order three",
)

# eta=L=0 is the exact, not merely jet-level, fixed-seed factorization.
p_cosmetic = p0 + p1 * u
q_cosmetic = sp.Rational(3, 2) * u * p_cosmetic
chi_polynomial = x - p_cosmetic * x**2
same(u**2 - x + p_cosmetic * x**2, u**2 - chi_polynomial,
     "L=0 exact A factorization")
same(
    u**3 - u - sp.Rational(3, 2) * u * x + q_cosmetic * x**2,
    u**3 - u - sp.Rational(3, 2) * u * chi_polynomial,
    "L=0 exact C factorization",
)

# A separate rational Mobius change also factors through the fixed seed and
# keeps its original monic cubic recovery of u.
lam = sp.symbols("lam")
chi_mobius = x / (1 + lam * x)
A_mobius = u**2 - chi_mobius
C_mobius = u**3 - u - sp.Rational(3, 2) * u * chi_mobius
D_mobius = sp.factor(jac(A_mobius, C_mobius, x, u))
same(D_mobius, (1 + sp.Rational(3, 2) * chi_mobius) / (1 + lam * x) ** 2,
     "Mobius change different factorization")
same(
    (u**3 + (2 - 3 * A) * u + 2 * C).subs(
        {A: A_mobius, C: C_mobius}, simultaneous=True
    ),
    0,
    "Mobius family retains fixed monic recovery cubic",
)


# -------------------------------------------------------------------------
# 6. Frozen packet and optimization-safe source check.
# -------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
check(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "source contains no optimization-sensitive assert",
)

semantic = {
    "classification": (
        "affine p,q exhaustive: p1; constant p generic; d=a; "
        "a=d=-4e2; p=0 q-linear/q-constant/q-zero"
    ),
    "elimination": "independent 4x4 Sylvester determinant; pR=pM2+LM+UL2",
    "different": "R_u_on_seed=L*D",
    "exception": "a=d=-4e2; valuation witness proves hidden u nonintegral",
    "x_only": "u,w integral; x can pole; two nonconstant boundary differents",
    "boundary_scope": "conditional bracket identity only; no Keller composite",
    "gauge": "unique quadratic jet; L=0 exact factor; Mobius exact rational factor",
    "status": "INDEPENDENTLY_VERIFIED_FINITE_EXACT;JC2_OPEN",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_RUSSELL_ALTERNATIVE_QUADRATIC_SEED_RECOVERY_INDEPENDENT_AUDIT")
print("status=INDEPENDENTLY_VERIFIED_FINITE_EXACT;SCOUT_SCOPE;JC2_OPEN")
print("elimination=4x4_SYLVESTER;pR=pM^2+LM+UL^2;R_u=L*D_on_seed")
print("affine_classification=EXHAUSTIVE_OVER_CHAR0")
print("exception=a=d=-4e^2;hidden_u_nonintegral_by_t_adic_valuation")
print("hidden_u_boundary=A_constant;C_affine;D=16e^2*v_generically_nonzero")
print("x_only=u_and_w_integral;x_poles;boundary_D=(3X+2)/2,(3X-2)/2")
print("boundary_scope=NONCONSTANT_RESIDUES+CONDITIONAL_BRACKET;NO_JC_CONCLUSION")
print("gauge=UNIQUE_QUADRATIC_JET;L0_EXACT_FACTOR;MOBIUS_EXACT_FACTOR")
print(f"semantic_sha256={semantic_sha}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
