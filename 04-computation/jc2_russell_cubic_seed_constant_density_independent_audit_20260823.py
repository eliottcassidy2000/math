#!/usr/bin/env python3
"""Independent finite-exact audit of the cubic Russell-seed scout.

This implementation deliberately avoids the primary script's ``solve``-based
classification.  It derives the arm action by Taylor columns, solves the raw
affine rows by full-rank coefficient matrices, rebuilds the section-free
degree cascade from the universal Jacobian, and reads pole orders directly
from numerator/denominator valuations.

Status: FINITE-EXACT / INDEPENDENT AUDIT.  The audited cells are empty; no
Keller map is constructed and JC(2) remains open.
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
    """Optimization-safe truth gate."""

    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    difference = left - right
    if isinstance(difference, sp.MatrixBase):
        check(all(sp.cancel(sp.factor(v)) == 0 for v in difference), label)
    else:
        check(sp.cancel(sp.factor(difference)) == 0, label)


def jacobian(left: sp.Expr, right: sp.Expr, first: sp.Symbol,
             second: sp.Symbol) -> sp.Expr:
    matrix = sp.Matrix(
        [
            [sp.diff(left, first), sp.diff(left, second)],
            [sp.diff(right, first), sp.diff(right, second)],
        ]
    )
    return sp.expand(matrix.det())


def coefficient_system(polynomial: sp.Expr, variable: sp.Symbol,
                       unknowns: list[sp.Symbol]) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the full coefficient system ``matrix*unknowns=rhs``."""

    poly = sp.Poly(sp.expand(polynomial), variable)
    equations = [
        poly.coeff_monomial(variable**degree)
        for degree in range(poly.degree() + 1)
    ]
    return sp.linear_eq_to_matrix(equations, unknowns)


def t_order_and_lead(expression: sp.Expr, t: sp.Symbol) -> tuple[int, sp.Expr]:
    """Exact order and leading coefficient at t=0 in an EX coefficient field."""

    numerator, denominator = sp.fraction(sp.cancel(expression))
    num_poly = sp.Poly(numerator, t, domain="EX")
    den_poly = sp.Poly(denominator, t, domain="EX")
    num_order = min(monomial[0] for monomial, _ in num_poly.terms())
    den_order = min(monomial[0] for monomial, _ in den_poly.terms())
    lead = sp.cancel(
        num_poly.coeff_monomial(t**num_order)
        / den_poly.coeff_monomial(t**den_order)
    )
    return num_order - den_order, lead


x, u = sp.symbols("x u")
A0 = u**2 - x
C0 = u**3 - u - sp.Rational(3, 2) * u * x
D0 = jacobian(A0, C0, x, u)
same(D0, 1 + sp.Rational(3, 2) * x, "normalized seed density")


# -------------------------------------------------------------------------
# 1. Taylor-column derivation of the filtered source action.
# -------------------------------------------------------------------------
Xi, Eta = sp.symbols("Xi Eta")
expected_arm_matrix = sp.Matrix(
    [
        [-1, 2 * u],
        [-sp.Rational(3, 2) * u, 3 * u**2 - 1],
    ]
)

for grade in (2, 3, 7):
    X_grade = x + Xi * x**grade
    U_grade = u + Eta * x**grade
    A_grade = sp.expand(A0.subs({x: X_grade, u: U_grade}, simultaneous=True))
    C_grade = sp.expand(C0.subs({x: X_grade, u: U_grade}, simultaneous=True))
    grade_vector = sp.Matrix(
        [A_grade.coeff(x, grade), C_grade.coeff(x, grade)]
    )
    derived_matrix = grade_vector.jacobian([Xi, Eta])
    same(derived_matrix, expected_arm_matrix,
         f"Taylor columns reproduce arm matrix at grade {grade}")

M = expected_arm_matrix
same(M.det(), 1, "arm matrix is unimodular over k[u]")
same(
    M.inv(),
    sp.Matrix(
        [
            [3 * u**2 - 1, -2 * u],
            [sp.Rational(3, 2) * u, -1],
        ]
    ),
    "polynomial inverse of arm matrix",
)

xi2, eta2, xi3, eta3 = sp.symbols("xi2 eta2 xi3 eta3")
X23 = x + xi2 * x**2 + xi3 * x**3
U23 = u + eta2 * x**2 + eta3 * x**3
A23 = sp.expand(A0.subs({x: X23, u: U23}, simultaneous=True))
C23 = sp.expand(C0.subs({x: X23, u: U23}, simultaneous=True))
same(A23.coeff(x, 2), (M * sp.Matrix([xi2, eta2]))[0],
     "quadratic target action")
same(C23.coeff(x, 2), (M * sp.Matrix([xi2, eta2]))[1],
     "quadratic target action in second coordinate")
same(A23.coeff(x, 3), (M * sp.Matrix([xi3, eta3]))[0],
     "fresh cubic target action")
same(
    C23.coeff(x, 3),
    (M * sp.Matrix([xi3, eta3]))[1] - sp.Rational(3, 2) * eta2,
    "quadratic-to-cubic carry",
)

rhoA, rhoC = sp.symbols("rhoA rhoC")
same(M * (M.inv() * sp.Matrix([rhoA, rhoC])),
     sp.Matrix([rhoA, rhoC]), "every fresh residual is movable")


# -------------------------------------------------------------------------
# 2. Determinant cocycle and a logically decisive orbit witness.
# -------------------------------------------------------------------------
xi = sp.Function("xi")(u)
eta = sp.Function("eta")(u)
X_phi = x + xi * x**3
U_phi = u + eta * x**3
A_phi = sp.expand(A0.subs({x: X_phi, u: U_phi}, simultaneous=True))
C_phi = sp.expand(C0.subs({x: X_phi, u: U_phi}, simultaneous=True))
J_phi = jacobian(X_phi, U_phi, x, u)
same(
    jacobian(A_phi, C_phi, x, u),
    D0.subs(x, X_phi) * J_phi,
    "full symbolic determinant cocycle",
)
same(
    J_phi,
    1 + 3 * xi * x**2 + sp.diff(eta, u) * x**3
    + 3 * (xi * sp.diff(eta, u) - sp.diff(xi, u) * eta) * x**5,
    "source determinant through its top cubic-source bucket",
)

chi_exact = sp.Rational(2, 3) * (sp.sqrt(1 + 3 * x) - 1)
A_darboux = A0.subs(x, chi_exact)
C_darboux = C0.subs(x, chi_exact)
same(jacobian(A_darboux, C_darboux, x, u), 1,
     "constant-density representative in the same formal I1 orbit")
same(sp.series(chi_exact - x, x, 0, 3).removeO(),
     -sp.Rational(3, 4) * x**2,
     "formal representative differs first inside I1")
check(D0 != 1, "identity representative in that orbit is nonconstant-density")


# -------------------------------------------------------------------------
# 3. Raw affine classification by invertible coefficient matrices.
# -------------------------------------------------------------------------
p0, p1, q0, q1, r0, r1, s0, s1 = sp.symbols(
    "p0 p1 q0 q1 r0 r1 s0 s1"
)
p = p0 + p1 * u
q = q0 + q1 * u
r = r0 + r1 * u
s = s0 + s1 * u
A_raw = A0 + p * x**2 + r * x**3
C_raw = C0 + q * x**2 + s * x**3
D_raw = jacobian(A_raw, C_raw, x, u)
check(sp.Poly(D_raw, x).degree() == 5, "universal raw cell has buckets E0..E5")

unknowns_E1 = [p0, p1, q0, q1]
matrix_E1, rhs_E1 = coefficient_system(D_raw.coeff(x, 1), u, unknowns_E1)
check(matrix_E1.shape == (4, 4), "E1 is a square affine coefficient system")
check(matrix_E1.det() != 0, "E1 coefficient system has full rank")
solution_E1 = matrix_E1.inv() * rhs_E1
same(solution_E1, sp.Matrix([sp.Rational(3, 4), 0, 0,
                             sp.Rational(9, 8)]),
     "unique raw E1 solution")
subs_E1 = dict(zip(unknowns_E1, solution_E1))

unknowns_E2 = [r0, r1, s0, s1]
E2_after_E1 = sp.expand(D_raw.subs(subs_E1)).coeff(x, 2)
matrix_E2, rhs_E2 = coefficient_system(E2_after_E1, u, unknowns_E2)
check(matrix_E2.shape == (4, 4), "E2 is a square affine coefficient system")
check(matrix_E2.det() != 0, "E2 coefficient system has full rank")
solution_E2 = matrix_E2.inv() * rhs_E2
same(solution_E2, sp.Matrix([-sp.Rational(9, 8), 0, 0,
                             -sp.Rational(27, 16)]),
     "unique raw E2 continuation")
subs_E2 = dict(zip(unknowns_E2, solution_E2))

raw_solution = {**subs_E1, **subs_E2}
A_candidate = sp.expand(A_raw.subs(raw_solution))
C_candidate = sp.expand(C_raw.subs(raw_solution))
D_candidate = sp.expand(D_raw.subs(raw_solution))
chi3 = x - sp.Rational(3, 4) * x**2 + sp.Rational(9, 8) * x**3
same(A_candidate, A0.subs(x, chi3), "candidate is a fixed-seed precomposition")
same(C_candidate, C0.subs(x, chi3),
     "second candidate coordinate is the same precomposition")
same(D_candidate, (1 + sp.Rational(3, 2) * chi3) * sp.diff(chi3, x),
     "candidate density factors by the one-variable chain rule")
same(D_candidate.coeff(x, 3), sp.Rational(135, 16),
     "raw affine obstruction occurs at E3")
check(D_candidate != 1, "raw affine constant-density cell is empty")


# -------------------------------------------------------------------------
# 4. Reversion coefficient, monic recovery, and the nonsquare wall.
# -------------------------------------------------------------------------
c2, c3, c4 = sp.symbols("c2 c3 c4")
chi4_unknown = x + c2 * x**2 + c3 * x**3 + c4 * x**4
density4 = sp.expand(
    (1 + sp.Rational(3, 2) * chi4_unknown) * sp.diff(chi4_unknown, x)
)
reversion_equations = [density4.coeff(x, degree) for degree in (1, 2, 3)]


def unique_linear_root(expression: sp.Expr, variable: sp.Symbol,
                       label: str) -> sp.Expr:
    """Solve one nondegenerate affine equation without ``solve``."""

    check(sp.diff(expression, variable, 2) == 0, f"{label} is affine")
    slope = sp.diff(expression, variable)
    check(slope != 0, f"{label} has nonzero slope")
    return sp.cancel(-expression.subs(variable, 0) / slope)


c2_value = unique_linear_root(reversion_equations[0], c2, "reversion row c2")
c3_value = unique_linear_root(
    reversion_equations[1].subs(c2, c2_value), c3, "reversion row c3"
)
c4_value = unique_linear_root(
    reversion_equations[2].subs({c2: c2_value, c3: c3_value}),
    c4,
    "reversion row c4",
)
solution_rev = sp.Matrix([c2_value, c3_value, c4_value])
same(solution_rev, sp.Matrix([-sp.Rational(3, 4), sp.Rational(9, 8),
                              -sp.Rational(135, 64)]),
     "independent density-ODE reversion coefficients")
same(chi_exact + sp.Rational(3, 4) * chi_exact**2, x,
     "integrated quadratic relation")
same((1 + sp.Rational(3, 2) * chi_exact) * sp.diff(chi_exact, x), 1,
     "exact Catalan/Hensel resummation")
same(sp.series(chi_exact, x, 0, 5).removeO(),
     chi3 - sp.Rational(135, 64) * x**4,
     "first omitted coefficient")
check(sp.Poly(1 + 3 * x, x).degree() == 1, "radicand has odd simple valuation")
check(sp.gcd(1 + 3 * x, sp.diff(1 + 3 * x, x)) == 1,
      "radicand divisor is squarefree")

A_t, C_t, Chi_t = sp.symbols("A_t C_t Chi_t")
u_monic = u**3 + (2 - 3 * A_t) * u + 2 * C_t
same(u_monic.subs({A_t: A_candidate, C_t: C_candidate}, simultaneous=True),
     0, "monic cubic recovers u integrally")
same((u**2 - A_t).subs(A_t, A_candidate), chi3,
     "chi3 recovers after u")
x_monic = x**3 - sp.Rational(2, 3) * x**2 + sp.Rational(8, 9) * x \
    - sp.Rational(8, 9) * Chi_t
same(x_monic.subs(Chi_t, chi3), 0, "monic cubic recovers x integrally")


# -------------------------------------------------------------------------
# 5. Full E1 fibre and section-independent affine-residual cascade.
# -------------------------------------------------------------------------
eta_fun = sp.Function("eta_fun")(u)
alpha_fun = sp.Function("alpha_fun")(u)
beta_fun = sp.Function("beta_fun")(u)
p_fibre = sp.Rational(3, 4) + 2 * u * eta_fun
q_fibre = sp.Rational(9, 8) * u + (3 * u**2 - 1) * eta_fun
r_fibre = alpha_fun
s_fibre = beta_fun - sp.Rational(3, 2) * eta_fun
A_fibre = A0 + p_fibre * x**2 + r_fibre * x**3
C_fibre = C0 + q_fibre * x**2 + s_fibre * x**3
D_fibre = jacobian(A_fibre, C_fibre, x, u)
same(D_fibre.coeff(x, 1), 0, "parameterized full E1 fibre")
check(sp.gcd(2 * u, 3 * u**2 - 1) == 1,
      "coprime syzygy coefficients make the E1 parameter exhaustive")
same(
    D_fibre.coeff(x, 2),
    sp.diff(eta_fun, u) - sp.Rational(27, 8)
    + 3 * (3 * u**2 - 1) * alpha_fun - 6 * u * beta_fun,
    "E2 on the full E1 fibre",
)
same(
    D_fibre.coeff(x, 3),
    -sp.diff(beta_fun, u) + sp.Rational(27, 16) + 9 * u * eta_fun
    + (12 * u**2 + 4) * eta_fun**2
    - sp.Rational(9, 2) * alpha_fun
    + sp.Rational(3, 2) * u * sp.diff(alpha_fun, u),
    "E3 on the full E1 fibre",
)

a0, a1, b0, b1, h0, h1, h2, h3, h4 = sp.symbols(
    "a0 a1 b0 b1 h0 h1 h2 h3 h4"
)
alpha_affine = a0 + a1 * u
beta_affine = b0 + b1 * u
eta_poly = h0 + h1 * u + h2 * u**2 + h3 * u**3 + h4 * u**4
E2_poly = sp.expand(
    D_fibre.coeff(x, 2).subs(
        {
            eta_fun: eta_poly,
            sp.diff(eta_fun, u): sp.diff(eta_poly, u),
            alpha_fun: alpha_affine,
            beta_fun: beta_affine,
        },
        simultaneous=True,
    )
)
E3_poly = sp.expand(
    D_fibre.coeff(x, 3).subs(
        {
            eta_fun: eta_poly,
            alpha_fun: alpha_affine,
            sp.diff(alpha_fun, u): a1,
            beta_fun: beta_affine,
            sp.diff(beta_fun, u): b1,
        },
        simultaneous=True,
    )
)

cascade = [(h4, 10), (h3, 8), (h2, 6), (h1, 4), (h0, 2)]
zeroed: dict[sp.Symbol, int] = {}
for coefficient, degree in cascade:
    current = sp.expand(E3_poly.subs(zeroed))
    same(sp.Poly(current, u).coeff_monomial(u**degree), 12 * coefficient**2,
         f"E3 descending degree gate for {coefficient}")
    zeroed[coefficient] = 0

E2_eta_zero = sp.expand(E2_poly.subs(zeroed))
matrix_residual, rhs_residual = coefficient_system(
    E2_eta_zero, u, [a0, a1, b0, b1]
)
check(matrix_residual.det() != 0, "eta=0 leaves a full-rank E2 system")
solution_residual = matrix_residual.inv() * rhs_residual
same(solution_residual,
     sp.Matrix([-sp.Rational(9, 8), 0, 0, -sp.Rational(27, 16)]),
     "unique affine residual after the cascade")
residual_subs = dict(zip([a0, a1, b0, b1], solution_residual))
same(E3_poly.subs({**zeroed, **residual_subs}), sp.Rational(135, 16),
     "section-independent residual contradiction")


# -------------------------------------------------------------------------
# 6. Pole-chart effectivity by exact t-valuations.
# -------------------------------------------------------------------------
e = sp.symbols("e", nonzero=True)
A_exceptional = A0 - 4 * e**2 * x**2
C_exceptional = C0 + (-4 * e**2 * u + e) * x**2
D_exceptional = jacobian(A_exceptional, C_exceptional, x, u)
expected_E1_exceptional = (
    sp.Rational(3, 2) + 8 * e**2 - 4 * e * u - 8 * e**2 * u**2
)
same(D_exceptional.coeff(x, 1), expected_E1_exceptional,
     "exceptional seed has immutable nonconstant E1")

t, v = sp.symbols("t v")
u_pole = 1 / t
x_pole = 1 / (2 * e * t) - 1 / (8 * e**2) - t / (4 * e) + v * t**2
A_boundary = sp.cancel(A_exceptional.subs({x: x_pole, u: u_pole}))
C_boundary = sp.cancel(C_exceptional.subs({x: x_pole, u: u_pole}))
check(A_boundary.is_polynomial(t, v), "A is regular on the hostile pole chart")
check(C_boundary.is_polynomial(t, v), "C is regular on the hostile pole chart")
same(jacobian(x_pole, u_pole, t, v), 1, "pole chart preserves the two-form")
same(sp.cancel(D_exceptional.subs({x: x_pole, u: u_pole})).subs(t, 0),
     16 * e**2 * v, "boundary density is nonconstant")

order_slope, lead_slope = t_order_and_lead(u_pole * x_pole**3, t)
order_constant, lead_constant = t_order_and_lead(x_pole**3, t)
check(order_slope == -4, "affine slope correction has pole order four")
same(lead_slope, 1 / (8 * e**3), "affine slope leading coefficient")
check(order_constant == -3, "affine constant correction has pole order three")
same(lead_constant, 1 / (8 * e**3), "affine constant leading coefficient")

R0, R1, S0, S1 = sp.symbols("R0 R1 S0 S1")
A_cubic = A_exceptional + (R0 + R1 * u) * x**3
C_cubic = C_exceptional + (S0 + S1 * u) * x**3
same(jacobian(A_cubic, C_cubic, x, u).coeff(x, 1),
     expected_E1_exceptional, "cubic additions cannot alter E1")


# -------------------------------------------------------------------------
# 7. Frozen semantic packet and optimized-run safety.
# -------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "source contains no optimization-sensitive assert")

semantic = {
    "status": "FINITE-EXACT/INDEPENDENT-AUDIT;JC2-OPEN",
    "arm_action": "Taylor columns give M=[[-1,2u],[-3u/2,3u2-1]],det=1",
    "cocycle": "J(FoPhi)=J(F)(Phi)J(Phi);Catalan representative proves predicate loss",
    "raw_affine": "full-rank E1/E2 systems;unique continuation;E3=135/16",
    "resummation": "c2=-3/4,c3=9/8,c4=-135/64;linear radicand nonsquare",
    "recovery": "u and x satisfy explicit monic cubics",
    "section_free": "E3 degree cascade kills eta;E2 unique;E3=135/16",
    "pole": "basis valuations -4 and -3 force all affine cubic corrections zero",
    "scope": "no rational residual classification;no new pole chart;no Keller map",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_RUSSELL_CUBIC_SEED_CONSTANT_DENSITY_INDEPENDENT_AUDIT")
print("status=FINITE-EXACT+INDEPENDENT-AUDIT;AFFINE_CELLS_EMPTY;JC2_OPEN")
print("arm_action=TAYLOR_DERIVED;det_M=1;fresh_residuals=ALL_MOVABLE")
print("determinant_cocycle=VERIFIED_SYMBOLICALLY")
print("quotient_predicate_witness=IDENTITY_NONCONSTANT;CATALAN_REPRESENTATIVE_CONSTANT")
print("raw_affine_matrices=E1_FULL_RANK;E2_FULL_RANK")
print("raw_affine_solution=p=3/4;q=9u/8;r=-9/8;s=-27u/16")
print("raw_affine_first_debt=E3=135/16")
print("reversion_coefficients=c2=-3/4;c3=9/8;c4=-135/64")
print("monic_recovery=u_CUBIC;x_CUBIC")
print("section_free_cascade=h4=h3=h2=h1=h0=0;E3=135/16")
print("pole_effectivity=ord_t(u*x^3)=-4;ord_t(x^3)=-3;AFFINE_CORRECTIONS_ZERO")
print("scope=NO_RATIONAL_CLASSIFICATION;NO_NEW_POLE_CHART;NO_KELLER_MAP")
print(f"semantic_sha256={semantic_sha}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
