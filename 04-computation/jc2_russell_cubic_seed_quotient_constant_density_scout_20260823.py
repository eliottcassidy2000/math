#!/usr/bin/env python3
"""Finite-exact cubic-seed quotient and constant-density scout.

Work in the normalized Russell arm coordinates ``(x,u)`` and start from

    A0 = u^2 - x,
    C0 = u^3 - u - (3/2) u x.

The script separates three questions which must not be conflated:

1. the full formal source-addition quotient preserving the arm and first
   normal jet;
2. the raw polynomial affine-coefficient cubic seed cell; and
3. effectivity of cubic corrections on the known nonintegral quadratic pole
   chart.

The full cubic residual is cosmetic: at every normal grade the source action
is the same determinant-one arm matrix.  Constancy of the Jacobian does not
descend through that quotient without retaining the source-determinant
cocycle.  In raw seed space the affine cubic cell contains no constant-
Jacobian point.  The unique point satisfying the first two nonconstant
Jacobian rows is an exact fixed-seed precomposition and has monic recovery.
Finally, every nonzero affine cubic correction destroys regularity on the
known hidden-u pole chart.

Status is FINITE-EXACT/SCOUT.  The script does not classify rational cubic
seeds, nonmonogenic multi-control orders, or polynomial normal depth six, and
it has no consequence for JC(2).
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
    """Optimization-safe exact check."""

    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    difference = left - right
    if isinstance(difference, sp.MatrixBase):
        check(
            all(sp.cancel(sp.factor(entry)) == 0 for entry in difference),
            label,
        )
        return
    check(sp.cancel(sp.factor(difference)) == 0, label)


def jac(
    left: sp.Expr,
    right: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
) -> sp.Expr:
    return sp.expand(
        sp.diff(left, first) * sp.diff(right, second)
        - sp.diff(left, second) * sp.diff(right, first)
    )


x, u = sp.symbols("x u")
A0 = u**2 - x
C0 = u**3 - u - sp.Rational(3, 2) * u * x
D0 = jac(A0, C0, x, u)
same(D0, 1 + sp.Rational(3, 2) * x, "fixed normalized seed Jacobian")
same(A0.subs({x: 0, u: 1}), A0.subs({x: 0, u: -1}),
     "arm collision A at u=+-1")
same(C0.subs({x: 0, u: 1}), C0.subs({x: 0, u: -1}),
     "arm collision C at u=+-1")


# -------------------------------------------------------------------------
# 1. Universal cubic seed and its six exact Jacobian buckets.
# -------------------------------------------------------------------------
p_fun = sp.Function("p")(u)
q_fun = sp.Function("q")(u)
r_fun = sp.Function("r")(u)
s_fun = sp.Function("s")(u)

A_universal = A0 + p_fun * x**2 + r_fun * x**3
C_universal = C0 + q_fun * x**2 + s_fun * x**3
D_universal = sp.expand(jac(A_universal, C_universal, x, u))

E1 = (
    sp.Rational(3, 2)
    + (6 * u**2 - 2) * p_fun
    - 4 * u * q_fun
)
E2 = (
    -sp.diff(q_fun, u)
    - 3 * p_fun
    + sp.Rational(3, 2) * u * sp.diff(p_fun, u)
    + (9 * u**2 - 3) * r_fun
    - 6 * u * s_fun
)
E3 = (
    -sp.diff(s_fun, u)
    + 2 * (p_fun * sp.diff(q_fun, u) - sp.diff(p_fun, u) * q_fun)
    - sp.Rational(9, 2) * r_fun
    + sp.Rational(3, 2) * u * sp.diff(r_fun, u)
)
E4 = (
    2 * p_fun * sp.diff(s_fun, u)
    + 3 * r_fun * sp.diff(q_fun, u)
    - 3 * sp.diff(p_fun, u) * s_fun
    - 2 * sp.diff(r_fun, u) * q_fun
)
E5 = 3 * (r_fun * sp.diff(s_fun, u) - sp.diff(r_fun, u) * s_fun)

same(
    D_universal,
    1 + E1 * x + E2 * x**2 + E3 * x**3 + E4 * x**4 + E5 * x**5,
    "universal cubic seed has exactly six Jacobian buckets",
)


# -------------------------------------------------------------------------
# 2. Full filtered cosmetic action and the quotient loss ledger.
# -------------------------------------------------------------------------
# I_0=x*k[u][[x]]^2 preserves the arm value.  The actual seed packet also
# fixes the first normal jet, so its cosmetic ideal is
# I_1=x^2*k[u][[x]]^2.  At every grade n>=2, the associated-graded action is
# the following arm derivative matrix.
arm_matrix = sp.Matrix(
    [
        [-1, 2 * u],
        [-sp.Rational(3, 2) * u, 3 * u**2 - 1],
    ]
)
check(arm_matrix.det() == 1, "arm addition matrix determinant one")
check(arm_matrix.rank() == 2, "arm addition matrix has full polynomial rank")
same(arm_matrix.inv() * arm_matrix, sp.eye(2), "arm matrix polynomial inverse")

xi2 = sp.Function("xi2")(u)
eta2 = sp.Function("eta2")(u)
xi3 = sp.Function("xi3")(u)
eta3 = sp.Function("eta3")(u)
X_control = x + xi2 * x**2 + xi3 * x**3
U_control = u + eta2 * x**2 + eta3 * x**3
A_control = sp.expand(U_control**2 - X_control)
C_control = sp.expand(
    U_control**3 - U_control - sp.Rational(3, 2) * U_control * X_control
)

control_quadratic = arm_matrix * sp.Matrix([xi2, eta2])
same(A_control.coeff(x, 2), control_quadratic[0], "quadratic A action")
same(C_control.coeff(x, 2), control_quadratic[1], "quadratic C action")

control_cubic = arm_matrix * sp.Matrix([xi3, eta3])
same(A_control.coeff(x, 3), control_cubic[0], "cubic A action")
same(
    C_control.coeff(x, 3),
    control_cubic[1] - sp.Rational(3, 2) * eta2,
    "cubic C action including the quadratic-section carry",
)

p_symbol, q_symbol = sp.symbols("p_symbol q_symbol")
eta_section = sp.factor(sp.Rational(3, 2) * p_symbol * u - q_symbol)
xi_section = sp.factor(2 * u * eta_section - p_symbol)
section_vector = arm_matrix * sp.Matrix([xi_section, eta_section])
same(section_vector[0], p_symbol, "unique quadratic section recovers p")
same(section_vector[1], q_symbol, "unique quadratic section recovers q")

# Relative to that chosen quadratic section, a displayed cubic residual is
# rho=(r,s+3eta/2).  It is not a quotient invariant: a cubic addition with
# coefficients M^{-1}rho realizes it exactly.
rho_A, rho_C = sp.symbols("rho_A rho_C")
cubic_killer = arm_matrix.inv() * sp.Matrix([rho_A, rho_C])
same(
    arm_matrix * cubic_killer,
    sp.Matrix([rho_A, rho_C]),
    "every cubic residual is a full cosmetic addition",
)
same(
    cubic_killer[0],
    (3 * u**2 - 1) * rho_A - 2 * u * rho_C,
    "explicit cubic xi coefficient",
)
same(
    cubic_killer[1],
    sp.Rational(3, 2) * u * rho_A - rho_C,
    "explicit cubic eta coefficient",
)

# The quotient forgets the next native operation.  A cubic source addition
# changes the determinant, so constant-J cannot be tested on the bare orbit.
X_hostile = x + x**3
U_hostile = u
A_hostile = A0.subs({x: X_hostile, u: U_hostile}, simultaneous=True)
C_hostile = C0.subs({x: X_hostile, u: U_hostile}, simultaneous=True)
D_hostile = jac(A_hostile, C_hostile, x, u)
same(
    D_hostile,
    (1 + sp.Rational(3, 2) * X_hostile) * (1 + 3 * x**2),
    "determinant cocycle under a cubic cosmetic addition",
)
check(D_hostile != D0, "constant-density data is lost by the full quotient")


# -------------------------------------------------------------------------
# 3. Raw affine p,q,r,s cell: exact empty constant-J classification.
# -------------------------------------------------------------------------
p0, p1, q0, q1, r0, r1, s0, s1 = sp.symbols(
    "p0 p1 q0 q1 r0 r1 s0 s1"
)
p_affine = p1 * u + p0
q_affine = q1 * u + q0
r_affine = r1 * u + r0
s_affine = s1 * u + s0
A_affine = A0 + p_affine * x**2 + r_affine * x**3
C_affine = C0 + q_affine * x**2 + s_affine * x**3
D_affine = sp.expand(jac(A_affine, C_affine, x, u))

E1_affine = sp.Poly(D_affine.coeff(x, 1), u)
E1_equations = [E1_affine.coeff_monomial(u**degree) for degree in range(4)]
E1_solutions = sp.solve(
    E1_equations, [p1, q1, q0, p0], dict=True, simplify=False
)
expected_E1_solution = {p0: sp.Rational(3, 4), p1: 0,
                        q0: 0, q1: sp.Rational(9, 8)}
check(E1_solutions == [expected_E1_solution],
      "raw affine E1 has the unique quadratic section")

D_after_E1 = sp.expand(D_affine.subs(expected_E1_solution))
E2_after_E1 = sp.Poly(D_after_E1.coeff(x, 2), u)
E2_equations = [E2_after_E1.coeff_monomial(u**degree) for degree in range(4)]
E2_solutions = sp.solve(
    E2_equations, [r1, s1, s0, r0], dict=True, simplify=False
)
expected_E2_solution = {
    r0: -sp.Rational(9, 8),
    r1: 0,
    s0: 0,
    s1: -sp.Rational(27, 16),
}
check(E2_solutions == [expected_E2_solution],
      "raw affine E2 has a unique cubic continuation")

raw_candidate_substitution = {**expected_E1_solution, **expected_E2_solution}
A_candidate = sp.expand(A_affine.subs(raw_candidate_substitution))
C_candidate = sp.expand(C_affine.subs(raw_candidate_substitution))
D_candidate = sp.factor(D_affine.subs(raw_candidate_substitution))

chi = x - sp.Rational(3, 4) * x**2 + sp.Rational(9, 8) * x**3
same(A_candidate, u**2 - chi, "unique two-row candidate A factors through chi")
same(
    C_candidate,
    u**3 - u - sp.Rational(3, 2) * u * chi,
    "unique two-row candidate C factors through chi",
)
expected_candidate_D = (
    1
    + sp.Rational(135, 16) * x**3
    - sp.Rational(405, 64) * x**4
    + sp.Rational(729, 128) * x**5
)
same(D_candidate, expected_candidate_D, "first surviving determinant debt")
same(sp.expand(D_candidate).coeff(x, 3), sp.Rational(135, 16),
     "raw affine E3 obstruction")
check(D_candidate != 1, "raw affine cubic constant-J cell is empty")

# The unique raw candidate is not accidental: it is exactly the cubic
# truncation of the one-variable density ODE
#
#     (1+3*chi/2) chi' = 1,       chi(0)=0.
#
# Its algebraic solution is the canonical Catalan/Hensel sheet from
# THM-3846, in these normalized coordinates.
chi_exact = sp.Rational(2, 3) * (sp.sqrt(1 + 3 * x) - 1)
same(
    chi_exact + sp.Rational(3, 4) * chi_exact**2,
    x,
    "exact density primitive quadratic equation",
)
same(
    (1 + sp.Rational(3, 2) * chi_exact) * sp.diff(chi_exact, x),
    1,
    "exact one-variable density ODE",
)
chi_exact_order4 = sp.series(chi_exact, x, 0, 5).removeO()
same(
    chi_exact_order4,
    chi - sp.Rational(135, 64) * x**4,
    "Catalan resummation through the first omitted coefficient",
)
chi_order4_density = sp.expand(
    (1 + sp.Rational(3, 2) * chi_exact_order4)
    * sp.diff(chi_exact_order4, x)
)
same(chi_order4_density.coeff(x, 3), 0,
     "first omitted coefficient cancels the 135/16 debt")
check(sp.Poly(1 + 3 * x, x).degree() == 1,
      "exact resummation radicand has odd simple degree")

# The same candidate lies on the old fixed-seed monic-recovery side, not on
# the nonintegral side of the frontier.
A_target, C_target, Chi_target = sp.symbols("A_target C_target Chi_target")
u_recovery = u**3 + (2 - 3 * A_target) * u + 2 * C_target
zero_on_candidate = u_recovery.subs(
    {A_target: A_candidate, C_target: C_candidate}, simultaneous=True
)
same(zero_on_candidate, 0, "candidate retains the fixed monic u cubic")
same(
    (u**2 - A_target).subs(A_target, A_candidate),
    chi,
    "candidate recovers chi after u",
)
x_recovery = (
    x**3
    - sp.Rational(2, 3) * x**2
    + sp.Rational(8, 9) * x
    - sp.Rational(8, 9) * Chi_target
)
same(x_recovery.subs(Chi_target, chi), 0, "candidate has monic x recovery")


# -------------------------------------------------------------------------
# 4. Section-independent affine residual cell.
# -------------------------------------------------------------------------
# E1=0 does not select eta=0.  Its full polynomial quadratic fibre is
# p=3/4+2u*eta, q=9u/8+(3u^2-1)eta.  Keep eta, then write the residual after
# the chosen quadratic section as alpha*x^3 and beta*x^3.
eta = sp.Function("eta")(u)
alpha = sp.Function("alpha")(u)
beta = sp.Function("beta")(u)
p_eta = sp.Rational(3, 4) + 2 * u * eta
q_eta = sp.Rational(9, 8) * u + (3 * u**2 - 1) * eta
r_eta = alpha
s_eta = beta - sp.Rational(3, 2) * eta
A_eta = A0 + p_eta * x**2 + r_eta * x**3
C_eta = C0 + q_eta * x**2 + s_eta * x**3
D_eta = sp.expand(jac(A_eta, C_eta, x, u))

same(D_eta.coeff(x, 1), 0, "full quadratic fibre cancels E1")
E2_eta_formula = (
    sp.diff(eta, u)
    - sp.Rational(27, 8)
    + 3 * (3 * u**2 - 1) * alpha
    - 6 * u * beta
)
same(D_eta.coeff(x, 2), E2_eta_formula,
     "section-independent E2 formula")
E3_eta_formula = (
    -sp.diff(beta, u)
    + sp.Rational(27, 16)
    + 9 * u * eta
    + (12 * u**2 + 4) * eta**2
    - sp.Rational(9, 2) * alpha
    + sp.Rational(3, 2) * u * sp.diff(alpha, u)
)
same(D_eta.coeff(x, 3), E3_eta_formula,
     "section-independent E3 formula")

# Let the relative cubic residuals be affine.  E2 integrates eta exhaustively:
# eta' has degree at most three, so every polynomial eta has degree <=4 plus
# one integration constant.
a0, a1, b0, b1, h0 = sp.symbols("a0 a1 b0 b1 h0")
alpha_affine = a1 * u + a0
beta_affine = b1 * u + b0
eta_solution = (
    h0
    + (sp.Rational(27, 8) + 3 * a0) * u
    + (sp.Rational(3, 2) * a1 + 3 * b0) * u**2
    + (-3 * a0 + 2 * b1) * u**3
    - sp.Rational(9, 4) * a1 * u**4
)
same(
    sp.diff(eta_solution, u),
    sp.Rational(27, 8)
    - 3 * (3 * u**2 - 1) * alpha_affine
    + 6 * u * beta_affine,
    "exhaustive E2 integration for affine relative residuals",
)

# The E3 degree mechanism is independent of the E2 parameterization.  For a
# polynomial eta=sum h_i u^i of degree <=4, the quadratic eta^2 term forces
# h4,h3,h2,h1,h0 successively to vanish.
h1, h2, h3, h4 = sp.symbols("h1 h2 h3 h4")
eta_generic = h0 + h1 * u + h2 * u**2 + h3 * u**3 + h4 * u**4
E3_degree_packet = sp.expand(
    (12 * u**2 + 4) * eta_generic**2
    + 9 * u * eta_generic
    - 3 * a1 * u
    - sp.Rational(9, 2) * a0
    - b1
    + sp.Rational(27, 16)
)
same(sp.Poly(E3_degree_packet, u).coeff_monomial(u**10), 12 * h4**2,
     "E3 degree cascade h4")
same(
    sp.Poly(E3_degree_packet.subs(h4, 0), u).coeff_monomial(u**8),
    12 * h3**2,
    "E3 degree cascade h3",
)
same(
    sp.Poly(E3_degree_packet.subs({h4: 0, h3: 0}), u).coeff_monomial(u**6),
    12 * h2**2,
    "E3 degree cascade h2",
)
same(
    sp.Poly(
        E3_degree_packet.subs({h4: 0, h3: 0, h2: 0}), u
    ).coeff_monomial(u**4),
    12 * h1**2,
    "E3 degree cascade h1",
)
same(
    sp.Poly(
        E3_degree_packet.subs({h4: 0, h3: 0, h2: 0, h1: 0}), u
    ).coeff_monomial(u**2),
    12 * h0**2,
    "E3 degree cascade h0",
)

eta_zero_equations = [
    sp.Poly(eta_solution, u).coeff_monomial(u**degree)
    for degree in range(5)
]
eta_zero_solutions = sp.solve(
    eta_zero_equations, [h0, a1, b1, b0, a0], dict=True, simplify=False
)
expected_eta_zero = {
    a0: -sp.Rational(9, 8),
    a1: 0,
    b0: 0,
    b1: -sp.Rational(27, 16),
    h0: 0,
}
check(eta_zero_solutions == [expected_eta_zero],
      "E2 plus eta=0 recovers unique raw continuation")
same(
    E3_degree_packet.subs(
        {**expected_eta_zero, h1: 0, h2: 0, h3: 0, h4: 0}
    ),
    sp.Rational(135, 16),
    "section-independent affine residual contradiction",
)


# -------------------------------------------------------------------------
# 5. The known nonintegral quadratic pole chart rejects cubic additions.
# -------------------------------------------------------------------------
e = sp.symbols("e", nonzero=True)
p_exceptional = -4 * e**2
q_exceptional = -4 * e**2 * u + e
A_exceptional = A0 + p_exceptional * x**2
C_exceptional = C0 + q_exceptional * x**2
D_exceptional = sp.factor(jac(A_exceptional, C_exceptional, x, u))

expected_exceptional_E1 = (
    sp.Rational(3, 2) + 8 * e**2 - 4 * e * u - 8 * e**2 * u**2
)
same(sp.expand(D_exceptional).coeff(x, 1), expected_exceptional_E1,
     "nonintegral exceptional seed already fails E1")
check(sp.Poly(expected_exceptional_E1, u).degree() == 2,
      "exceptional E1 is nonconstant for e nonzero")

t, v = sp.symbols("t v")
u_pole = 1 / t
x_pole = 1 / (2 * e * t) - 1 / (8 * e**2) - t / (4 * e) + v * t**2
A_pole = sp.cancel(A_exceptional.subs({u: u_pole, x: x_pole}))
C_pole = sp.cancel(C_exceptional.subs({u: u_pole, x: x_pole}))
check(A_pole.is_polynomial(t, v), "exceptional pole chart keeps A regular")
check(C_pole.is_polynomial(t, v), "exceptional pole chart keeps C regular")
same(jac(x_pole, u_pole, t, v), 1, "exceptional pole chart is symplectic")
D_exceptional_pole = sp.factor(
    sp.cancel(D_exceptional.subs({u: u_pole, x: x_pole}))
)
same(
    D_exceptional_pole.subs(t, 0),
    16 * e**2 * v,
    "exceptional boundary Jacobian remains nonconstant",
)

ra0, ra1, sa0, sa1 = sp.symbols("ra0 ra1 sa0 sa1")
delta_A = sp.cancel((ra1 * u_pole + ra0) * x_pole**3)
delta_C = sp.cancel((sa1 * u_pole + sa0) * x_pole**3)
same(sp.limit(t**4 * delta_A, t, 0), ra1 / (8 * e**3),
     "affine cubic A slope creates order-four pole")
same(
    sp.limit(t**3 * delta_A.subs(ra1, 0), t, 0),
    ra0 / (8 * e**3),
    "affine cubic A constant creates order-three pole",
)
same(sp.limit(t**4 * delta_C, t, 0), sa1 / (8 * e**3),
     "affine cubic C slope creates order-four pole")
same(
    sp.limit(t**3 * delta_C.subs(sa1, 0), t, 0),
    sa0 / (8 * e**3),
    "affine cubic C constant creates order-three pole",
)

# Cubic terms start at the x^2 Jacobian bucket, so they cannot repair the
# exceptional seed's nonzero x^1 bucket even before the pole test.
A_exceptional_cubic = A_exceptional + (ra1 * u + ra0) * x**3
C_exceptional_cubic = C_exceptional + (sa1 * u + sa0) * x**3
D_exceptional_cubic = jac(A_exceptional_cubic, C_exceptional_cubic, x, u)
same(
    D_exceptional_cubic.coeff(x, 1),
    expected_exceptional_E1,
    "cubic residuals cannot change the exceptional E1 debt",
)


# -------------------------------------------------------------------------
# 6. Frozen scope packet and optimization-safe source audit.
# -------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
check(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "source contains no optimization-sensitive assert",
)

semantic = {
    "normalization": "A0=u2-x;C0=u3-u-3ux/2;D0=1+3x/2",
    "addition_ideals": "I0=x*k[u][[x]]2;fixed-first-normal I1=x2*k[u][[x]]2",
    "graded_action": "M=[[-1,2u],[-3u/2,3u2-1]];det(M)=1 at every grade",
    "quotient": "chosen quadratic-section cubic residual is not invariant;M kills it",
    "lost_sidecar": "source determinant cocycle;constant-J does not descend",
    "raw_affine": "E1 unique p=3/4,q=9u/8;E2 unique r=-9/8,s=-27u/16;E3=135/16",
    "raw_candidate": "fixed seed after chi=x-3x2/4+9x3/8;u and x monically recover",
    "resummation": "chi=2(sqrt(1+3x)-1)/3;first omitted=-135x4/64;THM3846 nonsquare wall",
    "section_independent": "full polynomial quadratic fibre plus affine relative residuals is empty at E3",
    "nonintegral_hostile": "exceptional quadratic pole chart rejects every nonzero affine cubic correction",
    "boundary": "polynomial constant-J normal depth<=5 is separately closed by THM3856/3861/3867/3871",
    "status": "FINITE-EXACT/SCOUT;no Keller seed;no JC2 conclusion",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_RUSSELL_CUBIC_SEED_QUOTIENT_CONSTANT_DENSITY_SCOUT")
print("status=FINITE-EXACT+SCOUT;AFFINE_CUBIC_CELL_EMPTY;JC2_OPEN")
print("zero_jet_ideal=I0=x*k[u][[x]]^2;fixed_first_normal=I1=x^2*k[u][[x]]^2")
print("graded_cosmetic_action=M;det_M=1;rank_M=2;all_cubic_residuals_move")
print("quotient_verdict=NO_CUBIC_RESIDUAL_INVARIANT_WITHOUT_DETERMINANT_COCYCLE")
print("raw_affine_E1=p=3/4;q=9u/8")
print("raw_affine_E2=r=-9/8;s=-27u/16")
print("raw_affine_E3=135/16;constant_J=IMPOSSIBLE")
print("near_candidate=FIXED_SEED_PRECOMPOSED_BY_chi=x-3x^2/4+9x^3/8")
print("resummation=chi=2*(sqrt(1+3x)-1)/3;first_omitted=-135x^4/64")
print("resummation_scope=FORMAL_ALGEBRAIC_NOT_POLYNOMIAL_OR_RATIONAL;THM3846_SQUARE_WALL")
print("near_candidate_recovery=u_MONIC_CUBIC;chi_POLYNOMIAL;x_MONIC_CUBIC")
print("section_independent_affine_residual_cell=EMPTY_AT_E3")
print("known_nonintegral_chart=AFFINE_CUBIC_EFFECTIVITY_FORCES_r=s=0")
print("known_nonintegral_E1=3/2+8e^2-4eu-8e^2u^2;UNREPAIRED")
print("next=DEPTH6_OR_RATIONAL_DECAYING_RESIDUALS_OR_NONMONOGENIC_MULTI_CONTROL")
print(f"semantic_sha256={semantic_sha}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
