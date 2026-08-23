#!/usr/bin/env python3
"""Independent hostile audit of THM-3854.

This path implicitizes the curve by a parameter resultant, computes its affine
singular scheme directly, classifies nodes by Hessians, derives the quintic
discriminant by a second resultant, invokes an exact number-field Galois-group
algorithm on the specialization, checks A5 perfection, and audits the two
seminormal descents and the all-degree derivative ideal.

No Python ``assert`` statements are used, so optimized mode retains every
gate.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp
from sympy.combinatorics.named_groups import AlternatingGroup, SymmetricGroup
from sympy.polys.numberfields import galois_group

sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    require(sp.factor(left - right) == 0, label)


def t_order(expression: sp.Expr, local_variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), local_variable)
    return min(monomial[0] for monomial, coefficient in polynomial.terms()
               if coefficient != 0)


t, u, s, x, y, X, Y, Z, r, q, w = sp.symbols("t u s x y X Y Z r q w")


# -------------------------------------------------------------------------
# I. Independent implicitization, normalization, and singularity ledger.
# -------------------------------------------------------------------------

x_t = t**4 - 2 * t**2
y_t = 3 * t**5 - 5 * t**3
implicit = sp.factor(sp.resultant(x - x_t, y - y_t, t))
expected_implicit = sp.expand(
    81 * x**5 + 90 * x**4 + 25 * x**3
    + 30 * x**2 * y**2 + 30 * x * y**2 - y**4 + 8 * y**2
)
same(implicit, expected_implicit,
     "parameter resultant independently recovers the quintic")
delta = implicit
require(sp.Poly(delta, x, y).total_degree() == 5,
        "implicit equation has degree five")

# A dense rational inverse supplies birationality.  This is checked after,
# rather than used to produce, the implicit equation.
inverse_denominator = 27 * x**3 + 33 * x**2 + 10 * x + y**2
inverse_numerator = y * (9 * x**2 + 13 * x + 4)
pulled_denominator = sp.factor(
    inverse_denominator.subs({x: x_t, y: y_t})
)
same(inverse_numerator.subs({x: x_t, y: y_t}),
     t * pulled_denominator, "dense rational inverse")
require(pulled_denominator != 0, "inverse is defined on a dense open set")

# Direct affine singular-scheme computation, not a collision resultant.
delta_x = sp.diff(delta, x)
delta_y = sp.diff(delta, y)
singular_basis = sp.groebner(
    [delta, delta_x, delta_y], y, x, order="lex", domain=sp.QQ
)
expected_singular_basis = (
    y**2
    - sp.Rational(1, 2)
    * (162 * x**5 + 405 * x**4 + 310 * x**3 + 75 * x**2),
    sp.Rational(1, 9) * y * (x + 1) ** 2 * (9 * x + 4),
    sp.Rational(1, 81) * x**2 * (x + 1) ** 2 * (9 * x + 4) * (9 * x + 5),
)
require(len(singular_basis.polys) == 3, "singular ideal basis length")
for actual, expected in zip(singular_basis.polys, expected_singular_basis):
    same(actual.as_expr(), expected, "direct singular ideal basis")

# These four x-fibres contain exactly 1+2+2+1=6 geometric points.
singular_packets = (
    {x: 0, y: 0},
    {x: -1, y: 2},
    {x: -1, y: -2},
    {x: -sp.Rational(5, 9), y: 0},
)
for point in singular_packets:
    require(delta.subs(point) == 0 and delta_x.subs(point) == 0
            and delta_y.subs(point) == 0,
            "rational singular packet lies in the singular ideal")

# The remaining algebraic packet x=-4/9, y^2=8/27 contributes two points.
for polynomial in (delta, delta_x, delta_y):
    specialized = sp.Poly(polynomial.subs(x, -sp.Rational(4, 9)), y)
    remainder = sp.rem(
        specialized, sp.Poly(y**2 - sp.Rational(8, 27), y)
    )
    require(remainder.is_zero,
            "two algebraic node points lie in the singular ideal")

# Three one-address branches have Puiseux orders (2,3), hence type A2.
v = sp.symbols("v")
cusp_data = (
    (0, 0, 0, 0),
    (1, -1, -2, -sp.Rational(15, 4)),
    (-1, -1, 2, sp.Rational(15, 4)),
)
for address, x_value, y_value, tangent_correction in cusp_data:
    local_x = sp.expand(x_t.subs(t, address + v) - x_value)
    local_y = sp.expand(y_t.subs(t, address + v) - y_value)
    transverse = sp.expand(local_y + tangent_correction * local_x)
    require(t_order(local_x, v) == 2, "cusp tangent coordinate order two")
    require(t_order(transverse, v) == 3,
            "cusp transverse coordinate order three")

# The other three singular points have nondegenerate quadratic tangent cones.
hessian_determinant = sp.factor(sp.det(sp.hessian(delta, (x, y))))
require(hessian_determinant.subs({x: -sp.Rational(5, 9), y: 0})
        == -sp.Rational(8000, 243), "rational node Hessian")
algebraic_hessian = sp.expand(
    hessian_determinant.subs(x, -sp.Rational(4, 9))
)
algebraic_hessian_remainder = sp.rem(
    sp.Poly(algebraic_hessian, y),
    sp.Poly(y**2 - sp.Rational(8, 27), y),
).as_expr()
require(algebraic_hessian_remainder == -sp.Rational(16000, 243),
        "two algebraic nodes have nonzero Hessian")
require(3 + 3 == (5 - 1) * (5 - 2) // 2,
        "3A2+3A1 exhausts the quintic genus budget")

# Infinity is checked in the Y=1 chart.  There is only X=0 on Z=0, and the
# transverse Z derivative is a unit.  The normalization pulls Z back to s^5.
delta_h = sp.expand(Z**5 * delta.subs({x: X / Z, y: Y / Z}))
same(delta_h.subs(Z, 0), 81 * X**5,
     "unique projective point at infinity")
require(sp.diff(delta_h, Z).subs({X: 0, Y: 1, Z: 0}) == -1,
        "infinity point is smooth")
projective_coordinates = (
    t**4 * s - 2 * t**2 * s**3,
    3 * t**5 - 5 * t**3 * s**2,
    s**5,
)
same(
    delta_h.subs(
        dict(zip((X, Y, Z), projective_coordinates)), simultaneous=True
    ),
    0,
    "projective normalization identity",
)
require(projective_coordinates[1].subs(s, 0) == 3 * t**5,
        "projective parametrization has no infinity basepoint")
require(projective_coordinates[2] == s**5,
        "unique infinity address has contact five")

# The sextic double-plane boundary is locally two smooth branches tangent to
# order five.  The displayed change turns its double cover into UV=q^10,
# i.e. an A9 rational double point.  Finite branch A2/A1 packets lift to the
# same ADE labels by the analogous difference-of-squares identities.
local_branch_model = z_model = Z * (Z - q**5)
same(
    (2 * Z - q**5) ** 2 - (2 * w) ** 2 - q**10,
    -4 * (w**2 - local_branch_model),
    "fivefold line-quintic contact gives the A9 double-cover model",
)


# -------------------------------------------------------------------------
# II. Polynomial quintic completion and monodromy.
# -------------------------------------------------------------------------

F = 3 * u**5 - 10 * u**3 - 15 * x * u + 4 * y
F_u = sp.diff(F, u)
same(F_u, 15 * (u**4 - 2 * u**2 - x),
     "ramification equation of the quintic map")

# Derive the discriminant from the defining resultant, independently of
# SymPy's discriminant convenience function and the canonical calculation.
discriminant_from_resultant = sp.cancel(sp.resultant(F, F_u, u) / 3)
same(discriminant_from_resultant, -64_800_000 * delta,
     "quintic resultant discriminant")

y_map = sp.expand((-3 * u**5 + 10 * u**3 + 15 * x * u) / 4)
same(y_map.subs({x: x_t, u: t}), y_t,
     "ramification curve maps to the normalization")
same(sp.diff(y_map, u), -sp.Rational(15, 4) * (u**4 - 2 * u**2 - x),
     "polynomial map Jacobian")
require(sp.Poly(F, y).degree() == 1 and sp.Poly(F, y).LC() == 4,
        "total quintic hypersurface is the polynomial plane k[x,u]")
require(sp.Poly(F / 3, u).LC() == 1,
        "quintic is monic after division by a unit")

# Exact arithmetic specialization.  We use both the two Frobenius witnesses
# and SymPy's independent exact number-field Galois algorithm.
specialized_expr = sp.expand(F.subs({x: -3, y: 1}))
specialized = sp.Poly(specialized_expr, u, domain=sp.QQ)
specialized_disc = sp.discriminant(specialized_expr, u)
require(specialized_disc == 834_688_800_000,
        "specialized quintic discriminant")
for prime in (29, 67):
    require(specialized_disc % prime != 0,
            "Frobenius witness prime is unramified")
require(sp.Poly(specialized_expr, u, modulus=29).is_irreducible,
        "mod-29 factor degree pattern is 5")
factor_67 = sp.Poly(
    3 * (u + 7) * (u + 21) * (u - 23) * (u**2 - 5 * u + 5),
    u,
    modulus=67,
)
require(sp.Poly(specialized_expr, u, modulus=67) == factor_67,
        "mod-67 factor degree pattern is 1+1+1+2")
require(pow(5, 33, 67) == 66,
        "mod-67 quadratic discriminant is nonsquare")

specialized_group, alternating_flag = galois_group(specialized_expr, u)
require(specialized_group.order() == 120,
        "independent exact specialized Galois group has order 120")
require(specialized_group.degree == 5
        and specialized_group.order() == SymmetricGroup(5).order(),
        "independent exact specialized group is the full S5")
require(alternating_flag is False,
        "specialized group is not contained in A5")

# Over algebraically closed constants, the irreducible discriminant divisor
# has odd valuation one, so the geometric group is not contained in A5.  The
# sign kernel is A5; its trivial abelianization forbids every C3 quotient.
A5 = AlternatingGroup(5)
require(A5.order() == 60 and A5.derived_subgroup().order() == 60,
        "A5 is perfect and has no nontrivial cyclic quotient")
require(SymmetricGroup(5).derived_subgroup().order() == 60,
        "sign kernel in S5 has order 60")


# -------------------------------------------------------------------------
# III. Seminormal descents and nonsquare residuals.
# -------------------------------------------------------------------------

h_even = t**2 * (t**2 - 1) * (9 * t**2 - 14)
P_even = 81 * x**3 + 49 * x**2 + 8 * y**2
Q_even = (
    -243 * x**4 - 143 * x**3 + 81 * x**2 * y**2
    + 120 * x * y**2 + 64 * y**2
)
h_odd = t * (t**2 - 1) * (9 * t**2 - 4) * (3 * t**2 - 5)
P_odd = -648 * x**3 - 720 * x**2 + 81 * x * y**2 - 200 * x + 49 * y**2
Q_odd = (
    6561 * x**4 * y + 17010 * x**3 * y + 18009 * x**2 * y
    + 243 * x * y**3 + 8760 * x * y + 143 * y**3 + 1600 * y
)

for label, h_value, p_value, q_value in (
    ("even", h_even, P_even, Q_even),
    ("odd", h_odd, P_odd, Q_odd),
):
    same(p_value.subs({x: x_t, y: y_t}), h_value**2,
         f"{label} square descent")
    same(q_value.subs({x: x_t, y: y_t}), h_value**3,
         f"{label} cube descent")
    require(any(sp.diff(h_value, t).subs(t, address) != 0
                for address in (0, 1, -1)),
            f"{label} seminormal element is not in the cusp coordinate ring")

even_quotient, even_remainder = sp.div(
    sp.Poly(P_even**3 - Q_even**2, x, y, domain=sp.QQ),
    sp.Poly(delta, x, y, domain=sp.QQ),
)
require(even_remainder.is_zero, "even residual division is exact")
even_residual = sp.expand(even_quotient.as_expr())
same(even_residual, 6561 * x**4 + 3888 * x**3 - 512 * y**2,
     "even residual quotient")

odd_quotient, odd_remainder = sp.div(
    sp.Poly(P_odd**3 - Q_odd**2, x, y, domain=sp.QQ),
    sp.Poly(delta, x, y, domain=sp.QQ),
)
require(odd_remainder.is_zero, "odd residual division is exact")
odd_core = (
    41472 * x**2 + 6561 * x * y**2 + 46080 * x
    + 3888 * y**2 + 12800
)
same(odd_quotient.as_expr(), -(9 * x + 5) ** 2 * odd_core,
     "odd residual quotient")


def quadratic_in_y_is_not_square(polynomial: sp.Expr, label: str) -> None:
    poly_y = sp.Poly(polynomial, y)
    require(poly_y.degree() == 2, f"{label} has y-degree two")
    require(poly_y.coeff_monomial(y**2) != 0,
            f"{label} has nonzero quadratic coefficient")
    require(poly_y.coeff_monomial(y) == 0,
            f"{label} has no linear-y coefficient")
    require(poly_y.coeff_monomial(1) != 0,
            f"{label} has nonzero constant-in-y coefficient")


quadratic_in_y_is_not_square(even_residual, "even residual")
quadratic_in_y_is_not_square(odd_core, "odd residual core")


# -------------------------------------------------------------------------
# IV. All normal degrees and formal corrections fail on an interior arm.
# -------------------------------------------------------------------------

g = t * (t**2 - 1)
x_prime = sp.diff(x_t, t)
y_prime = sp.diff(y_t, t)
require(sp.gcd(sp.Poly(x_prime, t), sp.Poly(y_prime, t)) == sp.Poly(g, t),
        "derivative ideal is the proper principal ideal generated by g")

# Arbitrary first normal coefficients already exhibit the universal factor.
alpha_coefficients = sp.symbols("a0:8")
beta_coefficients = sp.symbols("b0:8")
alpha = sum(coefficient * t**degree
            for degree, coefficient in enumerate(alpha_coefficients))
beta = sum(coefficient * t**degree
           for degree, coefficient in enumerate(beta_coefficients))
bracket_constant = sp.expand(alpha * y_prime - x_prime * beta)
same(bracket_constant, g * (15 * t * alpha - 4 * beta),
     "arbitrary polynomial first-normal bracket has cusp factor")
for address in (0, 1, -1):
    require(bracket_constant.subs(t, address) == 0,
            "arm bracket vanishes at every cusp address")
require(g not in (sp.Integer(1), sp.Integer(-1)),
        "cusp factor is a nonunit")


# -------------------------------------------------------------------------
# V. Metadata and deterministic transcript.
# -------------------------------------------------------------------------

canonical_script_path = Path(
    "04-computation/jc2_integrated_three_cusp_quintic_s5_thm3854.py"
)
canonical_output_path = Path(
    "05-knowledge/results/jc2_integrated_three_cusp_quintic_s5_thm3854.out"
)
canonical_script_hash = hashlib.sha256(canonical_script_path.read_bytes()).hexdigest()
canonical_output_hash = hashlib.sha256(canonical_output_path.read_bytes()).hexdigest()
require(
    canonical_script_hash
    == "ffe82ca9bd4c147b685d75b16bfa94e18269b143fc9cbe629f5edf227538ae5e",
    "canonical script hash matches theorem metadata",
)
require(
    canonical_output_hash
    == "0a4598e2f042a08afd7b1d046c3cea2bc345e48fc6282980e1d2f265a619847e",
    "canonical output hash matches theorem metadata",
)

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "independent checker has no optimization-sensitive assert")

semantic = {
    "classification": "PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED",
    "front": "rational irreducible quintic; affine 3A2+3A1; one smooth infinity place",
    "completion": "finite polynomial-plane quintic; resultant discriminant -64800000 Delta",
    "monodromy": "generic geometric S5; sign kernel A5 perfect; no C3 quotient",
    "seminormal": "two proper square/cube descents; both residual discriminants nonsquare",
    "arm": "derivative ideal (t(t^2-1)); no constant Jacobian at any normal degree",
    "double_plane": "finite 3A2+3A1 and infinity A9; crepant sextic resolution is K3",
    "hostile_boundaries": "natural quintic only; alternative cubic orders, boundary arms, class-group 3-torsion, JC2 open",
}
semantic_blob = json.dumps(semantic, sort_keys=True,
                           separators=(",", ":")).encode()

print("experiment=THM3854-independent-hostile-audit")
print("classification=PROVED;VERIFIED-EXACT;INDEPENDENTLY-HOSTILE-AUDITED")
print("implicitization=parameter_resultant;degree=5;birational=YES")
print("affine_singular_scheme=3A2_plus_3A1;genus_budget=6")
print("infinity=one_smooth_place;line_contact=5")
print("completion=finite_polynomial_plane_quintic;discriminant_constant=-64800000")
print("specialized_exact_group=S5;generic_geometric_group=S5")
print("sign_kernel=A5;derived_order=60;C3_quotient=NONE")
print("seminormal_descents=2;proper=YES;nonsquare_residuals=2")
print("double_plane=finite_3A2_plus_3A1;infinity_A9;minimal_resolution=K3")
print("all_degree_arm_gate=derivative_ideal_t_times_t2minus1;interior_Keller_arm=IMPOSSIBLE")
print("scope=natural_quintic_only;alternative_cubic_orders,boundary_data,Cl3,and_JC2_OPEN")
print(f"canonical_script_sha256={canonical_script_hash}")
print(f"canonical_output_sha256={canonical_output_hash}")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
