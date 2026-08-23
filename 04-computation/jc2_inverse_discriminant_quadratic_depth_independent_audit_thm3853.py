#!/usr/bin/env python3
"""Independent hostile audit of THM-3853.

The canonical checker inserts the closed discriminant formula, eliminates the
target scalar through the C^5 row, and saturates it with an auxiliary variable
using Buchberger/grevlex.  This checker instead:

* derives the discriminant as a resultant;
* reverses the quadratic basis and coefficient-group order;
* uses degree-stratified scaling to normalize every nonzero target scalar to 1;
* substitutes the forced C^2 coefficient before elimination; and
* computes the remaining two unit ideals only after this dimension-reducing
  normalization, with reversed generator and variable order.

It also checks the two normalizations, the excluded reducible orientation,
lambda=0, and an active near miss.  No tracked file is modified.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(newline="\n")

import sympy as sp


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    require(sp.factor(left - right) == 0, label)


A, C, T, Z, MU = sp.symbols("A C T Z MU")
aa, bb, cc, dd, x = sp.symbols("aa bb cc dd x")


# Derive, rather than insert, the binary-cubic discriminant.
generic_cubic = aa * x**3 + bb * x**2 + cc * x + dd
resultant_quotient = sp.cancel(
    -sp.resultant(generic_cubic, sp.diff(generic_cubic, x), x) / aa
)
require(sp.denom(resultant_quotient) == 1,
        "resultant quotient is a polynomial")
generic_discriminant = sp.expand(resultant_quotient)
expected_discriminant = sp.expand(
    bb**2 * cc**2 - 4 * aa * cc**3 - 4 * bb**3 * dd
    - 27 * aa**2 * dd**2 + 18 * aa * bb * cc * dd
)
same(generic_discriminant, expected_discriminant,
     "resultant rederives the cubic discriminant")


def discriminant(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        generic_discriminant.subs({aa: a, bb: b, cc: c, dd: d})
    )


linear_packet = (A, C, 7 * A, -3 * A)
delta_zero = sp.factor(discriminant(*linear_packet))
expected_delta_zero = sp.factor(
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
)
same(delta_zero, expected_delta_zero, "base packet discriminant")


# Geometry in complementary (M,L) coordinates.  The gcd check supplies the
# four simple finite roots; the degree pair (5,4) supplies one infinity place.
D_C = sp.factor(delta_zero.subs({A: T, C: 1}))
D_SUM = sp.factor(delta_zero.subs({A: T, C: 1 - T}))
same(D_C, -T * (5 * T + 1) * (17 * T - 3) * (19 * T + 4),
     "C-orientation quartic")
same(D_SUM, -T * (4 * T + 1) * (15 * T + 4) * (20 * T - 3),
     "A+C-orientation quartic")
for label, polynomial in (("C", D_C), ("SUM", D_SUM)):
    require(sp.degree(polynomial, T) == 4, f"{label} quartic degree")
    require(sp.gcd(polynomial, sp.diff(polynomial, T)) == 1,
            f"{label} has four simple roots")
    require(sp.LC(sp.Poly(polynomial, T)) != 0,
            f"{label} integral normalization equation is genuinely monicable")

LAMBDA = sp.symbols("LAMBDA", nonzero=True)
delta_C_target = sp.expand(delta_zero + LAMBDA * C**5)
delta_SUM_target = sp.expand(delta_zero + LAMBDA * (A + C) ** 5)

# Localize at L and set T=M/L.  The localized equation is linear in L, while
# the displayed parametrization gives the birational finite normalization.
same(
    sp.expand(delta_C_target.subs({A: T * C}) / C**4),
    D_C + LAMBDA * C,
    "C localization equation",
)
same(
    delta_C_target.subs({A: -T * D_C / LAMBDA,
                         C: -D_C / LAMBDA}),
    0,
    "C normalization parametrization",
)
same(
    sp.expand(
        delta_SUM_target.subs({A: T * Z, C: (1 - T) * Z}) / Z**4
    ),
    D_SUM + LAMBDA * Z,
    "sum localization equation",
)
same(
    delta_SUM_target.subs({A: -T * D_SUM / LAMBDA,
                           C: -(1 - T) * D_SUM / LAMBDA}),
    0,
    "sum normalization parametrization",
)
require(sp.degree(-D_C / LAMBDA, T) == 4
        and sp.degree(-T * D_C / LAMBDA, T) == 5,
        "C normalization has degree pair L/M=(4,5)")
require(sp.degree(-D_SUM / LAMBDA, T) == 4
        and sp.degree(-T * D_SUM / LAMBDA, T) == 5
        and sp.degree(-(1 - T) * D_SUM / LAMBDA, T) == 5,
        "sum normalization has one degree-four line and degree-five coordinates")
require(delta_C_target.subs(C, 0) != 0,
        "C does not divide its target")
require(delta_SUM_target.subs({A: T, C: -T}) != 0,
        "A+C does not divide its target")

# The tempting third orientation is genuinely reducible for every lambda.
delta_A_target = sp.factor(delta_zero + LAMBDA * A**5)
require(sp.rem(delta_A_target, A, A) == 0,
        "L=A target retains the base factor A")
require(sp.factor(delta_A_target / A).is_polynomial(A, C),
        "L=A quotient is polynomial")


# Reverse both the monomial basis and coefficient-group order relative to the
# canonical checker.  v0,v1,v2 belong to q_d; the last triple belongs to q_a.
v = sp.symbols("v0:12")
quadratic_basis = (C**2, A * C, A**2)
qd = sum(v[index] * quadratic_basis[index] for index in range(3))
qc = sum(v[3 + index] * quadratic_basis[index] for index in range(3))
qb = sum(v[6 + index] * quadratic_basis[index] for index in range(3))
qa = sum(v[9 + index] * quadratic_basis[index] for index in range(3))
perturbed_packet = (A + qa, C + qb, 7 * A + qc, -3 * A + qd)
delta = discriminant(*perturbed_packet)
error = sp.Poly(sp.expand(delta - delta_zero), A, C)

for coefficient in perturbed_packet:
    same(coefficient.subs({A: 0, C: 0}), 0,
         "all perturbed coefficients remain in the origin maximal ideal")

# The base degree d coefficient is homogeneous of degree d-4 in the twelve
# perturbation parameters.  This gives an exact alternative to saturation:
# if q realizes nonzero lambda, then q/lambda realizes lambda=1.
degree_rows: dict[int, list[sp.Expr]] = {}
for degree in range(9):
    degree_rows[degree] = []
    for a_power in range(degree + 1):
        coefficient = sp.expand(
            error.coeff_monomial(A**a_power * C ** (degree - a_power))
        )
        if coefficient == 0:
            continue
        degree_rows[degree].append(coefficient)
        parameter_poly = sp.Poly(coefficient, *v)
        require(
            all(sum(monomial) == degree - 4
                for monomial, _ in parameter_poly.terms()),
            "degree row has the predicted parameter homogeneity",
        )
require(all(not degree_rows[degree] for degree in range(5)),
        "linear packet freezes all degrees below five")
require(all(degree_rows[degree] for degree in range(5, 9)),
        "every possible error degree five through eight is represented")

# A direct symbolic spot check of the scaling law, now derived coefficient by
# coefficient rather than inferred from the quartic formula.
for degree in range(5, 9):
    for coefficient in degree_rows[degree]:
        scaled = sp.expand(coefficient.subs({variable: MU * variable for variable in v}))
        same(scaled, MU ** (degree - 4) * coefficient,
             "quadratic-layer scaling law")


def normalized_target_system(line: sp.Expr) -> tuple[list[sp.Expr], sp.Expr]:
    """Set lambda=1 and eliminate the forced q_d(C^2) coefficient."""

    target = sp.Poly(sp.expand(line**5), A, C)
    c5_error = error.coeff_monomial(C**5)
    same(c5_error, -4 * v[0], "C^5 row reads the target scalar")
    forced = sp.Rational(-1, 4)
    equations: list[sp.Expr] = []
    for degree in range(5, 9):
        for a_power in range(degree + 1):
            monomial = A**a_power * C ** (degree - a_power)
            wanted = target.coeff_monomial(monomial) if degree == 5 else 0
            equation = sp.expand(
                (error.coeff_monomial(monomial) - wanted).subs(v[0], forced)
            )
            if equation != 0:
                equations.append(equation)
    return equations, forced


systems: dict[str, list[sp.Expr]] = {}
for label, line in (("C", C), ("SUM", A + C)):
    equations, forced_v0 = normalized_target_system(line)
    require(forced_v0 == -sp.Rational(1, 4),
            f"{label} lambda-one normalization")
    require(len(equations) == 29, f"{label} has 29 remaining equations")
    systems[label] = list(reversed(equations))

# Lambda zero must not be smuggled into the contradiction.
zero_parameters = {variable: 0 for variable in v}
for degree in range(5, 9):
    for coefficient in degree_rows[degree]:
        same(coefficient.subs(zero_parameters), 0,
             "lambda-zero unperturbed control survives")

# Active near miss: q_d=C^2 gives scalar -4 but does not clear the debt.
near_miss = dict(zero_parameters)
near_miss[v[0]] = 1
near_error = sp.factor((delta - delta_zero).subs(near_miss))
expected_near_error = sp.factor(
    -4 * C**5 + 126 * A**2 * C**3 + 162 * A**3 * C**2
    - 27 * A**2 * C**4
)
same(near_error, expected_near_error, "active q_d=C^2 debt")
require(sp.factor(near_error + 4 * C**5) != 0,
        "active near miss has non-target residual terms")


# Different exact elimination route: lambda has been normalized to one,
# v0 has been substituted, the equation order and variable order are reversed,
# leaving a different eleven-variable presentation from the canonical
# saturated thirteen-variable presentation.
remaining_variables = tuple(reversed(v[1:]))
unit_bases: dict[str, list[sp.Expr]] = {}
for label in ("C", "SUM"):
    basis = sp.groebner(systems[label], *remaining_variables, order="grevlex")
    unit_bases[label] = list(basis)
    require(unit_bases[label] == [sp.Integer(1)],
            f"{label} normalized inverse ideal is the unit ideal")


canonical_script = Path(
    "04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py"
)
canonical_output = Path(
    "05-knowledge/results/jc2_inverse_discriminant_quadratic_depth_thm3853.out"
)
canonical_script_hash = hashlib.sha256(canonical_script.read_bytes()).hexdigest()
canonical_output_hash = hashlib.sha256(canonical_output.read_bytes()).hexdigest()
require(
    canonical_script_hash
    == "2e5b4a240a4f81c793adc5cda1a3edc5da24056d5e49324e4617a1745be84bed",
    "canonical script hash matches theorem metadata",
)
require(
    canonical_output_hash
    == "f74d5652d24fa1b881d6da4ca0ac3cdcc5da144b42d5011e8c14b27f20c9c0a7",
    "canonical output hash matches theorem metadata",
)

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "independent checker has no optimization-sensitive assert")

semantic = {
    "classification": "PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED",
    "field": "algebraically closed characteristic zero",
    "linear_packet": ["A", "C", "7A", "-3A"],
    "quadratic_universe": "all four coefficients plus arbitrary homogeneous degree-two forms",
    "targets": ["Delta0+lambda*C^5", "Delta0+lambda*(A+C)^5"],
    "target_quantifier": "each target and every nonzero lambda",
    "geometry": "irreducible; affine normalization A1; four simple origins; one infinity place",
    "alternative_algebra": "normalize lambda=1 by graded scaling; force qd_C2=-1/4; eleven-variable unit ideals",
    "hostile_boundaries": ["lambda=0 survives", "L=A reducible", "q_d=C^2 leaves debt"],
    "scope": "quadratic-depth inverse-discriminant obstruction only; no cubic-cover or JC2 conclusion",
}
semantic_blob = json.dumps(semantic, sort_keys=True,
                           separators=(",", ":")).encode()

print("experiment=THM3853-independent-hostile-audit")
print("classification=PROVED;VERIFIED-EXACT;INDEPENDENTLY-HOSTILE-AUDITED")
print("discriminant_route=resultant_not_inserted_formula")
print("geometry=C_and_AplusC_affine_normalization_A1;four_simple_origins;one_infinity_place")
print("reducible_boundary=L_equals_A")
print("scalar_route=graded_scaling_normalizes_every_nonzero_lambda_to_one")
print("forced_normalized_qd_C2=-1/4")
print(f"C_equations={len(systems['C'])};C_normalized_basis={unit_bases['C']}")
print(f"SUM_equations={len(systems['SUM'])};SUM_normalized_basis={unit_bases['SUM']}")
print("hostile_controls=lambda_zero_survives;q_d_C2_has_residual_debt")
print("quantifiers=both_named_lines;every_nonzero_lambda;all_12_homogeneous_quadratic_parameters")
print("scope=no_quadratic_depth_realization_only;higher_degree_formal_other_packets_and_JC2_open")
print(f"canonical_script_sha256={canonical_script_hash}")
print(f"canonical_output_sha256={canonical_output_hash}")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
