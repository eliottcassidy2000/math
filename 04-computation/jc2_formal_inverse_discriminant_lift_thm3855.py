#!/usr/bin/env python3
"""Exact companion for THM-3855's formal inverse-discriminant lift."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


A, C = sp.symbols("A C")
a, b, c, d = sp.symbols("a b c d")


def discriminant(aa: sp.Expr, bb: sp.Expr, cc: sp.Expr, dd: sp.Expr) -> sp.Expr:
    return sp.expand(
        bb**2 * cc**2
        - 4 * aa * cc**3
        - 4 * bb**3 * dd
        - 27 * aa**2 * dd**2
        + 18 * aa * bb * cc * dd
    )


def homogeneous(expression: sp.Expr, degree: int) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), A, C)
    return sp.Add(
        *[
            coefficient * A**i * C**j
            for (i, j), coefficient in polynomial.terms()
            if i + j == degree
        ]
    )


base = (A, C, 7 * A, -3 * A)
Delta0 = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
zero(discriminant(*base) - Delta0, "four-ray discriminant")

D_symbolic = discriminant(a, b, c, d)
gradients = [
    sp.expand(sp.diff(D_symbolic, variable).subs(dict(zip((a, b, c, d), base))))
    for variable in (a, b, c, d)
]
expected_gradients = [
    -2 * A**2 * (929 * A + 189 * C),
    -2 * A * (189 * A**2 - 49 * A * C - 18 * C**2),
    -2 * A * (294 * A**2 + 27 * A * C - 7 * C**2),
    2 * (81 * A**3 + 63 * A**2 * C - 2 * C**3),
]
for index, (actual, expected) in enumerate(zip(gradients, expected_gradients)):
    zero(actual - expected, f"coefficient gradient {index}")

cubic_monomials = (A**3, A**2 * C, A * C**2, C**3)
gradient_matrix = sp.Matrix(
    [
        [sp.Poly(gradient, A, C).coeff_monomial(monomial) for gradient in gradients]
        for monomial in cubic_monomials
    ]
)
expected_matrix = sp.Matrix(
    [
        [-1858, -378, -588, 162],
        [-378, 98, -54, 126],
        [0, 36, 14, 0],
        [0, 0, 0, -4],
    ]
)
gate(gradient_matrix == expected_matrix, "gradient coefficient matrix")
gate(gradient_matrix.det() == 640000, "gradient determinant")
gradient_inverse = gradient_matrix.inv()
gate(gradient_matrix * gradient_inverse == sp.eye(4), "gradient right inverse")

groebner = sp.groebner(gradients, A, C, order="grevlex")
gate(
    [sp.Poly(polynomial, A, C).monic().as_expr() for polynomial in groebner.polys]
    == [A**3, A**2 * C, A * C**2, C**3],
    "gradient ideal is the cube of the maximal ideal",
)

# Stronger mechanism: the target is already formally right-equivalent to the
# four-line discriminant by a tangent-identity change of the base variables.
Delta_A = sp.expand(sp.diff(Delta0, A))
Delta_C = sp.expand(sp.diff(Delta0, C))
expected_Delta_A = -2 * (
    3230 * A**3 + 567 * A**2 * C - 49 * A * C**2 - 6 * C**3
)
expected_Delta_C = -2 * A * (189 * A**2 - 49 * A * C - 18 * C**2)
zero(Delta_A - expected_Delta_A, "base A-gradient")
zero(Delta_C - expected_Delta_C, "base C-gradient")
gate(sp.Poly(sp.gcd(Delta_A, Delta_C), A, C).monic().as_expr() == 1,
     "base gradients are coprime")
gate(
    sp.factor(sp.resultant(Delta_A, Delta_C, A)) == 36864000000 * C**9,
    "base-gradient resultant",
)

quadratic_monomials = (A**2, A * C, C**2)
quintic_monomials = tuple(A**i * C ** (5 - i) for i in range(6))
base_columns = [
    sp.expand(gradient * monomial)
    for gradient in (Delta_A, Delta_C)
    for monomial in quadratic_monomials
]
base_degree_five_matrix = sp.Matrix(
    [
        [sp.Poly(column, A, C).coeff_monomial(monomial) for column in base_columns]
        for monomial in quintic_monomials
    ]
)
gate(base_degree_five_matrix.det() == -36864000000,
     "base-gradient degree-five determinant")
base_groebner = sp.groebner([Delta_A, Delta_C], A, C, order="grevlex")
for monomial in quintic_monomials:
    zero(base_groebner.reduce(monomial)[1], "base Jacobian ideal contains m^5")

alpha_2 = (
    -sp.Rational(1722409941, 16000000) * A**2
    + sp.Rational(3586046429, 72000000) * A * C
    + sp.Rational(1, 12) * C**2
)
beta_2 = (
    sp.Rational(26492305283, 14400000) * A**2
    - sp.Rational(2460621541, 48000000) * A * C
    - sp.Rational(1211682143, 72000000) * C**2
)
zero(Delta_A * alpha_2 + Delta_C * beta_2 - C**5,
     "displayed tangent-identity quadratic correction")


def base_right_inverse(form: sp.Expr, degree: int) -> tuple[sp.Expr, sp.Expr]:
    """Deterministic lift of a homogeneous form through (Delta_A,Delta_C)."""

    source_degree = degree - 3
    source_monomials = tuple(
        A**i * C ** (source_degree - i) for i in range(source_degree + 1)
    )
    target_monomials = tuple(A**i * C ** (degree - i) for i in range(degree + 1))
    columns = [
        sp.expand(gradient * monomial)
        for gradient in (Delta_A, Delta_C)
        for monomial in source_monomials
    ]
    matrix = sp.Matrix(
        [
            [sp.Poly(column, A, C).coeff_monomial(monomial) for column in columns]
            for monomial in target_monomials
        ]
    )
    vector = sp.Matrix(
        [sp.Poly(sp.expand(form), A, C).coeff_monomial(monomial)
         for monomial in target_monomials]
    )
    _, pivots = matrix.rref()
    gate(len(pivots) == degree + 1, "base right-inverse full row rank")
    pivot_matrix = matrix[:, list(pivots)]
    pivot_solution = pivot_matrix.inv() * vector
    solution = [sp.S.Zero] * len(columns)
    for pivot, value in zip(pivots, pivot_solution):
        solution[pivot] = value
    width = len(source_monomials)
    alpha = sp.expand(sum(solution[i] * source_monomials[i] for i in range(width)))
    beta = sp.expand(
        sum(solution[width + i] * source_monomials[i] for i in range(width))
    )
    zero(Delta_A * alpha + Delta_C * beta - form,
         "base homogeneous right inverse")
    return alpha, beta


# A finite replay of the formal coordinate recursion.  The first step is
# unique and equals the displayed quadratic pair above.
P, Q = A, C
first_base_correction = None
for total_degree in range(5, 13):
    current = homogeneous(Delta0.subs({A: P, C: Q}, simultaneous=True), total_degree)
    wanted = C**5 if total_degree == 5 else sp.S.Zero
    correction = base_right_inverse(wanted - current, total_degree)
    P = sp.expand(P + correction[0])
    Q = sp.expand(Q + correction[1])
    zero(
        homogeneous(Delta0.subs({A: P, C: Q}, simultaneous=True), total_degree)
        - wanted,
        f"formal right-equivalence recursion through degree {total_degree}",
    )
    if total_degree == 5:
        first_base_correction = correction
gate(first_base_correction == (alpha_2, beta_2),
     "first base correction is the displayed unique solution")
base_error = sp.expand(Delta0.subs({A: P, C: Q}, simultaneous=True) - (Delta0 + C**5))
gate(
    all(homogeneous(base_error, degree) == 0 for degree in range(0, 13)),
    "truncated right-equivalence has no error through degree twelve",
)
gate(homogeneous(base_error, 13) != 0,
     "right-equivalence degree-thirteen hostile remains after truncation")

# The coefficient germ itself is formally rigid under tangent base changes
# and SL2 changes of the binary variables.  We use the infinitesimal action
# e=X*d/dY, f=Y*d/dX, h=X*d/dX-Y*d/dY.
X, Y = sp.symbols("X Y")
binary_form = A * X**3 + C * X**2 * Y + 7 * A * X * Y**2 - 3 * A * Y**3


def binary_coefficient_vector(form: sp.Expr) -> sp.Matrix:
    polynomial = sp.Poly(sp.expand(form), X, Y)
    return sp.Matrix(
        [polynomial.coeff_monomial(monomial)
         for monomial in (X**3, X**2 * Y, X * Y**2, Y**3)]
    )


v_A = sp.Matrix([1, 0, 7, -3])
v_C = sp.Matrix([0, 1, 0, 0])
sl2_e = binary_coefficient_vector(X * sp.diff(binary_form, Y))
sl2_f = binary_coefficient_vector(Y * sp.diff(binary_form, X))
sl2_h = binary_coefficient_vector(
    X * sp.diff(binary_form, X) - Y * sp.diff(binary_form, Y)
)
gate(sl2_e == sp.Matrix([C, 14 * A, -9 * A, 0]), "SL2 e action")
gate(sl2_f == sp.Matrix([0, 3 * A, 2 * C, 7 * A]), "SL2 f action")
gate(sl2_h == sp.Matrix([3 * A, C, -7 * A, 9 * A]), "SL2 h action")


def coefficient_quotient(vector: sp.Matrix) -> sp.Matrix:
    """Quotient coefficient space by the two base-coordinate directions."""

    return sp.Matrix([vector[2] - 7 * vector[0], vector[3] + 3 * vector[0]])


gate(coefficient_quotient(v_A) == sp.zeros(2, 1), "A base direction kernel")
gate(coefficient_quotient(v_C) == sp.zeros(2, 1), "C base direction kernel")
gauge_matrix = sp.Matrix.hstack(
    coefficient_quotient(sl2_e),
    coefficient_quotient(sl2_f),
    coefficient_quotient(sl2_h),
)
expected_gauge_matrix = sp.Matrix(
    [[-9 * A - 7 * C, 2 * C, -28 * A], [3 * C, 7 * A, 18 * A]]
)
gate(gauge_matrix == expected_gauge_matrix, "quotient SL2 action matrix")

gauge_minors = [
    sp.expand(gauge_matrix[:, [left, right]].det())
    for left, right in ((0, 1), (0, 2), (1, 2))
]
expected_gauge_minors = [
    -63 * A**2 - 49 * A * C - 6 * C**2,
    -162 * A**2 - 42 * A * C,
    196 * A**2 + 36 * A * C,
]
gate(gauge_minors == expected_gauge_minors, "quotient SL2 maximal minors")
minor_coefficient_matrix = sp.Matrix(
    [
        [sp.Poly(minor, A, C).coeff_monomial(monomial) for minor in gauge_minors]
        for monomial in (A**2, A * C, C**2)
    ]
)
gate(minor_coefficient_matrix.det() == -14400,
     "quotient SL2 minors span the square maximal ideal")
minor_groebner = sp.groebner(gauge_minors, A, C, order="grevlex")
gate(
    [sp.Poly(polynomial, A, C).monic().as_expr()
     for polynomial in minor_groebner.polys]
    == [A**2, A * C, C**2],
    "quotient SL2 Fitting ideal is (A,C)^2",
)

# Finite homogeneous controls for the all-degree Fitting-annihilator proof.
for total_degree in range(2, 9):
    source_monomials = tuple(
        A**i * C ** (total_degree - 1 - i) for i in range(total_degree)
    )
    target_monomials = tuple(
        A**i * C ** (total_degree - i) for i in range(total_degree + 1)
    )
    columns = [
        sp.expand(gauge_matrix[:, column] * monomial)
        for column in range(3)
        for monomial in source_monomials
    ]
    homogeneous_gauge_matrix = sp.Matrix(
        [
            [
                sp.Poly(column[component], A, C).coeff_monomial(monomial)
                for column in columns
            ]
            for component in range(2)
            for monomial in target_monomials
        ]
    )
    gate(
        homogeneous_gauge_matrix.rank() == 2 * (total_degree + 1),
        f"quotient SL2 action surjective in degree {total_degree}",
    )


def right_inverse(form: sp.Expr, degree: int) -> list[sp.Expr]:
    """Lift a homogeneous degree-`degree` form through the cubic gradients."""

    corrections = [sp.S.Zero] * 4
    for (i, j), coefficient in sp.Poly(sp.expand(form), A, C).terms():
        gate(i + j == degree, "right-inverse input is homogeneous")
        if i >= 3:
            basis_index = 0
            multiplier = A ** (i - 3) * C**j
        elif i == 2:
            basis_index = 1
            multiplier = C ** (j - 1)
        elif i == 1:
            basis_index = 2
            multiplier = C ** (j - 2)
        else:
            basis_index = 3
            multiplier = C ** (j - 3)
        for column in range(4):
            corrections[column] += (
                coefficient * gradient_inverse[column, basis_index] * multiplier
            )
    zero(
        sum(gradient * correction for gradient, correction in zip(gradients, corrections))
        - form,
        "homogeneous gradient right inverse",
    )
    return [sp.expand(correction) for correction in corrections]


# The one-place target of THM-3853.  The recursion is performed at lambda=1;
# homogeneity gives the displayed scalar first correction for arbitrary lambda.
target = sp.expand(Delta0 + C**5)
coefficients = list(base)
first_correction = None
for total_degree in range(5, 13):
    current = homogeneous(discriminant(*coefficients), total_degree)
    wanted = C**5 if total_degree == 5 else sp.S.Zero
    correction = right_inverse(wanted - current, total_degree)
    gate(
        all(
            homogeneous(entry, total_degree - 3) == entry
            for entry in correction
        ),
        f"correction degree {total_degree - 3}",
    )
    coefficients = [
        sp.expand(entry + update) for entry, update in zip(coefficients, correction)
    ]
    zero(
        homogeneous(discriminant(*coefficients), total_degree) - wanted,
        f"inverse-discriminant recursion through degree {total_degree}",
    )
    if total_degree == 5:
        first_correction = correction

expected_first = [
    sp.Rational(91449, 40000) * C**2,
    sp.Rational(151263, 40000) * C**2,
    -sp.Rational(194481, 20000) * C**2,
    -sp.Rational(1, 4) * C**2,
]
gate(first_correction == expected_first, "explicit first quadratic correction")

error = sp.expand(discriminant(*coefficients) - target)
gate(
    all(homogeneous(error, degree) == 0 for degree in range(0, 13)),
    "truncated lift has no error through degree twelve",
)
gate(homogeneous(error, 13) != 0, "degree-thirteen hostile remains after truncation")
gate(
    all(entry.subs({A: 0, C: 0}) == 0 for entry in coefficients),
    "all lifted index coefficients remain in the maximal ideal",
)

semantic = {
    "base_packet": "(A,C,7A,-3A)",
    "gradient_matrix_det": 640000,
    "gradient_ideal": "(A,C)^3",
    "base_gradient_resultant": "36864000000*C^9",
    "base_jacobian_ideal": "complete intersection containing (A,C)^5",
    "right_equivalence": "tangent-identity formal base automorphism",
    "coefficient_rigidity": "formal base automorphism times SL2 gauge",
    "gauge_fitting_ideal": "(A,C)^2",
    "gauge_minor_det": -14400,
    "target": "Delta0+C^5",
    "formal_lift": "all m-adic orders; corrections begin in degree 2",
    "finite_replay": "degrees 5 through 12; first residual degree 13",
    "index_gate": "all coefficients remain in (A,C), so no unit value",
    "scope": "formal completion only; polynomial algebraization and connected Keller surface open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate")
print("base=(A,C,7A,-3A);gradient_det=640000;gradient_ideal=(A,C)^3")
print("base_gradient_resultant=36864000000*C^9;base_jacobian_CI_contains=(A,C)^5")
print("right_equivalence=tangent_identity_formal_base_automorphism")
print("coefficient_rigidity=formal_base_change_times_SL2_gauge")
print("gauge_quotient_minors_generate=(A,C)^2;minor_coefficient_det=-14400")
print(
    "base_first_correction="
    "alpha2=(-1722409941/16000000)A^2+(3586046429/72000000)AC+(1/12)C^2;"
    "beta2=(26492305283/14400000)A^2-(2460621541/48000000)AC-"
    "(1211682143/72000000)C^2"
)
print("target=Delta0+C^5;formal_lift=ALL_ORDERS")
print("first_correction=C^2*(91449/40000,151263/40000,-194481/20000,-1/4)")
print("finite_replay=through_degree_12;first_truncation_residual_degree=13")
print("index=coefficients_remain_in_(A,C);unit_value=NONE")
print("scope=formal_only;polynomial_algebraization_and_connected_Keller_surface_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
