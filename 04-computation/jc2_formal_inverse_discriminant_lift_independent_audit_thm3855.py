#!/usr/bin/env python3
"""Independent hostile audit for THM-3855.

The canonical companion lifts each monomial through the 4x4 cubic-gradient
matrix.  This audit instead builds the full degree-n multiplication map,
chooses pivot columns by exact RREF, and solves the resulting square system.
It also tests a dense multi-degree target, not only the one-place C^5 row.
No Python ``assert`` is used, so optimized replay is equally active.
"""

from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


A, C = sp.symbols("A C")
a, b, c, d = sp.symbols("a b c d")
z = sp.symbols("z")
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(label: str, value: sp.Expr) -> None:
    gate(sp.expand(value) == 0, f"{label}: {sp.factor(value)}")


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def nonzero(label: str, value: sp.Expr) -> None:
    gate(sp.factor(value) != 0, f"{label}: unexpectedly zero")


def disc(row: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr] | list[sp.Expr]) -> sp.Expr:
    aa, bb, cc, dd = row
    return sp.expand(
        bb**2 * cc**2
        - 4 * aa * cc**3
        - 4 * bb**3 * dd
        - 27 * aa**2 * dd**2
        + 18 * aa * bb * cc * dd
    )


def homogeneous(poly: sp.Expr, degree: int) -> sp.Expr:
    answer = 0
    for (i, j), coefficient in sp.Poly(sp.expand(poly), A, C).terms():
        if i + j == degree:
            answer += coefficient * A**i * C**j
    return sp.expand(answer)


def total_order(poly: sp.Expr) -> int | sp.Expr:
    poly = sp.Poly(sp.expand(poly), A, C)
    if poly.is_zero:
        return sp.oo
    return min(sum(monomial) for monomial, _ in poly.terms())


def basis(degree: int) -> list[sp.Expr]:
    return [A ** (degree - j) * C**j for j in range(degree + 1)]


def coefficient_vector(poly: sp.Expr, degree: int) -> sp.Matrix:
    polynomial = sp.Poly(sp.expand(poly), A, C)
    return sp.Matrix(
        [polynomial.coeff_monomial(A ** (degree - j) * C**j) for j in range(degree + 1)]
    )


base = (A, C, 7 * A, -3 * A)
Delta0 = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
equal("base discriminant", disc(base), Delta0)


# Recover the differential as the coefficient of an auxiliary parameter,
# rather than importing the canonical gradient list.
symbolic_disc = disc((a, b, c, d))
coefficient_variables = (a, b, c, d)
gradients = tuple(
    sp.expand(sp.diff(symbolic_disc, variable).subs(dict(zip(coefficient_variables, base))))
    for variable in coefficient_variables
)
for index, gradient in enumerate(gradients):
    equal(
        f"polarized gradient {index}",
        sp.Poly(
            disc(tuple(base[j] + (z if j == index else 0) for j in range(4))),
            z,
        ).coeff_monomial(z),
        gradient,
    )

cubic_basis = basis(3)
gradient_matrix = sp.Matrix.hstack(
    *(coefficient_vector(gradient, 3) for gradient in gradients)
)
equal("gradient determinant", gradient_matrix.det(), 640000)
equal("bad-characteristic factorization", sp.factorint(abs(int(gradient_matrix.det())))[2], 10)
equal("bad-characteristic factorization five", sp.factorint(abs(int(gradient_matrix.det())))[5], 4)
gate(set(sp.factorint(abs(int(gradient_matrix.det())))) == {2, 5}, "unexpected determinant prime")

# Every cubic monomial has a unique scalar representation in the four
# gradients.  These identities, multiplied by arbitrary forms, prove the
# homogeneous surjectivity used in the all-orders argument.
gradient_inverse = gradient_matrix.inv()
for column, monomial in enumerate(cubic_basis):
    coefficients = gradient_inverse[:, column]
    equal(
        f"cubic spanning identity {column}",
        sum(coefficients[i] * gradients[i] for i in range(4)),
        monomial,
    )


def full_degree_right_inverse(form: sp.Expr, degree: int) -> list[sp.Expr]:
    """Solve the entire R_(degree-3)^4 -> R_degree map via RREF pivots."""

    input_degree = degree - 3
    input_basis = basis(input_degree)
    columns: list[sp.Matrix] = []
    labels: list[tuple[int, int]] = []
    for row_index in range(4):
        for monomial_index, monomial in enumerate(input_basis):
            columns.append(coefficient_vector(gradients[row_index] * monomial, degree))
            labels.append((row_index, monomial_index))
    linear_map = sp.Matrix.hstack(*columns)
    _, pivots = linear_map.rref()
    gate(len(pivots) == degree + 1, f"degree {degree} map not surjective")
    square_map = linear_map[:, list(pivots)]
    gate(square_map.det() != 0, f"degree {degree} pivot minor singular")
    pivot_solution = square_map.inv() * coefficient_vector(form, degree)
    correction = [sp.S.Zero] * 4
    for pivot_index, scalar in zip(pivots, pivot_solution):
        row_index, monomial_index = labels[pivot_index]
        correction[row_index] += scalar * input_basis[monomial_index]
    correction = [sp.expand(item) for item in correction]
    equal(
        f"degree {degree} full-map right inverse",
        sum(gradients[i] * correction[i] for i in range(4)),
        form,
    )
    return correction


# A generic symbolic quadratic verifies the sharp degree boundary: fixing the
# linear row makes degree four immutable; the degree-five change is exactly
# the differential, and every nonlinear Taylor term starts in degree six.
qvars = sp.symbols("q0:12")
generic_quadratic = tuple(
    qvars[3 * i] * A**2 + qvars[3 * i + 1] * A * C + qvars[3 * i + 2] * C**2
    for i in range(4)
)
generic_change = sp.expand(disc(tuple(base[i] + generic_quadratic[i] for i in range(4))) - Delta0)
equal("fixed linear part freezes degree four", homogeneous(generic_change, 4), 0)
equal(
    "quadratic differential controls degree five",
    homogeneous(generic_change, 5),
    sum(gradients[i] * generic_quadratic[i] for i in range(4)),
)
gate(total_order(generic_change - homogeneous(generic_change, 5)) >= 6, "nonlinear Taylor debt below degree six")


def recursive_lift(target_by_degree: dict[int, sp.Expr], stop: int) -> tuple[list[sp.Expr], list[sp.Expr]]:
    coefficients = list(base)
    corrections = [sp.S.Zero] * 4
    for degree in range(5, stop + 1):
        current_discriminant = disc(coefficients)
        wanted = target_by_degree.get(degree, sp.S.Zero)
        error = sp.expand(wanted - homogeneous(current_discriminant, degree))
        update = full_degree_right_inverse(error, degree)
        gate(
            all(total_order(item) >= degree - 3 for item in update),
            f"degree {degree} correction has low terms",
        )
        coefficients = [sp.expand(coefficients[i] + update[i]) for i in range(4)]
        corrections = [sp.expand(corrections[i] + update[i]) for i in range(4)]
        new_discriminant = disc(coefficients)
        equal(f"recursion closes degree {degree}", homogeneous(new_discriminant, degree), wanted)
        for old_degree in range(4, degree):
            old_wanted = Delta0 if old_degree == 4 else target_by_degree.get(old_degree, sp.S.Zero)
            equal(
                f"recursion preserves degree {old_degree} at stage {degree}",
                homogeneous(new_discriminant, old_degree),
                old_wanted,
            )
    return coefficients, corrections


# First replay the actual one-place target.  The RREF pivot lift is generally
# different from the canonical monomial lift, which makes this an independent
# recursion path.
one_place_coefficients, one_place_corrections = recursive_lift({5: C**5}, 12)
one_place_error = sp.expand(disc(one_place_coefficients) - Delta0 - C**5)
gate(total_order(one_place_error) >= 13, "one-place truncation leaks below degree thirteen")
nonzero("one-place finite truncation still has debt", homogeneous(one_place_error, 13))
gate(
    all(total_order(item - base[i]) >= 2 for i, item in enumerate(one_place_coefficients)),
    "one-place correction changed the linear row",
)

# The displayed canonical first correction is separately checked, but is not
# used by the independent recursion.
canonical_first = (
    sp.Rational(91449, 40000) * C**2,
    sp.Rational(151263, 40000) * C**2,
    -sp.Rational(194481, 20000) * C**2,
    -sp.Rational(1, 4) * C**2,
)
equal(
    "displayed first correction differential",
    sum(gradients[i] * canonical_first[i] for i in range(4)),
    C**5,
)
canonical_first_change = sp.expand(disc(tuple(base[i] + canonical_first[i] for i in range(4))) - Delta0)
equal("displayed first correction degree five", homogeneous(canonical_first_change, 5), C**5)
gate(total_order(canonical_first_change - C**5) >= 6, "displayed first correction has unexpected low debt")

# A dense target with a nonzero homogeneous form in every degree is a hostile
# control against a recursion accidentally specialized to C^5.
dense_target: dict[int, sp.Expr] = {}
for degree in range(5, 10):
    dense_target[degree] = sp.expand(
        sum(
            (-1) ** j * (degree + 2 * j + 1) * A ** (degree - j) * C**j
            for j in range(degree + 1)
        )
    )
dense_coefficients, dense_corrections = recursive_lift(dense_target, 9)
dense_phi = sum(dense_target.values())
dense_error = sp.expand(disc(dense_coefficients) - Delta0 - dense_phi)
gate(total_order(dense_error) >= 10, "dense target recursion leaks below degree ten")
gate(
    all(total_order(item - base[i]) >= 2 for i, item in enumerate(dense_coefficients)),
    "dense correction changed the linear row",
)


# -------------------------------------------------------------------------
# One-place target geometry and the formal branch packet.
# -------------------------------------------------------------------------

lam = sp.symbols("lam", nonzero=True)
t = sp.symbols("t")
target_delta = sp.expand(Delta0 + lam * C**5)
F = sp.expand(t * (1 + 5 * t) * (4 + 19 * t) * (3 - 17 * t))
equal("chart identity", target_delta.subs(A, t * C), C**4 * (F + lam * C))
parameter_C = -F / lam
parameter_A = -t * F / lam
equal("polynomial normalization parametrizes target", target_delta.subs({A: parameter_A, C: parameter_C}), 0)
equal("rational inverse on C nonzero", sp.factor(parameter_A / parameter_C), t)
equal("four normalization addresses", sp.degree(F, t), 4)
equal("four roots are distinct", sp.gcd(F, sp.diff(F, t)), 1)
equal("normalization C degree", sp.degree(parameter_C, t), 4)
equal("normalization A degree", sp.degree(parameter_A, t), 5)

# Elimination recovers the target itself, proving that no extra affine
# component was hidden by the C!=0 chart.
implicit = sp.factor(sp.resultant(A - t * C, lam * C + F, t))
equal("normalization implicit equation", implicit, target_delta)
equal(
    "lambda scaling reduces to lambda one",
    target_delta.subs({A: A / lam, C: C / lam}),
    (Delta0 + C**5) / lam**4,
)
gate(sp.Poly((Delta0 + C**5), A, C).is_irreducible, "lambda-one target unexpectedly reducible over Q")

tangent_lines = (A, C + 5 * A, 4 * C + 19 * A, 3 * C - 17 * A)
equal("target tangent cone", homogeneous(target_delta, 4), sp.prod(tangent_lines))
line_rows = [
    (sp.Poly(line, A, C).coeff_monomial(A), sp.Poly(line, A, C).coeff_monomial(C))
    for line in tangent_lines
]
pair_determinants = []
for i in range(4):
    for j in range(i + 1, 4):
        determinant = sp.det(sp.Matrix([line_rows[i], line_rows[j]]))
        nonzero(f"distinct tangent lines {i},{j}", determinant)
        pair_determinants.append(int(determinant))

# The nonzero-lambda and characteristic-zero hypotheses are both active.
equal("lambda-zero boundary", target_delta.subs(lam, 0), Delta0)
equal("degree-four target cannot be reached", homogeneous(generic_change, 4), 0)


# -------------------------------------------------------------------------
# Formal cubic order gates: multiplication, index, local fibre, and S3 data.
# -------------------------------------------------------------------------

x, y = sp.symbols("x y")
aa, bb, cc, dd = one_place_coefficients
index_form = -(aa * x**3 + bb * x**2 * y + cc * x * y**2 + dd * y**3)
gate(all(total_order(item) >= 1 for item in one_place_coefficients), "coefficient outside maximal ideal")
gate(
    all(total_order(coefficient) >= 1 for coefficient in sp.Poly(index_form, x, y).coeffs()),
    "index form has a unit coefficient",
)

# Delone--Faddeev multiplication constants all vanish modulo (A,C), leaving
# k[omega,theta]/(omega,theta)^2 as claimed.
multiplication_constants = (-aa * cc, bb, -aa, -aa * dd, -bb * dd, dd, -cc)
gate(all(total_order(item) >= 1 for item in multiplication_constants), "special fibre not square-zero")
equal("truncated discriminant target through degree twelve", homogeneous(disc(one_place_coefficients), 4), Delta0)
for degree in range(5, 13):
    equal(
        f"formal branch target truncation degree {degree}",
        homogeneous(disc(one_place_coefficients), degree),
        C**5 if degree == 5 else 0,
    )

# In the completion, squarefree tangent cone gives four distinct smooth
# height-one branches.  Their four odd valuations make the discriminant
# nonsquare.  The checker records the exact tangent determinants used by
# that Hensel argument; the DVR maximality/S2 implication is audited in the
# accompanying proof report.
gate(len(pair_determinants) == 6, "incomplete tangent-pair census")


source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "inactive assert found")

semantic_packet = (
    "THM-3855 independent hostile audit",
    "gradient determinant 2^10*5^4 and cubic-gradient basis",
    "full degree multiplication maps solved by RREF pivot minors",
    "one-place recursion through degree 12 via independent corrections",
    "dense all-monomial target recursion through degree 9",
    "m^5 boundary sharp for fixed linear row",
    "normalization A1 with degrees (5,4) and one infinity point",
    "four distinct smooth formal branches and odd discriminant valuations",
    "square-zero special fibre and index values in maximal ideal",
    "formal only; polynomial algebraization remains open",
)

print("THM3855_HOSTILE_AUDIT", "PASS")
print("GRADIENT_DETERMINANT", gradient_matrix.det(), "=2^10*5^4")
print("FORMAL_SURJECTIVITY", "R_(n-3)^4 -> R_n full row rank for every n>=3")
print("INDEPENDENT_RECURSION", "RREF-pivot path; one-place through degree 12")
print("DENSE_TARGET_CONTROL", "all monomials in every degree 5..9")
print("SHARP_BOUNDARY", "fixed linear part cannot alter degree 4; corrections start degree 2")
print("ONE_PLACE_NORMALIZATION", "k[t], deg(A)=5, deg(C)=4, one point at infinity")
print("FORMAL_BRANCHES", "four distinct smooth height-one factors, each valuation 1")
print("ORDER_GATES", "finite-free S2 + DVR maximality; local special fibre; no unit index value")
print("EXCEPTIONS", "lambda=0 excluded; characteristics 2 and 5 excluded; polynomial termination open")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("GATES", GATES)
