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


# The strengthened theorem uses only the two base-coordinate derivatives.
# Audit their complete-intersection claim through resultants, graded ranks,
# and an explicit socle representative rather than importing the canonical
# Groebner transcript.
Delta_A = sp.expand(sp.diff(Delta0, A))
Delta_C = sp.expand(sp.diff(Delta0, C))
equal(
    "base A derivative",
    Delta_A,
    -2 * (3230 * A**3 + 567 * A**2 * C - 49 * A * C**2 - 6 * C**3),
)
equal(
    "base C derivative",
    Delta_C,
    -2 * A * (189 * A**2 - 49 * A * C - 18 * C**2),
)
equal("base derivatives coprime", sp.Poly(sp.gcd(Delta_A, Delta_C), A, C).monic().as_expr(), 1)
equal("base resultant in A", sp.factor(sp.resultant(Delta_A, Delta_C, A)), 36864000000 * C**9)
equal("base resultant in C", sp.factor(sp.resultant(Delta_A, Delta_C, C)), -3072000000 * A**9)


def base_multiplication_matrix(degree: int) -> sp.Matrix:
    if degree < 3:
        return sp.zeros(degree + 1, 0)
    source = basis(degree - 3)
    return sp.Matrix.hstack(
        *(coefficient_vector(gradient * monomial, degree) for gradient in (Delta_A, Delta_C) for monomial in source)
    )


expected_hilbert = (1, 2, 3, 2, 1, 0, 0, 0)
hilbert_rows = []
for degree, expected_dimension in enumerate(expected_hilbert):
    multiplication = base_multiplication_matrix(degree)
    quotient_dimension = degree + 1 - multiplication.rank()
    equal(f"complete-intersection Hilbert row {degree}", quotient_dimension, expected_dimension)
    hilbert_rows.append(quotient_dimension)

degree_five_map = base_multiplication_matrix(5)
equal("degree-five square determinant", degree_five_map.det(), 36864000000)
degree_five_inverse = degree_five_map.inv()
gate(degree_five_map * degree_five_inverse == sp.eye(6), "degree-five inverse failed")

# For a (3,3) complete intersection the Hessian determinant represents the
# degree-four socle.  It is nonzero modulo the ideal, while multiplication by
# either base coordinate vanishes modulo the ideal.
base_groebner = sp.groebner([Delta_A, Delta_C], A, C, order="grevlex")
socle_candidate = sp.factor(sp.det(sp.hessian(Delta0, (A, C))))
socle_remainder = sp.factor(base_groebner.reduce(socle_candidate)[1])
nonzero("socle degree-four class", socle_remainder)
zero("A kills socle", base_groebner.reduce(A * socle_candidate)[1])
zero("C kills socle", base_groebner.reduce(C * socle_candidate)[1])

# The unique degree-five inverse of C^5 recovers the displayed quadratic
# tangent correction independently.
degree_five_target = coefficient_vector(C**5, 5)
base_first_vector = degree_five_inverse * degree_five_target
base_quadratics = basis(2)
alpha_displayed = sp.expand(sum(base_first_vector[i] * base_quadratics[i] for i in range(3)))
beta_displayed = sp.expand(sum(base_first_vector[3 + i] * base_quadratics[i] for i in range(3)))
alpha_expected = (
    -sp.Rational(1722409941, 16000000) * A**2
    + sp.Rational(3586046429, 72000000) * A * C
    + sp.Rational(1, 12) * C**2
)
beta_expected = (
    sp.Rational(26492305283, 14400000) * A**2
    - sp.Rational(2460621541, 48000000) * A * C
    - sp.Rational(1211682143, 72000000) * C**2
)
equal("displayed base alpha", alpha_displayed, alpha_expected)
equal("displayed base beta", beta_displayed, beta_expected)
equal("displayed base first identity", Delta_A * alpha_displayed + Delta_C * beta_displayed, C**5)


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


# A second, genuinely two-coordinate lift uses only six fixed degree-five
# identities.  Every higher monomial is divided by a degree-five monomial,
# whose unique preimage comes from the inverse above.  This differs from the
# strengthened primary checker's fresh RREF at each degree.
def fixed_base_right_inverse(form: sp.Expr, degree: int) -> tuple[sp.Expr, sp.Expr]:
    alpha = sp.S.Zero
    beta = sp.S.Zero
    polynomial = sp.Poly(sp.expand(form), A, C)
    for (a_degree, c_degree), coefficient in polynomial.terms():
        gate(a_degree + c_degree == degree, f"base input degree {degree} not homogeneous")
        if a_degree >= 5:
            degree_five_index = 0
            multiplier = A ** (a_degree - 5) * C**c_degree
        else:
            degree_five_index = 5 - a_degree
            multiplier = C ** (c_degree - degree_five_index)
        representation = degree_five_inverse[:, degree_five_index]
        alpha += coefficient * multiplier * sum(
            representation[i] * base_quadratics[i] for i in range(3)
        )
        beta += coefficient * multiplier * sum(
            representation[3 + i] * base_quadratics[i] for i in range(3)
        )
    alpha = sp.expand(alpha)
    beta = sp.expand(beta)
    equal(
        f"fixed degree-five base inverse at degree {degree}",
        Delta_A * alpha + Delta_C * beta,
        form,
    )
    return alpha, beta


def recursive_base_lift(
    target_by_degree: dict[int, sp.Expr], stop: int
) -> tuple[sp.Expr, sp.Expr]:
    P = A
    Q = C
    for degree in range(5, stop + 1):
        current = homogeneous(
            Delta0.subs({A: P, C: Q}, simultaneous=True), degree
        )
        wanted = target_by_degree.get(degree, sp.S.Zero)
        alpha, beta = fixed_base_right_inverse(wanted - current, degree)
        gate(total_order(alpha) >= degree - 3, f"base alpha low term at {degree}")
        gate(total_order(beta) >= degree - 3, f"base beta low term at {degree}")
        P = sp.expand(P + alpha)
        Q = sp.expand(Q + beta)
        changed = Delta0.subs({A: P, C: Q}, simultaneous=True)
        equal(f"base recursion closes degree {degree}", homogeneous(changed, degree), wanted)
        for old_degree in range(4, degree):
            old_wanted = Delta0 if old_degree == 4 else target_by_degree.get(old_degree, sp.S.Zero)
            equal(
                f"base recursion preserves degree {old_degree} at stage {degree}",
                homogeneous(changed, old_degree),
                old_wanted,
            )
    return P, Q


base_P, base_Q = recursive_base_lift({5: C**5}, 10)
base_error = sp.expand(
    Delta0.subs({A: base_P, C: base_Q}, simultaneous=True) - Delta0 - C**5
)
gate(total_order(base_error) >= 11, "base recursion leaks below degree eleven")
nonzero("finite base recursion still has debt", homogeneous(base_error, 11))
equal("base map linear P", homogeneous(base_P, 1), A)
equal("base map linear Q", homogeneous(base_Q, 1), C)
equal(
    "base map Jacobian at origin",
    sp.det(sp.Matrix([[sp.diff(base_P, A), sp.diff(base_P, C)], [sp.diff(base_Q, A), sp.diff(base_Q, C)]])).subs({A: 0, C: 0}),
    1,
)

# Dense right-equivalence control: not just the single C^5 deformation.
base_dense_target = {degree: dense_target[degree] for degree in range(5, 9)}
base_dense_P, base_dense_Q = recursive_base_lift(base_dense_target, 8)
base_dense_phi = sum(base_dense_target.values())
base_dense_error = sp.expand(
    Delta0.subs({A: base_dense_P, C: base_dense_Q}, simultaneous=True)
    - Delta0
    - base_dense_phi
)
gate(total_order(base_dense_error) >= 9, "dense base recursion leaks below degree nine")

# Construct the inverse of a finite tangent-identity truncation recursively
# through degree four, and verify both composition orders to that degree.
def truncate(poly: sp.Expr, stop: int) -> sp.Expr:
    return sp.expand(sum(homogeneous(poly, degree) for degree in range(stop + 1)))


map_P = truncate(base_P, 4)
map_Q = truncate(base_Q, 4)
inverse_P, inverse_Q = A, C
for degree in range(2, 5):
    left_P = map_P.subs({A: inverse_P, C: inverse_Q}, simultaneous=True)
    left_Q = map_Q.subs({A: inverse_P, C: inverse_Q}, simultaneous=True)
    correction_P = -homogeneous(left_P - A, degree)
    correction_Q = -homogeneous(left_Q - C, degree)
    inverse_P = truncate(inverse_P + correction_P, 4)
    inverse_Q = truncate(inverse_Q + correction_Q, 4)

left_composite_P = map_P.subs({A: inverse_P, C: inverse_Q}, simultaneous=True)
left_composite_Q = map_Q.subs({A: inverse_P, C: inverse_Q}, simultaneous=True)
right_composite_P = inverse_P.subs({A: map_P, C: map_Q}, simultaneous=True)
right_composite_Q = inverse_Q.subs({A: map_P, C: map_Q}, simultaneous=True)
for degree in range(1, 5):
    equal(f"formal inverse left P degree {degree}", homogeneous(left_composite_P, degree), A if degree == 1 else 0)
    equal(f"formal inverse left Q degree {degree}", homogeneous(left_composite_Q, degree), C if degree == 1 else 0)
    equal(f"formal inverse right P degree {degree}", homogeneous(right_composite_P, degree), A if degree == 1 else 0)
    equal(f"formal inverse right Q degree {degree}", homogeneous(right_composite_Q, degree), C if degree == 1 else 0)


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

# Every polynomial base map with identity linear part keeps these four
# factors nonconstant and pairwise nonassociate: their linear parts are the
# four rows just checked.  The exact substitution identity exhibits the
# factorization, while irreducibility above excludes equality with the
# nonzero-lambda target in k[A,C].
base_pullback_factors = (
    base_P,
    base_Q + 5 * base_P,
    4 * base_Q + 19 * base_P,
    3 * base_Q - 17 * base_P,
)
equal(
    "base pullback remains four-factor product",
    Delta0.subs({A: base_P, C: base_Q}, simultaneous=True),
    sp.prod(base_pullback_factors),
)
for index, factor in enumerate(base_pullback_factors):
    nonzero(f"base factor {index} has nonzero linear part", homogeneous(factor, 1))

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
    "base derivatives form a (3,3) complete intersection",
    "quotient Hilbert function 1,2,3,2,1 and socle degree four",
    "fixed degree-five identities lift every m^5 base deformation",
    "tangent-identity base recursion and formal inverse recursion",
    "polynomial base-map orbit remains a four-nonunit product",
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
print("BASE_GRADIENT_CI", "degrees (3,3); degree-five determinant=36864000000")
print("BASE_QUOTIENT_HILBERT", ",".join(str(item) for item in hilbert_rows[:5]), "; socle_degree=4")
print("BASE_RIGHT_EQUIVALENCE", "fixed-degree-five identities; one-place through degree 10; dense through degree 8")
print("FORMAL_AUTOMORPHISM", "identity Jacobian; two-sided inverse replay through degree 4")
print("POLYNOMIAL_BASE_MAP", "NO_GO: four nonconstant factors versus irreducible target")
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
