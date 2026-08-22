#!/usr/bin/env python3
"""Finite exact small-rational scout in short JC(2) Gauss charts.

The search is deliberately narrow and completely stated.  It does not claim
bounded emptiness outside the printed universe.  Determinant one is built in;
both row curls are solved over QQ before any collision test.  Every surviving
map is then given a polynomial inverse, and an exact rational point bank is a
redundant collision hostile.  A second finite control tests the first
correction equation for 1,000 signed-monomial five-word base triples.

Run with both

    python3 04-computation/jc2_small_rational_four_gauss_search.py
    python3 -O 04-computation/jc2_small_rational_four_gauss_search.py

The program never uses ``assert``.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

from flint import fmpq, fmpq_mat


Q = Fraction
Monomial = tuple[int, int]
Poly = dict[Monomial, Q]


def require(label: str, condition: bool, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(f"FAIL: {label}")
    checks.append(label)


def clean(poly: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in poly.items()
            if coefficient}


def add(left: Poly, right: Poly, scale: Q = Q(1)) -> Poly:
    answer = dict(left)
    for monomial, coefficient in right.items():
        new_coefficient = answer.get(monomial, Q(0)) + scale * coefficient
        if new_coefficient:
            answer[monomial] = new_coefficient
        else:
            answer.pop(monomial, None)
    return answer


def scalar(poly: Poly, value: Q) -> Poly:
    return clean({monomial: value * coefficient
                  for monomial, coefficient in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    answer: Poly = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            monomial = (i + k, j + ell)
            answer[monomial] = answer.get(monomial, Q(0)) + a * b
    return clean(answer)


def power(poly: Poly, exponent: int) -> Poly:
    answer = {(0, 0): Q(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        n //= 2
    return answer


def derivative(poly: Poly, variable: int) -> Poly:
    answer: Poly = {}
    for (i, j), coefficient in poly.items():
        exponent = (i, j)[variable]
        if not exponent:
            continue
        monomial = (i - 1, j) if variable == 0 else (i, j - 1)
        answer[monomial] = coefficient * exponent
    return answer


def derivative_monomial(monomial: Monomial, variable: int) -> Poly:
    exponent = monomial[variable]
    if not exponent:
        return {}
    lowered = ((monomial[0] - 1, monomial[1]) if variable == 0
               else (monomial[0], monomial[1] - 1))
    return {lowered: Q(exponent)}


def degree(poly: Poly) -> int:
    if not poly:
        return -1
    return max(i + j for i, j in poly)


def evaluate(poly: Poly, x_value: Q, y_value: Q) -> Q:
    return sum(coefficient * x_value**i * y_value**j
               for (i, j), coefficient in poly.items())


def monomial_basis(maximum_degree: int) -> list[Monomial]:
    return [(i, total - i)
            for total in range(maximum_degree + 1)
            for i in range(total + 1)]


def as_fmpq(value: Q) -> fmpq:
    return fmpq(value.numerator, value.denominator)


def as_fraction(value: fmpq) -> Q:
    return Q(int(value.numerator), int(value.denominator))


def rational_system_rows(
    operator_columns: list[Poly],
    right_hand_side: Poly,
    forced_zero_columns: tuple[int, ...] = (),
) -> list[list[Q]]:
    """Return the augmented coefficient rows for a polynomial identity."""

    monomials = sorted(
        set(right_hand_side).union(*(set(column) for column in operator_columns))
    )
    number_of_unknowns = len(operator_columns)
    rows = [
        [column.get(monomial, Q(0)) for column in operator_columns]
        + [right_hand_side.get(monomial, Q(0))]
        for monomial in monomials
    ]
    for column_index in forced_zero_columns:
        rows.append(
            [Q(int(index == column_index)) for index in range(number_of_unknowns)]
            + [Q(0)]
        )
    return rows


def modular_rank(rows: list[list[Q]], number_of_columns: int, prime: int) -> int:
    """Compute a matrix rank over one prime field, using exact reductions."""

    matrix = [
        [
            (value.numerator * pow(value.denominator, -1, prime)) % prime
            for value in row[:number_of_columns]
        ]
        for row in rows
    ]
    pivot_row = 0
    for column in range(number_of_columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [value * inverse % prime
                             for value in matrix[pivot_row]]
        for row in range(pivot_row + 1, len(matrix)):
            multiplier = matrix[row][column]
            if multiplier:
                matrix[row] = [
                    (left - multiplier * right) % prime
                    for left, right in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def full_rank_modular_certificate(
    operator_columns: list[Poly],
    right_hand_side: Poly,
    forced_zero_columns: tuple[int, ...] = (),
    prime: int = 1009,
) -> str | None:
    """Soundly certify zero kernel or inconsistency from full modular rank.

    Full coefficient-column rank modulo ``prime`` forces full rank over QQ.
    If the augmented modular rank is one larger, inconsistency over QQ follows.
    For a homogeneous system, full column rank proves that zero is its only
    rational solution.  No conclusion is drawn from a rank drop.
    """

    rows = rational_system_rows(
        operator_columns, right_hand_side, forced_zero_columns
    )
    number_of_unknowns = len(operator_columns)
    coefficient_rank = modular_rank(rows, number_of_unknowns, prime)
    if coefficient_rank != number_of_unknowns:
        return None
    if not right_hand_side:
        return "zero-only"
    augmented_rank = modular_rank(rows, number_of_unknowns + 1, prime)
    if augmented_rank == number_of_unknowns + 1:
        return "inconsistent"
    return None


def affine_solution_space(
    operator_columns: list[Poly],
    right_hand_side: Poly,
    forced_zero_columns: tuple[int, ...] = (),
) -> tuple[Poly, list[Poly]] | None:
    """Solve an exact polynomial-coefficient linear system over QQ.

    Return one particular solution and a basis for the homogeneous solution
    space, each encoded in the coefficient-column basis.  ``None`` means the
    system is inconsistent.
    """

    number_of_unknowns = len(operator_columns)
    rational_rows = rational_system_rows(
        operator_columns, right_hand_side, forced_zero_columns
    )
    reduced, _ = fmpq_mat(
        [[as_fmpq(value) for value in row] for row in rational_rows]
    ).rref()

    pivot_rows: list[tuple[int, int]] = []
    for row in range(reduced.nrows()):
        pivot = next(
            (column for column in range(number_of_unknowns)
             if reduced[row, column] != 0),
            None,
        )
        if pivot is None:
            if reduced[row, number_of_unknowns] != 0:
                return None
        else:
            pivot_rows.append((row, pivot))

    pivot_columns = {pivot for _, pivot in pivot_rows}
    free_columns = [column for column in range(number_of_unknowns)
                    if column not in pivot_columns]

    particular_vector = [Q(0) for _ in range(number_of_unknowns)]
    for row, pivot in pivot_rows:
        particular_vector[pivot] = as_fraction(reduced[row, number_of_unknowns])

    homogeneous_vectors: list[list[Q]] = []
    for free in free_columns:
        vector = [Q(0) for _ in range(number_of_unknowns)]
        vector[free] = Q(1)
        for row, pivot in pivot_rows:
            vector[pivot] = -as_fraction(reduced[row, free])
        homogeneous_vectors.append(vector)

    coefficient_labels = [{(index, 0): Q(1)}
                          for index in range(number_of_unknowns)]

    def vector_as_label_poly(vector: list[Q]) -> Poly:
        return clean({(index, 0): value for index, value in enumerate(vector)})

    # ``coefficient_labels`` exists only to make the encoding convention
    # explicit: exponent ``(i,0)`` here means coefficient-column i, not x^i.
    require_local = len(coefficient_labels) == number_of_unknowns
    if not require_local:
        raise RuntimeError("internal coefficient-label mismatch")
    return (vector_as_label_poly(particular_vector),
            [vector_as_label_poly(vector) for vector in homogeneous_vectors])


def coefficient_vector_to_polynomial(
    encoded: Poly, basis: list[Monomial]
) -> Poly:
    answer: Poly = {}
    for (index, dummy), coefficient in encoded.items():
        if dummy != 0:
            raise RuntimeError("bad encoded coefficient vector")
        answer[basis[index]] = coefficient
    return clean(answer)


def proportional(left: Poly, right: Poly) -> bool:
    if not left or not right:
        return left == right
    common = next(iter(set(left).intersection(right)), None)
    if common is None:
        return False
    ratio = left[common] / right[common]
    return left == scalar(right, ratio)


def coefficient_span_rank(entries: tuple[Poly, Poly, Poly, Poly]) -> int:
    support = sorted(set().union(*(set(entry) for entry in entries)) - {(0, 0)})
    if not support:
        return 0
    rows = [[as_fmpq(entry.get(monomial, Q(0))) for monomial in support]
            for entry in entries]
    return fmpq_mat(rows).rank()


def transform_signed_permutation(poly: Poly, action: str) -> Poly:
    answer: Poly = {}
    for (i, j), coefficient in poly.items():
        if action == "central":
            monomial, sign = (i, j), (-1) ** (i + j)
        elif action == "x_reflection":
            monomial, sign = (i, j), (-1) ** i
        elif action == "y_reflection":
            monomial, sign = (i, j), (-1) ** j
        elif action == "exchange":
            monomial, sign = (j, i), 1
        elif action == "negative_exchange":
            monomial, sign = (j, i), (-1) ** (i + j)
        else:
            raise RuntimeError(f"unknown action {action}")
        answer[monomial] = answer.get(monomial, Q(0)) + sign * coefficient
    return clean(answer)


def linear_span_coordinates(target: Poly, first: Poly, second: Poly) -> tuple[Q, Q] | None:
    support = sorted(set(target) | set(first) | set(second))
    pivot_pair: tuple[Monomial, Monomial] | None = None
    determinant = Q(0)
    for left, right in combinations(support, 2):
        determinant = (first.get(left, Q(0)) * second.get(right, Q(0))
                       - first.get(right, Q(0)) * second.get(left, Q(0)))
        if determinant:
            pivot_pair = (left, right)
            break
    if pivot_pair is None:
        return None
    left, right = pivot_pair
    a = (target.get(left, Q(0)) * second.get(right, Q(0))
         - target.get(right, Q(0)) * second.get(left, Q(0))) / determinant
    b = (first.get(left, Q(0)) * target.get(right, Q(0))
         - first.get(right, Q(0)) * target.get(left, Q(0))) / determinant
    if add(add(target, first, -a), second, -b):
        return None
    return a, b


def has_named_linear_involution(first: Poly, second: Poly) -> bool:
    for action in (
        "central", "x_reflection", "y_reflection", "exchange", "negative_exchange"
    ):
        row_one = linear_span_coordinates(
            transform_signed_permutation(first, action), first, second
        )
        row_two = linear_span_coordinates(
            transform_signed_permutation(second, action), first, second
        )
        if row_one is None or row_two is None:
            continue
        a, b = row_one
        c, d = row_two
        if (a * a + b * c == 1 and a * b + b * d == 0
                and c * a + d * c == 0 and c * b + d * d == 1):
            return True
    return False


def main() -> None:
    checks: list[str] = []

    zero: Poly = {}
    one: Poly = {(0, 0): Q(1)}
    x: Poly = {(1, 0): Q(1)}
    y: Poly = {(0, 1): Q(1)}
    x_squared = power(x, 2)
    xy = multiply(x, y)
    y_squared = power(y, 2)

    small_nonzero = (Q(-2), Q(-1), Q(-1, 2), Q(1, 2), Q(1), Q(2))
    small_with_zero = (Q(0),) + small_nonzero
    u_monomials = (x, y, x_squared, xy, y_squared)
    u_labels = ("x", "y", "x^2", "xy", "y^2")
    w_directions = {
        "y": y,
        "x": x,
        "x+y": add(x, y),
        "x-y": add(x, y, Q(-1)),
        "x+2y": add(x, y, Q(2)),
        "x-2y": add(x, y, Q(-2)),
        "2x+y": add(scalar(x, Q(2)), y),
        "2x-y": add(scalar(x, Q(2)), y, Q(-1)),
    }

    u_universe: list[tuple[str, Poly]] = []
    for support_size in (1, 2, 3):
        for support in combinations(range(len(u_monomials)), support_size):
            for coefficients in product(small_nonzero, repeat=support_size):
                polynomial: Poly = {}
                label_parts: list[str] = []
                for index, coefficient in zip(support, coefficients):
                    polynomial = add(polynomial, scalar(u_monomials[index], coefficient))
                    label_parts.append(f"{coefficient}*{u_labels[index]}")
                u_universe.append(("+".join(label_parts), polynomial))

    require("U universe count", len(u_universe) == 2550, checks)
    require("eight W directions", len(w_directions) == 8, checks)

    v_basis = monomial_basis(4)
    v_constant_column = v_basis.index((0, 0))
    first_row_histogram: Counter[tuple[str, str]] = Counter()
    first_row_survivors: list[tuple[Q, Poly, Poly]] = []
    symmetry_breaking_survivors = 0
    modular_certificates: Counter[str] = Counter()

    for w_label, w_poly in w_directions.items():
        for _, u_poly in u_universe:
            uw = multiply(u_poly, w_poly)
            delta = add(derivative(u_poly, 1), derivative(uw, 0), Q(-1))
            one_plus_uw = add(one, uw)
            operator_columns: list[Poly] = []
            for monomial in v_basis:
                column = add(
                    multiply(u_poly, derivative_monomial(monomial, 1)),
                    multiply(one_plus_uw, derivative_monomial(monomial, 0)),
                    Q(-1),
                )
                column = add(column, multiply(delta, {monomial: Q(1)}))
                operator_columns.append(column)

            certificate = full_rank_modular_certificate(
                operator_columns,
                derivative(w_poly, 0),
                forced_zero_columns=(v_constant_column,),
            )
            if certificate is not None:
                first_row_histogram[(w_label, certificate)] += 1
                modular_certificates[certificate] += 1
                continue

            solution = affine_solution_space(
                operator_columns,
                derivative(w_poly, 0),
                forced_zero_columns=(v_constant_column,),
            )
            if solution is None:
                first_row_histogram[(w_label, "inconsistent")] += 1
                continue
            particular_encoded, homogeneous_encoded = solution
            particular = coefficient_vector_to_polynomial(particular_encoded, v_basis)
            homogeneous = [coefficient_vector_to_polynomial(vector, v_basis)
                           for vector in homogeneous_encoded]
            require("homogeneous first-row systems have zero particular", not particular,
                    checks)
            nonzero_basis = [vector for vector in homogeneous if vector]
            if not nonzero_basis:
                first_row_histogram[(w_label, "zero-only")] += 1
                continue
            first_row_histogram[(w_label, f"dimension-{len(nonzero_basis)}")] += 1
            if delta:
                symmetry_breaking_survivors += 1

            # In the four surviving rows, U=alpha*(x+y^2/2) and the entire
            # degree-four V-kernel is Q*q, q=y+alpha*s^2/2.
            alpha = u_poly.get((1, 0), Q(0))
            shear_s = add(x, scalar(y_squared, Q(1, 2)))
            expected_u = scalar(shear_s, alpha)
            shear_q = add(y, scalar(power(shear_s, 2), alpha / 2))
            require("surviving U is a scaled shear coordinate", u_poly == expected_u,
                    checks)
            require("surviving first-row defect delta vanishes", not delta, checks)
            require("surviving V kernel has dimension one", len(nonzero_basis) == 1,
                    checks)
            require("surviving V kernel is Q*q", proportional(nonzero_basis[0], shear_q),
                    checks)
            first_row_survivors.append((alpha, shear_s, shear_q))

    expected_histogram = Counter({
        ("y", "zero-only"): 2546,
        ("y", "dimension-1"): 4,
        ("x", "inconsistent"): 2550,
        ("x+y", "inconsistent"): 2550,
        ("x-y", "inconsistent"): 2550,
        ("x+2y", "inconsistent"): 2550,
        ("x-2y", "inconsistent"): 2550,
        ("2x+y", "inconsistent"): 2550,
        ("2x-y", "inconsistent"): 2550,
    })
    require("first-row histogram", first_row_histogram == expected_histogram, checks)
    require("four first-row families", len(first_row_survivors) == 4, checks)
    require("no delta-nonzero first-row survivor", symmetry_breaking_survivors == 0,
            checks)

    z_basis = monomial_basis(8)
    z_constant_column = z_basis.index((0, 0))
    maps: list[tuple[Q, Q, Q, Poly, Poly]] = []
    second_row_dimensions: Counter[int] = Counter()

    for alpha, shear_s, shear_q in first_row_survivors:
        u_poly = scalar(shear_s, alpha)
        w_poly = y
        for t in small_nonzero:
            v_poly = scalar(shear_q, t)
            shear_p = add(shear_s, scalar(power(shear_q, 2), t / 2))
            a = add(one, multiply(u_poly, v_poly))
            b = add(v_poly, multiply(w_poly, a))
            operator_columns = [
                add(
                    multiply(a, derivative_monomial(monomial, 1)),
                    multiply(b, derivative_monomial(monomial, 0)),
                    Q(-1),
                )
                for monomial in z_basis
            ]
            solution = affine_solution_space(
                operator_columns,
                zero,
                forced_zero_columns=(z_constant_column,),
            )
            require("second-row system is consistent", solution is not None, checks)
            if solution is None:
                raise RuntimeError("unreachable")
            particular_encoded, homogeneous_encoded = solution
            particular = coefficient_vector_to_polynomial(particular_encoded, z_basis)
            homogeneous = [coefficient_vector_to_polynomial(vector, z_basis)
                           for vector in homogeneous_encoded]
            nonzero_basis = [vector for vector in homogeneous if vector]
            require("second-row particular is zero", not particular, checks)
            second_row_dimensions[len(nonzero_basis)] += 1
            require("second-row degree-eight kernel has dimension one",
                    len(nonzero_basis) == 1, checks)
            require("second-row kernel is Q*p", proportional(nonzero_basis[0], shear_p),
                    checks)

            for z_parameter in small_with_zero:
                shear_r = add(
                    shear_q,
                    scalar(power(shear_p, 2), z_parameter / 2),
                )
                maps.append((alpha, t, z_parameter, shear_p, shear_r))

    require("second-row dimension histogram", second_row_dimensions == Counter({1: 24}),
            checks)
    require("finite map count", len(maps) == 168, checks)

    span_histogram: Counter[int] = Counter()
    involution_intertwining = 0
    point_collision_hits = 0
    inverse_certificates = 0
    balanced_span_four_controls = 0
    point_bank = tuple(product(small_with_zero, repeat=2))

    for alpha, t, z_parameter, shear_p, shear_r in maps:
        shear_s = add(x, scalar(y_squared, Q(1, 2)))
        shear_q = add(y, scalar(power(shear_s, 2), alpha / 2))
        u_poly = scalar(shear_s, alpha)
        v_poly = scalar(shear_q, t)
        z_poly = scalar(shear_p, z_parameter)
        a = add(one, multiply(u_poly, v_poly))
        b = add(v_poly, multiply(y, a))
        c = add(u_poly, multiply(z_poly, a))
        d = add(add(one, multiply(u_poly, y)), multiply(z_poly, b))

        determinant = add(multiply(a, d), multiply(b, c), Q(-1))
        require("survivor determinant one", determinant == one, checks)
        require("survivor first curl", derivative(a, 1) == derivative(b, 0), checks)
        require("survivor second curl", derivative(c, 1) == derivative(d, 0), checks)
        span = coefficient_span_rank((a, b, c, d))
        span_histogram[span] += 1

        # The inverse is global, not merely a bounded collision test:
        # q=r-z*p^2/2; s=p-t*q^2/2; y=q-alpha*s^2/2; x=s-y^2/2.
        # Substitution in this triangular order cancels identically.
        inverse_q = add(shear_r, scalar(power(shear_p, 2), -z_parameter / 2))
        inverse_s = add(shear_p, scalar(power(inverse_q, 2), -t / 2))
        inverse_y = add(inverse_q, scalar(power(inverse_s, 2), -alpha / 2))
        inverse_x = add(inverse_s, scalar(power(inverse_y, 2), Q(-1, 2)))
        require("inverse recovers q", inverse_q == shear_q, checks)
        require("inverse recovers s", inverse_s == shear_s, checks)
        require("inverse recovers y", inverse_y == y, checks)
        require("inverse recovers x", inverse_x == x, checks)
        inverse_certificates += 1

        seen: dict[tuple[Q, Q], tuple[Q, Q]] = {}
        for point in point_bank:
            image = (
                evaluate(shear_p, point[0], point[1]),
                evaluate(shear_r, point[0], point[1]),
            )
            if image in seen and seen[image] != point:
                point_collision_hits += 1
            seen[image] = point

        if has_named_linear_involution(shear_p, shear_r):
            involution_intertwining += 1

        if z_parameter == 0:
            # Section 7 of the reflection proves absolute irreducibility of
            # the four displayed edges.  The exact rank computation here
            # verifies that all such balanced controls have full span.
            require("z=0 balanced control has coefficient span four", span == 4,
                    checks)
            balanced_span_four_controls += 1

    require("all maps have explicit inverse", inverse_certificates == 168, checks)
    require("no exact bounded-point collision", point_collision_hits == 0, checks)
    require("24 balanced full-span controls", balanced_span_four_controls == 24,
            checks)

    # A determinant-one hostile that dies at the first curl: U=W=x,V=0.
    hostile_u = x
    hostile_w = x
    hostile_v = zero
    hostile_a = add(one, multiply(hostile_u, hostile_v))
    hostile_b = add(hostile_v, multiply(hostile_w, hostile_a))
    hostile_first_curl = add(derivative(hostile_a, 1), derivative(hostile_b, 0), Q(-1))
    require("determinant-one hostile has first-curl defect -1",
            hostile_first_curl == {(0, 0): Q(-1)}, checks)

    # First live five-word orientation:
    # E_+(R) E_-(Z) E_+(V) E_-(U) E_+(W).  For the right three-word
    # prefix, epsilon=A_y-B_x and delta=U_y-(UW)_x.  The lower row after
    # E_-(Z) is closed iff
    #
    #     A Z_y-B Z_x+epsilon Z=-delta.
    #
    # This exact control takes U,V,W independently from the ten signed
    # nonconstant monomials +/-{x,y,x^2,xy,y^2} and solves the entire
    # degree-six Z space, including its constant.  The reflection proves
    # the stronger all-exponent, all-Z-degree monomial no-go.
    signed_monomials = tuple(
        scalar(monomial, sign)
        for monomial in u_monomials
        for sign in (Q(-1), Q(1))
    )
    require("ten signed monomial five-word parameters",
            len(signed_monomials) == 10, checks)
    depth_five_z_basis = monomial_basis(6)
    depth_five_histogram: Counter[str] = Counter()
    depth_five_systems = 0
    for u_five, v_five, w_five in product(signed_monomials, repeat=3):
        depth_five_systems += 1
        uv_five = multiply(u_five, v_five)
        a_five = add(one, uv_five)
        b_five = add(v_five, multiply(w_five, a_five))
        epsilon_five = add(
            derivative(a_five, 1), derivative(b_five, 0), Q(-1)
        )
        delta_five = add(
            derivative(u_five, 1),
            derivative(multiply(u_five, w_five), 0),
            Q(-1),
        )
        operator_columns = []
        for monomial in depth_five_z_basis:
            column = add(
                multiply(a_five, derivative_monomial(monomial, 1)),
                multiply(b_five, derivative_monomial(monomial, 0)),
                Q(-1),
            )
            column = add(column, multiply(epsilon_five, {monomial: Q(1)}))
            operator_columns.append(column)
        right_hand_side = scalar(delta_five, Q(-1))
        certificate = full_rank_modular_certificate(
            operator_columns, right_hand_side
        )
        if certificate == "inconsistent":
            depth_five_histogram["mod-1009-inconsistent"] += 1
            continue
        solution = affine_solution_space(operator_columns, right_hand_side)
        if solution is None:
            depth_five_histogram["exact-inconsistent"] += 1
        else:
            depth_five_histogram["survivor"] += 1

    require("one thousand signed-monomial five-word systems",
            depth_five_systems == 1000, checks)
    require("all signed-monomial systems classified",
            sum(depth_five_histogram.values()) == depth_five_systems, checks)
    require("signed-monomial five-word first stage is empty",
            depth_five_histogram.get("survivor", 0) == 0, checks)

    semantic_lines = [
        "status=FINITE-EXACT; no JC(2) conclusion outside the stated universe",
        "word=M=E_-(Z)E_+(V)E_-(U)E_+(W);det(M)=1",
        "small_nonzero={-2,-1,-1/2,1/2,1,2}",
        "W={y,x,x+y,x-y,x+2y,x-2y,2x+y,2x-y}; U=support<=3 in {x,y,x^2,xy,y^2}",
        "V=QQ[x,y]_(total_degree<=4),V(0)=0; curls solved over all QQ",
        "Z=QQ[x,y]_(total_degree<=8),Z(0)=0; curls solved over all QQ",
        f"raw_U={len(u_universe)};raw_UW={len(u_universe) * len(w_directions)}",
        f"full_rank_mod_1009_certificates={dict(sorted(modular_certificates.items()))}",
        "first_row_histogram=" + ",".join(
            f"{w}:{status}:{count}"
            for (w, status), count in sorted(first_row_histogram.items())
        ),
        ("delta_nonzero_nonunit_eligible_first_row_survivors="
         f"{symmetry_breaking_survivors}"),
        "survivor_alpha={-2,-1,1,2};s=x+y^2/2;q=y+alpha*s^2/2",
        "full_solution_spaces=V=QQ*q and Z=QQ*p after V(0)=Z(0)=0",
        "small_points:V=t*q,t in small_nonzero;p=s+t*q^2/2",
        "small_points:Z=z*p,z in small_nonzero union {0};r=q+z*p^2/2",
        f"maps={len(maps)};inverse_certificates={inverse_certificates}",
        f"coefficient_span_histogram={dict(sorted(span_histogram.items()))}",
        f"balanced_z0_span4_controls={balanced_span_four_controls}",
        f"named_linear_involution_intertwining_maps={involution_intertwining}",
        f"exact_point_bank={len(point_bank)};collision_hits={point_collision_hits}",
        "counterexample_survivors=0;mechanism=every survivor is four sequential shears",
        "hostile=U=x,W=x,V=0 has determinant one but first-curl defect -1",
        "depth5_word=E_+(R)E_-(Z)E_+(V)E_-(U)E_+(W)",
        ("depth5_base=U,V,W in +/-{x,y,x^2,xy,y^2};"
         "Z=QQ[x,y]_(total_degree<=6),constant allowed"),
        (f"depth5_first_correction_systems={depth_five_systems};"
         f"histogram={dict(sorted(depth_five_histogram.items()))}"),
    ]
    semantic_digest = sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()

    for line in semantic_lines:
        print(line)
    print(f"checks={len(checks)}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"source_sha256={source_digest}")


if __name__ == "__main__":
    main()
