#!/usr/bin/env python3
"""Independent exact verifier for the proposed THM-4376 finite-row lemma.

The candidate THM-4376 primary is neither imported nor read at runtime.  This
script starts from the already audited, stdlib-only THM-4366 referee's exact
row-ten P2-selected carrier, then independently:

* constructs every row-eleven P2/P3 projected-depth generator;
* computes the full exact ranks and left kernels by sparse rational RREF;
* appends the complete determinant-compatible degree-capped row eleven;
* compares full restricted annihilator rows with all admissible hierarchy rows;
* checks the new clock positions and constant active minors;
* eliminates the A equations and identifies the exact remaining C ideal; and
* checks that THM-4366's separate target/bracket residual kills that ideal.

Everything here is finite-row and relative to the declared source-normal
THM-4308/4315/4364/4366 universe.  In particular, JC(2) remains open.
"""

from __future__ import annotations

import contextlib
from fractions import Fraction as Q
import importlib.util
import io
from math import comb
from pathlib import Path
import sys
from typing import Any


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = (
    ROOT
    / "04-computation"
    / "jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_independent_referee_thm4366.py"
)
SPEC = importlib.util.spec_from_file_location("thm4366_clean_base", BASE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load audited THM-4366 clean-room base")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def dense_rank(matrix: list[list[Q]]) -> tuple[int, tuple[int, ...]]:
    """Exact rational rank and pivot columns."""

    if not matrix:
        return 0, ()
    rows = [[Q(value) for value in row] for row in matrix]
    nr, nc = len(rows), len(rows[0])
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(nc):
        selected = next(
            (row for row in range(pivot_row, nr) if rows[row][column]), None
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(nr):
            if row == pivot_row or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [
                rows[row][j] - scale * rows[pivot_row][j] for j in range(nc)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == nr:
            break
    return pivot_row, tuple(pivot_columns)


def independent_row_indices(matrix: list[list[Q]]) -> tuple[int, ...]:
    if not matrix:
        return ()
    _, pivots = dense_rank([list(column) for column in zip(*matrix)])
    return pivots


def determinant(matrix: list[list[Q]]) -> Q:
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "determinant square")
    rows = [[Q(value) for value in row] for row in matrix]
    result = Q(1)
    for column in range(size):
        selected = next((row for row in range(column, size) if rows[row][column]), None)
        if selected is None:
            return Q(0)
        if selected != column:
            rows[column], rows[selected] = rows[selected], rows[column]
            result = -result
        pivot = rows[column][column]
        result *= pivot
        for row in range(column + 1, size):
            if not rows[row][column]:
                continue
            scale = rows[row][column] / pivot
            for j in range(column + 1, size):
                rows[row][j] -= scale * rows[column][j]
    return result


def solve_square(matrix: list[list[Q]], rhs: list[Any]) -> list[Any]:
    """Solve a rational square matrix against Fraction or BASE.R values."""

    size = len(matrix)
    require(size == len(rhs) and all(len(row) == size for row in matrix), "solve square")
    rows = [[Q(value) for value in row] for row in matrix]
    values = list(rhs)
    for column in range(size):
        selected = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[selected] = rows[selected], rows[column]
        values[column], values[selected] = values[selected], values[column]
        pivot = rows[column][column]
        rows[column] = [value / pivot for value in rows[column]]
        values[column] = values[column] / pivot
        for row in range(size):
            if row == column or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [
                rows[row][j] - scale * rows[column][j] for j in range(size)
            ]
            values[row] = values[row] - scale * values[column]
    return values


def sparse_rref(
    sparse_rows: list[dict[int, Q]], number_columns: int
) -> tuple[list[dict[int, Q]], tuple[int, ...]]:
    rows = [
        {column: Q(value) for column, value in row.items() if value}
        for row in sparse_rows
    ]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(number_columns):
        selected = next(
            (row for row in range(pivot_row, len(rows)) if rows[row].get(column)),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = {
            index: value / scale for index, value in rows[pivot_row].items()
        }
        pivot = rows[pivot_row]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            scale = rows[row].get(column, Q(0))
            if not scale:
                continue
            updated = dict(rows[row])
            for index, value in pivot.items():
                new_value = updated.get(index, Q(0)) - scale * value
                if new_value:
                    updated[index] = new_value
                else:
                    updated.pop(index, None)
            rows[row] = updated
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, tuple(pivots)


def left_kernel(
    generator_columns: list[dict[int, Q]], number_coordinates: int
) -> tuple[int, list[dict[int, Q]]]:
    reduced, pivots = sparse_rref(generator_columns, number_coordinates)
    pivot_set = set(pivots)
    basis: list[dict[int, Q]] = []
    for free in range(number_coordinates):
        if free in pivot_set:
            continue
        vector = {free: Q(1)}
        for row, pivot in enumerate(pivots):
            value = reduced[row].get(free, Q(0))
            if value:
                vector[pivot] = -value
        basis.append(vector)
    return len(pivots), basis


def modular_rank(
    sparse_rows: list[dict[int, Q]], number_columns: int, prime: int
) -> int:
    rows: list[dict[int, int]] = []
    for row in sparse_rows:
        converted: dict[int, int] = {}
        for column, value in row.items():
            require(value.denominator % prime != 0, f"mod {prime} denominator")
            residue = value.numerator * pow(value.denominator, -1, prime) % prime
            if residue:
                converted[column] = residue
        rows.append(converted)
    pivot_row = 0
    for column in range(number_columns):
        selected = next(
            (row for row in range(pivot_row, len(rows)) if rows[row].get(column)),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = {
            index: value * inverse % prime
            for index, value in rows[pivot_row].items()
            if value * inverse % prime
        }
        pivot = rows[pivot_row]
        for row in range(pivot_row + 1, len(rows)):
            scale = rows[row].get(column, 0)
            if not scale:
                continue
            updated = dict(rows[row])
            for index, value in pivot.items():
                new_value = (updated.get(index, 0) - scale * value) % prime
                if new_value:
                    updated[index] = new_value
                else:
                    updated.pop(index, None)
            rows[row] = updated
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def depth_universe(
    maximum_row: int, depth: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]], list[dict[int, Q]]]:
    """Rebuild pi_m(P_d) in an order different from the candidate primary."""

    coordinates = [
        (row, degree)
        for row in range(maximum_row + 1)
        for degree in range(row + depth + 1)
    ]
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    monomials: list[tuple[int, int, int, int]] = []
    columns: list[dict[int, Q]] = []
    # Deliberately order e,c,b,a rather than the primary's b,a,e,c loop.
    for exponent_y in range(maximum_row // 2 + 1):
        for exponent_p in range(maximum_row + 1):
            for exponent_u in range(depth + 1):
                for exponent_x in range(depth - exponent_u + 1):
                    start_row = exponent_u + exponent_p + 2 * exponent_y
                    if start_row > maximum_row:
                        continue
                    packet_order = exponent_p + exponent_y
                    start_degree = exponent_x + 2 * exponent_u + exponent_y
                    column: dict[int, Q] = {}
                    for packet_index in range(packet_order + 1):
                        row = start_row + packet_index
                        degree = start_degree + 2 * packet_index
                        if row <= maximum_row:
                            coordinate = (row, degree)
                            require(coordinate in coordinate_index, "packet coordinate inside depth cap")
                            index = coordinate_index[coordinate]
                            column[index] = column.get(index, Q(0)) + Q(
                                comb(packet_order, packet_index)
                            )
                    monomials.append(
                        (exponent_x, exponent_u, exponent_p, exponent_y)
                    )
                    columns.append(column)
    require(len(set(monomials)) == len(monomials), "monomial enumeration unique")
    return coordinates, monomials, columns


def determinant_operator(maximum_row: int) -> list[list[Q]]:
    number_a = maximum_row + 2
    number_c = maximum_row + 3
    rows = [
        [Q(0) for _ in range(number_a + number_c)]
        for _ in range(maximum_row + 4)
    ]
    for degree in range(number_a):
        rows[degree][degree] += Q(3 * maximum_row, 4)
        rows[degree + 2][degree] += Q(3 * maximum_row, 8)
    for degree in range(number_c):
        rows[degree + 1][number_a + degree] += Q(maximum_row, 2)
    return rows


def canonical_particular(
    arows: dict[int, dict[int, Any]],
    crows: dict[int, dict[int, Any]],
    row: int,
) -> tuple[dict[int, Any], dict[int, Any]]:
    operator = determinant_operator(row)
    rank, pivots = dense_rank(operator)
    require(rank == row + 4 and len(pivots) == row + 4, f"row {row} operator surjective")
    bvalue = BASE.brow(arows, crows, row)
    rhs = [-BASE.xcoef(bvalue, degree) for degree in range(row + 4)]
    square = [[operator[i][column] for column in pivots] for i in range(row + 4)]
    solution = solve_square(square, rhs)
    values = [BASE.R() for _ in range(2 * row + 5)]
    for column, value in zip(pivots, solution):
        values[column] = value
    avalue = BASE.xp({degree: values[degree] for degree in range(row + 2)})
    cvalue = BASE.xp(
        {degree: values[row + 2 + degree] for degree in range(row + 3)}
    )
    lhs = BASE.xmany(
        [
            BASE.xscale(BASE.xmul(BASE.xder(arows[0]), cvalue), row),
            BASE.xscale(BASE.xmul(BASE.xder(crows[0]), avalue), -row),
            bvalue,
        ]
    )
    require(not lhs, f"row {row} canonical particular determinant identity")
    return avalue, cvalue


def poly_sub(left: dict[int, Any], right: dict[int, Any]) -> dict[int, Any]:
    return BASE.xadd(left, BASE.xneg(right))


def capture_thm4366_nonzero_branch() -> tuple[dict[str, Any], tuple[Any, Any]]:
    captured: dict[str, Any] = {}

    def tracer(frame, event, arg):
        if frame.f_code is BASE.branch_phi_nonzero.__code__ and event == "return":
            captured.update(frame.f_locals)
        return tracer

    BASE.CHECKS = 0
    old_trace = sys.gettrace()
    sys.settrace(tracer)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            result = BASE.branch_phi_nonzero()
    finally:
        sys.settrace(old_trace)
    require(bool(captured), "capture audited THM-4366 nonzero branch")
    require(isinstance(result, tuple) and len(result) == 2, "THM-4366 return ledger")
    return captured, result


ZERO_MONOMIAL = (0,) * BASE.NVAR


def rational_constant(value: Any) -> Q:
    value = BASE.rr(value)
    require(
        all(monomial == ZERO_MONOMIAL for monomial in value.terms),
        "extension direction coefficient is rational",
    )
    return value.terms.get(ZERO_MONOMIAL, Q(0))


AffineRow = tuple[list[Q], Any]


def affine_jet(
    constant_rows: dict[int, dict[int, Any]],
    directions: list[dict[int, dict[int, Any]]],
    coordinates: list[tuple[int, int]],
) -> list[AffineRow]:
    result: list[AffineRow] = []
    for row, degree in coordinates:
        constant = BASE.xcoef(constant_rows[row], degree)
        coefficients = [
            rational_constant(BASE.xcoef(direction.get(row, {}), degree))
            for direction in directions
        ]
        result.append((coefficients, constant))
    return result


def evaluate_left(left: dict[int, Q], jet: list[AffineRow]) -> AffineRow:
    number_variables = len(jet[0][0])
    coefficients = [Q(0) for _ in range(number_variables)]
    constant = BASE.R()
    for coordinate, weight in left.items():
        row_coefficients, row_constant = jet[coordinate]
        constant += row_constant * weight
        for variable, value in enumerate(row_coefficients):
            coefficients[variable] += weight * value
    return coefficients, constant


def affine_reduction(equations: list[AffineRow]) -> tuple[list[AffineRow], tuple[int, ...]]:
    if not equations:
        return [], ()
    rows = [([Q(value) for value in coefficients], constant) for coefficients, constant in equations]
    number_variables = len(rows[0][0])
    pivots: list[int] = []
    pivot_row = 0
    for column in range(number_variables):
        selected = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][0][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        coefficients, constant = rows[pivot_row]
        scale = coefficients[column]
        coefficients = [value / scale for value in coefficients]
        constant = constant / scale
        rows[pivot_row] = (coefficients, constant)
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][0][column]:
                continue
            scale = rows[row][0][column]
            updated_coefficients = [
                rows[row][0][j] - scale * coefficients[j]
                for j in range(number_variables)
            ]
            updated_constant = rows[row][1] - scale * constant
            rows[row] = (updated_coefficients, updated_constant)
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, tuple(pivots)


def affine_ranks(equations: list[AffineRow]) -> tuple[int, int, tuple[int, ...]]:
    reduced, pivots = affine_reduction(equations)
    coefficient_rank = len(pivots)
    inconsistent = any(
        not any(coefficients) and bool(constant) for coefficients, constant in reduced
    )
    return coefficient_rank, coefficient_rank + int(inconsistent), pivots


def independent_affine_rows(equations: list[AffineRow]) -> tuple[int, ...]:
    selected: list[int] = []
    current_rank = 0
    for index in range(len(equations)):
        trial = selected + [index]
        new_rank = affine_ranks([equations[j] for j in trial])[1]
        if new_rank > current_rank:
            selected.append(index)
            current_rank = new_rank
    return tuple(selected)


def hierarchy_labels(maximum_row: int, depth: int) -> list[tuple[int, int]]:
    return [
        (intercept, order)
        for intercept in range(2, 2 * maximum_row + 1)
        for order in range((intercept + 2) // 3)
        if maximum_row + order >= intercept
        and depth <= maximum_row + order - intercept
    ]


def hierarchy_left(
    coordinates: list[tuple[int, int]],
    maximum_row: int,
    intercept: int,
    order: int,
) -> dict[int, Q]:
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    start = (intercept + 1) // 2
    result: dict[int, Q] = {}
    for row in range(start, maximum_row + 1):
        coordinate = (row, 2 * row - intercept)
        if coordinate not in coordinate_index:
            continue
        value = Q((-1) ** (row - start) * comb(maximum_row + order - row, order))
        if value:
            result[coordinate_index[coordinate]] = value
    return result


def substitute_pivots(
    equations: list[AffineRow], pivot_columns: tuple[int, ...]
) -> tuple[dict[int, AffineRow], tuple[int, ...]]:
    coefficient_matrix = [coefficients for coefficients, _ in equations]
    rows = independent_row_indices(coefficient_matrix)
    require(len(rows) == len(pivot_columns), "pivot equation count")
    square = [
        [coefficient_matrix[row][column] for column in pivot_columns] for row in rows
    ]
    require(determinant(square) != 0, "pivot equation minor")
    free_columns = tuple(
        column
        for column in range(len(coefficient_matrix[0]))
        if column not in pivot_columns
    )
    constants = solve_square(square, [-equations[row][1] for row in rows])
    free_coefficients: dict[int, list[Q]] = {}
    for free in free_columns:
        free_coefficients[free] = solve_square(
            square, [-coefficient_matrix[row][free] for row in rows]
        )
    substitutions: dict[int, AffineRow] = {}
    for index, pivot in enumerate(pivot_columns):
        coefficients = [Q(0) for _ in free_columns]
        for free_index, free in enumerate(free_columns):
            coefficients[free_index] = free_coefficients[free][index]
        substitutions[pivot] = (coefficients, constants[index])
    return substitutions, free_columns


def apply_substitutions(
    equation: AffineRow,
    substitutions: dict[int, AffineRow],
    free_columns: tuple[int, ...],
) -> AffineRow:
    coefficients, constant = equation
    result_coefficients = [coefficients[column] for column in free_columns]
    result_constant = constant
    for pivot, (replacement_coefficients, replacement_constant) in substitutions.items():
        scale = coefficients[pivot]
        result_constant += scale * replacement_constant
        for index, value in enumerate(replacement_coefficients):
            result_coefficients[index] += scale * value
    return result_coefficients, result_constant


def polynomial_trim(poly: dict[int, Q]) -> dict[int, Q]:
    return {degree: Q(value) for degree, value in poly.items() if value}


def polynomial_divmod(
    dividend: dict[int, Q], divisor: dict[int, Q]
) -> tuple[dict[int, Q], dict[int, Q]]:
    remainder = polynomial_trim(dividend)
    divisor = polynomial_trim(divisor)
    quotient: dict[int, Q] = {}
    divisor_degree = max(divisor)
    divisor_lead = divisor[divisor_degree]
    while remainder and max(remainder) >= divisor_degree:
        degree = max(remainder) - divisor_degree
        scale = remainder[max(remainder)] / divisor_lead
        quotient[degree] = quotient.get(degree, Q(0)) + scale
        for index, value in divisor.items():
            target = index + degree
            remainder[target] = remainder.get(target, Q(0)) - scale * value
            if not remainder[target]:
                remainder.pop(target)
    return polynomial_trim(quotient), polynomial_trim(remainder)


def polynomial_gcd(left: dict[int, Q], right: dict[int, Q]) -> dict[int, Q]:
    left, right = polynomial_trim(left), polynomial_trim(right)
    while right:
        _, remainder = polynomial_divmod(left, right)
        left, right = right, remainder
    lead = left[max(left)]
    return {degree: value / lead for degree, value in left.items()}


def polynomial_derivative(poly: dict[int, Q]) -> dict[int, Q]:
    return {degree - 1: degree * value for degree, value in poly.items() if degree}


def main() -> None:
    captured, inherited_result = capture_thm4366_nonzero_branch()
    arows = dict(captured["aa"])
    crows = dict(captured["cc"])
    qpoly, rpoly = inherited_result
    require(captured["qpoly"] == qpoly and captured["rpoly"] == rpoly, "captured Q/R")
    require(set(arows) == set(range(11)) and set(crows) == set(range(11)), "rows zero through ten")

    expected_q = (
        BASE.P**6 * BASE.R.c(373891487235896675830078125)
        - BASE.P**4 * BASE.R.c(15097287707154073014589440000000)
        - BASE.P**2 * BASE.R.c(101452811911656563438652405841920000)
        - BASE.R.c(434321509795518334240224474125496745984)
    )
    expected_r = (
        BASE.P**8 * BASE.R.c(6846329377771290182382546697998046875)
        - BASE.P**6 * BASE.R.c(713835723041306505264998768716800000000000)
        - BASE.P**4 * BASE.R.c(2754991513504883058403611855707575418880000000)
        - BASE.P**2 * BASE.R.c(31916203206707002973657986739896008646412206080000)
        - BASE.R.c(156854967149983010817497418504735580308619473018945536)
    )
    require(qpoly == expected_q and rpoly == expected_r, "literal Q(Phi^2), R(Phi^2)")

    operator11 = determinant_operator(11)
    operator_rank, operator_pivots = dense_rank(operator11)
    require((len(operator11), len(operator11[0]), operator_rank) == (15, 27, 15), "row11 determinant operator")
    tangent_vectors: list[list[Q]] = []
    for power in range(12):
        vector = [Q(0) for _ in range(27)]
        vector[power + 1] = Q(1, 2)
        vector[13 + power] = Q(-3, 4)
        vector[13 + power + 2] = Q(-3, 8)
        tangent_vectors.append(vector)
        require(
            all(
                sum(operator11[row][column] * vector[column] for column in range(27)) == 0
                for row in range(15)
            ),
            f"theta11 tangent {power} in determinant kernel",
        )
    require(dense_rank(tangent_vectors)[0] == 12, "theta11 tangent independence")
    require(27 - operator_rank == 12, "theta11 tangent completeness")

    base_a11, base_c11 = canonical_particular(arows, crows, 11)
    b11 = BASE.brow(arows, crows, 11)
    constant_a = dict(arows)
    constant_c = dict(crows)
    constant_a[11] = base_a11
    constant_c[11] = base_c11
    directions_a: list[dict[int, dict[int, Any]]] = []
    directions_c: list[dict[int, dict[int, Any]]] = []

    for power in range(8):
        monomial = BASE.xp({power: 1})
        tangent_a10 = BASE.xmul(monomial, BASE.xder(arows[0]))
        tangent_c10 = BASE.xmul(monomial, BASE.xder(crows[0]))
        varied_a = dict(arows)
        varied_c = dict(crows)
        varied_a[10] = BASE.xadd(arows[10], tangent_a10)
        varied_c[10] = BASE.xadd(crows[10], tangent_c10)
        varied_a11, varied_c11 = canonical_particular(varied_a, varied_c, 11)
        tangent_a11 = poly_sub(varied_a11, base_a11)
        tangent_c11 = poly_sub(varied_c11, base_c11)
        directions_a.append({10: tangent_a10, 11: tangent_a11})
        directions_c.append({10: tangent_c10, 11: tangent_c11})
        delta_b = poly_sub(BASE.brow(varied_a, varied_c, 11), b11)
        linearized = BASE.xmany(
            [
                BASE.xscale(BASE.xmul(BASE.xder(arows[0]), tangent_c11), 11),
                BASE.xscale(BASE.xmul(BASE.xder(crows[0]), tangent_a11), -11),
                delta_b,
            ]
        )
        require(not linearized, f"theta10_{power} determinant-compatible continuation")

    for power in range(12):
        monomial = BASE.xp({power: 1})
        tangent_a11 = BASE.xmul(monomial, BASE.xder(arows[0]))
        tangent_c11 = BASE.xmul(monomial, BASE.xder(crows[0]))
        directions_a.append({11: tangent_a11})
        directions_c.append({11: tangent_c11})
        linearized = BASE.xadd(
            BASE.xmul(BASE.xder(arows[0]), tangent_c11),
            BASE.xneg(BASE.xmul(BASE.xder(crows[0]), tangent_a11)),
        )
        require(not linearized, f"theta11_{power} homogeneous determinant tangent")

    require(len(directions_a) == len(directions_c) == 20, "twenty affine coordinates")
    require(max(base_a11, default=-1) <= 12 and max(base_c11, default=-1) <= 13, "row11 base degree caps")
    require(
        all(max(direction.get(11, {}), default=-1) <= 12 for direction in directions_a),
        "row11 A direction degree caps",
    )
    require(
        all(max(direction.get(11, {}), default=-1) <= 13 for direction in directions_c),
        "row11 C direction degree caps",
    )

    acoordinates, amonomials, acolumns = depth_universe(11, 2)
    ccoordinates, cmonomials, ccolumns = depth_universe(11, 3)
    arank, anull = left_kernel(acolumns, len(acoordinates))
    crank, cnull = left_kernel(ccolumns, len(ccoordinates))
    require((len(acoordinates), len(amonomials), arank, len(anull)) == (102, 228, 77, 25), "complete P2 universe")
    require((len(ccoordinates), len(cmonomials), crank, len(cnull)) == (114, 361, 94, 20), "complete P3 universe")
    for prime in (101, 1009):
        require(modular_rank(acolumns, len(acoordinates), prime) == 77, f"P2 modular rank {prime}")
        require(modular_rank(ccolumns, len(ccoordinates), prime) == 94, f"P3 modular rank {prime}")
    for index, left in enumerate(anull):
        for column, generator in enumerate(acolumns):
            require(
                sum(left.get(coordinate, Q(0)) * value for coordinate, value in generator.items()) == 0,
                f"P2 left kernel {index}/{column}",
            )
    for index, left in enumerate(cnull):
        for column, generator in enumerate(ccolumns):
            require(
                sum(left.get(coordinate, Q(0)) * value for coordinate, value in generator.items()) == 0,
                f"P3 left kernel {index}/{column}",
            )

    ajet = affine_jet(constant_a, directions_a, acoordinates)
    cjet = affine_jet(constant_c, directions_c, ccoordinates)
    aequations = [evaluate_left(left, ajet) for left in anull]
    cequations = [evaluate_left(left, cjet) for left in cnull]
    aranks = affine_ranks(aequations)
    cranks = affine_ranks(cequations)
    jointranks = affine_ranks(aequations + cequations)
    require(aranks == (4, 4, (6, 7, 18, 19)), "restricted A ranks and pivots")
    require(cranks == (3, 4, (7, 18, 19)), "restricted C ranks and pivots")
    require(jointranks == (4, 5, (6, 7, 18, 19)), "restricted joint ranks and pivots")

    alabels = hierarchy_labels(11, 2)
    clabels = hierarchy_labels(11, 3)
    ahierarchy: list[AffineRow] = []
    chierarchy: list[AffineRow] = []
    for label in alabels:
        left = hierarchy_left(acoordinates, 11, *label)
        for column, generator in enumerate(acolumns):
            require(
                sum(left.get(coordinate, Q(0)) * value for coordinate, value in generator.items()) == 0,
                f"A hierarchy {label}/{column}",
            )
        ahierarchy.append(evaluate_left(left, ajet))
    for label in clabels:
        left = hierarchy_left(ccoordinates, 11, *label)
        for column, generator in enumerate(ccolumns):
            require(
                sum(left.get(coordinate, Q(0)) * value for coordinate, value in generator.items()) == 0,
                f"C hierarchy {label}/{column}",
            )
        chierarchy.append(evaluate_left(left, cjet))

    require(len(alabels) == 24 and len(clabels) == 19, "complete admissible hierarchy counts")
    require(affine_ranks(ahierarchy)[:2] == (4, 4), "A hierarchy restricted affine rank")
    require(affine_ranks(chierarchy)[:2] == (3, 4), "C hierarchy restricted affine rank")
    require(affine_ranks(aequations + ahierarchy)[:2] == (4, 4), "A restricted row-space equality")
    require(affine_ranks(cequations + chierarchy)[:2] == (3, 4), "C restricted row-space equality")
    abasis = independent_affine_rows(ahierarchy)
    cbasis = independent_affine_rows(chierarchy)
    require([alabels[index] for index in abasis] == [(10, 1), (11, 2), (12, 3), (13, 4)], "A hierarchy affine basis labels")
    require([clabels[index] for index in cbasis] == [(9, 1), (10, 2), (10, 3), (11, 3)], "C hierarchy affine basis labels")

    anew = sorted(set(alabels) - set(hierarchy_labels(10, 2)))
    cnew = sorted(set(clabels) - set(hierarchy_labels(10, 3)))
    require(anew == [(9, 0), (10, 1), (11, 2), (12, 3), (13, 4)], "new A clock positions")
    require(cnew == [(8, 0), (9, 1), (10, 2), (11, 3)], "new C clock positions")
    amap = dict(zip(alabels, ahierarchy))
    cmap = dict(zip(clabels, chierarchy))
    azero = [label for label in anew if not any(amap[label][0]) and not amap[label][1]]
    czero = [label for label in cnew if not any(cmap[label][0]) and not cmap[label][1]]
    require(azero == [(9, 0)] and czero == [(8, 0)], "unique clock silences")

    aactive = [amap[label] for label in anew if label not in azero]
    cactive = [cmap[label] for label in cnew if label not in czero]
    aminor_columns = (6, 7, 18, 19)
    cminor_columns = (7, 18, 19)
    aminor = [[row[0][column] for column in aminor_columns] for row in aactive]
    cminor = [[row[0][column] for column in cminor_columns] for row in cactive]
    adet = determinant(aminor)
    cdet = determinant(cminor)
    require(adet == Q(5, 4), "active A clock determinant 5/4")
    require(cdet == Q(27, 128), "active C clock determinant 27/128")
    require(affine_ranks(aequations + aactive)[:2] == (4, 4), "new A positions span full affine equations")
    require(dense_rank([row[0] for row in cequations + cactive])[0] == 3, "new C positions span full coefficient equations")

    # Check the two silence mechanisms directly, without relying on the
    # deliberately different coordinate lengths at rows ten and eleven.
    require(
        BASE.xcoef(constant_a[11], 13) == 0
        and all(BASE.xcoef(direction.get(11, {}), 13) == 0 for direction in directions_a),
        "A new diagonal coefficient x13 is degree-capped",
    )
    require(
        BASE.xcoef(constant_c[11], 14) == 0
        and all(BASE.xcoef(direction.get(11, {}), 14) == 0 for direction in directions_c),
        "C new diagonal coefficient x14 is degree-capped",
    )
    # Removing the final row term from each silent hierarchy leaves the
    # inherited row-ten trace, so zero of the full affine row checks it too.
    require(not any(amap[(9, 0)][0]) and not amap[(9, 0)][1], "A inherited trace silence")
    require(not any(cmap[(8, 0)][0]) and not cmap[(8, 0)][1], "C inherited trace silence")

    asubstitutions, remaining = substitute_pivots(aequations, aranks[2])
    require(remaining == tuple(list(range(6)) + list(range(8, 18))), "sixteen free coordinates after A")
    areduced = [apply_substitutions(equation, asubstitutions, remaining) for equation in aequations]
    require(all(not any(coefficients) and not constant for coefficients, constant in areduced), "A solve satisfies all A equations")
    creduced = [apply_substitutions(equation, asubstitutions, remaining) for equation in cequations]
    require(all(not any(coefficients) for coefficients, _ in creduced), "C coefficients collapse after A")
    nonzero_c = [(index, constant) for index, (_, constant) in enumerate(creduced) if constant]
    require([index for index, _ in nonzero_c] == [12, 17], "canonical C null-basis zero-based residual indices")
    expected_c_residuals = {
        12: -qpoly * BASE.R.mon(p=-5) / 15381515239492820800781250,
        17: -qpoly * BASE.R.mon(p=-5) / 5127171746497606933593750,
    }
    for index, value in nonzero_c:
        require(value == expected_c_residuals[index], f"exact C residual ratio {index}")
        require(BASE.zero_mod(value, qpoly), f"C residual {index} vanishes modulo Q")

    chierarchy_reduced = [
        apply_substitutions(equation, asubstitutions, remaining)
        for equation in chierarchy
    ]
    nonzero_chierarchy = [
        (label, constant)
        for label, (coefficients, constant) in zip(clabels, chierarchy_reduced)
        if any(coefficients) or constant
    ]
    require(len(nonzero_chierarchy) == 1 and nonzero_chierarchy[0][0] == (10, 3), "sole C hierarchy obstruction after A")
    require(
        not any(chierarchy_reduced[clabels.index((10, 3))][0])
        and nonzero_chierarchy[0][1]
        == qpoly * BASE.R.mon(p=-5) / 15381515239492820800781250,
        "exact inherited C hierarchy ratio",
    )

    q_y = {
        3: Q(373891487235896675830078125),
        2: Q(-15097287707154073014589440000000),
        1: Q(-101452811911656563438652405841920000),
        0: Q(-434321509795518334240224474125496745984),
    }
    r_y = {
        4: Q(6846329377771290182382546697998046875),
        3: Q(-713835723041306505264998768716800000000000),
        2: Q(-2754991513504883058403611855707575418880000000),
        1: Q(-31916203206707002973657986739896008646412206080000),
        0: Q(-156854967149983010817497418504735580308619473018945536),
    }
    require(q_y[0] != 0, "Q has no zero root")
    require(polynomial_gcd(q_y, polynomial_derivative(q_y)) == {0: Q(1)}, "Q squarefree")
    require(polynomial_gcd(q_y, r_y) == {0: Q(1)}, "Q and bracket R coprime")
    q_phi = {2 * degree: value for degree, value in q_y.items()}
    require(polynomial_gcd(q_phi, polynomial_derivative(q_phi)) == {0: Q(1)}, "Q(Phi^2) has six simple geometric points")

    bracket_residuals = captured["res11"]
    require([index for index, value in enumerate(bracket_residuals) if value] == [8, 9], "THM-4366 bracket residual support")
    require(
        bracket_residuals[9]
        == qpoly * BASE.R.mon(p=-5) * Q(8, 84598333817210514404296875),
        "bracket Q residual",
    )
    require(
        bracket_residuals[8]
        == -rpoly * BASE.R.mon(p=-6) / 71525718603423911567894897460937500,
        "bracket R residual",
    )

    # The exact Q multiples above are units times Q in Q[Phi,Phi^-1].
    # Therefore compatibility is iff Q=0, and the constant A minor gives a
    # rank-four system at every one of the six reduced geometric points.
    require(adet != 0 and len(remaining) == 16, "counterfactual fibre dimension sixteen")

    print("THM4376_CLEANROOM")
    print("INHERITED_THM4366_STDLIB_CHECKS", BASE.CHECKS)
    print("DETERMINANT_OPERATOR", 15, 27, operator_rank, "KERNEL", 12)
    print("P2_UNIVERSE", len(acoordinates), len(acolumns), arank, len(anull), "MOD", 77, 77)
    print("P3_UNIVERSE", len(ccoordinates), len(ccolumns), crank, len(cnull), "MOD", 94, 94)
    print("RESTRICTED_RANKS", "A", aranks[:2], "C", cranks[:2], "JOINT", jointranks[:2])
    print("HIERARCHY_COUNTS", len(alabels), len(clabels), "ROW_SPACES_EQUAL", True)
    print("NEW_CLOCK_A", anew, "SILENT", azero, "DET", adet)
    print("NEW_CLOCK_C", cnew, "SILENT", czero, "DET", cdet)
    print("C_AFTER_A", len(remaining), "COEFF_RANK", 0, "RESIDUAL_INDICES_ZERO_BASED", [index for index, _ in nonzero_c])
    print("C_AFTER_A_IDEAL", "Q(Phi^2)", "FIBRE_DIMENSION", 16)
    print("BRACKET_SIDECAR", "gcd(Q,R)=1", "ACTUAL_ROW11_BRACKET_FIBRE_EMPTY")
    print("CAUSAL_SCOPE", "COUNTERFACTUAL_DETERMINANT_ROW_ONLY_PROJECTED_DEPTH")
    print("JC(2)_STATUS", "OPEN")
    print(f"CHECKS={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()
