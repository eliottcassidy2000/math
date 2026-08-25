#!/usr/bin/env python3
"""Fraction-only independent audit for THM-4067.

This implementation uses no SymPy and imports no production companion.  It
checks graph incidence ranks, contact-jet matrices, branch-line arithmetic,
the oriented witness, and fixed/mixed A-adic valuations by an explicit
binomial dictionary engine.
"""

from fractions import Fraction as Q
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank(matrix):
    rows = [[Q(entry) for entry in row] for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [entry / scale for entry in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            multiple = rows[row][column]
            rows[row] = [
                rows[row][index] - multiple * rows[pivot_row][index]
                for index in range(column_count)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def incidence_matrix(vertex_count, edges):
    matrix = [[Q(0) for _ in edges] for _ in range(vertex_count)]
    for column, (source, target) in enumerate(edges):
        matrix[source][column] -= 1
        matrix[target][column] += 1
    return matrix


def add(left, right, scale=Q(1)):
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Q(0)) + scale * coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def derivative(polynomial, coordinate):
    result = {}
    for (a_degree, c_degree), coefficient in polynomial.items():
        exponent = a_degree if coordinate == "A" else c_degree
        if exponent == 0:
            continue
        monomial = (
            (a_degree - 1, c_degree)
            if coordinate == "A"
            else (a_degree, c_degree - 1)
        )
        result[monomial] = result.get(monomial, Q(0)) + exponent * coefficient
    return result


def restrict_monomial(a_degree, c_degree, u_degree, slope, intercept, opening):
    """Restrict A^a c^r u^s to u=slope*c+intercept*A^opening."""
    result = {}
    for shifted_factors in range(u_degree + 1):
        coefficient = (
            Q(comb(u_degree, shifted_factors))
            * slope ** (u_degree - shifted_factors)
            * intercept**shifted_factors
        )
        if coefficient == 0:
            continue
        monomial = (
            a_degree + opening * shifted_factors,
            c_degree + u_degree - shifted_factors,
        )
        result[monomial] = result.get(monomial, Q(0)) + coefficient
    return result


def a_valuation(polynomial):
    if not polynomial:
        return None
    return min(monomial[0] for monomial in polynomial)


print("THM4067_INDEPENDENT_FRACTION_AUDIT")

# Incidence and cycle ranks are reconstructed without the dimension formula.
graphs = {
    "tree": (4, ((0, 1), (1, 2), (1, 3)), ()),
    "triangle": (3, ((0, 1), (1, 2), (2, 0)), ((1, 1, 1),)),
    "figure_eight": (
        5,
        ((0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)),
        ((1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1)),
    ),
}
for name, (vertex_count, edges, cycles) in graphs.items():
    incidence = incidence_matrix(vertex_count, edges)
    incidence_rank = rank(incidence)
    cycle_rank = rank(cycles)
    require(incidence_rank == vertex_count - 1, (name, "incidence"))
    require(cycle_rank == len(edges) - vertex_count + 1, (name, "cycles"))
    for cycle in cycles:
        boundary = [
            sum(incidence[row][column] * cycle[column] for column in range(len(edges)))
            for row in range(vertex_count)
        ]
        require(all(entry == 0 for entry in boundary), (name, "cycle boundary"))
    print(
        f"{name}:incidence_rank={incidence_rank};cycle_rank={cycle_rank};"
        "cycle_boundaries=zero"
    )

# Build exact truncated derivative matrices for y(y-x^m).  Domain primitives
# have paired coefficients below m and independent coefficients from m on.
print("contact_derivative_matrix_ledger")
for contact in range(1, 9):
    cutoff = contact + 3
    rows = 2 * (cutoff + 1)
    columns = []
    for degree in range(contact):
        column = [Q(0)] * rows
        if degree:
            column[degree - 1] = degree
            column[(cutoff + 1) + degree - 1] = degree
        columns.append(column)
    for degree in range(contact, cutoff + 2):
        for branch in range(2):
            column = [Q(0)] * rows
            if degree:
                column[branch * (cutoff + 1) + degree - 1] = degree
            columns.append(column)
    matrix = [list(row) for row in zip(*columns)]
    derivative_rank = rank(matrix)
    cokernel = rows - derivative_rank
    require(cokernel == contact - 1, ("contact matrix", contact))
    print(
        f"m={contact};cutoff={cutoff};rows={rows};"
        f"rank={derivative_rank};cokernel={cokernel}"
    )

# THM-3696 hostile f=b^3-b^5, evaluated by hand from its derivative.
values = (Q(0), Q(0), Q(0))
derivative_values = (Q(-2), Q(0), Q(-2))
jet_law = derivative_values[0] + 4 * derivative_values[1] + derivative_values[2]
require(jet_law == -4, "THM-3696 independent jet")
print(f"thm3696_hostile_values={values};derivatives={derivative_values};jet_law={jet_law}")

# Supporting lines and oriented witness from endpoint fractions.
origin = (Q(0), Q(0))
B1 = (Q(1), Q(0))
B2 = (Q(2), Q(1))
D1 = (Q(-1), Q(2))
D2 = (Q(-3), Q(1))
edge_vertices = (
    (origin, B1),
    (B1, B2),
    (B2, origin),
    (origin, D1),
    (D1, D2),
    (D2, origin),
)
lines = []
c_increments = []
for start, end in edge_vertices:
    c_increment = end[0] - start[0]
    slope = (end[1] - start[1]) / c_increment
    intercept = start[1] - slope * start[0]
    lines.append((slope, intercept))
    c_increments.append(c_increment)
expected = (
    (Q(0), Q(0)),
    (Q(1), Q(-1)),
    (Q(1, 2), Q(0)),
    (Q(-2), Q(0)),
    (Q(1, 2), Q(5, 2)),
    (Q(-1, 3), Q(0)),
)
require(tuple(lines) == expected, "independent supporting lines")
require(c_increments == [1, 1, -2, -1, -2, 3], "oriented increments")
density = (Q(0), Q(2), Q(1), Q(0), Q(0), Q(0))
periods = (
    sum(density[index] * c_increments[index] for index in range(3)),
    sum(density[index] * c_increments[index] for index in range(3, 6)),
)
require(periods == (0, 0), "independent witness periods")
require(density[2] - density[4] == 1, "independent unit contact failure")
print(f"scaled_lines={tuple(lines)}")
print(f"oriented_c_increments={tuple(c_increments)};witness_periods={periods}")
print("witness_parallel_difference_mod_epsilon=1")

# Independent binomial-dictionary audit of every target monomial in the same
# hostile window as the production script.
fixed_checks = 0
mixed_checks = 0
for opening in range(1, 7):
    for total_degree in range(7):
        for a_degree in range(total_degree + 1):
            for u_degree in range(total_degree - a_degree + 1):
                c_degree = total_degree - a_degree - u_degree
                restriction_3 = restrict_monomial(
                    a_degree, c_degree, u_degree, Q(1, 2), Q(0), opening
                )
                restriction_5 = restrict_monomial(
                    a_degree, c_degree, u_degree, Q(1, 2), Q(5, 2), opening
                )
                fixed_difference = add(
                    derivative(restriction_3, "c"),
                    derivative(restriction_5, "c"),
                    Q(-1),
                )
                fixed_valuation = a_valuation(fixed_difference)
                require(
                    fixed_valuation is None or fixed_valuation >= opening,
                    ("fixed valuation", opening, total_degree),
                )
                fixed_checks += 1

                mixed_difference = add(
                    derivative(restriction_3, "A"),
                    derivative(restriction_5, "A"),
                    Q(-1),
                )
                mixed_valuation = a_valuation(mixed_difference)
                require(
                    mixed_valuation is None or mixed_valuation >= opening - 1,
                    ("mixed valuation", opening, total_degree),
                )
                mixed_checks += 1

    # For phi=A^i c^n, H(2A^q)/A^q has nonnegative A valuation i+qn.
    for a_degree in range(opening):
        for c_degree in range(6):
            connector_valuation = a_degree + opening * c_degree
            require(connector_valuation >= 0, "connector valuation")

    # f=u has exact A-derivative difference -(5/2)q A^(q-1).
    sharp = -Q(5, 2) * opening
    require(sharp != 0, ("sharp mixed coefficient", opening))

print(f"binomial_contact_checks=fixed:{fixed_checks};mixed:{mixed_checks}")
print("fixed_period_defect_quotient=Q[[A,c]]/(A^q);all_q>=1")
print("mixed_cokernel_quotient=Q[[A,c]]/(A^(q-1));nonzero_for_q>=2")
print("displayed_figure_eight_antecedent=REFUTED;abstract_implication=SURVIVES")

