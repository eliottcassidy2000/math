#!/usr/bin/env python3
"""Standalone split-prime hostile audit for the reserved THM-4060 candidate.

This file imports no repository code.  It reconstructs Q_alpha, all three
source jets, the affine mixed-pair tangent, an independently generated exact
two-form bank, and the leading affine triangle directly over F_137.
"""

from hashlib import sha256
from math import comb
from time import perf_counter

from flint import nmod_mat


PRIME = 137
ROOTS = (44, 82, 92, 134)
POINTS = (-1, 0, 1)
BRANCH_S = (2, 0, -2)
CUTOFFS = tuple(range(4, 13))
CARRIER_K = tuple(range(13))
MAX_RUNTIME_SECONDS = 120
EXPECTED_TABLE_SHA256 = "6833a2c709ed77ab45f3e02cd3bdda198e649d5c9f546740b9ccc602962ff0fe"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inverse(value):
    value %= PRIME
    require(value != 0, ("division by zero", value))
    return pow(value, -1, PRIME)


def rational(numerator, denominator=1):
    return numerator % PRIME * inverse(denominator) % PRIME


def add(left, right):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = (answer.get(monomial, 0) + coefficient) % PRIME
        if not answer[monomial]:
            del answer[monomial]
    return answer


def scale(value, scalar):
    scalar %= PRIME
    return {
        monomial: coefficient * scalar % PRIME
        for monomial, coefficient in value.items()
        if coefficient * scalar % PRIME
    }


def multiply(left, right, cutoff):
    answer = {}
    for (i, j), a_value in left.items():
        for (u, v), b_value in right.items():
            if i + j + u + v <= cutoff:
                key = (i + u, j + v)
                answer[key] = (
                    answer.get(key, 0) + a_value * b_value
                ) % PRIME
    return {key: value for key, value in answer.items() if value}


def power(value, exponent, cutoff):
    answer = {(0, 0): 1}
    factor = value
    while exponent:
        if exponent & 1:
            answer = multiply(answer, factor, cutoff)
        factor = multiply(factor, factor, cutoff)
        exponent //= 2
    return answer


def differentiate(value, axis):
    answer = {}
    for (i, j), coefficient in value.items():
        exponent = i if axis == 0 else j
        if exponent:
            key = (i - (axis == 0), j - (axis == 1))
            answer[key] = coefficient * exponent % PRIME
    return answer


def series_inverse(value, cutoff):
    constant = value.get((0, 0), 0)
    require(constant != 0, "nonunit power series")
    constant_inverse = inverse(constant)
    remainder = add(scale(value, constant_inverse), {(0, 0): -1 % PRIME})
    answer = {}
    term = {(0, 0): 1}
    for degree in range(cutoff + 1):
        answer = add(answer, scale(term, -1 if degree % 2 else 1))
        term = multiply(term, remainder, cutoff)
    answer = scale(answer, constant_inverse)
    require(
        multiply(value, answer, cutoff) == {(0, 0): 1},
        "series inverse residual",
    )
    return answer


def shift_second(value, amount, cutoff):
    return {
        (i, j + amount): coefficient
        for (i, j), coefficient in value.items()
        if i + j + amount <= cutoff
    }


def univariate_add(left, right):
    answer = dict(left)
    for degree, coefficient in right.items():
        answer[degree] = (answer.get(degree, 0) + coefficient) % PRIME
        if not answer[degree]:
            del answer[degree]
    return answer


def univariate_scale(value, scalar):
    scalar %= PRIME
    return {
        degree: coefficient * scalar % PRIME
        for degree, coefficient in value.items()
        if coefficient * scalar % PRIME
    }


def univariate_shift(value, amount):
    return {degree + amount: coefficient for degree, coefficient in value.items()}


def polynomial_value(polynomial, value):
    return sum(coefficient * pow(value, degree, PRIME)
               for degree, coefficient in polynomial.items()) % PRIME


def polynomial_derivative_value(polynomial, value):
    return sum(
        degree * coefficient * pow(value, degree - 1, PRIME)
        for degree, coefficient in polynomial.items()
        if degree
    ) % PRIME


def expand_around(polynomial, point, cutoff):
    answer = {}
    for degree, coefficient in polynomial.items():
        for local_degree in range(min(degree, cutoff) + 1):
            key = (local_degree, 0)
            answer[key] = (
                answer.get(key, 0)
                + coefficient * comb(degree, local_degree)
                * pow(point, degree - local_degree, PRIME)
            ) % PRIME
    return {key: value for key, value in answer.items() if value}


Q_BASE = {
    0: rational(-3, 4),
    1: 1,
    2: rational(-27, 4),
    3: -2 % PRIME,
    4: rational(9, 2),
    5: 1,
}
P_BASE = {2: 1, 4: -2 % PRIME, 6: 1}
Q6 = univariate_add(Q_BASE, univariate_scale(P_BASE, rational(-259, 36)))
R1 = univariate_add(P_BASE, univariate_scale(univariate_shift(P_BASE, 2), -1))
R2 = univariate_add(
    univariate_scale(P_BASE, 4),
    univariate_scale(univariate_shift(P_BASE, 1), -9),
)


def q_polynomial(root):
    theta = (
        rational(520, 9) * root * root
        - rational(1688, 81) * root
        - rational(5717, 729)
    ) % PRIME
    return univariate_add(
        Q6,
        univariate_add(univariate_scale(R1, theta), univariate_scale(R2, root)),
    )


def rho_reduction(root):
    return (
        rational(-2073506706944, 1678822119)
        + rational(372679949312, 62178597) * root
        - rational(184159683584, 6908733) * root * root
        - rational(73442787328, 2302911) * root * root * root
    ) % PRIME


def quartic_value(root):
    return (
        72783360 * root ** 4
        - 77822208 * root ** 3
        - 28419741 * root ** 2
        + 7849770 * root
        - 1276420
    ) % PRIME


def compose_univariate(polynomial, value, cutoff):
    answer = {}
    for degree, coefficient in polynomial.items():
        answer = add(answer, scale(power(value, degree, cutoff), coefficient))
    return answer


def graph_branches(root, cutoff=6):
    qalpha = q_polynomial(root)
    variable_c = {(0, 1): 1}
    target_a = {(0, 0): rational(-3, 4), (1, 0): 1}
    answer = []
    for base_s in BRANCH_S:
        correction = {}
        for iteration in range(5):
            s_value = add({(0, 0): base_s % PRIME}, correction)
            equation = add(
                add(
                    scale(s_value, 3),
                    multiply(target_a, power(s_value, 3, cutoff), cutoff),
                ),
                scale(variable_c, -1),
            )
            derivative = add(
                {(0, 0): 3},
                scale(
                    multiply(target_a, power(s_value, 2, cutoff), cutoff),
                    3,
                ),
            )
            correction = add(
                correction,
                scale(
                    multiply(equation, series_inverse(derivative, cutoff), cutoff),
                    -1,
                ),
            )
            require(iteration < 5, "Newton iteration sentinel")
        s_value = add({(0, 0): base_s % PRIME}, correction)
        equation = add(
            add(
                scale(s_value, 3),
                multiply(target_a, power(s_value, 3, cutoff), cutoff),
            ),
            scale(variable_c, -1),
        )
        require(not equation, ("implicit branch residual", root, base_s))
        denominator = add(
            {(0, 0): 1},
            multiply(target_a, power(s_value, 2, cutoff), cutoff),
        )
        source_x = multiply(s_value, series_inverse(denominator, cutoff), cutoff)
        source_q = multiply(target_a, power(denominator, 2, cutoff), cutoff)
        graph = add(
            source_q,
            scale(compose_univariate(qalpha, source_x, cutoff), -1),
        )
        answer.append(graph)
    return answer


def substitute_c(polynomial, c_series, cutoff):
    answer = {}
    for (a_degree, c_degree), coefficient in polynomial.items():
        term = scale(power(c_series, c_degree, cutoff), coefficient)
        for (degree, dummy_degree), value in term.items():
            require(dummy_degree == 0, "c substitution acquired a c-variable")
            if degree + a_degree <= cutoff:
                key = (degree + a_degree, 0)
                answer[key] = (answer.get(key, 0) + value) % PRIME
    return {key: value for key, value in answer.items() if value}


def implicit_c_root(difference, cutoff):
    linear = difference.get((0, 1), 0)
    require(linear != 0, "pairwise graph slope collision")
    answer = {}
    for degree in range(1, cutoff + 1):
        coefficient = substitute_c(difference, answer, cutoff).get((degree, 0), 0)
        if coefficient:
            answer[(degree, 0)] = -coefficient * inverse(linear) % PRIME
    require(not substitute_c(difference, answer, cutoff), "pairwise root residual")
    return answer


def triangle_and_carrier_audit(root):
    cutoff = 6
    branches = graph_branches(root, cutoff)
    slopes = (rational(3, 4), rational(-1, 3), rational(-3, 4))
    for branch, slope in zip(branches, slopes):
        require(branch.get((1, 0), 0) == 1, ("graph A slope", root))
        require(branch.get((0, 1), 0) == slope, ("graph c slope", root))

    vertices = {}
    for i, j in ((0, 1), (1, 2), (2, 0)):
        c_value = implicit_c_root(add(branches[i], scale(branches[j], -1)), cutoff)
        u_left = substitute_c(branches[i], c_value, cutoff)
        u_right = substitute_c(branches[j], c_value, cutoff)
        require(u_left == u_right, ("pairwise graph value", root, i, j))
        vertices[(i, j)] = (c_value, u_left)

    vertex_order = (vertices[(2, 0)], vertices[(0, 1)], vertices[(1, 2)])
    for coordinate in (0, 1):
        for degree in range(1, 5):
            values = [
                vertex[coordinate].get((degree, 0), 0)
                for vertex in vertex_order
            ]
            require(
                len(set(values)) == 1,
                ("common triangle jet", root, coordinate, degree, values),
            )

    c_edges = []
    u_edges = []
    cyclic = vertex_order[1:] + vertex_order[:1]
    for start, end in zip(vertex_order, cyclic):
        c_edges.append(
            (end[0].get((5, 0), 0) - start[0].get((5, 0), 0)) % PRIME
        )
        u_edges.append(
            (end[1].get((5, 0), 0) - start[1].get((5, 0), 0)) % PRIME
        )

    rho = rho_reduction(root)
    delta = c_edges[2]
    require(rho != 0 and delta != 0, ("zero rho or delta", root, rho, delta))
    require(delta == rational(26, 15) * rho % PRIME, ("delta/rho", root))
    require(
        c_edges == [
            delta * rational(5, 13) % PRIME,
            delta * rational(-18, 13) % PRIME,
            delta,
        ],
        ("leading c edges", root, c_edges),
    )
    require(
        u_edges == [slope * edge % PRIME for slope, edge in zip(slopes, c_edges)],
        ("leading u edges", root, u_edges),
    )

    current_u = 0
    area = 0
    for c_edge, slope in zip(c_edges, slopes):
        area = (
            area + current_u * c_edge + rational(1, 2) * slope * c_edge * c_edge
        ) % PRIME
        current_u = (current_u + slope * c_edge) % PRIME
    require(current_u == 0 and sum(c_edges) % PRIME == 0, ("triangle closure", root))
    require(area == rational(-15, 52) * delta * delta % PRIME,
            ("triangle area", root, area))
    require(area * inverse(delta) % PRIME == rational(-1, 2) * rho % PRIME,
            ("normalized u period", root, area))

    carrier_responses = []
    for k_value in CARRIER_K:
        # Pi(A^k u) starts with area*A^(k+10).  The owner operation
        # -4*d/dA therefore has normalized leading response below.
        response = (
            -4 * (k_value + 10) * area * inverse(delta)
        ) % PRIME
        expected = 2 * (k_value + 10) * rho % PRIME
        require(response == expected,
                ("carrier response", root, k_value, response, expected))
        carrier_responses.append(response)
    return rho, delta, area, tuple(carrier_responses)


def source_branches(root, cutoff):
    jet_cutoff = cutoff + 1
    qalpha = q_polynomial(root)
    answer = {}
    for point in POINTS:
        qjet = add(expand_around(qalpha, point, jet_cutoff), {(0, 1): 1})
        xjet = {(0, 0): point % PRIME, (1, 0): 1}
        denominator = add(
            {(0, 0): 1},
            multiply(power(xjet, 2, jet_cutoff), qjet, jet_cutoff),
        )
        denominator_inverse = series_inverse(denominator, jet_cutoff)
        ajet = multiply(
            qjet,
            power(denominator_inverse, 2, jet_cutoff),
            jet_cutoff,
        )
        Ajet = add(ajet, {(0, 0): rational(3, 4)})
        cjet = multiply(
            multiply(xjet, denominator, jet_cutoff),
            add(denominator, {(0, 0): 2}),
            jet_cutoff,
        )
        require(Ajet.get((0, 0), 0) == 0, ("A does not vanish", root, point))
        require(cjet.get((0, 0), 0) == 0, ("c does not vanish", root, point))

        a_x = differentiate(ajet, 0)
        a_u = differentiate(ajet, 1)
        c_x = differentiate(cjet, 0)
        c_u = differentiate(cjet, 1)
        jacobian = add(
            multiply(a_x, c_u, cutoff),
            scale(multiply(a_u, c_x, cutoff), -1),
        )
        require(jacobian == {(0, 0): -3 % PRIME},
                ("Jac(A,c) identity", root, point, cutoff, jacobian))

        A_powers = [power(Ajet, degree, jet_cutoff)
                    for degree in range(cutoff + 2)]
        c_powers = [power(cjet, degree, jet_cutoff)
                    for degree in range(cutoff + 2)]
        ac = {}
        axac = {}
        cxac = {}
        for a_degree in range(cutoff + 2):
            for c_degree in range(cutoff + 2 - a_degree):
                product = multiply(
                    A_powers[a_degree], c_powers[c_degree], jet_cutoff,
                )
                ac[(a_degree, c_degree)] = product
                axac[(a_degree, c_degree)] = multiply(a_x, product, cutoff)
                cxac[(a_degree, c_degree)] = multiply(c_x, product, cutoff)
        answer[point] = {"ac": ac, "axac": axac, "cxac": cxac}
    return answer


def retained_rows(cutoff):
    return tuple(
        (stable_degree, point, source_degree)
        for stable_degree in range(cutoff + 1)
        for source_degree in range(cutoff + 1 - stable_degree)
        for point in POINTS
    )


def make_column(branches, rows, cutoff, descriptors):
    values = {}
    for point in POINTS:
        density = {}
        for scalar, family, a_degree, c_degree, stable_degree in descriptors:
            require(a_degree >= 0 and c_degree >= 0 and stable_degree >= 0,
                    ("negative descriptor", descriptors))
            term = shift_second(
                branches[point][family][(a_degree, c_degree)],
                stable_degree,
                cutoff,
            )
            density = add(density, scale(term, scalar))
        values[point] = density
    return [
        values[point].get((source_degree, stable_degree), 0)
        for stable_degree, point, source_degree in rows
    ]


def build_banks(root, cutoff):
    rows = retained_rows(cutoff)
    branches = source_branches(root, cutoff)
    fixed_columns = []
    mixed_columns = []
    exact_columns = []
    column_identity_count = 0

    for a_degree in range(cutoff + 2):
        for c_degree in range(cutoff + 2 - a_degree):
            for u_degree in range(cutoff + 2 - a_degree - c_degree):
                f_terms = []
                g_terms = []
                U_terms = []
                V_terms = []
                W_terms = []

                if a_degree:
                    f_terms.append((12 * a_degree, "ac", a_degree - 1,
                                    c_degree, u_degree))
                    V_terms.append((-3 * a_degree, "ac", a_degree - 1,
                                    c_degree, u_degree))
                    W_terms.append((a_degree, "axac", a_degree - 1,
                                    c_degree, u_degree))
                if c_degree:
                    g_terms.append((-3 * c_degree, "ac", a_degree,
                                    c_degree - 1, u_degree))
                    U_terms.append((3 * c_degree, "ac", a_degree,
                                    c_degree - 1, u_degree))
                    W_terms.append((c_degree, "cxac", a_degree,
                                    c_degree - 1, u_degree))
                if u_degree:
                    f_terms.append((4 * u_degree, "cxac", a_degree,
                                    c_degree, u_degree - 1))
                    g_terms.append((u_degree, "axac", a_degree,
                                    c_degree, u_degree - 1))
                    U_terms.append((-u_degree, "axac", a_degree,
                                    c_degree, u_degree - 1))
                    V_terms.append((-u_degree, "cxac", a_degree,
                                    c_degree, u_degree - 1))

                f_column = make_column(branches, rows, cutoff, f_terms)
                g_column = make_column(branches, rows, cutoff, g_terms)
                U_column = make_column(branches, rows, cutoff, U_terms)
                V_column = make_column(branches, rows, cutoff, V_terms)

                require(
                    g_column == [(-value) % PRIME for value in U_column],
                    ("g versus U exact-form identity", root, cutoff,
                     a_degree, c_degree, u_degree),
                )
                require(
                    f_column == [(-4 * value) % PRIME for value in V_column],
                    ("f versus V exact-form identity", root, cutoff,
                     a_degree, c_degree, u_degree),
                )
                column_identity_count += 2

                if a_degree + c_degree + u_degree:
                    fixed_columns.append(g_column)
                    mixed_columns.extend((f_column, g_column))
                if c_degree + u_degree:
                    exact_columns.append(U_column)
                if a_degree + u_degree:
                    exact_columns.append(V_column)
                if a_degree + c_degree:
                    exact_columns.append(
                        make_column(branches, rows, cutoff, W_terms)
                    )

    monomial_count = comb(cutoff + 4, 3)
    require(len(rows) == 3 * comb(cutoff + 2, 2),
            ("row census", root, cutoff, len(rows)))
    require(len(fixed_columns) == monomial_count - 1,
            ("fixed column census", root, cutoff, len(fixed_columns)))
    require(len(mixed_columns) == 2 * (monomial_count - 1),
            ("mixed column census", root, cutoff, len(mixed_columns)))
    require(len(exact_columns) == 3 * (monomial_count - (cutoff + 2)),
            ("exact column census", root, cutoff, len(exact_columns)))
    require(column_identity_count == 2 * monomial_count,
            ("column identity census", root, cutoff, column_identity_count))
    return rows, fixed_columns, mixed_columns, exact_columns


def matrix_from_columns(columns, row_count):
    require(columns, "empty matrix bank")
    require(all(len(column) == row_count for column in columns),
            "ragged matrix bank")
    entries = [
        columns[column][row]
        for row in range(row_count)
        for column in range(len(columns))
    ]
    return nmod_mat(row_count, len(columns), entries, PRIME)


def stable_target(rows, degree, scalar=1):
    return [
        scalar % PRIME if stable == degree and source == 0 else 0
        for stable, point, source in rows
    ]


def left_relations(matrix, row_count, rank):
    relation_matrix, nullity = matrix.transpose().nullspace()
    require(nullity == row_count - rank,
            ("left nullity", nullity, row_count, rank))
    return relation_matrix, nullity


def belongs_to_image(target, relation_matrix, nullity):
    for relation in range(nullity):
        response = sum(
            int(relation_matrix[row, relation]) * target[row]
            for row in range(len(target))
        ) % PRIME
        if response:
            return False
    return True


def cutoff_audit(root, cutoff):
    rows, fixed_columns, mixed_columns, exact_columns = build_banks(root, cutoff)
    row_count = len(rows)
    fixed_matrix = matrix_from_columns(fixed_columns, row_count)
    mixed_matrix = matrix_from_columns(mixed_columns, row_count)
    exact_matrix = matrix_from_columns(exact_columns, row_count)
    fixed_rank = fixed_matrix.rank()
    mixed_rank = mixed_matrix.rank()
    exact_rank = exact_matrix.rank()
    require(fixed_rank == row_count - (cutoff + 1),
            ("fixed rank", root, cutoff, fixed_rank, row_count))
    require(mixed_rank == row_count - 4,
            ("mixed rank", root, cutoff, mixed_rank, row_count))
    require(exact_rank == row_count - 4,
            ("exact rank", root, cutoff, exact_rank, row_count))

    mixed_relations, mixed_nullity = left_relations(
        mixed_matrix, row_count, mixed_rank,
    )
    stable_memberships = tuple(
        belongs_to_image(stable_target(rows, degree), mixed_relations, mixed_nullity)
        for degree in range(cutoff + 1)
    )
    require(all(stable_memberships),
            ("stable RHS outside mixed image", root, cutoff, stable_memberships))
    constant_membership = stable_memberships[0]
    simple_zero_rhs_membership = belongs_to_image(
        stable_target(rows, 1, -24), mixed_relations, mixed_nullity,
    )
    require(constant_membership and simple_zero_rhs_membership,
            ("constant or -24u membership", root, cutoff))

    return (
        root,
        cutoff,
        row_count,
        len(fixed_columns),
        fixed_rank,
        row_count - fixed_rank,
        len(mixed_columns),
        mixed_rank,
        row_count - mixed_rank,
        len(exact_columns),
        exact_rank,
        row_count - exact_rank,
        int(constant_membership),
        int(simple_zero_rhs_membership),
        sum(stable_memberships),
    )


def main():
    started = perf_counter()
    require(all(quartic_value(root) == 0 for root in ROOTS), "bad split roots")
    require(tuple(rho_reduction(root) for root in ROOTS) == (8, 85, 12, 135),
            "rho reduction mismatch")
    for root in ROOTS:
        polynomial = q_polynomial(root)
        values = tuple(polynomial_value(polynomial, point) for point in POINTS)
        derivatives = tuple(
            polynomial_derivative_value(polynomial, point) for point in POINTS
        )
        require(values == tuple(value % PRIME for value in (-3, rational(-3, 4), -3)),
                ("Q values", root, values))
        require(derivatives == tuple(value % PRIME for value in
                                     (rational(-9, 2), 1, rational(9, 2))),
                ("Q derivatives", root, derivatives))

    triangle_rows = []
    for root in ROOTS:
        rho, delta, area, responses = triangle_and_carrier_audit(root)
        triangle_rows.append((root, rho, delta, area) + responses)

    rank_rows = []
    for root in ROOTS:
        for cutoff in CUTOFFS:
            require(perf_counter() - started < MAX_RUNTIME_SECONDS,
                    ("runtime gate before case", root, cutoff))
            rank_rows.append(cutoff_audit(root, cutoff))
            require(perf_counter() - started < MAX_RUNTIME_SECONDS,
                    ("runtime gate after case", root, cutoff))

    serialization = "\n".join(
        ",".join(str(value) for value in row)
        for row in triangle_rows + rank_rows
    )
    table_hash = sha256(serialization.encode()).hexdigest()
    require(table_hash == EXPECTED_TABLE_SHA256,
            ("deterministic table hash", table_hash, EXPECTED_TABLE_SHA256))

    print("THM-4060 MIXED/EXACT -- STANDALONE FINITE-FIELD HOSTILE AUDIT")
    print("field=F_137;split_roots=(44,82,92,134);cutoffs=4..12")
    print("Q_values=(-3,-3/4,-3);Q_derivatives=(-9/2,1,9/2);all_roots=True")
    print("source=affine_H(u)=u;coordinates=(A=a+3/4,c,u);Jac(A,c)=-3")
    print("mixed_form=-4*df^dc+dA^dg;D=g_c-4*f_A")
    print("exact_form=d(U*dA+V*dc+W*du);mixed_to_exact_columns=(g=-U,f=-V/4)")
    print("owner_period(D)=-4*d_A(Pi(f));f=A^k*u;normalize_by=delta*A^(k+9)")
    for row in triangle_rows:
        root, rho, delta, area = row[:4]
        print(
            f"triangle;alpha={root};rho={rho};delta={delta};area={area};"
            "common_vertex_through_A4=1;carrier_k=0..12;"
            "normalized_owner_response=2*(k+10)*rho"
        )
    for row in rank_rows:
        (
            root,
            cutoff,
            row_count,
            fixed_columns,
            fixed_rank,
            fixed_cokernel,
            mixed_columns,
            mixed_rank,
            mixed_cokernel,
            exact_columns,
            exact_rank,
            exact_cokernel,
            constant_membership,
            simple_zero_membership,
            stable_membership_count,
        ) = row
        print(
            f"rank;alpha={root};N={cutoff};rows={row_count};"
            f"fixed_cols={fixed_columns};fixed_rank={fixed_rank};"
            f"fixed_coker={fixed_cokernel};mixed_cols={mixed_columns};"
            f"mixed_rank={mixed_rank};mixed_coker={mixed_cokernel};"
            f"exact_cols={exact_columns};exact_rank={exact_rank};"
            f"exact_coker={exact_cokernel};constant_in_mixed={constant_membership};"
            f"minus24u_in_mixed={simple_zero_membership};"
            f"pure_u_memberships={stable_membership_count}/{cutoff + 1}"
        )
    print("all_mixed_columns_identified_with_exact_UV_columns=1")
    print("mixed_image_equals_exact_form_image=1")
    print("runtime_gate_seconds=120;runtime_gate=PASS")
    print("table_sha256=" + table_hash)
    print("scope=finite_field_hostile_audit_not_characteristic_zero_proof")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
