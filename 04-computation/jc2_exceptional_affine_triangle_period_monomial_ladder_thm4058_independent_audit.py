#!/usr/bin/env python3
"""Independent split-prime audit of the THM-4058 period and obstruction ladder."""

from hashlib import sha256
from math import comb


PRIME = 137
ROOTS = (44, 82, 92, 134)
POINTS = (-1, 0, 1)
BRANCH_S = (2, 0, -2)
RHO = {44: 8, 82: 85, 92: 12, 134: 135}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inv(value):
    return pow(value % PRIME, -1, PRIME)


def rat(numerator, denominator=1):
    return numerator % PRIME * inv(denominator) % PRIME


def add(left, right):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = (answer.get(monomial, 0) + coefficient) % PRIME
        if not answer[monomial]:
            del answer[monomial]
    return answer


def scale(value, scalar):
    scalar %= PRIME
    return {monomial: coefficient * scalar % PRIME
            for monomial, coefficient in value.items()
            if coefficient * scalar % PRIME}


def mul(left, right, cutoff):
    answer = {}
    for (i, j), a in left.items():
        for (u, v), b in right.items():
            if i + j + u + v <= cutoff:
                key = (i + u, j + v)
                answer[key] = (answer.get(key, 0) + a * b) % PRIME
    return {key: value for key, value in answer.items() if value}


def power(value, exponent, cutoff):
    answer = {(0, 0): 1}
    while exponent:
        if exponent & 1:
            answer = mul(answer, value, cutoff)
        value = mul(value, value, cutoff)
        exponent //= 2
    return answer


def diff(value, axis):
    answer = {}
    for (i, j), coefficient in value.items():
        exponent = i if axis == 0 else j
        if exponent:
            key = (i - (axis == 0), j - (axis == 1))
            answer[key] = coefficient * exponent % PRIME
    return answer


def series_inverse(value, cutoff):
    constant = value.get((0, 0), 0)
    require(constant != 0, "nonunit series")
    inverse_constant = inv(constant)
    remainder = add(scale(value, inverse_constant), {(0, 0): -1 % PRIME})
    answer = {}
    term = {(0, 0): 1}
    for degree in range(cutoff + 1):
        answer = add(answer, scale(term, -1 if degree % 2 else 1))
        term = mul(term, remainder, cutoff)
    return scale(answer, inverse_constant)


def uadd(left, right):
    answer = dict(left)
    for degree, coefficient in right.items():
        answer[degree] = (answer.get(degree, 0) + coefficient) % PRIME
        if not answer[degree]:
            del answer[degree]
    return answer


def uscale(value, scalar):
    scalar %= PRIME
    return {degree: coefficient * scalar % PRIME
            for degree, coefficient in value.items()
            if coefficient * scalar % PRIME}


def ushift(value, amount):
    return {degree + amount: coefficient for degree, coefficient in value.items()}


def around(value, point, cutoff):
    answer = {}
    for degree, coefficient in value.items():
        for local_degree in range(min(degree, cutoff) + 1):
            key = (local_degree, 0)
            answer[key] = (
                answer.get(key, 0)
                + coefficient * comb(degree, local_degree)
                * point ** (degree - local_degree)
            ) % PRIME
    return {key: coefficient for key, coefficient in answer.items() if coefficient}


Q_BASE = {
    0: rat(-3, 4), 1: 1, 2: rat(-27, 4), 3: -2 % PRIME,
    4: rat(9, 2), 5: 1,
}
P_BASE = {2: 1, 4: -2 % PRIME, 6: 1}
Q6 = uadd(Q_BASE, uscale(P_BASE, rat(-259, 36)))
R1 = uadd(P_BASE, uscale(ushift(P_BASE, 2), -1))
R2 = uadd(uscale(P_BASE, 4), uscale(ushift(P_BASE, 1), -9))


def q_polynomial(root):
    theta = (
        rat(520, 9) * root * root
        - rat(1688, 81) * root
        - rat(5717, 729)
    ) % PRIME
    return uadd(Q6, uadd(uscale(R1, theta), uscale(R2, root)))


def compose_univariate(polynomial, value, cutoff):
    answer = {}
    for degree, coefficient in polynomial.items():
        answer = add(answer, scale(power(value, degree, cutoff), coefficient))
    return answer


def graph_branches(root, cutoff=6):
    qalpha = q_polynomial(root)
    variable_a = {(1, 0): 1}
    variable_c = {(0, 1): 1}
    a = {(0, 0): rat(-3, 4), (1, 0): 1}
    answer = []
    for base_s in BRANCH_S:
        correction = {}
        for _ in range(5):
            s = add({(0, 0): base_s % PRIME}, correction)
            equation = add(
                add(scale(s, 3), mul(a, power(s, 3, cutoff), cutoff)),
                scale(variable_c, -1),
            )
            derivative = add(
                {(0, 0): 3},
                scale(mul(a, power(s, 2, cutoff), cutoff), 3),
            )
            correction = add(
                correction,
                scale(mul(equation, series_inverse(derivative, cutoff), cutoff), -1),
            )
        s = add({(0, 0): base_s % PRIME}, correction)
        equation = add(
            add(scale(s, 3), mul(a, power(s, 3, cutoff), cutoff)),
            scale(variable_c, -1),
        )
        require(not equation, ("implicit branch residual", root, base_s))
        denominator = add({(0, 0): 1}, mul(a, power(s, 2, cutoff), cutoff))
        x = mul(s, series_inverse(denominator, cutoff), cutoff)
        q = mul(a, power(denominator, 2, cutoff), cutoff)
        graph = add(q, scale(compose_univariate(qalpha, x, cutoff), -1))
        answer.append(graph)
    require(variable_a == {(1, 0): 1}, "unused A-variable sentinel")
    return answer


def substitute_c(polynomial, c_series, cutoff):
    answer = {}
    for (a_degree, c_degree), coefficient in polynomial.items():
        term = scale(power(c_series, c_degree, cutoff), coefficient)
        for (degree, dummy), value in term.items():
            require(dummy == 0, "univariate substitution acquired c")
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
            answer[(degree, 0)] = -coefficient * inv(linear) % PRIME
    require(not substitute_c(difference, answer, cutoff), "pairwise root residual")
    return answer


def triangle_audit(root):
    cutoff = 6
    branches = graph_branches(root, cutoff)
    slopes = (rat(3, 4), rat(-1, 3), rat(-3, 4))
    for branch, slope in zip(branches, slopes):
        require(branch.get((1, 0), 0) == 1, ("A slope", root))
        require(branch.get((0, 1), 0) == slope, ("c slope", root))

    vertices = {}
    for i, j in ((0, 1), (1, 2), (2, 0)):
        c_value = implicit_c_root(add(branches[i], scale(branches[j], -1)), cutoff)
        w_left = substitute_c(branches[i], c_value, cutoff)
        w_right = substitute_c(branches[j], c_value, cutoff)
        require(w_left == w_right, ("intersection value", root, i, j))
        vertices[(i, j)] = (c_value, w_left)

    vertex_order = (vertices[(2, 0)], vertices[(0, 1)], vertices[(1, 2)])
    for coordinate in (0, 1):
        for degree in range(1, 5):
            values = [vertex[coordinate].get((degree, 0), 0)
                      for vertex in vertex_order]
            require(len(set(values)) == 1, ("common vertex jet", root, coordinate, degree))

    c_edges = []
    w_edges = []
    for start, end in zip(vertex_order, vertex_order[1:] + vertex_order[:1]):
        c_edges.append((end[0].get((5, 0), 0) - start[0].get((5, 0), 0)) % PRIME)
        w_edges.append((end[1].get((5, 0), 0) - start[1].get((5, 0), 0)) % PRIME)
    delta = c_edges[2]
    require(delta == rat(26, 15) * RHO[root] % PRIME, ("delta/rho", root, delta))
    expected_c = [delta * rat(5, 13) % PRIME,
                  delta * rat(-18, 13) % PRIME, delta]
    require(c_edges == expected_c, ("c-edge vector", root, c_edges))
    require(w_edges == [slope * edge % PRIME for slope, edge in zip(slopes, c_edges)],
            ("w-edge vector", root, w_edges))

    current_w = 0
    area = 0
    for c_edge, slope in zip(c_edges, slopes):
        area = (area + current_w * c_edge + rat(1, 2) * slope * c_edge * c_edge) % PRIME
        current_w = (current_w + slope * c_edge) % PRIME
    require(current_w == 0 and sum(c_edges) % PRIME == 0, ("polygon closure", root))
    require(area == rat(-15, 52) * delta * delta % PRIME, ("triangle area", root, area))
    require(area * inv(delta) % PRIME == rat(-1, 2) * RHO[root] % PRIME,
            ("reduced period", root))
    return delta, area


def fixed_a_columns(root, cutoff):
    rows = tuple(
        (stable, point, source)
        for stable in range(cutoff + 1)
        for source in range(cutoff + 1 - stable)
        for point in POINTS
    )
    qalpha = q_polynomial(root)
    local_cutoff = cutoff + 2
    packets = {}
    for point in POINTS:
        qjet = add(around(qalpha, point, local_cutoff), {(0, 1): 1})
        xjet = {(0, 0): point % PRIME, (1, 0): 1}
        denominator = add(
            {(0, 0): 1},
            mul(mul(xjet, xjet, local_cutoff), qjet, local_cutoff),
        )
        denominator_inverse = series_inverse(denominator, local_cutoff)
        ajet = mul(
            qjet,
            mul(denominator_inverse, denominator_inverse, local_cutoff),
            local_cutoff,
        )
        alocal = add(ajet, {(0, 0): rat(3, 4)})
        cjet = mul(
            mul(xjet, denominator, local_cutoff),
            add(denominator, {(0, 0): 2}),
            local_cutoff,
        )
        packets[point] = (
            [power(alocal, degree, cutoff + 1) for degree in range(cutoff + 2)],
            [power(cjet, degree, cutoff + 1) for degree in range(cutoff + 2)],
            diff(ajet, 0),
        )

    columns = []
    for a_degree in range(cutoff + 2):
        for c_degree in range(cutoff + 2 - a_degree):
            for w_degree in range(cutoff + 2 - a_degree - c_degree):
                if a_degree + c_degree + w_degree == 0:
                    continue
                branch_values = {}
                for point in POINTS:
                    a_powers, c_powers, a_x = packets[point]
                    density = {}
                    if c_degree:
                        derivative_c = mul(
                            a_powers[a_degree], c_powers[c_degree - 1], cutoff,
                        )
                        derivative_c = mul(
                            derivative_c, {(0, w_degree): 1}, cutoff,
                        )
                        density = add(density, scale(derivative_c, -3 * c_degree))
                    if w_degree:
                        derivative_w = mul(
                            a_powers[a_degree], c_powers[c_degree], cutoff,
                        )
                        derivative_w = mul(
                            derivative_w, {(0, w_degree - 1): 1}, cutoff,
                        )
                        density = add(
                            density,
                            scale(mul(a_x, derivative_w, cutoff), w_degree),
                        )
                    branch_values[point] = density
                columns.append([
                    branch_values[point].get((source, stable), 0)
                    for stable, point, source in rows
                ])
    require(len(columns) == comb(cutoff + 4, 3) - 1, ("column census", cutoff))
    return rows, columns


def rref(matrix):
    rows = [[entry % PRIME for entry in row] for row in matrix]
    pivots = []
    pivot_row = 0
    for column in range(len(rows[0])):
        selected = next((index for index in range(pivot_row, len(rows))
                         if rows[index][column]), None)
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scalar = inv(rows[pivot_row][column])
        rows[pivot_row] = [scalar * entry % PRIME for entry in rows[pivot_row]]
        for index in range(len(rows)):
            if index != pivot_row and rows[index][column]:
                scalar = rows[index][column]
                rows[index] = [
                    (entry - scalar * pivot) % PRIME
                    for entry, pivot in zip(rows[index], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivots


def relation_basis(columns):
    echelon, pivots = rref(columns)
    free = [column for column in range(len(columns[0])) if column not in pivots]
    answer = []
    for column in free:
        vector = [0] * len(columns[0])
        vector[column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -echelon[row][column] % PRIME
        answer.append(vector)
    return len(pivots), answer


def inverse_derivative_coefficients(r, cutoff):
    return {
        r * index: (index, (-1 if index % 2 else 1)
                    * comb((r + 1) * index, index) % PRIME)
        for index in range(cutoff // r + 1)
    }


def response_polynomials(rows, relations, r, cutoff):
    coefficients = inverse_derivative_coefficients(r, cutoff)
    answer = []
    for relation in relations:
        polynomial = {}
        for index, (stable, _, source) in enumerate(rows):
            if source == 0 and stable in coefficients:
                gamma_degree, coefficient = coefficients[stable]
                polynomial[gamma_degree] = (
                    polynomial.get(gamma_degree, 0)
                    + coefficient * relation[index]
                ) % PRIME
        polynomial = {degree: value for degree, value in polynomial.items() if value}
        if polynomial:
            answer.append(polynomial)
    return answer


quartic = lambda root: (
    72783360 * root ** 4 - 77822208 * root ** 3 - 28419741 * root ** 2
    + 7849770 * root - 1276420
) % PRIME
require(all(quartic(root) == 0 for root in ROOTS), "split roots")

table = []
triangle_rows = []
for root in ROOTS:
    delta, area = triangle_audit(root)
    triangle_rows.append((root, delta, area))
    cache = {}
    for cutoff in range(4, 16):
        rows, columns = fixed_a_columns(root, cutoff)
        rank, relations = relation_basis(columns)
        require(rank == len(rows) - (cutoff + 1), ("rank", root, cutoff, rank))
        require(len(relations) == cutoff + 1, ("relation dimension", root, cutoff))
        cache[cutoff] = (rows, relations, rank)
    for m in range(2, 13):
        r = m - 1
        pass_cutoff = r + 3
        fail_cutoff = r + 4
        pass_rows, pass_relations, pass_rank = cache[pass_cutoff]
        fail_rows, fail_relations, fail_rank = cache[fail_cutoff]
        pass_responses = response_polynomials(
            pass_rows, pass_relations, r, pass_cutoff,
        )
        fail_responses = response_polynomials(
            fail_rows, fail_relations, r, fail_cutoff,
        )
        expected = RHO[root] * comb(m, 2) % PRIME
        require(not pass_responses, ("premature response", root, m, pass_responses))
        require(fail_responses == [{1: expected}],
                ("first response", root, m, fail_responses, expected))
        require(fail_responses != [{1: (expected + 1) % PRIME}], "hostile response")
        table.append((root, m, pass_cutoff, pass_rank, fail_cutoff, fail_rank, expected))

serialization = "\n".join(
    ",".join(str(value) for value in row)
    for row in triangle_rows + table
)
table_hash = sha256(serialization.encode()).hexdigest()
require(table_hash == "1641250c558fbdeea9f2e2ccff1a515761d4af596f95ae5bce4587d12edc72c9",
        ("table hash", table_hash))

print("EXCEPTIONAL AFFINE TRIANGLE PERIOD / MONOMIAL LADDER -- INDEPENDENT MOD-137 AUDIT")
print("split_roots=(44,82,92,134);all_good=True")
print("branch_slopes=(3/4,-1/3,-3/4);pairwise_vertices_common_through_A4=True")
print("rho_reductions=(8,85,12,135);delta=(26/15)rho")
print("c_edge_ratios=(5/13,-18/13,1);area=(-15/52)delta^2")
print("reduced_period_w=-rho/2")
print("m_range=2..12;pass_cutoff=m+2;fail_cutoff=m+3")
print("all_pass_responses_zero=True;all_fail_responses=binom(m,2)*gamma*rho")
print("rank_formula=rows-(N+1);inverse_derivative_gamma_polynomials_checked=True")
print("table_sha256=" + table_hash)
print("scope=independent_split-prime_hostile_audit_not_characteristic-zero_proof")
print("RESULT=PASS")
