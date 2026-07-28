#!/usr/bin/env python3
"""Exact audit for the tree-incidence A/D Weyl clutch in THM-2770.

The companion checks reduced-incidence unimodularity on every recursive
parent array through eight vertices and performs the four-vertex star/path
factor, fan, Farey-flank, and balanced-coefficient calculations exactly.
It uses explicit exceptions and no truth-bearing Python assertions.
"""

from itertools import permutations, product
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def parent_arrays(vertex_count):
    if vertex_count == 1:
        yield (-1,)
        return
    for tail in product(*(range(child) for child in range(1, vertex_count))):
        yield (-1,) + tail


def matmul(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0]) if right else 0
    require(all(len(row) == middle for row in left),
            "matrix dimensions stopped matching")
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(middle))
              for j in range(columns))
        for i in range(rows)
    )


def reduced_incidence_and_paths(parents):
    vertex_count = len(parents)
    edge_count = vertex_count - 1
    incidence = []
    for child in range(1, vertex_count):
        row = [0] * edge_count
        row[child - 1] = 1
        parent = parents[child]
        if parent:
            row[parent - 1] = -1
        incidence.append(tuple(row))

    paths = []
    for vertex in range(1, vertex_count):
        row = [0] * edge_count
        current = vertex
        while current:
            row[current - 1] = 1
            current = parents[current]
        paths.append(tuple(row))
    return tuple(incidence), tuple(paths)


def identity(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def add_term(polynomial, exponent, coefficient):
    if coefficient == 0:
        return
    polynomial[exponent] = polynomial.get(exponent, 0) + coefficient
    if polynomial[exponent] == 0:
        del polynomial[exponent]


def multiply(left, right):
    variable_count = len(next(iter(left)))
    out = {}
    for first_exponent, first_coefficient in left.items():
        for second_exponent, second_coefficient in right.items():
            exponent = tuple(
                first_exponent[index] + second_exponent[index]
                for index in range(variable_count)
            )
            add_term(out, exponent, first_coefficient * second_coefficient)
    return out


def linear_form(coefficients):
    variable_count = len(coefficients)
    out = {}
    for index, coefficient in enumerate(coefficients):
        if coefficient:
            exponent = tuple(int(index == other) for other in range(variable_count))
            out[exponent] = coefficient
    return out


def product_of_forms(forms):
    variable_count = len(forms[0])
    out = {tuple([0] * variable_count): 1}
    for form in forms:
        out = multiply(out, linear_form(form))
    return out


def edge_forms(kind):
    if kind == "star":
        vertex_forms = (
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (-1, 1, 0), (-1, 0, 1), (0, -1, 1),
        )
    elif kind == "path":
        vertex_forms = (
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 1, 0), (0, 1, 1), (1, 1, 1),
        )
    else:
        raise RuntimeError("unknown four-vertex tree")
    d3_forms = (
        (-1, 1, 0), (1, 1, 0),
        (-1, 0, 1), (1, 0, 1),
        (0, -1, 1), (0, 1, 1),
    )
    return vertex_forms + d3_forms


def canonical_normal(vector):
    common = 0
    for value in vector:
        common = gcd(common, abs(value))
    require(common > 0, "zero hyperplane normal appeared")
    reduced = tuple(value // common for value in vector)
    first = next(value for value in reduced if value)
    if first < 0:
        reduced = tuple(-value for value in reduced)
    return reduced


def direct_vertex_polynomial(parents):
    vertex_count = len(parents)
    forms = []
    for first in range(vertex_count):
        for second in range(first + 1, vertex_count):
            form = [0] * vertex_count
            form[first] = -1
            form[second] = 1
            forms.append(tuple(form))
    for first in range(1, vertex_count):
        for second in range(first + 1, vertex_count):
            y_first = [0] * vertex_count
            y_first[first] = 1
            y_first[parents[first]] = -1
            y_second = [0] * vertex_count
            y_second[second] = 1
            y_second[parents[second]] = -1
            forms.append(tuple(y_second[i] - y_first[i]
                               for i in range(vertex_count)))
            forms.append(tuple(y_second[i] + y_first[i]
                               for i in range(vertex_count)))
    return product_of_forms(tuple(forms))


def translated_balanced_contributions(polynomial):
    rows = []
    for exponent, coefficient in sorted(polynomial.items()):
        if min(exponent) < 3:
            continue
        # y_i=x_i-x_0.  The total excess over x_i^3 is three, so every
        # contribution has the common sign (-1)^3.
        contribution = -coefficient
        for value in exponent:
            contribution *= comb(value, 3)
        if contribution:
            rows.append((exponent, coefficient, contribution))
    return tuple(rows)


def determinant_three(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def main():
    incidence_controls = 0
    for vertex_count in range(2, 9):
        for parents in parent_arrays(vertex_count):
            incidence, paths = reduced_incidence_and_paths(parents)
            unit = identity(vertex_count - 1)
            require(matmul(incidence, paths) == unit,
                    "reduced incidence lost its path-sum inverse")
            require(matmul(paths, incidence) == unit,
                    "path-sum matrix lost its incidence inverse")
            require(all(incidence[i][i] == 1 for i in range(vertex_count - 1))
                    and all(incidence[i][j] == 0
                            for i in range(vertex_count - 1)
                            for j in range(i + 1, vertex_count - 1)),
                    "reduced incidence stopped being unit lower triangular")
            incidence_controls += 1

    for vertex_count in range(2, 13):
        a_roots = vertex_count * (vertex_count - 1) // 2
        d_roots = (vertex_count - 1) * (vertex_count - 2)
        expected = (vertex_count - 1) * (3 * vertex_count - 4) // 2
        require(a_roots + d_roots == expected,
                "A/D positive-root degree identity changed")

    b3 = set()
    for index in range(3):
        vector = [0, 0, 0]
        vector[index] = 1
        b3.add(canonical_normal(tuple(vector)))
    for first in range(3):
        for second in range(first + 1, 3):
            for sign in (-1, 1):
                vector = [0, 0, 0]
                vector[first] = 1
                vector[second] = sign
                b3.add(canonical_normal(tuple(vector)))
    diagonal = canonical_normal((1, 1, 1))
    star_support = {canonical_normal(form) for form in edge_forms("star")}
    path_support = {canonical_normal(form) for form in edge_forms("path")}
    require(star_support == b3, "star clutch stopped supporting B3")
    require(path_support == b3 | {diagonal},
            "path clutch stopped being B3 plus one diagonal")

    # The half-Hadamard D3=A3 map.  Pair differences of the four displayed
    # numerator rows, divided by two, are exactly the six D3 normals.
    even_sign_rows = (
        (1, 1, 1), (1, -1, -1),
        (-1, 1, -1), (-1, -1, 1),
    )
    a3_images = set()
    for first in range(4):
        for second in range(first + 1, 4):
            numerator = tuple(even_sign_rows[first][i]
                              - even_sign_rows[second][i] for i in range(3))
            require(all(value % 2 == 0 for value in numerator),
                    "half-Hadamard root stopped being integral")
            a3_images.add(canonical_normal(tuple(value // 2 for value in numerator)))
    d3 = {canonical_normal(form) for form in edge_forms("star")[6:]}
    require(a3_images == d3 and len(d3) == 6,
            "D3=A3 root-hyperplane identification changed")

    cut_chambers = 0
    cut_controls = 0
    for absolute_order in permutations(range(3)):
        largest = absolute_order[-1]
        for signs in product((-1, 1), repeat=3):
            cut = (signs[absolute_order[0]] == signs[absolute_order[1]]
                   and signs[largest] == -signs[absolute_order[0]])
            if not cut:
                continue
            cut_chambers += 1
            magnitudes = [0, 0, 0]
            for rank, coordinate in enumerate(absolute_order, start=1):
                magnitudes[coordinate] = rank
            point = tuple(signs[i] * magnitudes[i] for i in range(3))
            require(sum(point) == 0,
                    "a claimed cut B3 chamber missed the diagonal")
            require(all(value and abs(point[i]) != abs(point[j])
                        for i, value in enumerate(point)
                        for j in range(i + 1, 3)),
                    "diagonal test point hit a B3 wall")
            cut_controls += 1
    require(cut_chambers == 12 and cut_controls == 12,
            "the diagonal stopped cutting exactly twelve B3 chambers")

    e_a = (1, 0, 0)
    e_b = (0, 1, 0)
    e_c = (0, 0, 1)
    farey = (1, 0, 1)
    farey_determinants = (
        determinant_three((e_a, e_b, farey)),
        determinant_three((e_c, e_b, farey)),
    )
    require(tuple(abs(value) for value in farey_determinants) == (1, 1),
            "the path diagonal stopped giving a unimodular Farey split")

    star_edge_polynomial = product_of_forms(edge_forms("star"))
    path_edge_polynomial = product_of_forms(edge_forms("path"))
    star_rows = translated_balanced_contributions(star_edge_polynomial)
    path_rows = translated_balanced_contributions(path_edge_polynomial)
    star_balanced = sum(row[2] for row in star_rows)
    path_balanced = sum(row[2] for row in path_rows)
    require(star_balanced == 120 and path_balanced == 0,
            "four-vertex balanced coefficients changed")
    require(sorted(row[2] for row in star_rows) == [-40] * 3 + [40] * 6,
            "star coefficient contribution table changed")
    require(sorted(row[2] for row in path_rows) == [-80, -40, -20, 20, 40, 80],
            "path cancellation table changed")

    central_by_shape = {"star": [], "path": []}
    for parents in parent_arrays(4):
        degree = [0] * 4
        for child in range(1, 4):
            degree[child] += 1
            degree[parents[child]] += 1
        shape = "star" if max(degree) == 3 else "path"
        coefficient = direct_vertex_polynomial(parents).get((3, 3, 3, 3), 0)
        central_by_shape[shape].append(coefficient)
    require(sorted(abs(value) for value in central_by_shape["star"]) == [120, 120]
            and central_by_shape["path"] == [0, 0, 0, 0],
            "the labelled four-vertex central census changed")

    path_labels = (0, 3, 1, 2)
    path_differences = tuple(
        abs(path_labels[index + 1] - path_labels[index]) for index in range(3)
    )
    require(path_differences == (3, 2, 1),
            "the graceful P4 stopping control changed")

    print("TREE INCIDENCE A/D WEYL CLUTCH AUDIT")
    print(f"reduced_incidence_parent_arrays_n2_to_n8={incidence_controls}")
    print("reduced_incidence=unimodular path_sum_inverse=exact")
    print("degree=|A_(n-1)^+|+|D_(n-1)^+|=(n-1)(3n-4)/2")
    print("D3_to_A3_half_hadamard_root_hyperplanes=6")
    print("star_support=B3 hyperplanes=9 chambers=48")
    print("path_support=B3_plus_diagonal hyperplanes=10")
    print(f"path_diagonal_cut_chambers={cut_chambers} path_chambers={48 + cut_chambers}")
    print(f"farey_flank_determinants={farey_determinants}")
    print(f"star_balanced_contributions={len(star_rows)} coefficient={star_balanced}")
    print(f"path_balanced_contributions={len(path_rows)} coefficient={path_balanced}")
    print("four_vertex_recursive_shapes=2_star_coeff_abs120,4_path_coeff0")
    print("P4_graceful_control=labels(0,3,1,2)_differences(3,2,1)")
    print("SCOPE: Weyl/fan/coefficient theorem; no graceful-tree closure")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
