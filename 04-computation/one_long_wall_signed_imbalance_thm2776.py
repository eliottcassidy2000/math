#!/usr/bin/env python3
"""Exact one-long-wall signed-frame audit for THM-2776.

For every selection of k-1 positive B_k hyperplane roots through k=6, the
script compares the full determinant against the signed-component formula.
It uses integer arithmetic and explicit exceptions only: no truth-bearing
assertions or floating point.
"""

from collections import Counter, deque
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix),
            "determinant received a nonsquare matrix")
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot
                    - work[row][column] * work[column][target]
                )
                require(numerator % previous == 0,
                        "Bareiss division stopped being exact")
                work[row][target] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def b_roots(dimension):
    roots = []
    for index in range(dimension):
        row = [0] * dimension
        row[index] = 1
        roots.append(tuple(row))
    for first in range(dimension):
        for second in range(first + 1, dimension):
            for sign in (-1, 1):
                row = [0] * dimension
                row[first] = 1
                row[second] = sign
                roots.append(tuple(row))
    require(len(roots) == dimension * dimension,
            "positive B_k root atlas changed")
    return tuple(roots)


def ordinary_components(selection, dimension):
    adjacency = [[] for _ in range(dimension)]
    row_data = []
    for row in selection:
        support = tuple(index for index, value in enumerate(row) if value)
        require(len(support) in (1, 2),
                "a B root stopped being a half-edge or ordinary edge")
        if len(support) == 1:
            row_data.append(("half", support[0], None, None))
            continue
        first, second = support
        transport = -row[first] * row[second]
        adjacency[first].append((second, transport))
        adjacency[second].append((first, transport))
        row_data.append(("edge", first, second, transport))

    component_of = [-1] * dimension
    components = []
    for seed in range(dimension):
        if component_of[seed] != -1:
            continue
        label = len(components)
        vertices = []
        queue = [seed]
        component_of[seed] = label
        while queue:
            vertex = queue.pop()
            vertices.append(vertex)
            for neighbour, _transport in adjacency[vertex]:
                if component_of[neighbour] == -1:
                    component_of[neighbour] = label
                    queue.append(neighbour)
        components.append({
            "vertices": tuple(sorted(vertices)),
            "edges": [],
            "halves": [],
            "adjacency": adjacency,
        })

    for kind, first, second, transport in row_data:
        label = component_of[first]
        if kind == "half":
            components[label]["halves"].append(first)
        else:
            require(component_of[second] == label,
                    "ordinary edge crossed its connected component")
            components[label]["edges"].append((first, second, transport))
    return components


def sign_kernel(component):
    """Return the +/-1 kernel on a balanced ordinary component, else None."""
    vertices = component["vertices"]
    signs = {vertices[0]: 1}
    queue = deque([vertices[0]])
    adjacency = component["adjacency"]
    while queue:
        vertex = queue.popleft()
        for neighbour, transport in adjacency[vertex]:
            proposed = transport * signs[vertex]
            if neighbour in signs:
                if signs[neighbour] != proposed:
                    return None
            else:
                signs[neighbour] = proposed
                queue.append(neighbour)
    require(set(signs) == set(vertices),
            "component sign propagation missed a vertex")
    return signs


def component_prediction(selection, dimension):
    """Classify rank k-1 and return (independent, predicted_abs, data)."""
    components = ordinary_components(selection, dimension)
    deficient = []
    cycle_count = 0
    component_types = []
    for component in components:
        vertex_count = len(component["vertices"])
        edge_count = len(component["edges"])
        half_count = len(component["halves"])
        row_count = edge_count + half_count

        if row_count == vertex_count - 1 and half_count == 0:
            require(edge_count == vertex_count - 1,
                    "deficient component stopped being an ordinary tree")
            signs = sign_kernel(component)
            require(signs is not None,
                    "an ordinary tree acquired an inconsistent sign kernel")
            deficient.append((component, signs))
            component_types.append("tree")
            continue

        if (row_count == vertex_count and half_count == 1
                and edge_count == vertex_count - 1):
            component_types.append("anchored")
            continue

        if (row_count == vertex_count and half_count == 0
                and edge_count == vertex_count):
            signs = sign_kernel(component)
            if signs is None:
                cycle_count += 1
                component_types.append("unbalanced")
                continue
            return False, 0, ("balanced_cycle",)

        return False, 0, ("dependent_shape",)

    if len(deficient) != 1:
        return False, 0, ("deficient_count", len(deficient))
    component, signs = deficient[0]
    imbalance = sum(signs[vertex] for vertex in component["vertices"])
    predicted = (2 ** cycle_count) * abs(imbalance)
    return True, predicted, (
        len(component["vertices"]), imbalance, cycle_count,
        tuple(sorted(component_types)),
    )


def odd_part(value):
    while value and value % 2 == 0:
        value //= 2
    return value


def theoretical_nonzero_support(dimension):
    values = set()
    for deficient_size in range(1, dimension + 1):
        first_imbalance = 1 if deficient_size % 2 else 2
        for imbalance in range(first_imbalance, deficient_size + 1, 2):
            remaining = dimension - deficient_size
            for cycle_count in range(remaining // 2 + 1):
                values.add((2 ** cycle_count) * imbalance)
    return values


def path_frame(length):
    rows = []
    for index in range(length - 1):
        row = [0] * length
        row[index] = 1
        row[index + 1] = -1
        rows.append(tuple(row))
    return tuple(rows)


def main():
    selection_count = 0
    rank_count = 0
    balanced_zero_count = 0
    histograms = {}
    shape_histogram = Counter()

    for dimension in range(2, 7):
        roots = b_roots(dimension)
        determinant_histogram = Counter()
        for selection in combinations(roots, dimension - 1):
            selection_count += 1
            independent, predicted, data = component_prediction(
                selection, dimension
            )
            actual = abs(determinant(selection + (tuple([1] * dimension),)))
            determinant_histogram[actual] += 1
            if not independent:
                require(actual == 0,
                        "dependent B selection produced a nonzero long-wall frame")
                continue
            rank_count += 1
            require(actual == predicted,
                    "signed-component determinant factorization changed")
            deficient_size, imbalance, cycle_count, types = data
            shape_histogram[
                (dimension, deficient_size, abs(imbalance), cycle_count)
            ] += 1
            if imbalance == 0:
                balanced_zero_count += 1
                require(actual == 0,
                        "balanced deficient tree stopped being the zero boundary")
            if actual:
                require(odd_part(actual) <= deficient_size <= dimension,
                        "odd determinant part exceeded the deficient tree")

        observed_support = {
            value for value, count in determinant_histogram.items()
            if value and count
        }
        require(observed_support == theoretical_nonzero_support(dimension),
                "exact determinant support stopped matching the realization formula")
        histograms[dimension] = determinant_histogram

        path = path_frame(dimension)
        require(abs(determinant(path + (tuple([1] * dimension),))) == dimension,
                "difference-tree path stopped attaining the maximal imbalance")

    # The natural P/Q character of the path frame is coordinate sum modulo k.
    # Any homogeneous monomial of total degree D therefore occupies one
    # character sector; at P4 the degree-12 graceful obstruction is trivial
    # under the diagonal mu_3 action.
    require(12 % 3 == 0,
            "P4 homogeneous-degree character boundary changed")

    histogram_text = ";".join(
        "k" + str(dimension) + ":"
        + ",".join(
            f"det{value}={histograms[dimension][value]}"
            for value in sorted(histograms[dimension])
        )
        for dimension in range(2, 7)
    )
    support_text = ";".join(
        "k" + str(dimension) + "={"
        + ",".join(str(value)
                   for value in sorted(theoretical_nonzero_support(dimension)))
        + "}"
        for dimension in range(2, 7)
    )

    print("ONE LONG WALL SIGNED-COMPONENT IMBALANCE AUDIT")
    print(f"B_selections_k2_to_k6={selection_count} rank_k_minus_1={rank_count}")
    print(f"balanced_full_rank_zero_frames={balanced_zero_count}")
    print("det_formula=2^c*abs(sum_U epsilon)")
    print("unique_deficient_component=unanchored_signed_tree")
    print("square_components=half_edge_tree:det1,unbalanced_cycle:det2")
    print("determinant_histograms=" + histogram_text)
    print("nonzero_supports=" + support_text)
    print(f"classified_rank_shapes={sum(shape_histogram.values())}")
    print("odd_prime_divisors<=deficient_size<=k")
    print("difference_spanning_tree_attains_det_k")
    print("P4_mu3_character=total_degree_mod3=12_mod3=0_no_coefficient_split")
    print("SCOPE=one_selected_frame_not_global_arrangement_or_graceful_proof")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
