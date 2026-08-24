#!/usr/bin/env python3
"""Exact finite analogue of THM-3990's componentwise Poisson obstruction.

Universe
--------
* every labeled simple graph on 1 <= n <= 5 vertices (1,099 graphs);
* both Laplacian signs L and -L;
* the integral Smith invariants of every simple graph on n <= 3, including
  the triangle hostile (1,3,0);
* four additional disconnected rationally edge-weighted real/rational
  controls, kept separate from the integral sidecar.

For each graph, over Q, the checker verifies

    rank(L) = n-c,
    rank([1 | L]) = n-c+1,
    rank([1_C1 | ... | 1_Cc | L]) = n,

where c is the number of connected components.  It also verifies directly
that component sums annihilate L, that a target with equal component averages
is in im[1|L], and that the hostile target 1_C1 (when c > 1) has nonzero total
sum but is not in im[1|L], while it is in the componentwise repair image.

This is an exact rational companion, not a proof of the smooth theorem.  The
smooth proof uses Green's identity, the kernel description for the
Laplace--Beltrami operator, the Fredholm alternative, and elliptic regularity.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd


GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise AssertionError(label)


def rank(matrix):
    """Exact row rank over Q; matrix is a list of rows."""
    if not matrix:
        return 0
    a = [[Fraction(x) for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    r = 0
    for col in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        pivot_value = a[r][col]
        a[r] = [x / pivot_value for x in a[r]]
        for i in range(rows):
            if i != r and a[i][col]:
                factor = a[i][col]
                a[i] = [a[i][j] - factor * a[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def determinant(matrix):
    """Exact determinant for the tiny integer Smith controls."""
    n = len(matrix)
    if n == 0:
        return 1
    if n == 1:
        return matrix[0][0]
    total = 0
    for j, value in enumerate(matrix[0]):
        minor = [row[:j] + row[j + 1:] for row in matrix[1:]]
        total += (-1) ** j * value * determinant(minor)
    return total


def smith_invariants_by_minors(matrix):
    """Smith diagonal from determinantal divisors; intended for n <= 3."""
    n = len(matrix)
    r = rank(matrix)
    previous_divisor = 1
    invariants = []
    for k in range(1, r + 1):
        divisor = 0
        for row_indices in combinations(range(n), k):
            for col_indices in combinations(range(n), k):
                minor = [[int(matrix[i][j]) for j in col_indices]
                         for i in row_indices]
                divisor = gcd(divisor, abs(determinant(minor)))
        require(divisor > 0 and divisor % previous_divisor == 0,
                f"Smith determinantal divisor k={k}")
        invariants.append(divisor // previous_divisor)
        previous_divisor = divisor
    invariants.extend([0] * (n - r))
    return tuple(invariants)


def columns_to_rows(columns, n):
    return [[column[i] for column in columns] for i in range(n)]


def append_column(matrix, column):
    return [row + [column[i]] for i, row in enumerate(matrix)]


def mat_vec(matrix, vector):
    return [sum((row[j] * vector[j] for j in range(len(vector))), Fraction(0))
            for row in matrix]


def in_column_space(matrix, target):
    return rank(matrix) == rank(append_column(matrix, target))


def components(n, weighted_edges):
    adjacency = [[] for _ in range(n)]
    for u, v, weight in weighted_edges:
        if weight:
            adjacency[u].append(v)
            adjacency[v].append(u)
    answer = []
    unseen = set(range(n))
    while unseen:
        root = min(unseen)
        stack = [root]
        unseen.remove(root)
        component = []
        while stack:
            u = stack.pop()
            component.append(u)
            for v in adjacency[u]:
                if v in unseen:
                    unseen.remove(v)
                    stack.append(v)
        answer.append(tuple(sorted(component)))
    return tuple(answer)


def laplacian(n, weighted_edges, sign=1):
    matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for u, v, weight in weighted_edges:
        w = Fraction(weight)
        matrix[u][u] += w
        matrix[v][v] += w
        matrix[u][v] -= w
        matrix[v][u] -= w
    return [[Fraction(sign) * x for x in row] for row in matrix]


def indicator(n, component):
    support = set(component)
    return [Fraction(int(i in support)) for i in range(n)]


def component_averages(vector, comps):
    return tuple(sum((vector[i] for i in comp), Fraction(0)) / len(comp)
                 for comp in comps)


def audit_graph(n, weighted_edges, tag):
    comps = components(n, weighted_edges)
    c = len(comps)
    component_columns = [indicator(n, comp) for comp in comps]
    ones = [Fraction(1) for _ in range(n)]

    for sign in (1, -1):
        L = laplacian(n, weighted_edges, sign)
        scalar_matrix = columns_to_rows([ones] + [
            [L[i][j] for i in range(n)] for j in range(n)
        ], n)
        repaired_matrix = columns_to_rows(component_columns + [
            [L[i][j] for i in range(n)] for j in range(n)
        ], n)

        require(rank(L) == n - c, f"{tag}: rank L, sign={sign}")
        require(rank(scalar_matrix) == n - c + 1,
                f"{tag}: scalar image rank, sign={sign}")
        require(rank(repaired_matrix) == n,
                f"{tag}: repaired image rank, sign={sign}")
        require(n - rank(scalar_matrix) == c - 1,
                f"{tag}: quotient dimension, sign={sign}")

        # Every component-sum functional annihilates every Laplacian column.
        for comp in comps:
            for j in range(n):
                require(sum((L[i][j] for i in comp), Fraction(0)) == 0,
                        f"{tag}: component sum annihilates L, sign={sign}")

        # Positive control: mu*1 + Lz has the same normalized average on all
        # components and lies in the scalar image, for a nonconstant exact z.
        mu = Fraction(7, 3)
        z = [Fraction((i + 1) * (i + 2), i + 3) for i in range(n)]
        Lz = mat_vec(L, z)
        control = [mu + value for value in Lz]
        require(component_averages(control, comps) == (mu,) * c,
                f"{tag}: equal-average control, sign={sign}")
        require(in_column_space(scalar_matrix, control),
                f"{tag}: scalar control solvable, sign={sign}")

        if c > 1:
            # Hostile: total sum is nonzero, but normalized component means
            # are (1,0,...,0), so no single scalar can repair it.
            hostile = component_columns[0]
            hostile_means = component_averages(hostile, comps)
            require(sum(hostile, Fraction(0)) > 0,
                    f"{tag}: hostile total sum nonzero, sign={sign}")
            require(hostile_means[0] == 1 and all(x == 0 for x in hostile_means[1:]),
                    f"{tag}: hostile component means, sign={sign}")
            require(not in_column_space(scalar_matrix, hostile),
                    f"{tag}: scalar hostile rejected, sign={sign}")
            require(in_column_space(repaired_matrix, hostile),
                    f"{tag}: componentwise hostile repaired, sign={sign}")

            # Canonical mixed-sign hostile.  The coefficient -(n+1) makes its
            # total sum strictly negative for every disconnected n-vertex
            # graph, so a nonzero global sum cannot hide the sign conflict.
            mixed_hostile = [
                Fraction(1) if i in comps[0] else Fraction(-(n + 1))
                for i in range(n)
            ]
            mixed_means = component_averages(mixed_hostile, comps)
            require(mixed_means[0] == 1 and
                    all(x == -(n + 1) for x in mixed_means[1:]),
                    f"{tag}: mixed-sign component means, sign={sign}")
            require(sum(mixed_hostile, Fraction(0)) < 0,
                    f"{tag}: mixed-sign total sum nonzero, sign={sign}")
            require(not in_column_space(scalar_matrix, mixed_hostile),
                    f"{tag}: mixed-sign scalar hostile rejected, sign={sign}")
            require(in_column_space(repaired_matrix, mixed_hostile),
                    f"{tag}: mixed-sign componentwise repair, sign={sign}")

    return c


def simple_graph_edges(n, mask):
    possible = list(combinations(range(n), 2))
    return [(u, v, Fraction(1)) for bit, (u, v) in enumerate(possible)
            if (mask >> bit) & 1]


def main():
    graph_count = 0
    connected_count = 0
    disconnected_count = 0
    max_components = 0

    for n in range(1, 6):
        edge_count = n * (n - 1) // 2
        for mask in range(1 << edge_count):
            tag = f"simple/n={n}/mask={mask}"
            c = audit_graph(n, simple_graph_edges(n, mask), tag)
            graph_count += 1
            connected_count += int(c == 1)
            disconnected_count += int(c > 1)
            max_components = max(max_components, c)

    weighted_suite = (
        (4, ((0, 1, Fraction(1, 2)), (2, 3, Fraction(2, 3)))),
        (5, ((0, 1, Fraction(3, 5)), (1, 2, Fraction(7, 11)))),
        (6, ((0, 1, Fraction(2, 7)), (1, 2, Fraction(5, 13)),
             (3, 4, Fraction(11, 17)))),
        (6, ((0, 1, Fraction(1, 3)), (1, 2, Fraction(4, 9)),
             (2, 0, Fraction(5, 8)), (3, 4, Fraction(7, 10)),
             (4, 5, Fraction(2, 11)))),
    )
    weighted_components = []
    for index, (n, raw_edges) in enumerate(weighted_suite):
        edges = [tuple(edge) for edge in raw_edges]
        weighted_components.append(audit_graph(n, edges, f"weighted/{index}"))

    # Integral sidecar.  Determinantal divisors audit every simple graph through
    # three vertices, so the triangle really is the smallest *simple* torsion
    # hostile.  (Parallel-edge multigraphs are deliberately outside this
    # universe and can already have torsion on two vertices.)
    torsion_hostiles = []
    triangle_snf = None
    for n in range(1, 4):
        edge_count = n * (n - 1) // 2
        for mask in range(1 << edge_count):
            L = laplacian(n, simple_graph_edges(n, mask), sign=1)
            snf = smith_invariants_by_minors(L)
            if any(value > 1 for value in snf):
                torsion_hostiles.append((n, mask, snf))
            if n == 3 and mask == 7:
                triangle_snf = snf

    require(torsion_hostiles == [(3, 7, (1, 3, 0))],
            "triangle is the smallest simple integral torsion hostile")
    require(triangle_snf == (1, 3, 0), "triangle Smith form")

    triangle = laplacian(3, simple_graph_edges(3, 7), sign=1)
    hostile = [Fraction(1), Fraction(-1), Fraction(0)]
    rational_solution = [x / 3 for x in hostile]
    require(sum(hostile, Fraction(0)) == 0, "triangle hostile real quotient vanishes")
    require(mat_vec(triangle, rational_solution) == hostile,
            "triangle hostile lies in rational image")
    require(mat_vec(triangle, hostile) == [3 * x for x in hostile],
            "three times triangle hostile lies in integral image")
    rows_mod_three = tuple(tuple(int(x) % 3 for x in row) for row in triangle)
    require(len(set(rows_mod_three)) == 1,
            "triangle image coordinates are congruent modulo three")
    hostile_mod_three = tuple(int(x) % 3 for x in hostile)
    require(len(set(hostile_mod_three)) == 3,
            "triangle hostile is not in the integral image")

    require(graph_count == 1099, "simple graph universe cardinality")
    require(connected_count + disconnected_count == graph_count,
            "connected/disconnected partition")
    require(max_components == 5, "empty five-vertex graph boundary")
    require(tuple(weighted_components) == (2, 3, 3, 2),
            "weighted-suite component counts")

    manifest = "|".join((
        "THM-3990 finite exact analogue",
        "simple labeled graphs n=1..5",
        f"graphs={graph_count}",
        f"connected={connected_count}",
        f"disconnected={disconnected_count}",
        "signs=+L,-L",
        f"triangle_snf={triangle_snf}",
        "triangle_hostile_order=3",
        "mixed_sign_hostile=+1 on C1, -(n+1) on all other components",
        f"weighted_component_counts={tuple(weighted_components)}",
        f"gates={GATES}",
        "normalized component averages; scalar diagonal; componentwise repair",
    ))

    print("THM-3990 componentwise Laplacian exact audit: PASS")
    print("arithmetic=fractions.Fraction (exact Q)")
    print("simple_universe=labeled loopless unweighted graphs, 1<=n<=5")
    print(f"simple_graphs={graph_count}")
    print(f"connected={connected_count}")
    print(f"disconnected={disconnected_count}")
    print("laplacian_signs=+L,-L")
    print("integral_simple_universe=all labeled simple graphs, 1<=n<=3")
    print(f"triangle_snf={triangle_snf}")
    print("triangle_hostile=(1,-1,0): rationally solvable, integral order=3")
    print(f"rational_weighted_controls={len(weighted_suite)}")
    print("rational_weighted_scope=positive-edge extension; not an integral claim")
    print(f"weighted_component_counts={tuple(weighted_components)}")
    print("positive_control=mu*1+Lz has equal normalized component averages")
    print("hostile_control=1_C1 has nonzero total sum but unequal component averages")
    print("hostile_result=single scalar rejected; componentwise constants repair")
    print("mixed_sign_hostile=+1 on C1, -(n+1) elsewhere; nonzero total sum")
    print("mixed_sign_result=single scalar rejected; componentwise constants repair")
    print(f"gates={GATES}")
    print(f"semantic_sha256={sha256(manifest.encode('utf-8')).hexdigest()}")


if __name__ == "__main__":
    main()
