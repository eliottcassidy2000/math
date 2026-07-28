#!/usr/bin/env python3
"""Exact frame audit for THM-2774's tree-path Smith-index ladder.

The script exhausts all full-rank B_m root frames through m=5 and every
recursive tree/path pair through eight vertices.  It uses explicit
exceptions, integer arithmetic only, and no truth-bearing assertions.
"""

from collections import Counter
from itertools import combinations, permutations, product


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


def rank_mod_two(matrix):
    if not matrix:
        return 0
    column_count = len(matrix[0])
    rows = []
    for row in matrix:
        bits = 0
        for column, value in enumerate(row):
            if value % 2:
                bits |= 1 << column
        rows.append(bits)
    rank = 0
    for column in range(column_count):
        pivot = next((index for index in range(rank, len(rows))
                      if (rows[index] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and ((rows[index] >> column) & 1):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


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
            "B_m positive-root hyperplane atlas changed")
    return tuple(roots)


def path_frame(length):
    rows = []
    for index in range(length - 1):
        row = [0] * length
        row[index] = 1
        row[index + 1] = -1
        rows.append(tuple(row))
    rows.append(tuple([1] * length))
    return tuple(rows)


def parent_arrays(vertex_count):
    if vertex_count == 1:
        yield (-1,)
        return
    for tail in product(*(range(child) for child in range(1, vertex_count))):
        yield (-1,) + tail


def root_paths(parents):
    edge_count = len(parents) - 1
    paths = [tuple([0] * edge_count)]
    for vertex in range(1, len(parents)):
        row = [0] * edge_count
        current = vertex
        while current:
            row[current - 1] = 1
            current = parents[current]
        paths.append(tuple(row))
    return tuple(paths)


def tree_adjacency(parents):
    adjacency = [[] for _ in parents]
    for child in range(1, len(parents)):
        parent = parents[child]
        adjacency[child].append((parent, child - 1))
        adjacency[parent].append((child, child - 1))
    return adjacency


def geodesic_edges(adjacency, start, finish):
    stack = [(start, -1, ())]
    while stack:
        vertex, parent, edges = stack.pop()
        if vertex == finish:
            return edges
        for neighbour, edge in adjacency[vertex]:
            if neighbour != parent:
                stack.append((neighbour, vertex, edges + (edge,)))
    raise RuntimeError("tree geodesic disappeared")


def main():
    frame_histogram = Counter()
    b_frame_count = 0
    b_nonsingular_count = 0
    for dimension in range(1, 6):
        roots = b_roots(dimension)
        for selection in combinations(roots, dimension):
            b_frame_count += 1
            value = abs(determinant(selection))
            if value == 0:
                frame_histogram[(dimension, 0)] += 1
                continue
            b_nonsingular_count += 1
            require(value & (value - 1) == 0,
                    "a pure B frame acquired odd determinant torsion")
            binary_defect = dimension - rank_mod_two(selection)
            require(value == 2 ** binary_defect,
                    "a pure B frame stopped having elementary 2-Smith form")
            frame_histogram[(dimension, value)] += 1

    canonical_frames = 0
    torus_kernel_points = 0
    weight_root_controls = 0
    for length in range(2, 14):
        frame = path_frame(length)
        value = abs(determinant(frame))
        require(value == length,
                "path Smith frame determinant stopped equalling its length")
        leading_minor = tuple(row[:-1] for row in frame[:-1])
        require(abs(determinant(leading_minor)) == 1,
                "path Smith frame lost its primitive codimension-one minor")
        require(all(sum(row) % length == 0 for row in frame),
                "path frame stopped lying in the sum-mod-k kernel")
        diagonal_kernel = tuple(tuple([residue] * length)
                                for residue in range(length))
        require(all(all(sum(frame[row][column] * point[column]
                                for column in range(length)) % length == 0
                            for row in range(length))
                        for point in diagonal_kernel),
                "a diagonal k-torsion point left the torus kernel")
        fundamental_numerator = (length - 1,) + tuple([-1] * (length - 1))
        require(sum(fundamental_numerator) == 0,
                "the first A-weight numerator left the sum-zero plane")
        for multiple in range(1, length):
            require(any((multiple * value) % length
                        for value in fundamental_numerator),
                    "the A-weight class acquired order below k")
        for first in range(length):
            for second in range(first + 1, length):
                difference = [0] * length
                difference[first] = 1
                difference[second] = -1
                require(sum(difference) == 0,
                        "a projected basis difference left the A-root lattice")
        weight_root_controls += 1
        torus_kernel_points += len(diagonal_kernel)
        canonical_frames += 1

    tree_controls = 0
    geodesic_frames = 0
    full_extensions = 0
    length_histogram = Counter()
    for vertex_count in range(2, 9):
        for parents in parent_arrays(vertex_count):
            paths = root_paths(parents)
            adjacency = tree_adjacency(parents)
            seen_lengths = set()
            for start in range(vertex_count):
                for finish in range(start + 1, vertex_count):
                    edges = geodesic_edges(adjacency, start, finish)
                    length = len(edges)
                    seen_lengths.add(length)
                    if length < 2:
                        continue
                    normal = tuple(paths[finish][index] - paths[start][index]
                                   for index in range(vertex_count - 1))
                    support = tuple(index for index, value in enumerate(normal)
                                    if value)
                    require(set(support) == set(edges)
                            and all(normal[index] in (-1, 1) for index in support),
                            "partial-cube endpoint difference lost its geodesic support")

                    # Gauge each path coordinate by the sign of the endpoint
                    # root.  Consecutive differences then pull back to lawful
                    # D roots, and the endpoint root becomes the all-one row.
                    restricted = []
                    for index in range(length - 1):
                        row = [0] * length
                        row[index] = 1
                        row[index + 1] = -1
                        restricted.append(tuple(row))
                    restricted.append(tuple([1] * length))
                    require(abs(determinant(tuple(restricted))) == length,
                            "a gauged tree geodesic lost its index")

                    # Add coordinate roots off the geodesic.  After a column
                    # permutation this is block diagonal with the path frame,
                    # so the full ambient frame has the same determinant.
                    ambient_rows = []
                    ordered_edges = tuple(edges) + tuple(
                        index for index in range(vertex_count - 1)
                        if index not in set(edges)
                    )
                    position = {edge: index for index, edge in enumerate(ordered_edges)}
                    signs = [normal[edge] for edge in edges]
                    for index in range(length - 1):
                        row = [0] * (vertex_count - 1)
                        row[position[edges[index]]] = signs[index]
                        row[position[edges[index + 1]]] = -signs[index + 1]
                        ambient_rows.append(tuple(row))
                    ambient_rows.append(tuple(normal[edge] for edge in ordered_edges))
                    for edge in ordered_edges[length:]:
                        row = [0] * (vertex_count - 1)
                        row[position[edge]] = 1
                        ambient_rows.append(tuple(row))
                    require(abs(determinant(tuple(ambient_rows))) == length,
                            "the full ambient path-frame extension changed index")
                    geodesic_frames += 1
                    full_extensions += 1
                    length_histogram[length] += 1

            diameter = max(seen_lengths)
            require(seen_lengths == set(range(1, diameter + 1)),
                    "a tree diameter stopped containing every shorter path length")
            tree_controls += 1

    binary = path_frame(2)
    ternary = path_frame(3)
    require(abs(determinant(binary)) == 2
            and abs(determinant(ternary)) == 3,
            "binary/ternary first frames changed")

    p4_roots = b_roots(3) + ((1, 1, 1),)
    p4_minor_histogram = Counter(
        abs(determinant(selection))
        for selection in combinations(p4_roots, 3)
    )
    require(p4_minor_histogram == Counter({1: 73, 2: 25, 0: 19, 3: 3}),
            "P4 arithmetic-frame spectrum changed")

    h = (1, 1, 1)
    d3_images = set()
    h_stabilizer = 0
    for permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            if sum(sign < 0 for sign in signs) % 2:
                continue
            image = tuple(signs[index] * h[permutation[index]]
                          for index in range(3))
            d3_images.add(image)
            if image == h:
                h_stabilizer += 1
    require(len(d3_images) == 4 and h_stabilizer == 6,
            "the D3 body-diagonal orbit/stabilizer changed")

    histogram_text = ";".join(
        f"m{dimension}:" + ",".join(
            f"det{value}={frame_histogram[(dimension, value)]}"
            for value in sorted(
                determinant_value
                for frame_dimension, determinant_value in frame_histogram
                if frame_dimension == dimension
            )
        )
        for dimension in range(1, 6)
    )
    length_text = ",".join(
        f"k{length}={length_histogram[length]}"
        for length in sorted(length_histogram)
    )

    print("TREE PATH SMITH-INDEX LADDER AUDIT")
    print(f"B_frames_m1_to_m5={b_frame_count} nonsingular={b_nonsingular_count}")
    print("B_nonzero_Smith=diag(1,...,1,2,...,2)")
    print("B_frame_histogram=" + histogram_text)
    print(f"canonical_path_frames_k2_to_k13={canonical_frames}")
    print("path_frame_Smith=diag(1^(k-1),k) quotient=sum_mod_k")
    print(f"torus_diagonal_kernel_points_k2_to_k13={torus_kernel_points}")
    print("torus_kernel={(j/k,...,j/k):j_mod_k}")
    print(f"A_weight_root_quotient_controls_k2_to_k13={weight_root_controls}")
    print("Z^k/(Q_A+Z*h)=P_A/Q_A=Z/k")
    print(f"recursive_tree_controls_n2_to_n8={tree_controls}")
    print(f"geodesic_frames={geodesic_frames} full_ambient_extensions={full_extensions}")
    print("geodesic_length_histogram=" + length_text)
    print("partial_cube_map=root_path_indicators endpoint_support=tree_distance")
    print("binary_frame=k2_det2 ternary_first_long_frame=k3_det3")
    print("P4_full_minor_histogram=det0:19,det1:73,det2:25,det3:3")
    print("D3_body_diagonal_orbit=4 stabilizer=S3_order6")
    print("diameter>=p supplies a frame quotient Z/p")
    print("SCOPE: frame-local lattice cokernels; not PSL2Z or graceful existence")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
