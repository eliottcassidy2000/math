#!/usr/bin/env python3
"""Exact graph edge-sum discriminant and graceful-gauge audit for THM-2761.

For a simple graph, every pair of edges is adjacent or disjoint.  Adjacent
edge-sum differences reduce to vertex differences, with exponent given by
vertex codegree; disjoint pairs retain the additive-energy factor.  On a
bipartite graph, a sign gauge converts edge sums into oriented label
differences.  Squaring those differences requires a second mirror factor to
detect collisions up to sign.  All checks are exact and use explicit
exceptions rather than Python assertions.
"""

from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def all_edges(n_vertices):
    return tuple(combinations(range(n_vertices), 2))


def graph_from_mask(n_vertices, mask):
    possible = all_edges(n_vertices)
    return tuple(
        edge for index, edge in enumerate(possible) if (mask >> index) & 1
    )


def discriminant_from_values(values):
    out = 1
    for i, j in combinations(range(len(values)), 2):
        out *= (values[i] - values[j]) ** 2
    return out


def edge_sums(edges, vertex_values):
    return tuple(vertex_values[u] + vertex_values[v] for u, v in edges)


def neighborhoods(n_vertices, edges):
    out = [set() for _ in range(n_vertices)]
    for u, v in edges:
        out[u].add(v)
        out[v].add(u)
    return tuple(frozenset(row) for row in out)


def codegree_table(n_vertices, edges):
    neighbor_sets = neighborhoods(n_vertices, edges)
    return {
        (u, v): len(neighbor_sets[u] & neighbor_sets[v])
        for u, v in combinations(range(n_vertices), 2)
    }


def factor_inventory(n_vertices, edges):
    adjacent_counts = {
        pair: 0 for pair in combinations(range(n_vertices), 2)
    }
    disjoint_pairs = []
    for first_index, second_index in combinations(range(len(edges)), 2):
        first = set(edges[first_index])
        second = set(edges[second_index])
        common = first & second
        if common:
            require(len(common) == 1,
                    "distinct simple edges acquired two common endpoints")
            outer = tuple(sorted((
                next(iter(first - common)), next(iter(second - common))
            )))
            adjacent_counts[outer] += 1
        else:
            disjoint_pairs.append((first_index, second_index))
    return adjacent_counts, tuple(disjoint_pairs)


def factorized_discriminant(n_vertices, edges, vertex_values):
    codegrees = codegree_table(n_vertices, edges)
    out = 1
    for (u, v), codegree in codegrees.items():
        out *= (vertex_values[u] - vertex_values[v]) ** (2 * codegree)
    sums = edge_sums(edges, vertex_values)
    for first_index, second_index in combinations(range(len(edges)), 2):
        if set(edges[first_index]).isdisjoint(edges[second_index]):
            out *= (sums[first_index] - sums[second_index]) ** 2
    return out


def bipartition(n_vertices, edges):
    neighbor_sets = neighborhoods(n_vertices, edges)
    colors = [None] * n_vertices
    for start in range(n_vertices):
        if colors[start] is not None:
            continue
        colors[start] = 0
        stack = [start]
        while stack:
            vertex = stack.pop()
            for neighbor in neighbor_sets[vertex]:
                proposed = 1 - colors[vertex]
                if colors[neighbor] is None:
                    colors[neighbor] = proposed
                    stack.append(neighbor)
                elif colors[neighbor] != proposed:
                    return None
    return tuple(colors)


def connected(n_vertices, edges):
    if n_vertices == 0:
        return True
    neighbor_sets = neighborhoods(n_vertices, edges)
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in neighbor_sets[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == n_vertices


def oriented_differences(edges, labels, colors):
    out = []
    for u, v in edges:
        if colors[u] == 0 and colors[v] == 1:
            out.append(labels[u] - labels[v])
        elif colors[u] == 1 and colors[v] == 0:
            out.append(labels[v] - labels[u])
        else:
            raise RuntimeError("an edge failed to cross the bipartition")
    return tuple(out)


def mirror_factor(differences):
    out = 1
    for i, j in combinations(range(len(differences)), 2):
        out *= (differences[i] + differences[j]) ** 2
    return out


def main():
    graph_count = 0
    edge_pair_instances = 0
    adjacent_instances = 0
    disjoint_instances = 0
    value_patterns = (
        lambda n: tuple(index ** 3 + 2 * index + 1 for index in range(n)),
        lambda n: tuple(((-1) ** index) * (index + 2) for index in range(n)),
    )

    for n_vertices in range(2, 7):
        possible_count = len(all_edges(n_vertices))
        for mask in range(1 << possible_count):
            graph_count += 1
            edges = graph_from_mask(n_vertices, mask)
            adjacent_counts, disjoint_pairs = factor_inventory(
                n_vertices, edges
            )
            codegrees = codegree_table(n_vertices, edges)
            require(adjacent_counts == codegrees,
                    "adjacent edge pairs stopped being counted by codegrees")
            pair_count = len(edges) * (len(edges) - 1) // 2
            require(sum(adjacent_counts.values()) + len(disjoint_pairs)
                    == pair_count,
                    "the adjacent/disjoint edge-pair partition changed")
            edge_pair_instances += pair_count
            adjacent_instances += sum(adjacent_counts.values())
            disjoint_instances += len(disjoint_pairs)

            for make_values in value_patterns:
                values = make_values(n_vertices)
                lhs = discriminant_from_values(edge_sums(edges, values))
                rhs = factorized_discriminant(n_vertices, edges, values)
                require(lhs == rhs,
                        "the graph edge-sum discriminant factorization failed")

    require(graph_count == 33866
            and edge_pair_instances == 871926
            and adjacent_instances == 499398
            and disjoint_instances == 372528,
            "the all-graphs census changed")

    # K4 recovers the quartic pair-sum theorem: every codegree is two and the
    # three disjoint edge pairs produce the opposite-matching factor.
    k4_edges = all_edges(4)
    k4_values = (0, 1, 3, 7)
    k4_codegrees = codegree_table(4, k4_edges)
    require(set(k4_codegrees.values()) == {2},
            "K4 codegree stopped being two")
    k4_discriminant = discriminant_from_values(edge_sums(k4_edges, k4_values))
    vertex_discriminant = discriminant_from_values(k4_values)
    opposite_differences = tuple(
        k4_values[a] + k4_values[b] - k4_values[c] - k4_values[d]
        for (a, b), (c, d) in (
            ((0, 1), (2, 3)),
            ((0, 2), (1, 3)),
            ((0, 3), (1, 2)),
        )
    )
    opposite_square = 1
    for value in opposite_differences:
        opposite_square *= value * value
    require(k4_discriminant
            == vertex_discriminant ** 2 * opposite_square,
            "the K4 specialization stopped recovering THM-2758")

    # Exhaustive graceful criterion on connected bipartite graphs through four
    # vertices, over every injection into {0,...,m}.
    graceful_labelings = 0
    injection_controls = 0
    bipartite_graphs = 0
    for n_vertices in range(2, 5):
        possible_count = len(all_edges(n_vertices))
        for mask in range(1 << possible_count):
            edges = graph_from_mask(n_vertices, mask)
            if not edges or not connected(n_vertices, edges):
                continue
            colors = bipartition(n_vertices, edges)
            if colors is None:
                continue
            edge_count = len(edges)
            if n_vertices > edge_count + 1:
                continue
            bipartite_graphs += 1
            for labels in permutations(range(edge_count + 1), n_vertices):
                injection_controls += 1
                differences = oriented_differences(edges, labels, colors)
                signed_discriminant = discriminant_from_values(differences)
                mirror = mirror_factor(differences)
                squared_discriminant = discriminant_from_values(
                    tuple(value * value for value in differences)
                )
                require(squared_discriminant
                        == signed_discriminant * mirror,
                        "the absolute-value mirror factorization failed")
                graceful = sorted(abs(value) for value in differences) == list(
                    range(1, edge_count + 1)
                )
                require(graceful == (squared_discriminant != 0),
                        "the bounded-injection graceful iff changed")
                if graceful:
                    graceful_labelings += 1

    # Sharp hostile: oriented differences are distinct but collide in absolute
    # value.  The second factor, not the first discriminant, detects it.
    path_edges = ((0, 1), (1, 2))
    path_colors = (0, 1, 0)
    hostile_labels = (0, 1, 2)
    hostile_differences = oriented_differences(
        path_edges, hostile_labels, path_colors
    )
    hostile_signed = discriminant_from_values(hostile_differences)
    hostile_mirror = mirror_factor(hostile_differences)
    hostile_squared = discriminant_from_values(
        tuple(value * value for value in hostile_differences)
    )
    require(hostile_differences == (-1, 1)
            and hostile_signed == 4
            and hostile_mirror == 0
            and hostile_squared == 0,
            "the oriented-versus-absolute hostile changed")

    positive_labels = (0, 2, 1)
    positive_differences = oriented_differences(
        path_edges, positive_labels, path_colors
    )
    require(sorted(abs(value) for value in positive_differences) == [1, 2]
            and discriminant_from_values(
                tuple(value * value for value in positive_differences)
            ) != 0,
            "the graceful P3 positive control changed")

    print("GRAPH EDGE-SUM DISCRIMINANT / GRACEFUL SIGN-GAUGE AUDIT")
    print("factorization=codegree_vertex_factors*disjoint_edge_energy_factors")
    print(f"graph_census_n2_to_n6={graph_count}")
    print(
        f"edge_pair_instances={edge_pair_instances} "
        f"adjacent={adjacent_instances} disjoint={disjoint_instances}"
    )
    print("numeric_factorization_patterns=2_per_graph all_exact")
    print("K4_specialization=codegree2+3_opposite_pairs=THM2758")
    print("bipartite_gauge=x_A=label x_B=-label edge_sum=oriented_difference")
    print("disc(edge_difference_squares)=disc(oriented)*mirror_plus_factor")
    print(
        f"graceful_census_graphs={bipartite_graphs} "
        f"injections={injection_controls} graceful={graceful_labelings}"
    )
    print("hostile=P3 labels(0,1,2) oriented=(-1,1) disc=4 mirror=0")
    print("positive=P3 labels(0,2,1) absolute_differences=(2,1)")
    print("SCOPE: universal graph identity and graceful criterion; no existence theorem")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
