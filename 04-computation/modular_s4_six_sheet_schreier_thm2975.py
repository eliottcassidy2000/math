"""Exact companion for THM-2975.

Compare the marked C2*C3 -> S4 actions on the six edges of K4 and the six
oriented Hamilton cycles.  Audit their discriminant-subgroup restriction,
Bass--Serre incidence shadows and orbifold signatures, local V4 frame,
chiral A4 octahedron, partial-cube/graceful/tournament boundaries, and
finite-kernel hostiles.
"""

from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(g, h):
    return tuple(g[h[i]] for i in range(len(g)))


def inverse(g):
    result = [None] * len(g)
    for i, j in enumerate(g):
        result[j] = i
    return tuple(result)


def parity(g):
    return sum(
        g[i] > g[j]
        for i in range(len(g))
        for j in range(i + 1, len(g))
    ) % 2


def rotate_to_zero(cycle):
    j = cycle.index(0)
    return cycle[j:] + cycle[:j]


def relabel_edge(g, edge):
    return tuple(sorted((g[edge[0]], g[edge[1]])))


def relabel_cycle(g, cycle):
    return rotate_to_zero(tuple(g[i] for i in cycle))


def induced_permutation(objects, relabel, g):
    lookup = {x: i for i, x in enumerate(objects)}
    return tuple(lookup[relabel(g, x)] for x in objects)


def cycles_of(p):
    result = []
    seen = set()
    for i in range(len(p)):
        if i in seen:
            continue
        cycle = []
        j = i
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = p[j]
        result.append(tuple(cycle))
    return tuple(sorted(result, key=lambda cycle: (len(cycle), cycle)))


def permutation_order(p):
    identity = tuple(range(len(p)))
    power = identity
    for order in range(1, 100):
        power = compose(p, power)
        if power == identity:
            return order
    raise RuntimeError("permutation order exceeded exact search bound")


def undirected_generator_edges(p, include_loops=False):
    edges = set()
    for i, j in enumerate(p):
        if i == j and not include_loops:
            continue
        edges.add(tuple(sorted((i, j))))
    return edges


def simple_connected(n, edges):
    if n == 0:
        return True
    seen = {0}
    frontier = [0]
    while frontier:
        i = frontier.pop()
        for edge in edges:
            if i not in edge:
                continue
            j = edge[1] if edge[0] == i else edge[0]
            if j not in seen:
                seen.add(j)
                frontier.append(j)
    return len(seen) == n


def bipartite(n, edges):
    color = {}
    for root in range(n):
        if root in color:
            continue
        color[root] = 0
        frontier = [root]
        while frontier:
            i = frontier.pop()
            for edge in edges:
                if i not in edge:
                    continue
                j = edge[1] if edge[0] == i else edge[0]
                if j in color:
                    if color[j] == color[i]:
                        return False
                else:
                    color[j] = 1 - color[i]
                    frontier.append(j)
    return True


def degree_sequence(n, edges):
    return tuple(sorted(sum(i in edge for edge in edges) for i in range(n)))


def graceful_labeling(n, edges):
    """Return one standard graceful labelling, or None."""
    q = len(edges)
    for labels in permutations(range(q + 1), n):
        differences = {abs(labels[i] - labels[j]) for i, j in edges}
        if differences == set(range(1, q + 1)):
            return labels
    return None


def permutation_group(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    frontier = [identity]
    while frontier:
        g = frontier.pop()
        for h in generators:
            for product in (compose(g, h), compose(h, g)):
                if product not in seen:
                    seen.add(product)
                    frontier.append(product)
    return seen


def edge_orbits(group, n):
    pairs = set(combinations(range(n), 2))
    orbits = []
    while pairs:
        seed = min(pairs)
        orbit = {
            tuple(sorted((g[seed[0]], g[seed[1]])))
            for g in group
        }
        pairs -= orbit
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits, key=lambda orbit: (len(orbit), orbit)))


def directed_pair_orbits(group, n):
    pairs = {(i, j) for i in range(n) for j in range(n) if i != j}
    orbits = []
    while pairs:
        seed = min(pairs)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        pairs -= orbit
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits, key=lambda orbit: (len(orbit), orbit)))


def has_invariant_tournament(group, n):
    for i, j in combinations(range(n), 2):
        if any(g[i] == j and g[j] == i for g in group):
            return False
    return True


def orbit_partition(p):
    return cycles_of(p)


def incidence_matrix(left_orbits, right_orbits):
    return tuple(
        tuple(len(set(left) & set(right)) for right in right_orbits)
        for left in left_orbits
    )


def simple_incidence_edges(left_orbits, right_orbits):
    matrix = incidence_matrix(left_orbits, right_orbits)
    return {
        (i, len(left_orbits) + j)
        for i, row in enumerate(matrix)
        for j, value in enumerate(row)
        if value
    }


S4 = list(permutations(range(4)))
A4 = [g for g in S4 if parity(g) == 0]

# THM-2968's marked modular quotient.
a = (0, 1, 3, 2)  # (2 3), order two
b = (2, 0, 1, 3)  # (0 2 1), order three
c = compose(compose(a, b), a)

require(compose(a, a) == tuple(range(4)), "a is not an involution")
require(
    compose(compose(b, b), b) == tuple(range(4)),
    "b is not order three",
)
require(
    compose(compose(c, c), c) == tuple(range(4)),
    "aba is not order three",
)
require(permutation_group((a, b)) == set(S4), "a,b do not generate S4")
require(permutation_group((b, c)) == set(A4), "b,aba do not generate A4")
require(
    permutation_order(compose(a, b)) == 4,
    "the parabolic word ab did not collapse to order four",
)
require(
    permutation_order(compose(b, c)) == 2,
    "the free C3*C3 word bc did not collapse to order two",
)

edges = ((0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2))
orientation_pairs = (
    ((0, 2, 1, 3), (0, 3, 1, 2)),
    ((0, 3, 2, 1), (0, 1, 2, 3)),
    ((0, 1, 3, 2), (0, 2, 3, 1)),
)
orientations = tuple(obj for pair in orientation_pairs for obj in pair)

actions = {}
for name, objects, relabel in (
    ("edge", edges, relabel_edge),
    ("orientation", orientations, relabel_cycle),
):
    pa = induced_permutation(objects, relabel, a)
    pb = induced_permutation(objects, relabel, b)
    pc = induced_permutation(objects, relabel, c)
    require(
        len(permutation_group((pa, pb))) == 24,
        f"{name}: S4 image changed",
    )
    require(
        len(permutation_group((pb, pc))) == 12,
        f"{name}: A4 image changed",
    )
    actions[name] = (pa, pb, pc)

edge_a, edge_b, edge_c = actions["edge"]
or_a, or_b, or_c = actions["orientation"]

global_pair_flip = (1, 0, 3, 2, 5, 4)
for g in S4:
    edge_g = induced_permutation(edges, relabel_edge, g)
    orientation_g = induced_permutation(orientations, relabel_cycle, g)
    expected_orientation = (
        edge_g
        if parity(g) == 0
        else compose(global_pair_flip, edge_g)
    )
    require(
        orientation_g == expected_orientation,
        f"pointwise central twist failed at {g}",
    )

require(
    edge_b == or_b and edge_c == or_c,
    "A4 restrictions are not identical",
)
common_a4_group = permutation_group((edge_b, edge_c))

# Marked descent of the common A4-set.  Its equivariant automorphism group is
# exactly the three-pair flip C2.  Relative to the chosen odd generator a,
# there are exactly two semilinear involutions, namely the edge and
# orientation descents.
identity6 = tuple(range(6))
sheet_permutations = tuple(permutations(range(6)))
equivariant_automorphisms = {
    u
    for u in sheet_permutations
    if all(compose(u, g) == compose(g, u) for g in common_a4_group)
}
require(
    equivariant_automorphisms == {identity6, global_pair_flip},
    "Aut_A4(X) is not the expected C2",
)
require(
    compose(global_pair_flip, edge_a) == compose(edge_a, global_pair_flip),
    "the quotient C2 acts nontrivially on the descent coefficient C2",
)
semilinear_descents = {
    u
    for u in sheet_permutations
    if compose(u, u) == identity6
    and compose(compose(u, edge_b), inverse(u)) == edge_c
    and compose(compose(u, edge_c), inverse(u)) == edge_b
}
require(
    semilinear_descents == {edge_a, or_a},
    "the marked A4-set does not have exactly the two expected descents",
)

for name, (pa, pb, pc) in actions.items():
    original_edges = (
        undirected_generator_edges(pa) | undirected_generator_edges(pb)
    )
    plus_edges = (
        undirected_generator_edges(pb) | undirected_generator_edges(pc)
    )
    require(
        simple_connected(6, original_edges),
        f"{name}: original graph disconnected",
    )
    require(
        not bipartite(6, original_edges),
        f"{name}: original graph unexpectedly bipartite",
    )
    expected_plus = {
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 5),
        (3, 4),
    }
    require(
        plus_edges == expected_plus,
        f"{name}: common A4 simple graph changed",
    )
    require(
        not bipartite(6, plus_edges),
        f"{name}: A4 graph unexpectedly bipartite",
    )

# The normal V4 has three two-point orbits, not a free six-point torsor.
v4_sheet = {
    (0, 1, 2, 3),
    (1, 0, 3, 2),
    (2, 3, 0, 1),
    (3, 2, 1, 0),
}
for name, objects, relabel in (
    ("edge", edges, relabel_edge),
    ("orientation", orientations, relabel_cycle),
):
    v4_action = {
        induced_permutation(objects, relabel, g)
        for g in v4_sheet
    }
    point_orbits = []
    unseen = set(range(6))
    while unseen:
        i = min(unseen)
        orbit = {g[i] for g in v4_action}
        unseen -= orbit
        point_orbits.append(tuple(sorted(orbit)))
    require(
        tuple(sorted(point_orbits)) == ((0, 1), (2, 3), (4, 5)),
        f"{name}: V4 fibre blocks changed",
    )
    require(
        all(sum(g[i] == i for g in v4_action) == 2 for i in range(6)),
        f"{name}: V4 point stabilizer is not C2",
    )

edge_s4_group = permutation_group((edge_a, edge_b))
or_s4_group = permutation_group((or_a, or_b))
require(
    not has_invariant_tournament(edge_s4_group, 6),
    "edge S4 tournament unexpectedly exists",
)
require(
    not has_invariant_tournament(or_s4_group, 6),
    "orientation S4 tournament unexpectedly exists",
)
require(
    not has_invariant_tournament(common_a4_group, 6),
    "common A4 tournament unexpectedly exists",
)

s4_directed_orbits = directed_pair_orbits(edge_s4_group, 6)
a4_directed_orbits = directed_pair_orbits(common_a4_group, 6)
require(
    tuple(len(orbit) for orbit in s4_directed_orbits) == (6, 24),
    "S4 directed orbital sizes changed",
)
require(
    tuple(len(orbit) for orbit in a4_directed_orbits) == (6, 12, 12),
    "A4 directed orbital sizes changed",
)
require(
    all(
        (orbit[0][1], orbit[0][0]) in orbit
        for orbit in s4_directed_orbits
    ),
    "an S4 orbital is not self-paired",
)
require(
    (a4_directed_orbits[0][0][1], a4_directed_orbits[0][0][0])
    in a4_directed_orbits[0],
    "A4 opposite orbital is not self-paired",
)
require(
    all(
        (orbit[0][1], orbit[0][0]) not in orbit
        for orbit in a4_directed_orbits[1:]
    ),
    "an A4 incident orbital is self-paired",
)
require(
    {(j, i) for i, j in a4_directed_orbits[1]}
    == set(a4_directed_orbits[2]),
    "A4 incident orbitals are not inverse",
)
a4_chiral_arcs = set(a4_directed_orbits[1])
require(
    all(sum(i == v for i, _ in a4_chiral_arcs) == 2 for v in range(6)),
    "A4 chiral octahedron outdegree changed",
)
require(
    all(sum(j == v for _, j in a4_chiral_arcs) == 2 for v in range(6)),
    "A4 chiral octahedron indegree changed",
)
directed_triangles = {
    tuple(sorted((i, j, k)))
    for i, j, k in permutations(range(6), 3)
    if (
        (i, j) in a4_chiral_arcs
        and (j, k) in a4_chiral_arcs
        and (k, i) in a4_chiral_arcs
    )
}
require(
    len(directed_triangles) == 8,
    "A4 chiral octahedron does not orient every triangular face cyclically",
)

edge_graph = (
    undirected_generator_edges(edge_a) | undirected_generator_edges(edge_b)
)
or_graph = (
    undirected_generator_edges(or_a) | undirected_generator_edges(or_b)
)
common_graph = (
    undirected_generator_edges(edge_b)
    | undirected_generator_edges(edge_c)
)

edge_incidence = incidence_matrix(
    orbit_partition(edge_a), orbit_partition(edge_b)
)
or_incidence = incidence_matrix(
    orbit_partition(or_a), orbit_partition(or_b)
)
plus_incidence = incidence_matrix(
    orbit_partition(edge_b), orbit_partition(edge_c)
)
require(
    edge_incidence == ((1, 0), (0, 1), (1, 1), (1, 1)),
    "edge Bass--Serre quotient changed",
)
require(
    or_incidence == ((1, 1), (2, 0), (0, 2)),
    "orientation Bass--Serre quotient changed",
)
require(
    plus_incidence == ((1, 2), (2, 1)),
    "A4 Bass--Serre quotient changed",
)
edge_incidence_simple = simple_incidence_edges(
    orbit_partition(edge_a), orbit_partition(edge_b)
)
or_incidence_simple = simple_incidence_edges(
    orbit_partition(or_a), orbit_partition(or_b)
)
plus_incidence_simple = simple_incidence_edges(
    orbit_partition(edge_b), orbit_partition(edge_c)
)
require(
    bipartite(6, edge_incidence_simple),
    "edge incidence shadow is not bipartite",
)
require(
    bipartite(5, or_incidence_simple),
    "orientation incidence shadow is not bipartite",
)
require(
    bipartite(4, plus_incidence_simple),
    "A4 incidence shadow is not bipartite",
)
require(
    degree_sequence(6, edge_incidence_simple) == (1, 1, 2, 2, 3, 3),
    "edge incidence shadow changed",
)
require(
    degree_sequence(5, or_incidence_simple) == (1, 1, 2, 2, 2),
    "orientation incidence shadow is not P5",
)
require(
    degree_sequence(4, plus_incidence_simple) == (2, 2, 2, 2),
    "A4 incidence shadow is not C4",
)

# The weighted common C4 has abstract regular V4 symmetry if the two C3
# factor colors may be exchanged.
factor_subgroups = (
    frozenset(permutation_group((b,))),
    frozenset(permutation_group((c,))),
)
factor_pair_normalizer = {
    g
    for g in S4
    if {
        frozenset(
            compose(compose(g, h), inverse(g))
            for h in subgroup
        )
        for subgroup in factor_subgroups
    }
    == set(factor_subgroups)
}
require(
    len(factor_pair_normalizer) == 4,
    "factor-pair normalizer is not V4",
)
require(
    all(
        compose(g, h) == compose(h, g)
        for g in factor_pair_normalizer
        for h in factor_pair_normalizer
    ),
    "factor-pair normalizer is not abelian",
)

common_point_stabilizer = {
    g for g in A4 if relabel_edge(g, edges[0]) == edges[0]
}
require(
    len(common_point_stabilizer) == 2,
    "common A4 point stabilizer is not C2",
)
require(
    common_point_stabilizer <= factor_pair_normalizer,
    "common C2 left the factor-pair normalizer",
)

# Local D8/C2=V4 extension frame around one common A4 sheet.
d = (1, 0, 3, 2)  # (0 1)(2 3)
require(
    common_point_stabilizer == {tuple(range(4)), d},
    "common stabilizer generator changed",
)
normal_v4 = v4_sheet
edge_stabilizer = {
    g for g in S4 if relabel_edge(g, edges[0]) == edges[0]
}
orientation_stabilizer = {
    g
    for g in S4
    if relabel_cycle(g, orientations[0]) == orientations[0]
}
centralizer_d = {
    g for g in S4 if compose(g, d) == compose(d, g)
}
require(
    len(normal_v4)
    == len(edge_stabilizer)
    == len(orientation_stabilizer)
    == 4,
    "order-four subgroup census changed",
)
require(len(centralizer_d) == 8, "common centralizer is not D8")
require(
    edge_stabilizer == factor_pair_normalizer,
    "factor-pair normalizer is not the selected edge V4",
)
require(
    edge_stabilizer & set(A4)
    == orientation_stabilizer & set(A4)
    == common_point_stabilizer,
    "A4 stabilizers are not the common C2",
)
require(
    normal_v4 <= centralizer_d
    and edge_stabilizer <= centralizer_d
    and orientation_stabilizer <= centralizer_d,
    "an order-four extension left D8",
)
require(
    normal_v4 & edge_stabilizer
    == normal_v4 & orientation_stabilizer
    == edge_stabilizer & orientation_stabilizer
    == common_point_stabilizer,
    "the three D8 intermediate lines do not meet in the common C2",
)
require(
    any(permutation_order(g) == 4 for g in orientation_stabilizer),
    "orientation stabilizer is not C4",
)
require(
    all(
        permutation_order(g) <= 2
        for g in normal_v4 | edge_stabilizer
    ),
    "a V4 extension acquired order four",
)

# Reconstruct the full S4 action from each semilinear descent and verify the
# actual base-sheet stabilizer.  Cycle types alone would not identify the two
# H^1 classes.
a4_sheet_action = {
    g: induced_permutation(edges, relabel_edge, g)
    for g in A4
}


def reconstructed_s4_action(descent, g):
    if parity(g) == 0:
        return a4_sheet_action[g]
    even_part = compose(g, a)  # g=(g a)a and g a lies in A4
    require(even_part in a4_sheet_action, "odd element did not split over A4")
    return compose(a4_sheet_action[even_part], descent)


for descent, objects, relabel, expected_stabilizer in (
    (edge_a, edges, relabel_edge, edge_stabilizer),
    (or_a, orientations, relabel_cycle, orientation_stabilizer),
):
    require(
        all(
            reconstructed_s4_action(descent, g)
            == induced_permutation(objects, relabel, g)
            for g in S4
        ),
        "semilinear descent failed to reconstruct the marked S4 action",
    )
    reconstructed_stabilizer = {
        g for g in S4 if reconstructed_s4_action(descent, g)[0] == 0
    }
    require(
        reconstructed_stabilizer == expected_stabilizer,
        "descent class did not recover its claimed stabilizer",
    )

factor_nodes = tuple(
    frozenset(orbit)
    for orbit in orbit_partition(edge_b) + orbit_partition(edge_c)
)
factor_node_lookup = {node: i for i, node in enumerate(factor_nodes)}
require(
    len(factor_node_lookup) == 4,
    "factor orbit nodes are not four distinct sets",
)


def induced_factor_node_action(objects, relabel, g):
    pg = induced_permutation(objects, relabel, g)
    result = []
    for node in factor_nodes:
        image = frozenset(pg[i] for i in node)
        require(
            image in factor_node_lookup,
            "normalizer did not preserve factor nodes",
        )
        result.append(factor_node_lookup[image])
    return tuple(result)


edge_factor_v4 = {
    induced_factor_node_action(edges, relabel_edge, g)
    for g in factor_pair_normalizer
}
or_factor_v4 = {
    induced_factor_node_action(orientations, relabel_cycle, g)
    for g in factor_pair_normalizer
}

four_vertices = range(4)
weighted_edges = {
    tuple(sorted((i, 2 + j))): plus_incidence[i][j]
    for i in range(2)
    for j in range(2)
}
weighted_aut = set()
for p4 in permutations(four_vertices):
    transported = {
        tuple(sorted((p4[i], p4[j]))): weight
        for (i, j), weight in weighted_edges.items()
    }
    if transported == weighted_edges:
        weighted_aut.add(p4)
require(
    len(weighted_aut) == 4,
    "weighted common C4 automorphism group is not V4",
)
require(
    all(len({g[i] for g in weighted_aut}) == 4 for i in four_vertices),
    "abstract V4 does not act regularly",
)
factor_action_table = tuple(
    sorted(
        (
            g,
            induced_factor_node_action(edges, relabel_edge, g),
            induced_factor_node_action(
                orientations, relabel_cycle, g
            ),
        )
        for g in factor_pair_normalizer
    )
)
require(
    len(edge_factor_v4) == 2,
    "edge factor-frame action did not have C2 image",
)
require(
    or_factor_v4 == weighted_aut,
    "orientation factor-frame action did not realize abstract V4",
)
require(
    all(len({g[i] for g in or_factor_v4}) == 4 for i in four_vertices),
    "orientation factor V4 is not regular",
)

edge_graceful = graceful_labeling(6, tuple(sorted(edge_graph)))
or_graceful = graceful_labeling(6, tuple(sorted(or_graph)))
common_graceful = graceful_labeling(6, tuple(sorted(common_graph)))

# Modular-orbifold signatures of the two index-six sheet stabilizers.  The
# standard parabolic is T=a^{-1}b=ab because a is an involution.  A fixed
# point of a (respectively b) is an elliptic point of order two
# (respectively three), and the T-cycle lengths are the cusp widths.
edge_parabolic = compose(edge_a, edge_b)
or_parabolic = compose(or_a, or_b)
edge_cusp_widths = tuple(sorted(len(cycle) for cycle in cycles_of(edge_parabolic)))
or_cusp_widths = tuple(sorted(len(cycle) for cycle in cycles_of(or_parabolic)))
edge_e2 = sum(edge_a[i] == i for i in range(6))
or_e2 = sum(or_a[i] == i for i in range(6))
edge_e3 = sum(edge_b[i] == i for i in range(6))
or_e3 = sum(or_b[i] == i for i in range(6))


def twelve_times_genus(index, e2, e3, cusp_count):
    """Return 12g from the finite-index PSL2(Z) signature formula."""

    return 12 + index - 3 * e2 - 4 * e3 - 6 * cusp_count


require(edge_cusp_widths == (2, 4), "edge cusp widths are not (2,4)")
require(or_cusp_widths == (1, 1, 4), "orientation cusp widths are not (1,1,4)")
require((edge_e2, edge_e3) == (2, 0), "edge elliptic census changed")
require((or_e2, or_e3) == (0, 0), "orientation elliptic census changed")
require(
    twelve_times_genus(6, edge_e2, edge_e3, len(edge_cusp_widths)) == 0,
    "edge stabilizer genus is not zero",
)
require(
    twelve_times_genus(6, or_e2, or_e3, len(or_cusp_widths)) == 0,
    "orientation stabilizer genus is not zero",
)

# The weighted factor-orbit incidence graphs are the Bass--Serre
# graphs-of-groups.  Their cycle ranks retain the free factors erased by the
# simple P5/C4 shadows.
edge_incidence_beta = 6 - len(edge_incidence) - len(edge_incidence[0]) + 1
or_incidence_beta = 6 - len(or_incidence) - len(or_incidence[0]) + 1
plus_incidence_beta = 6 - len(plus_incidence) - len(plus_incidence[0]) + 1
require(
    (edge_incidence_beta, or_incidence_beta, plus_incidence_beta) == (1, 2, 3),
    "Bass--Serre incidence cycle ranks changed",
)

print("modular S4 six-sheet Schreier probe")
print(f"sheet_generators=a:{a};b:{b};aba:{c}")
print(f"edge_cycles=a:{cycles_of(edge_a)};b:{cycles_of(edge_b)}")
print(f"orientation_cycles=a:{cycles_of(or_a)};b:{cycles_of(or_b)}")
print(
    "edge_loopless="
    f"edges:{len(edge_graph)};degrees:{degree_sequence(6, edge_graph)};"
    "partial_cube:NO"
)
print(
    "orientation_loopless="
    f"edges:{len(or_graph)};degrees:{degree_sequence(6, or_graph)};"
    "partial_cube:NO"
)
print(
    "A4_generators=b,aba;common_graph=K2,4_plus_2K2;"
    f"edges:{len(common_graph)};degrees:{degree_sequence(6, common_graph)};"
    "partial_cube:NO"
)
print("V4_orbits=((0,1),(2,3),(4,5));point_stabilizer_order=2;torsor:NO")
print(
    "local_extension_frame=C_S4((01)(23))/C2=D8/C2=V4;"
    "lines=normal_V4,edge_V4,orientation_C4"
)
print(
    "binary_descent_h1=Aut_A4(X):C2;semilinear_involutions:2;"
    "classes:edge_V4,orientation_C4;difference:three_pair_flip"
)
print(
    "Bass_Serre_incidence="
    f"edge:{edge_incidence};orientation:{or_incidence};A4:{plus_incidence}"
)
print(
    "incidence_simple_shadows=edge:C4_plus_two_opposite_leaves;"
    "orientation:P5;A4:C4;partial_cubes:YES"
)
print(
    "modular_signatures="
    f"edge:e2={edge_e2},e3={edge_e3},cusps={edge_cusp_widths},g=0;"
    f"orientation:e2={or_e2},e3={or_e3},cusps={or_cusp_widths},g=0"
)
print(
    "Bass_Serre_graph_of_groups="
    f"edge:C2*C2*Z,beta={edge_incidence_beta};"
    f"orientation:F2,beta={or_incidence_beta};"
    f"A4_stabilizer:F3,beta={plus_incidence_beta}"
)
print(
    "A4_weighted_C4_abstract_aut="
    f"V4_regular:{len(weighted_aut)};"
    f"physical_factor_pair_normalizer:{len(factor_pair_normalizer)};"
    f"common_stabilizer:{len(common_point_stabilizer)}"
)
print(
    f"factor_frame_actions={factor_action_table};"
    "edge_image:C2;orientation_image:V4_regular"
)
print(
    "Farey_hostile=(ab)^4 collapses although ab is parabolic infinite order;"
    "cover_hostile=(b*aba)^2 collapses in C3*C3"
)
print(
    f"graceful=edge:{edge_graceful};orientation:{or_graceful};"
    f"A4:{common_graceful}"
)
print("invariant_tournament=S4_edge:NO;S4_orientation:NO;A4_common:NO")
print(
    "pair_orbits_S4_edge="
    f"{tuple(len(orbit) for orbit in edge_orbits(edge_s4_group, 6))}"
)
print(
    "pair_orbits_A4_common="
    f"{tuple(len(orbit) for orbit in edge_orbits(common_a4_group, 6))}"
)
print(
    "directed_orbitals="
    f"S4:{tuple(len(orbit) for orbit in s4_directed_orbits)};"
    f"A4:{tuple(len(orbit) for orbit in a4_directed_orbits)}"
)
print(
    "A4_chiral_partial_tournament="
    "outdegree2;indegree2;antipodal_ties3;"
    f"directed_triangles:{len(directed_triangles)}"
)
print("all_exact_controls=PASS")
