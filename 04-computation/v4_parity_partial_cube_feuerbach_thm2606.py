#!/usr/bin/env python3
"""Exact companion for THM-2606 (dependency-free, integer only)."""

from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


V = ((0, 0), (0, 1), (1, 0), (1, 1))
ZERO = (0, 0)
NONZERO = tuple(x for x in V if x != ZERO)


def add(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def dot(u, x):
    return (u[0] * x[0]) ^ (u[1] * x[1])


def mat_apply(matrix, x):
    return (
        (matrix[0][0] * x[0]) ^ (matrix[0][1] * x[1]),
        (matrix[1][0] * x[0]) ^ (matrix[1][1] * x[1]),
    )


def determinant(matrix):
    return (matrix[0][0] * matrix[1][1]) ^ (
        matrix[0][1] * matrix[1][0]
    )


def graph_image(edges, permutation):
    return frozenset(
        tuple(sorted((permutation[i], permutation[j]))) for i, j in edges
    )


def degrees(n, edges):
    return tuple(sum(vertex in edge for edge in edges) for vertex in range(n))


def connected(n, edges):
    if n == 0:
        return True
    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for edge in edges:
            if vertex not in edge:
                continue
            other = edge[0] if edge[1] == vertex else edge[1]
            if other not in seen:
                seen.add(other)
                frontier.append(other)
    return len(seen) == n


def bipartite(n, edges):
    adjacency = [set() for _ in range(n)]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    colors = {}
    for start in range(n):
        if start in colors:
            continue
        colors[start] = 0
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            for other in adjacency[vertex]:
                if other not in colors:
                    colors[other] = 1 - colors[vertex]
                    frontier.append(other)
                elif colors[other] == colors[vertex]:
                    return False
    return True


def graceful(labels, edges):
    return {
        abs(labels[i] - labels[j]) for i, j in edges
    } == set(range(1, len(edges) + 1))


def canonical_cycle_edges(order):
    return frozenset(
        tuple(sorted((order[i], order[(i + 1) % len(order)])))
        for i in range(len(order))
    )


GL = []
for bits in product(range(2), repeat=4):
    matrix = (bits[:2], bits[2:])
    if determinant(matrix):
        GL.append(matrix)
require(len(GL) == 6, "GL(2,F2) order changed")

GL_PERMS = tuple(
    tuple(V.index(mat_apply(matrix, x)) for x in V) for matrix in GL
)
TRANSLATIONS = tuple(
    tuple(V.index(add(x, shift)) for x in V) for shift in V
)
AGL_PERMS = tuple(
    dict.fromkeys(
        tuple(V.index(add(mat_apply(matrix, x), shift)) for x in V)
        for matrix in GL
        for shift in V
    )
)
require(len(AGL_PERMS) == 24, "AGL(2,F2) order changed")

POINT_EDGES = tuple(combinations(range(4), 2))
PARITY_GRAPHS = []
PARITY_KERNELS = []
for functional in NONZERO:
    kernel = tuple(x for x in NONZERO if dot(functional, x) == 0)
    edges = frozenset(
        (i, j)
        for i, j in POINT_EDGES
        if dot(functional, add(V[i], V[j])) == 1
    )
    require(len(kernel) == 1, "nonzero parity kernel is not a line")
    require(len(edges) == 4 and degrees(4, edges) == (2, 2, 2, 2),
            "parity graph is not C4")
    require(connected(4, edges) and bipartite(4, edges),
            "parity graph is not a connected partial cube control")
    missing = frozenset(POINT_EDGES) - edges
    expected = frozenset(
        tuple(sorted((i, V.index(add(V[i], kernel[0]))))) for i in range(4)
    )
    require(missing == expected, "omitted direction is not the diagonal matching")
    PARITY_GRAPHS.append(edges)
    PARITY_KERNELS.append(kernel[0])
require(len(set(PARITY_GRAPHS)) == 3, "parity channel count changed")
require(set(PARITY_KERNELS) == set(NONZERO), "kernel-direction bijection failed")

# The GL action on nonzero covectors is the full S3 action.
parity_signatures = tuple(tuple(dot(u, x) for x in V) for u in NONZERO)
parity_action = set()
for matrix in GL:
    image = []
    for u in NONZERO:
        signature = tuple(dot(u, mat_apply(matrix, x)) for x in V)
        image.append(parity_signatures.index(signature))
    parity_action.add(tuple(image))
require(len(parity_action) == 6, "GL action on parity channels is not S3")

# Explicit marked generators: S has order two, T order three.
S = ((0, 1), (1, 0))
T = ((0, 1), (1, 1))


def induced_parity_permutation(matrix):
    result = []
    for u in NONZERO:
        signature = tuple(dot(u, mat_apply(matrix, x)) for x in V)
        result.append(parity_signatures.index(signature))
    return tuple(result)


S_ACTION = induced_parity_permutation(S)
T_ACTION = induced_parity_permutation(T)
require(sorted(S_ACTION) == [0, 1, 2] and sorted(T_ACTION) == [0, 1, 2],
        "marked parity action is not a permutation")
require(sum(S_ACTION[i] == i for i in range(3)) == 1,
        "order-two element does not fix exactly one parity")
require(sum(T_ACTION[i] == i for i in range(3)) == 0,
        "order-three element does not cycle the parities")

# Full-affine and pointed-linear invariant graph classifications.
all_point_graphs = []
for mask in range(1 << len(POINT_EDGES)):
    graph = frozenset(
        POINT_EDGES[i] for i in range(len(POINT_EDGES)) if (mask >> i) & 1
    )
    all_point_graphs.append(graph)
affine_invariant = tuple(
    graph
    for graph in all_point_graphs
    if all(graph_image(graph, permutation) == graph for permutation in AGL_PERMS)
)
linear_invariant = tuple(
    graph
    for graph in all_point_graphs
    if all(graph_image(graph, permutation) == graph for permutation in GL_PERMS)
)
require(sorted(map(len, affine_invariant)) == [0, 6],
        "full-affine invariant graph classification changed")
require(sorted(map(len, linear_invariant)) == [0, 3, 3, 6],
        "pointed-linear invariant graph classification changed")
linear_connected_bipartite = tuple(
    graph for graph in linear_invariant if connected(4, graph) and bipartite(4, graph)
)
require(len(linear_connected_bipartite) == 1, "pointed invariant partial cube not unique")
STAR = linear_connected_bipartite[0]
require(degrees(4, STAR) == (3, 1, 1, 1), "unique pointed graph is not K1,3")

# No translation-invariant tournament on the four affine points.
fixed_point_tournaments = 0
for mask in range(1 << len(POINT_EDGES)):
    good = True
    for permutation in TRANSLATIONS:
        for edge_index, (i, j) in enumerate(POINT_EDGES):
            direction = (mask >> edge_index) & 1
            ii, jj = permutation[i], permutation[j]
            if ii < jj:
                image_index = POINT_EDGES.index((ii, jj))
                image_direction = direction
            else:
                image_index = POINT_EDGES.index((jj, ii))
                image_direction = 1 - direction
            if ((mask >> image_index) & 1) != image_direction:
                good = False
                break
        if not good:
            break
    fixed_point_tournaments += int(good)
require(fixed_point_tournaments == 0, "translation-invariant tournament appeared")

# Graceful positive controls.
K4 = frozenset(POINT_EDGES)
C4 = frozenset(((0, 1), (1, 2), (2, 3), (0, 3)))
require(graceful((0, 1, 4, 6), K4), "K4 graceful control failed")
require(graceful((0, 4, 1, 2), C4), "C4 graceful control failed")
require(graceful((0, 1, 2, 3), STAR), "K1,3 graceful control failed")

# Six unordered point-pairs and the octahedral/disjoint association scheme.
EDGE_VERTICES = POINT_EDGES
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGE_VERTICES)}
EDGE_ACTIONS = tuple(
    tuple(
        EDGE_INDEX[tuple(sorted((permutation[i], permutation[j])))]
        for i, j in EDGE_VERTICES
    )
    for permutation in permutations(range(4))
)
SIX_PAIRS = tuple(combinations(range(6), 2))
incident_pairs = frozenset(
    (i, j)
    for i, j in SIX_PAIRS
    if set(EDGE_VERTICES[i]) & set(EDGE_VERTICES[j])
)
disjoint_pairs = frozenset(SIX_PAIRS) - incident_pairs
require(len(incident_pairs) == 12 and len(disjoint_pairs) == 3,
        "six-edge association-scheme orbit sizes changed")
require(degrees(6, incident_pairs) == (4,) * 6,
        "intersection graph is not octahedral")
require(degrees(6, disjoint_pairs) == (1,) * 6,
        "disjointness graph is not a perfect matching")

fixed_six_tournaments = 0
for mask in range(1 << len(SIX_PAIRS)):
    good = True
    for permutation in EDGE_ACTIONS:
        for pair_index, (i, j) in enumerate(SIX_PAIRS):
            direction = (mask >> pair_index) & 1
            ii, jj = permutation[i], permutation[j]
            if ii < jj:
                image_index = SIX_PAIRS.index((ii, jj))
                image_direction = direction
            else:
                image_index = SIX_PAIRS.index((jj, ii))
                image_direction = 1 - direction
            if ((mask >> image_index) & 1) != image_direction:
                good = False
                break
        if not good:
            break
    fixed_six_tournaments += int(good)
require(fixed_six_tournaments == 0, "S4-invariant six-edge tournament appeared")

# Four origin-indexed equatorial C6 graphs.
opposite = {
    i: EDGE_INDEX[tuple(sorted(set(range(4)) - set(edge)))]
    for i, edge in enumerate(EDGE_VERTICES)
}
origin_cycles = set()
for origin in range(4):
    graph = set()
    for middle in range(4):
        if middle == origin:
            continue
        for last in range(4):
            if last in (origin, middle):
                continue
            first_edge = EDGE_INDEX[tuple(sorted((origin, middle)))]
            second_edge = EDGE_INDEX[tuple(sorted((middle, last)))]
            graph.add(tuple(sorted((first_edge, second_edge))))
    graph = frozenset(graph)
    require(len(graph) == 6 and degrees(6, graph) == (2,) * 6,
            "origin graph is not 2-regular")
    require(connected(6, graph) and bipartite(6, graph),
            "origin graph is not a C6 partial cube control")
    origin_cycles.add(graph)
require(len(origin_cycles) == 4, "origin-indexed C6 count changed")

all_opposition_cycles = set()
for order in permutations(range(6)):
    if all(order[i + 3] == opposite[order[i]] for i in range(3)):
        all_opposition_cycles.add(canonical_cycle_edges(order))
require(all_opposition_cycles == origin_cycles,
        "origin-to-opposition-compatible-C6 bijection failed")

# C6 is not graceful: direct exhaustive control, matching the parity proof.
one_c6 = next(iter(origin_cycles))
c6_graceful_count = 0
for labels in permutations(range(7), 6):
    c6_graceful_count += int(graceful(labels, one_c6))
require(c6_graceful_count == 0, "C6 unexpectedly became graceful")

# Feuerbach sign labels: sign triples modulo global reversal form V4.
SIGNS = tuple(product((-1, 1), repeat=3))


def sign_class(signs):
    return signs if signs[0] == 1 else tuple(-entry for entry in signs)


SIGN_CLASSES = tuple(sorted({sign_class(signs) for signs in SIGNS}))
require(len(SIGN_CLASSES) == 4, "projective sign quotient is not V4-sized")
IDENTITY_SIGN = (1, 1, 1)
require(IDENTITY_SIGN in SIGN_CLASSES, "all-positive Feuerbach class missing")
coordinate_actions = set()
for permutation in permutations(range(3)):
    action = tuple(
        SIGN_CLASSES.index(sign_class(tuple(signs[i] for i in permutation)))
        for signs in SIGN_CLASSES
    )
    coordinate_actions.add(action)
require(len(coordinate_actions) == 6, "side-permutation action is not S3")
identity_index = SIGN_CLASSES.index(IDENTITY_SIGN)
require(all(action[identity_index] == identity_index for action in coordinate_actions),
        "side permutations do not fix the incenter sign")

# Exact integer controls for the four-body/resolvent/Jacobi formulas.
body_controls = 0
ordered_controls = 0
for x1, x2, x3 in product(range(-4, 5), repeat=3):
    x4 = -(x1 + x2 + x3)
    roots = (x1, x2, x3, x4)
    p = sum(roots[i] * roots[j] for i, j in combinations(range(4), 2))
    e3 = sum(
        roots[i] * roots[j] * roots[k]
        for i, j, k in combinations(range(4), 3)
    )
    q = -e3
    r = x1 * x2 * x3 * x4
    s_values = (x1 + x2, x1 + x3, x1 + x4)
    u_values = tuple(value * value for value in s_values)
    sum_u = sum(u_values)
    pair_u = sum(u_values[i] * u_values[j] for i, j in combinations(range(3), 2))
    product_u = u_values[0] * u_values[1] * u_values[2]
    require(sum_u == -2 * p, "resolvent U^2 coefficient failed")
    require(pair_u == p * p - 4 * r, "resolvent U coefficient failed")
    require(product_u == q * q, "resolvent constant coefficient failed")
    require(s_values[0] * s_values[1] * s_values[2] == -q,
            "pairing orientation product failed")
    X, Y, Z = (2 * value for value in s_values)
    require(X * X + Y * Y + Z * Z == 4 * sum(value * value for value in roots),
            "Hadamard/Jacobi norm identity failed")
    body_controls += 1
    if x1 < x2 < x3 < x4 and q != 0:
        require(u_values[0] > u_values[1] > u_values[2] > 0,
                "ordered real-root resolvent channel order failed")
        ordered_controls += 1
require(ordered_controls > 0, "ordered real-root hostile/positive controls empty")

print("THM-2606 AFFINE V4 / PARTIAL-CUBE / FEUERBACH AUDIT")
print(f"GL2 order={len(GL)}, AGL2 order={len(AGL_PERMS)}, parity channels={len(PARITY_GRAPHS)}")
print(f"marked parity actions: C2={S_ACTION}, C3={T_ACTION}")
print(f"AGL-invariant graph edge counts={sorted(map(len, affine_invariant))}")
print(f"pointed GL-invariant graph edge counts={sorted(map(len, linear_invariant))}")
print(f"translation-invariant tournaments on X={fixed_point_tournaments}")
print("graceful controls: K4=PASS, C4=PASS, K1,3=PASS")
print(f"six-edge orbitals: incident={len(incident_pairs)}, disjoint={len(disjoint_pairs)}")
print(f"S4-invariant tournaments on six edges={fixed_six_tournaments}")
print(f"opposition-compatible C6 graphs={len(all_opposition_cycles)}=origins")
print(f"graceful labelings of C6={c6_graceful_count}")
print(f"Feuerbach sign classes={len(SIGN_CLASSES)}, side action order={len(coordinate_actions)}")
print(f"four-body integer controls={body_controls}, ordered controls={ordered_controls}")
print("FAILED CHECKS: NONE")
