#!/usr/bin/env python3
"""Finite exact referee for THM-3157's pointed edge-hexagon alignment.

This file is deliberately self-contained: it uses only permutation arithmetic
on four sheet labels, their six unordered pairs, and all 2^15 tournaments on
those pairs.  Every truth-bearing check goes through ``require`` so ``-O``
executes the same verification as an ordinary run.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def compose(p, q):
    """Return p after q."""
    require(len(p) == len(q), "composition degree mismatch")
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    ans = [None] * len(p)
    for i, j in enumerate(p):
        ans[j] = i
    require(all(x is not None for x in ans), "nonpermutation inverse")
    return tuple(ans)


def power(p, exponent):
    require(exponent >= 0, "negative permutation exponent")
    ans = tuple(range(len(p)))
    base = p
    while exponent:
        if exponent & 1:
            ans = compose(base, ans)
        base = compose(base, base)
        exponent //= 2
    return ans


def conjugate(g, p):
    return compose(compose(g, p), inverse(g))


def parity(p):
    inversions = sum(
        p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))
    )
    return -1 if inversions % 2 else 1


def cycle_type(p):
    seen = set()
    lengths = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        lengths.append(length)
    return tuple(sorted(lengths))


X = tuple(range(4))
EDGES = tuple(combinations(X, 2))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
ID4 = tuple(range(4))
ID6 = tuple(range(6))
ID3 = tuple(range(3))


def edge_action(point_perm):
    return tuple(
        EDGE_INDEX[tuple(sorted((point_perm[a], point_perm[b])))]
        for a, b in EDGES
    )


def edge_complement(i):
    complement = tuple(x for x in X if x not in EDGES[i])
    require(len(complement) == 2, "edge complement cardinality")
    return EDGE_INDEX[complement]


C = tuple(edge_complement(i) for i in range(6))
require(cycle_type(C) == (2, 2, 2), "complement cycle type")

MATCHINGS = []
seen_edges = set()
for i in range(6):
    if i in seen_edges:
        continue
    block = tuple(sorted((i, C[i])))
    MATCHINGS.append(block)
    seen_edges.update(block)
MATCHINGS = tuple(MATCHINGS)
require(MATCHINGS == ((0, 5), (1, 4), (2, 3)), "matching order")


def matching_action(point_perm):
    ep = edge_action(point_perm)
    image = []
    for block in MATCHINGS:
        target_edge = ep[block[0]]
        targets = [j for j, target in enumerate(MATCHINGS) if target_edge in target]
        require(len(targets) == 1, "matching image not unique")
        require(ep[block[1]] in MATCHINGS[targets[0]], "opposition not preserved")
        image.append(targets[0])
    return tuple(image)


S4 = tuple(permutations(X))
A4 = tuple(p for p in S4 if parity(p) == 1)
require(len(S4) == 24 and len(A4) == 12, "S4/A4 census")

S4_EDGE = tuple(edge_action(p) for p in S4)
require(len(set(S4_EDGE)) == 24, "S4 edge action not faithful")
require(all(compose(C, g) == compose(g, C) for g in S4_EDGE), "C not central")
require(C not in S4_EDGE, "complement unexpectedly induced by S4")
OCTAHEDRAL = {compose(power(C, bit), g) for bit in (0, 1) for g in S4_EDGE}
require(len(OCTAHEDRAL) == 48, "octahedral extension order")

edge_cycle_types = Counter(cycle_type(g) for g in S4_EDGE)
require(
    edge_cycle_types
    == Counter({(1, 1, 1, 1, 1, 1): 1, (1, 1, 2, 2): 9, (3, 3): 8, (2, 4): 6}),
    "S4 edge cycle-type atlas",
)

matching_actions = Counter(matching_action(p) for p in S4)
require(len(matching_actions) == 6, "resolvent image is not S3")
require(set(matching_actions.values()) == {4}, "resolvent fibres are not V4 cosets")
V4 = tuple(p for p in S4 if matching_action(p) == ID3)
require(len(V4) == 4 and all(parity(p) == 1 for p in V4), "V4 kernel")
require(
    Counter(matching_action(p) for p in A4)
    == Counter({ID3: 4, (1, 2, 0): 4, (2, 0, 1): 4}),
    "A4/C3 resolvent quotient",
)

# The chosen orientation of the three matching roots.
BAR_R = (1, 2, 0)
R_LIFTS = {}
for owner in X:
    candidates = [
        p
        for p in A4
        if p[owner] == owner and matching_action(p) == BAR_R
    ]
    require(len(candidates) == 1, f"nonunique pointed C3 lift at owner {owner}")
    R_LIFTS[owner] = candidates[0]

require(R_LIFTS[0] == (0, 2, 3, 1), "displayed owner-zero lift")


def pointed_data(owner, orientation):
    require(owner in X, "invalid owner")
    require(orientation in (1, -1), "invalid orientation")
    r_point = R_LIFTS[owner]
    if orientation == -1:
        r_point = inverse(r_point)
    r_edge = edge_action(r_point)
    h_edge = compose(C, inverse(r_edge))
    require(power(r_edge, 3) == ID6, "R order")
    require(power(h_edge, 2) == r_edge, "h^2=R")
    require(power(h_edge, 3) == C, "h^3=C")
    require(power(h_edge, 6) == ID6, "h order divides six")
    require(all(power(h_edge, k) != ID6 for k in range(1, 6)), "h order not six")
    require(cycle_type(h_edge) == (6,), "h is not a six-cycle")
    return r_point, r_edge, h_edge


example_h = pointed_data(0, 1)[2]
example_cycle = []
cursor = EDGE_INDEX[(0, 1)]
for _ in range(6):
    example_cycle.append(EDGES[cursor])
    cursor = example_h[cursor]
require(cursor == EDGE_INDEX[(0, 1)], "displayed hexagon not closed")
require(
    tuple(example_cycle)
    == ((0, 1), (1, 2), (0, 2), (2, 3), (0, 3), (1, 3)),
    "displayed hexagon changed",
)

VERTEX_PAIRS = tuple(combinations(range(6), 2))
PAIR_INDEX = {pair: i for i, pair in enumerate(VERTEX_PAIRS)}


def beats(mask, a, b):
    require(a != b, "tournament loop")
    lo, hi = sorted((a, b))
    bit = (mask >> PAIR_INDEX[(lo, hi)]) & 1
    return bool(bit) if (a, b) == (lo, hi) else not bool(bit)


def mask_from_arcs(arcs):
    arc_set = set(arcs)
    require(all(a != b for a, b in arc_set), "loop in arc set")
    mask = 0
    for pair_index, (a, b) in enumerate(VERTEX_PAIRS):
        forward = (a, b) in arc_set
        backward = (b, a) in arc_set
        require(forward != backward, f"missing or doubled tournament pair {(a, b)}")
        if forward:
            mask |= 1 << pair_index
    return mask


def invariant(mask, perm):
    return all(
        beats(mask, a, b) == beats(mask, perm[a], perm[b])
        for a, b in VERTEX_PAIRS
    )


def push_mask(mask, perm):
    arcs = []
    for a, b in VERTEX_PAIRS:
        if beats(mask, a, b):
            arcs.append((perm[a], perm[b]))
        else:
            arcs.append((perm[b], perm[a]))
    return mask_from_arcs(arcs)


def selected_tournament(owner, orientation):
    _, r_edge, h_edge = pointed_data(owner, orientation)
    incident = tuple(i for i, edge in enumerate(EDGES) if owner in edge)
    nonincident = tuple(i for i in range(6) if i not in incident)
    require(len(incident) == len(nonincident) == 3, "owner star split")
    arcs = set()

    # Both R-triangles follow the selected resolvent orientation.
    for a in incident + nonincident:
        arcs.add((a, r_edge[a]))

    # Each opposite pair is oriented from the owner star to its complement.
    for a in incident:
        arcs.add((a, C[a]))

    # The remaining six cross-pairs follow the h-hexagon.
    for a in range(6):
        arcs.add((a, h_edge[a]))

    require(len(arcs) == 15, "three tournament prescriptions overlap")
    return mask_from_arcs(arcs)


selected = {}
for owner in X:
    for orientation in (1, -1):
        mask = selected_tournament(owner, orientation)
        selected[(owner, orientation)] = mask
        _, r_edge, h_edge = pointed_data(owner, orientation)
        require(invariant(mask, r_edge), "selected tournament not C3 invariant")
        require(
            tuple(sorted(sum(beats(mask, i, j) for j in range(6) if i != j) for i in range(6)))
            == (2, 2, 2, 3, 3, 3),
            "selected score sequence",
        )
        require(
            {frozenset((i, h_edge[i])) for i in range(6)}
            == {frozenset((i, inverse(h_edge)[i])) for i in range(6)},
            "hexagon reversal changes unoriented cycle",
        )
require(len(set(selected.values())) == 8, "pointed/oriented tournament family not faithful")

# Full 32/16/4/1 census at the displayed owner and orientation.  The same
# check is repeated at all eight pointed/oriented choices below.
def tournament_census(owner, orientation):
    _, r_edge, h_edge = pointed_data(owner, orientation)
    incident = tuple(i for i, edge in enumerate(EDGES) if owner in edge)
    all_masks = range(1 << len(VERTEX_PAIRS))
    r_invariant = [mask for mask in all_masks if invariant(mask, r_edge)]
    opposite = [
        mask for mask in r_invariant if all(beats(mask, i, C[i]) for i in incident)
    ]
    triangles = [
        mask
        for mask in opposite
        if all(beats(mask, i, r_edge[i]) for i in range(6))
    ]
    hexagon = [
        mask
        for mask in triangles
        if all(beats(mask, i, h_edge[i]) for i in range(6))
    ]
    return tuple(map(len, (r_invariant, opposite, triangles, hexagon))), hexagon


for owner in X:
    for orientation in (1, -1):
        census, hexagon = tournament_census(owner, orientation)
        require(census == (32, 16, 4, 1), "32/16/4/1 census")
        require(hexagon == [selected[(owner, orientation)]], "wrong unique tournament")

selected_zero = selected[(0, 1)]
AUT = tuple(
    p for p in permutations(range(6)) if invariant(selected_zero, p)
)
_, R0_EDGE, _ = pointed_data(0, 1)
EXPECTED_AUT = {ID6, R0_EDGE, power(R0_EDGE, 2)}
require(len(AUT) == 3 and set(AUT) == EXPECTED_AUT, "automorphism group is not C3")

# The abstract tournament is inherited from THM-2597's tile mask 873.  What
# is new here is its relative placement on the six quartic edge slots: the
# displayed isomorphism carries the pointed lift R_o to the tournament's
# order-three automorphism.  It does not carry THM-2597's bicycle-support
# rotation to that automorphism (THM-2597 proves those are different).
THM2597_TILES = (
    (5, 0),
    (4, 0),
    (3, 0),
    (2, 0),
    (5, 1),
    (4, 1),
    (3, 1),
    (5, 2),
    (4, 2),
    (5, 3),
)
THM2597_TILE_EDGES = tuple((b, a) for a, b in THM2597_TILES)
thm2597_tile_index = {edge: i for i, edge in enumerate(THM2597_TILE_EDGES)}
thm2597_arcs = []
for a, b in VERTEX_PAIRS:
    tile = thm2597_tile_index.get((a, b))
    if tile is not None and 873 & (1 << tile):
        thm2597_arcs.append((a, b))
    else:
        thm2597_arcs.append((b, a))
THM2597_MASK873 = mask_from_arcs(thm2597_arcs)
EDGE_TO_THM2597 = (1, 3, 5, 4, 2, 0)
THM2597_AUT_GENERATOR = (2, 3, 4, 5, 0, 1)
require(
    push_mask(selected_zero, EDGE_TO_THM2597) == THM2597_MASK873,
    "selected tournament is not THM2597 mask 873",
)
require(
    conjugate(EDGE_TO_THM2597, R0_EDGE) == THM2597_AUT_GENERATOR,
    "pointed C3 does not align with THM2597 tournament automorphism",
)

# A4 itself cannot preserve a tournament: each opposite pair is swapped by
# an element of its V4 kernel.  Check both the local obstruction and the full
# tournament enumeration.
for i, j in MATCHINGS:
    swappers = [p for p in V4 if edge_action(p)[i] == j and edge_action(p)[j] == i]
    require(len(swappers) == 2, "opposite pair V4 swapper census")
require(
    not any(all(edge_action(p)[i] == j for i, j in MATCHINGS) for p in V4),
    "one V4 element unexpectedly swaps all three opposite pairs",
)
require(
    sum(
        all(invariant(mask, edge_action(p)) for p in A4)
        for mask in range(1 << len(VERTEX_PAIRS))
    )
    == 0,
    "unexpected unpointed A4-invariant tournament",
)

# Even relabelings preserve the resolvent orientation; odd relabelings reverse
# it.  This checks the R lift, h cycle, and selected relative tournament at all
# 96 (sheet permutation, owner) pairs.
reversal_checks = 0
for tau in S4:
    tau_edge = edge_action(tau)
    for owner in X:
        target_owner = tau[owner]
        target_orientation = parity(tau)
        r_source, _, h_source = pointed_data(owner, 1)
        r_target, _, h_target = pointed_data(target_owner, target_orientation)
        require(conjugate(tau, r_source) == r_target, "point-lift reversal law")
        require(conjugate(tau_edge, h_source) == h_target, "hexagon reversal law")
        require(
            push_mask(selected[(owner, 1)], tau_edge)
            == selected[(target_owner, target_orientation)],
            "tournament reversal law",
        )
        reversal_checks += 1
require(reversal_checks == 96, "reversal census")

# In particular, odd elements fixing the owner reverse the oriented hexagon.
odd_owner_fixers = []
for tau in S4:
    if parity(tau) == -1 and tau[0] == 0:
        tau_edge = edge_action(tau)
        require(
            conjugate(tau_edge, example_h) == inverse(example_h),
            "fixed-owner odd reversal",
        )
        odd_owner_fixers.append(tau)
require(len(odd_owner_fixers) == 3, "fixed-owner odd-reversal census")

# ---------------------------------------------------------------------------
# The six-cycle is the C6 abelianization shadow of C2*C3, not the A4 mod-three
# quotient.  In additive C6 notation s maps to 3, r maps to 2, and therefore
# s r^-1 maps to 1.  The corresponding edge actions are c, R, and h.
# ---------------------------------------------------------------------------
C6 = tuple(range(6))
s_exp = 3
r_exp = 2
sr_inverse_exp = (s_exp - r_exp) % 6
require(sr_inverse_exp == 1, "C6 image of s r^-1")
require(len({power(example_h, exponent) for exponent in C6}) == 6, "C6 action")
require(power(example_h, s_exp) == C, "abelianized s image")
require(power(example_h, r_exp) == R0_EDGE, "abelianized r image")
require(
    compose(power(example_h, s_exp), inverse(power(example_h, r_exp)))
    == example_h,
    "abelianized s r^-1 image",
)
require(compose(C, R0_EDGE) == compose(R0_EDGE, C), "abelian shadow not abelian")


def generated_group(generators):
    degree = len(generators[0])
    identity = tuple(range(degree))
    group = {identity}
    frontier = [identity]
    while frontier:
        g = frontier.pop()
        for generator in generators:
            for product in (compose(g, generator), compose(generator, g)):
                if product not in group:
                    group.add(product)
                    frontier.append(product)
    return group


S_POINT = (1, 0, 3, 2)
R_POINT = R_LIFTS[0]
require(len(generated_group((S_POINT, R_POINT))) == 12, "mod-three image not A4")
point_commutator = compose(
    compose(compose(S_POINT, R_POINT), inverse(S_POINT)), inverse(R_POINT)
)
require(point_commutator != ID4, "A4 generators unexpectedly commute")
require(
    compose(S_POINT, inverse(R_POINT)) != compose(inverse(R_POINT), S_POINT),
    "A4 forgets modular word order",
)


def matmul2(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def det2(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def inverse_sl2(a):
    require(det2(a) == 1, "matrix not in SL2")
    return ((a[1][1], -a[0][1]), (-a[1][0], a[0][0]))


def projectively_equal(a, b):
    return a == b or a == tuple(tuple(-entry for entry in row) for row in b)


S_MATRIX = ((0, -1), (1, 0))
R_MATRIX = ((0, -1), (1, 1))
require(det2(S_MATRIX) == det2(R_MATRIX) == 1, "modular generator determinant")
word_left = matmul2(S_MATRIX, inverse_sl2(R_MATRIX))
word_right = matmul2(inverse_sl2(R_MATRIX), S_MATRIX)
require(
    not projectively_equal(word_left, word_right),
    "Farey ancestry hostile collapsed in PSL2(Z)",
)
require(word_left == ((1, 0), (1, 1)), "left hostile matrix")
require(word_right == ((1, -1), (0, 1)), "right hostile matrix")

# ---------------------------------------------------------------------------
# Johnson/octahedral and partial-cube typing.  J(4,2) joins two edge slots
# when their two-subsets meet in one point.  The subdivided K4 instead joins
# singleton and two-subset vertices by inclusion and embeds isometrically in
# Q4 through characteristic vectors.
# ---------------------------------------------------------------------------
J_EDGES = tuple(
    (i, j)
    for i, j in combinations(range(6), 2)
    if len(set(EDGES[i]).intersection(EDGES[j])) == 1
)
require(len(J_EDGES) == 12, "J(4,2) edge census")
J_NEIGHBORS = {
    i: {j if i == a else a for a, j in J_EDGES if i in (a, j)}
    for i in range(6)
}
require(all(len(J_NEIGHBORS[i]) == 4 for i in range(6)), "octahedron degree")
for i in range(6):
    nonneighbors = {j for j in range(6) if j != i and j not in J_NEIGHBORS[i]}
    require(nonneighbors == {C[i]}, "complement not the Johnson antipode")
require(
    sum(
        all(tuple(sorted(pair)) in J_EDGES for pair in combinations(triple, 2))
        for triple in combinations(range(6), 3)
    )
    == 8,
    "octahedral triangular-face census",
)
hexagon_pairs = {
    tuple(sorted((i, example_h[i])))
    for i in range(6)
}
require(len(hexagon_pairs) == 6, "hexagon edge census")
require(hexagon_pairs.issubset(set(J_EDGES)), "h is not a Johnson Hamilton cycle")

# THM-2606 already identifies four owner-selected antipode-compatible C6
# subgraphs C_o.  Verify that h adds precisely an orientation to C_o, rather
# than producing a fifth underlying hexagon.
OWNER_CYCLES = {}
for owner in X:
    incident = {i for i, edge in enumerate(EDGES) if owner in edge}
    owner_cycle = {
        (i, j)
        for i, j in J_EDGES
        if (i in incident) != (j in incident)
    }
    require(len(owner_cycle) == 6, "THM2606 owner C6 edge census")
    require(
        all(sum(vertex in edge for edge in owner_cycle) == 2 for vertex in range(6)),
        "THM2606 owner C6 degree",
    )
    OWNER_CYCLES[owner] = frozenset(owner_cycle)
    for orientation in (1, -1):
        h_owner = pointed_data(owner, orientation)[2]
        h_pairs = frozenset(tuple(sorted((i, h_owner[i]))) for i in range(6))
        require(h_pairs == OWNER_CYCLES[owner], "h does not orient THM2606 C_o")
require(len(set(OWNER_CYCLES.values())) == 4, "THM2606 four-owner C6 census")

# THM-2632's Farey residue C6 is the Cayley graph of S3 by two adjacent
# transpositions.  It is abstractly a hexagon, but no modular left action is
# a one-step six-cycle; that cyclic action belongs to abelianization and needs
# a noncanonical gluing (THM-2641) unless the owner and C3 orientation are
# supplied as above.
S3 = tuple(permutations(range(3)))
FAREY_L = (0, 2, 1)
FAREY_R = (1, 0, 2)
FAREY_EDGES = frozenset(
    tuple(sorted((g, compose(g, generator))))
    for g in S3
    for generator in (FAREY_L, FAREY_R)
)
require(len(FAREY_EDGES) == 6, "Farey C6 edge census")
FAREY_NEIGHBORS = {g: set() for g in S3}
for g, k in FAREY_EDGES:
    FAREY_NEIGHBORS[g].add(k)
    FAREY_NEIGHBORS[k].add(g)
require(all(len(FAREY_NEIGHBORS[g]) == 2 for g in S3), "Farey C6 degree")
farey_seen = {S3[0]}
farey_frontier = [S3[0]]
while farey_frontier:
    farey_here = farey_frontier.pop()
    for farey_there in FAREY_NEIGHBORS[farey_here]:
        if farey_there not in farey_seen:
            farey_seen.add(farey_there)
            farey_frontier.append(farey_there)
require(len(farey_seen) == 6, "Farey C6 connectedness")
require(all(cycle_type(g) in ((1, 1, 1), (1, 2), (3,)) for g in S3), "S3 orders")
require(not any(cycle_type(g) == (6,) for g in S3), "S3 has an order-six element")
S3_INDEX = {g: i for i, g in enumerate(S3)}
FAREY_LEFT_ACTIONS = tuple(
    tuple(S3_INDEX[compose(g, state)] for state in S3)
    for g in S3
)
require(
    not any(cycle_type(action) == (6,) for action in FAREY_LEFT_ACTIONS),
    "Farey modular left action contains a one-step six-cycle",
)

owner_zero_edges = OWNER_CYCLES[0]
farey_isomorphisms = 0
for images in permutations(range(6)):
    image = {g: images[i] for i, g in enumerate(S3)}
    if {
        tuple(sorted((image[g], image[k])))
        for g, k in FAREY_EDGES
    } == set(owner_zero_edges):
        farey_isomorphisms += 1
require(farey_isomorphisms == 12, "noncanonical Farey/owner-C6 gluing census")

INCIDENCE_VERTICES = tuple(frozenset((i,)) for i in X) + tuple(
    frozenset(edge) for edge in EDGES
)
INCIDENCE_INDEX = {vertex: i for i, vertex in enumerate(INCIDENCE_VERTICES)}
INCIDENCE_EDGES = tuple(
    (i, j)
    for i, j in combinations(range(len(INCIDENCE_VERTICES)), 2)
    if len(INCIDENCE_VERTICES[i]) != len(INCIDENCE_VERTICES[j])
    and (
        INCIDENCE_VERTICES[i].issubset(INCIDENCE_VERTICES[j])
        or INCIDENCE_VERTICES[j].issubset(INCIDENCE_VERTICES[i])
    )
)
require(len(INCIDENCE_VERTICES) == 10 and len(INCIDENCE_EDGES) == 12, "subdivided K4 census")
INCIDENCE_NEIGHBORS = {i: set() for i in range(10)}
for i, j in INCIDENCE_EDGES:
    INCIDENCE_NEIGHBORS[i].add(j)
    INCIDENCE_NEIGHBORS[j].add(i)


def graph_distances(source, neighbors):
    distances = {source: 0}
    frontier = [source]
    while frontier:
        here = frontier.pop(0)
        for there in neighbors[here]:
            if there not in distances:
                distances[there] = distances[here] + 1
                frontier.append(there)
    return distances


for i, vertex in enumerate(INCIDENCE_VERTICES):
    distances = graph_distances(i, INCIDENCE_NEIGHBORS)
    require(len(distances) == 10, "incidence graph disconnected")
    for j, other in enumerate(INCIDENCE_VERTICES):
        require(
            distances[j] == len(vertex.symmetric_difference(other)),
            "Q4 embedding is not isometric",
        )

theta_classes = Counter()
for i, j in INCIDENCE_EDGES:
    difference = INCIDENCE_VERTICES[i].symmetric_difference(INCIDENCE_VERTICES[j])
    require(len(difference) == 1, "incidence edge changes multiple cube coordinates")
    theta_classes[next(iter(difference))] += 1
require(theta_classes == Counter({0: 3, 1: 3, 2: 3, 3: 3}), "Q4 Theta classes")
for i in range(6):
    j = example_h[i]
    require(
        len(set(EDGES[i]).symmetric_difference(EDGES[j])) == 2,
        "Johnson h-step Hamming distance",
    )
    left = INCIDENCE_INDEX[frozenset(EDGES[i])]
    right = INCIDENCE_INDEX[frozenset(EDGES[j])]
    require(tuple(sorted((left, right))) not in INCIDENCE_EDGES, "h is a partial-cube edge")

# The complement operator gives only a Z/2-graded representation: the three
# opposite-pair sums and differences are its +1 and -1 eigenspaces.  All
# operations generated here preserve the grading, so no Lie-superalgebra or
# supergroup structure is being asserted.
def act_on_vector(perm, vector):
    answer = [0] * len(vector)
    for i, coefficient in enumerate(vector):
        answer[perm[i]] += coefficient
    return tuple(answer)


PLUS_BASIS = []
MINUS_BASIS = []
for i, j in MATCHINGS:
    plus = tuple(int(k == i) + int(k == j) for k in range(6))
    minus = tuple(int(k == i) - int(k == j) for k in range(6))
    require(act_on_vector(C, plus) == plus, "plus eigenspace")
    require(act_on_vector(C, minus) == tuple(-x for x in minus), "minus eigenspace")
    PLUS_BASIS.append(plus)
    MINUS_BASIS.append(minus)
require(len(PLUS_BASIS) == len(MINUS_BASIS) == 3, "complement eigenspace dimensions")
for g in tuple(S4_EDGE) + (C, example_h):
    require(
        all(act_on_vector(C, act_on_vector(g, v)) == act_on_vector(g, v) for v in PLUS_BASIS),
        "plus grading not preserved",
    )
    require(
        all(
            act_on_vector(C, act_on_vector(g, v))
            == tuple(-x for x in act_on_vector(g, v))
            for v in MINUS_BASIS
        ),
        "minus grading not preserved",
    )

semantic = {
    "edges": EDGES,
    "matchings": MATCHINGS,
    "owner_zero_R": R_LIFTS[0],
    "owner_zero_hexagon": tuple(example_cycle),
    "tournament_census": (32, 16, 4, 1),
    "scores": (3, 3, 3, 2, 2, 2),
    "aut_order": len(AUT),
    "thm2597_alignment": (EDGE_TO_THM2597, THM2597_AUT_GENERATOR),
    "pointed_oriented_tournaments": len(set(selected.values())),
    "reversal_checks": reversal_checks,
    "abelianization": (s_exp, r_exp, sr_inverse_exp),
    "farey_hostile": (word_left, word_right),
    "johnson_edges": J_EDGES,
    "owner_cycles": tuple(sorted(tuple(sorted(cycle)) for cycle in OWNER_CYCLES.values())),
    "farey_isomorphisms": farey_isomorphisms,
    "incidence_theta_classes": tuple(sorted(theta_classes.items())),
    "complement_eigenspaces": (len(PLUS_BASIS), len(MINUS_BASIS)),
}
semantic_sha = sha256(repr(semantic).encode("ascii")).hexdigest()

print("THM3157_POINTED_RESOLVENT_EDGE_HEXAGON_ABELIANIZATION_REFEREE")
print("S4_ACTIONS=24 A4_ACTIONS=12 V4_KERNEL=4 RESOLVENT_ACTIONS=6")
print("S4_EDGE_CYCLE_TYPES=1^6:1,1^2_2^2:9,3^2:8,2_4:6")
print("COMPLEMENT_TYPE=2^3 COMPLEMENT_IN_S4=0 OCTAHEDRAL_EXTENSION=48")
print("POINTED_C3_LIFTS=4 ORIENTED_POINTED_FRAMES=8")
print("OWNER0_R=(1 2 3)")
print("OWNER0_H=01>12>02>23>03>13>01")
print("HEXAGON_IDENTITIES=h^2=R,h^3=c,h^6=1")
print("TOURNAMENT_CENSUS=32,16,4,1")
print("TOURNAMENT_SCORES=3,3,3,2,2,2 AUT=C3")
print("TOURNAMENT_TYPE=INHERITED_THM2597_MASK873 RELATIVE_ALIGNMENT=NEW")
print("UNPOINTED_A4_INVARIANT_TOURNAMENTS=0")
print("REVERSAL_CHECKS=96 FIXED_OWNER_ODD_REVERSALS=3")
print("MODULAR_SHADOW=C6_ABELIANIZATION s=h^3 r=h^2 sr^-1=h")
print("NONABELIAN_CONTRAST=A4_ORDER12 COMMUTATOR_NONTRIVIAL=1")
print("FAREY_HOSTILE=sr^-1!=r^-1s SAME_C6_IMAGE=h")
print("INHERITED_C6=THM2606_OWNER_C_o COUNT=4 h_ORIENTS_C_o=1")
print("FAREY_C6_GLUE=ABSTRACT_ISOMORPHISMS_12 ACTION_COMPATIBLE_ONE_STEP=0")
print("JOHNSON=J(4,2) V=6 E=12 TRIANGLES=8 COMPLEMENT=ANTIPODE")
print("PARTIAL_CUBE=SUBDIVIDED_K4_IN_Q4 V=10 E=12 THETA=3,3,3,3")
print("H_TYPE=JOHNSON_HAMILTON_CYCLE Q4_STEP_DISTANCE=2")
print("COMPLEMENT_EIGENSPACES=3+3 Z2_GRADED_REPRESENTATION_ONLY=1")
print("DESCENT=unoriented_hexagon_only OWNER_AND_C3_ORIENTATION_REQUIRED=1")
print("SCOPE=LOCAL_CONDITIONAL NO_KELLER_OR_LRC_CLOSURE=1")
print(f"SEMANTIC_SHA256={semantic_sha}")
