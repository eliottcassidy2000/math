#!/usr/bin/env python3
"""Exact finite controls for the THM-3157 six-state/action meta-pattern."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    ans = [None] * len(p)
    for i, j in enumerate(p):
        ans[j] = i
    require(all(x is not None for x in ans), "inverse input")
    return tuple(ans)


def power(p, exponent):
    require(exponent >= 0, "negative power")
    ans = tuple(range(len(p)))
    while exponent:
        if exponent & 1:
            ans = compose(p, ans)
        p = compose(p, p)
        exponent //= 2
    return ans


def parity(p):
    inversions = sum(
        p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))
    )
    return -1 if inversions % 2 else 1


def cycle_type(p):
    seen = set()
    answer = []
    for i in range(len(p)):
        if i in seen:
            continue
        here = i
        length = 0
        while here not in seen:
            seen.add(here)
            length += 1
            here = p[here]
        answer.append(length)
    return tuple(sorted(answer))


def generated_group(identity, generators, multiply):
    group = {identity}
    frontier = [identity]
    while frontier:
        g = frontier.pop()
        for generator in generators:
            for candidate in (multiply(g, generator), multiply(generator, g)):
                if candidate not in group:
                    group.add(candidate)
                    frontier.append(candidate)
    return group


# ---------------------------------------------------------------------------
# Two nonisomorphic regular six-state quotients of C2*C3 and their joint
# eighteen-state fibre product over the common binary/sign coordinate.
# ---------------------------------------------------------------------------
S3 = tuple(permutations(range(3)))
ID3 = (0, 1, 2)
S3_S = (1, 0, 2)
S3_R = (1, 2, 0)
C6_S = 3
C6_R = 2


def pair_multiply(left, right):
    return compose(left[0], right[0]), (left[1] + right[1]) % 6


JOINT = generated_group((ID3, 0), ((S3_S, C6_S), (S3_R, C6_R)), pair_multiply)
FIBRE_PRODUCT = {
    (sigma, residue)
    for sigma in S3
    for residue in range(6)
    if parity(sigma) == (-1 if residue % 2 else 1)
}
require(len(JOINT) == 18 and JOINT == FIBRE_PRODUCT, "S3 x_C2 C6 image")
require({sigma for sigma, _ in JOINT} == set(S3), "joint S3 projection")
require({residue for _, residue in JOINT} == set(range(6)), "joint C6 projection")
S3_PROJECTION_KERNEL = {
    (sigma, residue) for sigma, residue in JOINT if sigma == ID3
}
C6_PROJECTION_KERNEL = {
    (sigma, residue) for sigma, residue in JOINT if residue == 0
}
require(
    S3_PROJECTION_KERNEL == {(ID3, 0), (ID3, 2), (ID3, 4)},
    "S3-kernel C3",
)
require(
    C6_PROJECTION_KERNEL == {(ID3, 0), (S3_R, 0), (power(S3_R, 2), 0)},
    "C6-kernel C3",
)


def pair_inverse(value):
    return inverse(value[0]), (-value[1]) % 6


joint_s = (S3_S, C6_S)
joint_r = (S3_R, C6_R)
joint_identity = (ID3, 0)
require(pair_multiply(joint_s, joint_s) == joint_identity, "binary factor order")
require(joint_s != joint_identity, "binary factor injectivity")
joint_r_squared = pair_multiply(joint_r, joint_r)
require(joint_r_squared != joint_identity, "ternary factor injectivity")
require(pair_multiply(joint_r_squared, joint_r) == joint_identity, "ternary factor order")

# Since both free factors retain their exact orders, the joint kernel is
# torsion-free.  Its quotient Bass--Serre graph has one edge per joint-image
# element, one binary vertex per C2 coset, and one ternary vertex per C3 coset.
JOINT_BS_EDGES = len(JOINT)
JOINT_BS_BINARY_VERTICES = len(JOINT) // 2
JOINT_BS_TERNARY_VERTICES = len(JOINT) // 3
JOINT_FREE_KERNEL_RANK = (
    JOINT_BS_EDGES - JOINT_BS_BINARY_VERTICES - JOINT_BS_TERNARY_VERTICES + 1
)
require(
    (JOINT_BS_EDGES, JOINT_BS_BINARY_VERTICES, JOINT_BS_TERNARY_VERTICES)
    == (18, 9, 6)
    and JOINT_FREE_KERNEL_RANK == 4,
    "joint Bass--Serre free-kernel rank",
)

commutator = pair_multiply(
    pair_multiply(pair_multiply(joint_s, joint_r), pair_inverse(joint_s)),
    pair_inverse(joint_r),
)
sr_squared = pair_multiply(pair_multiply(joint_s, joint_r), pair_multiply(joint_s, joint_r))
require(commutator[0] != ID3 and commutator[1] == 0, "commutator hostile")
require(sr_squared == (ID3, 4), "(sr)^2 hostile")

S3_INDEX = {g: i for i, g in enumerate(S3)}
S3_REGULAR_ACTIONS = tuple(
    tuple(S3_INDEX[compose(g, state)] for state in S3)
    for g in S3
)
C6_REGULAR_ACTIONS = tuple(
    tuple((state + shift) % 6 for state in range(6))
    for shift in range(6)
)
S3_REGULAR_TYPES = Counter(cycle_type(action) for action in S3_REGULAR_ACTIONS)
C6_REGULAR_TYPES = Counter(cycle_type(action) for action in C6_REGULAR_ACTIONS)
require(
    S3_REGULAR_TYPES == Counter({(1, 1, 1, 1, 1, 1): 1, (2, 2, 2): 3, (3, 3): 2}),
    "regular S3 cycle inventory",
)
require(
    C6_REGULAR_TYPES
    == Counter({(1, 1, 1, 1, 1, 1): 1, (6,): 2, (3, 3): 2, (2, 2, 2): 1}),
    "regular C6 cycle inventory",
)

# ---------------------------------------------------------------------------
# Oriented pointed quartic frames form S4/C3.  Their natural graph is the
# cube K4,4 minus the owner matching, and orientation-forgetting is its
# antipodal quotient to the four owner vertices.
# ---------------------------------------------------------------------------
X = tuple(range(4))
S4 = tuple(permutations(X))
A4 = tuple(g for g in S4 if parity(g) == 1)
ID4 = tuple(range(4))
V4 = {
    g
    for g in A4
    if cycle_type(g) in ((1, 1, 1, 1), (2, 2))
}
require(len(V4) == 4 and ID4 in V4, "canonical quartic V4")
EDGES = tuple(combinations(X, 2))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
ID6 = tuple(range(6))


def edge_action(g):
    return tuple(
        EDGE_INDEX[tuple(sorted((g[a], g[b])))]
        for a, b in EDGES
    )


C = tuple(
    EDGE_INDEX[tuple(x for x in X if x not in edge)]
    for edge in EDGES
)
MATCHINGS = ((0, 5), (1, 4), (2, 3))


def matching_action(g):
    ge = edge_action(g)
    return tuple(
        next(j for j, block in enumerate(MATCHINGS) if ge[source[0]] in block)
        for source in MATCHINGS
    )


BAR_R = (1, 2, 0)
R_LIFTS = {}
for owner in X:
    candidates = [
        g for g in A4 if g[owner] == owner and matching_action(g) == BAR_R
    ]
    require(len(candidates) == 1, "pointed C3 lift")
    R_LIFTS[owner] = candidates[0]

FRAMES = tuple((owner, orientation) for owner in X for orientation in (1, -1))


def act_frame(g, frame):
    owner, orientation = frame
    return g[owner], parity(g) * orientation


base_frame = (0, 1)
frame_orbit = {act_frame(g, base_frame) for g in S4}
frame_stabilizer = {g for g in S4 if act_frame(g, base_frame) == base_frame}
require(frame_orbit == set(FRAMES), "frame transitivity")
require(
    frame_stabilizer
    == {ID4, R_LIFTS[0], power(R_LIFTS[0], 2)},
    "frame stabilizer is not C3",
)
require(
    {g for g in S4 if all(act_frame(g, frame) == frame for frame in FRAMES)} == {ID4},
    "frame action not faithful",
)

FRAME_INDEX = {frame: i for i, frame in enumerate(FRAMES)}
S4_FRAME_ACTIONS = tuple(
    tuple(FRAME_INDEX[act_frame(g, frame)] for frame in FRAMES)
    for g in S4
)
require(len(set(S4_FRAME_ACTIONS)) == 24, "S4/C3 action census")

FRAME_R = {
    frame: (R_LIFTS[frame[0]] if frame[1] == 1 else inverse(R_LIFTS[frame[0]]))
    for frame in FRAMES
}
ALL_THREE_CYCLES = {g for g in S4 if cycle_type(g) == (1, 3)}
require(set(FRAME_R.values()) == ALL_THREE_CYCLES, "frames are not the 3-cycle conjugacy class")
require(
    all(FRAME_R[(owner, -1)] == inverse(FRAME_R[(owner, 1)]) for owner in X),
    "frame antipodes are not inverse 3-cycles",
)

FRAME_EDGES = frozenset(
    tuple(sorted((FRAME_INDEX[left], FRAME_INDEX[right])))
    for left, right in combinations(FRAMES, 2)
    if left[0] != right[0] and left[1] == -right[1]
)
require(len(FRAME_EDGES) == 12, "frame cube edge census")
FRAME_NEIGHBORS = {i: set() for i in range(8)}
for i, j in FRAME_EDGES:
    FRAME_NEIGHBORS[i].add(j)
    FRAME_NEIGHBORS[j].add(i)
require(all(len(FRAME_NEIGHBORS[i]) == 3 for i in range(8)), "frame cube degree")
for left, right in combinations(FRAMES, 2):
    product_is_v4 = cycle_type(compose(FRAME_R[left], FRAME_R[right])) == (2, 2)
    is_adjacent = tuple(sorted((FRAME_INDEX[left], FRAME_INDEX[right]))) in FRAME_EDGES
    require(product_is_v4 == is_adjacent, "intrinsic 3-cycle cube adjacency")

CUBE_PLUS = {
    0: (0, 0, 0),
    1: (0, 1, 1),
    2: (1, 0, 1),
    3: (1, 1, 0),
}


def cube_coordinate(frame):
    owner, orientation = frame
    plus = CUBE_PLUS[owner]
    return plus if orientation == 1 else tuple(1 - bit for bit in plus)


def graph_distances(source, neighbors):
    distance = {source: 0}
    frontier = [source]
    while frontier:
        here = frontier.pop(0)
        for there in neighbors[here]:
            if there not in distance:
                distance[there] = distance[here] + 1
                frontier.append(there)
    return distance


for frame in FRAMES:
    i = FRAME_INDEX[frame]
    distances = graph_distances(i, FRAME_NEIGHBORS)
    require(len(distances) == 8, "frame cube connectedness")
    for other in FRAMES:
        hamming = sum(
            left != right
            for left, right in zip(cube_coordinate(frame), cube_coordinate(other))
        )
        require(distances[FRAME_INDEX[other]] == hamming, "frame Q3 embedding")


def frame_antipode(frame):
    return frame[0], -frame[1]


ANTIPODE_ACTION = tuple(FRAME_INDEX[frame_antipode(frame)] for frame in FRAMES)
require(cycle_type(ANTIPODE_ACTION) == (2, 2, 2, 2), "frame antipode")
require(
    all(
        act_frame(g, frame_antipode(frame)) == frame_antipode(act_frame(g, frame))
        for g in S4
        for frame in FRAMES
    ),
    "frame antipode not central",
)
require(ANTIPODE_ACTION not in S4_FRAME_ACTIONS, "frame antipode induced by S4")
FRAME_OCTAHEDRAL = {
    compose(power(ANTIPODE_ACTION, bit), action)
    for bit in (0, 1)
    for action in S4_FRAME_ACTIONS
}
require(len(FRAME_OCTAHEDRAL) == 48, "frame cube full extension")

# The canonical V4 sheet translations and the orientation antipode form the
# regular cube-translation group.  No base is needed for the action; choosing
# a base frame identifies the frame torsor with this group.
TRANSLATION_ACTIONS = {
    (v, bit): tuple(
        FRAME_INDEX[
            frame_antipode(act_frame(v, frame)) if bit else act_frame(v, frame)
        ]
        for frame in FRAMES
    )
    for v in V4
    for bit in (0, 1)
}
translation_values = set(TRANSLATION_ACTIONS.values())
require(len(translation_values) == 8, "frame translation order")
require(
    tuple(range(8)) in translation_values
    and all(compose(left, right) in translation_values for left in translation_values for right in translation_values),
    "frame translations do not form a group",
)
require(
    all(
        {action[FRAME_INDEX[source]] for action in TRANSLATION_ACTIONS.values()}
        == set(range(8))
        for source in FRAMES
    ),
    "frame translations are not regular",
)
NONZERO_V4 = tuple(sorted(v for v in V4 if v != ID4))
AXIS_ACTIONS = tuple(TRANSLATION_ACTIONS[(v, 1)] for v in NONZERO_V4)
require(all(cycle_type(axis) == (2, 2, 2, 2) for axis in AXIS_ACTIONS), "cube axes")
AXIS_EDGES = frozenset(
    tuple(sorted((i, axis[i])))
    for axis in AXIS_ACTIONS
    for i in range(8)
)
require(AXIS_EDGES == FRAME_EDGES, "V4 x antipode Cayley edge law")
require(
    generated_group(tuple(range(8)), AXIS_ACTIONS, compose) == translation_values,
    "cube axes do not generate the translation group",
)
axis_product = tuple(range(8))
for axis in AXIS_ACTIONS:
    axis_product = compose(axis, axis_product)
require(axis_product == ANTIPODE_ACTION, "cube axis product is not antipode")
for i, left in enumerate(NONZERO_V4):
    for right in NONZERO_V4[i + 1 :]:
        require(
            compose(TRANSLATION_ACTIONS[(left, 1)], TRANSLATION_ACTIONS[(right, 1)])
            == TRANSLATION_ACTIONS[(compose(left, right), 0)],
            "two-axis V4 product",
        )

# Exhaust the abstract cube automorphisms and identify the realized group.
CUBE_AUTOMORPHISMS = {
    p
    for p in permutations(range(8))
    if frozenset(tuple(sorted((p[i], p[j]))) for i, j in FRAME_EDGES)
    == FRAME_EDGES
}
require(len(CUBE_AUTOMORPHISMS) == 48, "abstract cube automorphism census")
require(FRAME_OCTAHEDRAL == CUBE_AUTOMORPHISMS, "S4 x antipode is not Aut(Q3)")
require(
    all(
        compose(compose(g, t), inverse(g)) in translation_values
        for g in CUBE_AUTOMORPHISMS
        for t in translation_values
    ),
    "cube translation group is not normal",
)

# The cube alone has two S4 complements to its central antipode.  The actual
# sheet action is the sign-twisted one; the plain owner action is the hostile
# alternative, and the two complements meet exactly in A4.
PLAIN_S4_ACTIONS = {
    tuple(FRAME_INDEX[(g[owner], orientation)] for owner, orientation in FRAMES)
    for g in S4
}
A4_FRAME_ACTIONS = {
    tuple(FRAME_INDEX[act_frame(g, frame)] for frame in FRAMES)
    for g in A4
}
sheet_action_set = set(S4_FRAME_ACTIONS)
require(len(PLAIN_S4_ACTIONS) == 24, "plain S4 complement order")
require(
    tuple(range(8)) in PLAIN_S4_ACTIONS
    and all(compose(left, right) in PLAIN_S4_ACTIONS for left in PLAIN_S4_ACTIONS for right in PLAIN_S4_ACTIONS)
    and tuple(range(8)) in sheet_action_set
    and all(compose(left, right) in sheet_action_set for left in sheet_action_set for right in sheet_action_set),
    "S4 complement subgroup law",
)
require(PLAIN_S4_ACTIONS != sheet_action_set, "two S4 complements collapsed")
require(
    PLAIN_S4_ACTIONS & sheet_action_set == A4_FRAME_ACTIONS
    and len(A4_FRAME_ACTIONS) == 12,
    "S4-complement intersection is not A4",
)
require(
    ANTIPODE_ACTION not in PLAIN_S4_ACTIONS
    and ANTIPODE_ACTION not in sheet_action_set,
    "antipode lies in an S4 complement",
)
PLAIN_EXTENSION = {
    compose(power(ANTIPODE_ACTION, bit), action)
    for bit in (0, 1)
    for action in PLAIN_S4_ACTIONS
}
require(PLAIN_EXTENSION == CUBE_AUTOMORPHISMS, "plain S4 complement extension")

# The sign cocycle is the exact obstruction to choosing one oriented frame
# over each owner S4-equivariantly.  It splits into two sections over A4.
def section_equivariant(signs, group):
    return all(
        act_frame(g, (owner, signs[owner])) == (g[owner], signs[g[owner]])
        for g in group
        for owner in X
    )


all_sections = tuple(dict(zip(X, signs)) for signs in product((1, -1), repeat=4))
S4_SECTIONS = tuple(signs for signs in all_sections if section_equivariant(signs, S4))
A4_SECTIONS = tuple(signs for signs in all_sections if section_equivariant(signs, A4))
require(len(S4_SECTIONS) == 0, "unexpected S4 frame section")
require(
    A4_SECTIONS == ({0: 1, 1: 1, 2: 1, 3: 1}, {0: -1, 1: -1, 2: -1, 3: -1}),
    "A4 orientation sections",
)
for g in S4:
    for k in S4:
        require(parity(compose(g, k)) == parity(g) * parity(k), "sign cocycle law")

odd_base_fixers = [g for g in S4 if g[0] == 0 and parity(g) == -1]
require(len(odd_base_fixers) == 3, "odd owner-fixer hostile census")
require(all(act_frame(g, base_frame) == (0, -1) for g in odd_base_fixers), "owner hostile")

# The THM-3157 directed hexagon and tournament selectors are equivariant
# bijections from the eight frames.  Orientation forgetting is exactly the
# frame antipode quotient.
def frame_h(frame):
    owner, orientation = frame
    r = R_LIFTS[owner] if orientation == 1 else inverse(R_LIFTS[owner])
    r_edge = edge_action(r)
    return compose(C, inverse(r_edge))


def frame_tournament(frame):
    owner, orientation = frame
    r = R_LIFTS[owner] if orientation == 1 else inverse(R_LIFTS[owner])
    r_edge = edge_action(r)
    h = frame_h(frame)
    incident = {i for i, edge in enumerate(EDGES) if owner in edge}
    arcs = set()
    for i in range(6):
        arcs.add((i, r_edge[i]))
        arcs.add((i, h[i]))
    for i in incident:
        arcs.add((i, C[i]))
    require(len(arcs) == 15, "frame tournament pair split")
    require(
        all(((i, j) in arcs) != ((j, i) in arcs) for i, j in combinations(range(6), 2)),
        "frame tournament is not a complete orientation",
    )
    return frozenset(arcs)


def conjugate(g, p):
    return compose(compose(g, p), inverse(g))


def push_arcs(g, arcs):
    return frozenset((g[i], g[j]) for i, j in arcs)


HEXAGONS = {frame: frame_h(frame) for frame in FRAMES}
TOURNAMENTS = {frame: frame_tournament(frame) for frame in FRAMES}
require(all(cycle_type(h) == (6,) for h in HEXAGONS.values()), "frame hexagon cycle type")
require(len(set(HEXAGONS.values())) == 8, "frame hexagon selector not injective")
require(len(set(TOURNAMENTS.values())) == 8, "frame tournament selector not injective")
for g in S4:
    ge = edge_action(g)
    for frame in FRAMES:
        target = act_frame(g, frame)
        require(conjugate(ge, HEXAGONS[frame]) == HEXAGONS[target], "hexagon equivariance")
        require(push_arcs(ge, TOURNAMENTS[frame]) == TOURNAMENTS[target], "tournament equivariance")

UNDIRECTED_HEXAGONS = {
    frame: frozenset(tuple(sorted((i, h[i]))) for i in range(6))
    for frame, h in HEXAGONS.items()
}
require(len(set(UNDIRECTED_HEXAGONS.values())) == 4, "owner hexagon quotient")
require(
    all(UNDIRECTED_HEXAGONS[frame] == UNDIRECTED_HEXAGONS[frame_antipode(frame)] for frame in FRAMES),
    "orientation forgetting is not antipodal quotient",
)

semantic = {
    "joint_order": len(JOINT),
    "joint_elements": tuple(sorted(JOINT)),
    "joint_S3_projection_kernel": tuple(sorted(S3_PROJECTION_KERNEL)),
    "joint_C6_projection_kernel": tuple(sorted(C6_PROJECTION_KERNEL)),
    "joint_Bass_Serre_census": (
        JOINT_BS_EDGES,
        JOINT_BS_BINARY_VERTICES,
        JOINT_BS_TERNARY_VERTICES,
        JOINT_FREE_KERNEL_RANK,
    ),
    "regular_S3_types": tuple(sorted(S3_REGULAR_TYPES.items())),
    "regular_C6_types": tuple(sorted(C6_REGULAR_TYPES.items())),
    "frame_stabilizer": tuple(sorted(frame_stabilizer)),
    "frame_three_cycles": tuple((frame, FRAME_R[frame]) for frame in FRAMES),
    "frame_edges": tuple(sorted(FRAME_EDGES)),
    "cube_coordinates": tuple((frame, cube_coordinate(frame)) for frame in FRAMES),
    "frame_translations": tuple(
        (key, action) for key, action in sorted(TRANSLATION_ACTIONS.items())
    ),
    "frame_axes": AXIS_ACTIONS,
    "cube_automorphism_order": len(CUBE_AUTOMORPHISMS),
    "S4_complement_intersection": tuple(sorted(PLAIN_S4_ACTIONS & sheet_action_set)),
    "S4_sections": len(S4_SECTIONS),
    "A4_sections": len(A4_SECTIONS),
    "hexagons": tuple((frame, HEXAGONS[frame]) for frame in FRAMES),
    "tournaments": tuple((frame, tuple(sorted(TOURNAMENTS[frame]))) for frame in FRAMES),
}
semantic_sha = sha256(repr(semantic).encode("ascii")).hexdigest()

print("THM3157_SIX_STATE_ACTION_FRAME_CUBE_META_REFEREE")
print("REGULAR_S3_TYPES=1^6:1,2^3:3,3^2:2 ORDER6=0")
print("REGULAR_C6_TYPES=1^6:1,2^3:1,3^2:2,6:2")
print(
    "JOINT_QUOTIENT=S3_x_C2_C6 ORDER=18 "
    "PROJECTION_KERNELS_IN_IMAGE=C3,C3 FREE_SOURCE_KERNEL_RANK=4"
)
print("HOSTILE_C6_BLIND=commutator S3_NONZERO=1")
print("HOSTILE_S3_BLIND=(sr)^2 C6_RESIDUE=4")
print("ORIENTED_POINTED_FRAMES=S4/C3 COUNT=8 STABILIZER=3")
print("FRAME_CLASS=ALL_8_THREE_CYCLES ADJACENCY=PRODUCT_IN_V4_NONZERO")
print("FRAME_GRAPH=Q3 EDGES=12 ANTIPODAL_OWNER_QUOTIENT=K4")
print("FRAME_ANTIPODE=CENTRAL_NOT_IN_S4 EXTENSION_ORDER=48")
print("FRAME_TRANSLATIONS=V4_x_C2 ORDER=8 REGULAR=1 EDGE_AXES=(v,a) AXIS_PRODUCT=a")
print("FRAME_AUT=Aut(Q3)=48=S4_x_C2 TWO_S4_COMPLEMENTS=1 INTERSECTION=A4")
print("ORIENTATION_COCYCLE=sign S4_SECTIONS=0 A4_SECTIONS=2")
print("SELECTOR_MAPS=8_HEXAGONS_8_TOURNAMENTS S4_EQUIVARIANT=1")
print("FORGET_ORIENTATION=4_OWNER_C_o NO_PHYSICAL_ACTION=1")
print(f"SEMANTIC_SHA256={semantic_sha}")
