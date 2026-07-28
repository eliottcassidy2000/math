#!/usr/bin/env python3
"""Exact finite referee for THM-2780.

The script uses only integer/set arithmetic and explicit exceptions.  It
checks that the four marked D3 determinant-weight colourings are precisely
the four-point V4 torsor, classifies their S4/A4 stabilizers, identifies the
three nonzero even-sign inertia words on marks and root lines, and computes
the complete orbit-symmetrized pair spectrum.
"""

from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


POINTS = tuple(range(4))
EDGES = tuple(combinations(POINTS, 2))
EDGE_PAIRS = tuple(combinations(EDGES, 2))
S4 = tuple(permutations(POINTS))


def parity(perm):
    inversions = sum(
        1
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
        if perm[i] > perm[j]
    )
    return inversions % 2


A4 = tuple(perm for perm in S4 if parity(perm) == 0)


def act_point(perm, point):
    return perm[point]


def act_edge(perm, edge):
    return tuple(sorted((perm[edge[0]], perm[edge[1]])))


def are_opposite(edge_a, edge_b):
    return set(edge_a).isdisjoint(edge_b)


def marked_weight(mark, edge_a, edge_b):
    if are_opposite(edge_a, edge_b):
        return 2
    if mark not in edge_a and mark not in edge_b:
        return 3
    return 1


def colouring(mark):
    return tuple(marked_weight(mark, edge_a, edge_b) for edge_a, edge_b in EDGE_PAIRS)


COLOURINGS = tuple(colouring(mark) for mark in POINTS)
require(len(set(COLOURINGS)) == 4, "the four marked colourings must be distinct")


def colouring_stabilizer(group, mark):
    result = []
    for perm in group:
        good = True
        for edge_a, edge_b in EDGE_PAIRS:
            lhs = marked_weight(mark, act_edge(perm, edge_a), act_edge(perm, edge_b))
            rhs = marked_weight(mark, edge_a, edge_b)
            if lhs != rhs:
                good = False
                break
        if good:
            result.append(perm)
    return tuple(result)


for mark in POINTS:
    weight_histogram = {
        value: colouring(mark).count(value)
        for value in (1, 2, 3)
    }
    require(weight_histogram == {1: 9, 2: 3, 3: 3}, "wrong marked spectrum")

    weight_three_pairs = tuple(
        pair
        for pair in EDGE_PAIRS
        if marked_weight(mark, pair[0], pair[1]) == 3
    )
    triangle_edges = tuple(edge for edge in EDGES if mark not in edge)
    require(
        set(edge for pair in weight_three_pairs for edge in pair) == set(triangle_edges),
        "weight-three support must be the opposite triangle",
    )
    omitted = set(POINTS) - set(vertex for edge in triangle_edges for vertex in edge)
    require(omitted == {mark}, "weight-three support must reconstruct the mark")

    s4_stabilizer = colouring_stabilizer(S4, mark)
    a4_stabilizer = colouring_stabilizer(A4, mark)
    require(len(s4_stabilizer) == 6, "S4 colouring stabilizer must have order six")
    require(len(a4_stabilizer) == 3, "A4 colouring stabilizer must have order three")
    require(
        {perm for perm in S4 if perm[mark] == mark} == set(s4_stabilizer),
        "S4 colouring stabilizer must equal the point stabilizer",
    )
    require(
        {perm for perm in A4 if perm[mark] == mark} == set(a4_stabilizer),
        "A4 colouring stabilizer must equal the point stabilizer",
    )

for perm in S4:
    for mark in POINTS:
        for edge_a, edge_b in EDGE_PAIRS:
            require(
                marked_weight(
                    act_point(perm, mark),
                    act_edge(perm, edge_a),
                    act_edge(perm, edge_b),
                )
                == marked_weight(mark, edge_a, edge_b),
                "marked weighting must be S4-equivariant",
            )


def compose(left, right):
    return tuple(left[right[i]] for i in POINTS)


IDENTITY = POINTS
V4 = tuple(
    perm
    for perm in A4
    if compose(perm, perm) == IDENTITY
)
require(len(V4) == 4, "V4 must have four elements")

for source in POINTS:
    images = {perm[source] for perm in V4}
    require(images == set(POINTS), "V4 must act transitively on marks")
    for target in POINTS:
        count = sum(1 for perm in V4 if perm[source] == target)
        require(count == 1, "V4 action on marks must be simply transitive")

for mark in POINTS:
    require(
        sum(1 for perm in V4 if perm in colouring_stabilizer(V4, mark)) == 1,
        "one colouring must have trivial V4 stabilizer",
    )


H_STATES = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)
ROOT_LINES = {
    "a12": (1, -1, 0),
    "b12": (1, 1, 0),
    "a13": (1, 0, -1),
    "b13": (1, 0, 1),
    "a23": (0, 1, -1),
    "b23": (0, 1, 1),
}
ROOT_EDGE = {
    "a12": (1, 2),
    "b12": (0, 3),
    "a13": (1, 3),
    "b13": (0, 2),
    "a23": (2, 3),
    "b23": (0, 1),
}
EXPECTED_FIXED = {
    "110": {"a12", "b12"},
    "101": {"a13", "b13"},
    "011": {"a23", "b23"},
}
SIGN_WORDS = {
    "000": (1, 1, 1),
    "110": (-1, -1, 1),
    "101": (-1, 1, -1),
    "011": (1, -1, -1),
}


def canonical_line(vector):
    for entry in vector:
        if entry:
            if entry < 0:
                return tuple(-value for value in vector)
            return tuple(vector)
    raise RuntimeError("zero vector is not a root line")


LINE_TO_NAME = {canonical_line(vector): name for name, vector in ROOT_LINES.items()}


def determinant3(row_a, row_b, row_c):
    return (
        row_a[0] * (row_b[1] * row_c[2] - row_b[2] * row_c[1])
        - row_a[1] * (row_b[0] * row_c[2] - row_b[2] * row_c[0])
        + row_a[2] * (row_b[0] * row_c[1] - row_b[1] * row_c[0])
    )


def mod_two_right_kernel(rows):
    kernel = []
    for vector in (
        (x_0, x_1, x_2)
        for x_0 in (0, 1)
        for x_1 in (0, 1)
        for x_2 in (0, 1)
    ):
        if all(sum(row[i] * vector[i] for i in range(3)) % 2 == 0 for row in rows):
            kernel.append(vector)
    return tuple(kernel)


WORD_VECTOR = {
    word: tuple(0 if sign == 1 else 1 for sign in signs)
    for word, signs in SIGN_WORDS.items()
}
PAIR_TO_WORD = {
    frozenset(names): word
    for word, names in EXPECTED_FIXED.items()
}

root_names = tuple(ROOT_LINES)
frame_kernel_histogram = {1: 0, 2: 0, 3: 0}
for mark, h_state in enumerate(H_STATES):
    for name_a, name_b in combinations(root_names, 2):
        edge_a = ROOT_EDGE[name_a]
        edge_b = ROOT_EDGE[name_b]
        expected_weight = marked_weight(mark, edge_a, edge_b)
        determinant_weight = abs(
            determinant3(ROOT_LINES[name_a], ROOT_LINES[name_b], h_state)
        )
        require(
            determinant_weight == expected_weight,
            "root determinant and K4 marked weight disagree",
        )
        kernel = mod_two_right_kernel(
            (ROOT_LINES[name_a], ROOT_LINES[name_b], h_state)
        )
        if expected_weight == 2:
            word = PAIR_TO_WORD.get(frozenset((name_a, name_b)))
            require(word is not None, "weight-two frame is not an inertia pair")
            require(
                set(kernel) == {(0, 0, 0), WORD_VECTOR[word]},
                "weight-two mod-two kernel is not its inertia word",
            )
        else:
            require(kernel == ((0, 0, 0),), "odd frame must be invertible mod two")
        frame_kernel_histogram[expected_weight] += 1

require(
    frame_kernel_histogram == {1: 36, 2: 12, 3: 12},
    "wrong four-mark frame census",
)


def sign_action_on_states(signs):
    action = []
    for state in H_STATES:
        image = tuple(signs[i] * state[i] for i in range(3))
        require(image in H_STATES, "even sign word left the half-Hadamard states")
        action.append(H_STATES.index(image))
    return tuple(action)


def sign_action_on_lines(signs):
    action = {}
    for name, root in ROOT_LINES.items():
        image = canonical_line(tuple(signs[i] * root[i] for i in range(3)))
        require(image in LINE_TO_NAME, "even sign word left the D3 root lines")
        action[name] = LINE_TO_NAME[image]
    return action


INERTIA_PERMS = {}
for word, signs in SIGN_WORDS.items():
    state_action = sign_action_on_states(signs)
    line_action = sign_action_on_lines(signs)
    INERTIA_PERMS[word] = state_action

    for name, edge in ROOT_EDGE.items():
        require(
            ROOT_EDGE[line_action[name]] == act_edge(state_action, edge),
            "root-line and K4-edge actions disagree",
        )

    fixed_marks = {mark for mark in POINTS if state_action[mark] == mark}
    fixed_lines = {name for name in ROOT_LINES if line_action[name] == name}
    if word == "000":
        require(fixed_marks == set(POINTS), "zero inertia must fix every mark")
        require(fixed_lines == set(ROOT_LINES), "zero inertia must fix every line")
    else:
        require(not fixed_marks, "nonzero V4 inertia must fix no mark")
        require(fixed_lines == EXPECTED_FIXED[word], "wrong fixed root-line pair")
        fixed_edges = tuple(ROOT_EDGE[name] for name in sorted(fixed_lines))
        require(are_opposite(*fixed_edges), "fixed root lines must be a weight-two pair")
        for mark in POINTS:
            require(
                marked_weight(mark, fixed_edges[0], fixed_edges[1]) == 2,
                "inertia-fixed pair must have mark-independent weight two",
            )
        cycle_lengths = []
        unseen = set(POINTS)
        while unseen:
            start = min(unseen)
            current = start
            length = 0
            while current in unseen:
                unseen.remove(current)
                length += 1
                current = state_action[current]
            cycle_lengths.append(length)
        require(sorted(cycle_lengths) == [2, 2], "nonzero inertia must pair the marks")


profile_histogram = {}
for edge_a, edge_b in EDGE_PAIRS:
    profile = tuple(sorted(marked_weight(mark, edge_a, edge_b) for mark in POINTS))
    profile_histogram[profile] = profile_histogram.get(profile, 0) + 1
    if are_opposite(edge_a, edge_b):
        require(profile == (2, 2, 2, 2), "wrong symmetrized opposite profile")
    else:
        require(profile == (1, 1, 1, 3), "wrong symmetrized adjacent profile")
require(
    profile_histogram == {(2, 2, 2, 2): 3, (1, 1, 1, 3): 12},
    "wrong complete symmetrized spectrum",
)

row_110 = INERTIA_PERMS["110"]
require(row_110 == (3, 2, 1, 0), "110 must act as (h1 h4)(h2 h3)")

print("MARKED D3 WEIGHT-TORSOR DESCENT AUDIT")
print("marked_colourings=4 each_spectrum=1^9,2^3,3^3")
print("weight3_triangle_reconstructs_mark=YES")
print("S4_stabilizer=S3_order6 A4_stabilizer=C3_order3")
print("V4_action_on_colourings=simply_transitive")
print("generic_marked_colouring_fields=S4:S3_index4,A4:C3_index4")
print("inertia_000=marks_fixed4")
print("inertia_110=mark_cycles=(h1 h4)(h2 h3),fixed_lines=a12,b12")
print("inertia_101=fixed_lines=a13,b13")
print("inertia_011=fixed_lines=a23,b23")
print("nonzero_inertia_mark_fixed_points=0")
print("mod2_frame_nulls=weight1:0,weight2:110/101/011,weight3:0")
print("A2_index3_frames=mod2_unimodular")
print("orbit_symmetrized_opposite_profile=2,2,2,2 count=3")
print("orbit_symmetrized_adjacent_profile=1,1,1,3 count=12")
print("SCOPE=local_inertia_avatar_not_global_section_or_JC2")
print("FAILED CHECKS: NONE")
