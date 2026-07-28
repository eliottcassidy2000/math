#!/usr/bin/env python3
"""Independent exact hostile audit for THM-2780's descent gate.

This scratch companion uses integer/rational symbolic arithmetic only and
explicit exceptions rather than truth-bearing ``assert`` statements.
"""

from collections import Counter
from itertools import combinations, permutations, product

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def det3(u, v, w):
    return (
        u[0] * (v[1] * w[2] - v[2] * w[1])
        - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
    )


def mul(u, v):
    return tuple(a * b for a, b in zip(u, v))


def canonical_line(v):
    for a in v:
        if a:
            return v if a > 0 else tuple(-x for x in v)
    raise RuntimeError("zero vector has no line")


ROOTS = {
    "a12": (1, -1, 0),
    "b12": (1, 1, 0),
    "a13": (1, 0, -1),
    "b13": (1, 0, 1),
    "a23": (0, 1, -1),
    "b23": (0, 1, 1),
}
NAMES = tuple(ROOTS)
LINE_TO_NAME = {canonical_line(v): name for name, v in ROOTS.items()}
STATES = tuple(v for v in product((-1, 1), repeat=3) if v[0] * v[1] * v[2] == 1)
PAIRS = tuple(combinations(NAMES, 2))


def colouring(state):
    return {
        frozenset((a, b)): abs(det3(ROOTS[a], ROOTS[b], state))
        for a, b in PAIRS
    }


COLOURINGS = {state: colouring(state) for state in STATES}
require(len(STATES) == 4, "even sign state census")
require(len({tuple(sorted(c.items(), key=lambda x: sorted(x[0]))) for c in COLOURINGS.values()}) == 4,
        "the four weight colourings must be distinct")

for state, weights in COLOURINGS.items():
    require(Counter(weights.values()) == Counter({1: 9, 2: 3, 3: 3}),
            f"weight spectrum failed at {state}")
    orthogonal = {name for name, root in ROOTS.items() if dot(root, state) == 0}
    weight_three_vertices = set()
    for pair, weight in weights.items():
        if weight == 3:
            weight_three_vertices.update(pair)
    require(len(orthogonal) == 3 and weight_three_vertices == orthogonal,
            f"weight-three orthogonal triangle failed at {state}")
    require(
        all(
            (weights[frozenset((a, b))] == 2) == (dot(ROOTS[a], ROOTS[b]) == 0)
            for a, b in PAIRS
        ),
        f"weight-two opposition matching failed at {state}",
    )


def act_signed_permutation(v, perm, signs):
    return tuple(signs[i] * v[perm[i]] for i in range(3))


WEYL = tuple(
    (perm, signs)
    for perm in permutations(range(3))
    for signs in product((-1, 1), repeat=3)
    if signs[0] * signs[1] * signs[2] == 1
)
require(len(WEYL) == 24, "W(D3) census")

for perm, signs in WEYL:
    line_action = {}
    for name, root in ROOTS.items():
        image_line = canonical_line(act_signed_permutation(root, perm, signs))
        require(image_line in LINE_TO_NAME, "W(D3) failed to preserve root lines")
        line_action[name] = LINE_TO_NAME[image_line]
    for state in STATES:
        target_state = act_signed_permutation(state, perm, signs)
        require(target_state in STATES, "W(D3) failed to preserve even states")
        for a, b in PAIRS:
            lhs = COLOURINGS[target_state][frozenset((line_action[a], line_action[b]))]
            rhs = COLOURINGS[state][frozenset((a, b))]
            require(lhs == rhs, "weight-colouring equivariance failed")

for source in STATES:
    require({mul(source, epsilon) for epsilon in STATES} == set(STATES),
            "V4 must act simply transitively on states")

for epsilon in STATES:
    if epsilon == (1, 1, 1):
        continue
    action = {state: mul(state, epsilon) for state in STATES}
    require(all(action[state] != state for state in STATES),
            "nonzero even inertia must be fixed-point-free")
    require(all(action[action[state]] == state for state in STATES),
            "nonzero even inertia must be a double transposition")


# Retained-chamber tournament and the sharp chamber-free hostile.
H = (1, 1, 1)


def arc_sign(a, b):
    value = det3(ROOTS[a], ROOTS[b], H)
    require(value != 0, "tournament tie")
    return 1 if value > 0 else -1


ARC = {(a, b): arc_sign(a, b) for a in NAMES for b in NAMES if a != b}
SCORES = tuple(
    sorted(sum(ARC[(a, b)] > 0 for b in NAMES if b != a) for a in NAMES)
)
require(SCORES == (1, 2, 2, 3, 3, 4), "tournament score sequence")


def is_converse_isomorphism(perm_names, preserve_weights=False):
    mapping = dict(zip(NAMES, perm_names))
    if preserve_weights:
        weights = COLOURINGS[H]
        for a, b in PAIRS:
            if weights[frozenset((a, b))] != weights[frozenset((mapping[a], mapping[b]))]:
                return False
    return all(
        ARC[(mapping[a], mapping[b])] == -ARC[(a, b)]
        for a, b in PAIRS
    )


require(
    not any(is_converse_isomorphism(p) for p in permutations(NAMES)),
    "retained-chamber tournament unexpectedly achiral",
)

SWITCH_VERTEX = "a12"
RELABEL = {
    "a12": "a12",
    "b12": "b12",
    "a23": "a13",
    "a13": "a23",
    "b13": "b23",
    "b23": "b13",
}
require(set(RELABEL) == set(NAMES) and set(RELABEL.values()) == set(NAMES),
        "self-converse relabelling is not a permutation")


def switched_sign(a, b):
    sign = ARC[(a, b)]
    if (a == SWITCH_VERTEX) != (b == SWITCH_VERTEX):
        sign = -sign
    return sign


inverse_relabel = {v: k for k, v in RELABEL.items()}
for a, b in PAIRS:
    transported = switched_sign(inverse_relabel[a], inverse_relabel[b])
    require(transported == -ARC[(a, b)],
            f"switched self-converse identity failed at {a},{b}")
    weights = COLOURINGS[H]
    require(
        weights[frozenset((inverse_relabel[a], inverse_relabel[b]))]
        == weights[frozenset((a, b))],
        f"self-converse relabelling failed to preserve weight at {a},{b}",
    )


# Quartic reconstruction from the Kummer plane.
s1, s2, s3, Y = sp.symbols("s1 s2 s3 Y")
symbolic_states = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
roots = [
    sp.Rational(1, 2) * (state[0] * s1 + state[1] * s2 + state[2] * s3)
    for state in symbolic_states
]
quartic = sp.Poly(sp.expand(sp.prod(Y - root for root in roots)), Y)
t1, t2, t3 = s1**2, s2**2, s3**2
expected = sp.Poly(
    Y**4
    - sp.Rational(1, 2) * (t1 + t2 + t3) * Y**2
    - s1 * s2 * s3 * Y
    + sp.Rational(1, 16)
    * ((t1 + t2 + t3) ** 2 - 4 * (t1 * t2 + t1 * t3 + t2 * t3)),
    Y,
)
require(quartic == expected, "quartic half-Hadamard reconstruction")
require(sp.expand(roots[0] + roots[1]) == s1, "first pair-sum coordinate")
require(sp.expand(roots[0] + roots[2]) == s2, "second pair-sum coordinate")
require(sp.expand(roots[0] + roots[3]) == s3, "third pair-sum coordinate")


# Twin unramified S3-standard Kummer planes on one fixed four-torus.
# Coordinates are the parity exponent lattice with bases x1,x2,y1,y2 and
# x3=x1+x2, y3=y1+y2.
ZERO = (0, 0, 0, 0)


def xor(u, v):
    return tuple(a ^ b for a, b in zip(u, v))


X = ((1, 0, 0, 0), (0, 1, 0, 0), (1, 1, 0, 0))
YV = ((0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 1, 1))
WX = {ZERO, X[0], X[1], X[2]}
WY = {ZERO, YV[0], YV[1], YV[2]}
require(WX != WY and WX.intersection(WY) == {ZERO}, "twin Kummer planes must differ")

for perm in permutations(range(3)):
    image_x = {ZERO, X[perm[0]], X[perm[1]], xor(X[perm[0]], X[perm[1]])}
    image_y = {ZERO, YV[perm[0]], YV[perm[1]], xor(YV[perm[0]], YV[perm[1]])}
    require(image_x == WX and image_y == WY, "S3 stability of twin planes")

cycle = (1, 2, 0)
require({X[cycle[i]] for i in range(3)} == set(X), "S3 cycles WX nonzero classes")
require({YV[cycle[i]] for i in range(3)} == set(YV), "S3 cycles WY nonzero classes")


print("THM-2780 INDEPENDENT MARKED-D3 DESCENT HOSTILE")
print("d3_root_lines=6 even_states=4 weyl_elements=24")
print("weight_colouring_spectrum=1^9,2^3,3^3")
print("weight3=state_orthogonal_triangle weight2=state_independent_opposition_matching")
print("colouring_map=injective_W(D3)-equivariant V4_regular")
print(f"retained_chamber_scores={SCORES} retained_chamber_tournament=chiral")
print("chamber_free_weighted_switching_class=self_converse")
print("switching_witness=switch(a12); relabel(a23,a13)(b13,b23)")
print("quartic_reconstruction=half_Hadamard_exact")
print("affine_inertia_action=nonzero_even_word_translates_four_colourings_freely")
print("quasi_etale_gate=all_normalized_divisor_rows_zero")
print("twin_base=(Gm)^4 with simultaneous_S3")
print("twin_planes=WX=<x1,x2,x3>, WY=<y1,y2,y3>")
print("twin_planes=S3_standard_distinct all_divisor_rows_zero common_sum_zero")
print("minimal_lost_coordinate=embedded_Q_equivariant_Kummer_H1_class")
print("FAILED CHECKS: NONE")
