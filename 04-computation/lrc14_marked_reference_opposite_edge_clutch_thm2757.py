#!/usr/bin/env python3
"""Exact companion for THM-2757.

The calculation is finite and integral.  It verifies the marked K4
star/opposite transform, its Hadamard identities, a formal two-channel clutch
factorization, the complete F_13 charge criterion, and the sharp S4/A4
equal-corolla collision.  The (12,2) specialization is deliberately only an
algebraic control after MISTAKE-313: no physical wing-diagonal operator is
claimed.  No optimized-away assertions are used.
"""

from itertools import product


P = 13
ALPHA = 12
BETA = 2
INV2 = pow(2, -1, P)
INV3 = pow(3, -1, P)
OMEGA = 3


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in transpose(b)] for row in a]


def matvec(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def modv(v):
    return tuple(x % P for x in v)


def edge_table(v):
    v0, v1, v2, v3 = v
    star = (v0 + v1, v0 + v2, v0 + v3)
    opposite = (v2 + v3, v1 + v3, v1 + v2)
    return modv(star), modv(opposite)


def clutch(star, opposite):
    return tuple((ALPHA * x + BETA * y) % P for x, y in zip(star, opposite))


def q3(v):
    mean = sum(v) * INV3 % P
    return tuple((x - mean) % P for x in v)


def fourier3(v, k):
    return INV3 * sum(v[i] * pow(OMEGA, (-k * i) % 3, P) for i in range(3)) % P


def permute_values(v, permutation):
    # All permutations used below are involutions, so inverse=itself.
    return tuple(v[permutation[i]] for i in range(4))


def matching_action(permutation):
    matchings = (
        (frozenset((0, 1)), frozenset((2, 3))),
        (frozenset((0, 2)), frozenset((1, 3))),
        (frozenset((0, 3)), frozenset((1, 2))),
    )
    canonical = {frozenset(m): i for i, m in enumerate(matchings)}
    image = []
    for matching in matchings:
        moved = frozenset(
            frozenset((permutation[min(edge)], permutation[max(edge)])) for edge in matching
        )
        image.append(canonical[moved])
    return tuple(image)


D = [
    [1, 1, -1, -1],
    [1, -1, 1, -1],
    [1, -1, -1, 1],
]
I3 = [[int(i == j) for j in range(3)] for i in range(3)]
I4 = [[int(i == j) for j in range(4)] for i in range(4)]
J4 = [[1 for _ in range(4)] for _ in range(4)]

require(matmul(D, transpose(D)) == [[4 * x for x in row] for row in I3], "D D^T")
require(
    matmul(transpose(D), D)
    == [[4 * I4[i][j] - J4[i][j] for j in range(4)] for i in range(4)],
    "D^T D",
)

charged = 0
flat = 0
for v in product(range(P), repeat=4):
    star, opposite = edge_table(v)
    total = sum(v) % P
    require(all((x + y) % P == total for x, y in zip(star, opposite)), "matching total")
    difference = modv(x - y for x, y in zip(star, opposite))
    require(difference == modv(matvec(D, v)), "opposite-edge transform")
    target = clutch(star, opposite)
    rhs = tuple(
        ((ALPHA + BETA) * total + (ALPHA - BETA) * d) % P
        for d in difference
    )
    require(tuple(2 * x % P for x in target) == rhs, "integral clutch factorization")
    is_flat = len(set(target)) == 1
    require(is_flat == (v[1] == v[2] == v[3]), "complete charge criterion")
    flat += int(is_flat)
    charged += int(not is_flat)

require(flat == P * P, "flat census")
require(charged == P**4 - P**2, "charged census")

s01 = (1, 0, 2, 3)
d01_23 = (1, 0, 3, 2)
for x, a in product(range(P), repeat=2):
    base = (x, a, a, a)
    star0, opposite0 = edge_table(base)
    require(len(set(clutch(star0, opposite0))) == 1, "unmoved invariant corolla")

    moved_s = permute_values(base, s01)
    moved_d = permute_values(base, d01_23)
    require(moved_s == moved_d, "equal-corolla S4/A4 collision")
    star, opposite = edge_table(moved_s)
    delta = (x - a) % P
    require(modv(u - w for u, w in zip(star, opposite)) == (delta, -delta % P, -delta % P), "swapped D-profile")
    target = clutch(star, opposite)
    require(q3(target) == modv((11 * delta, delta, delta)), "charged coordinate profile")
    require(fourier3(target, 1) == 12 * delta % P, "chi profile")
    require(fourier3(target, 2) == 12 * delta % P, "chi2 profile")
    rectangles = (
        (star[0] + opposite[1] - star[1] - opposite[0]) % P,
        (star[0] + opposite[2] - star[2] - opposite[0]) % P,
        (star[1] + opposite[2] - star[2] - opposite[1]) % P,
    )
    require(rectangles == (2 * delta % P, 2 * delta % P, 0), "mixed rectangles")
    require((target[0] - target[1]) % P == 5 * rectangles[0] % P, "rectangle clutch")

require(matching_action(s01) != matching_action(d01_23), "matching quotient separation")
require(matching_action(d01_23) == (0, 1, 2), "V4 matching kernel")

print("THM2757 exact marked-reference opposite-edge clutch transgression")
print("field=F13 formal_alpha=12 formal_beta=2 omega=3")
print("DDt=4I3")
print("DtD=4I4-J4")
print("matching_total_constant=1")
print("integral_factorization=2T-(alpha+beta)S-(alpha-beta)D=0")
print(f"flat_vertex_tables={flat}")
print(f"charged_vertex_tables={charged}")
print("charge_iff_three_arms_not_equal=1")
print("marked_swap_D_profile=delta*(1,-1,-1)")
print("marked_swap_Q3_target=delta*(11,1,1)")
print("marked_swap_chi1=12*delta")
print("marked_swap_chi2=12*delta")
print("mixed_rectangles=delta*(2,2,0)")
print("target_difference=5*mixed_rectangle")
print(f"s01_matching_action={matching_action(s01)}")
print(f"d01_23_matching_action={matching_action(d01_23)}")
print("equal_corolla_edge_table_collision=1")
print("formal_12_2_specialization_only=1")
print("physical_wing_diagonal_from_THM2751=0")
print("physical_LRC_four_point_carrier=0")
print("lrc14_closed=0")
