#!/usr/bin/env python3
"""Exact finite referee for THM-3067.

The twelve directed edges of the affine plane F_2^2 form a regular A_4
bitorsor.  The right moves are the order-two edge reversal and the order-three
rotation of the three directions at a fixed tail.  This script checks the
tetrahedral quotient of C_2*C_3, its information-loss boundaries, and the
sharp involution-role census in S_4.
"""

from hashlib import sha256
from itertools import combinations, permutations


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


V = tuple(range(4))
DIRECTIONS = (1, 2, 3)

# In the bit encoding 1=(1,0), 2=(0,1), 3=(1,1), this is an order-three
# linear map cycling the three nonzero vectors.
RHO = (0, 2, 3, 1)


def rho_power(x, exponent):
    for _ in range(exponent % 3):
        x = RHO[x]
    return x


FLAGS = tuple((x, d) for x in V for d in DIRECTIONS)
FLAG_INDEX = {flag: index for index, flag in enumerate(FLAGS)}


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in range(len(right)))


def perm_power(perm, exponent):
    result = tuple(range(len(perm)))
    for _ in range(exponent):
        result = compose(perm, result)
    return result


def permutation_of(action):
    return tuple(FLAG_INDEX[action(flag)] for flag in FLAGS)


def closure(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = [identity]
    while queue:
        current = queue.pop()
        for generator in generators:
            nxt = compose(generator, current)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return frozenset(seen)


def edge_reversal(flag):
    x, direction = flag
    return (x ^ direction, direction)


def direction_rotation(flag):
    x, direction = flag
    return (x, RHO[direction])


S = permutation_of(edge_reversal)
R = permutation_of(direction_rotation)
IDENTITY_12 = tuple(range(12))
SR = compose(R, S)  # first S, then R: right multiplication by S R

require(perm_power(S, 2) == IDENTITY_12, "S is not an involution")
require(perm_power(R, 3) == IDENTITY_12, "R is not order three")
require(perm_power(SR, 3) == IDENTITY_12, "the tetrahedral relation failed")

RIGHT_GROUP = closure((S, R))
require(len(RIGHT_GROUP) == 12, ("right group order", len(RIGHT_GROUP)))
base_index = FLAG_INDEX[(0, 1)]
require({g[base_index] for g in RIGHT_GROUP} == set(range(12)), "orbit not full")
require(sum(g[base_index] == base_index for g in RIGHT_GROUP) == 1,
        "right action not regular")


# The abstract affine group A_4=F_2^2 semidirect C_3.
AFFINE = tuple((translation, power) for translation in V for power in range(3))


def affine_product(left, right):
    a, k = left
    b, ell = right
    return (a ^ rho_power(b, k), (k + ell) % 3)


def flag_of_affine(element):
    translation, power = element
    return (translation, rho_power(1, power))


S_AFFINE = (1, 0)
R_AFFINE = (0, 1)
for element in AFFINE:
    require(flag_of_affine(affine_product(element, S_AFFINE)) ==
            edge_reversal(flag_of_affine(element)), ("right S", element))
    require(flag_of_affine(affine_product(element, R_AFFINE)) ==
            direction_rotation(flag_of_affine(element)), ("right R", element))


def left_affine(element, flag):
    translation, power = element
    x, direction = flag
    return (translation ^ rho_power(x, power), rho_power(direction, power))


left_permutations = {permutation_of(lambda flag, g=g: left_affine(g, flag))
                     for g in AFFINE}
require(len(left_permutations) == 12, "left affine action is not faithful")
for left in left_permutations:
    require(compose(left, S) == compose(S, left), "left/right S do not commute")
    require(compose(left, R) == compose(R, left), "left/right R do not commute")


# The direction quotient has image C_3 and kernel V_4.
direction_actions = {}
for group_element in RIGHT_GROUP:
    image = []
    for direction in DIRECTIONS:
        images = {FLAGS[group_element[FLAG_INDEX[(x, direction)]]][1] for x in V}
        require(len(images) == 1, ("direction depends on point", direction, images))
        image.append(next(iter(images)))
    direction_actions[group_element] = tuple(image)
direction_image = set(direction_actions.values())
direction_identity = DIRECTIONS
direction_kernel = {g for g, image in direction_actions.items()
                    if image == direction_identity}
require(len(direction_image) == 3, ("direction image", len(direction_image)))
require(len(direction_kernel) == 4, ("origin kernel", len(direction_kernel)))


# Quotienting by S gives the six unoriented edges of K_4, but R does not
# descend: two representatives of one S-orbit can rotate to different edges.
def edge_of(flag):
    x, direction = flag
    return tuple(sorted((x, x ^ direction)))


edges = {edge_of(flag) for flag in FLAGS}
require(edges == set(combinations(V, 2)), ("edge census", edges))
require(all(edge_of(flag) == edge_of(edge_reversal(flag)) for flag in FLAGS),
        "S does not reverse one edge")
hostile_flag = (0, 1)
hostile_reverse = edge_reversal(hostile_flag)
rotated_edges = (edge_of(direction_rotation(hostile_flag)),
                 edge_of(direction_rotation(hostile_reverse)))
require(rotated_edges[0] != rotated_edges[1], ("R unexpectedly descends", rotated_edges))

# Every geometric edge stabilizer in the left A_4 action contains the
# nonzero translation by that edge's direction, which swaps the endpoints.
for x, y in edges:
    direction = x ^ y
    require(left_affine((direction, 0), (x, direction)) == (y, direction),
            ("edge orientation unexpectedly invariant", x, y))


# An odd linear reflection fixes S but reverses the C_3 orientation.
REFLECTION = (0, 2, 1, 3)  # swap the two basis vectors


def reflected(flag):
    x, direction = flag
    return (REFLECTION[x], REFLECTION[direction])


J = permutation_of(reflected)
require(perm_power(J, 2) == IDENTITY_12, "J is not an involution")
require(compose(J, compose(S, J)) == S, "J does not fix S")
require(compose(J, compose(R, J)) == perm_power(R, 2), "J does not invert R")


# Fixed-R census of every involution in S_4.  This shows that the modular
# involution has three genuinely different jobs: inner reflection, cross
# reflection, or V_4 translation.
def point_compose(left, right):
    return tuple(left[right[index]] for index in V)


def point_closure(generators):
    identity = V
    seen = {identity}
    queue = [identity]
    while queue:
        current = queue.pop()
        for generator in generators:
            nxt = point_compose(generator, current)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


def cycle_lengths(perm):
    seen = set()
    lengths = []
    for start in V:
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = perm[current]
        if length > 1:
            lengths.append(length)
    return tuple(sorted(lengths))


def permutation_sign(perm):
    inversions = sum(perm[i] > perm[j]
                     for i in range(len(perm))
                     for j in range(i + 1, len(perm)))
    return -1 if inversions % 2 else 1


R_POINT = RHO
fixed_point = 0
census = {"identity": [], "inner_transposition": [],
          "cross_transposition": [], "double_transposition": []}
for perm in permutations(V):
    perm = tuple(perm)
    if point_compose(perm, perm) != V:
        continue
    cycles = cycle_lengths(perm)
    if not cycles:
        kind = "identity"
    elif cycles == (2,):
        kind = ("cross_transposition" if perm[fixed_point] != fixed_point
                else "inner_transposition")
    elif cycles == (2, 2):
        kind = "double_transposition"
    else:
        raise RuntimeError(("unexpected involution", perm, cycles))
    census[kind].append(len(point_closure((perm, R_POINT))))

require(census == {
    "identity": [3],
    "inner_transposition": [6, 6, 6],
    "cross_transposition": [24, 24, 24],
    "double_transposition": [12, 12, 12],
}, ("involution census", census))

# The chosen modular S is a nonzero affine translation and therefore a double
# transposition on the four points.
translation_one = tuple(x ^ 1 for x in V)
require(cycle_lengths(translation_one) == (2, 2), "chosen S is not V4-type")
require(len(point_closure((translation_one, R_POINT))) == 12,
        "chosen S,R do not generate A4")

affine_point_permutations = {
    tuple(translation ^ rho_power(x, power) for x in V)
    for translation, power in AFFINE
}
require(len(affine_point_permutations) == 12, "affine point action is not faithful")
require(all(permutation_sign(perm) == 1 for perm in affine_point_permutations),
        "the affine subgroup is not A4")
FULL_AFFINE = closure(tuple(left_permutations) + (J,))
require(len(FULL_AFFINE) == 24, ("full affine order", len(FULL_AFFINE)))


semantic = repr((
    FLAGS,
    sorted(RIGHT_GROUP),
    sorted(direction_image),
    sorted(edges),
    rotated_edges,
    census,
    sorted(FULL_AFFINE),
)).encode("ascii")
semantic_sha256 = sha256(semantic).hexdigest()

print("THM3067_TETRAHEDRAL_MODULAR_FLAG_BITORSOR")
print("flags=12;unoriented_edges=6;left_affine=12;right_modular_quotient=12")
print("relations=S^2=R^3=(SR)^3=1;right_action_regular=1;group=A4")
print("direction_image=3;origin_kernel=4;quotient=C3")
print("edge_reversal_orbits=6;R_descends_after_forgetting_orientation=0")
print(f"no_descent_hostile={hostile_flag}/{hostile_reverse}->{rotated_edges}")
print("odd_reflection=J^2=1,JSJ=S,JRJ=R^-1;full_affine_group=S4")
print("involution_census=id:1x3,inner_transposition:3x6,cross_transposition:3x24,double_transposition:3x12")
print("invariant_unoriented_edge_orientation=0;stabilizer_translation_swaps_endpoints=1")
print(f"semantic_sha256={semantic_sha256}")
print("scope=finite_chosen_orientation_quotient;no_physical_LRC_or_JC_intertwiner")
