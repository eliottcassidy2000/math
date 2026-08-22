#!/usr/bin/env python3
"""Assertion-independent exact companion for THM-3240.

The script checks the mod-13 Heisenberg action on the THM-2763 imbalance
and full exact-address residue quotients.  It uses an explicit THM-2350
owner-pivot normal form only to audit raw representatives; the action proof
itself is coordinate-free after a target axis and a Bezout section are fixed.
"""

from ast import Assert, Constant, parse
from collections import deque
from hashlib import sha256
from itertools import product
from pathlib import Path


P = 13


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "THM-2350": (
        ROOT / "01-canon/theorems/THM-2350-owner-pivot-dual-dipole-normal-form.md",
        "368fe20345351c1c78135a2e0347c98619282ed38e738087b17b14f375371fa0",
    ),
    "THM-2763": (
        ROOT / "01-canon/theorems/THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction.md",
        "1b4a64bc3dd7f307ded417b079ed85b42278f8bfb9c9bdb30dcd9b4f304df2a8",
    ),
    "THM-2779": (
        ROOT / "01-canon/theorems/THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md",
        "21ac7ec9b19b8ed0abe4763e0b7e13ebc1e5eb776168c8e0088540f29daabccb",
    ),
    "THM-3228": (
        ROOT / "01-canon/theorems/THM-3228-four-jet-heisenberg-minimal-faithful-permutation-carrier-gate.md",
        "4dd68b56fbf2b895c5ac4328fa79167992982f06940b8f3c9ebe6738d74b2dfb",
    ),
}


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency_name, (dependency_path, expected_hash) in DEPENDENCIES.items():
    actual_hash = lf_hash(dependency_path)
    require(actual_hash == expected_hash,
            ("dependency hash", dependency_name, actual_hash, expected_hash))

syntax_tree = parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, Assert) for node in __import__("ast").walk(syntax_tree)),
        "assert statements are optimization-sensitive")
require(not any(isinstance(node, Constant) and isinstance(node.value, float)
                for node in __import__("ast").walk(syntax_tree)),
        "floating literals are forbidden")


def h_mul(left, right):
    x, y, central = left
    xp, yp, centralp = right
    return ((x + xp) % P,
            (y + yp) % P,
            (central + centralp - y * xp) % P)


def h_inv(element):
    x, y, central = element
    return (-x % P, -y % P, (-central - x * y) % P)


def h_commutator(left, right):
    return h_mul(h_mul(h_mul(left, right), h_inv(left)), h_inv(right))


def act_pair(element, state):
    """Standard action on the chosen target coordinate and imbalance."""
    x, y, central = element
    target, imbalance = state
    return ((target + x) % P,
            (imbalance + central - y * target) % P)


def act_delta(element, state):
    target, spectator, imbalance = state
    target_out, imbalance_out = act_pair(element, (target, imbalance))
    return (target_out, spectator, imbalance_out)


def act_full(element, state):
    target, spectator, source_harmonic, target_harmonic = state
    imbalance = (source_harmonic - target_harmonic) % P
    target_out, imbalance_out = act_pair(element, (target, imbalance))
    source_out = (target_harmonic + imbalance_out) % P
    return (target_out, spectator, source_out, target_harmonic)


H = tuple(product(range(P), repeat=3))
PAIR_STATES = tuple(product(range(P), repeat=2))
DELTA_STATES = tuple(product(range(P), repeat=3))
FULL_STATES = tuple(product(range(P), repeat=4))

# Exhaust the coefficient form of rho(h)rho(h')=rho(hh').  The target offset,
# imbalance constant, and imbalance coefficient of the input target determine
# the affine map, so this checks all states without a redundant third loop.
composition_pairs = 0
for left in H:
    x, y, central = left
    for right in H:
        xp, yp, centralp = right
        product_element = h_mul(left, right)
        composed_signature = (
            (x + xp) % P,
            (central + centralp - y * xp) % P,
            -(y + yp) % P,
        )
        product_signature = (
            product_element[0],
            product_element[2],
            -product_element[1] % P,
        )
        require(composed_signature == product_signature,
                ("composition sign", left, right,
                 composed_signature, product_signature))
        composition_pairs += 1

permutation_signatures = {
    tuple(act_pair(element, state) for state in PAIR_STATES)
    for element in H
}
require(len(permutation_signatures) == P ** 3,
        ("faithful permutation image", len(permutation_signatures)))

base_orbit = {act_pair(element, (0, 0)) for element in H}
base_stabilizer = {
    element for element in H if act_pair(element, (0, 0)) == (0, 0)
}
expected_stabilizer = {(0, y, 0) for y in range(P)}
require(len(base_orbit) == P * P, ("base orbit", len(base_orbit)))
require(base_stabilizer == expected_stabilizer,
        ("base stabilizer", base_stabilizer))

X = (0, -1 % P, 0)
Y = (1, 0, 0)
Z = (0, 0, 1)
require(h_commutator(X, Y) == Z,
        ("generator commutator", h_commutator(X, Y), Z))
conjugated_stabilizer = {
    h_mul(h_mul(Y, element), h_inv(Y)) for element in base_stabilizer
}
require(conjugated_stabilizer != base_stabilizer,
        "stabilizer was falsely normal")
require(base_stabilizer.intersection({(0, 0, z) for z in range(P)})
        == {(0, 0, 0)}, "stabilizer met the center")
require(all(act_pair((0, 0, z), state) != state
            for z in range(1, P) for state in PAIR_STATES),
        "nontrivial center fixed a carrier point")


def orbit_sizes(states, generators, action):
    unseen = set(states)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        queue = deque([seed])
        while queue:
            current = queue.popleft()
            for generator in generators:
                child = action(generator, current)
                if child not in orbit:
                    orbit.add(child)
                    queue.append(child)
        unseen.difference_update(orbit)
        sizes.append(len(orbit))
    return tuple(sorted(sizes))


delta_orbits = orbit_sizes(DELTA_STATES, (X, Y, Z), act_delta)
full_orbits = orbit_sizes(FULL_STATES, (X, Y, Z), act_full)
require(delta_orbits == (P * P,) * P,
        ("G_delta orbit census", delta_orbits))
require(full_orbits == (P * P,) * (P * P),
        ("G_full orbit census", len(full_orbits), set(full_orbits)))

# Explicit THM-2350 normal form for the fixed THM-2763 word modulo thirteen.
# The first six word entries are one and the owner/target entries are zero.
W = (1, 1, 1, 1, 1, 1, 0, 0, 0)
E0 = (1, 0, 0, 0, 0, 0, 0, 0, 0)
EA = (0, 0, 0, 0, 0, 0, 0, 1, 0)
EB = (0, 0, 0, 0, 0, 0, 0, 0, 1)
L_BASIS = (
    (0, 0, 0, 0, 0, 0, -1 % P, 0, 0),
    (1, -1 % P, 0, 0, 0, 0, 0, 0, 0),
    (1, 0, -1 % P, 0, 0, 0, 0, 0, 0),
    (1, 0, 0, -1 % P, 0, 0, 0, 0, 0),
    (1, 0, 0, 0, -1 % P, 0, 0, -1 % P, 0),
    (1, 0, 0, 0, 0, -1 % P, 0, 0, -1 % P),
)


def vector_add(*vectors):
    return tuple(sum(vector[index] for vector in vectors) % P
                 for index in range(len(vectors[0])))


def vector_scale(scalar, vector):
    return tuple(scalar * value % P for value in vector)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right)) % P


def matrix_rank(rows):
    work = [list(row) for row in rows if any(value % P for value in row)]
    rank = 0
    column = 0
    while rank < len(work) and column < len(work[0]):
        pivot = next((row for row in range(rank, len(work))
                      if work[row][column] % P), None)
        if pivot is None:
            column += 1
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column] % P, P - 2, P)
        work[rank] = [inverse * value % P for value in work[rank]]
        for row in range(len(work)):
            if row != rank and work[row][column] % P:
                factor = work[row][column] % P
                work[row] = [
                    (value - factor * pivot_value) % P
                    for value, pivot_value in zip(work[row], work[rank])
                ]
        rank += 1
        column += 1
    return rank


require(dot(E0, W) == 1 and dot(EA, W) == dot(EB, W) == 0,
        "typed section/target axes")
require(all(dot(row, W) == 0 for row in L_BASIS), "L not in K")
require(matrix_rank(L_BASIS) == 6, "L rank")
require(matrix_rank(L_BASIS + (EA, EB)) == 8, "K=L+targets")

raw_action_checks = 0
for element in H:
    x, y, central = element
    for target, imbalance in PAIR_STATES:
        raw = vector_add(vector_scale(target, EA),
                         vector_scale(-imbalance, E0))
        require((dot(raw, W) + imbalance) % P == 0,
                ("raw input constraint", element, target, imbalance))
        raw_out = vector_add(
            raw,
            vector_scale(x, EA),
            vector_scale(y * target - central, E0),
        )
        target_out, imbalance_out = act_pair(
            element, (target, imbalance))
        require((dot(raw_out, W) + imbalance_out) % P == 0,
                ("raw output constraint", element, target, imbalance))
        sharp_out = vector_add(raw_out, vector_scale(imbalance_out, E0))
        expected_sharp = vector_scale(target_out, EA)
        require(sharp_out == expected_sharp,
                ("raw sharp coordinate", element, target, imbalance,
                 sharp_out, expected_sharp))
        raw_action_checks += 1

# Replacing either chosen lift by an L-vector changes the raw output only by
# an L-vector.  A genuinely different Bezout section is not the same action,
# although transporting both through their coordinate charts conjugates them.
for lam in L_BASIS:
    require(matrix_rank(L_BASIS + (lam,)) == 6, "L lift membership")
    for coefficient in range(P):
        require(matrix_rank(L_BASIS + (vector_scale(coefficient, lam),)) == 6,
                "scaled L lift membership")


def section_change(state):
    target, spectator, imbalance = state
    return ((target - imbalance) % P, spectator, imbalance)


def section_change_inverse(state):
    target, spectator, imbalance = state
    return ((target + imbalance) % P, spectator, imbalance)


def act_alternative_section(element, state):
    """Closed form in the z'=e0-ea Bezout chart, independent of transport."""
    x, y, central = element
    target, spectator, imbalance = state
    old_target = (target + imbalance) % P
    imbalance_out = (imbalance + central - y * old_target) % P
    target_out = (target + x - central + y * old_target) % P
    return (target_out, spectator, imbalance_out)


section_conjugacy_checks = 0
for element in H:
    for state in DELTA_STATES:
        transported = section_change(
            act_delta(element, section_change_inverse(state)))
        require(act_alternative_section(element, state) == transported,
                ("section conjugacy", element, state))
        section_conjugacy_checks += 1
require(act_alternative_section(Z, (0, 0, 0))
        != act_delta(Z, (0, 0, 0)),
        "different Bezout sections falsely gave the same action")

# Forgetting imbalance leaves only the x-translation.  The y-line and center
# are killed, so the old carrier-free THM-2334 target cannot see this clutch.
projected_signatures = {
    tuple((act_delta(element, state)[0], act_delta(element, state)[1])
          for state in DELTA_STATES)
    for element in H
}
require(len(projected_signatures) == P,
        ("old target projection image", len(projected_signatures)))
require(all(act_delta((0, 0, central), state)[:2] == state[:2]
            for central in range(P) for state in DELTA_STATES),
        "center survived the carrier-forgetting projection")

# Exact integral-lift hostile from THM-2788: the physical odometer is not the
# carry-suppressed low shift used by H_13.
def low_shift(state):
    low, high = state
    return ((low + 1) % P, high)


def odometer(state):
    low, high = state
    value = (low + P * high + 1) % (P * P)
    return (value % P, value // P)


def iterate(function, state, count):
    for _ in range(count):
        state = function(state)
    return state


carry_wall_count = sum(odometer(state) != low_shift(state)
                       for state in PAIR_STATES)
require(carry_wall_count == P, ("carry wall", carry_wall_count))
require(all(iterate(low_shift, state, P) == state for state in PAIR_STATES),
        "suppressed shift did not reset")
require(all(iterate(odometer, state, P)
            == (state[0], (state[1] + 1) % P)
            for state in PAIR_STATES),
        "odometer p-th power was not the central shift")

print("EXACT-ADDRESS HEISENBERG CLUTCH EXACT AUDIT")
print("p=13;H_order=2197;minimal_orbit=169;composition_pairs="
      + str(composition_pairs))
print("G_delta_size=2197;orbit_census=13x169;"
      "G_full_size=28561;orbit_census=169x169")
print("faithful_permutation_image=2197;base_stabilizer=(0,y,0):13;"
      "stabilizer_nonnormal=TRUE;center_fixed_points=0")
print("generators=X:(s,d)->(s,d+s);Y:(s,d)->(s+1,d);"
      "Z:(s,d)->(s,d+1);[X,Y]=Z")
print("raw_constraint_checks=" + str(raw_action_checks)
      + ";section_conjugacy_checks=" + str(section_conjugacy_checks))
print("chosen_axis=a;typed_bezout=e0;L_rank=6;K_rank=8;"
      "alternative_section=CONJUGATE_NOT_IDENTICAL")
print("carrier_forgetting_projection_image=13;center_killed=TRUE")
print("integral_lift_hostile=O_vs_Y;carry_wall_states=13;"
      "Y^13=I;O^13=Z")
print("EXTENDED_FACTOR_CURRENT=NOT_CONSTRUCTED;"
      "ALL_91_UNIT_COVARIANCE=NOT_CLAIMED;"
      "OWNER_OR_ROOT_MAP=NOT_CLAIMED")
print("FAILED_CHECKS=NONE")
