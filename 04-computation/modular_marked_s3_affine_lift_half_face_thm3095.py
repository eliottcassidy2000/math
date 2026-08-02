#!/usr/bin/env python3
"""Exact marked S3 affine-lift and co-occurrence-frame census for THM-3095."""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STORED = ROOT / "05-knowledge/results/modular_marked_s3_affine_lift_half_face_thm3095.out"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def identity(degree):
    return tuple(range(degree))


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def power(permutation, exponent):
    value = identity(len(permutation))
    base = permutation
    while exponent:
        if exponent & 1:
            value = compose(base, value)
        base = compose(base, base)
        exponent >>= 1
    return value


def permutation_order(permutation):
    for exponent in range(1, 25):
        if power(permutation, exponent) == identity(len(permutation)):
            return exponent
    raise RuntimeError(("order bound", permutation))


def generated_subgroup(generators):
    one = identity(len(generators[0]))
    seen = {one}
    queue = deque([one])
    while queue:
        value = queue.popleft()
        for generator in generators:
            target = compose(generator, value)
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return frozenset(seen)


def transposition(degree, first, second):
    value = list(range(degree))
    value[first], value[second] = value[second], value[first]
    return tuple(value)


def conjugate(element, gauge):
    return compose(compose(gauge, element), inverse(gauge))


def fixed_points(permutation):
    return tuple(index for index, target in enumerate(permutation) if index == target)


def cycle_string(permutation):
    unseen = set(range(len(permutation)))
    cycles = []
    while unseen:
        seed = min(unseen)
        cycle = []
        value = seed
        while value in unseen:
            unseen.remove(value)
            cycle.append(value)
            value = permutation[value]
        if len(cycle) > 1:
            cycles.append("(" + "".join(str(entry) for entry in cycle) + ")")
    return "".join(cycles) or "1"


S4 = tuple(permutations(range(4)))
ONE4 = identity(4)

MATCHINGS = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def act_on_matching(permutation, matching):
    return frozenset(
        frozenset(permutation[vertex] for vertex in edge) for edge in matching
    )


def matching_quotient(permutation):
    return tuple(
        MATCHINGS.index(act_on_matching(permutation, matching))
        for matching in MATCHINGS
    )


def matching_translation(matching):
    value = list(range(4))
    for edge in matching:
        first, second = sorted(edge)
        value[first] = second
        value[second] = first
    return tuple(value)


DIRECTIONS = tuple(matching_translation(matching) for matching in MATCHINGS)
V4 = frozenset((ONE4, *DIRECTIONS))


def fixed_direction(quotient_involution):
    fixed = fixed_points(quotient_involution)
    require(len(fixed) == 1, ("quotient fixed direction", quotient_involution))
    return DIRECTIONS[fixed[0]]


def local_cyclic_orientation(three_cycle):
    owner = fixed_points(three_cycle)
    require(len(owner) == 1, ("three-cycle owner", three_cycle))
    support = sorted(set(range(4)) - set(owner))
    seed = min(support)
    orientation = (seed, three_cycle[seed], power(three_cycle, 2)[seed])
    require(set(orientation) == set(support), ("orientation", three_cycle))
    return orientation


def direction_orientation(three_cycle):
    quotient_cycle = matching_quotient(three_cycle)
    require(
        permutation_order(quotient_cycle) == 3 and not fixed_points(quotient_cycle),
        ("direction orientation", three_cycle, quotient_cycle),
    )
    return (0, quotient_cycle[0], power(quotient_cycle, 2)[0])


def local_orientation_in_directions(three_cycle):
    owner = fixed_points(three_cycle)[0]
    orientation = []
    for point in local_cyclic_orientation(three_cycle):
        directions = [
            index
            for index, direction in enumerate(DIRECTIONS)
            if direction[owner] == point
        ]
        require(len(directions) == 1, ("point-direction identification", owner, point))
        orientation.append(directions[0])
    return tuple(orientation)


def cyclic_rotations(triple):
    require(len(triple) == 3, ("cyclic triple", triple))
    return {triple[index:] + triple[:index] for index in range(3)}


QUOTIENT = frozenset(matching_quotient(permutation) for permutation in S4)
KERNEL = frozenset(
    permutation for permutation in S4 if matching_quotient(permutation) == identity(3)
)
require(len(S4) == 24, "S4 order")
require(len(QUOTIENT) == 6, "matching quotient order")
require(KERNEL == V4, ("matching kernel", KERNEL))

QUOTIENT_PAIRS = tuple(
    (binary, ternary)
    for binary in QUOTIENT
    for ternary in QUOTIENT
    if permutation_order(binary) == 2
    and permutation_order(ternary) == 3
    and len(generated_subgroup((binary, ternary))) == 6
)
require(len(QUOTIENT_PAIRS) == 6, "marked quotient pairs")

LIFTS = []
for binary in S4:
    if permutation_order(binary) != 2 or permutation_order(matching_quotient(binary)) != 2:
        continue
    for ternary in S4:
        if permutation_order(ternary) != 3 or permutation_order(matching_quotient(ternary)) != 3:
            continue
        quotient_pair = (matching_quotient(binary), matching_quotient(ternary))
        require(quotient_pair in QUOTIENT_PAIRS, ("quotient pair", quotient_pair))
        word = compose(binary, ternary)
        quotient_word = compose(*quotient_pair)
        direction = fixed_direction(quotient_word)
        group_order = len(generated_subgroup((binary, ternary)))
        require(group_order in (6, 24), ("lift group", group_order))
        branch = 0 if group_order == 6 else 1
        owner = fixed_points(ternary)
        require(len(owner) == 1, ("owner", ternary))
        frame = (
            owner[0],
            direction,
            direction_orientation(ternary),
            branch,
        )
        require(
            power(word, 2) == (ONE4 if branch == 0 else direction),
            ("half-face formula", binary, ternary, word, direction, branch),
        )
        require(
            permutation_order(word) == (2 if branch == 0 else 4),
            ("word order", word, branch),
        )
        LIFTS.append((binary, ternary, quotient_pair, frame))

require(len(LIFTS) == 48, ("lift count", len(LIFTS)))
require(Counter(row[3][3] for row in LIFTS) == {0: 24, 1: 24}, "branch count")
require(len({row[3] for row in LIFTS}) == 48, "quadruple bijection")

# The explicit inverse to the frame map.  If c fixes o, d(o)=b, and
# c cycles a=c^{-1}(b), b, e=c(b), then the split and full binary lifts are
# (b e) and (o a), respectively.
for binary, ternary, quotient_pair, frame in LIFTS:
    owner, direction, orientation, branch = frame
    partner = direction[owner]
    previous = inverse(ternary)[partner]
    following = ternary[partner]
    split_binary = transposition(4, partner, following)
    full_binary = transposition(4, owner, previous)
    reconstructed = split_binary if branch == 0 else full_binary
    require(binary == reconstructed, ("explicit reconstruction", frame, binary, reconstructed))
    require(
        direction_orientation(ternary)
        in cyclic_rotations(local_orientation_in_directions(ternary)),
        ("local/direction orientation reconstruction", ternary),
    )
    require(direction_orientation(ternary) == orientation, "direction orientation reconstruction")
    require(matching_quotient(reconstructed) == quotient_pair[0], "binary quotient")
    require(
        len(generated_subgroup((split_binary, ternary))) == 6
        and len(generated_subgroup((full_binary, ternary))) == 24,
        ("split/full group", frame),
    )

FULL = tuple(row for row in LIFTS if row[3][3] == 1)
SPLIT = tuple(row for row in LIFTS if row[3][3] == 0)
require(len(FULL) == len(SPLIT) == 24, "full/split census")

FULL_TRIPLES = {(owner, direction, orientation) for _, _, _, (owner, direction, orientation, _) in FULL}
FLAGS = {(owner, direction) for owner, direction, _ in FULL_TRIPLES}
require(len(FULL_TRIPLES) == 24, "full co-occurrence triples")
require(len(FLAGS) == 12, "tetrahedral flags")
require(
    Counter((owner, direction) for owner, direction, _ in FULL_TRIPLES)
    == {flag: 2 for flag in FLAGS},
    "two orientations per flag",
)
require(Counter(owner for owner, _, _ in FULL_TRIPLES) == {index: 6 for index in range(4)}, "owner census")
require(Counter(direction for _, direction, _ in FULL_TRIPLES) == {direction: 8 for direction in DIRECTIONS}, "direction census")

# For full S4 markings, the matching quotient forgets exactly the owner.
FULL_QUOTIENT_FIBRES = defaultdict(list)
ALL_QUOTIENT_FIBRES = defaultdict(list)
for row in LIFTS:
    ALL_QUOTIENT_FIBRES[row[2]].append(row)
for row in FULL:
    FULL_QUOTIENT_FIBRES[row[2]].append(row)
require(set(FULL_QUOTIENT_FIBRES) == set(QUOTIENT_PAIRS), "full quotient keys")
require(set(ALL_QUOTIENT_FIBRES) == set(QUOTIENT_PAIRS), "all quotient keys")
for quotient_pair, fibre in FULL_QUOTIENT_FIBRES.items():
    require(len(fibre) == 4, ("full quotient fibre", quotient_pair))
    require({row[3][0] for row in fibre} == set(range(4)), ("owner fibre", quotient_pair))
    require(len({row[3][1] for row in fibre}) == 1, ("direction fibre", quotient_pair))
    require(len({row[3][2] for row in fibre}) == 1, ("orientation fibre", quotient_pair))
for quotient_pair, fibre in ALL_QUOTIENT_FIBRES.items():
    require(len(fibre) == 8, ("all quotient fibre", quotient_pair))
    require(Counter(row[3][0] for row in fibre) == {index: 2 for index in range(4)}, "owner/branch fibre")
    require(Counter(row[3][3] for row in fibre) == {0: 4, 1: 4}, "branch fibre")

# The quotient marking is equivalently one direction and one cyclic orientation.
DIRECTION_ORIENTATION_TO_QUOTIENT = {}
for _, _, quotient_pair, (owner, direction, orientation, branch) in LIFTS:
    key = (direction, quotient_pair[1])
    if key in DIRECTION_ORIENTATION_TO_QUOTIENT:
        require(DIRECTION_ORIENTATION_TO_QUOTIENT[key] == quotient_pair, ("quotient reconstruction", key))
    DIRECTION_ORIENTATION_TO_QUOTIENT[key] = quotient_pair
require(len(DIRECTION_ORIENTATION_TO_QUOTIENT) == 6, "direction-orientation quotient")

# Simultaneous S4 conjugation transports every coordinate equivariantly.
FRAME_LOOKUP = {(binary, ternary): frame for binary, ternary, _, frame in LIFTS}
for gauge in S4:
    for binary, ternary, _, frame in LIFTS:
        transported_pair = (conjugate(binary, gauge), conjugate(ternary, gauge))
        require(transported_pair in FRAME_LOOKUP, ("transported pair", transported_pair))
        transported = FRAME_LOOKUP[transported_pair]
        owner, direction, orientation, branch = frame
        transported_orientation = tuple(
            DIRECTIONS.index(conjugate(DIRECTIONS[index], gauge))
            for index in orientation
        )
        require(transported[0] == gauge[owner], ("owner equivariance", frame, gauge))
        require(transported[1] == conjugate(direction, gauge), ("direction equivariance", frame, gauge))
        require(
            transported[2] in cyclic_rotations(transported_orientation),
            ("orientation equivariance", frame, gauge),
        )
        require(transported[3] == branch, ("branch equivariance", frame, gauge))


def elementary_quartic(roots):
    require(sum(roots) == 0, ("depressed roots", roots))
    polynomial = [Fraction(1)]
    for root in roots:
        target = [Fraction(0)] * (len(polynomial) + 1)
        for degree, coefficient in enumerate(polynomial):
            target[degree] -= coefficient * root
            target[degree + 1] += coefficient
        polynomial = target
    require(polynomial[3] == 0 and polynomial[4] == 1, ("depressed polynomial", polynomial))
    return polynomial[2], polynomial[1], polynomial[0]


def resolvent_value(p_value, q_value, r_value, u_value):
    return (
        u_value**3
        + 2 * p_value * u_value**2
        + (p_value**2 - 4 * r_value) * u_value
        - q_value**2
    )


def multiply_polynomials(left, right):
    target = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            target[left_degree + right_degree] += left_coefficient * right_coefficient
    return target


ROOT_CONTROLS = (
    (Fraction(0), Fraction(1), Fraction(2), Fraction(-3)),
    (Fraction(1), Fraction(2), Fraction(4), Fraction(-7)),
    (Fraction(-5), Fraction(-1), Fraction(2), Fraction(4)),
    (Fraction(-8), Fraction(1), Fraction(3), Fraction(4)),
)

for roots in ROOT_CONTROLS:
    p_value, q_value, r_value = elementary_quartic(roots)
    for direction in DIRECTIONS:
        pair_sums = []
        for owner in range(4):
            signed_sum = roots[owner] + roots[direction[owner]]
            pair_sums.append(signed_sum)
            u_value = signed_sum**2
            require(
                resolvent_value(p_value, q_value, r_value, u_value) == 0,
                ("matching resolvent root", roots, direction, owner),
            )
            if signed_sum != 0:
                a_value = (p_value + u_value + q_value / signed_sum) / 2
                b_value = (p_value + u_value - q_value / signed_sum) / 2
                first = [a_value, -signed_sum, Fraction(1)]
                second = [b_value, signed_sum, Fraction(1)]
                require(
                    multiply_polynomials(first, second)
                    == [r_value, q_value, p_value, Fraction(0), Fraction(1)],
                    ("signed quadratic blocks", roots, direction, owner),
                )
        require(
            Counter(pair_sums) == {pair_sums[0]: 2, -pair_sums[0]: 2}
            or pair_sums[0] == 0,
            ("owner sign fibre", roots, direction, pair_sums),
        )

# Sharp same-quotient split/full and origin-sign controls.
C_REP = (0, 2, 3, 1)  # (1 2 3)
S_SPLIT_REP = transposition(4, 2, 3)
S_FULL_REP = transposition(4, 0, 1)
W_SPLIT_REP = compose(S_SPLIT_REP, C_REP)
W_FULL_REP = compose(S_FULL_REP, C_REP)
require(
    matching_quotient(S_SPLIT_REP) == matching_quotient(S_FULL_REP)
    and power(W_SPLIT_REP, 2) == ONE4
    and power(W_FULL_REP, 2) == DIRECTIONS[1],
    "same-quotient branch hostile",
)

ROOTS_REP = ROOT_CONTROLS[0]
H_REP = DIRECTIONS[0]
S_FULL_CONJ = conjugate(S_FULL_REP, H_REP)
C_CONJ = conjugate(C_REP, H_REP)
require(
    matching_quotient(S_FULL_CONJ) == matching_quotient(S_FULL_REP)
    and matching_quotient(C_CONJ) == matching_quotient(C_REP)
    and power(compose(S_FULL_CONJ, C_CONJ), 2) == power(W_FULL_REP, 2),
    "translation-conjugate full lifts",
)
selected_direction = power(W_FULL_REP, 2)
positive_owner = fixed_points(C_REP)[0]
negative_owner = fixed_points(C_CONJ)[0]
require((positive_owner, negative_owner) == (0, 1), "owner hostile labels")
positive_sum = ROOTS_REP[positive_owner] + ROOTS_REP[selected_direction[positive_owner]]
negative_sum = ROOTS_REP[negative_owner] + ROOTS_REP[selected_direction[negative_owner]]
require(positive_sum == 2 and negative_sum == -2, ("signed origin hostile", positive_sum, negative_sum))

SEMANTIC_ROWS = tuple(
    sorted(
        (
            cycle_string(binary),
            cycle_string(ternary),
            cycle_string(compose(binary, ternary)),
            owner,
            cycle_string(direction),
            orientation,
            branch,
        )
        for binary, ternary, _, (owner, direction, orientation, branch) in LIFTS
    )
)
SEMANTIC_DIGEST = hashlib.sha256(repr(SEMANTIC_ROWS).encode("ascii")).hexdigest()

LINES = [
    "THM3095_MARKED_AFFINE_LIFT_FRAME",
    f"s4_elements={len(S4)} matching_quotient={len(QUOTIENT)} kernel={len(KERNEL)}",
    f"marked_s3_pairs={len(QUOTIENT_PAIRS)} exact_affine_lifts={len(LIFTS)}",
    f"split_lifts={len(SPLIT)} full_lifts={len(FULL)}",
    f"lift_quadruples={len({row[3] for row in LIFTS})}",
    f"full_cooccurrence_triples={len(FULL_TRIPLES)} oriented_flags=4*3*2",
    f"unoriented_flags={len(FLAGS)} orientations_per_flag=2",
    "full_quotient_fibres=6x4_owners all_quotient_fibres=6x4_ownersx2_branches",
    "half_face_formula=w2=eta*d",
    "explicit_completions=split_(b,c(b))_full_(o,c^-1(b))",
    "full_direction_counts=" + ",".join(str(Counter(row[3][1] for row in FULL)[direction]) for direction in DIRECTIONS),
    f"quartic_root_controls={len(ROOT_CONTROLS)} signed_origin_hostile=2_to_-2",
    "representative_split=" + cycle_string(S_SPLIT_REP) + "+" + cycle_string(C_REP) + " word=" + cycle_string(W_SPLIT_REP),
    "representative_full=" + cycle_string(S_FULL_REP) + "+" + cycle_string(C_REP) + " word=" + cycle_string(W_FULL_REP) + " square=" + cycle_string(power(W_FULL_REP, 2)),
    f"semantic_digest={SEMANTIC_DIGEST}",
    "ALL_CHECKS_PASS",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-stored", action="store_true")
    arguments = parser.parse_args()
    payload = "\n".join(LINES) + "\n"
    if arguments.check_stored:
        stored = STORED.read_bytes().replace(b"\r\n", b"\n").decode("utf-8")
        require(stored == payload, "stored transcript changed")
    print(payload, end="")


if __name__ == "__main__":
    main()
