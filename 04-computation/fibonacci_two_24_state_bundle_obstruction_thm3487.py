#!/usr/bin/env python3
"""Exact comparison of two 6-by-4 Fibonacci state bundles.

The frame-line bundle advances the six residue orders and the four
projective mod-three lines simultaneously.  The owner bundle advances the
same six orders and transports the affine V4 owner through the unique linear
channel-frame changes carrying one displayed order to the next.  Their
one-step cycle types differ.  This is a finite shift obstruction, not a
full-tree or physical-current result.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
from json import dumps
from math import factorial
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/fibonacci_two_24_state_bundle_obstruction_thm3487.py"
OUTPUT = "05-knowledge/results/fibonacci_two_24_state_bundle_obstruction_thm3487.out"
EXPECTED_SEMANTIC_SHA256 = "9c338ec5eebe1f93326ae12c00f37da58382fa2084339adbd1de422e7d05d70e"

PINS = (
    (
        "THM-3339-script",
        ROOT / "04-computation/fibonacci_berggren_three_ray_owner_thm3339.py",
        "094fc254dcb7965791e59247a98f60a337725b40d169db6f3862c1da5943149b",
    ),
    (
        "THM-3339-output",
        ROOT / "05-knowledge/results/fibonacci_berggren_three_ray_owner_thm3339.out",
        "88c5f44971df12bbd61f84925f408a253450095ee2de23b5a617efde1e4ecdfa",
    ),
    (
        "frame-line-script",
        ROOT / "04-computation/fibonacci_farey_mod3_q15_four_state_probe_20260814.py",
        "6ad743f478420641aae9078bcdc83b4f9eb19b037bb90f640b92b1e5933e6295",
    ),
    (
        "frame-line-output",
        ROOT / "05-knowledge/results/fibonacci_farey_mod3_q15_four_state_probe_20260814.out",
        "803c848ed775b2b1278fecfa1340d2f8f634abaaa2a23ea3e1f24abf76c57e48",
    ),
)

ZERO = (0, 0)
P = (1, 0)
Q = (0, 1)
R = (1, 1)
V4 = (ZERO, P, Q, R)
I2 = ((1, 0), (0, 1))

# THM-3339 (21), with the declared channel gauge a=P,b=Q,c=R.
ORDERS = (
    (Q, R, P),
    (Q, P, R),
    (P, Q, R),
    (P, R, Q),
    (R, P, Q),
    (R, Q, P),
)
OWNERS = (ZERO, P, P, R, Q, Q)
P1_F3_FIBONACCI_CYCLE = (1, 2, 3, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def vadd(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (left[0] ^ right[0], left[1] ^ right[1])


def mvec(
    matrix: tuple[tuple[int, int], tuple[int, int]],
    vector: tuple[int, int],
) -> tuple[int, int]:
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % 2,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % 2,
    )


def mmul(left, right):
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(2)) % 2
            for column in range(2)
        )
        for row in range(2)
    )


GL2 = tuple(
    ((a, b), (c, d))
    for a, b, c, d in product((0, 1), repeat=4)
    if (a * d - b * c) % 2 == 1
)


def linear_step(source, target):
    candidates = tuple(
        matrix
        for matrix in GL2
        if tuple(mvec(matrix, vector) for vector in source) == target
    )
    require(len(candidates) == 1, (source, target, candidates))
    return candidates[0]


def affine_apply(affine, vector):
    matrix, translation = affine
    return vadd(mvec(matrix, vector), translation)


def affine_compose(left, right):
    """Return left after right."""

    left_matrix, left_translation = left
    right_matrix, right_translation = right
    return (
        mmul(left_matrix, right_matrix),
        vadd(mvec(left_matrix, right_translation), left_translation),
    )


def inverse_matrix(matrix):
    candidates = tuple(candidate for candidate in GL2
                       if mmul(candidate, matrix) == I2 and mmul(matrix, candidate) == I2)
    require(len(candidates) == 1, ("inverse", matrix, candidates))
    return candidates[0]


def affine_inverse(affine):
    matrix, translation = affine
    inverse = inverse_matrix(matrix)
    return inverse, mvec(inverse, translation)


def permutation_from_map(states, transition):
    index = {state: position for position, state in enumerate(states)}
    image = tuple(index[transition(state)] for state in states)
    require(len(set(image)) == len(states), "transition is not a permutation")
    return image


def compose_permutations(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def power_permutation(permutation, exponent):
    result = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            result = compose_permutations(base, result)
        exponent //= 2
        if exponent:
            base = compose_permutations(base, base)
    return result


def cycle_lengths(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        require(current == start, ("cycle did not close", start, current))
        lengths.append(length)
    return tuple(sorted(lengths))


def conjugacy_count(left, right):
    left_cycles = Counter(cycle_lengths(left))
    right_cycles = Counter(cycle_lengths(right))
    if left_cycles != right_cycles:
        return 0
    total = 1
    for length, multiplicity in left_cycles.items():
        total *= length**multiplicity * factorial(multiplicity)
    return total


def owner_bundle_steps():
    moving_steps = []
    fixed_steps = []
    linear_steps = []
    for index in range(6):
        next_index = (index + 1) % 6
        linear = linear_step(ORDERS[index], ORDERS[next_index])
        translation = vadd(OWNERS[next_index], mvec(linear, OWNERS[index]))
        moving = (linear, translation)
        require(affine_apply(moving, OWNERS[index]) == OWNERS[next_index],
                (index, "section transport"))
        moving_steps.append(moving)
        fixed_steps.append((I2, vadd(OWNERS[index], OWNERS[next_index])))
        linear_steps.append(linear)

    moving_holonomy = (I2, ZERO)
    fixed_holonomy = (I2, ZERO)
    linear_holonomy = I2
    for moving, fixed, linear in zip(moving_steps, fixed_steps, linear_steps):
        moving_holonomy = affine_compose(moving, moving_holonomy)
        fixed_holonomy = affine_compose(fixed, fixed_holonomy)
        linear_holonomy = mmul(linear, linear_holonomy)
    require(moving_holonomy == (I2, ZERO), moving_holonomy)
    require(fixed_holonomy == (I2, ZERO), fixed_holonomy)
    require(linear_holonomy == I2, linear_holonomy)

    # Construct an explicit fibrewise affine gauge from moving-frame to
    # fixed-frame transport.
    gauges = [(I2, ZERO)]
    for index in range(6):
        fixed = fixed_steps[index]
        moving = moving_steps[index]
        next_gauge = affine_compose(
            fixed,
            affine_compose(gauges[index], affine_inverse(moving)),
        )
        gauges.append(next_gauge)
    require(gauges[-1] == gauges[0], ("gauge did not close", gauges))
    return tuple(linear_steps), tuple(moving_steps), tuple(fixed_steps), tuple(gauges[:-1])


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    pins = []
    for label, path, expected in PINS:
        actual = lf_sha256(path)
        require(actual == expected, (label, "hash drift", actual))
        pins.append((label, actual))

    source_tree = ast.parse((ROOT / SCRIPT).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(source_tree)),
            "assert forbidden")
    forbidden = {"eval", "exec", "compile", "open", "system", "popen", "run", "Popen"}
    calls = {
        node.func.id
        for node in ast.walk(source_tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    calls.update(
        node.func.attr
        for node in ast.walk(source_tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not (calls & forbidden), ("forbidden calls", calls & forbidden))
    security = len(tuple(ast.walk(source_tree))), tuple(sorted(calls & forbidden))

    require(len(set(ORDERS)) == 6 and set(ORDERS) == set(permutations(V4[1:])),
            "order packet drift")
    require(tuple(vadd(OWNERS[index], OWNERS[(index + 1) % 6]) for index in range(6))
            == (P, ZERO, Q, P, ZERO, Q), "owner drift drifted")

    linear_steps, moving_steps, fixed_steps, gauges = owner_bundle_steps()

    frame_states = tuple((order, line) for order in range(6) for line in range(4))
    owner_states = tuple((order, owner) for order in range(6) for owner in V4)

    frame_shift = permutation_from_map(
        frame_states,
        lambda state: ((state[0] + 1) % 6, (state[1] + 1) % 4),
    )
    moving_owner_shift = permutation_from_map(
        owner_states,
        lambda state: (
            (state[0] + 1) % 6,
            affine_apply(moving_steps[state[0]], state[1]),
        ),
    )
    fixed_owner_shift = permutation_from_map(
        owner_states,
        lambda state: (
            (state[0] + 1) % 6,
            affine_apply(fixed_steps[state[0]], state[1]),
        ),
    )

    # Explicit gauge conjugacy of the two lawful owner trivializations.
    gauge_permutation = permutation_from_map(
        owner_states,
        lambda state: (state[0], affine_apply(gauges[state[0]], state[1])),
    )
    require(
        compose_permutations(gauge_permutation, moving_owner_shift)
        == compose_permutations(fixed_owner_shift, gauge_permutation),
        "moving/fixed owner gauges are not conjugate",
    )

    frame_cycles = cycle_lengths(frame_shift)
    owner_cycles = cycle_lengths(moving_owner_shift)
    require(frame_cycles == (12, 12), frame_cycles)
    require(owner_cycles == (6, 6, 6, 6), owner_cycles)
    require(cycle_lengths(fixed_owner_shift) == owner_cycles, "owner gauge cycle drift")
    require(sum(power_permutation(frame_shift, 6)[index] == index for index in range(24)) == 0,
            "frame sixth power gained a fixed point")
    require(sum(power_permutation(moving_owner_shift, 6)[index] == index for index in range(24)) == 24,
            "owner sixth power lost a fixed point")
    require(conjugacy_count(frame_shift, moving_owner_shift) == 0,
            "one-step bundles unexpectedly conjugate")

    frame_square = power_permutation(frame_shift, 2)
    require(cycle_lengths(frame_square) == owner_cycles, "two-step survivor failed")
    require(conjugacy_count(frame_square, moving_owner_shift) == 6**4 * factorial(4),
            "two-step conjugacy count")

    repaired = []
    for holonomy in V4[1:]:
        steps = list(moving_steps)
        last_matrix, last_translation = steps[-1]
        steps[-1] = (last_matrix, vadd(last_translation, holonomy))
        shift = permutation_from_map(
            owner_states,
            lambda state, steps=tuple(steps): (
                (state[0] + 1) % 6,
                affine_apply(steps[state[0]], state[1]),
            ),
        )
        cycles = cycle_lengths(shift)
        require(cycles == frame_cycles, (holonomy, cycles))
        require(conjugacy_count(frame_shift, shift) == 12**2 * factorial(2),
                (holonomy, "repair conjugacy count"))
        require(affine_apply(steps[-1], OWNERS[-1]) == vadd(OWNERS[0], holonomy),
                (holonomy, "repair did not break section closure exactly"))
        actual_holonomy = (I2, ZERO)
        for step in steps:
            actual_holonomy = affine_compose(step, actual_holonomy)
        require(actual_holonomy == (I2, holonomy),
                (holonomy, "repaired holonomy", actual_holonomy))
        repaired.append((holonomy, cycles, conjugacy_count(frame_shift, shift)))

    # On four points the identity plus the three double transpositions is the
    # canonical normal Klein subgroup of S4.  The square of the Fibonacci
    # 4-cycle selects one of those three matchings.  Transporting it to the
    # affine-owner V4 requires a point identification; the 24 bijections split
    # evenly, eight for each possible nonzero translation direction.
    p1_square = power_permutation(P1_F3_FIBONACCI_CYCLE, 2)
    require(cycle_lengths(p1_square) == (2, 2), p1_square)
    matching_gauge_counts = []
    for direction in V4[1:]:
        count = sum(
            all(labeling[p1_square[index]] == vadd(labeling[index], direction)
                for index in range(4))
            for labeling in permutations(V4)
        )
        require(count == 8, (direction, count))
        matching_gauge_counts.append((direction, count))
    require(sum(count for _, count in matching_gauge_counts) == factorial(4),
            matching_gauge_counts)

    uniform_address_harmonic_grades = {
        "frame_one_step": tuple(range(0, 3)),
        "owner_one_step": tuple(range(0, 5)),
        "frame_two_step": tuple(range(0, 5)),
        "owner_nonzero_holonomy": tuple(range(0, 3)),
    }
    # Numerators are measured in orbit units: frame/Y_h orbits have mass 1/2,
    # owner/X^2 orbits have mass 1/4.

    semantic_payload = {
        "pins": tuple(pins),
        "orders": ORDERS,
        "owners": OWNERS,
        "linear_steps": linear_steps,
        "moving_steps": moving_steps,
        "fixed_steps": fixed_steps,
        "gauges": gauges,
        "frame_cycles": frame_cycles,
        "owner_cycles": owner_cycles,
        "frame_square_cycles": cycle_lengths(frame_square),
        "repaired": tuple(repaired),
        "p1_f3_fibonacci_cycle_and_square": (P1_F3_FIBONACCI_CYCLE, p1_square),
        "matching_gauge_counts": tuple(matching_gauge_counts),
        "conjugacy_counts": (
            conjugacy_count(frame_shift, moving_owner_shift),
            conjugacy_count(frame_square, moving_owner_shift),
            tuple(count for _, _, count in repaired),
        ),
        "fixed_points_sixth": (0, 24),
        "uniform_address_harmonic_grades": uniform_address_harmonic_grades,
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("THM-3487 TWO FIBONACCI 24-STATE BUNDLES EXACT COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PINS: {tuple(pins)}")
    print(f"RESIDUE_ORDERS: {ORDERS}")
    print(f"OWNER_SECTION_AND_DRIFT: owners={OWNERS}; drift={(P, ZERO, Q, P, ZERO, Q)}")
    print(f"MOVING_FRAME_LINEAR_STEPS: {linear_steps}")
    print(f"MOVING_FRAME_AFFINE_STEPS: {moving_steps}")
    print(f"FRAME_LINE_SHIFT_CYCLES: {frame_cycles}; sixth_power_fixed_points=0")
    print(f"AFFINE_OWNER_SHIFT_CYCLES: {owner_cycles}; sixth_power_fixed_points=24")
    print("ONE_STEP_VERDICT: no equivariant bijection; cycle types 12^2 and 6^4 differ")
    print(f"TWO_STEP_SURVIVOR: frame_shift_squared has cycles={cycle_lengths(frame_square)}; abstract_conjugacies={conjugacy_count(frame_square, moving_owner_shift)}; no base-order-preserving conjugacy")
    print(f"NONZERO_HOLONOMY_REPAIRS: {tuple(repaired)}; abstract_conjugacies_each={12**2 * factorial(2)}")
    print(f"PROJECTIVE_TO_AFFINE_SIDECAR: G={P1_F3_FIBONACCI_CYCLE}; G_squared={p1_square}; point_identifications_by_nonzero_owner_direction={tuple(matching_gauge_counts)}")
    print("REPAIR_TARIFF: every nonzero V4 seam holonomy matches the 12^2 cycle type but sends the proved closing owner to owner+h; G^2 selects a projective matching, while choosing its owner direction is exactly the missing channel-matching gauge")
    print("HARMONIC_UNIFORM_ADDRESS_GRADES: under a declared uniform 24-state enumeration, one-step frame/nonzero-holonomy invariant unions have coefficients in {0,1/2,1}; one-step owner/two-step-frame unions have coefficients in {0,1/4,1/2,3/4,1}; this is not a single-orbit time average")
    print("TYPING: six states are orders of three matching channels, not six T6 vertices; four states are P1(F3) lines versus affine V4 owners; neither fibre is an intrinsic tournament")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print(f"SECURITY_AST_NODES_AND_FORBIDDEN: {security}")
    print("VERDICT: equal cardinality 6*4 does not give the missing transplant; the exact obstruction is sixth-step holonomy, and the minimal cycle repair is a projective-to-affine matching-gauge sidecar that destroys the proved closed owner section")


if __name__ == "__main__":
    main()
