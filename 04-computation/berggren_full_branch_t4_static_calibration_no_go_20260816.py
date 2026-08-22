#!/usr/bin/env python3
"""Exact full-branch calibration test for the Berggren/T4 four-state collision.

The mod-three projective line and the affine V4 owner set both have four
points.  This probe asks whether a single point calibration, or the more
general six-frame family of calibrations, can intertwine every Berggren
branch letter while retaining the letter's true mod-two action on the three
nonzero V4 directions.  It also records the surviving B-only transplant.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/berggren_full_branch_t4_static_calibration_no_go_20260816.py"
OUTPUT = "05-knowledge/results/berggren_full_branch_t4_static_calibration_no_go_20260816.out"

PINS = (
    (
        "THM-3339",
        ROOT / "01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md",
        "1e4aa8cd9d6cc9bf342328a4af0e5db8cd7eefb51eea08460d03f2c6410cee51",
    ),
    (
        "P1(F3)-four-state-output",
        ROOT / "05-knowledge/results/fibonacci_farey_mod3_q15_four_state_probe_20260814.out",
        "803c848ed775b2b1278fecfa1340d2f8f634abaaa2a23ea3e1f24abf76c57e48",
    ),
    (
        "THM-3487",
        ROOT / "01-canon/theorems/THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction.md",
        "a01da997e3f9d9ad8af42b43792309bd887cea86bd0a9aaed1b6cb89014fb2ad",
    ),
)

BRANCH = {
    "A": ((0, 1), (-1, 2)),
    "B": ((0, 1), (1, 2)),
    "C": ((1, 0), (2, 1)),
}
P1_F3 = ((1, 0), (0, 1), (1, 1), (1, 2))
ZERO = (0, 0)
P = (1, 0)
Q = (0, 1)
R = (1, 1)
V4 = (ZERO, P, Q, R)
I2 = ((1, 0), (0, 1))
J2 = ((0, 1), (1, 0))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def mvec(matrix, vector, modulus):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(2)) % modulus
        for row in range(2)
    )


def mmul(left, right, modulus):
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(2))
            % modulus
            for column in range(2)
        )
        for row in range(2)
    )


def vadd(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def same_projective_line(left, right):
    return (left[0] * right[1] - left[1] * right[0]) % 3 == 0


def projective_permutation(matrix):
    result = []
    for vector in P1_F3:
        image = mvec(matrix, vector, 3)
        matches = tuple(
            index for index, representative in enumerate(P1_F3)
            if same_projective_line(image, representative)
        )
        require(len(matches) == 1, (matrix, vector, image, matches))
        result.append(matches[0])
    return tuple(result)


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    return tuple(permutation.index(index) for index in range(len(permutation)))


def conjugate(calibration, action):
    return compose(compose(calibration, action), inverse(calibration))


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
        require(current == start, (start, current, permutation))
        lengths.append(length)
    return tuple(sorted(lengths))


def affine_permutation(linear, translation):
    return tuple(
        V4.index(vadd(mvec(linear, vector, 2), translation))
        for vector in V4
    )


GL2_F2 = tuple(
    ((a, b), (c, d))
    for a, b, c, d in product((0, 1), repeat=4)
    if (a * d - b * c) % 2 == 1
)
FRAMES = tuple(permutations((P, Q, R)))


def frame_permutation(linear):
    return tuple(
        FRAMES.index(tuple(mvec(linear, direction, 2) for direction in frame))
        for frame in FRAMES
    )


def product_permutation(left, right):
    """Diagonal action on the Cartesian product of two indexed sets."""
    right_size = len(right)
    return tuple(
        left[left_index] * right_size + right[right_index]
        for left_index in range(len(left))
        for right_index in range(right_size)
    )


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = deque((identity,))
    while queue:
        element = queue.popleft()
        for generator in generators:
            candidate = compose(generator, element)
            if candidate not in seen:
                seen.add(candidate)
                queue.append(candidate)
    return frozenset(seen)


def family_count_for_one_letter(source_line, target_owner, linear):
    """Count frame-preserving families f_pi satisfying one letter equation.

    The base frame action is induced by ``linear``.  A family consists of one
    point bijection P1(F3)->V4 for each of the six frames.  Propagation around
    every base orbit gives an exact count without searching 24^6 families.
    """
    base = frame_permutation(linear)
    all_calibrations = tuple(permutations(range(4)))
    unseen = set(range(6))
    total = 1
    orbit_data = []
    while unseen:
        start = min(unseen)
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            unseen.remove(current)
            current = base[current]
        require(current == start, (base, orbit, current))
        orbit_length = len(orbit)
        source_return = tuple(range(4))
        target_return = tuple(range(4))
        for _ in range(orbit_length):
            source_return = compose(source_line, source_return)
            target_return = compose(target_owner, target_return)
        seeds = sum(
            conjugate(calibration, source_return) == target_return
            for calibration in all_calibrations
        )
        orbit_data.append((orbit_length, seeds))
        total *= seeds
    return total, tuple(orbit_data)


def main() -> None:
    checked_pins = []
    for label, path, expected in PINS:
        actual = lf_sha256(path)
        require(actual == expected, (label, expected, actual))
        checked_pins.append((label, actual))

    projective = {letter: projective_permutation(matrix) for letter, matrix in BRANCH.items()}
    require(projective == {
        "A": (1, 3, 2, 0),
        "B": (1, 3, 0, 2),
        "C": (3, 1, 0, 2),
    }, projective)

    linear = {
        letter: tuple(tuple(entry % 2 for entry in row) for row in matrix)
        for letter, matrix in BRANCH.items()
    }
    require(linear == {"A": J2, "B": J2, "C": I2}, linear)

    all_calibrations = tuple(permutations(range(4)))
    all_affine = {
        affine_permutation(matrix, translation)
        for matrix in GL2_F2 for translation in V4
    }
    require(len(GL2_F2) == 6, len(GL2_F2))
    require(all_affine == set(all_calibrations), len(all_affine))

    affine_types_by_required_linear = {
        letter: tuple(
            (translation, cycle_lengths(affine_permutation(linear[letter], translation)))
            for translation in V4
        )
        for letter in BRANCH
    }
    source_types = {letter: cycle_lengths(action) for letter, action in projective.items()}

    admissible_targets = {
        letter: tuple(affine_permutation(linear[letter], translation) for translation in V4)
        for letter in BRANCH
    }
    static_counts = {
        letter: sum(
            conjugate(calibration, projective[letter]) in admissible_targets[letter]
            for calibration in all_calibrations
        )
        for letter in BRANCH
    }
    require(static_counts == {"A": 0, "B": 8, "C": 0}, static_counts)

    subset_counts = {}
    for size in range(1, 4):
        for subset in combinations(BRANCH, size):
            subset_counts["".join(subset)] = sum(
                all(
                    conjugate(calibration, projective[letter]) in admissible_targets[letter]
                    for letter in subset
                )
                for calibration in all_calibrations
            )
    require(subset_counts == {
        "A": 0, "B": 8, "C": 0,
        "AB": 0, "AC": 0, "BC": 0, "ABC": 0,
    }, subset_counts)

    fixed_lifts = {
        "A": affine_permutation(J2, P),
        "B": affine_permutation(J2, P),
        "C": affine_permutation(I2, P),
    }
    fixed_static_counts = {
        letter: sum(
            conjugate(calibration, projective[letter]) == fixed_lifts[letter]
            for calibration in all_calibrations
        )
        for letter in BRANCH
    }
    require(fixed_static_counts == {"A": 0, "B": 4, "C": 0}, fixed_static_counts)
    require(tuple(range(4)) in tuple(
        calibration for calibration in all_calibrations
        if conjugate(calibration, projective["B"]) == fixed_lifts["B"]
    ), "the pinned identity B-calibration disappeared")

    source_24 = {}
    target_24_by_translation = {}
    for letter in BRANCH:
        base = frame_permutation(linear[letter])
        source_24[letter] = product_permutation(base, projective[letter])
        target_24_by_translation[letter] = tuple(
            product_permutation(base, affine_permutation(linear[letter], translation))
            for translation in V4
        )

    source_24_types = {letter: cycle_lengths(action) for letter, action in source_24.items()}
    target_24_types = {
        letter: tuple(cycle_lengths(action) for action in actions)
        for letter, actions in target_24_by_translation.items()
    }
    require(source_24_types == {
        "A": (2, 2, 2, 6, 6, 6),
        "B": (4, 4, 4, 4, 4, 4),
        "C": (1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3),
    }, source_24_types)
    require(set(target_24_types["A"]) == {
        (2,) * 12, (4,) * 6,
    }, target_24_types["A"])
    require(set(target_24_types["B"]) == {
        (2,) * 12, (4,) * 6,
    }, target_24_types["B"])
    require(set(target_24_types["C"]) == {
        (1,) * 24, (2,) * 12,
    }, target_24_types["C"])

    family_counts = {
        letter: tuple(
            family_count_for_one_letter(
                projective[letter],
                affine_permutation(linear[letter], translation),
                linear[letter],
            )[0]
            for translation in V4
        )
        for letter in BRANCH
    }
    require(family_counts == {
        "A": (0, 0, 0, 0),
        "B": (0, 512, 512, 0),
        "C": (0, 0, 0, 0),
    }, family_counts)

    source_group = generated_group(tuple(projective.values()))
    fixed_target_group = generated_group(tuple(fixed_lifts.values()))
    require(len(source_group) == 24, len(source_group))
    require(len(fixed_target_group) == 8, len(fixed_target_group))

    # With the identity point calibration, recover each mod-three branch
    # permutation's unique affine decomposition.  A and C necessarily use an
    # order-three linear part, revealing exactly what the true mod-two action
    # would have to sacrifice.
    affine_decomposition = {}
    for letter, action in projective.items():
        candidates = tuple(
            (matrix, translation)
            for matrix in GL2_F2 for translation in V4
            if affine_permutation(matrix, translation) == action
        )
        require(len(candidates) == 1, (letter, candidates))
        matrix, translation = candidates[0]
        affine_decomposition[letter] = (
            matrix,
            translation,
            cycle_lengths(tuple(V4.index(mvec(matrix, vector, 2)) for vector in V4)),
        )

    semantic = {
        "projective": projective,
        "linear": linear,
        "source_types": source_types,
        "static_counts": static_counts,
        "fixed_static_counts": fixed_static_counts,
        "subset_counts": subset_counts,
        "source_24_types": source_24_types,
        "target_24_types": target_24_types,
        "family_counts": family_counts,
        "group_orders": (len(source_group), len(fixed_target_group)),
        "affine_decomposition": affine_decomposition,
    }
    semantic_sha256 = sha256(
        dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("BERGGREN FULL-BRANCH / T4 STATIC CALIBRATION EXACT PROBE")
    print("STATUS: VERIFIED-EXACT STRUCTURAL SIDECAR; NO LRC OR JACOBIAN CONSEQUENCE")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PINS: {tuple(checked_pins)}")
    print(f"P1_BRANCH_PERMUTATIONS: {projective}")
    print(f"TRUE_MOD2_LINEAR_PARTS: {linear}")
    print(f"P1_CYCLE_TYPES: {source_types}")
    print(f"AFFINE_TYPES_BY_TRUE_LINEAR_PART: {affine_types_by_required_linear}")
    print(f"STATIC_CALIBRATIONS_ALLOWING_TRANSLATION: {static_counts}")
    print(f"STATIC_SUBSET_COUNTS: {subset_counts}")
    print(f"THM3339_FIXED_LIFT_STATIC_COUNTS: {fixed_static_counts}")
    print("B_POSITIVE_CONTROL: under x0,x1,x2,x3 <-> 0,p,q,r, the projective B action is exactly x |-> Jx+p; four static gauges conjugate to this pinned lift, while eight work if either lawful four-cycle translation is allowed")
    print(f"SOURCE_24_CYCLE_TYPES: {source_24_types}")
    print(f"TARGET_24_CYCLE_TYPES_BY_TRANSLATION: {target_24_types}")
    print(f"BASE_FRAME_PRESERVING_FAMILY_COUNTS_BY_TRANSLATION_0_P_Q_R: {family_counts}")
    print("B_24_STATE_POSITIVE_CONTROL: each lawful four-cycle translation supports 8^3=512 calibration families, one independent return-compatible seed on each of the three swap-orbits of the six frame states")
    print(f"IMAGE_ORDER_OBSTRUCTION: projective_branch_image={len(source_group)}=|S4|; pinned_affine_branch_image={len(fixed_target_group)}=|V4 semidirect <J>|")
    print(f"UNCONSTRAINED_AFFINENESS_HOSTILE: AGL(2,2)=S4={len(all_affine)} permutations, so every four-point action is affine if its true mod-two linear part is forgotten")
    print(f"IDENTITY_GAUGE_AFFINE_DECOMPOSITIONS: {affine_decomposition}")
    print("MISSING_SIDECAR: the A and C projective three-cycles force an order-three affine linear part; translations over the prescribed J/I direction actions cannot supply that 3-primary mode")
    print("T4_TYPING: both 6*4 carriers index all 24 labelled transitive T4s, but equal atlas cardinality does not identify their branch connections; the B submonoid transplants, the full <A,B,C> monoid does not")
    print(f"SEMANTIC_SHA256: {semantic_sha256}")
    print("VERDICT: no static point calibration and no six-frame base-preserving calibration intertwines the full ternary Berggren action with affine owner transport while preserving the true mod-two direction quotient; B is the exact maximal one-letter positive control")


if __name__ == "__main__":
    main()
