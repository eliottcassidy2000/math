#!/usr/bin/env python3
"""Exact probe for THM-3509's fraction/window/K4/Fibonacci bridge.

The script uses only the Python standard library.  It checks the two
Berggren parity trees, Stern--Brocot mediant linearization on four-term
recurrence windows, the primitive harmonic-face converse, the octahedral
edge completion, the Euclid current and parity content, the Fibonacci
unit-Cassini classification, and hostile quotient controls.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd, isqrt
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY_HASHES = {
    "01-canon/theorems/THM-2596-modular-free-factor-farey-gram-owner-cocycle.md":
        "6adabe6c52dcd15bd0b73d4349470e6a923b1776d0741b7ecda0dbc407873114",
    "01-canon/theorems/THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md":
        "5400f79449a57276a9127a742420ca82986aa6d631b97c9969c73724b427f9b6",
    "01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md":
        "68d2b2d69316a07e3939471ffc753196b090a261b21526acc1f7cd15d5cc0391",
    "01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md":
        "1e4aa8cd9d6cc9bf342328a4af0e5db8cd7eefb51eea08460d03f2c6410cee51",
    "01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md":
        "ba24accf81d123d76ceee2ea7332d394d9ff23d6c9a0d47c7c76ab5b3ad9d446",
    "01-canon/theorems/THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order.md":
        "2f55187055b8158f5a99fd154df14beec5c2b6e7a0a4f65bc5995d495be7d058",
}
RELATED_HASHES = {
    "01-canon/theorems/THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary.md":
        "2bea649cefd56e52c47f92e759e52ec32efda31a1a0b17486a5c0b8bb0f2fa8a",
}
EXPECTED_SEMANTIC_SHA256 = "7ca802f1f8706658b62fb330504544b9c46d473db6481429cb011b078695a46e"

VERTICES = tuple(range(4))
EDGES = tuple(combinations(VERTICES, 2))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
MATCHING_INDEX_PAIRS = (
    (EDGE_INDEX[(0, 1)], EDGE_INDEX[(2, 3)]),
    (EDGE_INDEX[(0, 2)], EDGE_INDEX[(1, 3)]),
    (EDGE_INDEX[(0, 3)], EDGE_INDEX[(1, 2)]),
)


def normalized_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected_hash in DEPENDENCY_HASHES.items():
    actual_hash = normalized_sha256(REPO_ROOT / dependency)
    require(actual_hash == expected_hash, f"dependency hash changed: {dependency}")
for related, expected_hash in RELATED_HASHES.items():
    actual_hash = normalized_sha256(REPO_ROOT / related)
    require(actual_hash == expected_hash, f"related source hash changed: {related}")


def gcd3(a: int, b: int, c: int) -> int:
    return gcd(gcd(a, b), c)


def window(m: int, n: int) -> tuple[int, int, int, int]:
    """The additive-recurrence window attached to m/n."""
    return (n - m, m, n, n + m)


def fraction_from_seed(seed: tuple[int, int]) -> tuple[int, int]:
    x, y = seed
    return (y, x + y)


def edge_weights(values: tuple[int, int, int, int]) -> tuple[int, ...]:
    return tuple(values[left] * values[right] for left, right in EDGES)


def harmonic_face(weights: tuple[int, ...]) -> tuple[int, int, int]:
    return (
        weights[EDGE_INDEX[(0, 1)]],
        weights[EDGE_INDEX[(0, 2)]],
        weights[EDGE_INDEX[(1, 2)]],
    )


def euclid_triple(m: int, n: int) -> tuple[int, int, int]:
    return (n * n - m * m, 2 * m * n, n * n + m * m)


def triple_from_face(face: tuple[int, int, int]) -> tuple[int, int, int]:
    u, v, z = face
    return (u + v, 2 * z, v + 2 * z - u)


def primitive_unordered(triple: tuple[int, int, int]) -> tuple[int, int, int]:
    common = gcd3(*triple)
    a, b, c = (value // common for value in triple)
    return (min(a, b), max(a, b), c)


def branch(seed: tuple[int, int], letter: str) -> tuple[int, int]:
    x, y = seed
    if letter == "A":
        return (x, x + y)
    if letter == "B":
        return (x + 2 * y, x + y)
    if letter == "C":
        return (x + 2 * y, y)
    raise RuntimeError(f"unknown branch: {letter}")


def keller_seed_update(seed: tuple[int, int]) -> tuple[int, int]:
    """THM-3506's conditional face update in x=e-m,y=m coordinates."""
    x, y = seed
    return (4 * x + 4 * y, 3 * x + y)


def parent_step(seed: tuple[int, int]) -> tuple[tuple[int, int], str] | None:
    x, y = seed
    if x == y:
        require(seed == (1, 1), "nonprimitive 1/2 seam reached")
        return None
    if x == 2 * y:
        require(seed == (2, 1), "nonprimitive 1/3 seam reached")
        return None
    if y > x:
        parent, letter = (x, y - x), "A"
    elif y < x < 2 * y:
        parent, letter = (2 * y - x, x - y), "B"
    elif x > 2 * y:
        parent, letter = (x - 2 * y, y), "C"
    else:
        raise RuntimeError(f"unclassified positive seed: {seed}")
    require(min(parent) > 0, "Berggren parent left the positive cone")
    require(sum(parent) < sum(seed), "Berggren inverse did not descend")
    require(gcd(*parent) == 1, "Berggren inverse lost primitivity")
    require(parent[0] % 2 == x % 2, "Berggren inverse changed the parity lane")
    require(branch(parent, letter) == seed, "Berggren branch inverse failed")
    return parent, letter


def descend(seed: tuple[int, int]) -> tuple[tuple[int, int], tuple[str, ...]]:
    current = seed
    child_to_parent: list[str] = []
    while True:
        step = parent_step(current)
        if step is None:
            root = current
            break
        current, letter = step
        child_to_parent.append(letter)
    word = tuple(reversed(child_to_parent))
    rebuilt = root
    for letter in word:
        rebuilt = branch(rebuilt, letter)
    require(rebuilt == seed, "Berggren word did not rebuild its seed")
    return root, word


def parity_pair(seed: tuple[int, int]) -> tuple[int, int]:
    """Primitive seed form of r -> (1-r)/(1+r)."""
    x, y = seed
    if x % 2:
        return (2 * y, x)
    require(y % 2 == 1, "coprime even-x seed did not have odd y")
    return (y, x // 2)


def fib_pair_fast(n: int) -> tuple[int, int]:
    if n == 0:
        return (0, 1)
    a, b = fib_pair_fast(n // 2)
    c = a * (2 * b - a)
    d = a * a + b * b
    if n % 2:
        return (d, c + d)
    return (c, d)


def is_rational_square(value: Fraction) -> bool:
    if value <= 0:
        return False
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    return (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    )


semantic_rows: list[tuple[object, ...]] = []

# All reduced fractions through a complete denominator box, including the
# two parity trees, the harmonic face, the full octahedral lift, and the
# global leg-swap pairing.
FRACTION_DENOMINATOR_BOUND = 240
fraction_count = 0
odd_lane_count = 0
even_lane_count = 0
max_descent_depth = 0
for n in range(2, FRACTION_DENOMINATOR_BOUND + 1):
    for m in range(1, n):
        if gcd(m, n) != 1:
            continue
        fraction_count += 1
        x, y = n - m, m
        values = window(m, n)
        require(values == (x, y, x + y, x + 2 * y), "window recurrence changed")
        require(gcd(x, y) == 1, "reduced fraction did not give a primitive seed")

        weights = edge_weights(values)
        if (x, y) == (1, 1):
            require(len(set(values)) < 4 and len(set(weights)) < 6, "root tie disappeared")
        else:
            require(len(set(values)) == 4, "nonroot K4 vertex comparison acquired a tie")
            require(len(set(weights)) == 6, "nonroot octahedral edge comparison acquired a tie")
        u, v, z = harmonic_face(weights)
        require((u, v, z) == (x * y, x * (x + y), y * (x + y)), "harmonic face changed")
        require(v * z == u * (v + z), "harmonic equation failed")
        require(gcd3(u, v, z) == 1, "harmonic face was not primitive")
        require((gcd(u, v), gcd(u, z), gcd(v, z)) == (x, y, x + y), "gcd inverse changed")
        require(u * n == v * m, "harmonic face did not recover m/n")
        require(u * n == z * (n - m), "harmonic face did not recover the complement")

        w01, w02, w03, w12, w13, w23 = weights
        require(w03 == u + v, "first opposite-face extension changed")
        require(w13 == 2 * z - u, "second opposite-face extension changed")
        require(w23 == v + 2 * z == w03 + w13, "third opposite-face extension changed")
        require(w01 * w23 == w02 * w13 == w03 * w12, "octahedral antipodal products changed")

        raw_triple = euclid_triple(m, n)
        require(triple_from_face((u, v, z)) == raw_triple, "linear harmonic-to-Pythagorean map changed")
        a, b, c = raw_triple
        require(c * c == a * a + b * b, "Pythagorean consequence failed")
        expected_content = 1 if x % 2 else 2
        require(gcd3(a, b, c) == expected_content, "Euclid parity content changed")

        root, word = descend((x, y))
        expected_root = (1, 1) if x % 2 else (2, 1)
        require(root == expected_root, "Berggren parity tree reached the wrong root")
        max_descent_depth = max(max_descent_depth, len(word))
        if x % 2:
            odd_lane_count += 1
        else:
            even_lane_count += 1
            next_x, next_y = keller_seed_update((x, y))
            require(next_x % 2 == 0 and next_y % 2 == 1, "Keller sidecar left the even-x parity lane")
            require(gcd(next_x, next_y) == 1, "Keller sidecar lost primitivity on the even-x lane")

        paired = parity_pair((x, y))
        require(parity_pair(paired) == (x, y), "reduced parity involution failed")
        paired_m, paired_n = fraction_from_seed(paired)
        require(paired_m * (n + m) == (n - m) * paired_n, "fraction involution formula changed")
        paired_triple = euclid_triple(paired_m, paired_n)
        require(primitive_unordered(paired_triple) == primitive_unordered(raw_triple), "leg-swap pair lost its triple")
        require(gcd3(*paired_triple) == (2 if expected_content == 1 else 1), "parity involution did not swap content lanes")

        branch_swap = {"A": "C", "B": "B", "C": "A"}
        for letter in "ABC":
            require(
                parity_pair(branch((x, y), letter)) == branch(paired, branch_swap[letter]),
                "parity involution did not conjugate the outer branches",
            )

        semantic_rows.append(((m, n), values, (u, v, z), raw_triple, expected_root, word))

require(odd_lane_count + even_lane_count == fraction_count, "fraction parity census changed")

# The linear face-to-current formula has the harmonic equation as its exact
# Pythagorean defect, even away from the harmonic locus.
for u in range(1, 20):
    for v in range(1, 20):
        for z in range(1, 20):
            a, b, c = triple_from_face((u, v, z))
            require(c * c - a * a - b * b == 4 * (v * z - u * (v + z)), "Pythagorean defect identity changed")

# Independent forward generation of both ternary trees.
BERGGREN_DEPTH = 8
level = {(1, 1), (2, 1)}
all_berggren_nodes = set(level)
berggren_level_counts = [len(level)]
for depth in range(1, BERGGREN_DEPTH + 1):
    next_level = {branch(seed, letter) for seed in level for letter in "ABC"}
    require(len(next_level) == 2 * 3**depth, "Berggren level collision or omission")
    require(all_berggren_nodes.isdisjoint(next_level), "Berggren tree repeated an ancestor")
    for seed in next_level:
        root, word = descend(seed)
        require(len(word) == depth, "Berggren generated depth disagreed with descent")
        require(root == ((1, 1) if seed[0] % 2 else (2, 1)), "generated node crossed parity trees")
    all_berggren_nodes.update(next_level)
    level = next_level
    berggren_level_counts.append(len(level))

# Stern--Brocot mediants become componentwise addition of four-vertex
# recurrence windows, including the two boundary windows.
STERN_BROCOT_DEPTH = 12
stern_brocot_nodes: set[tuple[int, int]] = set()


def visit_stern_brocot(
    left: tuple[int, int], right: tuple[int, int], remaining_depth: int
) -> None:
    lm, ln = left
    rm, rn = right
    require(lm * rn - ln * rm == -1, "Stern--Brocot flank lost unimodularity")
    middle = (lm + rm, ln + rn)
    require(gcd(*middle) == 1, "Stern--Brocot mediant was not reduced")
    expected_window = tuple(a + b for a, b in zip(window(*left), window(*right)))
    require(window(*middle) == expected_window, "window did not linearize the mediant")
    require(middle not in stern_brocot_nodes, "Stern--Brocot node repeated")
    stern_brocot_nodes.add(middle)
    if remaining_depth > 1:
        visit_stern_brocot(left, middle, remaining_depth - 1)
        visit_stern_brocot(middle, right, remaining_depth - 1)


visit_stern_brocot((0, 1), (1, 1), STERN_BROCOT_DEPTH)
require(len(stern_brocot_nodes) == 2**STERN_BROCOT_DEPTH - 1, "Stern--Brocot census changed")

# All Farey-neighbor pairs in a denominator box, with the mixed polarization
# which the quadratic six-edge carrier needs after window addition.
FAREY_DENOMINATOR_BOUND = 64
farey_fractions = [
    (m, n)
    for n in range(2, FAREY_DENOMINATOR_BOUND + 1)
    for m in range(1, n)
    if gcd(m, n) == 1
]
farey_edge_count = 0
for left, right in combinations(farey_fractions, 2):
    lm, ln = left
    rm, rn = right
    determinant = lm * rn - ln * rm
    if abs(determinant) != 1:
        continue
    farey_edge_count += 1
    middle = (lm + rm, ln + rn)
    left_window, right_window = window(*left), window(*right)
    middle_window = window(*middle)
    require(middle_window == tuple(a + b for a, b in zip(left_window, right_window)), "Farey window addition failed")
    left_seed = (ln - lm, lm)
    right_seed = (rn - rm, rm)
    seed_determinant = left_seed[0] * right_seed[1] - left_seed[1] * right_seed[0]
    require(seed_determinant == -determinant, "seed/Farey determinant gauge changed")
    left_edges = edge_weights(left_window)
    right_edges = edge_weights(right_window)
    middle_edges = edge_weights(middle_window)
    for index, (i, j) in enumerate(EDGES):
        cross = left_window[i] * right_window[j] + left_window[j] * right_window[i]
        require(middle_edges[index] == left_edges[index] + right_edges[index] + cross, "edge polarization failed")

# Independent exhaustive converse for primitive positive harmonic triples.
HARMONIC_EDGE_BOUND = 600
primitive_harmonic_count = 0
for v in range(1, HARMONIC_EDGE_BOUND + 1):
    for z in range(1, HARMONIC_EDGE_BOUND + 1):
        numerator, denominator = v * z, v + z
        if numerator % denominator:
            continue
        u = numerator // denominator
        if gcd3(u, v, z) != 1:
            continue
        primitive_harmonic_count += 1
        common = gcd(v, z)
        x, y = v // common, z // common
        require(gcd(x, y) == 1, "harmonic converse parameters were not coprime")
        require(common == x + y, "harmonic converse common factor changed")
        require((u, v, z) == (x * y, x * (x + y), y * (x + y)), "harmonic converse failed")
        require((gcd(u, v), gcd(u, z), gcd(v, z)) == (x, y, x + y), "harmonic converse gcd decoder failed")

require(4 * 4 == 2 * (4 + 4), "scaled harmonic hostile was not harmonic")
require(gcd3(2, 4, 4) == 2, "scaled harmonic hostile became primitive")

# The equal-antipodal-product toric equations need a square-class sidecar over
# Q.  Six edge weights all equal to 2 satisfy the equations but force e0^2=2.
equal_two_edges = (2,) * 6
require(
    equal_two_edges[0] * equal_two_edges[5]
    == equal_two_edges[1] * equal_two_edges[4]
    == equal_two_edges[2] * equal_two_edges[3],
    "toric square-class hostile left the antipodal-product locus",
)
forced_e0_square = Fraction(equal_two_edges[0] * equal_two_edges[1], equal_two_edges[3])
require(forced_e0_square == 2 and not is_rational_square(forced_e0_square), "toric square-class hostile failed")

# Fibonacci recurrence and an independent fast-doubling implementation.
FIBONACCI_MAX_K = 120
fibonacci = [0, 1]
for _ in range(2, 2 * FIBONACCI_MAX_K + 3):
    fibonacci.append(fibonacci[-1] + fibonacci[-2])
for index, value in enumerate(fibonacci):
    require(fib_pair_fast(index)[0] == value, "fast-doubling Fibonacci path disagreed")

golden_rows = []
for k in range(2, FIBONACCI_MAX_K + 1):
    m, n = fibonacci[k], fibonacci[k + 1]
    values = window(m, n)
    expected_values = tuple(fibonacci[k - 1 + offset] for offset in range(4))
    require(values == expected_values, "golden window ceased to be consecutive Fibonacci")
    weights = edge_weights(values)
    u, v, z = harmonic_face(weights)
    cassini = u + v - z
    require(cassini == (-1) ** k, "Cassini matching orientation changed")
    expected_triple = (
        fibonacci[k - 1] * fibonacci[k + 2],
        2 * fibonacci[k] * fibonacci[k + 1],
        fibonacci[2 * k + 1],
    )
    require(euclid_triple(m, n) == expected_triple, "closed Fibonacci triple formula changed")
    expected_content = 2 if k % 3 == 1 else 1
    require(gcd3(*expected_triple) == expected_content, "Fibonacci mod-three content lane changed")

    if k < FIBONACCI_MAX_K:
        next_determinant = m * fibonacci[k + 2] - n * fibonacci[k + 1]
        require(next_determinant == (-1) ** (k + 1), "golden Farey orientation changed")
    if k >= 3:
        w01, w02, w03, w12, w13, w23 = weights
        require(w01 < w02 < min(w03, w12) < max(w03, w12) < w13 < w23, "golden T6 order changed")
        matching_gaps = tuple(abs(weights[left] - weights[right]) for left, right in MATCHING_INDEX_PAIRS)
        require(matching_gaps[2] == 1 and min(matching_gaps[:2]) > 1, "Cassini matching was not the unique unit gap")
        require(len(set(values)) == 4 and len(set(weights)) == 6, "tie entered the golden tournament range")
    golden_rows.append((k, values, (u, v, z), cassini, expected_triple, expected_content))

for k in range(3, FIBONACCI_MAX_K):
    previous = window(fibonacci[k - 1], fibonacci[k])
    current = window(fibonacci[k], fibonacci[k + 1])
    following = window(fibonacci[k + 1], fibonacci[k + 2])
    require(following == tuple(a + b for a, b in zip(previous, current)), "golden Stern--Brocot mediant recurrence changed")

UNIT_BOX = 600
unit_pairs = {
    (x, y)
    for x in range(1, UNIT_BOX + 1)
    for y in range(1, UNIT_BOX + 1)
    if gcd(x, y) == 1 and abs(x * x + x * y - y * y) == 1
}
expected_unit_pairs = set()
for k in range(2, len(fibonacci)):
    pair = (fibonacci[k - 1], fibonacci[k])
    if max(pair) > UNIT_BOX:
        break
    expected_unit_pairs.add(pair)
require(unit_pairs == expected_unit_pairs, "bounded hostile scan found a non-Fibonacci unit-Cassini seed")

# The U-spine is the first-coordinate-one parabolic lane and the telescoping
# primitive harmonic family (t,t+1,t(t+1)).
U_SPINE_BOUND = 240
u_seed = (1, 1)
for t in range(1, U_SPINE_BOUND + 1):
    if t > 1:
        u_seed = branch(u_seed, "A")
    require(u_seed == (1, t), "U-spine branch index changed")
    m, n = fraction_from_seed(u_seed)
    require((m, n) == (t, t + 1), "U-spine fraction changed")
    require(window(m, n) == (1, t, t + 1, 2 * t + 1), "U-spine window changed")
    require(harmonic_face(edge_weights(window(m, n))) == (t, t + 1, t * (t + 1)), "U-spine harmonic triple changed")
    require(m - n == -1, "U-spine left the fixed-cusp Farey fan")

# THM-3506's exact/conditional odd face orbit becomes an even-x harmonic
# face orbit.  The first two entries are full verified packets; (271,99) is
# the exact next exposed pair, while further iteration requires renewal.
KELLER_FACE_PAIRS = ((7, 3), (43, 15), (271, 99))
KELLER_PRIMITIVE_TRIPLES = ((20, 21, 29), (812, 645, 1037), (31820, 26829, 41621))
keller_seeds = []
for (e, m), expected_triple in zip(KELLER_FACE_PAIRS, KELLER_PRIMITIVE_TRIPLES):
    require(e % 2 == m % 2 == 1 and gcd(e, m) == 1, "THM-3506 pair lost its odd coprime sidecar")
    x, y = e - m, m
    require(x % 2 == 0 and y % 2 == 1 and gcd(x, y) == 1, "THM-3506 pair left the even-x tree")
    face = (m * (e - m), e * (e - m), e * m)
    require(face == harmonic_face(edge_weights((x, y, x + y, x + 2 * y))), "THM-3506 harmonic face chart changed")
    raw_triple = triple_from_face(face)
    require(raw_triple == (e * e - m * m, 2 * e * m, e * e + m * m), "THM-3506 raw current chart changed")
    require(gcd3(*raw_triple) == 2, "THM-3506 raw current left content two")
    require(tuple(value // 2 for value in raw_triple) == expected_triple, "THM-3506 primitive triple changed")
    keller_seeds.append((x, y))

require(keller_seed_update(keller_seeds[0]) == keller_seeds[1], "first verified Keller seed update changed")
require(keller_seed_update(keller_seeds[1]) == keller_seeds[2], "second exact Keller seed update changed")
require(keller_seed_update(keller_seeds[2]) == (1084, 615), "next conditional Keller seed changed")

source_root, source_word = descend(keller_seeds[0])
target_root, target_word = descend(keller_seeds[1])
require(source_root == target_root == (2, 1), "Keller sidecar left the even Berggren root")
require(source_word == ("B",), "first Keller ancestry word changed")
require(target_word == ("A",) * 6 + ("B",), "second Keller ancestry word changed")
require(target_word[:len(source_word)] != source_word, "Keller update became a descendant branch word")

branch_matrices = {
    "A": ((1, 0), (1, 1)),
    "B": ((1, 2), (1, 1)),
    "C": ((1, 2), (0, 1)),
}
keller_matrix = ((4, 4), (3, 1))


def det2(matrix: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


require(tuple(det2(branch_matrices[letter]) for letter in "ABC") == (1, -1, 1), "branch determinant signs changed")
require(det2(keller_matrix) == -8, "Keller seed determinant changed")
for branch_determinant in (1, -1):
    scalar_square = Fraction(det2(keller_matrix), branch_determinant)
    require(not is_rational_square(scalar_square), "Keller matrix became projectively unimodular over Q")

# Tournament and quotient hostiles.
root_window = window(1, 2)
root_edges = edge_weights(root_window)
require(len(set(root_window)) < 4 and len(set(root_edges)) < 6, "tied root became a tournament")

golden_3 = window(fibonacci[3], fibonacci[4])
golden_4 = window(fibonacci[4], fibonacci[5])
require(golden_3 == (1, 2, 3, 5) and golden_4 == (2, 3, 5, 8), "orientation hostile windows changed")
require(harmonic_face(edge_weights(golden_3))[0] + harmonic_face(edge_weights(golden_3))[1] - harmonic_face(edge_weights(golden_3))[2] == -1, "negative Cassini hostile changed")
require(harmonic_face(edge_weights(golden_4))[0] + harmonic_face(edge_weights(golden_4))[1] - harmonic_face(edge_weights(golden_4))[2] == 1, "positive Cassini hostile changed")

induced_edge_permutations = set()
for permutation in permutations(VERTICES):
    image = []
    for left, right in EDGES:
        moved = tuple(sorted((permutation[left], permutation[right])))
        image.append(EDGE_INDEX[moved])
    induced_edge_permutations.add(tuple(image))
require(len(induced_edge_permutations) == 24, "K4 edge action lost faithfulness")
isolated_swap = list(range(6))
isolated_swap[EDGE_INDEX[(0, 3)]], isolated_swap[EDGE_INDEX[(1, 2)]] = (
    isolated_swap[EDGE_INDEX[(1, 2)]],
    isolated_swap[EDGE_INDEX[(0, 3)]],
)
require(tuple(isolated_swap) not in induced_edge_permutations, "isolated Cassini swap became an S4 relabeling")

# Lost-ancestry, parity, word-order, and nonlinear-mediant controls.
half_triple = euclid_triple(1, 2)
third_triple = euclid_triple(1, 3)
require(primitive_unordered(half_triple) == primitive_unordered(third_triple) == (3, 4, 5), "1/2 versus 1/3 ancestry hostile changed")
require(gcd3(*half_triple) == 1 and gcd3(*third_triple) == 2, "1/2 versus 1/3 parity hostile changed")

root_seed = (1, 1)
seed_ab = branch(branch(root_seed, "A"), "B")
seed_ba = branch(branch(root_seed, "B"), "A")
require(fraction_from_seed(seed_ab) == (3, 8), "AB hostile changed")
require(fraction_from_seed(seed_ba) == (5, 8), "BA hostile changed")
require(seed_ab != seed_ba, "branch abelianization recovered ancestry")

left_window = window(1, 2)
right_window = window(2, 3)
mediant_window = window(3, 5)
require(mediant_window == tuple(a + b for a, b in zip(left_window, right_window)), "golden mediant window hostile changed")
require(edge_weights(mediant_window)[0] != edge_weights(left_window)[0] + edge_weights(right_window)[0], "quadratic edge map became additive")

semantic_rows.extend(("golden", *row) for row in golden_rows)
semantic_rows.extend(
    [
        ("counts", fraction_count, odd_lane_count, even_lane_count, max_descent_depth),
        ("berggren_levels", tuple(berggren_level_counts)),
        ("stern_brocot", len(stern_brocot_nodes)),
        ("farey_edges", farey_edge_count),
        ("harmonic_converse", primitive_harmonic_count),
        ("unit_pairs", tuple(sorted(unit_pairs))),
        ("keller_sidecar", KELLER_FACE_PAIRS, tuple(keller_seeds), KELLER_PRIMITIVE_TRIPLES, source_word, target_word),
        ("hostiles", half_triple, third_triple, seed_ab, seed_ba, tuple(isolated_swap)),
    ]
)
semantic_digest = hashlib.sha256(repr(semantic_rows).encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "PIN_AFTER_FIRST_RUN":
    require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic ledger changed")

print("== THM-3509 reduced-fraction harmonic K4/Fibonacci probe ==")
print(f"dependency LF hashes={len(DEPENDENCY_HASHES)}/{len(DEPENDENCY_HASHES)}")
print(f"related LF hashes={len(RELATED_HASHES)}/{len(RELATED_HASHES)}")
print(
    f"reduced fractions 0<m<n<={FRACTION_DENOMINATOR_BOUND}: "
    f"total={fraction_count}; primitive-raw odd-x lane={odd_lane_count}; "
    f"content-two even-x lane={even_lane_count}; max descent={max_descent_depth}"
)
print(f"two Berggren trees through depth {BERGGREN_DEPTH}: level counts={berggren_level_counts}; nodes={len(all_berggren_nodes)}")
print(f"Stern-Brocot depth {STERN_BROCOT_DEPTH}: nodes={len(stern_brocot_nodes)}; componentwise window mediants exact")
print(f"Farey denominator box {FAREY_DENOMINATOR_BOUND}: edges={farey_edge_count}; edge polarization exact")
print(f"primitive harmonic converse with v,z<={HARMONIC_EDGE_BOUND}: triples={primitive_harmonic_count}")
print("harmonic face: 1/u=1/v+1/z; octahedral completion=(u+v,2z-u,v+2z)")
print("Pythagorean current: (a,b,c)=(u+v,2z,v+2z-u); raw content=1 iff gcd(u,v) is odd")
print(f"Fibonacci rows k=2..{FIBONACCI_MAX_K}: unit Cassini, closed triples, mod-3 content, T6 boundary exact")
print(f"unit-Cassini hostile box x,y<={UNIT_BOX}: solutions={len(unit_pairs)}; all consecutive Fibonacci seeds")
print(f"U-spine t=1..{U_SPINE_BOUND}: face=(t,t+1,t(t+1)); fixed-cusp determinant=-1")
print("THM-3506 sidecar: (7,3),(43,15),(271,99) -> even-x harmonic faces and halved primitive triples")
print("THM-3506 update: N=((4,4),(3,1)), det=-8; branch determinants are +/-1, with no rational projective repair")
print("hostile ancestry/parity: 1/2 -> (3,4,5), 1/3 -> (8,6,10)/2 -> leg-swapped (3,4,5)")
print("hostile word order: chronological AB -> 3/8, BA -> 5/8")
print("hostile orientation: (1,2,3,5) has Cassini -1; (2,3,5,8) has +1; same transitive T4 order")
print("hostile toric lift: six weights 2 satisfy antipodal products but require irrational vertex square root")
print("hostile operation: Stern-Brocot mediants add K4 vertices; six edge products require mixed polarization")
print(f"semantic ledger sha256={semantic_digest}")
print("scope: arithmetic/representation theorem only; no physical LRC current, row exclusion, or JC transfer")
print("all exact checks passed")
