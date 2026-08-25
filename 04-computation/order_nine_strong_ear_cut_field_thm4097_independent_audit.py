#!/usr/bin/env python3
"""Independent exact audit for the THM-4097 strong-ear spectrum.

This referee does not import the primary THM-4097 program.  It imports only
the older validated canonical-augmentation universe through order eight,
then independently implements:

* strong connectivity and Held--Karp Hamiltonian-path counting;
* subset dynamic programming for path starts, path ends, and the one-slot
  repair matrix Q;
* the complement-even cut metric and complement-odd vertex field; and
* the complete image of every nonconstant ear from every strong order-eight
  isomorphism representative, retained separately by cut cardinality.

For a parent T, let s_i/e_i count Hamiltonian paths starting/ending at i and
let Q_ij count vertex permutations with i immediately followed by j and all
other adjacent steps valid.  Put

    w_ij = (Q_ij + Q_ji)/2,
    h_i  = s_i - e_i + (sum_j Q_ji - sum_j Q_ij)/2.

If a new vertex x beats exactly the cut S, the audited identity is

    H(T+x_S) = H(T) + cut_w(S) + sum_(i in S) h_i.

The older generator supplies one representative of every order-eight
isomorphism class.  Moon's 1966 induced strong-subtournament theorem is the
external theorem that turns the computed ear image into the full strong
order-nine support; this script audits the finite ear image, not that cited
theorem.

Finite universe: all 6,880 order-eight representatives, all 6,008 strong
representatives, and all 254 nonconstant cuts of each strong parent, for
1,526,032 exact parent/cut rows.  Optimization-safe ``require`` gates remain
active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import struct
import sys
from itertools import combinations, permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
GENERATOR_PATH = HERE / "strong_H_spectrum_m8_isoclass_monad_s5.py"
GENERATOR_SPEC = importlib.util.spec_from_file_location(
    "thm4097_independent_order8_generator", GENERATOR_PATH
)
if GENERATOR_SPEC is None or GENERATOR_SPEC.loader is None:
    raise RuntimeError(f"cannot load validated order-eight generator: {GENERATOR_PATH}")
order8_generator = importlib.util.module_from_spec(GENERATOR_SPEC)
GENERATOR_SPEC.loader.exec_module(order8_generator)
sys.stdout.reconfigure(newline="\n")


def require(condition: bool, message: str) -> None:
    """Optimization-safe assertion."""
    if not condition:
        raise RuntimeError(message)


def has_arc(adjacency: list[int], tail: int, tip: int) -> bool:
    return bool((adjacency[tail] >> tip) & 1)


def require_tournament(adjacency: list[int]) -> None:
    n = len(adjacency)
    for vertex, row in enumerate(adjacency):
        require(((row >> vertex) & 1) == 0, f"loop at vertex {vertex}")
        require(row >> n == 0, f"out-of-range adjacency bit at vertex {vertex}")
    for left, right in combinations(range(n), 2):
        require(
            has_arc(adjacency, left, right) != has_arc(adjacency, right, left),
            f"pair {left},{right} is not oriented exactly once",
        )


def reverse_rows(adjacency: list[int]) -> list[int]:
    n = len(adjacency)
    reverse = [0] * n
    for tail in range(n):
        row = adjacency[tail]
        while row:
            bit = row & -row
            tip = bit.bit_length() - 1
            reverse[tip] |= 1 << tail
            row -= bit
    return reverse


def reaches_all(adjacency: list[int], root: int) -> bool:
    full = (1 << len(adjacency)) - 1
    seen = 1 << root
    frontier = seen
    while frontier:
        next_frontier = 0
        active = frontier
        while active:
            bit = active & -active
            vertex = bit.bit_length() - 1
            next_frontier |= adjacency[vertex]
            active -= bit
        next_frontier &= ~seen
        seen |= next_frontier
        frontier = next_frontier
    return seen == full


def is_strong(adjacency: list[int]) -> bool:
    return reaches_all(adjacency, 0) and reaches_all(reverse_rows(adjacency), 0)


def hamiltonian_path_count(adjacency: list[int]) -> int:
    """Independent Held--Karp count."""
    n = len(adjacency)
    size = 1 << n
    endings = [[0] * n for _ in range(size)]
    for vertex in range(n):
        endings[1 << vertex][vertex] = 1
    for mask in range(size):
        for tail, count in enumerate(endings[mask]):
            if count == 0:
                continue
            available = adjacency[tail] & ~mask
            while available:
                bit = available & -available
                tip = bit.bit_length() - 1
                endings[mask | bit][tip] += count
                available -= bit
    return sum(endings[-1])


def add_vertex(adjacency: list[int], signature: int) -> list[int]:
    """Add x, with signature bit i equal to one exactly when x -> i."""
    n = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(n):
        if (signature >> vertex) & 1:
            child[n] |= 1 << vertex
        else:
            child[vertex] |= 1 << n
    return child


def flip_arc(adjacency: list[int], left: int, right: int) -> list[int]:
    flipped = list(adjacency)
    flipped[left] ^= 1 << right
    flipped[right] ^= 1 << left
    return flipped


def boundary_state(adjacency: list[int]) -> tuple[list[int], list[int], list[list[int]]]:
    """Subset-DP starts, ends, and one-slot repair counts Q."""
    n = len(adjacency)
    size = 1 << n
    full = size - 1
    ends_dp = [[0] * n for _ in range(size)]
    starts_dp = [[0] * n for _ in range(size)]
    for vertex in range(n):
        ends_dp[1 << vertex][vertex] = 1
        starts_dp[1 << vertex][vertex] = 1

    for mask in range(1, size):
        for vertex in range(n):
            vertex_bit = 1 << vertex
            if not (mask & vertex_bit):
                continue
            previous = mask ^ vertex_bit
            if previous == 0:
                continue

            ending_total = 0
            candidates = previous
            while candidates:
                bit = candidates & -candidates
                predecessor = bit.bit_length() - 1
                if adjacency[predecessor] & vertex_bit:
                    ending_total += ends_dp[previous][predecessor]
                candidates -= bit
            ends_dp[mask][vertex] = ending_total

            starting_total = 0
            candidates = previous
            while candidates:
                bit = candidates & -candidates
                successor = bit.bit_length() - 1
                if adjacency[vertex] & bit:
                    starting_total += starts_dp[previous][successor]
                candidates -= bit
            starts_dp[mask][vertex] = starting_total

    # If U contains a but not b, concatenate a path on U ending at a with a
    # path on V\U starting at b.  The a,b slot is the sole unconstrained slot.
    repair = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        complement = full ^ mask
        left_vertices = mask
        while left_vertices:
            left_bit = left_vertices & -left_vertices
            left = left_bit.bit_length() - 1
            left_vertices -= left_bit
            left_count = ends_dp[mask][left]
            if left_count == 0:
                continue
            right_vertices = complement
            while right_vertices:
                right_bit = right_vertices & -right_vertices
                right = right_bit.bit_length() - 1
                right_vertices -= right_bit
                right_count = starts_dp[complement][right]
                if right_count:
                    repair[left][right] += left_count * right_count

    return starts_dp[full], ends_dp[full], repair


def brute_boundary_state(
    adjacency: list[int],
) -> tuple[list[int], list[int], list[list[int]]]:
    """Permutation implementation used only on hostile controls."""
    n = len(adjacency)
    starts = [0] * n
    ends = [0] * n
    repair = [[0] * n for _ in range(n)]
    for order in permutations(range(n)):
        invalid = [
            index
            for index in range(n - 1)
            if not has_arc(adjacency, order[index], order[index + 1])
        ]
        if not invalid:
            starts[order[0]] += 1
            ends[order[-1]] += 1
        if len(invalid) <= 1:
            for index in range(n - 1):
                if not invalid or invalid == [index]:
                    repair[order[index]][order[index + 1]] += 1
    return starts, ends, repair


def cut_field(
    starts: list[int], ends: list[int], repair: list[list[int]]
) -> tuple[list[list[int]], list[int]]:
    n = len(starts)
    weights = [[0] * n for _ in range(n)]
    field: list[int] = []
    for left, right in combinations(range(n), 2):
        symmetric = repair[left][right] + repair[right][left]
        require(symmetric % 2 == 0, f"nonintegral cut weight at {left},{right}")
        weights[left][right] = weights[right][left] = symmetric // 2
    for vertex in range(n):
        column = sum(repair[other][vertex] for other in range(n))
        row = sum(repair[vertex])
        require((column - row) % 2 == 0, f"nonintegral field at {vertex}")
        field.append(starts[vertex] - ends[vertex] + (column - row) // 2)
    require(sum(field) == 0, "vertex field is not zero-sum")
    return weights, field


def insertion_formula(
    starts: list[int], ends: list[int], repair: list[list[int]], signature: int
) -> int:
    n = len(starts)
    total = 0
    for vertex in range(n):
        if (signature >> vertex) & 1:
            total += starts[vertex]
        else:
            total += ends[vertex]
            for beaten in range(n):
                if (signature >> beaten) & 1:
                    total += repair[vertex][beaten]
    return total


def cut_field_formula(
    parent_h: int,
    weights: list[list[int]],
    field: list[int],
    signature: int,
) -> tuple[int, int, int]:
    n = len(field)
    cut = sum(
        weights[left][right]
        for left, right in combinations(range(n), 2)
        if ((signature >> left) & 1) != ((signature >> right) & 1)
    )
    charge = sum(field[vertex] for vertex in range(n) if (signature >> vertex) & 1)
    return parent_h + cut + charge, cut, charge


def directed_triangle_count(adjacency: list[int]) -> int:
    count = 0
    for a, b, c in combinations(range(len(adjacency)), 3):
        first = has_arc(adjacency, a, b) and has_arc(adjacency, b, c) and has_arc(adjacency, c, a)
        second = has_arc(adjacency, a, c) and has_arc(adjacency, c, b) and has_arc(adjacency, b, a)
        count += int(first or second)
    return count


def odd_runs(values: set[int]) -> list[tuple[int, int, int]]:
    ordered = sorted(values)
    require(bool(ordered), "empty ear spectrum")
    runs: list[tuple[int, int, int]] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 2:
            previous = value
            continue
        runs.append((start, previous, (previous - start) // 2 + 1))
        start = previous = value
    runs.append((start, previous, (previous - start) // 2 + 1))
    return runs


def semantic_digest_update(
    digest: "hashlib._Hash", representative_index: int, adjacency: list[int], rows: list[int]
) -> None:
    """Stable binary schema: P|rep-index|8 rows, then R|sig|H for sig=1..254."""
    digest.update(b"P")
    digest.update(struct.pack("<I8B", representative_index, *adjacency))
    for signature, value in enumerate(rows, start=1):
        digest.update(b"R")
        digest.update(struct.pack("<HI", signature, value))


def carrier_digest_update(
    digest: "hashlib._Hash",
    representative_index: int,
    starts: list[int],
    ends: list[int],
    repair: list[list[int]],
    weights: list[list[int]],
    field: list[int],
) -> None:
    """Stable signed-64-bit digest of the complete parent boundary carrier."""
    digest.update(struct.pack("<I", representative_index))
    payload = starts + ends + [entry for row in repair for entry in row]
    payload += [weights[left][right] for left, right in combinations(range(8), 2)]
    payload += field
    for value in payload:
        digest.update(struct.pack("<q", value))


def main() -> None:
    expected_counts = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    representatives, counts = order8_generator.generate(8)
    require(counts == expected_counts, f"canonical universe mismatch: {counts}")

    spectrum: set[int] = set()
    by_weight = {weight: set() for weight in range(1, 8)}
    target_multiplicity = {613: 0, 623: 0}
    target_parent_support = {613: set(), 623: set()}
    strong_parents = 0
    ear_rows = 0
    pair_integrality_gates = 0
    field_integrality_gates = 0
    direct_child_controls = 0
    semantic_digest = hashlib.sha256()
    carrier_digest = hashlib.sha256()

    direct_order8_613: list[tuple[int, tuple[int, ...]]] = []

    for representative_index, adjacency in enumerate(representatives[8]):
        require_tournament(adjacency)
        if not is_strong(adjacency):
            continue
        strong_parents += 1
        parent_h = hamiltonian_path_count(adjacency)
        if parent_h == 613:
            direct_order8_613.append((representative_index, tuple(adjacency)))

        starts, ends, repair = boundary_state(adjacency)
        require(sum(starts) == parent_h, f"start total mismatch at parent {representative_index}")
        require(sum(ends) == parent_h, f"end total mismatch at parent {representative_index}")
        weights, field = cut_field(starts, ends, repair)
        pair_integrality_gates += 28
        field_integrality_gates += 8
        carrier_digest_update(
            carrier_digest,
            representative_index,
            starts,
            ends,
            repair,
            weights,
            field,
        )

        row_values: list[int] = []
        for signature in range(1, 255):
            via_incidence = insertion_formula(starts, ends, repair, signature)
            via_cut, _, _ = cut_field_formula(parent_h, weights, field, signature)
            require(
                via_incidence == via_cut,
                f"cut-field identity failed at parent={representative_index}, sig={signature}",
            )
            require(via_cut % 2 == 1, "Rédei parity failed in ear image")
            require(via_cut >= parent_h, "vertex-insertion monotonicity failed")
            row_values.append(via_cut)
            spectrum.add(via_cut)
            by_weight[signature.bit_count()].add(via_cut)
            ear_rows += 1
            if via_cut in target_multiplicity:
                target_multiplicity[via_cut] += 1
                target_parent_support[via_cut].add(representative_index)

        # One direct nine-vertex Held--Karp control per strong parent.  The
        # signature schedule is deterministic and always nonconstant.
        control_signature = representative_index % 254 + 1
        control_child = add_vertex(adjacency, control_signature)
        require(is_strong(control_child), "nonconstant ear lost strongness")
        require(
            hamiltonian_path_count(control_child) == row_values[control_signature - 1],
            f"direct child count failed at parent={representative_index}",
        )
        direct_child_controls += 1
        semantic_digest_update(semantic_digest, representative_index, adjacency, row_values)

    require(strong_parents == 6008, f"strong parent count changed: {strong_parents}")
    require(ear_rows == 1_526_032, f"ear universe changed: {ear_rows}")
    require(pair_integrality_gates == 168_224, "pair integrality universe changed")
    require(field_integrality_gates == 48_064, "field integrality universe changed")
    require(direct_child_controls == 6008, "direct child control universe changed")

    require(len(spectrum) == 1482, f"spectrum support changed: {len(spectrum)}")
    require(min(spectrum) == 75 and max(spectrum) == 3357, "spectrum endpoints changed")
    runs = odd_runs(spectrum)
    longest = max(runs, key=lambda row: row[2])
    require(longest == (85, 2881, 1399), f"solid interval changed: {longest}")
    low_holes = [value for value in range(75, 85, 2) if value not in spectrum]
    require(low_holes == [77, 79, 83], f"low holes changed: {low_holes}")
    require(target_multiplicity == {613: 674, 623: 780}, "target row multiplicities changed")

    expected_weight_sizes = {1: 962, 2: 1294, 3: 1467, 4: 1478, 5: 1467, 6: 1294, 7: 962}
    for weight in range(1, 8):
        require(
            len(by_weight[weight]) == expected_weight_sizes[weight],
            f"cut-cardinality support size changed at weight {weight}",
        )
        require(
            by_weight[weight] == by_weight[8 - weight],
            f"converse cut-cardinality symmetry failed at weight {weight}",
        )
        require(by_weight[weight] != spectrum, f"weight {weight} became an exact singleton basis")
    weight_four_misses = spectrum - by_weight[4]
    require(weight_four_misses == {89, 93, 105, 125}, "weight-four boundary changed")
    require(weight_four_misses <= by_weight[3], "weight-three repairs changed")
    require(by_weight[3] | by_weight[4] == spectrum, "central two-weight basis changed")

    spectrum_weight_payload = "ALL " + ",".join(map(str, sorted(spectrum))) + "\n"
    spectrum_weight_payload += "\n".join(
        f"W{weight} " + ",".join(map(str, sorted(by_weight[weight])))
        for weight in range(1, 8)
    )
    spectrum_weight_digest = hashlib.sha256(spectrum_weight_payload.encode("utf-8")).hexdigest()
    require(
        spectrum_weight_digest
        == "de2b67a2c9bc0a349d33f7c9e996508a53116fbfe1e4764edc2e918005acd736",
        "primary spectrum/weight semantic digest changed",
    )

    expected_order8_613 = [
        (4546, (30, 28, 184, 240, 96, 195, 135, 19)),
        (4572, (158, 28, 56, 240, 224, 67, 7, 102)),
    ]
    require(direct_order8_613 == expected_order8_613, "direct H=613 witnesses changed")

    # Full carrier hostile and displayed 623/735 complementary pair.
    parent_91 = [126, 60, 248, 112, 32, 64, 146, 59]
    require_tournament(parent_91)
    require(is_strong(parent_91), "displayed H=91 parent is not strong")
    require(hamiltonian_path_count(parent_91) == 91, "displayed parent H changed")
    starts, ends, repair = boundary_state(parent_91)
    require(
        (starts, ends, repair) == brute_boundary_state(parent_91),
        "subset and permutation boundary states disagree",
    )
    weights, field = cut_field(starts, ends, repair)
    require(field == [288, 84, 121, -50, -282, -263, -31, 133], "displayed field changed")
    signature_623 = 105
    signature_735 = 255 ^ signature_623
    value_623, cut_623, charge_623 = cut_field_formula(91, weights, field, signature_623)
    value_735, cut_735, charge_735 = cut_field_formula(91, weights, field, signature_735)
    require((cut_623, charge_623, value_623) == (588, -56, 623), "623 carrier changed")
    require((cut_735, charge_735, value_735) == (588, 56, 735), "735 carrier changed")
    child_623 = add_vertex(parent_91, signature_623)
    expected_child_623 = [126, 316, 504, 112, 288, 64, 146, 315, 105]
    require(child_623 == expected_child_623, "displayed child adjacency changed")
    require(is_strong(child_623), "displayed H=623 child is not strong")
    require(hamiltonian_path_count(child_623) == 623, "displayed child H changed")
    child_735 = add_vertex(parent_91, signature_735)
    require(is_strong(child_735), "complementary H=735 child is not strong")
    require(hamiltonian_path_count(child_735) == 735, "complementary child H changed")

    # Arc-flip parity explains w-integrality on the displayed carrier.
    for left, right in combinations(range(8), 2):
        flipped_h = hamiltonian_path_count(flip_arc(parent_91, left, right))
        if has_arc(parent_91, left, right):
            created_minus_destroyed = repair[right][left] - repair[left][right]
        else:
            created_minus_destroyed = repair[left][right] - repair[right][left]
        require(flipped_h - 91 == created_minus_destroyed, "arc-flip response identity failed")
        require(created_minus_destroyed % 2 == 0, "arc-flip Rédei parity failed")

    # THM-4093 firewall in both directions at p=7.
    require(directed_triangle_count(child_623) == 16, "H=623 triangle count changed")
    reverse_firewall = [14, 60, 24, 48, 33, 5]
    require_tournament(reverse_firewall)
    require(is_strong(reverse_firewall), "reverse firewall is not strong")
    require(hamiltonian_path_count(reverse_firewall) == 33, "reverse firewall H changed")
    require(directed_triangle_count(reverse_firewall) == 7, "reverse firewall c3 changed")

    spectrum_digest = hashlib.sha256(
        ",".join(str(value) for value in sorted(spectrum)).encode("ascii")
    ).hexdigest()

    print("THM-4097 independent strong-ear cut-field audit: PASS")
    print(
        "finite universe: order8 iso reps=6880; strong parents=6008; "
        "nonconstant cuts/parent=254; ear rows=1526032"
    )
    print(
        "carrier gates: pair-integrality=168224; field-integrality=48064; "
        "cut-field identities=1526032; direct child DP controls=6008"
    )
    print("order9 strong-ear support: count=1482; min=75; max=3357")
    print("longest solid odd interval: [85,2881], count=1399; low holes=[77,79,83]")
    print(
        "cut-cardinality spectra E1..E7 sizes: 962,1294,1467,1478,1467,1294,962; "
        "converse symmetry PASS"
    )
    print("central basis: E4 misses exactly [89,93,105,125]; E3 repairs all; E3 union E4 exact")
    print(
        "target parent/cut multiplicities: H613=674 across "
        f"{len(target_parent_support[613])} parents; H623=780 across "
        f"{len(target_parent_support[623])} parents"
    )
    print("direct strong order8 H613 reps: 4546,4572")
    print("displayed complement pair: Hparent=91; cut=588; field=-56/+56; H=623/735")
    print("c3 firewalls at p=7: (H,c3)=(623,16) and (33,7)")
    print(f"semantic row digest sha256={semantic_digest.hexdigest()}")
    print(f"spectrum/weight semantic digest sha256={spectrum_weight_digest}")
    print(f"boundary carrier digest sha256={carrier_digest.hexdigest()}")
    print(f"sorted support digest sha256={spectrum_digest}")


if __name__ == "__main__":
    main()
