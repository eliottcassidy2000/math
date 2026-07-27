#!/usr/bin/env python3
"""Exact referee for THM-2458.

The script deliberately uses only integer and Fraction arithmetic.  The
18,900 count is the deterministic repeated-root recursion traversal count;
canonicalizing the four unordered two-root masks gives 15,120 distinct
clean covers.  All theorem-level signature tests are performed on sets of
supports, so the traversal multiplicity is immaterial.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import gcd


P = 13
ROOTS = frozenset(range(P))
GATE = frozenset((0, 1))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def ap(start, step, length):
    return frozenset((start + j * step) % P for j in range(length))


def chord_step(pair):
    require(len(pair) == 2, "a chord must have two distinct roots")
    a, b = sorted(pair)
    d = (b - a) % P
    return min(d, P - d)


def repeated_root_pairings(items):
    """The historical deterministic traversal used for the 18,900 census.

    When ``items`` contains one repeated root, identical copies are not
    globally canonicalized.  This is why the traversal count is larger than
    the number of distinct unordered mask families.  Equal choices at one
    recursion level are suppressed, and a self-pair is forbidden.
    """

    items = tuple(sorted(items))
    if not items:
        yield ()
        return
    a = items[0]
    previous_b = None
    for j in range(1, len(items)):
        b = items[j]
        if b == a or b == previous_b:
            continue
        previous_b = b
        rest = items[1:j] + items[j + 1 :]
        for tail in repeated_root_pairings(rest):
            pair = tuple(sorted((a, b)))
            yield tuple(sorted((pair,) + tail))


def cover_key(guard, pairs):
    return (
        tuple(sorted(guard)),
        tuple(sorted(tuple(sorted(pair)) for pair in pairs)),
    )


def enumerate_clean_cover_records():
    """Enumerate the exact fixed-gate clean-cover universe.

    The guard is a four-term AP of unsigned step 1,...,6.  The other four
    dangers are genuine two-root masks disjoint from the gate/guard union.
    Total incidence is fourteen and every root is covered, hence there is
    exactly one double root.
    """

    records = []
    distinct = {}
    overlap_records = 0
    disjoint_records = 0
    for guard_step in range(1, 7):
        for guard_start in range(P):
            guard = ap(guard_start, guard_step, 4)
            base_overlap = GATE & guard
            if len(base_overlap) > 1:
                continue
            complement = sorted(ROOTS - (GATE | guard))
            if len(base_overlap) == 1:
                targets = (tuple(complement),)
            else:
                targets = tuple(tuple(sorted(complement + [extra])) for extra in complement)
            for target in targets:
                for raw_pairs in repeated_root_pairings(target):
                    pairs = tuple(frozenset(pair) for pair in raw_pairs)
                    masks = (GATE, guard) + pairs
                    counts = [sum(root in mask for mask in masks) for root in ROOTS]
                    require(sorted(counts) == [1] * 12 + [2], "not an excess-one cover")
                    require(set().union(*masks) == set(ROOTS), "cover misses a root")
                    require(
                        all(pair.isdisjoint(GATE | guard) for pair in pairs),
                        "low pair meets the gate/guard base",
                    )
                    require(all(len(pair) == 2 for pair in pairs), "degenerate pair")
                    records.append((guard_step, guard, pairs))
                    distinct.setdefault(
                        cover_key(guard, pairs),
                        (guard_step, guard, pairs),
                    )
                    if base_overlap:
                        overlap_records += 1
                    else:
                        disjoint_records += 1
    return records, tuple(distinct.values()), overlap_records, disjoint_records


def atom_support(qstar, remaining_roles, bits):
    support = set(ROOTS - qstar)
    for role, danger in zip(remaining_roles, bits):
        support &= set(role if danger else ROOTS - role)
    return frozenset(support)


def build_signature_atlases(records):
    guard_atlas = defaultdict(set)
    ordinary_atlas = defaultdict(set)
    zero_atom_checks = 0
    two_atom_checks = 0
    high_atom_checks = 0

    for guard_step, guard, pairs in records:
        pair_steps = tuple(chord_step(pair) for pair in pairs)
        for q_index, qstar in enumerate(pairs):
            ordinary_indices = tuple(i for i in range(4) if i != q_index)
            ordinary_roles = tuple(pairs[i] for i in ordinary_indices)
            remaining_roles = (guard,) + ordinary_roles

            for bits in product((0, 1), repeat=4):
                support = atom_support(qstar, remaining_roles, bits)
                danger_count = sum(bits)
                if danger_count == 0:
                    require(support <= GATE, "zero-danger support escapes the factored gate")
                    zero_atom_checks += 1
                elif danger_count >= 3:
                    require(not support, "three-danger clean atom is nonempty")
                    high_atom_checks += 1
                elif danger_count == 2:
                    require(len(support) <= 1, "two-danger atom has more than one root")
                    require(support.isdisjoint(GATE), "two-danger atom meets the gate")
                    two_atom_checks += 1

            guard_support = atom_support(qstar, remaining_roles, (1, 0, 0, 0))
            guard_signature = (
                pair_steps[q_index],
                guard_step,
            ) + tuple(sorted(pair_steps[i] for i in ordinary_indices))
            guard_atlas[guard_signature].add(guard_support)

            for local_index, danger_index in enumerate(ordinary_indices):
                bits = [0, 0, 0, 0]
                bits[local_index + 1] = 1
                support = atom_support(qstar, remaining_roles, tuple(bits))
                safe_indices = tuple(i for i in ordinary_indices if i != danger_index)
                signature = (
                    pair_steps[q_index],
                    guard_step,
                    pair_steps[danger_index],
                ) + tuple(sorted(pair_steps[i] for i in safe_indices))
                ordinary_atlas[signature].add(support)

    return (
        guard_atlas,
        ordinary_atlas,
        zero_atom_checks,
        two_atom_checks,
        high_atom_checks,
    )


def ordinary_one_danger_size_census(distinct_records):
    """Count labelled ordinary-danger atoms on canonical mask families."""

    census = {1: 0, 2: 0}
    for _, guard, pairs in distinct_records:
        for q_index, qstar in enumerate(pairs):
            ordinary_indices = tuple(i for i in range(4) if i != q_index)
            ordinary_roles = tuple(pairs[i] for i in ordinary_indices)
            remaining_roles = (guard,) + ordinary_roles
            for local_index in range(3):
                bits = [0, 0, 0, 0]
                bits[local_index + 1] = 1
                support = atom_support(qstar, remaining_roles, tuple(bits))
                require(
                    len(support) in (1, 2),
                    "ordinary one-danger support is not one or two roots",
                )
                census[len(support)] += 1
    return census


def rational_span_certificate(supports):
    """Return None if 1 lies in the column span, else an exact dual row.

    The returned vector y satisfies y^T 1 != 0 and y^T 1_C = 0 for every
    support C.  Row operations are tracked exactly over Fraction.
    """

    supports = tuple(sorted(supports, key=lambda s: tuple(sorted(s))))
    column_count = len(supports)
    matrix = [
        [Fraction(int(root in supports[j])) for j in range(column_count)]
        for root in range(P)
    ]
    rhs = [Fraction(1) for _ in range(P)]
    transform = [
        [Fraction(int(i == j)) for j in range(P)]
        for i in range(P)
    ]

    row = 0
    for column in range(column_count):
        pivot = next((i for i in range(row, P) if matrix[i][column]), None)
        if pivot is None:
            continue
        if pivot != row:
            matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
            rhs[row], rhs[pivot] = rhs[pivot], rhs[row]
            transform[row], transform[pivot] = transform[pivot], transform[row]
        scale = matrix[row][column]
        matrix[row] = [value / scale for value in matrix[row]]
        rhs[row] /= scale
        transform[row] = [value / scale for value in transform[row]]
        for i in range(P):
            if i == row:
                continue
            scale = matrix[i][column]
            if not scale:
                continue
            matrix[i] = [
                matrix[i][j] - scale * matrix[row][j]
                for j in range(column_count)
            ]
            rhs[i] -= scale * rhs[row]
            transform[i] = [
                transform[i][j] - scale * transform[row][j]
                for j in range(P)
            ]
        row += 1
        if row == P:
            break

    for i in range(P):
        if not any(matrix[i]) and rhs[i]:
            certificate = tuple(transform[i])
            require(sum(certificate) != 0, "dual certificate misses the constant")
            for support in supports:
                require(
                    sum(certificate[root] for root in support) == 0,
                    "dual certificate misses a support",
                )
            return certificate
    return None


def determinant_and_solution(matrix, rhs):
    size = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    vector = [Fraction(value) for value in rhs]
    determinant = Fraction(1)
    for column in range(size):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        require(pivot is not None, "singular square system")
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            vector[column], vector[pivot] = vector[pivot], vector[column]
            determinant = -determinant
        scale = work[column][column]
        determinant *= scale
        work[column] = [value / scale for value in work[column]]
        vector[column] /= scale
        for i in range(size):
            if i == column:
                continue
            scale = work[i][column]
            if not scale:
                continue
            work[i] = [
                work[i][j] - scale * work[column][j]
                for j in range(size)
            ]
            vector[i] -= scale * vector[column]
    return determinant, tuple(vector)


ATLAS_ROWS = (
    (0, 9, 7, 3, 12),
    (1, 2, 9, 5, 7),
    (2, 3, 10, 6, 3),
    (3, 2, 10, 6, 7),
    (4, 10, 7, 2, 11),
    (5, 6, 3, 9, 11),
    (6, 2, 9, 2, 7),
    (7, 8, 2, 2, 6),
    (8, 2, 6, 9, 11),
    (9, 8, 3, 2, 7),
    (10, 9, 4, 3, 3),
    (11, 2, 6, 9, 5),
    (12, 8, 6, 2, 11),
)


def verify_explicit_atlas(guard_atlas):
    signature = (2, 5, 1, 3, 5)
    supports = guard_atlas[signature]
    require(len(supports) == P, "chosen signature does not have thirteen supports")

    incidence = [0] * P
    atlas_supports = set()
    for guard_start, q_start, a_start, b_start, c_start in ATLAS_ROWS:
        guard = ap(guard_start, 5, 4)
        qstar = ap(q_start, 2, 2)
        ordinary_a = ap(a_start, 1, 2)
        ordinary_b = ap(b_start, 3, 2)
        ordinary_c = ap(c_start, 5, 2)
        masks = (GATE, guard, qstar, ordinary_a, ordinary_b, ordinary_c)
        counts = [sum(root in mask for mask in masks) for root in ROOTS]
        require(sorted(counts) == [1] * 12 + [2], "atlas row is not excess-one")
        require(set().union(*masks) == set(ROOTS), "atlas row is not a cover")
        require(
            all(
                role.isdisjoint(GATE | guard)
                for role in (qstar, ordinary_a, ordinary_b, ordinary_c)
            ),
            "guard-danger row has a low role in the gate/guard base",
        )
        require(
            tuple(
                (
                    chord_step(qstar),
                    5,
                    chord_step(ordinary_a),
                    chord_step(ordinary_b),
                    chord_step(ordinary_c),
                )
            )
            == signature,
            "atlas row has the wrong signature",
        )
        atlas_supports.add(guard)
        for root in guard:
            incidence[root] += 1

    require(atlas_supports == supports, "displayed atlas misses a signature support")
    require(incidence == [4] * P, "guard translates are not four-regular")

    matrix = [
        [int(root in ap(start, 5, 4)) for start in range(P)]
        for root in range(P)
    ]
    determinant, solution = determinant_and_solution(matrix, [1] * P)
    require(determinant == 4, "unexpected guard-circulant determinant")
    require(solution == (Fraction(1, 4),) * P, "wrong uniform chart weights")
    return signature, tuple(incidence), determinant, solution


def verify_replica_table():
    shared_support = ap(6, 5, 4)
    require(shared_support == frozenset((3, 6, 8, 11)), "wrong shared chart")
    word = tuple(len(shared_support & ap(shift, 1, 2)) for shift in range(P))
    expected = (0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0)
    require(word == expected, "wrong two-root overlap word")
    table = tuple(
        tuple((2 if source == 0 else 0) + word[target] for target in range(P))
        for source in range(7)
    )
    require(
        tuple(table[source][0] for source in range(7)) == (2, 0, 0, 0, 0, 0, 0),
        "delta anchor fails",
    )
    for source in range(1, 7):
        for target in range(P):
            rectangle = (
                table[source][target]
                - table[source][0]
                - table[0][target]
                + table[0][0]
            )
            require(rectangle == 0, "mixed rectangle survives")

    # Scale every clean chart by 1/13.  The thirteen owner atlas charts have
    # weights 1/52, the seven shared charts have weights 1/13, and an empty
    # chart has weight 11/52.  These are a probability distribution.
    weights = (Fraction(1, 52),) * 13 + (Fraction(1, 13),) * 7 + (
        Fraction(11, 52),
    )
    require(sum(weights) == 1, "Boolean chart mixture is not normalized")
    require(max(Fraction(value, 13) for row in table for value in row) <= 1, "bad scale")
    return shared_support, word, table


def ordinary_phase_code(speed_residue, start):
    """Nearest-integer residues for a two-root danger mask.

    The first entry accompanies delta_u in (-1/2,-1/14); the second
    accompanies delta_u in (1/14,1/2).
    """

    return ((-speed_residue * start) % P, (-speed_residue * start - 1) % P)


def guard_phase_code(speed_residue, start):
    """Nearest-integer residues for a four-root guard-danger mask.

    The entries accompany delta_H in (-1/2,-1/7) and (1/7,1/2).
    """

    return (
        (-speed_residue * start - 1) % P,
        (-speed_residue * start - 2) % P,
    )


def atlas_phase_code(row):
    guard_start, q_start, a_start, b_start, c_start = row
    return (
        ordinary_phase_code(1, 0),
        ordinary_phase_code(7, q_start),
        guard_phase_code(8, guard_start),
        ordinary_phase_code(1, a_start),
        ordinary_phase_code(9, b_start),
        ordinary_phase_code(8, c_start),
    )


def phase_intervals(speed, start, step, threshold, is_guard):
    """Open y-intervals realizing one prescribed oriented root mask."""

    require(speed > 0 and gcd(speed, P) == 1, "bad physical speed")
    require((speed * step) % P == 1, "step is not the inverse speed residue")
    if is_guard:
        cases = (
            (1, Fraction(-1, 2), -threshold),
            (2, threshold, Fraction(1, 2)),
        )
    else:
        cases = (
            (0, Fraction(-1, 2), -threshold),
            (1, threshold, Fraction(1, 2)),
        )
    intervals = []
    for central_offset, delta_left, delta_right in cases:
        residue = (-speed * (start + central_offset * step)) % P
        for nearest_integer in range(residue, speed + 2, P):
            left = max(
                Fraction(0),
                (Fraction(nearest_integer) + delta_left) / speed,
            )
            right = min(
                Fraction(1),
                (Fraction(nearest_integer) + delta_right) / speed,
            )
            if left < right:
                intervals.append((left, right))
    return tuple(intervals)


def intersect_intervals(left_intervals, right_intervals):
    answer = []
    for left_a, right_a in left_intervals:
        for left_b, right_b in right_intervals:
            left = max(left_a, left_b)
            right = min(right_a, right_b)
            if left < right:
                answer.append((left, right))
    return tuple(answer)


def direct_mask(speed, y, threshold):
    mask = set()
    for root in ROOTS:
        value = Fraction(speed) * (y + root) / P
        floor_value = value.numerator // value.denominator
        distance = min(value - floor_value, floor_value + 1 - value)
        if distance < threshold:
            mask.add(root)
    return frozenset(mask)


def validate_phase_interval(speed, start, step, threshold, is_guard, intervals):
    desired = ap(start, step, 4 if is_guard else 2)
    for left, right in intervals:
        midpoint = (left + right) / 2
        require(
            direct_mask(speed, midpoint, threshold) == desired,
            "nearest-integer phase interval has the wrong physical mask",
        )


def bounded_phase_orbit_probe():
    """Exploratory exact search over the first three lifts of each residue.

    This is intentionally not an exhaustive fixed-speed theorem.  It exposes
    the phase-orbit debt in a reproducible finite bank.
    """

    residues = (1, 7, 8, 1, 9, 8)
    candidates = tuple(
        tuple(residue + P * lift for lift in range(3))
        for residue in residues
    )
    start_sets = (
        (0,),
        tuple(sorted({row[1] for row in ATLAS_ROWS})),
        tuple(range(P)),
        tuple(sorted({row[2] for row in ATLAS_ROWS})),
        tuple(sorted({row[3] for row in ATLAS_ROWS})),
        tuple(sorted({row[4] for row in ATLAS_ROWS})),
    )
    steps = (1, 2, 5, 1, 3, 5)
    thresholds = (
        Fraction(1, 14),
        Fraction(1, 14),
        Fraction(1, 7),
        Fraction(1, 14),
        Fraction(1, 14),
        Fraction(1, 14),
    )
    guard_flags = (False, False, True, False, False, False)

    interval_bank = {}
    for role in range(6):
        for speed in candidates[role]:
            for start in start_sets[role]:
                intervals = phase_intervals(
                    speed,
                    start,
                    steps[role],
                    thresholds[role],
                    guard_flags[role],
                )
                validate_phase_interval(
                    speed,
                    start,
                    steps[role],
                    thresholds[role],
                    guard_flags[role],
                    intervals,
                )
                interval_bank[(role, speed, start)] = intervals

    tuple_count = 0
    realized_pairs = 0
    maximum_rows = 0
    for speeds in product(*candidates):
        if speeds[0] == speeds[3] or speeds[2] == speeds[5]:
            continue
        tuple_count += 1
        realized_here = 0
        for row in ATLAS_ROWS:
            starts = (0, row[1], row[0], row[2], row[3], row[4])
            current = ((Fraction(0), Fraction(1)),)
            for role in range(6):
                current = intersect_intervals(
                    current,
                    interval_bank[(role, speeds[role], starts[role])],
                )
                if not current:
                    break
            if current:
                realized_pairs += 1
                realized_here += 1
        maximum_rows = max(maximum_rows, realized_here)

    require(tuple_count == 324, "unexpected bounded lift-bank size")
    require(realized_pairs == 0, "bounded probe unexpectedly realizes an atlas row")
    return tuple_count, realized_pairs, maximum_rows


def main():
    (
        records,
        distinct_records,
        overlap_records,
        disjoint_records,
    ) = enumerate_clean_cover_records()
    require(len(records) == 18900, "wrong clean-cover traversal census")
    require(len(distinct_records) == 15120, "wrong distinct clean-cover census")
    require(overlap_records == 3780, "wrong gate/guard-overlap traversal count")
    require(disjoint_records == 15120, "wrong pair/pair-overlap traversal count")

    (
        guard_atlas,
        ordinary_atlas,
        zero_checks,
        two_checks,
        high_checks,
    ) = build_signature_atlases(records)
    require(len(guard_atlas) == 1949, "wrong guard-danger signature census")
    require(len(ordinary_atlas) == 4430, "wrong ordinary-danger signature census")
    require(max(map(len, guard_atlas.values())) == 13, "wrong guard support maximum")
    require(max(map(len, ordinary_atlas.values())) == 21, "wrong ordinary support maximum")

    ordinary_duals = 0
    for supports in ordinary_atlas.values():
        certificate = rational_span_certificate(supports)
        require(certificate is not None, "ordinary-danger constant enters the span")
        ordinary_duals += 1
    guard_span_feasible = sum(
        rational_span_certificate(supports) is None
        for supports in guard_atlas.values()
    )
    require(ordinary_duals == 4430, "missing ordinary exact dual")
    require(guard_span_feasible == 215, "wrong guard span-feasible census")
    ordinary_sizes = ordinary_one_danger_size_census(distinct_records)
    require(
        ordinary_sizes == {1: 68040, 2: 113400},
        "wrong ordinary one-danger size census",
    )

    signature, incidence, determinant, solution = verify_explicit_atlas(guard_atlas)
    shared_support, word, table = verify_replica_table()
    phase_codes = tuple(atlas_phase_code(row) for row in ATLAS_ROWS)
    require(len(set(phase_codes)) == 13, "phase-code rows collide")
    probe_count, probe_pairs, probe_maximum = bounded_phase_orbit_probe()

    print("THM-2458 CLEAN-ROOT GUARD-DANGER AUDIT")
    print(f"cover_traversal_records={len(records)}")
    print(f"distinct_clean_mask_families={len(distinct_records)}")
    print(
        "traversal_split="
        f"gate_guard_double:{overlap_records},pair_pair_double:{disjoint_records}"
    )
    print(
        "support_checks="
        f"zero:{zero_checks},two:{two_checks},three_or_four:{high_checks};max_size=4"
    )
    print(
        "ordinary_danger="
        f"signatures:{len(ordinary_atlas)},exact_duals:{ordinary_duals},"
        f"span_feasible:0,max_family:{max(map(len, ordinary_atlas.values()))}"
    )
    print(
        "ordinary_one_danger_sizes="
        f"1:{ordinary_sizes[1]},2:{ordinary_sizes[2]}"
    )
    print(
        "guard_danger="
        f"signatures:{len(guard_atlas)},span_feasible:{guard_span_feasible},"
        f"max_family:{max(map(len, guard_atlas.values()))}"
    )
    print(f"hostile_signature={signature}")
    print(f"atlas_rows={len(ATLAS_ROWS)};root_incidence={incidence}")
    print(f"guard_circulant_det={determinant};unique_alpha={solution[0]}")
    print(f"shared_support={tuple(sorted(shared_support))}")
    print(f"overlap_word={word};word_sum={sum(word)}")
    print(f"replica_anchor={tuple(table[source][0] for source in range(7))}")
    print("fixed_residues=(1,7,8,1,9,8);inverse_steps=(1,2,5,1,3,5)")
    for index, code in enumerate(phase_codes):
        print(f"phase_code_{index:02d}={code}")
    print(
        "bounded_orbit_probe="
        f"lift_tuples:{probe_count},realized_row_pairs:{probe_pairs},"
        f"max_rows_per_tuple:{probe_maximum}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
