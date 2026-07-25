#!/usr/bin/env python3
"""Exact branch certificate for the scalar blocker profile (1,2,3).

At N=13^4, a depth-one blocker is an all-or-nothing root-fibre bit
indexed by a unit sign class modulo 13^3, while a depth-two blocker is
indexed modulo 13^2.  For each of the 1,014*78 mixed pairs, we certify the
five largest conditional capacities among all 13,182 unit sign classes.

The branch bound is exact integer arithmetic.  If h_q(r) in {0,1,2} is
the number of q-dangerous guard-safe roots above phase r, and A,B are the
two blocker-active phase sets, then

  C_q(A,B) = F_q-X_q(A)-X_q(B)+X_q(A intersection B)
           <= F_q-X_q(A)-X_q(B)+2*|A intersection B|.

It is also bounded by each one-blocker conditional capacity.  Candidates
are visited in decreasing order of the minimum of these three bounds.
Once the next bound is strictly below the fifth exact capacity, the exact
top-five list, including its label tie-break, is certified.

No Python ``assert`` is used, so all checks remain active under ``-O``.
"""

from collections import Counter, defaultdict
import hashlib
import heapq

import numpy as np


P = 13
N = P**4
Q = P**3


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign_classes(modulus):
    """Unit classes modulo sign, in increasing canonical representatives."""
    values = np.arange(1, (modulus + 1) // 2, dtype=np.int64)
    return values[values % P != 0]


def norm_numerators(values, modulus):
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def rows_to_bitsets(rows):
    packed = np.packbits(rows, axis=1, bitorder="little")
    return tuple(int.from_bytes(row.tobytes(), "little") for row in packed)


def capacity_on_phase_set(phase_set, at_least_one, at_least_two, index):
    return (
        (phase_set & at_least_one[index]).bit_count()
        + (phase_set & at_least_two[index]).bit_count()
    )


def run():
    phases = np.arange(1, Q, dtype=np.int64)
    phases = phases[phases % P != 0]
    unit_labels = sign_classes(N)
    depth_one_labels = sign_classes(Q)
    depth_two_labels = sign_classes(P**2)
    sheets = np.arange(P, dtype=np.int64)

    require(len(phases) == 2028, "bad primitive image-phase count")
    require(len(unit_labels) == 13182, "bad unit sign-class count")
    require(len(depth_one_labels) == 1014,
            "bad depth-one sign-class count")
    require(len(depth_two_labels) == 78,
            "bad depth-two sign-class count")
    require(N % 7 != 0 and N % 14 != 0 and Q % 14 != 0,
            "a torsion endpoint can occur")

    roots = phases[:, None] + Q * sheets[None, :]
    guard = 7 * norm_numerators(roots, N) > N
    guard_counts = guard.sum(axis=1)
    guard_histogram = tuple(
        (size, int(np.sum(guard_counts == size))) for size in (9, 10)
    )
    require(guard_histogram == ((9, 1450), (10, 578)),
            "bad guard-fibre histogram")
    require(int(guard_counts.sum()) == 18830,
            "bad primitive guard-safe universe size")
    guard_nine = rows_to_bitsets((guard_counts == 9)[None, :])[0]
    phase_full = (1 << len(phases)) - 1

    # A depth-one coefficient 13*u reduces to u on the image phase.
    depth_one_rows = (
        14
        * norm_numerators(
            depth_one_labels[:, None] * phases[None, :],
            Q,
        )
        < Q
    )
    depth_one_active = rows_to_bitsets(depth_one_rows)

    # A depth-two coefficient 13^2*v reduces to 13*v.
    depth_two_rows = (
        14
        * norm_numerators(
            (P * depth_two_labels[:, None]) * phases[None, :],
            Q,
        )
        < Q
    )
    depth_two_active = rows_to_bitsets(depth_two_rows)

    # Definition-level parity for every blocker class on every root.
    for start in range(0, len(depth_one_labels), 13):
        batch = depth_one_labels[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P * batch[:, None, None]) * roots[None, :, :],
                N,
            )
            < N
        )
        reduced = depth_one_rows[start : start + 13, :, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-one direct/image parity failed",
        )

    for start in range(0, len(depth_two_labels), 13):
        batch = depth_two_labels[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P**2 * batch[:, None, None]) * roots[None, :, :],
                N,
            )
            < N
        )
        reduced = depth_two_rows[start : start + 13, :, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-two direct/image parity failed",
        )

    # Unit phase-incidence bitsets h_q>=1 and h_q>=2.
    unit_at_least_one = []
    unit_at_least_two = []
    for start in range(0, len(unit_labels), 24):
        batch = unit_labels[start : start + 24]
        residues = (batch[:, None, None] * roots[None, :, :]) % N
        dangerous = 14 * np.minimum(residues, N - residues) < N
        counts = np.sum(dangerous & guard[None, :, :], axis=2)
        require(int(counts.max()) <= 2,
                "a unit phase has more than two active roots")
        unit_at_least_one.extend(rows_to_bitsets(counts >= 1))
        unit_at_least_two.extend(rows_to_bitsets(counts >= 2))
    unit_at_least_one = tuple(unit_at_least_one)
    unit_at_least_two = tuple(unit_at_least_two)

    unit_label_list = tuple(map(int, unit_labels))
    number_of_units = len(unit_label_list)
    full_capacities = np.fromiter(
        (
            unit_at_least_one[index].bit_count()
            + unit_at_least_two[index].bit_count()
            for index in range(number_of_units)
        ),
        dtype=np.int16,
        count=number_of_units,
    )

    # Exact one-depth-two-blocker capacities, used in every branch bound.
    depth_two_capacities = []
    for blocker in depth_two_active:
        complement = phase_full ^ blocker
        depth_two_capacities.append(
            np.fromiter(
                (
                    capacity_on_phase_set(
                        complement,
                        unit_at_least_one,
                        unit_at_least_two,
                        index,
                    )
                    for index in range(number_of_units)
                ),
                dtype=np.int16,
                count=number_of_units,
            )
        )
    depth_two_capacities = np.asarray(
        depth_two_capacities, dtype=np.int16
    )
    depth_two_removed = (
        full_capacities[None, :] - depth_two_capacities
    )

    pair_digest = hashlib.sha256()
    branch_digest = hashlib.sha256()
    pair_count = 0
    candidate_evaluations = 0
    maximum_candidates = 0
    candidate_histogram = Counter()
    best_margin = None
    best_record = None
    margin_minimizers = 0

    for depth_one_index, (left, left_mask) in enumerate(
        zip(map(int, depth_one_labels), depth_one_active)
    ):
        left_complement = phase_full ^ left_mask
        depth_one_capacities = np.fromiter(
            (
                capacity_on_phase_set(
                    left_complement,
                    unit_at_least_one,
                    unit_at_least_two,
                    index,
                )
                for index in range(number_of_units)
            ),
            dtype=np.int16,
            count=number_of_units,
        )
        depth_one_removed = full_capacities - depth_one_capacities

        for depth_two_index, (right, right_mask) in enumerate(
            zip(map(int, depth_two_labels), depth_two_active)
        ):
            pair_count += 1
            residual = left_complement & (phase_full ^ right_mask)
            residual_size = (
                10 * residual.bit_count()
                - (residual & guard_nine).bit_count()
            )
            phase_overlap = (left_mask & right_mask).bit_count()

            inclusion_bound = (
                full_capacities
                - depth_one_removed
                - depth_two_removed[depth_two_index]
                + 2 * phase_overlap
            )
            upper_bounds = np.minimum(
                np.minimum(
                    depth_one_capacities,
                    depth_two_capacities[depth_two_index],
                ),
                inclusion_bound,
            )
            order = np.argsort(-upper_bounds, kind="stable")

            top_five_heap = []
            evaluated = 0
            stopping_bound = None
            for raw_index in order:
                index = int(raw_index)
                bound = int(upper_bounds[index])
                if (
                    len(top_five_heap) == 5
                    and bound < top_five_heap[0][0]
                ):
                    stopping_bound = bound
                    break

                capacity = capacity_on_phase_set(
                    residual,
                    unit_at_least_one,
                    unit_at_least_two,
                    index,
                )
                require(capacity <= bound,
                        "an inclusion-exclusion branch bound failed")
                evaluated += 1
                item = (capacity, -unit_label_list[index])
                if len(top_five_heap) < 5:
                    heapq.heappush(top_five_heap, item)
                elif item > top_five_heap[0]:
                    heapq.heapreplace(top_five_heap, item)

            require(len(top_five_heap) == 5,
                    "branch ended before finding five capacities")
            require(stopping_bound is not None,
                    "branch did not stop before exhausting all units")
            top_five = tuple(
                sorted(
                    (
                        (capacity, -negative_label)
                        for capacity, negative_label in top_five_heap
                    ),
                    key=lambda entry: (-entry[0], entry[1]),
                )
            )
            require(stopping_bound < top_five[-1][0],
                    "strict branch stopping certificate failed")

            margin = residual_size - sum(
                capacity for capacity, _ in top_five
            )
            require(margin > 0,
                    "a mixed depth-(1,2,3) pair has nonpositive margin")

            record = (left, right, residual_size, top_five)
            pair_digest.update(repr((record, margin)).encode())
            branch_digest.update(
                repr(
                    (
                        left,
                        right,
                        phase_overlap,
                        evaluated,
                        stopping_bound,
                    )
                ).encode()
            )

            candidate_evaluations += evaluated
            maximum_candidates = max(maximum_candidates, evaluated)
            candidate_histogram[evaluated] += 1

            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_record = record
                margin_minimizers = 1
            elif margin == best_margin:
                margin_minimizers += 1

    require(pair_count == 1014 * 78 == 79092,
            "bad mixed shallow-pair count")
    require(candidate_evaluations == 940857,
            "candidate-evaluation total changed")
    require(maximum_candidates == 27,
            "maximum branch prefix changed")
    require(best_margin == 1608,
            "minimum mixed depth-(1,2,3) margin changed")
    require(margin_minimizers == 1,
            "minimum mixed margin is no longer unique")
    require(
        best_record
        == (
            799,
            46,
            13526,
            (
                (2604, 5193),
                (2472, 10386),
                (2292, 7773),
                (2288, 10388),
                (2262, 7775),
            ),
        ),
        "unique hostile row changed",
    )
    require(
        pair_digest.hexdigest()
        == "79b9b75f3732e47b43c2bba726906250ce0c1069b90d8952c33caa8a8364570f",
        "mixed pair-table digest changed",
    )
    require(
        branch_digest.hexdigest()
        == "2699c79c62d4c0e5805529b26dc4ce7ae5ed1f0146be25714b898ea42176802a",
        "mixed branch-trace digest changed",
    )
    require(
        tuple(sorted(candidate_histogram.items()))
        == (
            (5, 48),
            (6, 119),
            (7, 152),
            (8, 413),
            (9, 4110),
            (10, 10029),
            (11, 15964),
            (12, 19420),
            (13, 17035),
            (14, 9897),
            (15, 1530),
            (16, 218),
            (17, 66),
            (18, 9),
            (19, 6),
            (20, 6),
            (21, 11),
            (22, 12),
            (23, 17),
            (24, 8),
            (25, 12),
            (26, 6),
            (27, 4),
        ),
        "mixed candidate-prefix histogram changed",
    )

    # Independent direct-torsion reconstruction of every unit capacity on
    # the unique hostile row.
    all_residues = np.arange(N, dtype=np.int64)
    direct_universe = all_residues[
        (all_residues % P != 0)
        & (7 * norm_numerators(all_residues, N) > N)
    ]
    direct_residual = direct_universe[
        (
            14
            * norm_numerators(P * best_record[0] * direct_universe, N)
            >= N
        )
        & (
            14
            * norm_numerators(P**2 * best_record[1] * direct_universe, N)
            >= N
        )
    ]
    require(len(direct_universe) == 18830,
            "bad direct primitive guard-safe universe")
    require(len(direct_residual) == best_record[2],
            "direct hostile residual size changed")

    direct_capacities = []
    for start in range(0, len(unit_labels), 64):
        batch = unit_labels[start : start + 64]
        residues = (batch[:, None] * direct_residual[None, :]) % N
        direct_capacities.extend(
            map(
                int,
                np.sum(
                    14 * np.minimum(residues, N - residues) < N,
                    axis=1,
                ),
            )
        )
    direct_capacities = tuple(direct_capacities)

    hostile_left_index = tuple(map(int, depth_one_labels)).index(
        best_record[0]
    )
    hostile_right_index = tuple(map(int, depth_two_labels)).index(
        best_record[1]
    )
    hostile_residual = (
        (phase_full ^ depth_one_active[hostile_left_index])
        & (phase_full ^ depth_two_active[hostile_right_index])
    )
    phase_capacities = tuple(
        capacity_on_phase_set(
            hostile_residual,
            unit_at_least_one,
            unit_at_least_two,
            index,
        )
        for index in range(number_of_units)
    )
    require(direct_capacities == phase_capacities,
            "direct/phase parity failed on the hostile all-unit row")
    direct_top = tuple(
        sorted(
            zip(direct_capacities, unit_label_list),
            key=lambda entry: (-entry[0], entry[1]),
        )[:5]
    )
    require(direct_top == best_record[3],
            "direct hostile top five changed")

    # The deepest coefficient 13^3*w is safe on the primitive layer.
    for deepest_unit_part in range(1, 7):
        require(
            not np.any(
                14
                * norm_numerators(
                    P**3 * deepest_unit_part * direct_universe,
                    N,
                )
                < N
            ),
            "a depth-three blocker was active on the primitive layer",
        )

    # Lift-family partition and the exact THM-2200 sum law on the hostile
    # residual.  This also records the family sidecar behind the top five.
    families = defaultdict(list)
    for index, label in enumerate(unit_label_list):
        residue = label % Q
        base = min(residue, (-residue) % Q)
        families[base].append(index)
    require(len(families) == 1014, "bad lift-family count")
    require(set(map(len, families.values())) == {13},
            "a lift family does not contain thirteen sign classes")

    hostile_weight = best_record[2]
    base_danger_rows = (
        14
        * norm_numerators(
            depth_one_labels[:, None] * phases[None, :],
            Q,
        )
        < Q
    )
    base_danger = rows_to_bitsets(base_danger_rows)
    base_to_index = {
        int(label): index for index, label in enumerate(depth_one_labels)
    }
    family_sum_checks = 0
    for base, indices in families.items():
        total = sum(phase_capacities[index] for index in indices)
        danger = hostile_residual & base_danger[base_to_index[base]]
        weighted_danger = (
            10 * danger.bit_count()
            - (danger & guard_nine).bit_count()
        )
        require(total == 2 * hostile_weight - weighted_danger,
                "hostile lift-family sum identity failed")
        family_sum_checks += 1
    require(family_sum_checks == 1014,
            "bad hostile family-sum check count")

    top_family_bases = []
    top_family_rows = []
    for _, label in best_record[3]:
        residue = label % Q
        base = min(residue, (-residue) % Q)
        top_family_bases.append(base)
        top_family_rows.append(
            (
                base,
                tuple(
                    (
                        unit_label_list[index],
                        phase_capacities[index],
                    )
                    for index in families[base]
                ),
            )
        )
    require(
        tuple(top_family_bases) == (799, 599, 1015, 597, 1013),
        "hostile top-family bases changed",
    )

    print("scalar depth-(1,2,3) exact branch certificate")
    print("arithmetic: exact integers; no float; no assert dependence")
    print()
    print("depth-four primitive root geometry")
    print(
        f"  N={N} Q={Q} primitive_phases={len(phases)} "
        f"|U_N|={int(guard_counts.sum())}"
    )
    print(f"  guard_fibre_histogram={guard_histogram}")
    print(
        f"  unit_sign_classes={len(unit_labels)} "
        f"depth_one_sign_classes={len(depth_one_labels)} "
        f"depth_two_sign_classes={len(depth_two_labels)}"
    )
    print("  depth-one direct/image parity=1014/1014 PASS")
    print("  depth-two direct/image parity=78/78 PASS")
    print("  six depth-three sign classes empty on U_N PASS")
    print()
    print("mixed shallow-pair branch certificate")
    print(f"  mixed_pairs={pair_count}")
    print(
        "  bound=min(one-blocker capacities, "
        "F-X_A-X_B+2|A intersection B|)"
    )
    print(
        f"  candidate_evaluations={candidate_evaluations} "
        f"average={candidate_evaluations / pair_count:.12f} "
        f"maximum={maximum_candidates}"
    )
    print(
        f"  unique_worst_pair=({best_record[0]},{best_record[1]}) "
        f"residual_size={best_record[2]}"
    )
    print(f"  top_five_(capacity,label)={best_record[3]}")
    print(
        f"  minimum_capacity_margin={best_margin} "
        f"minimizers={margin_minimizers}"
    )
    print(f"  pair_table_digest={pair_digest.hexdigest()}")
    print(f"  branch_trace_digest={branch_digest.hexdigest()}")
    print(
        "  candidate_prefix_histogram="
        f"{tuple(sorted(candidate_histogram.items()))}"
    )
    print("  direct all-unit hostile-row parity=PASS")
    print()
    print("hostile lift-family sidecar")
    print(f"  top_family_bases={tuple(top_family_bases)}")
    for base, rows in top_family_rows:
        print(
            f"  base={base} family_sum="
            f"{sum(capacity for _, capacity in rows)} rows={rows}"
        )
    print(f"  hostile_family_sum_checks={family_sum_checks} PASS")
    print()
    print("CONCLUSION")
    print("  actual blocker valuation profile (1,2,3) is empty")
    print("  this certificate does not address profile (1,1,3)")


if __name__ == "__main__":
    run()
