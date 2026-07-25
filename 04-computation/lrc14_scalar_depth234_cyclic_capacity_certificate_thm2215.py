#!/usr/bin/env python3
"""Literal-integer cyclic-capacity certificate for scalar profile (2,3,4).

At N=13^5, group the 26,364 primitive image phases modulo 13^4 by
their oriented residue modulo 13^3.  Each of the 2,028 groups has thirteen
phases.  The array H[q,s] is the exact guard-safe dangerous-root capacity
of unit sign class q over group s; it is uint8 because H[q,s] <= 26.

A depth-two blocker is constant on these groups and has 1,014 possible
sign classes.  A depth-three blocker has 78 possible sign classes.  For
every one of the 1,014*78 mixed pairs, the exact top-five unit capacities
are certified with the overlap-prefix bound

  C_q <= min(C_q(complement A), C_q(complement B),
             F_q-X_q(A)-X_q(B)+2*13*|A intersection B|).

The large removal tables use uint16 matrix multiplication.  This is
literal integer arithmetic without wraparound: every removal is a sub-sum
of a full capacity smaller than 2^16.  No float and no Python assert is
used, so all checks remain active under -O.
"""

from collections import Counter
import hashlib
import heapq

import numpy as np


P = 13
N = P**5
Q = P**4
R2 = P**3
R3 = P**2


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign_classes(modulus):
    values = np.arange(1, (modulus + 1) // 2, dtype=np.int64)
    return values[values % P != 0]


def norm_numerators(values, modulus):
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def family_group_capacities(
    base,
    phase_grid,
    group_residues,
    guard_safe_by_label,
    shifts,
):
    inverse = pow(int(base), -1, P)
    product = int(base) * phase_grid
    quotient = product // Q
    remainder = product % Q
    first_present = 14 * remainder < 13 * Q
    second_present = 14 * remainder > Q
    first_label = (inverse * quotient) % P
    lifted_first = (
        first_label[:, :, None]
        + shifts[inverse][:, None, :]
    ) % P
    lifted_second = (lifted_first + inverse) % P
    counts = (
        np.take_along_axis(
            guard_safe_by_label, lifted_first, axis=2
        )
        & first_present[:, :, None]
    ).astype(np.uint8)
    counts += (
        np.take_along_axis(
            guard_safe_by_label, lifted_second, axis=2
        )
        & second_present[:, :, None]
    )
    return counts.sum(axis=1, dtype=np.uint8)


def run():
    group_residues = np.arange(1, R2, dtype=np.int64)
    group_residues = group_residues[group_residues % P != 0]
    within_group = np.arange(Q // R2, dtype=np.int64)
    phase_grid = (
        group_residues[:, None] + R2 * within_group[None, :]
    )
    root_indices = np.arange(P, dtype=np.int64)
    roots = (
        phase_grid[:, :, None]
        + Q * root_indices[None, None, :]
    )
    guard_by_root = 7 * norm_numerators(roots, N) > N

    require(len(group_residues) == 2028,
            "bad oriented residue count modulo 2197")
    require(phase_grid.size == 26364,
            "bad primitive image-phase count")
    require(int(guard_by_root.sum()) == 244810,
            "bad primitive guard-safe universe size")
    guard_histogram = tuple(
        (
            size,
            int(np.sum(guard_by_root.sum(axis=2) == size)),
        )
        for size in (9, 10)
    )
    require(guard_histogram == ((9, 18830), (10, 7534)),
            "bad guard-fibre histogram")
    require(
        N % 7 != 0
        and N % 14 != 0
        and Q % 14 != 0
        and R2 % 14 != 0
        and R3 % 14 != 0,
        "a torsion endpoint can occur",
    )

    guard_safe_by_label = np.empty_like(guard_by_root)
    for root_index in range(P):
        guard_safe_by_label[:, :, (-root_index) % P] = (
            guard_by_root[:, :, root_index]
        )

    lift_parameters = np.arange(P, dtype=np.int64)
    shifts = {
        inverse: (
            inverse
            * (group_residues[:, None] % P)
            * lift_parameters[None, :]
        ) % P
        for inverse in range(1, P)
    }

    base_labels = sign_classes(Q)
    unit_labels = np.empty(len(base_labels) * P, dtype=np.int64)
    capacities_by_group = np.empty(
        (len(unit_labels), len(group_residues)),
        dtype=np.uint8,
    )
    for family_index, base in enumerate(map(int, base_labels)):
        start = family_index * P
        capacities_by_group[start : start + P] = (
            family_group_capacities(
                base,
                phase_grid,
                group_residues,
                guard_safe_by_label,
                shifts,
            ).T
        )
        raw_lifts = base + Q * lift_parameters
        unit_labels[start : start + P] = np.minimum(
            raw_lifts, N - raw_lifts
        )

    require(
        np.array_equal(np.sort(unit_labels), sign_classes(N)),
        "lift families do not partition unit sign classes modulo 13^5",
    )
    require(int(capacities_by_group.max()) <= 26,
            "a group capacity exceeded 2*13")

    depth_two_labels = sign_classes(R2)
    depth_three_labels = sign_classes(R3)
    depth_two_active = (
        14
        * norm_numerators(
            depth_two_labels[:, None]
            * group_residues[None, :],
            R2,
        )
        < R2
    )
    depth_three_active = (
        14
        * norm_numerators(
            depth_three_labels[:, None]
            * (group_residues[None, :] % R3),
            R3,
        )
        < R3
    )
    require(len(depth_two_labels) == 1014,
            "bad depth-two blocker sign-class count")
    require(len(depth_three_labels) == 78,
            "bad depth-three blocker sign-class count")

    # Every blocker class is checked at definition level on every root.
    for start in range(0, len(depth_two_labels), 8):
        batch = depth_two_labels[start : start + 8]
        direct = (
            14
            * norm_numerators(
                (P**2 * batch[:, None, None, None])
                * roots[None, :, :, :],
                N,
            )
            < N
        )
        reduced = depth_two_active[start : start + 8, :, None, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-two blocker direct/group parity failed",
        )

    for start in range(0, len(depth_three_labels), 13):
        batch = depth_three_labels[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P**3 * batch[:, None, None, None])
                * roots[None, :, :, :],
                N,
            )
            < N
        )
        reduced = depth_three_active[start : start + 13, :, None, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-three blocker direct/group parity failed",
        )

    full_capacities = capacities_by_group.sum(
        axis=1, dtype=np.uint32
    )
    require(int(full_capacities.max()) < 2**16,
            "a full unit capacity exceeds uint16")

    # Every blocker removal is a sub-sum of a full capacity.  Since the
    # latter is <2^16, these uint16 dot products cannot wrap.
    integer_capacity_matrix = capacities_by_group.astype(np.uint16)
    depth_two_removed = (
        integer_capacity_matrix
        @ depth_two_active.T.astype(np.uint16)
    )
    depth_three_removed = (
        integer_capacity_matrix
        @ depth_three_active.T.astype(np.uint16)
    )
    del integer_capacity_matrix
    require(
        int(depth_two_removed.max()) <= int(full_capacities.max())
        and int(depth_three_removed.max()) <= int(full_capacities.max()),
        "a blocker removal exceeded its full capacity",
    )
    full_capacities = full_capacities.astype(np.int32)

    # Independent literal-integer hostile column checks for both tables.
    for index in (0, 1, 500, 1013):
        exact = capacities_by_group[
            :, depth_two_active[index]
        ].sum(axis=1, dtype=np.uint32)
        require(
            np.array_equal(
                exact,
                depth_two_removed[:, index].astype(np.uint32),
            ),
            "depth-two removal table parity failed",
        )
    for index in (0, 1, 46, 77):
        exact = capacities_by_group[
            :, depth_three_active[index]
        ].sum(axis=1, dtype=np.uint32)
        require(
            np.array_equal(
                exact,
                depth_three_removed[:, index].astype(np.uint32),
            ),
            "depth-three removal table parity failed",
        )

    guard_mass_by_group = guard_by_root.sum(
        axis=(1, 2), dtype=np.uint16
    )
    pair_digest = hashlib.sha256()
    branch_digest = hashlib.sha256()
    candidate_histogram = Counter()
    pair_count = 0
    candidate_evaluations = 0
    maximum_prefix = 0
    branch_expansions = 0
    best_margin = None
    best_record = None
    margin_minimizers = 0

    for depth_two_index, left in enumerate(
        map(int, depth_two_labels)
    ):
        left_removed = depth_two_removed[
            :, depth_two_index
        ].astype(np.int32)
        left_capacity = full_capacities - left_removed

        for depth_three_index, right_raw in enumerate(
            map(int, depth_three_labels)
        ):
            right = int(right_raw)
            pair_count += 1
            right_removed = depth_three_removed[
                :, depth_three_index
            ].astype(np.int32)
            right_capacity = full_capacities - right_removed
            active_union = (
                depth_two_active[depth_two_index]
                | depth_three_active[depth_three_index]
            )
            active_intersection = (
                depth_two_active[depth_two_index]
                & depth_three_active[depth_three_index]
            )
            residual_size = int(
                guard_mass_by_group[~active_union].sum(
                    dtype=np.uint32
                )
            )
            inclusion_bound = (
                full_capacities
                - left_removed
                - right_removed
                + 2
                * (Q // R2)
                * int(active_intersection.sum())
            )
            upper_bound = np.minimum(
                np.minimum(left_capacity, right_capacity),
                inclusion_bound,
            )

            prefix_size = 64
            while True:
                partition = np.argpartition(
                    -upper_bound, prefix_size
                )
                candidates = partition[:prefix_size]
                stopping_bound = int(
                    upper_bound[partition[prefix_size]]
                )
                ordered = candidates[
                    np.argsort(
                        -upper_bound[candidates], kind="stable"
                    )
                ]
                heap = []
                evaluated = 0
                for raw_index in ordered:
                    index = int(raw_index)
                    bound = int(upper_bound[index])
                    if len(heap) == 5 and bound < heap[0][0]:
                        break
                    exact = int(
                        full_capacities[index]
                        - left_removed[index]
                        - right_removed[index]
                        + capacities_by_group[
                            index, active_intersection
                        ].sum(dtype=np.uint32)
                    )
                    require(exact <= bound,
                            "overlap branch bound failed")
                    evaluated += 1
                    item = (exact, -int(unit_labels[index]))
                    if len(heap) < 5:
                        heapq.heappush(heap, item)
                    elif item > heap[0]:
                        heapq.heapreplace(heap, item)
                if (
                    len(heap) == 5
                    and stopping_bound < heap[0][0]
                ):
                    break
                prefix_size *= 2
                branch_expansions += 1
                require(prefix_size < len(unit_labels),
                        "branch prefix exhausted all units")

            top_five = tuple(
                sorted(
                    (
                        (capacity, -negative_label)
                        for capacity, negative_label in heap
                    ),
                    key=lambda entry: (-entry[0], entry[1]),
                )
            )
            margin = residual_size - sum(
                capacity for capacity, _ in top_five
            )
            require(margin > 0,
                    "a depth-(2,3,4) pair has nonpositive margin")
            record = (left, right, residual_size, top_five)
            pair_digest.update(repr((record, margin)).encode())
            branch_digest.update(
                repr(
                    (
                        left,
                        right,
                        int(active_intersection.sum()),
                        prefix_size,
                        evaluated,
                        stopping_bound,
                    )
                ).encode()
            )
            candidate_evaluations += evaluated
            maximum_prefix = max(maximum_prefix, prefix_size)
            candidate_histogram[evaluated] += 1

            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_record = record
                margin_minimizers = 1
            elif margin == best_margin:
                margin_minimizers += 1

    require(pair_count == 1014 * 78 == 79092,
            "bad mixed blocker-pair count")
    require(candidate_evaluations == 878784,
            "candidate-evaluation total changed")
    require(maximum_prefix == 64,
            "maximum branch prefix changed")
    require(branch_expansions == 0,
            "a branch prefix unexpectedly expanded")
    require(best_margin == 27440,
            "minimum depth-(2,3,4) margin changed")
    require(margin_minimizers == 1,
            "minimum depth-(2,3,4) margin is not unique")
    require(
        best_record
        == (
            844,
            1,
            175248,
            (
                (29664, 142637),
                (29646, 2198),
                (29646, 142635),
                (29628, 2196),
                (29224, 6),
            ),
        ),
        "unique hostile depth-(2,3,4) row changed",
    )
    require(
        pair_digest.hexdigest()
        == "f141ef307efc773fb90e6a20c534f4a2d754b31bd1f63658591d34551161af3f",
        "depth-(2,3,4) pair-table digest changed",
    )

    # Direct hostile residual and reported-capacity replay.
    direct_universe = roots[guard_by_root]
    direct_residual = direct_universe[
        (
            14
            * norm_numerators(
                P**2 * best_record[0] * direct_universe, N
            )
            >= N
        )
        & (
            14
            * norm_numerators(
                P**3 * best_record[1] * direct_universe, N
            )
            >= N
        )
    ]
    require(len(direct_residual) == best_record[2],
            "direct hostile residual size changed")
    direct_top = tuple(
        (
            int(
                np.sum(
                    14
                    * norm_numerators(label * direct_residual, N)
                    < N
                )
            ),
            label,
        )
        for _, label in best_record[3]
    )
    require(direct_top == best_record[3],
            "direct hostile top capacities changed")

    for deepest_unit_part in range(1, 7):
        require(
            not np.any(
                14
                * norm_numerators(
                    P**4
                    * deepest_unit_part
                    * direct_universe,
                    N,
                )
                < N
            ),
            "a depth-four blocker is active on the primitive layer",
        )

    # Exact lift-family sums on the hostile residual.
    hostile_left = tuple(map(int, depth_two_labels)).index(
        best_record[0]
    )
    hostile_right = tuple(map(int, depth_three_labels)).index(
        best_record[1]
    )
    hostile_allowed_groups = ~(
        depth_two_active[hostile_left]
        | depth_three_active[hostile_right]
    )
    hostile_phase_mask = np.broadcast_to(
        hostile_allowed_groups[:, None], phase_grid.shape
    )
    guard_counts = guard_by_root.sum(axis=2)
    family_sum_checks = 0
    for family_index, base in enumerate(map(int, base_labels)):
        family_total = int(
            capacities_by_group[
                family_index * P : (family_index + 1) * P,
                hostile_allowed_groups,
            ].sum(dtype=np.uint32)
        )
        danger = (
            14 * norm_numerators(base * phase_grid, Q) < Q
        ) & hostile_phase_mask
        weighted_danger = int(
            guard_counts[danger].sum(dtype=np.uint32)
        )
        require(
            family_total
            == 2 * best_record[2] - weighted_danger,
            "hostile thirteen-lift family-sum identity failed",
        )
        family_sum_checks += 1
    require(family_sum_checks == 13182,
            "bad hostile family-sum check count")

    print("scalar depth-(2,3,4) cyclic-capacity certificate")
    print("arithmetic: exact integers; no float; no assert dependence")
    print()
    print("primitive depth-five universe")
    print(
        f"  N={N} Q={Q} primitive_phases={phase_grid.size} "
        f"|U_N|={int(guard_by_root.sum())}"
    )
    print(f"  guard_fibre_histogram={guard_histogram}")
    print(
        f"  unit_sign_classes={len(unit_labels)} "
        f"lift_families={len(base_labels)} "
        f"depth_two_sign_classes={len(depth_two_labels)} "
        f"depth_three_sign_classes={len(depth_three_labels)}"
    )
    print("  all depth-two blocker direct/group parity=1014/1014 PASS")
    print("  all depth-three blocker direct/group parity=78/78 PASS")
    print("  all six depth-four sign classes empty on U_N PASS")
    print()
    print("complete depth-(2,3,4) branch certificate")
    print(f"  mixed_pairs={pair_count}")
    print(
        f"  candidate_evaluations={candidate_evaluations} "
        f"maximum_prefix={maximum_prefix} "
        f"expansions={branch_expansions}"
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
        "  candidate_evaluation_histogram="
        f"{tuple(sorted(candidate_histogram.items()))}"
    )
    print("  direct hostile residual/top-five replay=PASS")
    print(f"  hostile_family_sum_checks={family_sum_checks} PASS")
    print()
    print("CONCLUSION")
    print("  actual blocker valuation profile (2,3,4) is empty")


if __name__ == "__main__":
    run()
