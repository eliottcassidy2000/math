#!/usr/bin/env python3
"""Exact cyclic-capacity certificate for scalar profile (3,3,4).

The primitive 13^5 layer has 26,364 image phases modulo 13^4 and
171,366 unit sign classes.  A naive pair of phase bitsets for every unit
would use more than one gigabyte.  Here phases are grouped by their
oriented residue modulo 13^2.  Each group contains 13^2 phases, and every
depth-three blocker is constant on a group.  For each unit coefficient the
array H[q,s] stores its exact guard-safe dangerous-root capacity over group
s.  H is uint16 because 0 <= H[q,s] <= 2*13^2=338.

For active blocker group sets A,B, the branch bound is

  C_q <= min(C_q(complement A), C_q(complement B),
             F_q-X_q(A)-X_q(B)+2*13^2*|A intersection B|).

All arithmetic is integer arithmetic.  Checks remain active under -O.
"""

from collections import Counter
import hashlib
import heapq

import numpy as np


P = 13
N = P**5
Q = P**4
R = P**2


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
    """Return H[s,l] for q_l=base+lQ, exactly."""
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
    return counts.sum(axis=1, dtype=np.uint16)


def run():
    group_residues = np.arange(1, R, dtype=np.int64)
    group_residues = group_residues[group_residues % P != 0]
    within_group = np.arange(Q // R, dtype=np.int64)
    phase_grid = (
        group_residues[:, None] + R * within_group[None, :]
    )
    root_indices = np.arange(P, dtype=np.int64)
    roots = (
        phase_grid[:, :, None]
        + Q * root_indices[None, None, :]
    )
    guard_by_root = 7 * norm_numerators(roots, N) > N

    require(len(group_residues) == 156,
            "bad oriented residue count modulo 169")
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
        and R % 14 != 0,
        "a torsion endpoint can occur",
    )

    # Convert root index k to the reversed guard label -k.
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
        dtype=np.uint16,
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

    expected_unit_labels = sign_classes(N)
    require(
        np.array_equal(np.sort(unit_labels), expected_unit_labels),
        "lift families do not partition unit sign classes modulo 13^5",
    )
    require(int(capacities_by_group.max()) <= 338,
            "a group capacity exceeded 2*13^2")

    blocker_labels = sign_classes(R)
    blocker_active = (
        14
        * norm_numerators(
            blocker_labels[:, None] * group_residues[None, :],
            R,
        )
        < R
    )
    require(len(blocker_labels) == 78,
            "bad depth-three blocker sign-class count")

    # All-class definition-level blocker parity on every root.
    for start in range(0, len(blocker_labels), 13):
        batch = blocker_labels[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P**3 * batch[:, None, None, None])
                * roots[None, :, :, :],
                N,
            )
            < N
        )
        reduced = blocker_active[start : start + 13, :, None, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-three blocker direct/group parity failed",
        )

    full_capacities = capacities_by_group.sum(
        axis=1, dtype=np.uint32
    ).astype(np.int32)
    # True removal totals are below 2^16, so uint16 matrix multiplication
    # is literal exact integer arithmetic without wraparound.
    removed = (
        capacities_by_group
        @ blocker_active.T.astype(np.uint16)
    )
    require(int(removed.max()) < 2**16,
            "uint16 removal accumulation can wrap")
    removed = removed.astype(np.int32)

    guard_mass_by_group = guard_by_root.sum(
        axis=(1, 2), dtype=np.uint32
    )
    pair_digest = hashlib.sha256()
    branch_digest = hashlib.sha256()
    candidate_histogram = Counter()
    pair_count = 0
    diagonal_count = 0
    candidate_evaluations = 0
    maximum_prefix = 0
    branch_expansions = 0
    best_margin = None
    best_record = None
    margin_minimizers = 0

    for left_index, left in enumerate(map(int, blocker_labels)):
        left_capacity = full_capacities - removed[:, left_index]
        for right_index in range(left_index, len(blocker_labels)):
            right = int(blocker_labels[right_index])
            pair_count += 1
            diagonal_count += left_index == right_index
            active_union = (
                blocker_active[left_index]
                | blocker_active[right_index]
            )
            active_intersection = (
                blocker_active[left_index]
                & blocker_active[right_index]
            )
            residual_size = int(
                guard_mass_by_group[~active_union].sum(
                    dtype=np.uint32
                )
            )

            if left_index == right_index:
                candidate_indices = np.argpartition(
                    left_capacity, -5
                )[-5:]
                top_five = tuple(
                    sorted(
                        (
                            (
                                int(left_capacity[index]),
                                int(unit_labels[index]),
                            )
                            for index in candidate_indices
                        ),
                        key=lambda entry: (-entry[0], entry[1]),
                    )
                )
                evaluated = 5
                prefix_size = 5
                stopping_bound = top_five[-1][0] - 1
            else:
                right_capacity = (
                    full_capacities - removed[:, right_index]
                )
                inclusion_bound = (
                    full_capacities
                    - removed[:, left_index]
                    - removed[:, right_index]
                    + 2
                    * (Q // R)
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
                            - removed[index, left_index]
                            - removed[index, right_index]
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
                    "a depth-(3,3,4) pair has nonpositive margin")
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

    require(pair_count == 3081,
            "bad unordered blocker-pair count")
    require(diagonal_count == 78,
            "diagonal blocker pairs were lost")
    require(candidate_evaluations == 35974,
            "candidate-evaluation total changed")
    require(maximum_prefix == 64,
            "maximum branch prefix changed")
    require(branch_expansions == 0,
            "a branch prefix unexpectedly expanded")
    require(best_margin == 24022,
            "minimum depth-(3,3,4) margin changed")
    require(margin_minimizers == 1,
            "minimum depth-(3,3,4) margin is not unique")
    require(
        best_record
        == (
            1,
            84,
            188290,
            (
                (33928, 2198),
                (33918, 2196),
                (32518, 1098),
                (32496, 1099),
                (31408, 6),
            ),
        ),
        "unique hostile depth-(3,3,4) row changed",
    )
    require(
        pair_digest.hexdigest()
        == "97ae53fb2a59bf4bd8a34e18dfb919063ef502eed938fda9fb30c0c901fdea68",
        "depth-(3,3,4) pair-table digest changed",
    )
    require(
        branch_digest.hexdigest()
        == "3f92815e67a77af35d97d44cbf50cd5622b444b0d1e13917e8fd44bc549aac3f",
        "depth-(3,3,4) branch-trace digest changed",
    )
    require(
        tuple(sorted(candidate_histogram.items()))
        == (
            (5, 78),
            (8, 312),
            (9, 201),
            (10, 190),
            (11, 109),
            (12, 1220),
            (13, 82),
            (14, 861),
            (15, 28),
        ),
        "depth-(3,3,4) candidate histogram changed",
    )

    # Direct torsion replay of the hostile residual and its reported top
    # capacities, independent of the grouped-capacity representation.
    direct_universe = roots[guard_by_root]
    direct_residual = direct_universe[
        (
            14
            * norm_numerators(
                P**3 * best_record[0] * direct_universe, N
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
                    * norm_numerators(
                        label * direct_residual, N
                    )
                    < N
                )
            ),
            label,
        )
        for _, label in best_record[3]
    )
    require(direct_top == best_record[3],
            "direct hostile top capacities changed")

    # A depth-four blocker is a nonzero thirteenth root everywhere on the
    # primitive layer.  Check all six sign classes directly.
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

    # All lift families satisfy the exact family-sum identity on the
    # hostile residual.
    hostile_left = tuple(map(int, blocker_labels)).index(
        best_record[0]
    )
    hostile_right = tuple(map(int, blocker_labels)).index(
        best_record[1]
    )
    hostile_allowed_groups = ~(
        blocker_active[hostile_left]
        | blocker_active[hostile_right]
    )
    hostile_phase_mask = np.broadcast_to(
        hostile_allowed_groups[:, None], phase_grid.shape
    )
    hostile_guard_mass = best_record[2]
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
            == 2 * hostile_guard_mass - weighted_danger,
            "hostile thirteen-lift family-sum identity failed",
        )
        family_sum_checks += 1
    require(family_sum_checks == 13182,
            "bad hostile family-sum check count")

    print("scalar depth-(3,3,4) cyclic-capacity certificate")
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
        f"depth_three_sign_classes={len(blocker_labels)}"
    )
    print("  all depth-three blocker direct/group parity=78/78 PASS")
    print("  all six depth-four sign classes empty on U_N PASS")
    print()
    print("complete depth-(3,3,4) branch certificate")
    print(
        f"  unordered_pairs={pair_count} "
        f"including_diagonal={diagonal_count}"
    )
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
    print("  actual blocker valuation profile (3,3,4) is empty")


if __name__ == "__main__":
    run()
