#!/usr/bin/env python3
"""Exact audit for THM-2207: scalar blocker profile (1,2,3).

At N=13^4, write every primitive residue uniquely as

    z = r + k*13^3,   r a unit modulo 13^3,   k in F_13.

For a unit coefficient q, h_q(r) is the number (zero, one, or two) of
these roots which are simultaneously guard-safe and q-dangerous.  A
depth-one blocker and a depth-two blocker are fibre-constant; their safe
phase masks are denoted P_1(u) and P_2(v).  The exact conditional capacity
of q is therefore

    C_q(u,v) = sum_r P_1(u,r) P_2(v,r) h_q(r).

The complete scan covers all 1,014*78 typed shallow pairs and all 13,182
unit sign classes.  To keep the computation small while retaining an exact
certificate, it bounds C_q(u,v) by three independently valid quantities,
sorts those bounds, and evaluates exact bitset capacities until the next
bound lies below the fifth exact capacity.

The script also checks a simpler independent certificate.  Since the pair
residual is contained in each one-owner residual, its top-five capacity is
at most the smaller of the two one-owner top-five capacities.  This proves
all but one typed shallow pair; the sole exception is then evaluated
directly.

All arithmetic is exact.  ``require`` remains active under ``python -O``.
"""

import hashlib

import numpy as np


P = 13
N = P**4
Q = P**3
Q2 = P**2


def require(condition, message):
    """Raise on a failed exact check, including under optimized Python."""
    if not condition:
        raise RuntimeError(message)


def sign_classes(modulus):
    """Represent unit classes modulo sign by 1 <= a < modulus/2."""
    values = np.arange(1, (modulus + 1) // 2, dtype=np.int64)
    return values[values % P != 0]


def norm_numerators(values, modulus):
    """Return modulus times the circle norm of values/modulus."""
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def rows_to_bitsets(rows):
    """Pack Boolean matrix rows as little-endian Python-integer bitsets."""
    packed = np.packbits(rows, axis=1, bitorder="little")
    return tuple(int.from_bytes(row.tobytes(), "little") for row in packed)


def top_five_insert(top, item):
    """Maintain five entries ordered by decreasing capacity, then label."""
    top.append(item)
    top.sort(key=lambda entry: (-entry[0], entry[1]))
    if len(top) > 5:
        top.pop()


def direct_residual(universe, depth_one, depth_two):
    """Definition-level residual for one typed shallow pair."""
    first = P * depth_one
    second = P**2 * depth_two
    first_safe = (
        14 * norm_numerators(first * universe, N) >= N
    )
    second_safe = (
        14 * norm_numerators(second * universe, N) >= N
    )
    return universe[first_safe & second_safe]


def direct_capacity(universe, coefficient):
    """Definition-level strict unit-danger capacity on a residue set."""
    return int(
        np.sum(
            14 * norm_numerators(coefficient * universe, N) < N
        )
    )


def run():
    phases = np.arange(1, Q, dtype=np.int64)
    phases = phases[phases % P != 0]
    sheets = np.arange(P, dtype=np.int64)
    roots = phases[:, None] + Q * sheets[None, :]

    units = sign_classes(N)
    depth_one_units = sign_classes(Q)
    depth_two_units = sign_classes(Q2)
    deepest_units = sign_classes(P)

    require(len(phases) == 2028, "bad primitive image-phase count")
    require(len(units) == 13182, "bad depth-four unit sign-class count")
    require(
        len(depth_one_units) == 1014,
        "bad depth-one blocker sign-class count",
    )
    require(
        len(depth_two_units) == 78,
        "bad depth-two blocker sign-class count",
    )
    require(
        len(deepest_units) == 6,
        "bad depth-three blocker sign-class count",
    )
    require(
        all(modulus % 7 and modulus % 14 for modulus in (N, Q, Q2)),
        "a strict torsion endpoint appeared",
    )

    # The guard-safe roots above a primitive phase have nine or ten sheets.
    guard = 7 * norm_numerators(roots, N) > N
    guard_counts = guard.sum(axis=1)
    guard_histogram = tuple(
        (size, int(np.sum(guard_counts == size))) for size in (9, 10)
    )
    require(
        guard_histogram == ((9, 1450), (10, 578)),
        "bad depth-four guard-fibre histogram",
    )
    require(
        int(guard_counts.sum()) == 18830,
        "bad primitive guard-safe annulus size",
    )
    guard_nine = rows_to_bitsets((guard_counts == 9)[None, :])[0]
    phase_full = (1 << len(phases)) - 1

    # The unique depth-three blocker is safe on the entire primitive layer,
    # not merely on the selected guard-safe roots.  Its value at z/N is the
    # nonzero thirteenth root wz/13.  Check all six unit sign classes
    # directly in the original modulus.
    deepest_residues = (
        (P**3 * deepest_units[:, None, None]) * roots[None, :, :]
    ) % N
    require(
        np.all(
            14
            * np.minimum(deepest_residues, N - deepest_residues)
            >= N
        ),
        "a depth-three blocker was dangerous on a primitive root",
    )

    # Exact two-threshold unit signatures.  The capacity of q on a phase
    # mask R is popcount(R & {h_q>=1}) + popcount(R & {h_q>=2}).
    unit_counts = np.empty((len(units), len(phases)), dtype=np.uint8)
    for start in range(0, len(units), 32):
        batch = units[start : start + 32]
        residues = (
            batch[:, None, None] * roots[None, :, :]
        ) % N
        danger = 14 * np.minimum(residues, N - residues) < N
        counts = np.sum(danger & guard[None, :, :], axis=2)
        require(
            int(counts.max()) <= 2,
            "a unit root window had more than two sheets",
        )
        unit_counts[start : start + len(batch)] = counts

    unit_at_least_one = rows_to_bitsets(unit_counts >= 1)
    unit_at_least_two = rows_to_bitsets(unit_counts >= 2)
    full_capacities = unit_counts.sum(axis=1, dtype=np.uint16)

    # A depth-one coefficient 13u is dangerous above r exactly when ur/Q
    # is dangerous.  A depth-two coefficient 13^2 v depends only on r
    # modulo 13^2.
    depth_one_active_rows = (
        14
        * norm_numerators(
            depth_one_units[:, None] * phases[None, :],
            Q,
        )
        < Q
    )
    depth_two_active_rows = (
        14
        * norm_numerators(
            depth_two_units[:, None] * (phases[None, :] % Q2),
            Q2,
        )
        < Q2
    )
    depth_one_active = rows_to_bitsets(depth_one_active_rows)
    depth_two_active = rows_to_bitsets(depth_two_active_rows)

    # Direct/image hostile parity for every shallow sign class.
    for start in range(0, len(depth_one_units), 13):
        batch = depth_one_units[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P * batch[:, None, None]) * roots[None, :, :],
                N,
            )
            < N
        )
        reduced = depth_one_active_rows[start : start + 13, :, None]
        require(
            np.array_equal(
                direct, np.broadcast_to(reduced, direct.shape)
            ),
            "depth-one direct/image fibre parity failed",
        )

    for start in range(0, len(depth_two_units), 13):
        batch = depth_two_units[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P**2 * batch[:, None, None]) * roots[None, :, :],
                N,
            )
            < N
        )
        reduced = depth_two_active_rows[start : start + 13, :, None]
        require(
            np.array_equal(
                direct, np.broadcast_to(reduced, direct.shape)
            ),
            "depth-two direct/image fibre parity failed",
        )

    # C_q(not A) and C_q(not B), for every unit q and shallow owner.
    # These are two of the three exact upper bounds used below, and their
    # top-five sums supply the independent one-owner certificate.
    single_depth_one = np.empty(
        (len(units), len(depth_one_units)), dtype=np.uint16
    )
    for index, active in enumerate(depth_one_active_rows):
        single_depth_one[:, index] = (
            full_capacities
            - unit_counts[:, active].sum(axis=1, dtype=np.uint16)
        )

    single_depth_two = np.empty(
        (len(units), len(depth_two_units)), dtype=np.uint16
    )
    for index, active in enumerate(depth_two_active_rows):
        single_depth_two[:, index] = (
            full_capacities
            - unit_counts[:, active].sum(axis=1, dtype=np.uint16)
        )

    top_five_depth_one = np.partition(
        single_depth_one, -5, axis=0
    )[-5:, :].sum(axis=0, dtype=np.uint32)
    top_five_depth_two = np.partition(
        single_depth_two, -5, axis=0
    )[-5:, :].sum(axis=0, dtype=np.uint32)

    table_digest = hashlib.sha256()
    best_margin = None
    best_record = None
    margin_minimizers = 0
    pair_count = 0
    exact_candidate_evaluations = 0
    maximum_candidates = 0

    coarse_nonpositive = []
    coarse_minimum_positive = None
    coarse_minimum_positive_record = None

    # Keep the two rows used by the independent direct-residue controls.
    requested_records = {}
    requested_pairs = {(799, 46), (1098, 84)}

    full_capacities_i32 = full_capacities.astype(np.int32)

    for i, left in enumerate(map(int, depth_one_units)):
        cap_not_left = single_depth_one[:, i].astype(np.int32)
        loss_left = full_capacities_i32 - cap_not_left
        left_active = depth_one_active[i]

        for j, right in enumerate(map(int, depth_two_units)):
            pair_count += 1
            cap_not_right = single_depth_two[:, j].astype(np.int32)
            loss_right = full_capacities_i32 - cap_not_right
            right_active = depth_two_active[j]

            active_intersection = left_active & right_active
            intersection_phase_count = active_intersection.bit_count()
            residual = phase_full ^ (left_active | right_active)
            residual_phase_count = residual.bit_count()
            residual_size = (
                10 * residual_phase_count
                - (residual & guard_nine).bit_count()
            )

            # The exact capacity is
            #
            # F_q-X_q(A)-X_q(B)+X_q(A intersect B).
            #
            # Since h_q(r)<=2, the last term is bounded by twice the
            # number of phases in A intersect B.  The residual is also
            # contained in not-A and not-B.  Taking the minimum gives a
            # valid integer upper bound for every q.
            overlap_bound = (
                full_capacities_i32
                - loss_left
                - loss_right
                + 2 * intersection_phase_count
            )
            upper = np.minimum(
                np.minimum(cap_not_left, cap_not_right),
                overlap_bound,
            )

            # Branch and bound.  Include every label at the current cutoff
            # so ties cannot be lost.  If the fifth exact value is below
            # the cutoff, double the candidate window.
            window = 32
            evaluated = set()
            top_five = []
            while True:
                if window >= len(units):
                    cutoff = -1
                    candidates = np.arange(len(units))
                else:
                    cutoff = int(np.partition(upper, -window)[-window])
                    candidates = np.flatnonzero(upper >= cutoff)

                order = candidates[
                    np.lexsort((units[candidates], -upper[candidates]))
                ]
                for unit_index in order:
                    unit_index = int(unit_index)
                    if unit_index in evaluated:
                        continue
                    bound = int(upper[unit_index])
                    if len(top_five) == 5 and bound < top_five[-1][0]:
                        break

                    capacity = (
                        (residual & unit_at_least_one[unit_index]).bit_count()
                        + (
                            residual
                            & unit_at_least_two[unit_index]
                        ).bit_count()
                    )
                    require(
                        capacity <= bound,
                        "a branch-and-bound capacity exceeded its bound",
                    )
                    evaluated.add(unit_index)
                    top_five_insert(
                        top_five, (capacity, int(units[unit_index]))
                    )

                if len(top_five) == 5 and top_five[-1][0] >= cutoff:
                    break
                if window >= len(units):
                    break
                window = min(len(units), 2 * window)

            candidate_count = len(evaluated)
            exact_candidate_evaluations += candidate_count
            maximum_candidates = max(maximum_candidates, candidate_count)

            margin = residual_size - sum(
                capacity for capacity, _ in top_five
            )
            require(
                margin > 0,
                "a depth-(1,2,3) pair had nonpositive exact margin",
            )
            record = (left, right, residual_size, tuple(top_five))
            table_digest.update(repr((record, margin)).encode())

            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_record = record
                margin_minimizers = 1
            elif margin == best_margin:
                margin_minimizers += 1

            # Independent one-owner Ky-Fan certificate.
            coarse_upper = min(
                int(top_five_depth_one[i]),
                int(top_five_depth_two[j]),
            )
            coarse_margin = residual_size - coarse_upper
            if coarse_margin <= 0:
                coarse_nonpositive.append(
                    (
                        left,
                        right,
                        residual_size,
                        int(top_five_depth_one[i]),
                        int(top_five_depth_two[j]),
                        coarse_upper,
                        coarse_margin,
                        tuple(top_five),
                        margin,
                    )
                )
            elif (
                coarse_minimum_positive is None
                or coarse_margin < coarse_minimum_positive
            ):
                coarse_minimum_positive = coarse_margin
                coarse_minimum_positive_record = (
                    left,
                    right,
                    residual_size,
                    coarse_upper,
                )

            if (left, right) in requested_pairs:
                requested_records[(left, right)] = (record, margin)

    require(pair_count == 79092, "bad typed shallow-pair count")
    require(
        exact_candidate_evaluations == 940857,
        "bad exact candidate-evaluation count",
    )
    require(maximum_candidates == 27, "bad maximum candidate count")
    require(best_margin == 1608, "bad minimum exact margin")
    require(margin_minimizers == 1, "exact minimum was not unique")
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
        "bad unique minimum row",
    )
    require(
        table_digest.hexdigest()
        == "79b9b75f3732e47b43c2bba726906250ce0c1069b90d8952c33caa8a8364570f",
        "depth-(1,2,3) pair-table digest changed",
    )

    require(
        coarse_minimum_positive == 36,
        "bad minimum positive one-owner margin",
    )
    require(
        coarse_minimum_positive_record == (1, 1, 13310, 13274),
        "bad one-owner minimum-positive row",
    )
    require(
        coarse_nonpositive
        == [
            (
                1098,
                84,
                13310,
                13874,
                13328,
                13328,
                -18,
                (
                    (2396, 6),
                    (2330, 14275),
                    (2324, 14278),
                    (2296, 5),
                    (2288, 12),
                ),
                1676,
            )
        ],
        "bad sole one-owner exception",
    )

    # Independent definition-level controls on the full residue universe.
    all_residues = np.arange(N, dtype=np.int64)
    direct_universe = all_residues[
        (all_residues % P != 0)
        & (7 * norm_numerators(all_residues, N) > N)
    ]
    require(
        len(direct_universe) == 18830,
        "direct primitive guard-safe universe size changed",
    )

    for pair in ((799, 46), (1098, 84)):
        record, margin = requested_records[pair]
        direct = direct_residual(direct_universe, *pair)
        require(
            len(direct) == record[2],
            f"direct residual size changed for {pair}",
        )
        direct_rows = tuple(
            (direct_capacity(direct, label), label)
            for _, label in record[3]
        )
        require(
            direct_rows == record[3],
            f"direct top-five capacities changed for {pair}",
        )
        require(
            len(direct)
            - sum(capacity for capacity, _ in direct_rows)
            == margin,
            f"direct margin changed for {pair}",
        )

    average_milli = (
        1000 * exact_candidate_evaluations + pair_count // 2
    ) // pair_count
    require(average_milli == 11896, "bad rounded candidate average")

    print("THM-2207 scalar depth-(1,2,3) labelled guard-hole audit")
    print("arithmetic: exact integers; no assert dependence")
    print()
    print("depth-four primitive root geometry")
    print(
        f"  N={N} Q={Q} primitive_phases={len(phases)} "
        f"|U_N|={int(guard_counts.sum())}"
    )
    print(f"  guard_fibre_histogram={guard_histogram}")
    print(
        f"  unit_sign_classes={len(units)} "
        f"depth_one_sign_classes={len(depth_one_units)} "
        f"depth_two_sign_classes={len(depth_two_units)} "
        f"depth_three_sign_classes={len(deepest_units)}"
    )
    print("  shallow direct/image parity=1014+78 PASS")
    print("  depth-three primitive-layer safety=6/6 PASS")
    print()
    print("full labelled branch-and-bound scan")
    print(f"  typed_shallow_pairs={pair_count}")
    print(
        "  exact_candidate_evaluations="
        f"{exact_candidate_evaluations} "
        f"average={average_milli // 1000}.{average_milli % 1000:03d} "
        f"maximum={maximum_candidates}"
    )
    print(f"  unique_minimum_margin={best_margin}")
    print(
        "  unique_minimum_row="
        f"(depth_one={best_record[0]},depth_two={best_record[1]},"
        f"residual={best_record[2]})"
    )
    print(f"  unique_minimum_top_five={best_record[3]}")
    print(f"  pair_table_sha256={table_digest.hexdigest()}")
    print()
    print("independent one-owner Ky-Fan certificate")
    print("  coarse_positive_pairs=79091/79092")
    print(
        "  minimum_positive_coarse_margin=36 "
        "at (depth_one,depth_two)=(1,1)"
    )
    print(
        "  sole_coarse_exception=(1098,84) "
        "residual=13310 coarse_upper=13328 coarse_margin=-18"
    )
    exception_record, exception_margin = requested_records[(1098, 84)]
    print(f"  exception_exact_top_five={exception_record[3]}")
    print(f"  exception_exact_margin={exception_margin}")
    print()
    print("definition-level full-residue controls")
    print("  unique minimum row direct parity=PASS")
    print("  sole coarse exception direct parity=PASS")
    print()
    print("verdict: depth-(1,2,3) excluded")


if __name__ == "__main__":
    run()
