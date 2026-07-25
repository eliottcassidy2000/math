#!/usr/bin/env python3
"""Exact audit for THM-2204: scalar blocker profile (2,2,3).

At N=13^4, the two depth-two blocker masks are constant on the thirteen
roots above each primitive phase modulo Q=13^3.  For every one of their
3,081 unordered sign-class pairs, this script computes the five largest
conditional capacities of all 13,182 unit sign classes.  A unit capacity is
represented by two phase bitsets: h>=1 and h>=2, where h is the exact number
of roots that are simultaneously guard-safe and unit-dangerous.

The script also verifies the exact all-depth thirteen-lift family-sum law
for every one of the 3,081*1,014 pair/family cases.  Arithmetic is exact;
``require`` remains active under ``python -O``.
"""

from collections import defaultdict
import hashlib

import numpy as np


P = 13
N = P**4
Q = P**3


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


def direct_danger_count(universe, coefficient):
    """Definition-level strict danger count on a named residue universe."""
    return int(
        np.sum(
            14 * norm_numerators(coefficient * universe, N) < N
        )
    )


def run():
    phases = np.arange(1, Q, dtype=np.int64)
    phases = phases[phases % P != 0]
    units = sign_classes(N)
    base_units = sign_classes(Q)
    shallow_units = sign_classes(P**2)
    sheets = np.arange(P, dtype=np.int64)

    require(len(phases) == 2028, "bad primitive phase count")
    require(len(units) == 13182, "bad depth-four unit sign-class count")
    require(len(base_units) == 1014,
            "bad depth-three unit sign-class count")
    require(len(shallow_units) == 78,
            "bad depth-two blocker sign-class count")
    require(N % 7 != 0 and N % 14 != 0 and Q % 14 != 0,
            "a strict torsion endpoint appeared")

    # Direct root geometry over every primitive image phase.
    roots = phases[:, None] + Q * sheets[None, :]
    guard = 7 * norm_numerators(roots, N) > N
    guard_counts = guard.sum(axis=1)
    require(
        tuple(sorted(set(map(int, guard_counts)))) == (9, 10),
        "guard roots did not have nine/ten-sheet fibres",
    )
    guard_histogram = tuple(
        (size, int(np.sum(guard_counts == size))) for size in (9, 10)
    )
    require(guard_histogram == ((9, 1450), (10, 578)),
            "bad depth-four guard-fibre histogram")
    require(int(guard_counts.sum()) == 18830,
            "bad primitive guard-safe annulus size")
    guard_nine = rows_to_bitsets((guard_counts == 9)[None, :])[0]

    # A depth-two mask at N is a depth-one mask on the image phase r/Q.
    shallow_active_rows = (
        14
        * norm_numerators(
            (P * shallow_units[:, None]) * phases[None, :],
            Q,
        )
        < Q
    )
    shallow_active = rows_to_bitsets(shallow_active_rows)
    phase_full = (1 << len(phases)) - 1

    # Hostile all-class parity: direct depth-two evaluation on all thirteen
    # roots agrees with the reduced, fibre-constant image mask.
    for start in range(0, len(shallow_units), 13):
        shallow = shallow_units[start : start + 13]
        direct = (
            14
            * norm_numerators(
                (P**2 * shallow[:, None, None]) * roots[None, :, :],
                N,
            )
            < N
        )
        reduced = shallow_active_rows[start : start + 13, :, None]
        require(
            np.array_equal(direct, np.broadcast_to(reduced, direct.shape)),
            "depth-two direct/image fibre parity failed",
        )

    # For each unit q, h_q(r) is the number of roots that are both
    # guard-safe and q-dangerous.  The root window has one or two sheets, so
    # h belongs to {0,1,2}.  Its conditional capacity on a phase set P is
    # popcount(P & {h>=1}) + popcount(P & {h>=2}).
    unit_at_least_one = []
    unit_at_least_two = []
    for start in range(0, len(units), 24):
        batch = units[start : start + 24]
        residues = (
            batch[:, None, None] * roots[None, :, :]
        ) % N
        dangerous = 14 * np.minimum(residues, N - residues) < N
        counts = np.sum(dangerous & guard[None, :, :], axis=2)
        require(int(counts.max()) <= 2,
                "a unit root window had more than two sheets")
        unit_at_least_one.extend(rows_to_bitsets(counts >= 1))
        unit_at_least_two.extend(rows_to_bitsets(counts >= 2))
    unit_at_least_one = tuple(unit_at_least_one)
    unit_at_least_two = tuple(unit_at_least_two)
    require(
        len(unit_at_least_one) == len(unit_at_least_two) == len(units),
        "bad unit phase-bitset count",
    )

    # The 13 sign-class lifts over each depth-three sign class.
    families = defaultdict(list)
    unit_base = []
    for index, q in enumerate(map(int, units)):
        residue = q % Q
        base = min(residue, (-residue) % Q)
        families[base].append(index)
        unit_base.append(base)
    require(len(families) == 1014, "bad thirteen-lift family count")
    require(set(map(len, families.values())) == {13},
            "a lift family did not have thirteen sign classes")
    require(set(families) == set(map(int, base_units)),
            "lift-family bases did not match the depth-three sign classes")

    # Base danger masks enter the exact family-sum identity
    #
    #   sum_l C_l = 2 W - W_a,
    #
    # where W is residual sheet-time mass and W_a is the same mass over
    # base phases dangerous for a.
    base_danger = rows_to_bitsets(
        14
        * norm_numerators(base_units[:, None] * phases[None, :], Q)
        < Q
    )
    base_to_index = {int(a): i for i, a in enumerate(base_units)}

    table_digest = hashlib.sha256()
    best_margin = None
    best_record = None
    margin_minimizers = 0
    pair_count = 0
    diagonal_count = 0
    family_identity_checks = 0
    maximum_lift_spread = None
    maximum_lift_spread_record = None

    for i, left in enumerate(map(int, shallow_units)):
        for j in range(i, len(shallow_units)):
            right = int(shallow_units[j])
            pair_count += 1
            diagonal_count += i == j
            residual = phase_full ^ (
                shallow_active[i] | shallow_active[j]
            )
            residual_phase_count = residual.bit_count()
            residual_size = (
                10 * residual_phase_count
                - (residual & guard_nine).bit_count()
            )

            top_five = []
            family_sums = defaultdict(int)
            family_mins = {}
            family_maxs = {}
            for index, q in enumerate(map(int, units)):
                capacity = (
                    (residual & unit_at_least_one[index]).bit_count()
                    + (residual & unit_at_least_two[index]).bit_count()
                )
                top_five_insert(top_five, (capacity, q))
                base = unit_base[index]
                family_sums[base] += capacity
                family_mins[base] = min(
                    family_mins.get(base, capacity), capacity
                )
                family_maxs[base] = max(
                    family_maxs.get(base, capacity), capacity
                )

            margin = residual_size - sum(
                capacity for capacity, _ in top_five
            )
            require(margin > 0,
                    "a depth-(2,2,3) pair had nonpositive margin")
            record = (left, right, residual_size, tuple(top_five))
            table_digest.update(repr((record, margin)).encode())
            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_record = record
                margin_minimizers = 1
            elif margin == best_margin:
                margin_minimizers += 1

            # Verify the family-sum law for every pair/base combination.
            for base, total in family_sums.items():
                danger = base_danger[base_to_index[base]] & residual
                weighted_danger = (
                    10 * danger.bit_count()
                    - (danger & guard_nine).bit_count()
                )
                expected = 2 * residual_size - weighted_danger
                require(total == expected,
                        "thirteen-lift family-sum identity failed")
                family_identity_checks += 1

                spread = family_maxs[base] - family_mins[base]
                if (
                    maximum_lift_spread is None
                    or spread > maximum_lift_spread
                ):
                    maximum_lift_spread = spread
                    maximum_lift_spread_record = (
                        left,
                        right,
                        base,
                        family_mins[base],
                        family_maxs[base],
                    )

    require(pair_count == 3081, "bad unordered shallow-pair count")
    require(diagonal_count == 78, "diagonal shallow pairs were lost")
    require(family_identity_checks == 3124134,
            "bad family-sum identity check count")
    require(best_margin == 1830, "bad minimum depth-(2,2,3) margin")
    require(margin_minimizers == 1,
            "depth-(2,2,3) minimum was not unique")
    require(
        best_record
        == (
            11,
            79,
            14484,
            (
                (2614, 1858),
                (2614, 1860),
                (2504, 929),
                (2504, 930),
                (2418, 6),
            ),
        ),
        "bad unique worst depth-(2,2,3) row",
    )
    require(
        table_digest.hexdigest()
        == "b82b4f22fe2242197fc8ee587a99acd026120d52001d9093e8b465557cc57024",
        "depth-(2,2,3) table digest changed",
    )
    require(maximum_lift_spread == 2460,
            "bad maximum within-family lift spread")
    require(
        maximum_lift_spread_record == (4, 4, 1098, 0, 2460),
        "bad first maximum-spread witness",
    )

    # Independent direct-torsion hostile control for the unique worst row.
    all_residues = np.arange(N, dtype=np.int64)
    direct_universe = all_residues[
        (all_residues % P != 0)
        & (7 * norm_numerators(all_residues, N) > N)
    ]
    direct_residual = direct_universe[
        (14 * norm_numerators(P**2 * 11 * direct_universe, N) >= N)
        & (14 * norm_numerators(P**2 * 79 * direct_universe, N) >= N)
    ]
    require(len(direct_universe) == 18830,
            "direct hostile universe size changed")
    require(len(direct_residual) == best_record[2],
            "direct hostile residual size did not match phase audit")
    direct_hostile_capacities = tuple(
        (direct_danger_count(direct_residual, label), label)
        for _, label in best_record[3]
    )
    require(direct_hostile_capacities == best_record[3],
            "direct hostile capacities did not match phase audit")
    require(
        len(direct_residual)
        - sum(capacity for capacity, _ in direct_hostile_capacities)
        == best_margin,
        "direct hostile margin did not match phase audit",
    )

    # Explicit concentration witness.  The thirteen sign-class lifts above
    # base a=1098 have a fixed total but range from zero to 2460 on the
    # diagonal shallow pair (4,4).
    shallow_four_index = tuple(map(int, shallow_units)).index(4)
    concentration_residual = (
        phase_full ^ shallow_active[shallow_four_index]
    )
    concentration_indices = families[1098]
    concentration_rows = tuple(
        (
            int(units[index]),
            (concentration_residual & unit_at_least_one[index]).bit_count()
            + (concentration_residual & unit_at_least_two[index]).bit_count(),
        )
        for index in concentration_indices
    )
    require(
        concentration_rows
        == (
            (1098, 2460),
            (1099, 2454),
            (3295, 2456),
            (3296, 2452),
            (5492, 2454),
            (5493, 2452),
            (7689, 2450),
            (7690, 2456),
            (9886, 2456),
            (9887, 2450),
            (12083, 2452),
            (12084, 2452),
            (14280, 0),
        ),
        "bad lift-concentration witness",
    )
    require(sum(value for _, value in concentration_rows) == 29444,
            "bad lift-concentration family sum")

    print("THM-2204 scalar depth-(2,2,3) thirteen-lift capacity audit")
    print("arithmetic: exact integers; no floats; no assert dependence")
    print()
    print("depth-four primitive root geometry")
    print(
        f"  N={N} Q={Q} primitive_phases={len(phases)} "
        f"|U_N|={int(guard_counts.sum())}"
    )
    print(f"  guard_fibre_histogram={guard_histogram}")
    print(
        f"  unit_sign_classes={len(units)} "
        f"base_unit_sign_classes={len(base_units)} "
        f"depth_two_sign_classes={len(shallow_units)}"
    )
    print("  depth-two direct/image parity=78/78 PASS")
    print()
    print("complete depth-(2,2,3) capacity certificate")
    print(
        f"  unordered_shallow_pairs={pair_count} "
        f"including_diagonal={diagonal_count}"
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
    print(f"  full_table_digest={table_digest.hexdigest()}")
    print("  independent direct hostile row parity=PASS")
    print()
    print("thirteen-lift recursive law and obstruction")
    print(
        f"  lift_families={len(families)} size=13 "
        f"family_sum_checks={family_identity_checks} PASS"
    )
    print("  exact_identity=sum_l C_l=2W-W_a")
    print(
        "  maximum_within_family_spread=2460 "
        "first_witness=(shallow_pair=(4,4),base=1098,min=0,max=2460)"
    )
    print(f"  witness_(sign_label,capacity)={concentration_rows}")
    print("  witness_family_sum=29444")
    print(
        "  consequence: family sums alone do not control lifted "
        "top-five order statistics"
    )
    print()
    print("CONCLUSION")
    print("  actual blocker valuation profile (2,2,3) is empty")
    print("  profiles (1,1,3) and (1,2,3) remain open")
    print()
    print("reproduce:")
    print(
        "  python3 "
        "04-computation/"
        "lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py"
    )
    print(
        "  python3 -O "
        "04-computation/"
        "lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py"
    )


if __name__ == "__main__":
    run()
