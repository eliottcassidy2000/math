#!/usr/bin/env python3
"""Exact common-sheet bank for the semantic root-zero overlap clutch.

This exploratory companion extends the single-sheet probe at commit
``ad375e51e``.  It inserts the full E3 mask and delayed Q_(3,{1,2}) word,
then enumerates all 9 x 9 target sheets common to the source and target
semantic points on canonical rail 8.  Source and target coefficients are
recomputed independently.  No endpoint-current or LRC(14) conclusion is
asserted.
"""

from collections import Counter
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCY = COMP / "lrc14_root_zero_overlap_clutch_20260728.py"
require(
    hashlib.sha256(
        DEPENDENCY.read_bytes().replace(bytes((13, 10)), bytes((10,)))
    ).hexdigest()
    == "e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8",
    "root-zero clutch dependency changed",
)

import lrc14_root_zero_overlap_clutch_20260728 as clutch


relative = clutch.relative
old = relative.private.old
P = relative.P
T = relative.lift.m.T
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
SOURCE_CARRY = clutch.SOURCE_STATE[2]
TARGET_CARRY = clutch.TARGET_STATE[2]
SOURCE_ROOT = clutch.SOURCE_STATE[3]
TARGET_ROOT = clutch.TARGET_STATE[3]


def intersect_danger(pieces, module, speed, period, lo, hi):
    danger = module.make_comb(speed, period, lo, hi)
    return tuple(old.intersect_weighted_union(pieces, danger))


def intersect_safe(pieces, module, speed, period, lo, hi):
    danger = module.make_comb(speed, period, lo, hi)
    return tuple(old.intersect_weighted_union(pieces, clutch.complement(danger)))


def restrict_e3(pieces, module):
    """Insert the complete unshifted E3 source mask."""
    pieces = intersect_safe(pieces, module, 1, 91, -13, 13)
    for speed in relative.semantic.UNITS:
        pieces = intersect_safe(pieces, module, speed, 182, -13, 13)
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[0], 182, -13, 13
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[1], 182, -13, 13
    )
    return intersect_danger(
        pieces, module, relative.semantic.BLOCKERS[2], 182, -13, 13
    )


def restrict_sheet(pieces, module, s, t):
    """Insert the four shifted safe factors defining U_(s,t)."""
    pieces = intersect_safe(
        pieces, module, relative.semantic.UNITS[0], 182,
        -14 * s - 13, -14 * s + 13,
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[1], 182,
        14 * s - 13, 14 * s + 13,
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.UNITS[1], 182,
        -14 * t - 13, -14 * t + 13,
    )
    return intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[2], 182,
        14 * t - 13, 14 * t + 13,
    )


def build_q3_pair_prefixes(module):
    """Delayed sector-zero prefixes with the full Q_(3,{1,2}) word."""
    word = relative.private.prior.sector_words(module)[0]
    by_clock = []
    for ell in range(7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        qell = module.intersect_comb(qell, module.C1, 182, -13, 13)
        by_h = []
        for h in range(P):
            pair = []
            for kappa in range(2):
                digit = old.sat.intersect_interval(
                    qell,
                    (2 * h + kappa) * T // (2 * P),
                    (2 * h + kappa + 1) * T // (2 * P),
                )
                pair.append(module.make_prefix(digit))
            by_h.append(tuple(pair))
        by_clock.append(tuple(by_h))
    return tuple(by_clock)


def scale_profile(profile, scalar):
    return tuple((scalar * value) % P for value in profile)


def profile_scalar(profile):
    require(profile is not None, "missing reduced profile")
    require(profile[0] and profile[1:] == (0, 0, 0, 0, 0),
            "profile left the nonzero constant line")
    return profile[0]


def rank_mod_p(matrix):
    work = [list(row) for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work))
             if work[index][column] % P),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column] % P, -1, P)
        work[row] = [(inverse * value) % P for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column] % P:
                continue
            scalar = work[index][column] % P
            work[index] = [
                (value - scalar * pivot_value) % P
                for value, pivot_value in zip(work[index], work[row])
            ]
        row += 1
    return row


def det3_mod_p(matrix, row_indices, column_indices):
    values = [
        [matrix[row][column] % P for column in column_indices]
        for row in row_indices
    ]
    return (
        values[0][0] * (
            values[1][1] * values[2][2] - values[1][2] * values[2][1]
        )
        - values[0][1] * (
            values[1][0] * values[2][2] - values[1][2] * values[2][0]
        )
        + values[0][2] * (
            values[1][0] * values[2][1] - values[1][1] * values[2][0]
        )
    ) % P


def main():
    module, _prefixes, _, _, rails, present, _starts = (
        relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = build_q3_pair_prefixes(module)
    target_pullback = clutch.shift_weighted(rails[8][3], -clutch.SHIFT)
    rail_pairs = clutch.intersect_weighted_profiles(
        rails[8][3], target_pullback
    )
    require(
        rail_pairs and all(row[2] == row[3] for row in rail_pairs),
        "rail-eight source/target weights stopped agreeing",
    )

    base = []
    for ell in range(7):
        source, target = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        base.append((restrict_e3(source, module), restrict_e3(target, module)))

    caches = tuple({} for _ell in range(7))
    rows = []
    for s in COMMON_S:
        for t in COMMON_T:
            source_vector = []
            target_vector = []
            source_masses = []
            target_masses = []
            for ell, (source_base, target_base) in enumerate(base):
                source = restrict_sheet(source_base, module, s, t)
                target = restrict_sheet(target_base, module, s, t)
                source_masses.append(sum((b - a) * w for a, b, w in source))
                target_masses.append(sum((b - a) * w for a, b, w in target))
                source_values = relative.private.delayed_carry_pair(
                    source, pair_prefixes[ell][6], caches[ell]
                )
                target_values = relative.private.delayed_carry_pair(
                    target, pair_prefixes[ell][6], caches[ell]
                )
                source_vector.append(source_values[SOURCE_CARRY][1])
                target_vector.append(target_values[TARGET_CARRY][1])

            source_vector = tuple(source_vector)
            target_vector = tuple(target_vector)
            source_content = all(value % clutch.CONTENT == 0
                                 for value in source_vector)
            target_content = all(value % clutch.CONTENT == 0
                                 for value in target_vector)
            source_profile = None
            target_profile = None
            source_unit = target_unit = False
            covariance = False
            if source_content:
                source_profile = clutch.normalized_profile(
                    source_vector, SOURCE_ROOT
                )[1]
                source_unit = relative.private.is_unit(
                    source_vector, SOURCE_ROOT, clutch.CONTENT
                )
            if target_content:
                target_profile = clutch.normalized_profile(
                    target_vector, TARGET_ROOT
                )[1]
                target_unit = relative.private.is_unit(
                    target_vector, TARGET_ROOT, clutch.CONTENT
                )
            if source_profile is not None and target_profile is not None:
                covariance = (
                    target_profile == scale_profile(source_profile, 7)
                    and scale_profile(source_profile, pow(SOURCE_CARRY, -1, P))
                    == scale_profile(target_profile, pow(TARGET_CARRY, -1, P))
                )
            rows.append({
                "sheet": (s, t),
                "source_vector": source_vector,
                "target_vector": target_vector,
                "source_mass": tuple(source_masses),
                "target_mass": tuple(target_masses),
                "source_profile": source_profile,
                "target_profile": target_profile,
                "source_unit": source_unit,
                "target_unit": target_unit,
                "covariance": covariance,
            })

    require(
        len(rows) == len(COMMON_S) * len(COMMON_T),
        "common-sheet bank cardinality changed",
    )
    profile_pairs = Counter(
        (row["source_profile"], row["target_profile"]) for row in rows
    )
    unit_count = sum(
        row["source_unit"] and row["target_unit"] for row in rows
    )
    covariance_count = sum(row["covariance"] for row in rows)
    constant_tail_count = sum(
        len(set(row["source_vector"][1:])) == 1
        and len(set(row["target_vector"][1:])) == 1
        and row["source_vector"][0] == row["target_vector"][0] == 0
        for row in rows
    )
    carry_ratio = TARGET_CARRY * pow(SOURCE_CARRY, -1, P) % P
    carry_difference = (TARGET_CARRY - SOURCE_CARRY) % P
    require(carry_ratio == carry_difference == 7,
            "canonical carry gain/difference stopped equalling seven")
    gain_census = Counter()
    shear_reference_census = Counter()
    covariance_sheets = []
    scalar_rows = []
    for row in rows:
        source_scalar = profile_scalar(row["source_profile"])
        target_scalar = profile_scalar(row["target_profile"])
        gain = target_scalar * pow(source_scalar, -1, P) % P
        reference = (
            (target_scalar - source_scalar) * pow(carry_difference, -1, P)
        ) % P
        require(
            target_scalar == (source_scalar + carry_difference * reference) % P,
            "sheetwise unipotent shear identity failed",
        )
        gain_census[gain] += 1
        shear_reference_census[reference] += 1
        if gain == carry_ratio:
            covariance_sheets.append(row["sheet"])
        scalar_rows.append(
            (row["sheet"], source_scalar, target_scalar, gain, reference)
        )

    expected_profile_pairs = Counter({
        ((5, 0, 0, 0, 0, 0), (9, 0, 0, 0, 0, 0)): 10,
        ((12, 0, 0, 0, 0, 0), (3, 0, 0, 0, 0, 0)): 19,
        ((9, 0, 0, 0, 0, 0), (8, 0, 0, 0, 0, 0)): 2,
        ((4, 0, 0, 0, 0, 0), (12, 0, 0, 0, 0, 0)): 9,
        ((7, 0, 0, 0, 0, 0), (12, 0, 0, 0, 0, 0)): 2,
        ((8, 0, 0, 0, 0, 0), (6, 0, 0, 0, 0, 0)): 5,
        ((11, 0, 0, 0, 0, 0), (6, 0, 0, 0, 0, 0)): 2,
        ((7, 0, 0, 0, 0, 0), (5, 0, 0, 0, 0, 0)): 20,
        ((2, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0)): 12,
    })
    require(
        constant_tail_count == unit_count == 81
        and profile_pairs == expected_profile_pairs
        and gain_census == Counter({3: 9, 4: 5, 7: 22, 10: 41, 11: 4})
        and shear_reference_census
        == Counter({3: 11, 8: 29, 9: 25, 10: 2, 11: 14})
        and scalar_rows[5] == ((0, 8), 12, 3, 10, 8),
        "common-sheet exact census changed",
    )

    def matrix(index):
        table = {(s, t): row[index] for (s, t), *row in scalar_rows}
        # row[index] addresses the tail after the sheet: source=0, target=1,
        # gain=2, reference=3.
        return tuple(tuple(table[s, t] for t in COMMON_T) for s in COMMON_S)

    source_matrix = matrix(0)
    target_matrix = matrix(1)
    gain_matrix = matrix(2)
    reference_matrix = matrix(3)
    matrix_ranks = tuple(
        rank_mod_p(table) for table in (
            source_matrix, target_matrix, gain_matrix, reference_matrix
        )
    )
    column_blocks = (
        (3, 4, 5, 6, 7),
        (8, 9),
        (10, 11),
    )
    minor_witnesses = (
        ("source", (0, 1, 2), (0, 5, 7), 9),
        ("target", (0, 1, 2), (0, 5, 7), 4),
        ("gain", (0, 1, 2), (0, 5, 7), 5),
        ("reference", (0, 2, 3), (0, 5, 7), 6),
    )
    tables = (source_matrix, target_matrix, gain_matrix, reference_matrix)
    require(
        matrix_ranks == (3, 3, 3, 3)
        and all(
            all(
                all(row[0] == row[index] for index in range(1, 5))
                and row[5] == row[6]
                and row[7] == row[8]
                for row in table
            )
            for table in tables
        )
        and tuple(
            det3_mod_p(table, rows_index, columns_index)
            for table, (_name, rows_index, columns_index, _expected)
            in zip(tables, minor_witnesses)
        ) == tuple(item[3] for item in minor_witnesses)
        and tuple(
            tuple(COMMON_S[index] for index in rows_index)
            for _name, rows_index, _columns_index, _expected in minor_witnesses
        ) == ((0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 2, 3))
        and tuple(COMMON_T[index] for index in (0, 5, 7)) == (3, 8, 10),
        "three-block rank atlas changed",
    )

    print("LRC14 SEMANTIC ROOT-ZERO COMMON-SHEET BANK")
    print("status=FINITE-EXACT CANDIDATE AWAITING INDEPENDENT AUDIT")
    print(f"common_s={COMMON_S}")
    print(f"common_t={COMMON_T}")
    print(f"sheets={len(rows)}")
    print(f"constant_tail_pairs={constant_tail_count}")
    print(f"private_unit_pairs={unit_count}")
    print(f"carry_ratio={carry_ratio} carry_difference={carry_difference}")
    print(f"uniform_gain7_pairs={covariance_count}")
    print(f"profile_pair_census={tuple(sorted(profile_pairs.items(), key=str))}")
    print(f"gain_census={tuple(sorted(gain_census.items()))}")
    print(f"shear_reference_census={tuple(sorted(shear_reference_census.items()))}")
    print(f"gain7_sheets={tuple(covariance_sheets)}")
    print("first_gain7_failure=((0,8),source=12,target=3,gain=10,reference=8)")
    print(f"source_scalar_matrix={source_matrix}")
    print(f"target_scalar_matrix={target_matrix}")
    print(f"gain_matrix={gain_matrix}")
    print(f"shear_reference_matrix={reference_matrix}")
    print(f"matrix_ranks=(source,target,gain,reference):{matrix_ranks}")
    print(f"three_separable_target_blocks={column_blocks}")
    print(
        "nonzero_3x3_minors="
        "source[s=0,1,2;t=3,8,10]:9,"
        "target[s=0,1,2;t=3,8,10]:4,"
        "gain[s=0,1,2;t=3,8,10]:5,"
        "reference[s=0,2,3;t=3,8,10]:6"
    )
    print(
        "one_dimensional_hostile="
        f"7^12:{pow(7, 12, P)},7^13:{pow(7, 13, P)},7^2:{pow(7, 2, P)}"
    )
    print(
        "SCOPE: full E3 + delayed Q_(3,{1,2}) on canonical rail 8; "
        "the sheetwise U(7) reference is pair-derived, not a physical packet; "
        "no endpoint current, additive-odometer identification, or LRC14 conclusion"
    )


if __name__ == "__main__":
    main()
