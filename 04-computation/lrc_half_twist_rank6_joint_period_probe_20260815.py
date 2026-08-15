#!/usr/bin/env python3
"""Exact companion for the lower-base-free half-twist joint-period closure.

The proof core classifies the sector in which 8,9,10,12 do not divide Q and
the selected quotient orders have lcm Q.  Adding one breaker owner to scaled
rank-four/rank-five atoms then gives the unrestricted primitive half-twist
support.  Joint period and augmented primitivity are kept distinct family by
family when Q is odd.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

from lrc_zero_mode_cochain_rank6_support_probe_20260815 import (
    anchor_complement_record,
    danger_mask,
    exact_order_masks,
    half_block_count,
)


ROOT = Path(__file__).resolve().parents[1]
LOWER_BASES = (8, 9, 10, 12)
EXPECTED_SEMANTIC_DIGEST = "4ac644546e4f81631b5c12a404779397c5c63da67f3e4dd5c1b88da1fa8beda1"
PINNED = (
    (
        "THM-3414",
        ROOT / "01-canon/theorems/THM-3414-fixed-zero-six-owner-base-classification.md",
        "5568a4e93bc4478566335e2722c488c999797462eeb7c95af364b20dba41e998",
    ),
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def lower_base_free(value):
    return not any(value % base == 0 for base in LOWER_BASES)


def critical_orders(threshold, last):
    require(Fraction(last + 7, 7 * (last + 1)) < threshold, (threshold, last))
    return tuple(
        order
        for order in range(2, last + 1)
        if lower_base_free(order)
        and Fraction(half_block_count(order), order) >= threshold
    )


def complement_table(anchor, quota, last):
    orders = critical_orders(quota, last)
    return tuple(
        (order, anchor_complement_record(anchor, order)[4])
        for order in orders
    )


def core_score(prime, order, mask):
    size = mask.bit_count()
    if order % prime:
        return Fraction(6 * size, order)
    fixed = (prime - 1) // 2
    fixed_count = sum(
        1 for sheet in range(order)
        if sheet % prime == fixed and mask >> sheet & 1
    )
    if prime == 11:
        return Fraction(11 * (size + fixed_count), 2 * order)
    require(prime == 23, prime)
    return Fraction(23 * (size + fixed_count), 4 * order)


def maximal_core_score(prime, order):
    masks = exact_order_masks(order, order)
    require(masks, (prime, order, "empty full-order bank"))
    return max(core_score(prime, order, mask) for mask in masks)


def core_score_audit():
    eleven_rows = tuple((order, maximal_core_score(11, order)) for order in (11, 22, 33, 44, 55, 66))
    twenty_three_rows = tuple((order, maximal_core_score(23, order)) for order in (23, 46, 69))
    require(
        eleven_rows
        == ((11, Fraction(1)), (22, Fraction(1)), (33, Fraction(1)),
            (44, Fraction(3, 4)), (55, Fraction(4, 5)), (66, Fraction(5, 6))),
        eleven_rows,
    )
    require(
        twenty_three_rows
        == ((23, Fraction(1)), (46, Fraction(3, 4)), (69, Fraction(5, 6))),
        twenty_three_rows,
    )

    eleven_finite = tuple(
        (order, maximal_core_score(11, order))
        for order in range(2, 77)
        if lower_base_free(order) and exact_order_masks(order, order)
    )
    twenty_three_finite = tuple(
        (order, maximal_core_score(23, order))
        for order in range(2, 75)
        if lower_base_free(order) and exact_order_masks(order, order)
    )
    eleven_above = tuple(order for order, score in eleven_finite if score > 1)
    eleven_equal = tuple(order for order, score in eleven_finite if score == 1)
    twenty_three_above = tuple(order for order, score in twenty_three_finite if score > 1)
    twenty_three_equal = tuple(order for order, score in twenty_three_finite if score == 1)
    require(eleven_above == (3, 5, 15, 17, 23, 29), eleven_above)
    require(eleven_equal == (11, 22, 33), eleven_equal)
    require(twenty_three_above == (3, 5, 11, 15, 17, 22, 29), twenty_three_above)
    require(twenty_three_equal == (23,), twenty_three_equal)

    # The balanced-interval tails used in the proof.
    for order in range(77, 77 + 11 * 100, 11):
        h = half_block_count(order)
        require(11 * (h + (h + 10) // 11) < 2 * order, (11, order, h))
    for order in range(92, 92 + 23 * 100, 23):
        h = half_block_count(order)
        require(23 * (h + (h + 22) // 23) < 4 * order, (23, order, h))

    return (
        eleven_rows,
        twenty_three_rows,
        eleven_above,
        eleven_equal,
        twenty_three_above,
        twenty_three_equal,
    )


def equality_types(q, allowed_orders):
    grouped = set()
    for residue in range(1, 2 * q):
        if residue % q == 0:
            continue
        order = q // gcd(q, residue)
        if order not in allowed_orders:
            continue
        mask = danger_mask(q, residue, 1)
        if mask:
            grouped.add((mask, order))
    return tuple(sorted(grouped, key=lambda item: (item[1], item[0])))


def equality_cover_gate(q, allowed_orders):
    types = equality_types(q, allowed_orders)
    full = (1 << q) - 1
    covers = 0
    periods = []
    for packet in combinations(types, 6):
        joined = 0
        for mask, _ in packet:
            joined |= mask
        if joined == full:
            covers += 1
            periods.append(lcm(*(order for _, order in packet)))
    return len(types), covers, tuple(sorted(periods))


def witness_record(q, epsilon, residues):
    full = (1 << q) - 1
    joined = 0
    for residue in residues:
        joined |= danger_mask(q, residue, epsilon)
    require(joined == full, (q, epsilon, residues, joined, full))
    orders = tuple(sorted(q // gcd(q, residue) for residue in residues))
    period = lcm(*orders)
    modulus = q if epsilon == 0 else 2 * q
    return q, epsilon, residues, orders, period, gcd(modulus, *residues)


def lower_base_pullback_controls():
    atoms = {
        8: (1, 3, 5, 7),
        9: (1, 5, 6, 7),
        10: (1, 3, 4, 7, 9),
        12: (1, 5, 7, 8, 11),
    }
    rows = []
    for base, residues in atoms.items():
        for multiplier in range(1, 41):
            q = base * multiplier
            lifted = tuple(multiplier * residue for residue in residues)
            if multiplier > 1:
                lifted = (1,) + lifted
            record = witness_record(q, 1, lifted)
            require(len(lifted) <= 6 and record[4:] == (q, 1), (base, multiplier, record))
            rows.append((base, multiplier, len(lifted), record[4], record[5]))
    return tuple(rows)


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    branch_rows = (
        (15, Fraction(4, 25), complement_table(15, Fraction(4, 25), 50)),
        (17, Fraction(14, 85), complement_table(17, Fraction(14, 85), 39)),
        (29, Fraction(24, 145), complement_table(29, Fraction(24, 145), 37)),
    )
    expected_branch_rows = (
        (15, Fraction(4, 25), ((3, Fraction(4, 15)), (5, Fraction(2, 15)),
            (11, Fraction(8, 55)), (15, Fraction(2, 15)), (17, Fraction(12, 85)),
            (22, Fraction(8, 55)), (23, Fraction(16, 115)), (25, Fraction(2, 15)),
            (29, Fraction(4, 29)), (31, Fraction(4, 31)), (37, Fraction(24, 185)),
            (43, Fraction(28, 215)))),
        (17, Fraction(14, 85), ((3, Fraction(14, 51)), (5, Fraction(14, 85)),
            (11, Fraction(28, 187)), (15, Fraction(14, 85)), (17, Fraction(2, 17)),
            (22, Fraction(28, 187)), (23, Fraction(56, 391)), (29, Fraction(70, 493)))),
        (29, Fraction(24, 145), ((3, Fraction(8, 29)), (5, Fraction(24, 145)),
            (11, Fraction(48, 319)), (15, Fraction(24, 145)), (17, Fraction(72, 493)),
            (22, Fraction(48, 319)), (23, Fraction(96, 667)), (29, Fraction(4, 29)))),
    )
    require(branch_rows == expected_branch_rows, branch_rows)

    core_rows = core_score_audit()
    q33_gate = equality_cover_gate(33, (11, 22, 33))
    q66_gate = equality_cover_gate(66, (11, 22, 33))
    require(q33_gate == (26, 1, (11,)), q33_gate)
    require(q66_gate[0:2] == (36, 13), q66_gate)
    require(set(q66_gate[2]) <= {11, 22} and 66 not in q66_gate[2], q66_gate)

    mixed_gaps = tuple(
        (a, b, Fraction(344 * a * b - 207 * a - 121 * b + 253, 10879))
        for a in range(1, 6)
        for b in range(1, 6 - a + 1)
    )
    require(min(gap for _, _, gap in mixed_gaps) == Fraction(269, 10879), mixed_gaps)

    positives = tuple(
        witness_record(q, 1, residues)
        for q, residues in (
            (11, (1, 2, 3, 5, 7, 9)),
            (15, (1, 4, 6, 7, 8, 10)),
            (22, (1, 4, 9, 10, 13, 21)),
            (23, (1, 4, 5, 7, 9, 11)),
            (25, (1, 9, 10, 11, 19, 21)),
        )
    )
    require(all(row[4] == row[0] and row[5] == 1 for row in positives), positives)
    lower_pullbacks = lower_base_pullback_controls()
    q27_control = witness_record(27, 1, (1, 3, 15, 18, 21))
    require(q27_control[4:] == (27, 1), q27_control)
    fixed_zero_q15 = witness_record(15, 0, (1, 2, 3, 4, 5, 7))
    require(fixed_zero_q15[4:] == (15, 1), fixed_zero_q15)

    hostiles = (
        witness_record(33, 1, (3, 6, 9, 15, 21, 27)),
        witness_record(46, 1, (2, 8, 10, 14, 18, 22)),
    )
    require(hostiles[0][4:] == (11, 3), hostiles[0])
    require(hostiles[1][4:] == (23, 2), hostiles[1])
    odd_q_parity_hostile = witness_record(15, 1, (2, 4, 6, 8, 10, 14))
    require(odd_q_parity_hostile[4:] == (15, 2), odd_q_parity_hostile)

    semantic_surface = (
        branch_rows,
        core_rows,
        q33_gate,
        q66_gate,
        mixed_gaps,
        positives,
        lower_pullbacks,
        q27_control,
        fixed_zero_q15,
        hostiles,
        odd_q_parity_hostile,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("Half-twist rank-six lower-base-free joint-period closure")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=PROVED_ELEMENTARY_PLUS_VERIFIED_EXACT_PLUS_INDEPENDENTLY_AUDITED;lower_base_free_joint_period_closure;all_Q_primitive_half_cap6_support;no_LRC14_decrement")
    print(f"maximal_15_17_29_complement_rows={branch_rows}")
    print(f"weighted_11_23_core_scores={core_rows}")
    print(f"candidate_type_gates_Q33_Q66=(types,covers,joint_periods)={(q33_gate,q66_gate)}")
    print(f"mixed_11x23_CRT_gaps=(rows,min)={(mixed_gaps,min(gap for _,_,gap in mixed_gaps))}")
    print(f"positive_half_joint_period_atoms={positives}")
    print(f"lower_base_pullback_breaker_controls=(rows,first,last)={(len(lower_pullbacks),lower_pullbacks[0],lower_pullbacks[-1])}")
    print(f"Q27_primitive_five_cover_control={q27_control}")
    print(f"fixed_zero_lower_base_free_atom={fixed_zero_q15}")
    print(f"pullback_hostiles_Q33_Q46={hostiles}")
    print(f"odd_Q_joint_period_without_parity_breaker_hostile={odd_q_parity_hostile}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=if_8,9,10,12_do_not_divide_Q_then_half_cap6_cover_with_joint_period_Q_iff_Q_in_11,15,22,23,25;primitive_half_cap6_iff_lower_base_divisibility_or_one_of_five_points;fixed_zero_lower_base_free_iff_Q15;family_level_joint_period_needs_parity_sidecar_for_odd_Q;LRC14_open")


if __name__ == "__main__":
    main()
