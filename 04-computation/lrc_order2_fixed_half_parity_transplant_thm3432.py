#!/usr/bin/env python3
"""Exact controls for THM-3432's order-two parity transplant.

The proof is algebraic.  This independent finite bank checks the two parity
charts at several strict arc widths, the order/parity formula, the compatible
unit action, direct minimum-cover searches on the named boundary moduli, and
the augmented-primitivity defect at R=15.  It is standard-library only.
"""

from __future__ import annotations

import ast
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


MAX_IDENTITY_R = 96
THRESHOLDS = ((1, 14), (1, 7), (2, 7), (1, 3), (1, 2))
RANK_CONTROLS = {
    8: (1, 3, 5, 7),
    9: (1, 5, 6, 7),
    10: (1, 3, 4, 7, 9),
    11: (1, 2, 3, 5, 7, 9),
    12: (1, 5, 7, 8, 11),
    13: (1, 2, 3, 5, 7, 9, 11),
    15: (1, 4, 6, 7, 8, 10),
    23: (1, 4, 5, 7, 9, 11),
    25: (1, 9, 10, 11, 19, 21),
    29: (1, 5, 7, 8, 12, 13, 22),
}
EXPECTED_RANKS = {8: 4, 9: 4, 10: 5, 11: 6, 12: 5, 13: 7,
                  15: 6, 23: 6, 25: 6, 29: 7}


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def in_arc(phase, modulus, numerator=1, denominator=14):
    phase %= modulus
    return denominator * min(phase, modulus - phase) < numerator * modulus


def fixed_mask(modulus, residue, numerator=1, denominator=14):
    mask = 0
    for sheet in range(modulus):
        if in_arc(residue * sheet, modulus, numerator, denominator):
            mask |= 1 << sheet
    return mask


def half_mask(modulus, residue, numerator=1, denominator=14):
    ambient = 2 * modulus
    mask = 0
    for sheet in range(modulus):
        if in_arc(residue * (2 * sheet + 1), ambient, numerator, denominator):
            mask |= 1 << sheet
    return mask


def extract_parity(mask, modulus, parity):
    answer = 0
    for sheet in range(modulus):
        if mask & (1 << (2 * sheet + parity)):
            answer |= 1 << sheet
    return answer


def pullback_mask(mask, modulus, permutation):
    answer = 0
    for sheet in range(modulus):
        if mask & (1 << permutation(sheet)):
            answer |= 1 << sheet
    return answer


def additive_order(modulus, residue):
    return modulus // gcd(modulus, residue)


def quotient_order(modulus, residue):
    return modulus // gcd(modulus, residue)


def family_union(mask_function, modulus, residues):
    answer = 0
    for residue in residues:
        answer |= mask_function(modulus, residue)
    return answer


def unique_bank(mask_function, modulus, coefficient_modulus, excluded):
    by_mask = {}
    for residue in range(1, coefficient_modulus):
        if residue in excluded:
            continue
        mask = mask_function(modulus, residue)
        if mask:
            by_mask.setdefault(mask, residue)
    return tuple(sorted((mask, residue) for mask, residue in by_mask.items()))


def exact_cover(full, items, cap, initial=0):
    masks = tuple(row[0] for row in items)
    residues = tuple(row[1] for row in items)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    require(all(by_bit[bit] for bit in range(width) if not initial & (1 << bit)),
            (width, "unreachable"))

    @lru_cache(maxsize=None)
    def solve(state, slots):
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted(
            ((mask & missing).bit_count() for mask in masks if mask | state != state),
            reverse=True,
        )
        if not gains or sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing & (1 << bit))
        pivot = min(
            missing_bits,
            key=lambda bit: sum(1 for index in by_bit[bit]
                                if masks[index] | state != state),
        )
        candidates = sorted(
            (index for index in by_bit[pivot] if masks[index] | state != state),
            key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
        )
        for index in candidates:
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (residues[index],) + suffix
        return None

    return solve(initial, cap)


def minimum_cover(mask_function, modulus, coefficient_modulus, cap,
                  initial=0, excluded=()):
    full = (1 << modulus) - 1
    bank = unique_bank(mask_function, modulus, coefficient_modulus, set(excluded))
    for slots in range(cap + 1):
        witness = exact_cover(full, bank, slots, initial)
        if witness is not None:
            return slots, tuple(sorted(witness))
    return None, ()


def identity_audit():
    chart_rows = 0
    sheet_cells = 0
    order_rows = 0
    unit_rows = 0
    for numerator, denominator in THRESHOLDS:
        require(0 < 2 * numerator <= denominator, (numerator, denominator))
        for modulus in range(2, MAX_IDENTITY_R + 1):
            ambient = 2 * modulus
            even_sheets = sum(1 << (2 * j) for j in range(modulus))
            require(
                fixed_mask(ambient, modulus, numerator, denominator) == even_sheets,
                (modulus, numerator, denominator, "order_two"),
            )
            require(
                half_mask(modulus, modulus, numerator, denominator) == 0,
                (modulus, numerator, denominator, "empty_half"),
            )
            for residue in range(ambient):
                upstairs = fixed_mask(ambient, residue, numerator, denominator)
                downstairs_zero = fixed_mask(modulus, residue, numerator, denominator)
                downstairs_half = half_mask(modulus, residue, numerator, denominator)
                require(
                    extract_parity(upstairs, modulus, 0) == downstairs_zero,
                    (modulus, residue, numerator, denominator, "even_chart"),
                )
                require(
                    extract_parity(upstairs, modulus, 1) == downstairs_half,
                    (modulus, residue, numerator, denominator, "odd_chart"),
                )
                divisor = gcd(modulus, residue)
                m = modulus // divisor
                normalized = residue // divisor
                predicted = 2 * m // gcd(2, normalized)
                require(
                    additive_order(ambient, residue) == predicted,
                    (modulus, residue, predicted),
                )
                chart_rows += 1
                sheet_cells += 2 * modulus
                order_rows += 1

            if modulus <= 32:
                for unit in range(1, ambient):
                    if gcd(unit, ambient) != 1:
                        continue
                    affine_shift = (unit - 1) // 2
                    for residue in range(ambient):
                        zero = fixed_mask(ambient, residue, numerator, denominator)
                        zero_unit = fixed_mask(ambient, unit * residue, numerator, denominator)
                        require(
                            zero_unit
                            == pullback_mask(zero, ambient, lambda x, u=unit: u * x % ambient),
                            (modulus, unit, residue, "fixed_unit"),
                        )
                        half = half_mask(modulus, residue, numerator, denominator)
                        half_unit = half_mask(modulus, unit * residue,
                                              numerator, denominator)
                        require(
                            half_unit
                            == pullback_mask(
                                half,
                                modulus,
                                lambda x, u=unit, b=affine_shift: (u * x + b) % modulus,
                            ),
                            (modulus, unit, residue, "half_affine_unit"),
                        )
                        unit_rows += 1
    return chart_rows, sheet_cells, order_rows, unit_rows


def rank_and_control_audit():
    rows = []
    for modulus, named in sorted(RANK_CONTROLS.items()):
        full_half = (1 << modulus) - 1
        require(family_union(half_mask, modulus, named) == full_half,
                (modulus, named, "named_half"))

        half_rank, half_witness = minimum_cover(half_mask, modulus, 2 * modulus, 7,
                                                excluded=(modulus,))
        require(half_rank == EXPECTED_RANKS[modulus],
                (modulus, half_rank, half_witness))

        ambient = 2 * modulus
        order_two = fixed_mask(ambient, modulus)
        fixed_rank_extra, fixed_witness = minimum_cover(
            fixed_mask,
            ambient,
            ambient,
            7,
            initial=order_two,
            excluded=(modulus,),
        )
        require(fixed_rank_extra == half_rank,
                (modulus, fixed_rank_extra, half_rank, fixed_witness))
        fixed_named = (modulus,) + named
        require(family_union(fixed_mask, ambient, fixed_named) == (1 << ambient) - 1,
                (modulus, fixed_named, "named_fixed"))
        rows.append((modulus, half_rank, 1 + fixed_rank_extra,
                     half_witness, fixed_witness))
    return tuple(rows)


def primitive_and_multiplicity_audit():
    all_even = (2, 4, 6, 8, 10, 14)
    modulus = 15
    require(family_union(half_mask, modulus, all_even) == (1 << modulus) - 1,
            "R15 all-even half cover")
    require(gcd(modulus, *all_even) == 1, "joint quotient period")
    require(gcd(2 * modulus, *all_even) == 2, "half augmented parity defect")
    upstairs = (modulus,) + all_even
    require(family_union(fixed_mask, 2 * modulus, upstairs) == (1 << (2 * modulus)) - 1,
            "Q30 fixed transplant")
    require(gcd(2 * modulus, *upstairs) == 1, "upstairs primitive")

    multiplicity_rows = 0
    for modulus, residues in sorted(RANK_CONTROLS.items()):
        for sheet in range(modulus):
            even_up = sum(
                1 for residue in (modulus,) + residues
                if fixed_mask(2 * modulus, residue) & (1 << (2 * sheet))
            )
            odd_up = sum(
                1 for residue in (modulus,) + residues
                if fixed_mask(2 * modulus, residue) & (1 << (2 * sheet + 1))
            )
            zero_down = sum(
                1 for residue in residues
                if fixed_mask(modulus, residue) & (1 << sheet)
            )
            half_down = sum(
                1 for residue in residues
                if half_mask(modulus, residue) & (1 << sheet)
            )
            require(even_up == 1 + zero_down, (modulus, sheet, "even_mult"))
            require(odd_up == half_down, (modulus, sheet, "odd_mult"))
            multiplicity_rows += 1
    return all_even, upstairs, multiplicity_rows


def main():
    chart = identity_audit()
    ranks = rank_and_control_audit()
    all_even, upstairs, multiplicity_rows = primitive_and_multiplicity_audit()

    semantic = sha256(repr((chart, ranks, all_even, upstairs,
                            multiplicity_rows)).encode()).hexdigest()
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node")
    require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                    for node in ast.walk(tree)), "float literal")

    print("THM-3432 order-two fixed/half parity transplant -- exact controls")
    print("status=VERIFIED_EXACT_CONTROLS_FOR_ELEMENTARY_ALL_MODULUS_PROOF;no_LRC14_cut")
    print("thresholds=" + repr(THRESHOLDS))
    print("identity_universe=(max_R,chart_rows,sheet_cells,order_rows,unit_rows)="
          + repr((MAX_IDENTITY_R,) + chart))
    print("direct_minimum_rows=(R,half_rank,fixed_with_order2_rank,half_witness,fixed_extra_witness)="
          + repr(ranks))
    print("primitive_parity_hostile=(R,all_even_half,upstairs_fixed,gcd_half,gcd_upstairs)="
          + repr((15, all_even, upstairs, gcd(30, *all_even), gcd(30, *upstairs))))
    print("multiplicity_rows=" + str(multiplicity_rows))
    print("semantic_sha256=" + semantic)
    print("script_sha256_lf=" + lf_sha256(source))
    print("scope=literal_fixed_zero_with_named_order2_owner_iff_literal_half_twist;"
          "joint_period_preserved;augmented_parity_sidecar_required;"
          "not_all_fixed_covers;not_mixed_rank7;not_LRC14")


if __name__ == "__main__":
    main()
