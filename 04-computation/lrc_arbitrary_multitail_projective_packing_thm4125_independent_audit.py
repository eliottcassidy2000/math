#!/usr/bin/env python3
"""Clean-room integer phase-table audit for THM-4125."""

from hashlib import sha256
from itertools import product
import json
from math import comb, gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
MODULUS = 19
EXPECTED_SEMANTIC = "2b80075c494e496faa4247cbe22d59cb07af20b31cc7c426e82d7be8e5dc7192"


def require(condition, label):
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def digest(value):
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def distance(value):
    value %= MODULUS
    return min(value, MODULUS - value)


def phase_numerator(residues, parameter):
    return min((2,) + tuple(distance(9 * residue * parameter) for residue in residues))


def aligned_multipliers(count):
    return tuple(1 + MODULUS * index for index in range(count))


def main():
    require(min(distance(9 * value) for value in BODY) == 2,
            "integer body clearance")
    signed_pairs = tuple((residue, MODULUS - residue) for residue in range(1, 10))
    support_histogram = {}
    minimum_support = {size: 10 for size in range(1, 19)}
    checked = 0

    # Each signed pair is absent, represented by its positive element, by its
    # negative element, or by both. This enumeration is independent of the
    # inverse-pair construction in the primary implementation.
    for choices in product(range(4), repeat=9):
        if not any(choices):
            continue
        residues = []
        for (positive, negative), choice in zip(signed_pairs, choices):
            if choice in (1, 3):
                residues.append(positive)
            if choice in (2, 3):
                residues.append(negative)
        support = sum(choice != 0 for choice in choices)
        bad = tuple(
            parameter for parameter in range(MODULUS)
            if phase_numerator(tuple(residues), parameter) < 2
        )
        require(len(bad) == 1 + 2 * support, "brute signed-support bad count")
        require(phase_numerator(tuple(residues), 0) == 0,
                "brute zero parameter clearance")
        require(all(
            phase_numerator(tuple(residues), parameter) == 1
            for parameter in bad if parameter
        ), "brute nonzero bad clearance")
        require(all(
            phase_numerator(tuple(residues), parameter) == 2
            for parameter in range(MODULUS) if parameter not in bad
        ), "brute safe clearance")
        support_histogram[support] = support_histogram.get(support, 0) + 1
        size = len(residues)
        minimum_support[size] = min(minimum_support[size], support)
        checked += 1

    expected_histogram = tuple(
        (support, comb(9, support) * 3**support) for support in range(1, 10)
    )
    require(checked == 4**9 - 1 == (1 << 18) - 1,
            "clean-room residue subset count")
    require(tuple(sorted(support_histogram.items())) == expected_histogram,
            "clean-room support histogram")
    distinct_residue_fixed_optima = tuple(
        (size, minimum_support[size], 18 - 2 * minimum_support[size])
        for size in range(1, 19)
    )
    require(distinct_residue_fixed_optima == tuple(
        (size, (size + 1) // 2, 18 - 2 * ((size + 1) // 2))
        for size in range(1, 19)
    ), "clean-room distinct-residue optimum")
    distinct_residue_natural_optima = tuple(
        (
            size,
            minimum_support[size],
            18 if size >= 7 else 18 - 2 * minimum_support[size],
        )
        for size in range(1, 19)
    )
    require(distinct_residue_natural_optima == tuple(
        (
            size,
            (size + 1) // 2,
            18 if size >= 7 else 18 - 2 * ((size + 1) // 2),
        )
        for size in range(1, 19)
    ), "clean-room distinct-residue natural-threshold optimum")

    aligned_bad = tuple(
        parameter for parameter in range(MODULUS)
        if phase_numerator((1,), parameter) < 2
    )
    aligned_safe = tuple(
        parameter for parameter in range(MODULUS) if parameter not in aligned_bad
    )
    require(aligned_bad == (0, 2, 17) and len(aligned_safe) == 16,
            "clean-room aligned optimum")
    aligned_controls = []
    for count in (1, 2, 3, 5, 6, 7, 8, 9, 18, 19, 40):
        multipliers = aligned_multipliers(count)
        require(len(set(multipliers)) == count, "clean-room actual distinctness")
        natural_safe_count = 0
        for parameter in range(MODULUS):
            residues = tuple((multiplier * parameter) % MODULUS for multiplier in multipliers)
            value = min((2,) + tuple(distance(9 * residue) for residue in residues))
            expected = 0 if parameter == 0 else 1 if parameter in aligned_bad else 2
            require(value == expected, "clean-room aligned full tail table")
            natural_safe_count += value * (12 + count) >= MODULUS
        require(2 * (12 + count) > MODULUS, "clean-room LRC threshold")
        expected_natural_safe_count = 18 if count >= 7 else 16
        require(natural_safe_count == expected_natural_safe_count,
                "clean-room aligned natural-threshold phase transition")
        aligned_controls.append(
            (count, 11 + count, (2, 19), (1, 12 + count), 16, natural_safe_count)
        )

    require(all(
        phase_numerator((1, 0), parameter) == 0 for parameter in range(MODULUS)
    ), "clean-room zero residue hostile")
    fixed_safe_by_support = tuple(
        (support, 18 - 2 * support) for support in range(1, 10)
    )
    signed_distinct_optima = tuple(
        (count, 18 - 2 * count) for count in range(1, 10)
    )
    ledger = {
        "body": BODY,
        "phase": (9, 19),
        "body_clearance": (2, 19),
        "fixed_threshold": (1, 14),
        "signed_pairs": signed_pairs,
        "nonzero_residue_subsets_checked": checked,
        "support_histogram": expected_histogram,
        "fixed_safe_classes_by_signed_support": fixed_safe_by_support,
        "distinct_residue_fixed_optima": distinct_residue_fixed_optima,
        "natural_threshold_rule": (("r_le_6", "18-2s"), ("r_ge_7", 18)),
        "distinct_residue_natural_optima": distinct_residue_natural_optima,
        "signed_distinct_optima": signed_distinct_optima,
        "aligned_exact_bad": aligned_bad,
        "aligned_exact_safe": aligned_safe,
        "aligned_controls": tuple(aligned_controls),
        "zero_residue_hostile_safe_classes": 0,
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "shared frozen semantic digest")

    print("status=PASS")
    print("implementation=four-state signed-pair enumeration plus integer phase tables")
    print(f"nonzero_residue_subsets_checked={checked}")
    print(f"support_histogram={expected_histogram}")
    print(f"fixed_1_over_14_safe_classes_by_signed_support={fixed_safe_by_support}")
    print(f"distinct_residue_fixed_optima={distinct_residue_fixed_optima}")
    print("natural_threshold_rule=r<=6:18-2s;r>=7:18")
    print(f"distinct_residue_natural_optima={distinct_residue_natural_optima}")
    print(f"signed_distinct_optima={signed_distinct_optima}")
    print(f"aligned_exact_bad={aligned_bad};aligned_exact_safe={aligned_safe};fixed_density=16/19")
    print(f"aligned_arbitrary_tail_controls={tuple(aligned_controls)}")
    print("zero_multiplier_residue_hostile=safe_classes_0")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
