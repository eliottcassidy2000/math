#!/usr/bin/env python3
"""Fraction and projective-support referee for THM-4125."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import comb, gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
MODULUS = 19
PHASE = Fraction(9, MODULUS)
EXPECTED_SEMANTIC = "2b80075c494e496faa4247cbe22d59cb07af20b31cc7c426e82d7be8e5dc7192"


def require(condition, label):
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def digest(value):
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def gcd_all(values):
    return reduce(gcd, values)


def norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(values):
    return min(norm(PHASE * value) for value in values)


def distance_mod(value):
    value %= MODULUS
    return min(value, MODULUS - value)


def signed_class(residue):
    residue %= MODULUS
    require(residue != 0, "signed class requires a unit residue")
    return min(residue, MODULUS - residue)


def bad_residues(multiplier_residues):
    residues = tuple(residue % MODULUS for residue in multiplier_residues)
    if any(residue == 0 for residue in residues):
        return tuple(range(MODULUS))
    bad = {0}
    for residue in residues:
        inverse = pow((9 * residue) % MODULUS, -1, MODULUS)
        bad.update((inverse, (-inverse) % MODULUS))
    return tuple(sorted(bad))


def numerator(multiplier_residues, parameter_residue):
    tail = tuple(
        distance_mod(9 * residue * parameter_residue)
        for residue in multiplier_residues
    )
    return min((2,) + tail)


def representative(residue):
    return next(value for value in range(23, 42) if value % MODULUS == residue)


def aligned_multipliers(count):
    return tuple(1 + MODULUS * index for index in range(count))


def row(multipliers, parameter):
    require(parameter > max(BODY), "tail scale above the body")
    require(len(set(multipliers)) == len(multipliers), "distinct multipliers")
    values = tuple(sorted(BODY + tuple(c * parameter for c in multipliers)))
    require(len(values) == len(BODY) + len(multipliers), "distinct row speeds")
    return values


def main():
    require(clearance(BODY) == Fraction(2, MODULUS), "body clearance")
    require(gcd_all(BODY) == 1, "primitive body")

    support_histogram = {}
    minimum_support_by_residue_count = {size: 10 for size in range(1, 19)}
    checked = 0
    for mask in range(1, 1 << 18):
        residues = tuple(
            residue for residue in range(1, MODULUS)
            if mask & (1 << (residue - 1))
        )
        support = len({signed_class(residue) for residue in residues})
        predicted_bad = bad_residues(residues)
        brute_bad = tuple(
            parameter for parameter in range(MODULUS)
            if numerator(residues, parameter) < 2
        )
        require(brute_bad == predicted_bad, "complete projective bad-set formula")
        require(len(predicted_bad) == 1 + 2 * support,
                "signed-support union cardinality")
        for parameter in range(MODULUS):
            value = numerator(residues, parameter)
            if parameter == 0:
                require(value == 0, "zero parameter clearance")
            elif parameter in predicted_bad:
                require(value == 1, "nonzero bad clearance numerator")
            else:
                require(value == 2, "safe clearance numerator")
        support_histogram[support] = support_histogram.get(support, 0) + 1
        size = len(residues)
        minimum_support_by_residue_count[size] = min(
            minimum_support_by_residue_count[size], support
        )
        checked += 1

    require(checked == (1 << 18) - 1, "all nonempty nonzero residue subsets")
    expected_histogram = tuple(
        (support, comb(9, support) * 3**support) for support in range(1, 10)
    )
    require(tuple(sorted(support_histogram.items())) == expected_histogram,
            "signed-support subset histogram")
    distinct_residue_fixed_optima = tuple(
        (
            size,
            minimum_support_by_residue_count[size],
            18 - 2 * minimum_support_by_residue_count[size],
        )
        for size in range(1, 19)
    )
    require(distinct_residue_fixed_optima == tuple(
        (size, (size + 1) // 2, 18 - 2 * ((size + 1) // 2))
        for size in range(1, 19)
    ), "sharp pairwise-distinct-residue optimum")
    distinct_residue_natural_optima = tuple(
        (
            size,
            minimum_support_by_residue_count[size],
            18 if size >= 7 else 18 - 2 * minimum_support_by_residue_count[size],
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
    ), "sharp pairwise-distinct-residue natural-threshold optimum")

    aligned_bad = bad_residues((1,))
    aligned_safe = tuple(
        residue for residue in range(MODULUS) if residue not in aligned_bad
    )
    require(aligned_bad == (0, 2, 17) and len(aligned_safe) == 16,
            "sharp aligned signed class")
    aligned_controls = []
    for count in (1, 2, 3, 5, 6, 7, 8, 9, 18, 19, 40):
        multipliers = aligned_multipliers(count)
        require(len(set(multipliers)) == count, "aligned actual multipliers distinct")
        target = Fraction(1, 12 + count)
        natural_safe_count = 0
        for residue in range(MODULUS):
            parameter = representative(residue)
            values = row(multipliers, parameter)
            require(gcd_all(values) == 1, "aligned family remains primitive")
            actual = clearance(values)
            if residue == 0:
                require(actual == 0, "aligned zero-class clearance")
            elif residue in aligned_bad:
                require(actual == Fraction(1, MODULUS),
                        "aligned nonzero bad clearance")
            else:
                require(actual == Fraction(2, MODULUS),
                        "aligned safe exact clearance")
            natural_safe_count += actual >= target
        require(Fraction(2, MODULUS) > target, "lonely-runner threshold consequence")
        expected_natural_safe_count = 18 if count >= 7 else 16
        require(natural_safe_count == expected_natural_safe_count,
                "aligned natural-threshold phase transition")
        aligned_controls.append(
            (
                count,
                11 + count,
                (2, 19),
                (target.numerator, target.denominator),
                16,
                natural_safe_count,
            )
        )

    require(bad_residues((1, 19)) == tuple(range(MODULUS)),
            "zero multiplier residue hostile")
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
        "signed_pairs": tuple((residue, MODULUS - residue) for residue in range(1, 10)),
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
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=PASS")
    print("family=U_union_{c*t:c_in_C};distinct_positive_C;t>22;phase=9/19")
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
