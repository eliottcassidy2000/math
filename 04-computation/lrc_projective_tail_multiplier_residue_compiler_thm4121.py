#!/usr/bin/env python3
"""Fraction-exact projective residue and supplier referee for THM-4121."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
AP8 = tuple(range(1, 9))
D0 = (3, 4, 5, 6, 8, 10, 12)
MODULUS = 19
PHASE = Fraction(9, MODULUS)
EXPECTED_SEMANTIC = "62712d11643e5e93ba830a327226530739cd5f2bbace569c472751516166a089"


def gcd_all(values):
    return reduce(gcd, values)


def norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(values):
    return min(norm(PHASE * value) for value in values)


def row(c, t):
    assert c > 1 and t > max(BODY)
    values = tuple(sorted(BODY + (t, c * t)))
    assert len(values) == 13
    return values


def multiplier_representative(residue):
    if residue == 0:
        return MODULUS
    if residue == 1:
        return MODULUS + 1
    return residue


def parameter_representative(residue):
    return next(
        value for value in range(max(BODY) + 1, max(BODY) + 1 + MODULUS)
        if value % MODULUS == residue
    )


def predicted_bad_residues(c_residue):
    c_residue %= MODULUS
    if c_residue == 0:
        return tuple(range(MODULUS))
    first_inverse = pow(9, -1, MODULUS)
    second_inverse = pow((9 * c_residue) % MODULUS, -1, MODULUS)
    return tuple(sorted({
        0,
        first_inverse,
        (-first_inverse) % MODULUS,
        second_inverse,
        (-second_inverse) % MODULUS,
    }))


def supplier_case(values, c, t, k):
    dilated = tuple(k * value for value in values)
    assert gcd_all(dilated) == k
    assert tuple(value // gcd_all(dilated) for value in dilated) == values
    assert not set(AP8) <= set(dilated)
    assert not set(D0) <= set(dilated)

    odds = tuple(value for value in dilated if value % 2)
    if k % 2 == 0:
        assert not odds
        return "k_even:no_odd_tails"
    if t % 2 and c % 2:
        assert len(odds) == 4
        return "k_odd_t_odd_c_odd:four_odds"
    if t % 2:
        assert c % 2 == 0 and len(odds) == 3
        return "k_odd_t_odd_c_even:three_odds"

    assert odds == (k, 15 * k)
    even_block = tuple(value for value in dilated if not value % 2)
    assert gcd_all(even_block) == 2 * k
    forced_core = tuple(value // (2 * k) for value in even_block)
    assert 1 not in forced_core
    return "k_odd_t_even:forced_core_missing_unit"


def main():
    assert gcd_all(BODY) == 1
    assert clearance(BODY) == Fraction(2, MODULUS)
    assert pow(9, -1, MODULUS) == 17

    bad_by_multiplier = []
    safe_counts = []
    supplier_cases_seen = set()
    for c_residue in range(MODULUS):
        c = multiplier_representative(c_residue)
        actual_bad = []
        for t_residue in range(MODULUS):
            t = parameter_representative(t_residue)
            values = row(c, t)
            assert gcd_all(values) == 1
            value = clearance(values)
            if value < Fraction(2, MODULUS):
                actual_bad.append(t_residue)
            else:
                assert value == Fraction(2, MODULUS)

            numerator = int(value * MODULUS)
            if c_residue == 0 or t_residue == 0:
                assert numerator == 0
            elif value < Fraction(2, MODULUS):
                assert numerator == 1
            else:
                assert numerator == 2
            assert tuple(sorted(set(AP8) - set(values))) == (2, 3, 5, 7)
            assert tuple(sorted(set(D0) - set(values))) == (3, 5)
            for k in (1, 2, 3, 7, 19, 38):
                supplier_cases_seen.add(supplier_case(values, c, t, k))

        actual_bad = tuple(actual_bad)
        assert actual_bad == predicted_bad_residues(c_residue)
        bad_by_multiplier.append(actual_bad)
        safe_counts.append(MODULUS - len(actual_bad))

    bad_by_multiplier = tuple(bad_by_multiplier)
    safe_counts = tuple(safe_counts)
    assert safe_counts[0] == 0
    assert tuple(c for c in range(1, MODULUS) if safe_counts[c] == 16) == (1, 18)
    assert all(safe_counts[c] == 14 for c in range(2, 18))
    unavoidable = {0, pow(9, -1, MODULUS), -pow(9, -1, MODULUS) % MODULUS}
    assert all(unavoidable <= set(bad_by_multiplier[c]) for c in range(1, MODULUS))
    assert max(safe_counts) == 16
    assert bad_by_multiplier[3] == (0, 2, 7, 12, 17)
    assert supplier_cases_seen == {
        "k_even:no_odd_tails",
        "k_odd_t_even:forced_core_missing_unit",
        "k_odd_t_odd_c_even:three_odds",
        "k_odd_t_odd_c_odd:four_odds",
    }

    ledger = {
        "body": BODY,
        "body_clearance": (2, 19),
        "phase": (9, 19),
        "bad_residues_by_multiplier": bad_by_multiplier,
        "safe_counts_by_multiplier": safe_counts,
        "optimal_multiplier_residues": (1, 18),
        "sharp_max_safe_classes": 16,
        "clearance_rule_over_19": (
            "c_zero_or_t_zero:0",
            "nonzero_bad:1",
            "safe:2",
        ),
        "thm4119_multiplier": 3,
        "thm4119_bad_residues": bad_by_multiplier[3],
        "direct_missing_ap8": (2, 3, 5, 7),
        "direct_missing_d0": (3, 5),
        "supplier_cases": tuple(sorted(supplier_cases_seen)),
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    print("status=PASS")
    print("family=S_(c,t)=BODY_union_{t,ct};c>1;t>22;phase=9/19")
    print(f"body={BODY};body_clearance=2/19")
    print(f"bad_residues_by_c_mod19={bad_by_multiplier}")
    print(f"safe_counts_by_c_mod19={safe_counts}")
    print("classification=c=0:0/19;c=+/-1:16/19;all_other_nonzero_c:14/19")
    print("clearance_rule=c_zero_or_t_zero:0;nonzero_bad:1/19;safe:2/19")
    print("sharpness=every_nonzero_multiplier_loses_t=0,+/-17;maximum=16/19")
    print("supplier_exclusion=AP8_and_D0_direct_core_failure;AP7_parity_or_forced_core_missing_unit")
    print("dilation_controls=(1,2,3,7,19,38);normalization_and_exclusion=PASS")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
