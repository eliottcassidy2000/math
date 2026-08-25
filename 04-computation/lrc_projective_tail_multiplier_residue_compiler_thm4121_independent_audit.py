#!/usr/bin/env python3
"""Clean-room modular and exhaustive supplier audit for THM-4121."""

from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
AP8 = tuple(range(1, 9))
D0 = (3, 4, 5, 6, 8, 10, 12)
MODULUS = 19
EXPECTED_SEMANTIC = "62712d11643e5e93ba830a327226530739cd5f2bbace569c472751516166a089"


def gcd_tuple(values):
    return reduce(gcd, values)


def residue_distance(value):
    value %= MODULUS
    return min(value, MODULUS - value)


def family(c, t):
    assert c > 1 and t > 22
    values = tuple(sorted(BODY + (t, c * t)))
    assert len(values) == 13
    return values


def multiplier_representative(residue):
    return {0: 19, 1: 20}.get(residue, residue)


def parameter_representative(residue):
    return next(value for value in range(23, 42) if value % MODULUS == residue)


def divisors(value):
    return tuple(candidate for candidate in range(1, value + 1) if value % candidate == 0)


def ap7_attempts(values):
    attempts = []
    for tail_indices in combinations(range(len(values)), 2):
        tails = tuple(values[index] for index in tail_indices)
        if any(not value % 2 for value in tails):
            continue
        block = tuple(
            value for index, value in enumerate(values) if index not in tail_indices
        )
        if any(value % 2 for value in block):
            continue
        for q in divisors(gcd_tuple(block) // 2):
            scale = 2 * q
            if any(value % scale for value in block):
                continue
            core = tuple(value // scale for value in block)
            missing = tuple(sorted(set(range(1, 8)) - set(core)))
            attempts.append((tails, q, core, missing))
    return tuple(attempts)


def positive_supplier_controls():
    ap8 = tuple(range(1, 9)) + (94, 188, 376, 752, 1504)
    d0 = (3, 4, 5, 6, 8, 10, 12, 16, 32, 64, 128, 256, 512)
    ap7_core = tuple(range(1, 8)) + (85, 86, 91, 101)
    ap7 = tuple(sorted(tuple(2 * value for value in ap7_core) + (103, 105)))
    assert set(AP8) <= set(ap8)
    assert set(D0) <= set(d0)
    successful = [attempt for attempt in ap7_attempts(ap7) if not attempt[3]]
    assert len(successful) == 1
    return tuple(len(values) for values in (ap7, ap8, d0))


def main():
    body_numerators = tuple(residue_distance(9 * value) for value in BODY)
    assert min(body_numerators) == 2

    bad_by_multiplier = []
    safe_counts = []
    observed_odd_counts = set()
    for c_residue in range(MODULUS):
        c = multiplier_representative(c_residue)
        bad = []
        for t_residue in range(MODULUS):
            tail_numerators = (
                residue_distance(9 * t_residue),
                residue_distance(9 * c_residue * t_residue),
            )
            exact_numerator = min((2,) + tail_numerators)
            if exact_numerator < 2:
                bad.append(t_residue)
            if c_residue == 0 or t_residue == 0:
                assert exact_numerator == 0
            elif exact_numerator < 2:
                assert exact_numerator == 1
            else:
                assert exact_numerator == 2

            t = parameter_representative(t_residue)
            values = family(c, t)
            assert gcd_tuple(values) == 1
            assert tuple(sorted(set(AP8) - set(values))) == (2, 3, 5, 7)
            assert tuple(sorted(set(D0) - set(values))) == (3, 5)
            for k in (1, 2, 3, 7, 19, 38):
                dilated = tuple(k * value for value in values)
                assert tuple(value // gcd_tuple(dilated) for value in dilated) == values
                assert not set(AP8) <= set(dilated)
                assert not set(D0) <= set(dilated)
                attempts = ap7_attempts(dilated)
                assert not any(not attempt[3] for attempt in attempts)
                odd_count = sum(value % 2 for value in dilated)
                observed_odd_counts.add(odd_count)
                if k % 2 == 0:
                    assert odd_count == 0 and not attempts
                elif t % 2:
                    assert odd_count == 3 + (c % 2)
                    assert not attempts
                else:
                    assert odd_count == 2 and attempts
                    exact_scale = [attempt for attempt in attempts if attempt[1] == k]
                    assert len(exact_scale) == 1
                    assert exact_scale[0][0] == (k, 15 * k)
                    assert exact_scale[0][3] == (1,)

        bad = tuple(bad)
        bad_by_multiplier.append(bad)
        safe_counts.append(MODULUS - len(bad))

    bad_by_multiplier = tuple(bad_by_multiplier)
    safe_counts = tuple(safe_counts)
    assert bad_by_multiplier[0] == tuple(range(MODULUS))
    assert safe_counts[0] == 0
    assert tuple(c for c in range(1, MODULUS) if safe_counts[c] == 16) == (1, 18)
    assert all(safe_counts[c] == 14 for c in range(2, 18))
    assert all({0, 2, 17} <= set(bad_by_multiplier[c]) for c in range(1, MODULUS))
    assert max(safe_counts) == 16
    assert bad_by_multiplier[3] == (0, 2, 7, 12, 17)
    assert observed_odd_counts == {0, 2, 3, 4}

    supplier_cases = (
        "k_even:no_odd_tails",
        "k_odd_t_even:forced_core_missing_unit",
        "k_odd_t_odd_c_even:three_odds",
        "k_odd_t_odd_c_odd:four_odds",
    )
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
        "supplier_cases": supplier_cases,
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    control_sizes = positive_supplier_controls()
    print("status=PASS")
    print("implementation=integer-residue table plus exhaustive odd-tail/scale recognizer")
    print(f"body={BODY};body_clearance=2/19;phase=9/19")
    print(f"bad_residues_by_c_mod19={bad_by_multiplier}")
    print(f"safe_counts_by_c_mod19={safe_counts}")
    print("classification=c=0:0/19;c=+/-1:16/19;all_other_nonzero_c:14/19")
    print("clearance_rule=c_zero_or_t_zero:0;nonzero_bad:1/19;safe:2/19")
    print("sharpness=unavoidable_bad_residues_(0,2,17);maximum=16/19")
    print("supplier_exclusion=exhaustive_AP7_tail_and_scale_search;direct_AP8_and_D0_failure")
    print(f"positive_supplier_controls=PASS;sizes={control_sizes}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
