#!/usr/bin/env python3
"""Clean-room modular and exhaustive-shape audit for THM-4119."""

from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
EXPECTED_SEMANTIC = "b1d6101e71ed922a4f1a873dcd11d0f7c41be9bd1f2a3c1a6c53f8f07c61ac35"


def gcd_tuple(values):
    return reduce(gcd, values)


def residue_distance(multiplier, residue, modulus=19):
    value = (multiplier * residue) % modulus
    return min(value, modulus - value)


def family(t):
    assert t > 22
    return tuple(sorted(BODY + (t, 3 * t)))


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def ap7_attempts(values):
    values = tuple(values)
    attempts = []
    for indices in combinations(range(len(values)), 2):
        tails = tuple(values[i] for i in indices)
        if any(value % 2 == 0 for value in tails):
            continue
        block = tuple(values[i] for i in range(len(values)) if i not in indices)
        if any(value % 2 for value in block):
            continue
        half_gcd = gcd_tuple(block) // 2
        for q in divisors(half_gcd):
            scale = 2 * q
            if any(value % scale for value in block):
                continue
            core = tuple(value // scale for value in block)
            missing = tuple(sorted(set(range(1, 8)) - set(core)))
            attempts.append((tails, q, core, missing))
    return tuple(attempts)


def smallest_over_22(residue):
    return next(value for value in range(23, 42) if value % 19 == residue)


def recognize_positive_controls():
    ap8 = tuple(range(1, 9)) + (94, 188, 376, 752, 1504)
    d0 = (3, 4, 5, 6, 8, 10, 12, 16, 32, 64, 128, 256, 512)
    ap7_core = tuple(range(1, 8)) + (85, 86, 91, 101)
    ap7 = tuple(sorted(tuple(2 * value for value in ap7_core) + (103, 105)))
    assert set(range(1, 9)) <= set(ap8)
    assert set((3, 4, 5, 6, 8, 10, 12)) <= set(d0)
    attempts = ap7_attempts(ap7)
    assert len(attempts) == 1 and attempts[0][3] == ()
    return tuple(len(values) for values in (ap7, ap8, d0))


def main():
    body_numerators = tuple(residue_distance(9, value % 19) for value in BODY)
    assert min(body_numerators) == 2
    clearance_numerators = []
    structural = {}
    for residue in range(19):
        high = min(residue_distance(9, residue), residue_distance(8, residue))
        clearance_numerators.append(min(2, high))

        t = smallest_over_22(residue)
        values = family(t)
        assert len(values) == 13 and gcd_tuple(values) == 1
        assert tuple(sorted(set(range(1, 9)) - set(values))) == (2, 3, 5, 7)
        assert tuple(sorted(set((3, 4, 5, 6, 8, 10, 12)) - set(values))) == (3, 5)
        attempts = ap7_attempts(values)
        odd_count = sum(value % 2 for value in values)
        if t % 2:
            assert odd_count == 4 and not attempts
            reason = "not_two_odd_tails"
        else:
            assert odd_count == 2 and len(attempts) == 1
            tails, q, core, missing = attempts[0]
            assert tails == (1, 15) and q == 1 and missing == (1,)
            assert core[:9] == (2, 3, 4, 5, 6, 7, 8, 9, 11)
            reason = "missing_ap7_core"
        structural[residue] = (t, odd_count, reason)

        for k in (1, 2, 3, 7, 19, 38):
            dilated = tuple(k * value for value in values)
            assert gcd_tuple(dilated) == k
            assert tuple(value // k for value in dilated) == values
            attempts_k = ap7_attempts(dilated)
            if k % 2 == 0 or t % 2:
                assert not attempts_k
            else:
                assert attempts_k and all(attempt[3] for attempt in attempts_k)
                exact_scale = [attempt for attempt in attempts_k if attempt[1] == k]
                assert len(exact_scale) == 1 and exact_scale[0][3] == (1,)

    numerators = tuple(clearance_numerators)
    safe = tuple(r for r, value in enumerate(numerators) if value >= 2)
    excluded = tuple(r for r, value in enumerate(numerators) if value < 2)
    assert safe == (1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 18)
    assert excluded == (0, 2, 7, 12, 17)
    assert numerators == (0, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2)
    assert (2**45) % 19 == 18
    control_sizes = recognize_positive_controls()

    ledger = {
        "body": BODY,
        "body_clearance": (2, 19),
        "phase": (9, 19),
        "clearance_numerators_by_residue": numerators,
        "safe_residues": safe,
        "excluded_residues": excluded,
        "direct_missing_ap8": (2, 3, 5, 7),
        "direct_missing_d0": (3, 5),
        "canonical_residue": (2**45) % 19,
        "natural_density": (14, 19),
        "structural_controls": structural,
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    print("status=PASS")
    print("implementation=integer-residue phase table plus exhaustive odd-tail/scale recognizer")
    print(f"body={BODY};phase=9/19;body_clearance=2/19")
    print(f"clearance_numerators_over_19_by_residue={numerators}")
    print(f"safe_residues_mod19={safe};excluded={excluded};density=14/19")
    print("supplier_exclusion=AP8_missing_(2,3,5,7);D0_missing_(3,5);AP7_parity_or_missing_unit")
    print(f"canonical_t=2^45;residue={(2**45) % 19};clearance=2/19")
    print(f"positive_supplier_controls=PASS;sizes={control_sizes}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
