#!/usr/bin/env python3
"""Fraction-exact referee for the infinite THM-4119 family."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd


BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
AP8 = tuple(range(1, 9))
D0 = (3, 4, 5, 6, 8, 10, 12)
PHASE = Fraction(9, 19)
EXPECTED_SEMANTIC = "b1d6101e71ed922a4f1a873dcd11d0f7c41be9bd1f2a3c1a6c53f8f07c61ac35"


def gcd_all(values):
    return reduce(gcd, values)


def norm(x):
    residue = x % 1
    return min(residue, 1 - residue)


def clearance(values, theta=PHASE):
    return min(norm(v * theta) for v in values)


def row(t):
    assert t > max(BODY)
    return tuple(sorted(BODY + (t, 3 * t)))


def direct_missing(values, base):
    return tuple(sorted(set(base) - set(values)))


def ap7_failure(values):
    odd = tuple(v for v in values if v % 2)
    even = tuple(v for v in values if not (v % 2))
    if len(odd) != 2:
        return {"odd_count": len(odd), "reason": "not_two_odd_tails"}
    block_gcd = gcd_all(even)
    core = tuple(v // block_gcd for v in even)
    missing = tuple(sorted(set(range(1, 8)) - set(core)))
    return {
        "odd_count": 2,
        "block_gcd": block_gcd,
        "core": core,
        "missing": missing,
        "reason": "missing_ap7_core" if missing else "passed_core",
    }


def representative(residue):
    t = residue
    while t <= max(BODY):
        t += 19
    return t


def main():
    assert gcd_all(BODY) == 1
    body_clearance = clearance(BODY)
    assert body_clearance == Fraction(2, 19)

    clearance_by_residue = []
    structural = {}
    for residue in range(19):
        t = representative(residue)
        values = row(t)
        assert len(values) == 13 and gcd_all(values) == 1
        value = clearance(values)
        clearance_by_residue.append(value)
        assert direct_missing(values, AP8) == (2, 3, 5, 7)
        assert direct_missing(values, D0) == (3, 5)
        failure = ap7_failure(values)
        if t % 2:
            assert failure == {"odd_count": 4, "reason": "not_two_odd_tails"}
        else:
            assert failure["block_gcd"] == 2
            assert failure["missing"] == (1,)
        structural[residue] = (t, failure["odd_count"], failure["reason"])

        # Finite dilation controls mirror the symbolic proof in the theorem.
        for k in (1, 2, 3, 7, 19, 38):
            dilated = tuple(k * v for v in values)
            assert tuple(v // gcd_all(dilated) for v in dilated) == values
            assert direct_missing(dilated, AP8)
            assert direct_missing(dilated, D0)
            if k % 2 == 0 or t % 2:
                assert len(tuple(v for v in dilated if v % 2)) != 2
            else:
                dilated_failure = ap7_failure(dilated)
                assert dilated_failure["missing"] == (1,)

    numerators = tuple(int(value * 19) for value in clearance_by_residue)
    safe = tuple(r for r, numerator in enumerate(numerators) if numerator >= 2)
    excluded = tuple(r for r, numerator in enumerate(numerators) if numerator < 2)
    assert safe == (1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 18)
    assert excluded == (0, 2, 7, 12, 17)
    assert numerators[0] == 0
    assert all(numerators[r] == 1 for r in (2, 7, 12, 17))
    assert all(numerators[r] == 2 for r in safe)
    assert (2**45) % 19 == 18
    assert clearance(row(2**45)) == Fraction(2, 19)

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
    print("family=S_t=BODY_union_{t,3t};t>22")
    print(f"body={BODY};phase={PHASE};body_clearance={body_clearance}")
    print(f"clearance_numerators_over_19_by_residue={numerators}")
    print(f"safe_residues_mod19={safe};excluded={excluded};density=14/19")
    print("supplier_exclusion=AP8_missing_(2,3,5,7);D0_missing_(3,5);AP7_parity_or_missing_unit")
    print(f"canonical_t=2^45;residue={(2**45) % 19};clearance={clearance(row(2**45))}")
    print("dilation_controls=(1,2,3,7,19,38);normalization_and_exclusion=PASS")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
