#!/usr/bin/env python3
"""Clean-room hostile audit for THM-4117.

This implementation imports no primary code.  It enumerates every possible
choice of AP7 odd tails and every admissible scale divisor of the remaining
block, and it evaluates all phase claims by integer residues.
"""

from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import gcd


ROW = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22, 2**45, 3 * 2**45)
AP8 = tuple(range(1, 9))
D0 = (3, 4, 5, 6, 8, 10, 12)
EXPECTED_SEMANTIC = "73563bafba2bc5e8dd870aee18e32047a053eed8a9bce84ad33e7aaecceb0f88"


def gcd_all(values):
    return reduce(gcd, values)


def positive_divisors(n):
    small = []
    large = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    return tuple(small + large[::-1])


def primitive(values):
    g = gcd_all(values)
    return tuple(v // g for v in values)


def direct_supplier(values, base, outlier_count, first):
    values = set(values)
    base = set(base)
    missing = tuple(sorted(base - values))
    if missing:
        return False, missing, ()
    outliers = tuple(sorted(values - base))
    valid = (
        len(outliers) == outlier_count
        and outliers[0] >= first
        and all(y >= 2 * x for x, y in zip(outliers, outliers[1:]))
    )
    return valid, (), outliers


def ap7_witnesses(values):
    """Enumerate all tail choices and all possible q; do not infer q by gcd."""
    values = tuple(values)
    found = []
    attempts = []
    for tail_indices in combinations(range(len(values)), 2):
        tails = tuple(values[i] for i in tail_indices)
        if any(x % 2 == 0 for x in tails):
            continue
        block = tuple(values[i] for i in range(len(values))
                      if i not in tail_indices)
        if any(x % 2 for x in block):
            continue
        half_gcd = gcd_all(block) // 2
        for q in positive_divisors(half_gcd):
            scale = 2 * q
            if any(x % scale for x in block):
                continue
            core = tuple(x // scale for x in block)
            missing = tuple(sorted(set(range(1, 8)) - set(core)))
            outliers = tuple(sorted(set(core) - set(range(1, 8))))
            attempts.append((tails, q, core, missing, outliers))
            if missing or len(outliers) != 4:
                continue
            a, b, c, d = outliers
            adaptive = 84 <= a < b < c < d and (a % 2 == 0 or b % 2 == 0)
            parity_free = 47 <= a and b >= 2 * a and b < c < d
            if adaptive or parity_free:
                found.append((tails, q, core, outliers,
                              adaptive, parity_free))
    return tuple(found), tuple(attempts)


def norm_numerator(v, numerator, denominator):
    residue = (v * numerator) % denominator
    return min(residue, denominator - residue)


def clearance_numerator(values, numerator, denominator):
    return min(norm_numerator(v, numerator, denominator) for v in values)


def recognize_controls():
    ap8 = AP8 + (94, 188, 376, 752, 1504)
    d0 = D0 + (16, 32, 64, 128, 256, 512)
    ap7_core = tuple(range(1, 8)) + (85, 86, 91, 101)
    ap7 = tuple(sorted(tuple(2 * x for x in ap7_core) + (103, 105)))
    assert direct_supplier(ap8, AP8, 5, 94)[0]
    assert direct_supplier(d0, D0, 6, 16)[0]
    assert len(ap7_witnesses(ap7)[0]) == 1
    return ap7, ap8, d0


def main():
    assert len(ROW) == 13 and len(set(ROW)) == 13
    assert gcd_all(ROW) == 1
    assert sum(ROW) == 140737488355454 < 91**12
    assert primitive(tuple(101 * x for x in ROW)) == ROW

    ap8_ok, ap8_missing, _ = direct_supplier(ROW, AP8, 5, 94)
    d0_ok, d0_missing, _ = direct_supplier(ROW, D0, 6, 16)
    ap7_found, ap7_attempts = ap7_witnesses(ROW)
    assert not ap8_ok and ap8_missing == (2, 3, 5, 7)
    assert not d0_ok and d0_missing == (3, 5)
    assert not ap7_found
    assert len(ap7_attempts) == 1
    tails, q, core, missing, outliers = ap7_attempts[0]
    assert tails == (1, 15) and q == 1
    assert core == (2, 3, 4, 5, 6, 7, 8, 9, 11, 2**44, 3 * 2**44)
    assert missing == (1,)
    assert outliers == (8, 9, 11, 2**44, 3 * 2**44)

    physical_num = clearance_numerator(ROW, 9, 19)
    theta_num_38 = clearance_numerator(ROW, 18, 38)
    shifted_num_38 = clearance_numerator(ROW, 37, 38)
    anti_num = min(theta_num_38, shifted_num_38)
    assert physical_num == 2
    assert anti_num == 1
    firewall_nums = tuple(
        clearance_numerator(ROW, a, 112) for a in (4, 60, 5, 61)
    )
    assert firewall_nums == (4, 4, 2, 2)

    even = tuple(x for x in ROW if x % 2 == 0)
    divided = tuple(x // 2 for x in even)
    residues = tuple(x % 56 for x in divided)
    assert gcd_all(even) == 2
    assert residues == (2, 3, 4, 5, 6, 7, 8, 9, 11, 32, 40)

    controls = recognize_controls()
    ledger = {
        "row": ROW,
        "row_gcd": gcd_all(ROW),
        "row_sum": sum(ROW),
        "box_bound": 91**12,
        "odd_tails": tails,
        "ap7_q_candidates": tuple(x[1] for x in ap7_attempts),
        "forced_core": core,
        "forced_missing_ap7": missing,
        "ap8_missing": ap8_missing,
        "d0_missing": d0_missing,
        "physical_clearance": (physical_num, 19),
        "antipodal_clearance": (anti_num, 38),
        "firewall_numerators_over_112": firewall_nums,
        "divided_mod_56": residues,
        "positive_control_sizes": tuple(len(x) for x in controls),
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    print("status=PASS")
    print("implementation=clean-room exhaustive tail/scale recognizer plus integer-residue phase audit")
    print(f"row_gcd={ledger['row_gcd']};row_sum={ledger['row_sum']};box_bound={ledger['box_bound']}")
    print(f"AP7_witnesses={len(ap7_found)};tail_scale_attempts={len(ap7_attempts)};forced_tails={tails};forced_q={q};missing={missing}")
    print(f"AP8=False;missing={ap8_missing};D0=False;missing={d0_missing}")
    print(f"physical_clearance={physical_num}/19;antipodal_clearance={anti_num}/38;firewall_numerators_over_112={firewall_nums}")
    print(f"positive_controls=PASS;sizes={ledger['positive_control_sizes']}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
