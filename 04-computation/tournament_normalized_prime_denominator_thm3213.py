#!/usr/bin/env python3
"""Exact prime-denominator controls for THM-3213."""

from functools import lru_cache
from math import comb, factorial


def ordered_stirling_profile(n):
    row = [1]
    for m in range(1, n + 1):
        nxt = [0] * (m + 1)
        for k in range(1, m + 1):
            left = row[k] if k < len(row) else 0
            right = row[k - 1] if k - 1 < len(row) else 0
            nxt[k] = k * (left + right)
        row = nxt
    return row


@lru_cache(None)
def numerator_power(d):
    base = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (0, 1, 1): 1,
        (1, 1, 1): 3,
    }
    poly = {(0, 0, 0): 1}
    for _ in range(d):
        out = {}
        for a, va in poly.items():
            for b, vb in base.items():
                e = tuple(a[i] + b[i] for i in range(3))
                out[e] = out.get(e, 0) + va * vb
        poly = out
    return poly


def kernel(d, c1, c2, c3):
    answer = 0
    for q in range(min(c1, c2, c3) + 1):
        coeff = numerator_power(d).get((c1 - q, c2 - q, c3 - q), 0)
        answer += coeff * comb(d + q - 1, d - 1)
    return answer


def output_coordinate(n, d):
    profile = ordered_stirling_profile(n)
    answer = 0
    for c1 in range(1, n + 1):
        for c2 in range(1, n + 1):
            for c3 in range(1, n + 1):
                w = kernel(d, c1, c2, c3)
                if w:
                    answer += w * profile[c1] * profile[c2] * profile[c3]
    return answer, profile


primes = (5, 7, 11, 13, 17, 19)
for p in primes:
    residues = []
    for d in (1, 2, 3):
        value, profile = output_coordinate(p, d)
        expected = (3, 6, 6)[d - 1]
        if value % p != expected % p:
            raise RuntimeError((p, d, value % p, expected))
        residues.append(value % p)
    if profile[1] % p != 1 or any(x % p for x in profile[2:]):
        raise RuntimeError((p, profile))
    print(f"p={p};F_output_d1_d2_d3_mod_p={tuple(residues)}")

print("w_d(1,1,1)=" + str(tuple(kernel(d, 1, 1, 1) for d in range(1, 8))))
print("consequence=d=1,2,3 factorial-normalized reduced denominators contain p^3")


print('scope=sample_primes_supporting_the_all_prime_proof_for_d_1_2_3')
print('FAILED_CHECKS=NONE')
