#!/usr/bin/env python3
"""Exact independent finite controls for the nonprimitive repetition lemma."""

from math import gcd
import json


def require(condition, *detail):
    if not condition:
        raise RuntimeError(detail)


def carry(word):
    b = 0
    threes = 1
    for i in range(len(word) - 1, -1, -1):
        if word[i]:
            b += (1 << i) * threes
            threes *= 3
    return b


shapes = occurrences = crossing = positive = negative = 0
for q in range(4, 101):
    for p in range(1, q):
        d = gcd(p, q)
        if d == 1:
            continue
        shapes += 1
        p0, q0 = p // d, q // d
        a, b = 3**p0, 2**q0
        h = sum(b ** (d - 1 - j) * a**j for j in range(d))
        d0, clock = b - a, 2**q - 3**p
        require(h > 1 and gcd(h, 6) == 1, "repetition factor", p, q)
        require(clock == h * d0, "clock factorization", p, q)
        original = [((i + 1) * p) // q - (i * p) // q for i in range(q)]
        require(original == original[:q0] * d, "mechanical repetition", p, q)
        for i in range(q):
            j = (i + 1) % q
            if original[i] == original[j]:
                continue
            occurrences += 1
            crossing += int(i == q - 1)
            positive += int(clock > 0)
            negative += int(clock < 0)
            edited = original.copy()
            edited[i], edited[j] = edited[j], edited[i]
            original_edited_carry = carry(edited)
            if i == q - 1:
                base = original[1:] + original[:1]
                edited_at_cut = edited[1:] + edited[:1]
                edge = q - 2
            else:
                base, edited_at_cut, edge = original, edited, i
            require(base == base[:q0] * d, "rotated repetition", p, q, i)
            base_carry, edited_carry = carry(base), carry(edited_at_cut)
            require(base_carry == h * carry(base[:q0]), "carry factorization", p, q, i)
            exponent_three = sum(base[edge + 2 :])
            expected_difference = (
                (base[edge] - base[edge + 1]) * 2**edge * 3**exponent_three
            )
            require(
                edited_carry - base_carry == expected_difference,
                "exact local difference",
                p,
                q,
                i,
            )
            require(gcd(h, edited_carry) == 1, "rotated H coprimality", p, q, i)
            require(gcd(h, original_edited_carry) == 1, "original-cut H coprimality", p, q, i)
            actual = gcd(abs(clock), original_edited_carry)
            require(
                actual == gcd(abs(d0), original_edited_carry) <= abs(d0) < abs(clock),
                "exact gcd reduction and properness",
                p,
                q,
                i,
            )
print(
    json.dumps(
        {
            "status": "PASS",
            "maximum_q": 100,
            "nonprimitive_shapes": shapes,
            "actual_swap_occurrences": occurrences,
            "cut_crossing_occurrences": crossing,
            "positive_clock_occurrences": positive,
            "negative_clock_occurrences": negative,
        },
        indent=2,
        sort_keys=True,
    )
)
