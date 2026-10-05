#!/usr/bin/env python3
"""Exact endpoint-versus-prefix audit; no floating-point threshold decisions."""
from fractions import Fraction
from math import comb
import json

checks = 0


def check(test, label):
    global checks
    checks += 1
    if not test:
        raise ValueError(label)


def v2(x):
    return (x & -x).bit_length()-1


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(1, total-length+2):
            for tail in compositions(total-first, length-1):
                yield (first,) + tail


def survives(word):
    total = 0
    for j, a in enumerate(word, 1):
        total += a
        if 2**total > 3**j:
            return False
    return True


def main():
    counts = {0: 1}
    rows = []
    rotations = 0
    for k in range(1, 101):
        bound = (3**k).bit_length()-1
        next_counts = {}
        for old, count in counts.items():
            for new in range(old+1, bound+1):
                next_counts[new] = next_counts.get(new, 0)+count
        counts = next_counts
        prefix = sum((Fraction(count, 2**a) for a, count in counts.items()), Fraction())
        endpoint = sum((Fraction(comb(a-1, k-1), 2**a)
                        for a in range(k, bound+1)), Fraction())
        tail = Fraction(sum(comb(bound, j) for j in range(k, bound+1)), 2**bound)
        check(endpoint == tail, "negative-binomial endpoint identity")
        check(endpoint/k <= prefix <= endpoint, "cyclic ballot sandwich")
        if k <= 8:
            words = [w for a in range(k, bound+1) for w in compositions(a, k)]
            brute = sum((Fraction(1, 2**sum(w)) for w in words if survives(w)), Fraction())
            check(brute == prefix, "composition prefix enumeration")
            for w in words:
                check(any(survives(w[i:]+w[:i]) for i in range(k)), "a surviving rotation")
                rotations += 1
            good_end = good_prefix = 0
            for n in range(1, 2**(bound+1), 2):
                x, total, good = n, 0, True
                for j in range(1, k+1):
                    a = v2(3*x+1)
                    total += a
                    good = good and 2**total < 3**j
                    x = (3*x+1)//2**a
                good_end += 2**total < 3**k
                good_prefix += good
            check(Fraction(good_end, 2**bound) == endpoint, "dyadic endpoint census")
            check(Fraction(good_prefix, 2**bound) == prefix, "dyadic prefix census")
        if k <= 8 or k in (20, 50, 100):
            rows.append({"k": k, "endpoint": str(endpoint), "prefix": str(prefix),
                         "ratio_prefix_endpoint": str(prefix/endpoint)})
    check(not survives((2, 1)) and 2**3 < 3**2, "minimal endpoint-prefix hostile")
    print(json.dumps({"status": "PASS; coefficient-prefix events, not a Collatz proof",
                      "rows": rows, "rotation_words": rotations, "checks": checks,
                      "DP_depth": 100, "independent_dyadic_depth": 8,
                      "witness": {"word": [2, 1], "route": [9, 7, 11],
                                  "endpoint_probability_k2": "1/2",
                                  "prefix_probability_k2": "3/8"}}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
