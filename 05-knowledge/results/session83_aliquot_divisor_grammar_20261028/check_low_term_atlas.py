#!/usr/bin/env python3
"""Independent proper-divisor-sieve audit of the complete layers zero to nine."""

import argparse
import json
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def pairs(n):
    ds = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ds.update((d, n // d))
    ds.discard(n)
    return sorted(ds)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100000)
    limit = parser.parse_args().limit
    require(8128 <= limit <= 2000000, "limit must be in [8128,2000000]")
    sums, counts = [0] * (limit + 1), [0] * (limit + 1)
    for d in range(1, limit // 2 + 1):
        for n in range(2 * d, limit + 1, d):
            sums[n] += d
            counts[n] += 1
    # Primes are obtained independently from a proper-divisor condition:
    # every prime has exactly the single proper divisor1, and conversely.
    primes = [n for n in range(2, limit + 1) if counts[n] == 1]
    family = [6 * p for p in primes if 5 <= p <= limit // 6]
    expected_abundant = {k: [] for k in range(10)}
    expected_abundant[5] = [12, 18, 20]
    expected_abundant[7] = sorted([24, 40, 54, 56, 70, 88, 104] + family)
    expected_abundant[8] = [36, 100, 196]
    expected_abundant[9] = [48, 80, 112, 162, 176, 208, 272, 304, 368, 464]
    expected_perfect = {k: [] for k in range(10)}
    expected_perfect.update({3: [6], 5: [28], 9: [496]})
    rows = {}
    for k in range(10):
        sources = [n for n in range(1, limit + 1) if counts[n] == k]
        abundant = [n for n in sources if sums[n] > n]
        perfect = [n for n in sources if sums[n] == n]
        require(abundant == expected_abundant[k], f"abundance row {k}")
        require(perfect == expected_perfect[k], f"perfection row {k}")
        require(not any(n % 2 for n in abundant + perfect), f"odd nondeficient row {k}")
        rows[str(k)] = {"sources": len(sources), "abundant_count": len(abundant), "perfect": perfect}

    finite = [12,18,20,24,36,40,48,54,56,70,80,88,100,104,
              112,162,176,196,208,272,304,368,464]
    require(len(finite) == len(set(finite)) == 23, "finite exception count")
    require(not set(finite) & set(family), "disjoint finite/infinite pieces")
    actual = [n for n in range(1, limit + 1) if counts[n] <= 9 and sums[n] > n]
    require(actual == sorted(finite + family), "compressed abundance classification")

    pi = [0] * (limit + 1)
    for p in primes:
        pi[p] = 1
    for x in range(1, limit + 1):
        pi[x] += pi[x - 1]
    running = 0
    for x in range(limit + 1):
        if x > 0 and counts[x] <= 9 and sums[x] > x:
            running += 1
        require(running == max(pi[x // 6] - 2, 0) + sum(e <= x for e in finite), f"count at {x}")
        if x >= 464:
            require(running == pi[x // 6] + 21, f"simplified count at {x}")

    controls = {}
    for n in (1, 6, 28, 36, 100, 162, 196, 210, 496, 8128):
        ds = pairs(n)
        require(sum(ds) == sums[n] and len(ds) == counts[n], f"pair control {n}")
        controls[str(n)] = {"s": sums[n], "terms": counts[n], "excess": sums[n] - n}
    require(controls["210"] == {"s": 366, "terms": 15, "excess": 156}, "210 excluded")
    require(controls["8128"] == {"s": 8128, "terms": 13, "excess": 0}, "8128 outside atlas")
    print(json.dumps({"status": "passed", "source_universe": [1, limit],
                      "term_layers": rows, "abundant_total": len(actual),
                      "finite_abundant_exceptions": finite,
                      "all_counting_cutoffs": [0, limit], "controls": controls},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
