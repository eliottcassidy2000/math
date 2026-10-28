#!/usr/bin/env python3
"""Independent finite checks of elementary aliquot source-shape theorems."""

import argparse
import json
from math import isqrt, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    out = []
    p = 2
    while p * p <= n:
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        if a:
            out.append((p, a))
        p += 1 if p == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def divisor_pair_list(n):
    ds = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ds.update((d, n // d))
    ds.discard(n)
    return sorted(ds)


def classify(fs):
    if not fs:
        return "unit"
    if len(fs) == 1:
        return "prime_power"
    if len(fs) == 2 and all(a == 1 for _, a in fs):
        return "distinct_semiprime"
    return "higher_mixed"


def partitions(n, minimum=1):
    if n == 0:
        yield ()
        return
    for first in range(minimum, n + 1):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def multiplicative_partitions(n, minimum=2):
    if n >= minimum:
        yield (n,)
    for first in range(minimum, isqrt(n) + 1):
        if n % first == 0:
            for tail in multiplicative_partitions(n // first, first):
                yield (first,) + tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100000)
    args = parser.parse_args()
    require(38880 <= args.limit <= 2000000, "limit must lie in [38880, 2000000]")
    limit = args.limit
    sums, counts = [0] * (limit + 1), [0] * (limit + 1)
    for d in range(1, limit // 2 + 1):
        for n in range(2 * d, limit + 1, d):
            sums[n] += d
            counts[n] += 1

    census = {key: 0 for key in ("unit", "prime_power", "distinct_semiprime", "higher_mixed")}
    five_abundant, five_perfect, three_nondeficient = [], [], []
    even_mixed, odd_even_mixed = [], []
    for n in range(1, limit + 1):
        fs = factor(n)
        tau = prod(a + 1 for _, a in fs)
        sigma = prod(sum(p ** j for j in range(a + 1)) for p, a in fs)
        omega = sum(a for _, a in fs)
        shape = tuple(sorted(a for _, a in fs))
        cls = classify(fs)
        census[cls] += 1
        require(sigma - n == sums[n], f"sigma mismatch at {n}")
        require(tau - 1 == counts[n], f"tau mismatch at {n}")
        require((counts[n] % 2 == 0) == (isqrt(n) ** 2 == n), f"parity at {n}")
        if n > 1:
            require(omega <= counts[n] <= 2 ** omega - 1, f"complexity bound {n}")
            require(2 ** len(fs) - 1 <= counts[n], f"support bound {n}")
            require((counts[n] == omega) == (len(fs) == 1), f"lower equality {n}")
            require((counts[n] == 2 ** omega - 1) == all(a == 1 for _, a in fs), f"upper equality {n}")
        require((counts[n] == 3) == (shape in ((3,), (1, 1))), f"3-term class {n}")
        require((counts[n] == 5) == (shape in ((5,), (1, 2))), f"5-term class {n}")
        if cls == "higher_mixed":
            require(counts[n] >= 5 and sums[n] >= 16, f"mixed minimum {n}")
            if counts[n] % 2 == 0:
                even_mixed.append(n)
                if n % 2:
                    odd_even_mixed.append(n)
        if counts[n] == 5:
            if sums[n] > n:
                five_abundant.append(n)
            if sums[n] == n:
                five_perfect.append(n)
        if counts[n] == 3 and sums[n] >= n:
            three_nondeficient.append(n)
        if n <= 2000:
            ds = divisor_pair_list(n)
            require(len(ds) == counts[n] and sum(ds) == sums[n], f"paired divisors {n}")

    require(five_abundant == [12, 18, 20], "five-term abundance")
    require(five_perfect == [28], "five-term perfection")
    require(three_nondeficient == [6], "three-term nondeficiency")
    require(even_mixed[0] == 36 and odd_even_mixed[0] == 225, "mixed square minima")

    examples = {}
    for n in (6, 10, 15, 14, 12, 18, 20, 28, 24, 36, 225, 28800, 38880):
        ds, fs = divisor_pair_list(n), factor(n)
        require(len(ds) == counts[n] and sum(ds) == sums[n], f"example oracle {n}")
        examples[str(n)] = {
            "factorization": fs, "proper_divisors": ds, "s": sums[n],
            "terms": counts[n], "degree": sum(a for _, a in fs) - 1,
            "class": classify(fs), "excess": sums[n] - n,
        }
    require(examples["20"]["proper_divisors"] == [1, 2, 4, 5, 10], "20 hostile")
    require(sums[20] == 20 + 1 + 1 and 20 not in divisor_pair_list(20), "20 identity vs grammar")
    require(examples["28800"]["terms"] == examples["38880"]["terms"] == 71, "shape collision terms")
    require(examples["28800"]["degree"] == examples["38880"]["degree"] == 10, "shape collision degree")

    fibers = {}
    for target, expected in ((6, [6, 25]), (13, [27, 35]), (16, [12, 26])):
        bound = (target - 1) ** 2
        found = [n for n in range(1, bound + 1) if sums[n] == target]
        independent = [n for n in range(1, bound + 1) if sum(divisor_pair_list(n)) == target]
        require(found == independent == expected, f"complete F{target}")
        fibers[str(target)] = {"proved_source_bound": bound, "sources": found}

    # Independent additive partitions vs multiplicative partitions of tau.
    shape_layers = {}
    for terms in range(1, 32):
        multiplicative = {tuple(f - 1 for f in ps) for ps in multiplicative_partitions(terms + 1)}
        additive = {ps for omega in range(1, terms + 1) for ps in partitions(omega)
                    if prod(a + 1 for a in ps) == terms + 1}
        require(multiplicative == additive, f"multiplicative partition layer {terms}")
        require(any(len(ps) > 1 for ps in multiplicative) == (len(factor(terms + 1)) != 1 or factor(terms + 1)[0][1] != 1),
                f"mixed existence layer {terms}")
        shape_layers[str(terms)] = sorted(multiplicative)

    first_collisions = []
    for omega in range(1, 12):
        seen = {}
        for shape in partitions(omega):
            key = (len(shape), prod(a + 1 for a in shape))
            if key in seen:
                first_collisions.append({"Omega": omega, "omega": key[0], "tau": key[1],
                                         "shapes": [seen[key], shape]})
            else:
                seen[key] = shape
    require(first_collisions == [{"Omega": 11, "omega": 3, "tau": 72,
                                  "shapes": [(1, 5, 5), (2, 2, 7)]}], "summary collision minimum")

    print(json.dumps({
        "status": "passed", "source_universe": [1, limit],
        "independent_divisor_pair_universe": [1, 2000], "class_census": census,
        "five_term_abundant": five_abundant, "five_term_perfect": five_perfect,
        "three_term_nondeficient": three_nondeficient,
        "first_mixed_even_term_source": even_mixed[0],
        "first_odd_mixed_even_term_source": odd_even_mixed[0],
        "examples": examples, "complete_fibers": fibers, "term_count_shapes": shape_layers,
        "first_three_invariant_shape_collision": first_collisions,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
