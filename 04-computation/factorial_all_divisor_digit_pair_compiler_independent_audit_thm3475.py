#!/usr/bin/env python3
"""Independent formula-only referee for the THM-3475 divisor compiler."""

from fractions import Fraction
import hashlib
import json


START = 2501
END = 2600
EXPECTED_FILTER_COUNTS = (89, 78, 69, 60, 52, 44, 38)
EXPECTED_RESIDUALS = (
    2501, 2502,
    *range(2510, 2521),
    2528, 2529, 2530, 2538,
    *range(2564, 2579),
    *range(2586, 2591),
    2600,
)
EXPECTED_SURVIVORS = {
    2516: (503, 1006, 1509, 2012),
    2564: (466, 699, 1165, 1631, 1864, 2097, 2330),
    2571: (2056,),
    2576: tuple(range(103, 2473, 103)),
    2586: (
        47, 141, 188, 235, 282, 329,
        2209, 2256, 2303, 2350, 2397, 2444, 2491, 2538,
    ),
}
EXPECTED_SEMANTIC_SHA256 = (
    "e019eb61019620cfaffa8b5bb5769e8d171f08fb5ffeb2153f209c0128d42115"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def residual_universe():
    alive = list(range(START, END + 1))
    counts = []
    for shift in range(7):
        alive = [d for d in alive if not is_prime(d - shift)]
        counts.append(len(alive))
    return tuple(counts), tuple(alive)


def prime_divisors(n):
    answer = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            answer.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        answer.append(n)
    return tuple(answer)


def factorial_valuation(n, prime):
    answer = 0
    while n:
        n //= prime
        answer += n
    return answer


def binomial_valuation(n, j, prime):
    return (
        factorial_valuation(n, prime)
        - factorial_valuation(j, prime)
        - factorial_valuation(n - j, prime)
    )


def lower_hull(points):
    hull = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            left = (y1 - y0) * (x2 - x1)
            right = (y2 - y1) * (x1 - x0)
            if left < right:
                break
            hull.pop()
        hull.append(point)
    require(len(hull) >= 2, ("short hull", points))
    return tuple(hull)


def slope_ledger(hull):
    ledger = {}
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        slope = Fraction(y1 - y0, x1 - x0)
        require(slope not in ledger, ("repeated slope", slope, hull))
        ledger[slope] = x1 - x0
    return ledger


def raw_weight(n, j, prime):
    return binomial_valuation(n, j, prime) + factorial_valuation(2 * j, prime)


def g_hull(n, prime):
    return lower_hull(tuple((j, raw_weight(n, j, prime)) for j in range(n + 1)))


def f_suffix_hull(n, prime):
    if prime == 2:
        exponents = range(1, n, 2)
    else:
        exponents = range((prime - 1) // 2, n)
    return lower_hull(
        tuple((j, raw_weight(n - 1, j, prime)) for j in exponents)
    )


def audit_f_unit_sidecar(n, prime):
    if prime == 2:
        for j in range(1, n, 2):
            require((n - 1 - j) % 2 == 0, ("binary ell=1", n, j))
            require((2 * j + 2) % 2 == 0, ("binary rising", n, j))
        for j in range(0, n - 1, 2):
            require((n - 1 - j) % 2 == 1, ("binary cancellation binomial", n, j))
            require((2 * j + 1) % 2 == 1, ("binary cancellation rising", n, j))
        return

    half = (prime - 1) // 2
    for j in range(half, n):
        residue = j % prime
        if residue == half:
            require((2 * j + 1) % prime == 0, ("odd half anchor", n, prime, j))
        elif residue == prime - 1:
            require((n - 1 - j) % prime == 0, ("odd binomial anchor", n, prime, j))
            require((2 * j + 2) % prime == 0, ("odd rising anchor", n, prime, j))


def pair_barcode(n, prime):
    audit_f_unit_sidecar(n, prime)
    g_polygon = g_hull(n, prime)
    f_polygon = f_suffix_hull(n, prime)
    g_ledger = slope_ledger(g_polygon)
    f_ledger = slope_ledger(f_polygon)

    degrees = {0}
    blocks = []
    for slope in sorted(set(g_ledger) & set(f_ledger)):
        require(slope > 0, ("nonpositive common suffix slope", n, prime, slope))
        capacity = min(g_ledger[slope], f_ledger[slope])
        denominator = slope.denominator
        require(capacity % denominator == 0,
                ("nonintegral capacity", n, prime, slope, capacity))
        blocks.append((slope.numerator, denominator, capacity))
        degrees = {
            old + use
            for old in degrees
            for use in range(0, capacity + 1, denominator)
        }
    return tuple(sorted(degrees)), tuple(blocks), g_polygon, f_polygon


def main():
    filter_counts, residuals = residual_universe()
    require(filter_counts == EXPECTED_FILTER_COUNTS, filter_counts)
    require(residuals == EXPECTED_RESIDUALS, residuals)

    records = []
    for d in residuals:
        n = d - 1
        possible = set(range(1, n))
        trace = []
        for prime in prime_divisors(n):
            local, blocks, _, _ = pair_barcode(n, prime)
            possible &= set(local)
            trace.append((prime, blocks, tuple(sorted(possible))))
        records.append((d, tuple(trace), tuple(sorted(possible))))

    survivors = {d: possible for d, _, possible in records if possible}
    require(survivors == EXPECTED_SURVIVORS, survivors)
    require(len(records) - len(survivors) == 33, (len(records), survivors))

    pair_6, _, g_6, _ = pair_barcode(6, 3)
    g_6_ledger = slope_ledger(g_6)
    g_6_blocks = tuple(
        (slope.numerator, slope.denominator, g_6_ledger[slope])
        for slope in sorted(g_6_ledger)
    )
    require(g_6_blocks == ((2, 3, 3), (1, 1, 3)), g_6_blocks)
    require(pair_6 == (0, 3), pair_6)

    _, _, _, f_10 = pair_barcode(10, 2)
    require(all(exponent % 2 == 1 for exponent, _ in f_10), f_10)

    semantic_sha256 = digest(tuple(records))
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, semantic_sha256)

    print("THM-3475 ALL-DIVISOR DIGIT PAIR COMPILER INDEPENDENT AUDIT")
    print("implementation=fresh Legendre/Kummer weights, integer-cross-product hull, denominator DP; no theorem engine imports")
    print("universe=%s" % (residuals,))
    print("seven_exit_filter_counts=%s" % (filter_counts,))
    print("rows=38 pair_closed=33 survivor_rows=5")
    print("survivors=%s" % (survivors,))
    print("controls=(N,p)=(6,3) raw G blocks ((2/3,3),(1,3)) but pair degrees (0,3); (N,p)=(10,2) F suffix hull uses odd anchors only")
    print("sidecar=all odd half/right anchors and binary unit/cancellation congruences checked over the full universe")
    print("semantic_sha256=%s" % semantic_sha256)
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
