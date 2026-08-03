#!/usr/bin/env python3
"""Independent exact six-exit/first-flag audit for THM-3180.

The polynomial recurrence, primality/factorization filters, valuations, and
Newton lower hull are implemented here from scratch.  This artifact imports
neither the primary THM-3152 engine nor its hull code.
"""

import hashlib
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd

import gmpy2
from flint import fmpz_poly


BANK = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
WORK_CHUNKS = ((2001, 2050), (2051, 2100), (2101, 2150), (2151, 2200))
OLD_REPORT_CHUNKS = (
    (1001, 1100),
    (1101, 1300),
    (1301, 1500),
    (1501, 1750),
    (1751, 2000),
)
MAX_WORKERS = 4
V = fmpz_poly([0, 1])


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def prime(n):
    if n < 2:
        return False
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            return n == divisor
        divisor += 1
    return True


def distinct_prime_factor_count(n):
    count = 0
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            count += 1
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    return count + int(n > 1)


def prime_power(n):
    return distinct_prime_factor_count(n) == 1


def exit_predicates():
    return (
        lambda d: prime(d),
        lambda d: prime_power(d - 1),
        lambda d: prime(d - 2),
        lambda d: prime(d - 3),
        lambda d: prime(d - 4),
        lambda d: prime(d - 5),
    )


def exit_counts(start, end):
    alive = list(range(start, end + 1))
    counts = []
    for exits in exit_predicates():
        alive = [d for d in alive if not exits(d)]
        counts.append(len(alive))
    return tuple(counts), tuple(alive)


def six_exit_residual(d):
    return all(not exits(d) for exits in exit_predicates())


def seven_exit_residual(d):
    return six_exit_residual(d) and not prime(d - 6)


def moments(d):
    previous = fmpz_poly([1])
    current = fmpz_poly([d - 1, 2])
    d_to_k = d
    for k in range(1, d - 1):
        following = (
            fmpz_poly([d_to_k * (d - k - 1)])
            + 2 * (k + 1) * (2 * k + 1) * V * current
            + k * (k + 1) * (previous - 4 * d * V * previous)
        )
        previous, current = current, following
        d_to_k *= d
    return previous, current


def first_full_remainder(p, q, d):
    n = d - 2
    return (2 * n - 1) * (q - 2 * (n + 1) * (2 * n + 1) * V * p) + 2 * d * (n + 1) * p


def valuation(integer, rational_prime):
    integer = gmpy2.mpz(integer)
    if not integer:
        return None
    return int(gmpy2.remove(integer, rational_prime)[1])


def determinant_lower_profile(poly, rational_prime):
    hull = []
    zero_order = None
    for exponent in range(poly.degree() + 1):
        ordinate = valuation(poly[exponent], rational_prime)
        if ordinate is None:
            continue
        if zero_order is None:
            zero_order = exponent
        point = (exponent, ordinate)
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
    require(zero_order is not None, ("zero polynomial", rational_prime))

    capacities = {}
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        delta_x = x1 - x0
        delta_y = y1 - y0
        common = gcd(abs(delta_y), delta_x)
        slope = (delta_y // common, delta_x // common)
        require(slope not in capacities, ("repeated lower-hull slope", slope, hull))
        capacities[slope] = delta_x
    return zero_order, capacities


def allowed_degrees(rows, rational_prime):
    profiles = tuple(determinant_lower_profile(row, rational_prime) for row in rows)
    possibilities = set(range(min(zero for zero, _ in profiles) + 1))
    common_slopes = set(profiles[0][1])
    for _, capacities in profiles[1:]:
        common_slopes &= set(capacities)
    for slope in sorted(common_slopes):
        denominator = slope[1]
        capacity = min(capacities[slope] for _, capacities in profiles)
        require(capacity % denominator == 0, (slope, capacity, denominator))
        possibilities = {
            old + use
            for old in possibilities
            for use in range(0, capacity + 1, denominator)
        }
    return possibilities


def flag_trace(rows):
    possible = set(range(1, min(row.degree() for row in rows) + 1))
    trace = []
    for rational_prime in BANK:
        possible &= allowed_degrees(rows, rational_prime)
        trace.append((rational_prime, tuple(sorted(possible))))
        if not possible:
            break
    return possible, tuple(trace)


def positive_controls():
    nonzero_rows = (
        (V + 1) * (V**3 + 2 * V + 3),
        (V + 1) * (2 * V**4 + V**2 + 5),
        (V + 1) * (3 * V**2 + V + 7),
    )
    zero_rows = (
        V * (V**3 + 2 * V + 3),
        V * (2 * V**4 + V**2 + 5),
        V * (3 * V**2 + V + 7),
    )
    for label, rows in (("v+1", nonzero_rows), ("v", zero_rows)):
        possible = set(range(1, min(row.degree() for row in rows) + 1))
        for rational_prime in BANK:
            possible &= allowed_degrees(rows, rational_prime)
        require(1 in possible, (label, tuple(sorted(possible))))
    return "planted(v+1) and planted(v), degree 1 retained"


def scan_chunk(chunk):
    start, end = chunk
    rows = []
    residuals = []
    for d in range(start, end + 1):
        if not six_exit_residual(d):
            continue
        residuals.append(d)
        p, q = moments(d)
        r = first_full_remainder(p, q, d)
        require(
            (p.degree(), q.degree(), r.degree()) == (d - 2, d - 1, d - 3),
            (d, p.degree(), q.degree(), r.degree()),
        )
        possible, trace = flag_trace((p, q, r))
        require(not possible, (d, tuple(sorted(possible)), "six-exit survivor"))
        rows.append((d, trace))
    require(tuple(residuals) == tuple(d for d, _ in rows), (chunk, residuals, rows))
    return {
        "chunk": chunk,
        "residuals": tuple(residuals),
        "rows": tuple(rows),
        "semantic_trace_digest": digest(rows),
    }


def main():
    control_message = positive_controls()

    expected_censuses = {
        (3, 1000): (831, 642, 513, 390, 301, 211),
        (1001, 2000): (865, 725, 617, 511, 426, 341),
        (3, 2000): (1696, 1367, 1130, 901, 727, 552),
        (2001, 2200): (176, 149, 130, 111, 95, 79),
    }
    censuses = {}
    residual_lists = {}
    for bounds, expected in expected_censuses.items():
        counts, residuals = exit_counts(*bounds)
        require(counts == expected, (bounds, counts, expected))
        censuses[bounds] = counts
        residual_lists[bounds] = residuals

    old_report_progression = tuple(
        exit_counts(start, end)[0][2:] for start, end in OLD_REPORT_CHUNKS
    )
    require(
        old_report_progression
        == ((56, 45, 38, 30), (123, 100, 79, 58), (121, 101, 86, 72),
            (154, 125, 104, 83), (163, 140, 119, 98)),
        old_report_progression,
    )
    new_chunk_progression = tuple(
        exit_counts(start, end)[0][2:] for start, end in WORK_CHUNKS
    )
    require(
        new_chunk_progression
        == ((31, 25, 20, 15), (29, 24, 20, 16),
            (31, 26, 22, 18), (39, 36, 33, 30)),
        new_chunk_progression,
    )

    seven_exit_rows = tuple(
        d for d in residual_lists[(2001, 2200)] if seven_exit_residual(d)
    )
    require(len(seven_exit_rows) == 66, len(seven_exit_rows))
    first_unaudited = next(
        d for d in range(2201, 10000) if seven_exit_residual(d)
    )
    require(first_unaudited == 2201, first_unaudited)
    require(
        (2201, 2200, 2199, 2198, 2197, 2196, 2195)
        == (31 * 71, 2**3 * 5**2 * 11, 3 * 733, 2 * 7 * 157,
            13**3, 2**2 * 3**2 * 61, 5 * 439),
        "first unaudited factorization invoice",
    )

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        chunk_results = list(pool.map(scan_chunk, WORK_CHUNKS))

    residuals = tuple(d for result in chunk_results for d in result["residuals"])
    rows = tuple(row for result in chunk_results for row in result["rows"])
    require(residuals == residual_lists[(2001, 2200)], (residuals, residual_lists[(2001, 2200)]))
    require(len(rows) == 79, len(rows))
    require(all(not trace[-1][1] for _, trace in rows), "a row did not close")

    killers = tuple((d, trace[-1][0]) for d, trace in rows)
    killer_histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    trace_2009 = next(trace for d, trace in rows if d == 2009)
    require(trace_2009[-1][1] == (), (2009, trace_2009))

    print("THM-3180 SIX-EXIT FIRST-EUCLIDEAN-FLAG INDEPENDENT AUDIT")
    print("implementation=self-contained recurrence + trial division + determinant lower hull")
    print("exit_order=(d prime,d-1 prime power,d-2 prime,d-3 prime,d-4 prime,d-5 prime)")
    print("census_3_1000=%s" % (censuses[(3, 1000)],))
    print("census_1001_2000=%s" % (censuses[(1001, 2000)],))
    print("census_3_2000=%s" % (censuses[(3, 2000)],))
    print("old_report_progression_after_exits_3_to_6=%s" % (old_report_progression,))
    print("census_2001_2200=%s" % (censuses[(2001, 2200)],))
    print("new_chunk_progression_after_exits_3_to_6=%s" % (new_chunk_progression,))
    print("new_six_exit_residual_count=%d residual_digest=%s" % (len(residuals), digest(residuals)))
    print("related_thm3176_seven_exit_residual_count=%d" % len(seven_exit_rows))
    print("bank=%s" % (BANK,))
    print("deterministic_work_chunks=%s max_workers=%d" % (WORK_CHUNKS, MAX_WORKERS))
    print("chunk_semantic_trace_digests=%s" % (tuple(result["semantic_trace_digest"] for result in chunk_results),))
    print("global_semantic_trace_digest=%s" % digest(rows))
    print("killer_histogram=%s max_killer_prime=%d" % (killer_histogram, max(prime for _, prime in killers)))
    print("hostile_control_d_2009_trace=%s" % (trace_2009,))
    print("positive_controls=%s" % control_message)
    print("closed=79 survivors=0")
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=2198 is nonzero in at least one slot")
    print("first_unaudited=r=2199 d=2201 seven_exit_residual=True")
    print("first_unaudited_factorization=(31*71,2^3*5^2*11,3*733,2*7*157,13^3,2^2*3^2*61,5*439)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
