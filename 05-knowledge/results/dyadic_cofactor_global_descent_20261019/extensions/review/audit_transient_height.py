#!/usr/bin/env python3
"""Independent exact finite controls for transient-height scaling.

No floating-point evaluation of fractional powers is used. Explicit checks
remain active under python -O. Finite controls supplement the symbolic proof.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import random


def check(ok, details):
    if not ok:
        raise RuntimeError(details)


def F(n):
    c = n & -n
    return (c - 1) * (n // c) + 2 * c - 1


def orbit(n, length):
    a = [n]
    for _ in range(length):
        a.append(F(a[-1]))
    return a


EXCEPTIONS = {2, 4, 8, 12, 16}


def cube_root_floor(n):
    low, high = 0, 1 << ((n.bit_length() + 2) // 3)
    while high - low > 1:
        mid = (low + high) // 2
        if mid ** 3 <= n:
            low = mid
        else:
            high = mid
    return low


def integer_envelope(n):
    # This is at most the manuscript's real-power envelope, since r<=n^(1/3).
    r = cube_root_floor(n)
    return max(46, 3 * n + 3 * r * r + 6 * r + 1)


def local_excursion(n):
    c = n & -n
    q = n // c + 1
    check(F(n) <= 3 * n - 2, ("universal_first_step", n))
    if n in EXCEPTIONS:
        check(max(orbit(n, 15)) <= (22 if n < 12 else 46), ("exception", n))
        return "exception"
    if q == 2 * c:
        check(F(n) == n, ("fixed", n))
        return "fixed"
    if q > 2 * c:
        name, steps, ceiling = "immediate", 1, 3 * n - 2
    elif q & (q - 1):
        name, steps, ceiling = "nonpower_q", 2, 3 * n - 2
    elif q == 2:
        name, steps, ceiling = "power_start", 3, 3 * n - 2
    elif c == q:
        name, steps, ceiling = "first_resonance", 4, 3 * q * q - 2
    elif c == q * q:
        name, steps, ceiling = "second_resonance", 6 if q == 4 else 5, 3 * q ** 3 - 2
    elif c < q * q:
        name, steps, ceiling = "middle", 3, 3 * n - 2
    else:
        name, steps, ceiling = "high", 4 if c == 2 * q * q else 3, 3 * n - 2
    v = orbit(n, steps)
    check(v[-1] < n, ("checkpoint", n, name, v))
    check(max(v) <= ceiling, ("local_envelope", n, name, v, ceiling))
    if name in ("middle", "high"):
        check(3 * v[1] <= 5 * n, ("first_bound", n))
        check(v[2] <= 2 * n, ("second_bound", n))
    return name


def full_orbit_max(n, cache):
    # Follow the map directly, independently of the proved case decomposition.
    path, positions = [], {}
    x = n
    while x not in cache and x not in positions:
        positions[x] = len(path)
        path.append(x)
        x = F(x)
    if x in cache:
        height = cache[x]
    else:
        j = positions[x]
        height = max(path[j:])
        for v in path[j:]:
            cache[v] = height
        path = path[:j]
    for v in reversed(path):
        height = max(height, v)
        cache[v] = height
    return cache[n]


def main():
    source = Path(__file__).resolve().parents[1] / "TRANSIENT_HEIGHT_SCALING.md"
    result = {"source_sha256": sha256(source.read_bytes()).hexdigest()}
    counts = Counter()
    for n in range(2, 4_000_001, 2):
        counts[local_excursion(n)] += 1
    result["local_excursion_exhaustive_limit"] = 4_000_000
    result["local_excursion_cases"] = dict(counts)
    structured = 0
    for a in range(1, 181):
        c = 1 << a
        for t in range(1, 181):
            for odd in (1, 3, 5, 7, 11):
                n = c * ((1 << t) * odd - 1)
                local_excursion(n)
                structured += 1
    result["large_structured_excursion_cases"] = structured
    cache = {}
    for n in range(2, 1_000_001, 2):
        height = full_orbit_max(n, cache)
        check(height <= integer_envelope(n), ("full_height", n, height))
    result["full_orbit_exhaustive_limit"] = 1_000_000
    result["cached_distinct_nodes"] = len(cache)
    family = []
    for t in range(3, 257):
        d = 1 << t
        n = d * d * (d - 1)
        initial = orbit(n, 5)
        check(initial == [n, d * (d * (d + 1) - 1), d ** 3,
                          3 * d ** 3 - 2, 3 * d ** 3 // 2 + 2,
                          3 * d ** 3 // 4 + 4], ("resonant_chain", d))
        height = full_orbit_max(n, cache)
        check(height == 3 * d ** 3 - 2, ("exact_resonant_max", d, height))
        m = initial[-1]
        check(m < d ** 3, ("smaller_cube", d))
        # Four times the polynomial gap in the proposed exact-height proof.
        check(3 * d ** 3 - 12 * d * d - 24 * d - 60 > 0, ("strict_gap", d))
        check(integer_envelope(m) < 3 * d ** 3 - 2, ("tail_bound", d))
        family.append({"dyadic_exponent": t, "start_bits": n.bit_length()})
    result["exact_resonant_maximum_cases"] = len(family)
    result["largest_resonant_start_bits"] = family[-1]["start_bits"]
    rng = random.Random(20261019)
    large_full = 0
    for _ in range(300):
        a = rng.randrange(1, 1025)
        p = rng.getrandbits(rng.randrange(1, 1025)) | 1
        n = (1 << a) * p
        local_excursion(n)
        height = full_orbit_max(n, cache)
        check(height <= integer_envelope(n), ("large_full_height", n, height))
        large_full += 1
    result["random_large_full_orbit_cases"] = large_full
    result["status"] = "PASS"
    out = Path(__file__).with_name("transient_height_audit_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
