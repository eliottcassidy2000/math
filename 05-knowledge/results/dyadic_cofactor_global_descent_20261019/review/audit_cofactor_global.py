#!/usr/bin/env python3
"""Independent exact controls for the dyadic cofactor-map classification.

Numerical checks supplement, and do not replace, the proof.
No floating-point arithmetic is used.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import random
import sympy as s


def F(n):
    assert n > 0 and n % 2 == 0
    c = n & -n
    p = n // c
    return (c - 1) * p + 2 * c - 1


def orbit(n, steps):
    out = [n]
    for _ in range(steps):
        out.append(F(out[-1]))
    return out


def fixed(n):
    c = n & -n
    return n // c == 2 * c - 1


EXCEPTIONS = {2, 4, 8, 12, 16}


def promised_descent_bound(n):
    c = n & -n
    q = n // c + 1
    if q > 2 * c:
        return 1, "immediate"
    if q == 2 * c:
        assert F(n) == n
        return None, "fixed"
    if n in EXCEPTIONS:
        return None, "exception"
    if q & (q - 1):
        return 2, "nonpower_q"
    if q == 2:
        return 3, "power_start"
    d = q
    if c == d:
        return 4, "resonance_c_d"
    if c < d * d:
        return 3, "middle"
    if c == d * d:
        return (6 if d == 4 else 5), "resonance_c_d2"
    return (4 if c == 2 * d * d else 3), "large"


def validate_point(n):
    bound, branch = promised_descent_bound(n)
    if bound is not None:
        values = orbit(n, bound)
        assert values[-1] < n, (n, branch, values)
    return branch


def symbolic_controls():
    c, d, p, e, h, C = s.symbols("c d p e h C")
    fcp = lambda cofactor, odd: (cofactor - 1) * odd + 2 * cofactor - 1
    n1 = fcp(c, p)
    n2 = fcp(d, n1 / d)
    equations = {
        "nonpower_second_difference": d * (n2 - c * p)
        - ((c + d) * (2 * d - p - 2) + p + 1),
        "power_start_first": fcp(C, 1) - (3 * C - 2),
        "power_start_second": fcp(2, (3 * C - 2) / 2) - (3 * C / 2 + 2),
        "power_start_third": fcp(2, (3 * C / 2 + 2) / 2) - (3 * C / 4 + 4),
    }
    n = e * d * (d - 1)
    n1 = fcp(e * d, d - 1)
    n2 = fcp(d, n1 / d)
    n3 = fcp(e, n2 / e)
    equations.update({
        "mersenne_first": n1 - d * (e * (d + 1) - 1),
        "mersenne_second": n2 - (e * (d * d - 1) + d),
        "middle_third": n3 - (e * d * d - d * d + e + d - d / e),
        "middle_difference": n3 - n - (-d * (d - e) + e + d - d / e),
    })
    n = h * d * d * (d - 1)
    n1 = fcp(h * d * d, d - 1)
    n2 = fcp(d, n1 / d)
    n3 = fcp(d, n2 / d)
    n4_h2 = fcp(d, n3.subs(h, 2) / d)
    equations.update({
        "large_second": n2 - d * (h * (d * d - 1) + 1),
        "large_third_difference": n3 - n - (-h * (d - 1) + 3 * d - 2),
        "h2_third": n3.subs(h, 2) - d * (2 * d * (d - 1) + 1),
        "h2_fourth_difference": n4_h2 - n.subs(h, 2) + (2 * d - 1) * (d - 2),
        "resonance_c_d_first": fcp(d, d - 1) - d * d,
        "resonance_c_d_fourth_difference": 3 * d * d / 4 + 4 - d * (d - 1)
        - (-d * d / 4 + d + 4),
        "resonance_c_d2_second": fcp(d, fcp(d * d, d - 1) / d) - d ** 3,
        "resonance_c_d2_fifth_difference": 3 * d ** 3 / 4 + 4 - d * d * (d - 1)
        - (-d ** 3 / 4 + d * d + 4),
    })
    for name, expression in equations.items():
        assert s.cancel(expression) == 0, (name, expression)
    return sorted(equations)


def valuation_controls():
    checked = 0
    for a in range(1, 181):
        c = 1 << a
        for t in range(1, 181):
            d = 1 << t
            n = c * (d - 1)
            validate_point(n)
            checked += 1
            if d >= 4 and c > d:
                v = orbit(n, 4)
                assert (v[1] & -v[1]) == d
                if c < d * d:
                    assert (v[2] & -v[2]) == c // d
                elif c == d * d:
                    assert v[2] == d ** 3
                else:
                    assert (v[2] & -v[2]) == d
                    if c == 2 * d * d:
                        assert (v[3] & -v[3]) == d
            elif d >= 4 and c == d:
                assert F(n) == d * d
        for t in range(1, a + 1):
            d = 1 << t
            for k in (3, 5, 7, 11, (2 * c // d) - 1, (2 * c // d) + 1):
                if k < 3 or k % 2 == 0:
                    continue
                q = d * k
                n = c * (q - 1)
                validate_point(n)
                checked += 1
                if q < 2 * c:
                    assert (F(n) & -F(n)) == d
    for a in range(3, 601):
        C = 1 << a
        v = orbit(C, 3)
        assert v == [C, 3 * C - 2, 3 * C // 2 + 2, 3 * C // 4 + 4]
        assert (v[1] & -v[1]) == (v[2] & -v[2]) == 2
        checked += 1
    rng = random.Random(20261019)
    for _ in range(20000):
        a = rng.randrange(1, 2049)
        bits = rng.randrange(1, 2049)
        p = rng.getrandbits(bits) | 1
        validate_point((1 << a) * p)
        checked += 1
    return checked


def canonical_cycle(cycle):
    j = cycle.index(min(cycle))
    return tuple(cycle[j:] + cycle[:j])


def independent_cycle_scan(limit):
    # Follow actual F orbits without using the claimed descent bound or cases.
    cache = {}
    cycles = set()
    longest_new_path = 0
    for start in range(2, limit + 1, 2):
        if start in cache:
            continue
        path, positions = [], {}
        n = start
        while n not in cache and n not in positions:
            positions[n] = len(path)
            path.append(n)
            n = F(n)
        longest_new_path = max(longest_new_path, len(path))
        if n in positions:
            cyc = canonical_cycle(path[positions[n]:])
            cycles.add(cyc)
            assert len(cyc) == 1 and fixed(cyc[0]) or cyc in {
                (8, 22, 14, 10), (16, 46, 26)
            }, cyc
        else:
            cyc = cache[n]
        for node in path:
            cache[node] = cyc
    return {
        "largest_start": limit,
        "distinct_nodes_visited": len(cache),
        "largest_node_visited": max(cache),
        "longest_uncached_path": longest_new_path,
        "cycles": [list(c) for c in sorted(cycles)],
    }


def independent_aliquot(n):
    # Independent divisor enumeration, including sqrt correction.
    total = 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            total += d
            if d * d != n:
                total += n // d
        d += 1
    return total


def native_guard_controls(limit):
    checked = 0
    for n in range(2, limit + 1, 2):
        c = n & -n
        p = n // c
        if s.isprime(p):
            assert F(n) == independent_aliquot(n), n
            checked += 1
    assert F(120) == 120 and independent_aliquot(120) == 240
    assert F(8) == 22 and independent_aliquot(8) == 7
    assert F(16) == 46 and independent_aliquot(16) == 15
    return checked


def main():
    source = Path(__file__).resolve().parents[1] / "COFACTOR_MAP_GLOBAL_CLASSIFICATION.md"
    result = {"source_sha256": sha256(source.read_bytes()).hexdigest()}
    result["symbolic_identities"] = symbolic_controls()
    result["parametric_and_large_integer_checks"] = valuation_controls()
    counts = Counter()
    first_descent_counts = Counter()
    last_six_step_examples = []
    for n in range(2, 4_000_001, 2):
        counts[validate_point(n)] += 1
        if fixed(n) or n in EXCEPTIONS:
            continue
        v = n
        for j in range(1, 7):
            v = F(v)
            if v < n:
                first_descent_counts[j] += 1
                if j == 6:
                    last_six_step_examples.append(n)
                break
        else:
            raise AssertionError(("No six-step descent", n))
    result["exhaustive_limit"] = 4_000_000
    result["exhaustive_cases"] = dict(counts)
    result["least_descent_step_counts"] = dict(first_descent_counts)
    result["six_step_starts"] = last_six_step_examples
    assert orbit(48, 6) == [48, 76, 64, 190, 98, 52, 46]
    for n in EXCEPTIONS:
        assert min(orbit(n, 30)) == n
    result["independent_cycle_scan"] = independent_cycle_scan(500_000)
    result["native_guard_cases"] = native_guard_controls(30_000)
    out = Path(__file__).with_name("cofactor_global_audit_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
