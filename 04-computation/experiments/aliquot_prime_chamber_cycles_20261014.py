#!/usr/bin/env python3
"""Exact aliquot-cycle, odd-cycle-block, and dyadic prime-chamber checks.

No network or third-party dependencies. Python 3.10+.
The main census is finite: 1 <= a,b,c <= 256. It checks all cyclic
representatives, hence all 256**3 ordered templates. Nothing here proves
the unrestricted absence of aliquot 3-cycles.
"""
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        answer[n] = 1
    return answer


def sigma(n):
    result = 1
    for p, a in factor(n).items():
        result *= (p ** (a + 1) - 1) // (p - 1)
    return result


def valuation(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def canonical(cycle):
    return min(tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle)))


def finite_cycles(bound):
    sums = [0] * (bound + 1)
    for d in range(1, bound // 2 + 1):
        for n in range(2 * d, bound + 1, d):
            sums[n] += d
    seen = bytearray(bound + 1)
    cycles = []
    for start in range(2, bound + 1):
        if seen[start]:
            continue
        path, positions = [], {}
        n = start
        while 1 < n <= bound and not seen[n] and n not in positions:
            positions[n] = len(path)
            path.append(n)
            n = sums[n]
        if n in positions:
            cycles.append(canonical(path[positions[n]:]))
        for m in path:
            seen[m] = 1
    cycles.sort()
    histogram = dict(sorted(Counter(map(len, cycles)).items()))
    require(histogram == {1: 4, 2: 40, 5: 1, 28: 1}, "finite-cycle census")
    traces = {k: sum(len(c) for c in cycles if k % len(c) == 0)
              for k in (1, 2, 3, 4, 5, 9, 28)}
    return {"bound_all_members": bound, "histogram": histogram,
            "traces": traces, "cycles": cycles}


def cycle_packet(cycle):
    k = len(cycle)
    sigmas = [sigma(n) for n in cycle]
    require(all(sigmas[i] - cycle[i] == cycle[(i + 1) % k]
                for i in range(k)), "displayed cycle has a false edge")
    g = reduce(gcd, cycle)
    h = reduce(gcd, sigmas)
    stable = 1
    for p, a in factor(g).items():
        if all(valuation(n, p) == a for n in cycle):
            stable *= p ** a
    if k % 2:
        require(h in (g, 2 * g), "odd neighbour-sum gcd")
        require((2 * g) % sigma(stable) == 0, "stable-prime supplier capacity")
        if k > 1 and g > 1:
            require(stable != g, "Dickson's common-unitary-gcd obstruction")
    product_ratio = Fraction(1)
    for i, n in enumerate(cycle):
        product_ratio *= Fraction(cycle[(i + 1) % k], n)
    require(product_ratio == 1, "cycle ratio product")
    return {"cycle": cycle, "length": k, "gcd_members": g,
            "gcd_sigma": h, "stable_kernel": stable,
            "ramified_common_kernel": g // stable,
            "sign_word": "".join("A" if sigmas[i] > 2 * n else
                                 "D" if sigmas[i] < 2 * n else "P"
                                 for i, n in enumerate(cycle)),
            "v2_members": [valuation(n, 2) for n in cycle],
            "v2_sigma": [valuation(n, 2) for n in sigmas]}


def dyadic_numerator_denominator(a, b, c):
    x, y, z = 1 << a, 1 << b, 1 << c
    d = x * y + y * z + z * x - x - y - z + 1
    n = ((z - 1) * (y - 1) * (2 * x - 1)
         + (z - 1) * y * (2 * y - 1) + y * z * (2 * z - 1))
    return n, d


def integer_closure(a, b, c):
    x, y, z = 1 << a, 1 << b, 1 << c
    n, d = dyadic_numerator_denominator(a, b, c)
    p, rem = divmod(n, d)
    if rem:
        return None
    q, rem = divmod((x - 1) * p + 2 * x - 1, y)
    if rem:
        return None
    r, rem = divmod((y - 1) * q + 2 * y - 1, z)
    if rem:
        return None
    require((z - 1) * r + 2 * z - 1 == x * p, "closed template")
    return p, q, r


def independent_fraction_closure(exponents):
    # Independent affine composition in n, not the numerator formula in p.
    slope, intercept = Fraction(1), Fraction(0)
    for a in exponents:
        t = 1 << a
        alpha = Fraction(t - 1, t)
        slope, intercept = alpha * slope, alpha * intercept + 2 * t - 1
    n = intercept / (1 - slope)
    values = []
    for a in exponents:
        t = 1 << a
        values.append(n / t)
        n = Fraction(t - 1, t) * n + 2 * t - 1
    return tuple(values)


def chamber_census(bound):
    representatives = 0
    closures = []
    powers = [1 << a for a in range(bound + 1)]
    for a in range(1, bound + 1):
        x = powers[a]
        for b in range(a, bound + 1):
            y = powers[b]
            for c in range(a, bound + 1):
                if (a, b, c) > min((b, c, a), (c, a, b)):
                    continue
                representatives += 1
                z = powers[c]
                d = x * y + y * z + z * x - x - y - z + 1
                n = ((z - 1) * (y - 1) * (2 * x - 1)
                     + (z - 1) * y * (2 * y - 1) + y * z * (2 * z - 1))
                p, rem = divmod(n, d)
                if rem:
                    continue
                q, rem = divmod((x - 1) * p + 2 * x - 1, y)
                if rem:
                    continue
                r, rem = divmod((y - 1) * q + 2 * y - 1, z)
                if rem:
                    continue
                require((z - 1) * r + 2 * z - 1 == x * p, "census closure")
                closures.append({"exponents": [a, b, c], "cofactors": [p, q, r]})
    expected = (bound ** 3 + 2 * bound) // 3
    require(representatives == expected, "Burnside representative count")
    heterogeneous = [r for r in closures if len(set(r["exponents"])) > 1]
    require(heterogeneous == [{"exponents": [1, 1, 4], "cofactors": [23, 13, 1]}],
            "heterogeneous integral closures")
    require(len(closures) == bound + 1, "closure count")
    for r in closures:
        if len(set(r["exponents"])) == 1:
            p = (1 << (r["exponents"][0] + 1)) - 1
            require(r["cofactors"] == [p, p, p], "constant closure")
    return {"exponent_bound_inclusive": bound,
            "ordered_templates_covered": bound ** 3,
            "cyclic_representatives_checked": representatives,
            "integral_closures_up_to_rotation": len(closures),
            "constant_closures": bound, "heterogeneous_closures": heterogeneous}


def independent_checks():
    fraction_checks = 0
    for exponents in product(range(1, 13), repeat=3):
        values = independent_fraction_closure(exponents)
        n, d = dyadic_numerator_denominator(*exponents)
        require(values[0] == Fraction(n, d), "independent affine formula")
        integral = tuple(int(v) for v in values) if all(v.denominator == 1 for v in values) else None
        require(integral == integer_closure(*exponents), "independent integer gate")
        fraction_checks += 1
    divisor_checks = 0
    for a in range(1, 65):
        x = 1 << a
        q = (2 * x * x - 1) * (3 * x * x - 3 * x + 1)
        for c in range(a, 3 * a + 3):
            z = 1 << c
            d = (2 * x - 1) * z + (x - 1) ** 2
            n, check_d = dyadic_numerator_denominator(a, a, c)
            require(d == check_d, "repeated-minimum denominator")
            ell = 6*x**3 + 4*x*x*z - 10*x*x - 2*x*z + 5*x - 1
            require((2*x - 1)**2 * n + x*q == d*ell, "divisor identity")
            require((q % d == 0) == (integer_closure(a, a, c) is not None),
                    "divisor gate equivalence")
            if q % d == 0:
                require(c <= 3*a + 1, "proved exponent bound")
            divisor_checks += 1
    require(integer_closure(1, 1, 4) == (23, 13, 1), "virtual loop")
    require([sigma(n)-n for n in (46, 26, 16)] == [26, 16, 15], "native prime gate")
    # I+C_k and its alternating inverse numerator, including the k=1 loop.
    matrix_checks = 0
    for k in range(1, 32):
        p = [[int(j == (i + 1) % k) for j in range(k)] for i in range(k)]
        m = [[int(i == j) + p[i][j] for j in range(k)] for i in range(k)]
        alternating = [[0] * k for _ in range(k)]
        for power in range(k):
            for i in range(k):
                alternating[i][(i + power) % k] += (-1) ** power
        product_matrix = [[sum(m[i][l] * alternating[l][j] for l in range(k))
                           for j in range(k)] for i in range(k)]
        expected = [[2 * int(i == j) if k % 2 else 0 for j in range(k)] for i in range(k)]
        require(product_matrix == expected, "odd/even cyclic block")
        matrix_checks += 1
    return {"independent_fraction_templates": fraction_checks,
            "repeated_minimum_divisor_checks": divisor_checks,
            "cycle_matrix_orders_checked": matrix_checks,
            "virtual_loop": {"template": [46, 26, 16, 46],
                             "actual_first_steps": [46, 26, 16, 15],
                             "failure": "cofactor 1 is not prime; sigma(1)=1, not 2"}}


def main():
    displayed = [
        [6], [220, 284],
        [1264460, 1547860, 1727636, 1305184],
        [12496, 14288, 15472, 14536, 14264],
        [805984760, 1268997640, 1803863720, 2308845400, 3059220620,
         3367978564, 2525983930, 2301481286, 1611969514],
        [14316, 19116, 31704, 47616, 83328, 177792, 295488, 629072,
         589786, 294896, 358336, 418904, 366556, 274924, 275444, 243760,
         376736, 381028, 285778, 152990, 122410, 97946, 48976, 45946,
         22976, 22744, 19916, 17716],
    ]
    result = {
        "status": "exact finite checks; unrestricted three-cycle problem remains open",
        "source_date": "2026-10-14",
        "cycle_data_source": "https://djm.cc/sociable.txt (dated 11-IX-2026)",
        "displayed_cycles": [cycle_packet(c) for c in displayed],
        "finite_aliquot_graph": finite_cycles(1_000_000),
        "dyadic_chambers": chamber_census(256),
        "independent_checks": independent_checks(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
