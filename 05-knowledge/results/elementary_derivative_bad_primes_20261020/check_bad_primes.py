#!/usr/bin/env python3
"""Exact standard-library controls; no asserts, random sampling, or dependencies."""
import json
import math


def require(condition, context):
    if not condition:
        raise RuntimeError(context)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def factor(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d:
            d = 3 if d == 2 else d + 2
            continue
        exponent = 0
        while n % d == 0:
            n //= d
            exponent += 1
        out.append((d, exponent))
    if n > 1:
        out.append((n, 1))
    return out


def order_via_lagrange(c, p, r=1):
    modulus = p ** r
    d = (p - 1) * p ** (r - 1)
    prime_divisors = [q for q, _ in factor(p - 1)]
    if r > 1:
        prime_divisors.append(p)
    for q in prime_divisors:
        while d % q == 0 and pow(c, d // q, modulus) == 1:
            d //= q
    require(pow(c, d, modulus) == 1, ("order endpoint", c, p, r, d))
    return d


def order_by_multiplication(c, p):
    v = 1
    for d in range(1, p):
        v = v * c % p
        if v == 1:
            return d
    raise RuntimeError(("order missing", c, p))


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def quotient(c, ell):
    require(abs(c) >= 2, ("dilation", c))
    numerator, denominator = c ** ell + 1, c + 1
    require(numerator % denominator == 0, ("integral quotient", c, ell))
    return numerator // denominator


def check(c, ell, p, independent_order=True):
    require(is_prime(ell) and ell % 2 == 1, ("odd prime exponent", ell))
    require((c + 1) % ell != 0, ("exponent filter", c, ell))
    n = quotient(c, ell)
    require(n > 1 and n % 2 == 1, ("positive odd quotient", c, ell, n))
    require(math.gcd(n, c * (c + 1)) == 1, ("coprime", c, ell, n))
    require(n % ell == 1, ("exponent cannot divide", c, ell, n))
    require(is_prime(p) and n % p == 0, ("prime divisor", c, ell, n, p))
    require(p != ell and (p - 1) % (2 * ell) == 0, ("prime class", c, ell, p))
    d = order_via_lagrange(c, p)
    require(d == 2 * ell, ("exact order", c, ell, p, d))
    if independent_order:
        require(order_by_multiplication(c, p) == d,
                ("independent order", c, ell, p, d))
    require(pow(-c, odd_part(p - 1), p) == 1, ("negative odd subgroup", c, p))
    lifts = []
    for r in range(1, 5):
        dr = order_via_lagrange(c, p, r)
        require(dr % d == 0 and odd_part(dr) * 2 == dr,
                ("order valuation lifts", c, p, r, dr))
        require((p ** (r - 1)) % (dr // d) == 0,
                ("kernel order", c, p, r, dr))
        require(pow(-c, odd_part((p - 1) * p ** (r - 1)), p ** r) == 1,
                ("negative odd subgroup lift", c, p, r))
        lifts.append(dr)
    return {"c": c, "ell": ell, "quotient": n, "prime": p,
            "order": d, "prime_power_orders_r_1_to_4": lifts}


table_specs = [
    (2, 5, 11), (3, 3, 7), (4, 3, 13),
    (16, 3, 241), (256, 3, 97), (65536, 3, 193),
    (-2, 3, 7), (-3, 3, 13), (-4, 5, 11),
    (-16, 7, 29), (-256, 7, 29), (-65536, 7, 29),
]
table = [check(*spec) for spec in table_specs]
extra_factors = []
for c, factors in [(256, [97, 673]), (65536, [193, 22253377])]:
    require(math.prod(factors) == quotient(c, 3), ("factorization", c, factors))
    extra_factors.extend(check(c, 3, p) for p in factors)

# Every prime divisor, not just a selected witness, in a bounded exact universe.
all_factor_rows = []
for c in range(-12, 13):
    if c in (-1, 0, 1):
        continue
    seen = {}
    for ell in (3, 5, 7):
        if (c + 1) % ell == 0:
            continue
        n = quotient(c, ell)
        for p, multiplicity in factor(n):
            require(p not in seen, ("different exponents share prime", c, ell, p, seen))
            seen[p] = ell
            row = check(c, ell, p)
            row["factor_multiplicity"] = multiplicity
            all_factor_rows.append(row)

boundary_primes = [p for p in range(3, 500) if is_prime(p)]
for p in boundary_primes:
    require(order_via_lagrange(-1, p) == 2, ("minus-one boundary", p))
    require(order_via_lagrange(1, p) == 1, ("one boundary", p))
    require(pow(1, odd_part(p - 1), p) == 1, ("minus-one coset", p))
    require(pow(-1, odd_part(p - 1), p) != 1, ("one not in coset", p))

hostiles = []
for c, ell, p, expected in [(2, 3, 3, 2), (-4, 3, 3, 2), (3, 9, 7, 6)]:
    n = quotient(c, ell)
    require(n % p == 0, ("hostile divisor", c, ell, p))
    d = order_via_lagrange(c, p)
    require(d == expected and d != 2 * ell, ("hostile order", c, ell, p, d))
    hostiles.append({"c": c, "exponent": ell, "quotient": n, "prime": p, "order": d})

ring_controls = []
for c, p, r in [(2, 11, 2), (2, 11, 3), (-2, 7, 2), (-2, 7, 3),
                (4, 13, 2), (4, 13, 3), (-16, 29, 2)]:
    n = p ** r
    h = [u for u in range(1, n) if pow(u, odd_part(p - 1), p) == 1]
    require(len(h) == odd_part((p - 1) * p ** (r - 1)),
            ("odd unit subgroup size", p, r))
    require((-c) % n in h, ("bad dilation lift", c, p, r))
    remaining = set(range(1, n))
    s = set()
    pairs = 0
    while remaining:
        x = min(remaining)
        orbit = {u * x % n for u in h}
        negative = {(-y) % n for y in orbit}
        require(orbit.isdisjoint(negative), ("antipodal orbit", p, r, x))
        require((orbit | negative) <= remaining, ("orbit partition", p, r, x))
        s.update(orbit)
        remaining.difference_update(orbit | negative)
        pairs += 1
    require(len(s) == (n - 1) // 2, ("tournament size", c, p, r))
    require(all(((x in s) != ((-x) % n in s)) for x in range(1, n)),
            ("tournament complement", c, p, r))
    require(all(((x in s) != ((c * x) % n in s)) for x in range(1, n)),
            ("complete derivative on all strata", c, p, r))
    nonunits = sum(x % p == 0 for x in s)
    require(nonunits == (p ** (r - 1) - 1) // 2,
            ("nonunit strata included", c, p, r, nonunits))
    ring_controls.append({"c": c, "prime": p, "exponent": r, "modulus": n,
                          "odd_subgroup_size": len(h), "orbit_pairs": pairs,
                          "tournament_size": len(s), "nonunits_in_tournament": nonunits,
                          "derivative_size": n - 1})

print(json.dumps({
    "status": "PASS",
    "method": "Trial-division primality and factorization; Lagrange order reduction, checked independently by repeated multiplication",
    "table": table,
    "complete_factorization_controls": extra_factors,
    "bounded_universe": {
        "integers": "[-12,12] excluding -1,0,1",
        "exponents": [3, 5, 7],
        "filter": "ell does not divide c+1",
        "scope": "Every distinct prime divisor of every admissible quotient",
        "prime_divisor_controls": len(all_factor_rows),
    },
    "all_factor_rows": all_factor_rows,
    "boundary_primes_count": len(boundary_primes),
    "hostiles": hostiles,
    "prime_power_scope": "Every positive control checked for r=1,2,3,4",
    "ring_circulant_controls_including_nonunits": ring_controls,
}, indent=2))
