#!/usr/bin/env python3
"""Exact audit for THM-4056.

The script checks the divisor-phase bijection, its denominator strata and
Ramanujan transform, the LCM prime-power filtration, the finite
Duffin--Schaeffer first-moment compiler, and the opposite-parity
Pythagorean/Berggren subpacket.  All truth-bearing checks use ``check`` rather
than ``assert`` so optimized mode runs the same gates.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import gcd, lcm

from sympy import Poly, QQ, Rational, cyclotomic_poly, symbols


GATES = 0


def check(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


@lru_cache(maxsize=None)
def divisors(n: int) -> tuple[int, ...]:
    lo: list[int] = []
    hi: list[int] = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            lo.append(d)
            if d * d != n:
                hi.append(n // d)
        d += 1
    return tuple(lo + hi[::-1])


@lru_cache(maxsize=None)
def units(q: int) -> tuple[int, ...]:
    # The unique residue modulo 1 is the identity of the trivial unit group.
    if q == 1:
        return (0,)
    return tuple(a for a in range(q) if gcd(a, q) == 1)


def phi(q: int) -> int:
    return len(units(q))


def mobius(n: int) -> int:
    if n == 1:
        return 1
    p = 2
    parity = 0
    while p * p <= n:
        if n % p == 0:
            n //= p
            parity ^= 1
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        parity ^= 1
    return -1 if parity else 1


def ramanujan(q: int, k: int) -> int:
    return sum(d * mobius(q // d) for d in divisors(gcd(q, k)))


def encode(n: int, q: int, a: int) -> int:
    check(q >= 1 and n % q == 0, "encode divisor")
    check((q == 1 and a == 0) or (0 <= a < q and gcd(a, q) == 1),
          "encode unit")
    return ((n // q) * a) % n


def decode(n: int, x: int) -> tuple[int, int]:
    x %= n
    g = gcd(x, n)
    q = n // g
    a = 0 if q == 1 else x // g
    return q, a


def exact_denominator(n: int, x: int) -> int:
    return n // gcd(x % n, n)


def odd_part(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def prime_power(n: int) -> tuple[int, int] | None:
    if n < 2:
        return None
    p = next((d for d in range(2, n + 1) if n % d == 0), None)
    if p is None:
        return None
    e = 0
    m = n
    while m % p == 0:
        m //= p
        e += 1
    return (p, e) if m == 1 else None


def valuation(n: int, p: int) -> int:
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def lcm_upto(q: int) -> int:
    ans = 1
    for n in range(1, q + 1):
        ans = lcm(ans, n)
    return ans


def union_length(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    ordered = sorted(intervals)
    merged: list[list[Fraction]] = []
    for left, right in ordered:
        check(Fraction(0) <= left <= right <= Fraction(1), "interval domain")
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return sum((right - left for left, right in merged), Fraction(0))


print("THM-4056 exact divisor-phase / Duffin--Schaeffer clock audit")
print()

# 1. Bijection, inverse, strata, and unit-orbit symmetry.
tested_n = list(range(1, 161)) + [210, 420, 840, 2520, 27720]
compiler_images = 0
for n in tested_n:
    domain = [(q, a) for q in divisors(n) for a in units(q)]
    image = [encode(n, q, a) for q, a in domain]
    check(len(domain) == n, f"totient partition N={n}")
    check(len(set(image)) == n, f"compiler injective N={n}")
    check(set(image) == set(range(n)), f"compiler surjective N={n}")
    compiler_images += len(image)
    for q, a in domain:
        x = encode(n, q, a)
        check(decode(n, x) == (q, a), f"inverse N={n},q={q},a={a}")
        check(gcd(x, n) == n // q, f"gcd stratum N={n},q={q},a={a}")
    for q in divisors(n):
        stratum = {x for x in range(n) if exact_denominator(n, x) == q}
        check(len(stratum) == phi(q), f"stratum phi N={n},q={q}")
        representative = encode(n, q, units(q)[0])
        orbit = {(u * representative) % n for u in units(n)}
        check(orbit == stratum, f"unit orbit N={n},q={q}")

print(f"compiler universes: {len(tested_n)}")
print(f"compiled exact phases: {compiler_images}")
print("bijection/inverse/phi strata/unit orbits: PASS")
print()

# 2. Arbitrary weights and the exact Duffin--Schaeffer first-moment mass.
for n in list(range(1, 101)) + [420, 2520]:
    weight = {q: Fraction(q * q + 3, 2 * q + 1) for q in divisors(n)}
    left = sum((phi(q) * weight[q] for q in divisors(n)), Fraction(0))
    right = sum((weight[exact_denominator(n, x)] for x in range(n)), Fraction(0))
    check(left == right, f"weighted compiler N={n}")

    psi = {q: Fraction(1, 2 * (q + 1)) for q in divisors(n)}
    ds_left = sum(
        (Fraction(2 * phi(q), q) * psi[q] for q in divisors(n)),
        Fraction(0),
    )
    ds_right = sum(
        (Fraction(2, exact_denominator(n, x)) * psi[exact_denominator(n, x)]
         for x in range(n)),
        Fraction(0),
    )
    check(ds_left == ds_right, f"DS mass compiler N={n}")

print("arbitrary rational weights and DS multiplicity mass: PASS")
print()

# 3. Exact Ramanujan transform, verified in cyclotomic quotient rings.
X = symbols("X")
fourier_cases = 0
for n in list(range(1, 25)) + [30, 36, 42, 60]:
    cyclo = Poly(cyclotomic_poly(n, X), X, domain=QQ)
    weights = {q: Fraction(q * q + q + 1, q + 2) for q in divisors(n)}
    for k in range(n):
        direct = sum(
            (Rational(weights[exact_denominator(n, x)].numerator,
                      weights[exact_denominator(n, x)].denominator)
             * X ** ((k * x) % n)
             for x in range(n)),
            Rational(0),
        )
        predicted_fraction = sum(
            (weights[q] * ramanujan(q, k) for q in divisors(n)),
            Fraction(0),
        )
        predicted = Rational(predicted_fraction.numerator, predicted_fraction.denominator)
        remainder = Poly(direct - predicted, X, domain=QQ).rem(cyclo)
        check(remainder.is_zero, f"Ramanujan DFT N={n},k={k}")
        fourier_cases += 1

print(f"exact cyclotomic Ramanujan transforms: {fourier_cases}")
print("DFT formula sum_{d|N} W(d)c_d(k): PASS")
print()

# 4. Prefix filtering and prime-power growth of LCM clocks.
prefix_rows: list[tuple[int, int, int, int]] = []
previous = 1
for qmax in range(1, 17):
    n = lcm_upto(qmax)
    prefix_domain = {
        encode(n, q, a)
        for q in range(1, qmax + 1)
        for a in units(q)
    }
    filtered = {x for x in range(n) if exact_denominator(n, x) <= qmax}
    check(prefix_domain == filtered, f"DS prefix filter Q={qmax}")
    check(len(filtered) == sum(phi(q) for q in range(1, qmax + 1)),
          f"DS prefix count Q={qmax}")

    embedded = {((n // previous) * x) % n for x in range(previous)}
    new = set(range(n)) - embedded
    if n == previous:
        check(not new, f"no LCM jump Q={qmax}")
        check(prime_power(qmax) is None or qmax == 1,
              f"nonjump is not prime power Q={qmax}")
    else:
        pp = prime_power(qmax)
        check(pp is not None, f"LCM jump prime power Q={qmax}")
        p, e = pp
        check(n == p * previous, f"LCM jump ratio Q={qmax}")
        check(len(new) == n - previous, f"new phase count Q={qmax}")
        for x in new:
            d = exact_denominator(n, x)
            check(valuation(d, p) == e, f"new valuation depth Q={qmax},x={x}")
            check(previous % d != 0, f"new denominator Q={qmax},x={x}")
    prefix_rows.append((qmax, n, len(filtered), len(new)))
    previous = n

print("Q  L_Q  prefix-centres  new-clock-phases")
for row in prefix_rows:
    print("%2d %6d %14d %16d" % row)
print("prefix filter and prime-power jump law: PASS")
print()

# The clocks form a phase-preserving direct system by multiplication, not by
# ordinary residue reduction.  The phase 1/2 is 3 in C_6 and 6 in C_12;
# reducing 6 modulo 6 would incorrectly return the zero phase.
old_half = encode(6, 2, 1)
new_half = encode(12, 2, 1)
check((old_half, new_half) == (3, 6), "half-phase direct embedding")
check(decode(6, old_half) == decode(12, new_half) == (2, 1),
      "direct embedding preserves reduced phase")
check(new_half % 6 != old_half, "ordinary reduction hostile")
print("clock inclusions use x -> (N'/N)x; ordinary residue reduction fails: PASS")
print()

# 5. Correct circle reflection, and the reciprocal hostile.
for n in list(range(2, 101)) + [420]:
    for x in range(n):
        check(exact_denominator(n, -x) == exact_denominator(n, x),
              f"cyclic reflection N={n},x={x}")
q0, a0 = 5, 2
reciprocal_mod_one = Fraction(q0, a0) % 1
check(reciprocal_mod_one == Fraction(1, 2), "reciprocal hostile value")
check(reciprocal_mod_one.denominator != q0, "reciprocal changes DS denominator")
print("x -> -x preserves denominator/radius; a/q -> q/a need not: PASS")
print()

# 6. Opposite-parity phases are exactly primitive Pythagorean spinors.
masters = [6, 60, 420, 27720]
print("N  odd-part  primitive-Pythagorean-phases  formula")
for n in masters:
    eligible: list[tuple[int, int]] = []
    triples: set[tuple[int, int, int]] = set()
    for x in range(n):
        q, a = decode(n, x)
        if q > 1 and a > 0 and (q - a) % 2 == 1:
            eligible.append((a, q))
            triple = (q * q - a * a, 2 * q * a, q * q + a * a)
            check(triple[0] * triple[0] + triple[1] * triple[1]
                  == triple[2] * triple[2], f"Pythagorean N={n},a/q={a}/{q}")
            check(gcd(gcd(triple[0], triple[1]), triple[2]) == 1,
                  f"primitive triple N={n},a/q={a}/{q}")
            triples.add(triple)
    formula = n - (odd_part(n) + 1) // 2
    check(len(eligible) == formula, f"opposite parity count N={n}")
    check(len(triples) == len(eligible), f"spinor injective N={n}")
    print(f"{n:5d} {odd_part(n):9d} {len(eligible):29d} {formula:8d}")
print("Pythagorean/Berggren subpacket count: PASS")
print()

# 7. Hostile controls: divisor completion is not a prefix, and mass is not union.
n = lcm_upto(5)
extra_denominators = [d for d in divisors(n) if d > 5]
check(n == 60 and extra_denominators, "unfiltered LCM has extra denominators")
check(60 in extra_denominators, "largest extra denominator")

# q=2 has centre 1/2, radius 1/4; q=4 has centres 1/4,3/4,
# radius 1/8.  Layer mass is 1, but cross-layer overlap makes union 3/4.
intervals = [
    (Fraction(1, 4), Fraction(3, 4)),
    (Fraction(1, 8), Fraction(3, 8)),
    (Fraction(5, 8), Fraction(7, 8)),
]
formal_mass = sum((b - a for a, b in intervals), Fraction(0))
actual_union = union_length(intervals)
check(formal_mass == 1, "hostile formal mass")
check(actual_union == Fraction(3, 4), "hostile union mass")
check(formal_mass != actual_union, "mass is not coverage")
print(f"Q=5 unfiltered divisor denominators above Q: {extra_denominators}")
print(f"overlap hostile: formal mass={formal_mass}, union measure={actual_union}")
print("scope hostiles: PASS")
print()

milestones = [(q, lcm_upto(q)) for q in (3, 5, 7, 11)]
check([n for _, n in milestones] == masters, "6-60-420-27720 LCM milestones")
print(f"LCM odd-prime milestones: {milestones}")
print(f"TOTAL EXACT GATES: {GATES}")
print("ALL CHECKS PASSED")
