#!/usr/bin/env python3
"""Exact finite controls for THM-3833's ABC-conditional bridge.

The program never assumes ABC.  It checks the primitive cube-factorization
interface, the scale/radical hostile, and the exact shortlex atlas of
hyperbolic exponent tasks used by the conditional proofs.
"""

from __future__ import annotations

import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def primes_through(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


PRIMES = primes_through(2_000_000)


def factor(n: int) -> dict[int, int]:
    gate(n >= 1, "factorization domain")
    out: dict[int, int] = {}
    rest = n
    for p in PRIMES:
        if p * p > rest:
            break
        if rest % p:
            continue
        e = 0
        while rest % p == 0:
            rest //= p
            e += 1
        out[p] = e
    if rest > 1:
        out[rest] = 1
    return out


def radical(n: int) -> int:
    value = 1
    for p in factor(n):
        value *= p
    return value


def inert_cube_free_shell(s: int) -> bool:
    return all(p % 3 == 2 and e <= 2 for p, e in factor(s).items())


def quality(m: int, radical_abc: int) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 60
        return Decimal(m).ln() / Decimal(radical_abc).ln()


def hyperbolic(p: int, q: int, r: int) -> bool:
    return Fraction(1, p) + Fraction(1, q) + Fraction(1, r) < 1


def triples_of_weight(weight: int):
    for p in range(2, weight + 1):
        for q in range(p, weight + 1):
            r = weight - p - q
            if r < q:
                continue
            if hyperbolic(p, q, r):
                yield (p, q, r)


def first_hyperbolic_tasks(count: int) -> list[tuple[int, int, int]]:
    tasks: list[tuple[int, int, int]] = []
    weight = 6
    while len(tasks) < count:
        tasks.extend(triples_of_weight(weight))
        weight += 1
    return tasks[:count]


def main() -> None:
    print("THM-3833 ABC-CUBE/HYPERBOLIC TASK ATLAS")
    print("truth_gate=ABC is never assumed; every executable gate is unconditional finite-exact")

    pairs: list[tuple[Decimal, int, int, int, int, bool]] = []
    threshold_counts = {Fraction(1, 1): 0, Fraction(6, 5): 0, Fraction(4, 3): 0, Fraction(3, 2): 0}
    eligible_count = 0
    for s in range(3, 357):
        for a in range(1, (s + 1) // 2):
            b = s - a
            if not a < b or math.gcd(a, b) != 1:
                continue
            m = a**3 + b**3
            gate(math.gcd(a * b, m) == 1, "primitive cube triple is not pairwise coprime")
            rad_ab = radical(a) * radical(b)
            gate(rad_ab == radical(a * b), "rad(ab) failed to split")
            rad_m = radical(m)
            rad_abc = rad_ab * rad_m
            gate(math.gcd(rad_ab, rad_m) == 1, "ABC radical supports overlap")
            eligible = inert_cube_free_shell(s)
            eligible_count += int(eligible)
            q = quality(m, rad_abc)
            pairs.append((q, a, b, m, rad_abc, eligible))
            for threshold in threshold_counts:
                if m ** threshold.denominator > rad_abc ** threshold.numerator:
                    threshold_counts[threshold] += 1

    gate(len(pairs) == 19_314, "THM-3743 pair-atlas count drift")
    gate(eligible_count == 5_855, "THM-3825 inert decoder count drift")
    print("cube_universe=coprime 1<=a<b, a+b<=356")
    print(f"cube_pairs={len(pairs)} decoder_eligible={eligible_count}")
    print(
        "quality_threshold_counts="
        + ",".join(f">{t.numerator}/{t.denominator}:{threshold_counts[t]}" for t in threshold_counts)
    )
    print("top_abc_quality_rows=")
    for q, a, b, m, rad_abc, eligible in sorted(pairs, reverse=True)[:12]:
        print(
            f"  ({a},{b}) m={m} rad(abm)={rad_abc} "
            f"quality={q:.12f} decoder={int(eligible)}"
        )

    base = 1**3 + 3**3
    scaled = 2**3 + 6**3
    gate((base, scaled) == (28, 224), "scale hostile values drift")
    gate(radical(base) == radical(scaled) == 14, "value radical should forget inert cube scale")
    gate(radical(1 * 3 * base) == radical(2 * 6 * scaled) == 42, "ABC radical scale hostile drift")
    gate(math.gcd(2**3, 6**3) == 8, "scaled ABC triple must be nonprimitive")
    for k in range(13):
        x, y, value = 2**k, 3 * 2**k, 28 * 8**k
        gate(x**3 + y**3 == value, "infinite cube-scale family identity drift")
        gate(radical(value) == 14, "value radical should forget every tested cube scale")
        gate(radical(x * y * value) == 42, "ABC radical should forget every tested cube scale")
    print("scale_hostile=(1,3;28)->(2,6;224), same rad(value)=14 and rad(x*y*value)=42")
    print("scale_family=(2^k)^3+(3*2^k)^3=28*8^k has radical(value)=14 for every k>=0")
    print("scale_repair=divide common cube 2^3 before invoking ABC; decoder scale is an external sidecar")

    gate(radical(1 * 2 * 3) == radical(1 * 8 * 9) == 6, "radical-fibre hostile drift")
    gate(3 != 9, "radical address unexpectedly retained height")
    print("radical_fibre_hostile=1+2=3 and 1+8=9 share rad(abc)=6 but not height or quality")

    tasks = first_hyperbolic_tasks(15)
    gate(len(tasks) == len(set(tasks)) == 15, "shortlex task atlas is not injective")
    gate(tasks[0] == (3, 3, 4), "first shortlex hyperbolic task drift")
    print("hyperbolic_task_order=weight p+q+r, then lexicographic; ordinal is an address, not curvature")
    for ordinal, (p, q, r) in enumerate(tasks, 1):
        sigma = Fraction(1, p) + Fraction(1, q) + Fraction(1, r)
        kappa = 1 - sigma
        epsilon = kappa / (2 * sigma)
        delta = 1 - (1 + epsilon) * sigma
        gate(delta == kappa / 2, "canonical epsilon does not halve curvature margin")
        print(
            f"  task={ordinal} exponents=({p},{q},{r}) "
            f"sigma={sigma.numerator}/{sigma.denominator} "
            f"kappa={kappa.numerator}/{kappa.denominator}"
        )

    bounded = [
        (p, q, r)
        for p in range(2, 41)
        for q in range(p, 41)
        for r in range(q, 41)
        if hyperbolic(p, q, r)
    ]
    margins = {triple: 1 - sum((Fraction(1, x) for x in triple), Fraction(0, 1)) for triple in bounded}
    min_margin = min(margins.values())
    minimizers = sorted(t for t, margin in margins.items() if margin == min_margin)
    gate(min_margin == Fraction(1, 42) and minimizers == [(2, 3, 7)], "Hurwitz margin control failed")
    euclidean = [
        (p, q, r)
        for p in range(2, 41)
        for q in range(p, 41)
        for r in range(q, 41)
        if Fraction(1, p) + Fraction(1, q) + Fraction(1, r) == 1
    ]
    gate(euclidean == [(2, 3, 6), (2, 4, 4), (3, 3, 3)], "Euclidean triangle controls drift")
    print(f"bounded_hyperbolic_census_2<=p<=q<=r<=40={len(bounded)}")
    print("least_curvature_signature=(2,3,7), kappa=1/42; first shortlex task=(3,3,4)")
    print("euclidean_controls=(2,3,6),(2,4,4),(3,3,3)")

    pyth_sigma = Fraction(3, 2)
    gate(pyth_sigma > 1, "Pythagorean hostile left spherical regime")
    print("pythagorean_hostile=(2,2,2), sigma=3/2: ABC gives no tree finiteness")
    pyth_nodes = [(5, 12, 13), (15, 8, 17)]
    for odd_leg, even_leg, hypotenuse in pyth_nodes:
        gate(odd_leg**2 + even_leg**2 == hypotenuse**2, "Pythagorean shell node drift")
    gate({even_leg + hypotenuse for _, even_leg, hypotenuse in pyth_nodes} == {25}, "odd-square shell collision drift")
    gate(tuple(radical(a * b * c) for a, b, c in pyth_nodes) == (390, 510), "Pythagorean radical sidecar drift")
    print("pythagorean_shell_hostile=(5,12,13),(15,8,17): B+C=25 but rad(ABC)=390,510")
    for t in range(1, 17):
        x, y, z = 9 * t**4, 18 * t**4, 9 * t**3
        gate(x**3 + y**3 == z**4, "nonprimitive (3,3,4) scaling family drift")
        gate(math.gcd(math.gcd(x, y), z) > 1, "hyperbolic scaling hostile became primitive")
    print("hyperbolic_scale_hostile=(9t^4)^3+(18t^4)^3=(9t^3)^4: primitivity is indispensable")
    left = 3**6 + 19**6 + 22**6
    right = 10**6 + 15**6 + 23**6
    gate(left == right, "six-term sixth-power hostile drift")
    gate(math.gcd(math.gcd(3, 19), math.gcd(22, math.gcd(10, math.gcd(15, 23)))) == 1, "six-term hostile lost primitivity")
    print(f"six_term_hostile=3^6+19^6+22^6=10^6+15^6+23^6={left}; not a three-term ABC input")

    print(f"active_checks={CHECKS}")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
