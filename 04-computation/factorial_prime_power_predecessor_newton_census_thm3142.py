#!/usr/bin/env python3
"""Exact controls for THM-3142 and composite resonances in its polynomial pair.

No repository imports are used.  Exact integer coefficients drive every
Newton polygon; the finite gcd census uses the proved THM-3124 recurrence
over two explicitly declared prime fields.
"""

from fractions import Fraction
from math import comb, factorial, isqrt

import sympy as sp
from flint import nmod_poly


EXACT_DIVISOR_LIMIT = 220
ALL_PLACE_LIMIT = 160
GCD_D_MIN = 203
GCD_D_MAX = 1000
GCD_PRIMES = (1_000_003, 1_000_033)


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % q for q in range(3, isqrt(n) + 1, 2))


def factor(n: int):
    answer = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            answer.append((p, e))
        p += 1 if p == 2 else 2
    if n > 1:
        answer.append((n, 1))
    return answer


def is_prime_power_at(n: int, p: int) -> bool:
    while n % p == 0:
        n //= p
    return n == 1


def vp(a: int, p: int):
    if a == 0:
        return None
    answer = 0
    while a % p == 0:
        a //= p
        answer += 1
    return answer


def coefficients(n: int, d: int):
    """[v^j] L((d-t+v*t^2)^n), exactly over Z."""
    facts = [factorial(k) for k in range(2 * n + 1)]
    answer = []
    for j in range(n + 1):
        inner = sum(
            comb(n - j, ell)
            * d ** (n - j - ell)
            * (-1) ** ell
            * facts[2 * j + ell]
            for ell in range(n - j + 1)
        )
        answer.append(comb(n, j) * inner)
    return answer


def lower_hull(coeffs, p: int):
    hull = []
    for x, a in enumerate(coeffs):
        if a == 0:
            continue
        point = (x, vp(a, p))
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            old = Fraction(y1 - y0, x1 - x0)
            new = Fraction(y2 - y1, x2 - x1)
            if old >= new:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def slopes(hull):
    return tuple(
        Fraction(y1 - y0, x1 - x0)
        for (x0, y0), (x1, y1) in zip(hull, hull[1:])
    )


def composite(n: int) -> bool:
    return n >= 4 and not is_prime(n)


def exact_place_census():
    divisor_places = 0
    divisor_overlaps = 0
    even_single_slope = 0
    predecessor_places = 0
    predecessor_separations = 0
    predecessor_classification = True
    composites = 0
    examples = {}
    for d in range(4, EXACT_DIVISOR_LIMIT + 1):
        if not composite(d):
            continue
        composites += 1
        a = coefficients(d - 2, d)
        b = coefficients(d - 1, d)
        for p, _ in factor(d):
            divisor_places += 1
            sa = slopes(lower_hull(a, p))
            sb = slopes(lower_hull(b, p))
            overlap = set(sa) & set(sb)
            divisor_overlaps += bool(overlap)
            require(overlap, f"missing divisor overlap d,p={d,p}")
            if p == 2:
                require(
                    sa == (Fraction(1),) and sb == (Fraction(1),),
                    f"binary slope-one gate d={d}",
                )
                even_single_slope += 1
            if (d, p) in ((4, 2), (6, 3), (132, 3)):
                examples[(d, p)] = (sa, sb, tuple(sorted(overlap)))
        for p, _ in factor(d - 1):
            predecessor_places += 1
            sa = slopes(lower_hull(a, p))
            sb = slopes(lower_hull(b, p))
            separated = set(sa).isdisjoint(sb)
            expected = is_prime_power_at(d - 1, p)
            predecessor_separations += separated
            predecessor_classification &= separated == expected
    require(divisor_places == divisor_overlaps, "divisor-place census")
    require(predecessor_classification, "predecessor-place classification")
    return (
        composites,
        divisor_places,
        even_single_slope,
        predecessor_places,
        predecessor_separations,
        examples,
    )


def prime_power_predecessor_controls():
    cases = []
    for p, max_k in ((2, 7), (3, 4), (5, 3), (7, 2), (11, 2), (13, 1)):
        for k in range(1, max_k + 1):
            n = p**k
            d = n + 1
            a = coefficients(n - 1, d)
            b = coefficients(n, d)
            ha = lower_hull(a, p)
            hb = lower_hull(b, p)
            lead = vp(factorial(2 * n), p)
            require(hb == [(0, 0), (n, lead)], f"prime-power hull p,k={p,k}")
            require(
                set(slopes(ha)).isdisjoint(slopes(hb)),
                f"prime-power separation p,k={p,k}",
            )
            require(lead % p != 0, f"primitive endpoint p,k={p,k}")
            cases.append((p, k, n, str(Fraction(lead, n))))
    return cases


def all_small_places_classification():
    tested = 0
    separating = 0
    classification = True
    for d in range(4, ALL_PLACE_LIMIT + 1):
        if not composite(d):
            continue
        a = coefficients(d - 2, d)
        b = coefficients(d - 1, d)
        for p in range(2, 2 * d + 1):
            if not is_prime(p):
                continue
            tested += 1
            separated = set(slopes(lower_hull(a, p))).isdisjoint(
                slopes(lower_hull(b, p))
            )
            expected = (d - 1) % p == 0 and is_prime_power_at(d - 1, p)
            separating += separated
            classification &= separated == expected
    require(classification, "all-small-place classification")
    return tested, separating


def recurrence_pair_mod(d: int, modulus: int):
    """Return A_(d-2), A_(d-1) via the division-free moment recurrence."""
    x = nmod_poly([0, 1], modulus)
    old = nmod_poly([1], modulus)
    current = nmod_poly([d - 1, 2], modulus)
    d_power = d % modulus
    for n in range(1, d - 1):
        boundary = d_power * (d - (n + 1)) % modulus
        following = (
            nmod_poly([boundary], modulus)
            + x * current * (2 * (n + 1) * (2 * n + 1))
            + (old + x * old * (-4 * d)) * (n * (n + 1))
        )
        old, current = current, following
        d_power = d_power * d % modulus
    return old, current


def recurrence_direct_controls():
    for d in range(4, 16):
        exact_a = coefficients(d - 2, d)
        exact_b = coefficients(d - 1, d)
        for modulus in GCD_PRIMES:
            a, b = recurrence_pair_mod(d, modulus)
            require(
                list(a) == [x % modulus for x in exact_a],
                f"recurrence/direct A d,q={d,modulus}",
            )
            require(
                list(b) == [x % modulus for x in exact_b],
                f"recurrence/direct B d,q={d,modulus}",
            )
    return True


def modular_gcd_census():
    counts = []
    for modulus in GCD_PRIMES:
        require(is_prime(modulus), f"modulus not prime q={modulus}")
        checked = 0
        flags = []
        for d in range(GCD_D_MIN, GCD_D_MAX + 1):
            if not composite(d):
                continue
            a, b = recurrence_pair_mod(d, modulus)
            require(
                a.degree() == d - 2 and b.degree() == d - 1,
                f"leading degree lost d,q={d,modulus}",
            )
            checked += 1
            degree = a.gcd(b).degree()
            if degree:
                flags.append((d, degree))
        require(not flags, f"modular common factor q={modulus}: {flags}")
        counts.append((modulus, is_prime(modulus), checked, flags))
    return counts


def rational_hostiles():
    x = sp.symbols("x")
    rows = []
    for d in (4, 6, 8, 9, 10, 12, 15, 16, 18, 20):
        a = sp.Poly.from_list(list(reversed(coefficients(d - 2, d))), gens=x)
        b = sp.Poly.from_list(list(reversed(coefficients(d - 1, d))), gens=x)
        degree = sp.gcd(a, b).degree()
        require(degree == 0, f"rational gcd d={d}")
        divisor_overlap = all(
            set(slopes(lower_hull(list(reversed(a.all_coeffs())), p)))
            & set(slopes(lower_hull(list(reversed(b.all_coeffs())), p)))
            for p, _ in factor(d)
        )
        require(divisor_overlap, f"rational hostile lost overlap d={d}")
        rows.append((d, degree))
    return rows


def main():
    print("THM-3142 PRIME-POWER PREDECESSOR NEWTON SEPARATION -- EXACT CONTROLS")
    (
        composites,
        divisor_places,
        even_single_slope,
        predecessor_places,
        predecessor_separations,
        examples,
    ) = exact_place_census()
    print(
        f"exact_composite_universe=d_4_through_{EXACT_DIVISOR_LIMIT} "
        f"composites={composites} divisor_places={divisor_places} "
        f"all_divisor_places_overlap=True"
    )
    print(f"even_d_p2_single_slope_one_cases={even_single_slope}")
    print(
        f"predecessor_places={predecessor_places} separations={predecessor_separations} "
        "separates_iff_d_minus_1_is_power_of_p=True"
    )
    for key in sorted(examples):
        sa, sb, overlap = examples[key]
        print(f"hostile_d={key[0]} p={key[1]} A_slopes={sa} B_slopes={sb} overlap={overlap}")
    pp = prime_power_predecessor_controls()
    print(f"prime_power_predecessor_exact_cases={len(pp)} first={pp[0]} last={pp[-1]}")
    tested, separating = all_small_places_classification()
    print(
        f"all_prime_places_p_le_2d_d_le_{ALL_PLACE_LIMIT} tested={tested} "
        f"separating={separating} exactly_prime_power_predecessors=True"
    )
    print(f"recurrence_matches_direct_d_4_through_15={recurrence_direct_controls()}")
    for modulus, prime_flag, checked, flags in modular_gcd_census():
        print(
            f"modular_gcd_prime={modulus} isprime={prime_flag} "
            f"composite_d={GCD_D_MIN}..{GCD_D_MAX} "
            f"checked={checked} common_factor_flags={flags}"
        )
    print(f"rational_gcd_hostile_controls={rational_hostiles()}")


if __name__ == "__main__":
    main()
