#!/usr/bin/env python3
"""Exact companion for THM-2455.

The proof in THM-2455 is uniform.  This dependency-free companion supplies
independent finite audits of its valuation-box criterion, prime-ray tails,
additive analogue, induced-cover criterion, and twin-centre specialization.
All correctness checks use ``require`` so ``python3`` and ``python3 -O`` run
the same test program.
"""

from fractions import Fraction
from itertools import combinations, product
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime_sieve(limit):
    table = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        table[0] = 0
    if limit >= 1:
        table[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if table[p]:
            start = p * p
            table[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return table


def factorization(n):
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            exponent = 0
            while n % d == 0:
                n //= d
                exponent += 1
            factors.append((d, exponent))
        d = 3 if d == 2 else d + 2
    if n > 1:
        factors.append((n, 1))
    return factors


def least_prime_factor(n):
    require(n >= 2, "least prime factor called below two")
    if n % 2 == 0:
        return 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return d
        d += 2
    return n


def omega(n):
    return sum(exponent for _, exponent in factorization(n))


def direct_factor_pairs(z, strict):
    pairs = []
    for a in range(2, isqrt(z) + 1):
        if z % a == 0:
            b = z // a
            if not strict or a < b:
                pairs.append((a, b))
    return tuple(pairs)


def valuation_box_pairs(z, strict):
    """Generate factor pairs from the exponent box, not by trial division."""
    fac = factorization(z)
    pairs = set()
    for beta in product(*(range(exponent + 1) for _, exponent in fac)):
        a = 1
        for (p, _), exponent in zip(fac, beta):
            a *= p**exponent
        if a == 1 or a == z:
            continue
        b = z // a
        if a > b or (strict and a == b):
            continue
        pairs.add((a, b))
    return tuple(sorted(pairs))


def artificial_multiplicative_atom(holes, z, pairs):
    return z not in holes and bool(pairs) and all(
        a in holes or b in holes for a, b in pairs
    )


def nonunit_divisors(n):
    divisors = []
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            divisors.append(d)
            if d * d != n:
                divisors.append(n // d)
    divisors.append(n)
    return tuple(sorted(set(divisors)))


def cage(holes):
    return frozenset(
        m for m in holes if all(d in holes for d in nonunit_divisors(m))
    )


def additive_pairs(z, strict):
    pairs = []
    for a in range(1, z // 2 + 1):
        b = z - a
        if not strict or a < b:
            pairs.append((a, b))
    return tuple(pairs)


def artificial_additive_atom(holes, z, pairs):
    return z not in holes and bool(pairs) and all(
        a in holes or b in holes for a, b in pairs
    )


def proper_divisors(n):
    return tuple(d for d in range(2, n) if n % d == 0)


def induced_cover_formula(holes, x, z):
    if x in holes or z in holes or x >= z or z % x:
        return False
    q = z // x
    return all(x * d in holes for d in proper_divisors(q))


def induced_cover_direct(holes, x, z):
    if x in holes or z in holes or x >= z or z % x:
        return False
    for y in range(x + 1, z):
        if y not in holes and y % x == 0 and z % y == 0:
            return False
    return True


def bitmask_set(mask, universe):
    return frozenset(
        value for bit, value in enumerate(universe) if mask & (1 << bit)
    )


def main():
    # One sieve supports the tail phase and the complete twin-centre census.
    is_prime = prime_sieve(1_000_001)

    # ------------------------------------------------------------------
    # 1. Direct factor fibres versus the exponent-box criterion.
    # ------------------------------------------------------------------
    direct_weak_500 = {
        z: direct_factor_pairs(z, False) for z in range(2, 501)
    }
    direct_strict_500 = {
        z: direct_factor_pairs(z, True) for z in range(2, 501)
    }
    box_weak_500 = {
        z: valuation_box_pairs(z, False) for z in range(2, 501)
    }
    box_strict_500 = {
        z: valuation_box_pairs(z, True) for z in range(2, 501)
    }
    require(
        any(
            direct_weak_500[z] != direct_strict_500[z]
            for z in range(2, 501)
        ),
        "weak/strict diagonal control was not activated",
    )

    valuation_checks = 0
    hole_universe = tuple(range(2, 13))
    for mask in range(1 << len(hole_universe)):
        holes = bitmask_set(mask, hole_universe)
        for z in range(2, 501):
            weak_direct = artificial_multiplicative_atom(
                holes, z, direct_weak_500[z]
            )
            weak_box = artificial_multiplicative_atom(
                holes, z, box_weak_500[z]
            )
            require(weak_direct == weak_box, ("weak box mismatch", holes, z))
            valuation_checks += 1

            strict_direct = artificial_multiplicative_atom(
                holes, z, direct_strict_500[z]
            )
            strict_box = artificial_multiplicative_atom(
                holes, z, box_strict_500[z]
            )
            require(
                strict_direct == strict_box,
                ("strict box mismatch", holes, z),
            )
            valuation_checks += 1
    require(valuation_checks == 2_043_904, "valuation audit count drift")

    # ------------------------------------------------------------------
    # 2. Prime cages, disjoint rays, and the M^3/M^4 remote bounds.
    #
    # The 179-hole atlas contains every singleton and pair in 2..12 and
    # 113 lexicographically first triples.  Each of its 15,000 targets
    # receives one weak and one strict tail check.
    # ------------------------------------------------------------------
    tail_holes = []
    for size in (1, 2):
        tail_holes.extend(frozenset(row) for row in combinations(hole_universe, size))
    tail_holes.extend(
        frozenset(row)
        for row in list(combinations(hole_universe, 3))[:113]
    )
    require(len(tail_holes) == 179, "tail atlas size drift")
    require(len(set(tail_holes)) == 179, "duplicate tail hole set")
    require(frozenset({2, 3, 6}) in tail_holes, "missing ray hostile")
    require(frozenset({8, 12}) in tail_holes, "missing diagonal hostile")

    tail_limit = 15_001
    weak_tail_pairs = {
        z: direct_factor_pairs(z, False) for z in range(2, tail_limit + 1)
    }
    strict_tail_pairs = {
        z: direct_factor_pairs(z, True) for z in range(2, tail_limit + 1)
    }
    tail_checks = 0
    rational_prefix_checks = 0

    for holes in tail_holes:
        M = max(holes)
        c_h = cage(holes)
        has_prime_hole = any(is_prime[h] for h in holes)
        require(
            bool(c_h) == has_prime_hole,
            ("cage/prime equivalence mismatch", holes, c_h),
        )

        ray_owner = {}
        for m in c_h:
            for q in range(M + 1, tail_limit // m + 1):
                if not is_prime[q]:
                    continue
                z = m * q
                require(z not in ray_owner, ("ray collision", holes, z))
                ray_owner[z] = (m, q)
        ray_targets = frozenset(ray_owner)

        if not has_prime_hole:
            require(not c_h and not ray_targets, "prime-free cage survived")
            prime_free_bound = max(
                h * least_prime_factor(h) for h in holes
            )
        else:
            prime_free_bound = None

        for z in range(2, tail_limit + 1):
            weak = artificial_multiplicative_atom(
                holes, z, weak_tail_pairs[z]
            )
            strict = artificial_multiplicative_atom(
                holes, z, strict_tail_pairs[z]
            )
            on_ray = z in ray_targets

            require(
                not on_ray or weak,
                ("prime ray missed weak atom", holes, z, ray_owner.get(z)),
            )
            require(
                on_ray or not weak or z <= M**3,
                ("weak finite remainder exceeded M^3", holes, z),
            )
            if prime_free_bound is not None and weak:
                require(
                    z <= prime_free_bound and omega(z) >= 3,
                    ("prime-free weak bound", holes, z, prime_free_bound),
                )
            tail_checks += 1

            require(
                not on_ray or strict,
                ("prime ray missed strict atom", holes, z, ray_owner.get(z)),
            )
            require(
                on_ray or not strict or z <= M**4,
                ("strict finite remainder exceeded M^4", holes, z),
            )
            if prime_free_bound is not None and strict:
                require(
                    z <= prime_free_bound and omega(z) >= 3,
                    ("prime-free strict bound", holes, z, prime_free_bound),
                )
            tail_checks += 1

        # An exact rational prefix of the ray Dirichlet product.
        prefix_primes = tuple(
            q for q in range(M + 1, 98) if is_prime[q]
        )
        direct_mass = sum(
            (Fraction(1, m * q) for m in c_h for q in prefix_primes),
            Fraction(0),
        )
        product_mass = sum(
            (Fraction(1, m) for m in c_h), Fraction(0)
        ) * sum((Fraction(1, q) for q in prefix_primes), Fraction(0))
        require(direct_mass == product_mass, ("ray mass mismatch", holes))
        require(
            len({m * q for m in c_h for q in prefix_primes})
            == len(c_h) * len(prefix_primes),
            ("prefix ray collision", holes),
        )
        rational_prefix_checks += 1

    require(tail_checks == 5_370_000, "tail audit count drift")
    require(rational_prefix_checks == 179, "ray prefix count drift")

    # Sharp multiplicative controls.
    for p in (2, 3, 5, 7, 11):
        prime_free = frozenset({p * p})
        z = p**3
        bound = max(h * least_prime_factor(h) for h in prime_free)
        require(z == bound, ("prime-free sharp bound drift", p))
        require(
            artificial_multiplicative_atom(
                prime_free, z, direct_factor_pairs(z, False)
            ),
            ("prime-free weak sharpness failed", p),
        )
        require(
            artificial_multiplicative_atom(
                prime_free, z, direct_factor_pairs(z, True)
            ),
            ("prime-free strict sharpness failed", p),
        )

        prime_hole = frozenset({p})
        require(
            artificial_multiplicative_atom(
                prime_hole, p**3, direct_factor_pairs(p**3, False)
            )
            and p**3 not in {
                m * q
                for m in cage(prime_hole)
                for q in range(p + 1, p**3 + 1)
                if q < len(is_prime) and is_prime[q]
            },
            ("weak M^3 sharpness failed", p),
        )
        require(
            artificial_multiplicative_atom(
                prime_hole, p**4, direct_factor_pairs(p**4, True)
            )
            and not artificial_multiplicative_atom(
                prime_hole, p**4, direct_factor_pairs(p**4, False)
            ),
            ("strict M^4 sharpness/diagonal failed", p),
        )

    ray_hostile = frozenset({2, 3, 6})
    require(6 in cage(ray_hostile), "6 missing from hostile cage")
    for q in (7, 11, 13, 17, 19):
        z = 6 * q
        require(
            omega(z) == 3
            and artificial_multiplicative_atom(
                ray_hostile, z, direct_factor_pairs(z, False)
            )
            and artificial_multiplicative_atom(
                ray_hostile, z, direct_factor_pairs(z, True)
            ),
            ("non-semiprime ray hostile failed", q),
        )

    require(
        not artificial_multiplicative_atom(
            frozenset({8}), 16, direct_factor_pairs(16, False)
        )
        and artificial_multiplicative_atom(
            frozenset({8}), 16, direct_factor_pairs(16, True)
        ),
        "H={8}, z=16 diagonal hostile failed",
    )
    require(
        not artificial_multiplicative_atom(
            frozenset({8, 12}), 24, direct_factor_pairs(24, False)
        )
        and not artificial_multiplicative_atom(
            frozenset({8, 12}), 24, direct_factor_pairs(24, True)
        ),
        "H={8,12}, z=24 surviving-pair hostile failed",
    )

    # ------------------------------------------------------------------
    # 3. Additive finite-hole bounds.
    # ------------------------------------------------------------------
    additive_weak_pairs = {
        z: additive_pairs(z, False) for z in range(2, 160)
    }
    additive_strict_pairs = {
        z: additive_pairs(z, True) for z in range(2, 160)
    }
    additive_checks = 0
    additive_universe = tuple(range(1, 9))
    for mask in range(1, 1 << len(additive_universe)):
        holes = bitmask_set(mask, additive_universe)
        L = max(holes)
        for z in range(2, 160):
            weak = artificial_additive_atom(
                holes, z, additive_weak_pairs[z]
            )
            strict = artificial_additive_atom(
                holes, z, additive_strict_pairs[z]
            )
            require(
                not weak or z <= 2 * L + 1,
                ("additive weak tail", holes, z),
            )
            require(
                not strict or z <= 2 * L + 2,
                ("additive strict tail", holes, z),
            )
            additive_checks += 1
    require(additive_checks == 40_290, "additive audit count drift")

    for L in range(1, 9):
        interval_holes = frozenset(range(1, L + 1))
        require(
            artificial_additive_atom(
                interval_holes,
                2 * L + 1,
                additive_pairs(2 * L + 1, False),
            ),
            ("additive weak sharpness", L),
        )
        require(
            artificial_additive_atom(
                interval_holes,
                2 * L + 2,
                additive_pairs(2 * L + 2, True),
            ),
            ("additive strict sharpness", L),
        )

    # ------------------------------------------------------------------
    # 4. Induced divisibility covers.
    # ------------------------------------------------------------------
    comparable_pairs = tuple(
        (x, z)
        for z in range(2, 160)
        for x in range(2, z)
        if z % x == 0
    )
    require(len(comparable_pairs) == 513, "cover carrier count drift")
    cover_checks = 0
    cover_universe = tuple(range(2, 9))
    for mask in range(1, 1 << len(cover_universe)):
        holes = bitmask_set(mask, cover_universe)
        M = max(holes)
        for x, z in comparable_pairs:
            direct = induced_cover_direct(holes, x, z)
            formula = induced_cover_formula(holes, x, z)
            require(direct == formula, ("cover formula mismatch", holes, x, z))
            if formula and not is_prime[z // x]:
                require(
                    x * z <= M * M and z <= M * M // 2,
                    ("composite cover bound", holes, x, z),
                )
            cover_checks += 1

    cover_controls = (
        (frozenset({6}), 2, 18, True),
        (frozenset({6}), 2, 12, False),
        (frozenset({6}), 3, 18, False),
        (frozenset({6}), 2, 8, False),
        (frozenset({4, 6}), 2, 8, True),
        (frozenset({4, 6}), 2, 12, True),
        (frozenset({4, 6}), 3, 12, True),
        (frozenset({4, 6}), 2, 18, True),
        (frozenset({4, 6}), 3, 18, False),
        (frozenset({4, 6}), 2, 24, False),
        (frozenset({4, 6}), 3, 24, False),
        (frozenset({2}), 3, 12, False),
        (frozenset({8}), 4, 16, True),
        (frozenset({8}), 2, 16, False),
        (frozenset({8, 12}), 4, 24, True),
        (frozenset({8, 12}), 2, 24, False),
        (frozenset({8, 12}), 3, 24, False),
        (frozenset({4}), 2, 8, True),
        (frozenset({4}), 2, 12, False),
        (frozenset({9}), 3, 27, True),
        (frozenset({6, 9}), 3, 27, True),
        (frozenset({10, 15}), 5, 30, True),
    )
    require(len(cover_controls) == 22, "cover control count drift")
    for holes, x, z, expected in cover_controls:
        direct = induced_cover_direct(holes, x, z)
        formula = induced_cover_formula(holes, x, z)
        require(
            direct == formula == expected,
            ("explicit cover hostile failed", holes, x, z, expected),
        )
        cover_checks += 1
    require(cover_checks == 65_173, "cover audit count drift")
    require(
        induced_cover_formula(frozenset({6}), 2, 18)
        and 2 * 18 == 6**2
        and 18 == 6**2 // 2,
        "sharp induced-cover bound failed",
    )

    # ------------------------------------------------------------------
    # 5. Twin-centre specialization and THM-2433 recovery.
    # ------------------------------------------------------------------
    twin_centres = tuple(
        c
        for c in range(2, 1_000_001)
        if is_prime[c - 1] and is_prime[c + 1]
    )
    require(len(twin_centres) == 8_169, "twin-centre census drift")

    startup_holes = frozenset({4, 6})
    startup_weak = frozenset(
        z
        for z in range(2, 501)
        if artificial_multiplicative_atom(
            startup_holes, z, direct_factor_pairs(z, False)
        )
    )
    startup_strict = frozenset(
        z
        for z in range(2, 501)
        if artificial_multiplicative_atom(
            startup_holes, z, direct_factor_pairs(z, True)
        )
    )
    require(
        startup_weak == startup_strict == frozenset({8, 12}),
        ("THM-2433 atom recovery", startup_weak, startup_strict),
    )

    artificial_twin_centres = []
    for c in twin_centres:
        weak = artificial_multiplicative_atom(
            startup_holes, c, direct_factor_pairs(c, False)
        )
        strict = artificial_multiplicative_atom(
            startup_holes, c, direct_factor_pairs(c, True)
        )
        require(weak == strict, ("twin weak/strict drift", c))
        if weak:
            require(
                c // 2 in startup_holes
                and c // 3 in startup_holes
                and c <= 2 * max(startup_holes),
                ("twin-centre necessary conditions", c),
            )
            artificial_twin_centres.append(c)
    require(
        artificial_twin_centres == [12],
        ("startup twin-centre classification", artificial_twin_centres),
    )

    print("THM-2455 finite-hole relative-atom exact audit")
    print(f"valuation-box weak/strict comparisons: {valuation_checks}")
    print(f"prime-ray and M^3/M^4 tail checks: {tail_checks}")
    print(f"exact rational ray-prefix checks: {rational_prefix_checks}")
    print(f"additive artificial-atom checks: {additive_checks}")
    print(f"induced-cover comparisons and hostiles: {cover_checks}")
    print(
        "twin centres through 1,000,000: "
        f"{len(twin_centres)}; artificial for H={{4,6}}: "
        f"{artificial_twin_centres}"
    )
    print("sharp controls: p^3, p^4, 2L+1, 2L+2, and (2,18) PASS")
    print("hostiles: H={2,3,6}; H={8},z=16; H={8,12},z=24 PASS")
    print("PASS")


if __name__ == "__main__":
    main()
