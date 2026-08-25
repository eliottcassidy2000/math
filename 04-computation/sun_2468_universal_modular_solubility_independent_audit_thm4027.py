#!/usr/bin/env python3
"""Finite gates for a proof of universal local solubility of Sun's form.

The accompanying report supplies the Cauchy--Davenport, Hensel, 2-adic, and
CRT arguments.  This program exhausts precisely the finite odd-prime base
cases left below the analytic cutoff, checks hostile 2-adic levels, and
independently builds the complete profile modulo 99.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import ceil, comb, gcd


N = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)
ODD_PRIMES_BELOW_89 = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_log_p(k: int, p: int) -> int:
    e = 0
    while k >= p:
        k //= p
        e += 1
    return e


def prime_period(p: int, a: int, k: int) -> int:
    return p ** (a + floor_log_p(k, p))


def regular_triangular_image(p: int) -> set[int]:
    # 2w-1 is twice the formal derivative of C(w,2), and is valid also at
    # p=3,5,7.  The w-coordinate has period p modulo every odd p.
    return {comb(w, 2) % p for w in range(p) if (2 * w - 1) % p}


def binomial_image_mod_prime(p: int, k: int) -> set[int]:
    return {comb(x, k) % p for x in range(prime_period(p, 1, k))}


def addsets(left: set[int], right: set[int], modulus: int) -> set[int]:
    return {(a + b) % modulus for a in left for b in right}


def finite_regular_cover(p: int) -> tuple[tuple[int, ...], int]:
    images = (regular_triangular_image(p),) + tuple(
        binomial_image_mod_prime(p, k) for k in (4, 6, 8)
    )
    total = {0}
    for image in images:
        total = addsets(total, image, p)
    return tuple(len(image) for image in images), len(total)


def factor_prime_powers(q: int) -> tuple[tuple[int, int], ...]:
    answer = []
    candidate = 2
    remaining = q
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            exponent = 0
            while remaining % candidate == 0:
                remaining //= candidate
                exponent += 1
            answer.append((candidate, exponent))
        candidate += 1
    if remaining > 1:
        answer.append((remaining, 1))
    return tuple(answer)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def composite_period(q: int, k: int) -> int:
    answer = 1
    for p, a in factor_prime_powers(q):
        answer = lcm(answer, prime_period(p, a, k))
    return answer


def cyclic_convolution(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    q = len(left)
    require(len(right) == q, "convolution modulus")
    answer = [0] * q
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    answer[(i + j) % q] += a * b
    return tuple(answer)


def composite_profile(q: int) -> tuple[tuple[Fraction, ...], tuple[int, ...]]:
    total = (1,) + (0,) * (q - 1)
    periods = []
    universe = 1
    for k in DEGREES:
        orbit = composite_period(q, k)
        hist = [0] * q
        for x in range(orbit):
            hist[comb(x, k) % q] += 1
        total = cyclic_convolution(total, tuple(hist))
        periods.append(orbit)
        universe *= orbit
    require(sum(total) == universe, "profile mass")
    return tuple(Fraction(q * count, universe) for count in total), tuple(periods)


def main() -> None:
    finite_rows = []
    for p in ODD_PRIMES_BELOW_89:
        sizes, covered = finite_regular_cover(p)
        require(covered == p, f"regular base cover failed at p={p}")
        finite_rows.append((p, sizes, covered))

    # The proof for every p>=89 uses
    #   (p-1)/2 + ceil(p/4)+ceil(p/6)+ceil(p/8)-3 >= p.
    # The coarser real lower bound 25p/24-7/2 is already >=p for p>=84.
    require(Fraction(25 * 89, 24) - Fraction(7, 2) >= 89, "CD cutoff")
    cd_at_89 = (89 - 1) // 2 + ceil(89 / 4) + ceil(89 / 6) + ceil(89 / 8) - 3
    require(cd_at_89 >= 89, "integer CD cutoff")

    # Exact hostile controls for the triangular 2-adic induction: over its
    # period 2^(a+1), every residue mod 2^a occurs exactly twice.
    two_adic_rows = []
    for a in range(1, 15):
        q = 2**a
        hist = [0] * q
        for w in range(2 * q):
            hist[comb(w, 2) % q] += 1
        require(set(hist) == {2}, f"2-adic balance failed at a={a}")
        # This is the exact toggle identity used by the induction.
        for w in range(min(2 * q, 257)):
            delta = comb(w + 2 * q, 2) - comb(w, 2)
            require(delta // q % 2 == 1, "2-adic toggle lost oddness")
        two_adic_rows.append((a, q))

    # Independent composite-modulus profile, rather than multiplying stored
    # prime profiles.  It verifies both the target and the transcript's true
    # worst mod-99 class.
    profile99, periods99 = composite_profile(99)
    target99 = N % 99
    minimum99 = min(profile99)
    minimizers99 = tuple(i for i, value in enumerate(profile99) if value == minimum99)
    require(target99 == 53, "target mod 99")
    require(profile99[target99] == Fraction(544, 1089), "target factor mod 99")
    require(minimum99 == Fraction(496, 1089), "minimum factor mod 99")
    require(minimizers99 == (86,), "minimum residue mod 99")
    require(all(value > 0 for value in profile99), "mod-99 image hole")

    payload = "\n".join(
        [f"{p}:{','.join(map(str, sizes))}:{covered}" for p, sizes, covered in finite_rows]
        + [f"2^{a}:{q}" for a, q in two_adic_rows]
        + [f"99:{target99}:{profile99[target99]}:{minimum99}:{minimizers99}"]
    )
    print("SUN_2468_UNIVERSAL_LOCAL_SOLUBILITY_FINITE_GATES")
    print("finite_odd_prime_regular_covers")
    for p, sizes, covered in finite_rows:
        print(f"p={p} image_sizes={sizes} regular_sumset={covered}/{p}")
    print(f"cauchy_davenport_cutoff=89 integer_lower_bound_at_cutoff={cd_at_89}")
    print("cauchy_davenport_scope=ALL_PRIMES_P_GE_89_BY_25P_OVER_24_MINUS_7_OVER_2")
    print(f"two_adic_exact_balance_levels=1..{two_adic_rows[-1][0]}")
    print("two_adic_theorem_scope=ALL_LEVELS_BY_TOGGLE_INDUCTION")
    print(
        f"mod99 periods={periods99} target={target99} target_factor={profile99[target99]} "
        f"minimum={minimum99}@{minimizers99} image=99/99"
    )
    print(f"semantic_sha256={sha256(payload.encode('ascii')).hexdigest()}")
    print("theorem=EVERY_RESIDUE_IS_REPRESENTED_MODULO_EVERY_POSITIVE_MODULUS")
    print("proof_sidecars=CAUCHY_DAVENPORT_PLUS_REGULAR_W_HENSEL_PLUS_2_ADIC_TOGGLE_PLUS_CRT")
    print("global_integer_coverage_inference=NONE")
    print("PASS")


if __name__ == "__main__":
    main()
