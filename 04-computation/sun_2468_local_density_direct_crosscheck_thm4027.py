#!/usr/bin/env python3
"""Independent direct-tuple checks for selected Sun 2-4-6-8 factors.

Unlike ``local_density_audit.py``, this program performs no convolution.  It
iterates the full Cartesian product of the four certified binomial periods and
counts only the target congruence.  The larger moduli 25 and 49 are retained as
hostile controls (about 9.8M and 40.4M tuples respectively).
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import comb


N = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_log_p(k: int, p: int) -> int:
    e = 0
    while k >= p:
        k //= p
        e += 1
    return e


def values(p: int, a: int, k: int) -> tuple[int, ...]:
    q = p**a
    orbit = p ** (a + floor_log_p(k, p))
    return tuple(comb(x, k) % q for x in range(orbit))


def direct_target_count(p: int, a: int) -> tuple[int, int, Fraction]:
    q = p**a
    target = N % q
    v2, v4, v6, v8 = (values(p, a, k) for k in DEGREES)
    count = 0
    # Deliberately direct: no residue histograms and no pair-sum tables.
    for a2 in v2:
        for a4 in v4:
            for a6 in v6:
                partial = (a2 + a4 + a6) % q
                for a8 in v8:
                    count += (partial + a8) % q == target
    universe = len(v2) * len(v4) * len(v6) * len(v8)
    return count, universe, Fraction(q * count, universe)


def main() -> None:
    cases = (
        (3, 1), (3, 2),
        (5, 1), (5, 2),
        (7, 1), (7, 2),
        (11, 1), (13, 1), (17, 1), (19, 1), (23, 1),
        (29, 1), (31, 1), (43, 1),
    )
    expected = {
        (3, 1): Fraction(22, 27),
        (3, 2): Fraction(68, 81),
        (5, 1): Fraction(113, 125),
        (5, 2): Fraction(566, 625),
        (7, 1): Fraction(311, 343),
        (7, 2): Fraction(310, 343),
        (11, 1): Fraction(72, 121),
        (13, 1): Fraction(154, 169),
        (17, 1): Fraction(240, 289),
        (19, 1): Fraction(316, 361),
        (23, 1): Fraction(472, 529),
        (29, 1): Fraction(832, 841),
        (31, 1): Fraction(942, 961),
        (43, 1): Fraction(1880, 1849),
    }
    rows = []
    for p, a in cases:
        count, universe, local_factor = direct_target_count(p, a)
        require(local_factor == expected[(p, a)], f"factor drift at {p}^{a}")
        rows.append((p, a, N % (p**a), count, universe, local_factor))

    payload = "\n".join(
        f"{p},{a},{target},{count},{universe},{factor.numerator}/{factor.denominator}"
        for p, a, target, count, universe, factor in rows
    )
    print("SUN_2468_DIRECT_CARTESIAN_CROSSCHECK")
    print("aggregation=DIRECT_FOUR_NESTED_LOOPS_NO_CONVOLUTION")
    for p, a, target, count, universe, local_factor in rows:
        print(
            f"p={p} a={a} q={p**a} target={target} count={count} "
            f"universe={universe} factor={local_factor}"
        )
    print(f"semantic_sha256={sha256(payload.encode('ascii')).hexdigest()}")
    print("PASS")


if __name__ == "__main__":
    main()
