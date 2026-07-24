#!/usr/bin/env python3
"""Exact companion for THM-2159 (opus-2026-07-24 puzzle atlas).

The load-bearing checks use Python integers only:

* quotient-ring coefficient recurrence modulo y^n-a;
* lossless evaluation at X=a^K;
* equality with the least modular remainders;
* the leading-coefficient quotient; and
* Shunia's final natural-number division.

The universe is every admissible (a,n) with 3 <= a <= 12.  This is a
finite hostile companion, not the proof of the universal inequalities.
"""

from __future__ import annotations

import sys

sys.set_int_max_str_digits(1_000_000)


def nth_root_floor(a: int, n: int) -> int:
    lo, hi = 0, a + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**n <= a:
            lo = mid
        else:
            hi = mid
    return lo


def admissible(a: int, n: int) -> bool:
    if a <= 2 or n <= 1 or n > a.bit_length():
        return False
    d = nth_root_floor(a, n)
    return d**n != a


def next_coefficients(coeff: list[int], a: int) -> list[int]:
    """Coefficients of (1+y)F modulo y^n-a."""
    out = [0] * len(coeff)
    out[0] = coeff[0] + a * coeff[-1]
    for r in range(1, len(coeff)):
        out[r] = coeff[r] + coeff[r - 1]
    return out


def evaluate(coeff: list[int], x: int) -> int:
    value = 0
    for c in reversed(coeff):
        value = value * x + c
    return value


def check_pair(a: int, n: int) -> dict[str, int | bool]:
    k = 2 * a**n
    x = a**k
    modulus = x**n - a

    coeff = [1] + [0] * (n - 1)
    for _ in range(k):
        coeff = next_coefficients(coeff, a)
    coeff_k = coeff
    coeff_k1 = next_coefficients(coeff_k, a)

    f_k = evaluate(coeff_k, x)
    f_k1 = evaluate(coeff_k1, x)
    r_k = pow(x + 1, k, modulus)
    r_k1 = pow(x + 1, k + 1, modulus)
    root = nth_root_floor(a, n)

    checks = {
        "coeff_sum_below_base_k": sum(coeff_k) < x,
        "coeff_sum_below_base_k1": sum(coeff_k1) < x,
        "least_remainder_k": f_k == r_k and 0 <= f_k < modulus,
        "least_remainder_k1": f_k1 == r_k1 and 0 <= f_k1 < modulus,
        "leading_ratio": coeff_k[-2] // coeff_k[-1] == root,
        "final_formula": r_k1 // r_k - 1 == root,
    }
    if not all(checks.values()):
        raise AssertionError((a, n, checks))
    return {
        "a": a,
        "n": n,
        "K": k,
        "root": root,
        "X_digits": len(str(x)),
        **checks,
    }


def main() -> None:
    rows = []
    for a in range(3, 13):
        for n in range(2, a.bit_length() + 1):
            if admissible(a, n):
                rows.append(check_pair(a, n))

    print("THM-2159 SHUNIA FINITE-ROOT SPECTRAL EXTRACTION -- EXACT COMPANION")
    print(f"universe: all admissible (a,n), 3<=a<=12; rows={len(rows)}")
    print("arithmetic: Python integers only")
    for row in rows:
        print(
            f"a={row['a']:2d} n={row['n']} K={row['K']:6d} "
            f"root={row['root']} Xdigits={row['X_digits']:6d} PASS"
        )
    print("ALL EXACT CHECKS PASS")


if __name__ == "__main__":
    main()
