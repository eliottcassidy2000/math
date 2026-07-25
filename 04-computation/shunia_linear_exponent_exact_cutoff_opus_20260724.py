#!/usr/bin/env python3
"""Exact finite-cutoff companion for THM-2165.

Universe: every admissible pair (a,n) with 3 <= a <= 500.
Arithmetic: Python integers only.

For K=a+1 the script constructs (1+y)^K modulo y^n-a coefficient by
coefficient.  It verifies the leading-coordinate quotient, lossless
Kronecker evaluation at X=a^K, the least-remainder inequalities, and the
final natural-number quotient.  Selected boundary and hostile rows are also
checked against Python's independent modular exponentiation path.
"""

from __future__ import annotations

import hashlib
import sys
from collections import Counter

sys.set_int_max_str_digits(1_000_000)


MAX_A = 500
MODULAR_CONTROL_A = {
    3,
    4,
    5,
    8,
    14,
    30,
    63,
    64,
    127,
    128,
    255,
    256,
    499,
    500,
}
PREDECESSOR_FAILURES = ((8, 4), (14, 4))


def nth_root_floor(a: int, n: int) -> int:
    lo, hi = 0, a + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**n <= a:
            lo = mid
        else:
            hi = mid
    return lo


def next_coefficients(coeff: list[int], a: int) -> list[int]:
    """Coefficients of (1+y)F modulo y^n-a."""
    return [coeff[0] + a * coeff[-1]] + [
        coeff[r] + coeff[r - 1] for r in range(1, len(coeff))
    ]


def evaluate(coeff: list[int], x: int) -> int:
    value = 0
    for c in reversed(coeff):
        value = value * x + c
    return value


def require(condition: bool, context: object) -> None:
    if not condition:
        raise AssertionError(context)


def main() -> None:
    counts: Counter[int] = Counter()
    digest = hashlib.sha256()
    rows = 0
    modular_controls = 0
    coefficient_carry_rows: list[tuple[int, int, str]] = []

    for a in range(3, MAX_A + 1):
        for n in range(2, a.bit_length() + 1):
            root = nth_root_floor(a, n)
            if root**n == a:
                continue

            k = a + 1
            coeff = [1] + [0] * (n - 1)
            for _ in range(k):
                coeff = next_coefficients(coeff, a)
            coeff_next = next_coefficients(coeff, a)

            top = coeff[-1]
            top_next = coeff_next[-1]
            require((root + 1) * top < top_next, (a, n, "leading lower"))
            require(top_next < (root + 2) * top, (a, n, "leading upper"))

            x = a**k
            if sum(coeff) >= x:
                coefficient_carry_rows.append((a, n, "K"))
            if sum(coeff_next) >= x:
                coefficient_carry_rows.append((a, n, "K+1"))

            f = evaluate(coeff, x)
            f_next = evaluate(coeff_next, x)
            modulus = x**n - a
            require(0 < f < modulus, (a, n, "least remainder K"))
            require(0 < f_next < modulus, (a, n, "least remainder K+1"))
            require(f_next // f - 1 == root, (a, n, "final quotient"))

            if a in MODULAR_CONTROL_A:
                require(
                    pow(x + 1, k, modulus) == f,
                    (a, n, "modular control K"),
                )
                require(
                    pow(x + 1, k + 1, modulus) == f_next,
                    (a, n, "modular control K+1"),
                )
                modular_controls += 1

            counts[n] += 1
            rows += 1
            digest.update(
                f"{a},{n},{root},{top_next // top},{f_next // f}\n".encode()
            )

    print("THM-2165 SHUNIA LINEAR-EXPONENT EXACT CUTOFF")
    print(f"universe: every admissible (a,n), 3<=a<={MAX_A}")
    print(f"rows: {rows}")
    print(
        "counts by n: "
        + " ".join(f"n={n}:{counts[n]}" for n in sorted(counts))
    )
    print(f"independent modular-exponentiation controls: {modular_controls}")
    print(f"coefficient-sum >= X rows (allowed): {coefficient_carry_rows}")
    predecessor_rows = []
    for a, n in PREDECESSOR_FAILURES:
        k = a
        x = a**k
        modulus = x**n - a
        extracted = pow(x + 1, k + 1, modulus) // pow(x + 1, k, modulus) - 1
        root = nth_root_floor(a, n)
        require(extracted != root, (a, n, "K=a should fail"))
        predecessor_rows.append((a, n, root, extracted))
    print(f"K=a exact failures (root, extracted): {predecessor_rows}")
    print(f"row digest sha256: {digest.hexdigest()}")
    print("ALL EXACT INTEGER CHECKS PASS")


if __name__ == "__main__":
    main()
