#!/usr/bin/env python3
"""
a000568_core_tail_quotient_codex.py

Exact exploration of a second quotient for odd partitions:

    odd partition = (>1 core) + (tail of 1s)

If lambda = mu union 1^r, where mu has odd parts >= 3, let
  m = |mu|          (core mass)
  t = #parts(mu)    (core length)

Then the Burnside exponent splits exactly as

  e(lambda) = e(mu) + C(r,2) + r*t

and z(lambda) = z(mu) * r!.

So the full A000568 sum factors as

  a(n) = sum_{m<=n} sum_t B[m,t] * 2^{C(n-m,2) + (n-m)t} / (n-m)!

where

  B[m,t] = sum_{mu core of mass m, length t} 2^{e(mu)} / z(mu).

This compresses the outer odd-part partition sum to O(n^2) states (m,t).
The remaining hard object is the core kernel B.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import factorial, gcd


KNOWN = {
    0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304, 15: 31021002160355166848,
    16: 63530415842308265100288, 17: 244912778438520759443245824,
    18: 1783398846284777975419600287232, 19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}


def odd_partitions_core(max_n: int):
    """Yield all core partitions with odd parts >=3 and sum <= max_n."""
    yield Counter()

    def rec(remaining: int, max_part: int, current: Counter[int]) -> None:
        start = min(remaining, max_part)
        if start % 2 == 0:
            start -= 1
        for part in range(start, 2, -2):
            max_mult = remaining // part
            for mult in range(1, max_mult + 1):
                current[part] += mult
                yield Counter(current)
                yield from rec(remaining - mult * part, part - 2, current)
                current[part] -= mult
                if current[part] == 0:
                    del current[part]

    yield from rec(max_n, max_n if max_n % 2 else max_n - 1, Counter())


def exponent_from_counter(parts: Counter[int]) -> int:
    lengths = []
    for k, m in parts.items():
        lengths.extend([k] * m)
    e = 0
    for k in lengths:
        e += (k - 1) // 2
    for i, a in enumerate(lengths):
        for b in lengths[i + 1:]:
            e += gcd(a, b)
    return e


def z_from_counter(parts: Counter[int]) -> int:
    z = 1
    for k, m in parts.items():
        z *= (k ** m) * factorial(m)
    return z


def kernel_table(max_n: int) -> dict[tuple[int, int], Fraction]:
    table: dict[tuple[int, int], Fraction] = {}
    seen: set[tuple[tuple[int, int], ...]] = set()
    for core in odd_partitions_core(max_n):
        key = tuple(sorted(core.items()))
        if key in seen:
            continue
        seen.add(key)
        m = sum(k * v for k, v in core.items())
        t = sum(core.values())
        w = Fraction(1 << exponent_from_counter(core), z_from_counter(core))
        table[(m, t)] = table.get((m, t), Fraction(0)) + w
    return table


def a_from_kernel(n: int, kernel: dict[tuple[int, int], Fraction]) -> int:
    total = Fraction(0)
    for (m, t), weight in kernel.items():
        if m > n:
            continue
        r = n - m
        tail = Fraction(1 << ((r * (r - 1) // 2) + r * t), factorial(r))
        total += weight * tail
    assert total.denominator == 1
    return total.numerator


def count_states(max_n: int) -> list[tuple[int, int, int]]:
    # dp[m][t] = count of cores of mass m with t parts (odd parts >= 3)
    dp = [[0] * (max_n + 1) for _ in range(max_n + 1)]
    dp[0][0] = 1
    rows = []
    for part in range(3, max_n + 1, 2):
        for m in range(part, max_n + 1):
            for t in range(1, max_n + 1):
                dp[m][t] += dp[m - part][t - 1]
    # odd partitions of n
    odd = [0] * (max_n + 1)
    odd[0] = 1
    for part in range(1, max_n + 1, 2):
        for m in range(part, max_n + 1):
            odd[m] += odd[m - part]
    active = 0
    for n in [20, 30, 40, 50, 60, 80, 100]:
        active = sum(1 for m in range(n + 1) for t in range(n + 1) if dp[m][t])
        rows.append((n, odd[n], active))
    return rows


def main() -> None:
    max_verify = 20
    kernel = kernel_table(max_verify)

    print("A000568 core-tail quotient")
    print("=" * 72)
    print("Exact factorization:")
    print("  odd partition = (>1 core) + 1-tail")
    print("  a(n) = sum_{m,t} B[m,t] * 2^(C(r,2)+r*t) / r!,  r=n-m")
    print("")
    print("Kernel support through max_verify:")
    support_by_n = []
    for n in range(max_verify + 1):
        support = sum(1 for (m, t) in kernel if m <= n)
        support_by_n.append((n, support))
    for n, support in support_by_n[::5]:
        print(f"  n<={n:2d}: kernel states available = {support}")
    print("")
    print("Verification against known A000568 values:")
    all_ok = True
    for n in range(max_verify + 1):
        val = a_from_kernel(n, kernel)
        ok = val == KNOWN[n]
        all_ok &= ok
        print(f"  a({n:2d}) = {val} [{'OK' if ok else 'FAIL'}]")
    print("")
    print(f"All checks through n={max_verify}: {all_ok}")
    print("")
    print("Outer quotient compression:")
    print(f"{'n':>4} | {'odd partitions':>14} | {'(m,t) states':>12} | {'compression':>11}")
    print("-" * 56)
    for n, odd_count, active in count_states(100):
        ratio = odd_count / active if active else 0
        print(f"{n:4d} | {odd_count:14d} | {active:12d} | {ratio:11.2f}x")
    print("")
    print("Interpretation:")
    print("  1. This is an exact second quotient on odd partitions for the outer sum.")
    print("  2. The remaining hard object is the core kernel B[m,t].")
    print("  3. Compression is already strong at the outer layer:")
    print("     for n=100, 444793 odd partitions collapse to 834 active (m,t) states.")
    print("  4. This does not by itself solve the kernel-generation problem,")
    print("     but it isolates where the true complexity lives.")


if __name__ == "__main__":
    main()
