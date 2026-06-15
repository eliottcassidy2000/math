#!/usr/bin/env python3
"""
a000568_divisor_profile_codex.py

Exact A000568 via odd-divisor profiles.

Main idea:
  Burnside already reduces tournament counting to odd-part partitions.
  The remaining cross-term

      sum_{i<j} m_i m_j gcd(l_i, l_j)

  is still usually computed pairwise across distinct part sizes.  Replace that
  by the arithmetic identity

      gcd(a,b) = sum_{d|a, d|b} phi(d).

  Then each time we add m copies of an odd part k, its interaction with the
  previously chosen odd parts is:

      m * sum_{d|k} phi(d) * S_d

  where S_d = number of already-chosen parts divisible by d.

  So the whole Burnside exponent can be built incrementally from a divisor
  profile S instead of repeated pairwise gcd scans.
"""

from __future__ import annotations

from math import factorial, gcd, isqrt
from time import perf_counter


KNOWN = {
    0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
}


def phi_upto(n: int) -> list[int]:
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for k in range(p, n + 1, p):
                phi[k] -= phi[k] // p
    return phi


def odd_divisors(k: int) -> list[int]:
    out = []
    r = isqrt(k)
    for d in range(1, r + 1, 2):
        if k % d == 0:
            out.append(d)
            e = k // d
            if e != d and e % 2 == 1:
                out.append(e)
    out.sort()
    return out


def a000568_pairwise(n: int) -> int:
    """Baseline exact implementation with pairwise gcd across distinct sizes."""
    if n <= 1:
        return 1
    fact = [1] * (n + 1)
    for i in range(2, n + 1):
        fact[i] = fact[i - 1] * i
    nfact = fact[n]
    pow2 = [1] * (n * (n - 1) // 2 + 1)
    for i in range(1, len(pow2)):
        pow2[i] = pow2[i - 1] << 1
    gcd_tab = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            gcd_tab[i][j] = gcd(i, j)

    total = 0

    def rec(remaining: int, max_size: int, parts: list[tuple[int, int]]) -> None:
        nonlocal total
        if remaining == 0:
            c = 0
            denom = 1
            for s, m in parts:
                c += m * (s // 2)
                c += m * (m - 1) // 2 * s
                denom *= pow(s, m) * fact[m]
            for i, (s1, m1) in enumerate(parts):
                for s2, m2 in parts[i + 1:]:
                    c += m1 * m2 * gcd_tab[s1][s2]
            total += (nfact // denom) * pow2[c]
            return
        start = min(remaining, max_size)
        if start % 2 == 0:
            start -= 1
        for size in range(start, 0, -2):
            max_mult = remaining // size
            for mult in range(1, max_mult + 1):
                parts.append((size, mult))
                rec(remaining - size * mult, size - 2, parts)
                parts.pop()

    rec(n, n, [])
    return total // nfact


def a000568_divisor_profile(n: int) -> int:
    """Exact divisor-profile engine."""
    if n <= 1:
        return 1

    fact = [1] * (n + 1)
    for i in range(2, n + 1):
        fact[i] = fact[i - 1] * i
    nfact = fact[n]
    max_exp = n * (n - 1) // 2
    pow2 = [1] * (max_exp + 1)
    for i in range(1, max_exp + 1):
        pow2[i] = pow2[i - 1] << 1

    phi = phi_upto(n)
    odds = list(range(1, n + 1, 2))
    divs = {k: odd_divisors(k) for k in odds}
    part_weight = {
        k: sum(phi[d] for d in divs[k]) for k in odds
    }
    # profile[d] = number of already chosen odd cycles with length divisible by d
    profile = [0] * (n + 1)

    total = 0

    def rec(remaining: int, max_size: int, exp_so_far: int, denom_so_far: int) -> None:
        nonlocal total
        if remaining == 0:
            total += (nfact // denom_so_far) * pow2[exp_so_far]
            return
        start = min(remaining, max_size)
        if start % 2 == 0:
            start -= 1
        for size in range(start, 0, -2):
            ds = divs[size]
            base_cross = sum(phi[d] * profile[d] for d in ds)
            max_mult = remaining // size
            kpow = 1
            mfact = 1
            for mult in range(1, max_mult + 1):
                kpow *= size
                mfact *= mult
                # self = m*(k-1)/2 + C(m,2)*k
                self_term = mult * (size // 2) + (mult * (mult - 1) // 2) * size
                # cross = m * sum_{d|k} phi(d) * S_d(previous)
                # plus the extra same-size interaction already accounted for in self_term
                cross_term = mult * base_cross
                for d in ds:
                    profile[d] += mult
                rec(
                    remaining - mult * size,
                    size - 2,
                    exp_so_far + self_term + cross_term,
                    denom_so_far * kpow * mfact,
                )
                for d in ds:
                    profile[d] -= mult

    rec(n, n, 0, 1)
    return total // nfact


def benchmark() -> list[str]:
    lines: list[str] = []
    lines.append("A000568 via divisor profiles")
    lines.append("=" * 72)
    lines.append("Arithmetic delegation:")
    lines.append("  odd-cycle filter from tournament Burnside")
    lines.append("  gcd(a,b) = sum_{d|a,b} phi(d)")
    lines.append("  state = odd-divisor profile S_d instead of pairwise gcd loops")
    lines.append("")
    lines.append("Verification against known values:")
    for n in range(0, 15):
        val = a000568_divisor_profile(n)
        mark = "OK" if val == KNOWN[n] else f"FAIL expected {KNOWN[n]}"
        lines.append(f"  a({n:2d}) = {val} [{mark}]")
    lines.append("")
    lines.append(
        f"{'n':>4} | {'pairwise s':>10} | {'div-prof s':>10} | {'speedup':>8} | {'digits':>8}"
    )
    lines.append("-" * 60)
    for n in [20, 25, 30, 35, 40, 45, 50]:
        t0 = perf_counter()
        pairwise = a000568_pairwise(n)
        t1 = perf_counter() - t0
        t0 = perf_counter()
        prof = a000568_divisor_profile(n)
        t2 = perf_counter() - t0
        assert pairwise == prof
        speedup = t1 / t2 if t2 else 0.0
        lines.append(
            f"{n:4d} | {t1:10.4f} | {t2:10.4f} | {speedup:8.2f} | {len(str(prof)):8d}"
        )
    lines.append("")
    lines.append("Odd-part partition geometry (why Python sees little speedup):")
    lines.append(
        f"{'n':>4} | {'odd parts':>9} | {'avg distinct':>12} | {'avg pair gcd slots':>18} | {'avg odd-div slots':>17}"
    )
    lines.append("-" * 72)
    for n in [20, 30, 40, 50, 60, 80, 100]:
        count = 0
        distinct_total = 0
        pair_slots = 0
        div_slots = 0

        def geom(remaining: int, max_size: int, parts: list[tuple[int, int]]) -> None:
            nonlocal count, distinct_total, pair_slots, div_slots
            if remaining == 0:
                count += 1
                d = len(parts)
                distinct_total += d
                pair_slots += d * (d - 1) // 2
                div_slots += sum(len(odd_divisors(size)) for size, _ in parts)
                return
            start = min(remaining, max_size)
            if start % 2 == 0:
                start -= 1
            for size in range(start, 0, -2):
                max_mult = remaining // size
                for mult in range(1, max_mult + 1):
                    parts.append((size, mult))
                    geom(remaining - mult * size, size - 2, parts)
                    parts.pop()

        geom(n, n, [])
        lines.append(
            f"{n:4d} | {count:9d} | {distinct_total / count:12.3f} | "
            f"{pair_slots / count:18.3f} | {div_slots / count:17.3f}"
        )
    lines.append("")
    lines.append("Number-theory interpretation:")
    lines.append("  Burnside leaves only odd partitions.")
    lines.append("  The tournament-specific combinatorics then enters only through")
    lines.append("  odd divisibility profiles and totients, not through explicit")
    lines.append("  pairwise cycle-size interactions.")
    lines.append("  But the benchmark shows the current Python bottleneck is partition")
    lines.append("  enumeration itself, not gcd accumulation.")
    lines.append("")
    lines.append("Route forward:")
    lines.append("  1. modularize the same divisor-profile recurrence for CRT/GMP engines;")
    lines.append("  2. memoize low-level profile fibers where the remaining odd divisors are sparse;")
    lines.append("  3. search for a closed recurrence on compressed divisor strata;")
    lines.append("  4. do not expect a Python-only speedup unless partition generation is also compressed.")
    return lines


if __name__ == "__main__":
    for line in benchmark():
        print(line)
