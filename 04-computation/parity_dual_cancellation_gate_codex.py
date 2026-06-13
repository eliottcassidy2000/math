#!/usr/bin/env python3
"""
parity_dual_cancellation_gate_codex.py

codex-2026-06-12-P5

Scout the centered parity pattern behind two repo threads:

1. triangular power-anchor balancing:
     F_p(c,n) = c^p + sum_{j=1}^n ((c-j)^p - (c+j)^p),
   where only odd Faulhaber moments survive;
2. tournament forward-edge statistics:
     fwd(sigma) + fwd(sigma^rev) = n-1,
   where all odd centered moments vanish.

This does not prove a new theorem. It checks the parity-dual bridge and
packages the exact finite evidence for a fresh hypothesis/reflection packet.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import permutations
from math import comb


getcontext().prec = 90


def power_sum(n: int, r: int) -> int:
    return sum(k**r for k in range(1, n + 1))


def odd_support(p: int) -> list[int]:
    support = []
    for r in range(p + 1):
        coeff = ((-1) ** r - 1) * comb(p, r)
        if coeff != 0:
            support.append(r)
    return support


def balance_direct(p: int, n: int, c: int) -> int:
    total = c**p
    for j in range(1, n + 1):
        total += (c - j) ** p - (c + j) ** p
    return total


def balance_odd_expansion(p: int, n: int, c: int) -> int:
    total = c**p
    for r in range(1, p + 1, 2):
        total -= 2 * comb(p, r) * (c ** (p - r)) * power_sum(n, r)
    return total


def alpha_formula(p: int) -> Fraction:
    return Fraction((p - 1) * (p - 2), 12 * p)


def beta_formula(p: int) -> Fraction:
    num = -(p - 1) * (p - 2) * (2 * p * p - 4 * p - 1)
    den = 180 * p**3
    return Fraction(num, den)


def approx_anchor(p: int, n: int) -> Decimal:
    u = Decimal(n * (n + 1))
    alpha = Decimal(alpha_formula(p).numerator) / Decimal(alpha_formula(p).denominator)
    beta = Decimal(beta_formula(p).numerator) / Decimal(beta_formula(p).denominator)
    return (
        Decimal(p * n * n + (p - 1) * n)
        + alpha
        + beta / u
    )


def anchor_balance_decimal(p: int, n: int, a: Decimal) -> Decimal:
    total = Decimal(0)
    for j in range(n + 1):
        total += (a + Decimal(j)) ** p
    for j in range(1, n + 1):
        total -= (a + Decimal(n + j)) ** p
    return total


def real_anchor(p: int, n: int) -> Decimal:
    guess = approx_anchor(p, n)
    lo = guess - Decimal(2)
    hi = guess + Decimal(2)
    f_lo = anchor_balance_decimal(p, n, lo)
    f_hi = anchor_balance_decimal(p, n, hi)

    while f_lo > 0:
        lo -= Decimal(2)
        f_lo = anchor_balance_decimal(p, n, lo)
    while f_hi < 0:
        hi += Decimal(2)
        f_hi = anchor_balance_decimal(p, n, hi)

    for _ in range(140):
        mid = (lo + hi) / 2
        f_mid = anchor_balance_decimal(p, n, mid)
        if f_mid >= 0:
            hi = mid
        else:
            lo = mid
    return hi


def tournament_from_bits(n: int, bits: int) -> list[list[int]]:
    adj = [[0] * n for _ in range(n)]
    bit_index = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> bit_index) & 1:
                adj[i][j] = 1
                adj[j][i] = 0
            else:
                adj[i][j] = 0
                adj[j][i] = 1
            bit_index += 1
    return adj


def forward_counts(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    out = []
    for sigma in permutations(range(n)):
        count = 0
        for i in range(n - 1):
            if adj[sigma[i]][sigma[i + 1]]:
                count += 1
        out.append(count)
    return out


def centered_moment(values: list[int], center: Fraction, r: int) -> Fraction:
    total = Fraction(0, 1)
    for value in values:
        total += (Fraction(value, 1) - center) ** r
    return total / len(values)


def format_fraction(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def main() -> None:
    print("PARITY-DUAL CANCELLATION GATE SCOUT")
    print("codex-2026-06-12-P5")
    print()

    print("1. Anchor-side odd support")
    for p in range(1, 11):
        print(f"  p={p:2d} surviving_j_powers={odd_support(p)}")
    print("  claim: only odd powers of j survive in (c-j)^p-(c+j)^p.")
    print()

    print("2. Direct midpoint balance vs odd-moment expansion")
    samples = [(3, 4, 61), (4, 5, 123), (5, 6, 217), (8, 4, 290)]
    for p, n, c in samples:
        direct = balance_direct(p, n, c)
        odd_form = balance_odd_expansion(p, n, c)
        status = "PASS" if direct == odd_form else "FAIL"
        print(
            f"  p={p:2d} n={n:2d} c={c:3d} "
            f"direct={direct} odd_form={odd_form} [{status}]"
        )
    print()

    print("3. Real-anchor asymptotics")
    print("  columns: root-minus-approx, scaled by n^4, and scaled by u^2=(n(n+1))^2")
    for p in range(3, 9):
        print(f"  p={p}")
        for n in (25, 50, 100):
            root = real_anchor(p, n)
            approx = approx_anchor(p, n)
            err = root - approx
            scaled = err * Decimal(n**4)
            u_scaled = err * Decimal((n * (n + 1)) ** 2)
            print(
                "    "
                f"n={n:3d} "
                f"err={err:.18E} "
                f"scaled_n4={scaled:.18E} "
                f"scaled_u2={u_scaled:.18E}"
            )
    print()

    print("4. Tournament reversal side (exhaustive n=5)")
    n = 5
    mu = Fraction(n - 1, 2)
    odd_failures = 0
    m2_values: set[Fraction] = set()
    m4_values: set[Fraction] = set()
    sample_records = []
    for bits in range(1 << comb(n, 2)):
        counts = forward_counts(tournament_from_bits(n, bits))
        m1 = centered_moment(counts, mu, 1)
        m3 = centered_moment(counts, mu, 3)
        m5 = centered_moment(counts, mu, 5)
        if any(moment != 0 for moment in (m1, m3, m5)):
            odd_failures += 1
        m2 = centered_moment(counts, mu, 2)
        m4 = centered_moment(counts, mu, 4)
        m2_values.add(m2)
        m4_values.add(m4)
        if len(sample_records) < 5:
            sample_records.append((bits, m2, m4, counts[:10]))

    print(f"  tournaments_checked={1 << comb(n, 2)}")
    print(f"  odd_centered_moment_failures={odd_failures}")
    print(
        "  distinct_even_m2_values="
        f"{len(m2_values)} range=[{format_fraction(min(m2_values))}, {format_fraction(max(m2_values))}]"
    )
    print(
        "  distinct_even_m4_values="
        f"{len(m4_values)} range=[{format_fraction(min(m4_values))}, {format_fraction(max(m4_values))}]"
    )
    print("  sample_tournaments(bits, m2, m4, first_10_fwd_counts):")
    for bits, m2, m4, first10 in sample_records:
        print(
            "    "
            f"bits={bits:4d} "
            f"m2={format_fraction(m2):>5s} "
            f"m4={format_fraction(m4):>5s} "
            f"counts={first10}"
        )
    print()

    print("5. Scout conclusion")
    print("  Anchor gate: midpoint antisymmetrization kills even Faulhaber channels.")
    print("  Tournament gate: reversal symmetry kills odd centered forward moments.")
    print("  This is a parity-dual pair of centered cancellation mechanisms.")


if __name__ == "__main__":
    main()
