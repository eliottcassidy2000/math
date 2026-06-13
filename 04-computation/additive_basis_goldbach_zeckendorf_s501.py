#!/usr/bin/env python3
"""Compare additive-basis regimes around Goldbach, polygonal numbers, and Zeckendorf.

The point is not to prove any classical theorem.  It is to put the four
phenomena in one finite coordinate system:

  * Goldbach / Hardy-Littlewood: dense prime basis, many representations.
  * Helfgott / weak Goldbach: one extra prime variable smooths the local noise.
  * Fermat polygonal: polynomially sparse bases, bounded arity.
  * Zeckendorf: exponentially sparse basis, variable arity, unique normal form.
"""

from __future__ import annotations

from collections import Counter
from math import gcd, isclose, log


N = 5000


def sieve(limit: int) -> list[int]:
    is_prime = [True] * (limit + 1)
    is_prime[0:2] = [False, False]
    p = 2
    while p * p <= limit:
        if is_prime[p]:
            for q in range(p * p, limit + 1, p):
                is_prime[q] = False
        p += 1
    return [i for i, ok in enumerate(is_prime) if ok]


def prime_pair_counts(primes: list[int], limit: int) -> list[int]:
    counts = [0] * (limit + 1)
    for p in primes:
        for q in primes:
            s = p + q
            if s <= limit:
                counts[s] += 1
    return counts


def prime_triple_counts(primes: list[int], pair_counts: list[int], limit: int) -> list[int]:
    counts = [0] * (limit + 1)
    for n in range(limit + 1):
        total = 0
        for p in primes:
            if p > n:
                break
            total += pair_counts[n - p]
        counts[n] = total
    return counts


def factor_distinct(n: int) -> list[int]:
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(n)
    return out


def twin_prime_constant(primes: list[int], cutoff: int) -> float:
    c = 1.0
    for p in primes:
        if p == 2:
            continue
        if p > cutoff:
            break
        c *= 1.0 - 1.0 / ((p - 1) ** 2)
    return c


def goldbach_prediction(n: int, c2: float) -> float:
    local = 1.0
    for p in factor_distinct(n):
        if p > 2:
            local *= (p - 1) / (p - 2)
    return 2 * c2 * n / (log(n) ** 2) * local


def polygonal(g: int, k: int) -> int:
    return ((g - 2) * k * k - (g - 4) * k) // 2


def polygonal_terms(g: int, limit: int) -> list[int]:
    vals = []
    k = 0
    while True:
        v = polygonal(g, k)
        if v > limit:
            break
        vals.append(v)
        k += 1
    return vals


def min_terms_for_basis(vals: list[int], limit: int) -> list[int]:
    inf = 10**9
    dp = [inf] * (limit + 1)
    dp[0] = 0
    positive = [v for v in vals if v > 0]
    for n in range(1, limit + 1):
        dp[n] = 1 + min(dp[n - v] for v in positive if v <= n)
    return dp


def fibs_upto(limit: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] + fibs[-2] <= limit:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


def zeckendorf(n: int, fibs: list[int]) -> list[int]:
    out = []
    remaining = n
    for f in reversed(fibs):
        if f <= remaining:
            out.append(f)
            remaining -= f
    assert remaining == 0
    return list(reversed(out))


def repeated_sumset(residues: set[int], k: int, m: int) -> set[int]:
    sums = {0}
    for _ in range(k):
        sums = {(a + b) % m for a in sums for b in residues}
    return sums


def zeckendorf_residues_mod(m: int, steps: int) -> set[int]:
    # DP over legal no-adjacent Fibonacci digit strings.
    fibs = [1, 2]
    while len(fibs) < steps:
        fibs.append(fibs[-1] + fibs[-2])
    states = {(0, 0)}  # (previous digit used, residue)
    for f in fibs:
        nxt = set()
        for prev, r in states:
            nxt.add((0, r))
            if prev == 0:
                nxt.add((1, (r + f) % m))
        states = nxt
    return {r for _, r in states}


def residue_table(max_mod: int = 18) -> list[tuple[int, str, str, str, str, str, str]]:
    rows = []
    for m in range(2, max_mod + 1):
        units = {r for r in range(m) if gcd(r, m) == 1}
        u2 = repeated_sumset(units, 2, m)
        u3 = repeated_sumset(units, 3, m)
        tri = {polygonal(3, k) % m for k in range(m)}
        sq = {polygonal(4, k) % m for k in range(m)}
        pent = {polygonal(5, k) % m for k in range(m)}
        rows.append(
            (
                m,
                f"{len(u2)}/{m}",
                f"{len(u3)}/{m}",
                "Y" if len(repeated_sumset(tri, 3, m)) == m else "n",
                "Y" if len(repeated_sumset(sq, 4, m)) == m else "n",
                "Y" if len(repeated_sumset(pent, 5, m)) == m else "n",
                "Y" if len(zeckendorf_residues_mod(m, 6 * m + 8)) == m else "n",
            )
        )
    return rows


def main() -> None:
    primes = sieve(N)
    prime_set = set(primes)
    pair_counts = prime_pair_counts(primes, N)
    triple_counts = prime_triple_counts(primes, pair_counts, N)
    c2 = twin_prime_constant(primes, N)

    evens = list(range(4, N + 1, 2))
    odds = list(range(7, N + 1, 2))
    missing_even = [n for n in evens if pair_counts[n] == 0]
    missing_odd = [n for n in odds if triple_counts[n] == 0]

    print("=== Additive basis comparison S501 ===")
    print(f"N={N}, primes={len(primes)}, largest_prime={primes[-1]}")
    print()

    print("Goldbach / Helfgott finite audit")
    print(f"  binary Goldbach missing evens: {missing_even[:10]} count={len(missing_even)}")
    print(f"  ternary Goldbach missing odds: {missing_odd[:10]} count={len(missing_odd)}")
    sample_evens = [100, 210, 420, 1000, 2310, 4096, 5000]
    print("  even n: ordered prime-pair count vs Hardy-Littlewood shape")
    for n in sample_evens:
        pred = goldbach_prediction(n, c2)
        ratio = pair_counts[n] / pred if pred else float("nan")
        print(
            f"    n={n:4d} r2={pair_counts[n]:4d} HL~{pred:7.2f} "
            f"ratio={ratio:5.2f} local_primes={factor_distinct(n)}"
        )
    sample_odds = [101, 211, 421, 1001, 2311, 4097, 4999]
    print("  odd n: ordered prime-triple count, normalized by n^2/log^3 n")
    for n in sample_odds:
        scale = n * n / (log(n) ** 3)
        print(f"    n={n:4d} r3={triple_counts[n]:7d} r3/scale={triple_counts[n]/scale:6.2f}")
    print()

    print("Fermat polygonal finite audit")
    for g in range(3, 9):
        vals = polygonal_terms(g, N)
        dp = min_terms_for_basis(vals, N)
        max_terms = max(dp[1:])
        hist = Counter(dp[1:])
        tight = [i for i in range(1, N + 1) if dp[i] == max_terms][:8]
        print(
            f"  {g}-gonal: atoms<={N}: {len(vals)-1:3d}, "
            f"max_terms={max_terms}, hist={dict(sorted(hist.items()))}, first_tight={tight}"
        )
    print()

    print("Zeckendorf normal-form audit")
    fibs = fibs_upto(N)
    digit_counts = []
    top_gaps = []
    for n in range(1, N + 1):
        z = zeckendorf(n, fibs)
        digit_counts.append(len(z))
        idx = [fibs.index(v) for v in z]
        if len(idx) > 1:
            top_gaps.extend(b - a for a, b in zip(idx, idx[1:]))
    hist = Counter(digit_counts)
    print(f"  Fibonacci atoms<={N}: {len(fibs)} ({fibs})")
    print(f"  digit_count_hist={dict(sorted(hist.items()))}")
    print(f"  max_digits={max(digit_counts)}, avg_digits={sum(digit_counts)/len(digit_counts):.3f}")
    print(f"  all adjacent gaps >=2: {all(g >= 2 for g in top_gaps)}")
    print()

    print("Residue smoothing table")
    print("  m | U+U | U+U+U | tri^3 | sq^4 | pent^5 | Zeck")
    for row in residue_table():
        print(f" {row[0]:2d} | {row[1]:>4s} | {row[2]:>5s} | {row[3]:>5s} | {row[4]:>4s} | {row[5]:>6s} | {row[6]:>4s}")
    print()

    print("Synthesis")
    print("  primes: dense atoms, many reps; ternary adds smoothing and enables Helfgott.")
    print("  polygonals: sparse polynomial atoms, but bounded arity kills local obstructions.")
    print("  Zeckendorf: exponentially sparse atoms; uniqueness replaces averaging.")
    print("  common object: representation hypergraph; the deep choice is smoothing vs normal form.")


if __name__ == "__main__":
    main()
