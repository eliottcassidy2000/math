"""
trienerment_iso_count.py — Burnside counting for t(r)ienerment isomorphism classes.

A t(r)ienerment on n vertices is a complete graph where each edge {u,v} is assigned
one of three states: directed u->v, directed v->u, or bidirectional (tie) u<->v.

The number of unlabeled isomorphism classes T(n) = A007025.
The distribution f(n,k) counts classes with exactly k ties.
f(n,0) = A000568 (tournament isomorphism classes).

Usage:
  python3 trienerment_iso_count.py
  python3 trienerment_iso_count.py 8   # compute through n=8 (slower)
"""
import sys
from math import gcd, factorial
from collections import defaultdict, Counter


def partitions(n):
    """Yield all integer partitions of n as sorted tuples (largest first)."""
    def helper(n, max_val):
        if n == 0:
            yield ()
            return
        for k in range(min(n, max_val), 0, -1):
            for rest in helper(n - k, k):
                yield (k,) + rest
    yield from helper(n, n)


def cycle_count(partition):
    """Number of permutations in S_n with the given cycle type."""
    n = sum(partition)
    c = Counter(partition)
    result = factorial(n)
    for l, cnt in c.items():
        result //= (l ** cnt * factorial(cnt))
    return result


def poly_for_partition(partition):
    """
    Generating polynomial for fixed t(r)ienerments under a permutation with
    given cycle type. Coefficient of x^k = number of fixed t(r)ienerments
    with exactly k ties.

    The polynomial is:
      prod_{even cycles l} x^{l/2}          (Case A: l/2 forced-tie pairs per even cycle)
    x prod_{cycles l} (2 + x^l)^{floor((l-1)/2)}   (Case B within-cycle: floor((l-1)/2) pairs,
                                                      each covering l unordered pairs)
    x prod_{pairs of cycles (l_a,l_b)} (2 + x^{lcm(l_a,l_b)})^{gcd(l_a,l_b)}
                                                     (Case B cross-cycle)
    """
    parts = list(partition)

    def poly_mult(p1, p2):
        result = {}
        for d1, c1 in p1.items():
            for d2, c2 in p2.items():
                result[d1+d2] = result.get(d1+d2, 0) + c1*c2
        return result

    poly = {0: 1}

    # Case A: one forced-tie block of l/2 pairs per even-length cycle
    for l in parts:
        if l % 2 == 0:
            poly = poly_mult(poly, {l//2: 1})

    # Case B (within cycle): floor((l-1)/2) free pairs per cycle, each covering l pairs
    for l in parts:
        num_pairs = (l - 1) // 2
        factor = {0: 2, l: 1}   # 2 directed choices or 1 all-tie choice
        for _ in range(num_pairs):
            poly = poly_mult(poly, factor)

    # Case B (cross-cycle): gcd(l_a, l_b) pairs for each pair of cycles
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            la, lb = parts[i], parts[j]
            g = gcd(la, lb)
            lc = la * lb // g   # lcm
            factor = {0: 2, lc: 1}
            for _ in range(g):
                poly = poly_mult(poly, factor)

    return poly


def f_distribution(n):
    """Return dict {k: f(n,k)} where f(n,k) = number of iso classes with exactly k ties."""
    max_edges = n * (n - 1) // 2
    coeff_sums = defaultdict(int)
    for part in partitions(n):
        cnt = cycle_count(part)
        poly = poly_for_partition(part)
        for deg, coef in poly.items():
            coeff_sums[deg] += cnt * coef
    nfact = factorial(n)
    return {k: coeff_sums.get(k, 0) // nfact for k in range(max_edges + 1)}


def trienerment_count(n):
    return sum(f_distribution(n).values())


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 7

    print("T(r)ienerment isomorphism class counts")
    print("=" * 55)
    print(f"  n   T(n)=A007025    f(n,0)=A000568")
    for n in range(1, max_n + 1):
        dist = f_distribution(n)
        T = sum(dist.values())
        print(f"  {n}   {T:<14}  {dist.get(0, 0)}")

    print()
    print("Distribution f(n,k) — iso classes with exactly k ties")
    print("=" * 55)
    for n in range(2, min(max_n + 1, 8)):
        dist = f_distribution(n)
        M = n * (n - 1) // 2
        vals = [dist[k] for k in range(M + 1)]
        print(f"  n={n}: {vals}")

    print()
    print("Stabilized tail values f(n, C(n,2)-j) for small j:")
    print(f"  {'n':>3}  {'C(n,2)':>6}  k=M  k=M-1  k=M-2  k=M-3  k=M-4")
    for n in range(2, min(max_n + 1, 8)):
        dist = f_distribution(n)
        M = n * (n - 1) // 2
        row = [dist.get(M - j, 0) for j in range(5)]
        print(f"  {n:>3}  {M:>6}  {row[0]:4}  {row[1]:5}  {row[2]:5}  {row[3]:5}  {row[4]:5}")
