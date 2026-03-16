#!/usr/bin/env python3
"""
self_comp_tournaments_s90dl.py — Self-converse tournament enumeration
opus-2026-03-16-S90dl

A tournament T on n vertices is SELF-CONVERSE (= self-complementary) if
T is isomorphic to T^op (all arcs reversed).

Counts:
  SC(n) = number of self-converse isomorphism classes (OEIS A002785)
  T(n)  = total tournament isomorphism classes (OEIS A000568)
  paired = (T(n) - SC(n)) / 2   (non-self-converse classes come in pairs)

Formula for SC(n) [from A002785, Andrew Howroyd]:
  The self-converse count uses Burnside on the wreath product.
  We sum over partitions p of m = n//2 with ALL PARTS ODD.
  Each such partition p = (p_1,...,p_k) of m defines a cycle type
  2*p = (2p_1,...,2p_k) in S_{2m} (or S_{2m+1} for odd n).

  SC(n) = (1/n!) * sum_p [ permcount(2*p) * 2^edges(p) * extra(n,p) ]

  where:
    permcount(v) = (sum v)! / prod_i (v_i^{m_i} * m_i!)
                   [conjugacy class size for cycle type v in S_{sum v}]
    edges(p) = sum_i p_i + 2 * sum_{i<j} gcd(p_i, p_j)
    extra(n,p) = 1             if n is even
               = n * 2^len(p)  if n is odd

  Only partitions with all odd parts contribute (even-length cycles
  cannot pair with the converse involution).

Formula for T(n) [A000568]:
  T(n) = (1/n!) * sum over partitions lambda of n INTO ODD PARTS:
         |C_lambda| * 2^{edges(lambda)}
  where edges(lambda) = sum_i floor(l_i/2) + sum_{i<j} gcd(l_i, l_j)
  (Permutations with even-length cycles fix zero tournaments.)
"""

from math import gcd, factorial
from collections import Counter
import time, sys, os

# ── partition generators and helpers ──────────────────────────────────
def odd_partitions(n, max_part=None):
    """Partitions of n into odd parts only."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    top = min(n, max_part)
    if top % 2 == 0:
        top -= 1
    for first in range(top, 0, -2):
        for rest in odd_partitions(n - first, first):
            yield (first,) + rest

def permcount(v):
    """Size of conjugacy class in S_N with cycle type v."""
    counts = Counter(v)
    denom = 1
    for l, m in counts.items():
        denom *= (l ** m) * factorial(m)
    return factorial(sum(v)) // denom

# ── T(n): total tournaments (A000568) ────────────────────────────────
def edges_T(parts):
    """Orbit count on arcs for a permutation with cycle type `parts`.
    Only valid when all parts are odd.
    edges = sum_i floor(l_i/2) + sum_{i<j} gcd(l_i, l_j)
    """
    s = sum(p // 2 for p in parts)
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            s += gcd(parts[i], parts[j])
    return s

def T(n):
    """Number of non-isomorphic tournaments on n vertices (A000568).
    Only partitions into odd parts contribute (even cycles fix 0 tournaments).
    """
    if n <= 1:
        return 1
    total = 0
    for parts in odd_partitions(n):
        total += permcount(parts) * (1 << edges_T(parts))
    return total // factorial(n)

# ── SC(n): self-converse tournaments (A002785) ───────────────────────
def SC(n):
    """Number of self-converse tournament isomorphism classes on n vertices."""
    if n <= 1:
        return 1
    m = n // 2
    total = 0
    for p in odd_partitions(m):
        # cycle type in S_n is 2*p (each part doubled)
        doubled = tuple(2 * x for x in p)
        # edges exponent: sum p_i + 2 * sum_{i<j} gcd(p_i, p_j)
        e = sum(p)
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                e += 2 * gcd(p[i], p[j])
        # extra factor for odd n
        if n % 2 == 1:
            extra = n * (1 << len(p))
        else:
            extra = 1
        total += permcount(doubled) * (1 << e) * extra
    return total // factorial(n)

# ── known values for validation ──────────────────────────────────────
KNOWN_SC = {
    1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176,
    9: 2752, 10: 8784, 11: 279968, 12: 1492288,
    13: 95458560, 14: 872687552
}
KNOWN_T = {
    0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456,
    8: 6880, 9: 191536, 10: 9733056, 11: 903753248,
    12: 154108311168
}

# ── main computation ─────────────────────────────────────────────────
def main():
    NMAX = 30
    out_lines = []

    def pr(s=""):
        print(s)
        out_lines.append(s)

    pr("Self-converse (self-complementary) tournament enumeration")
    pr("OEIS A002785 (SC) and A000568 (total)")
    pr("=" * 78)
    pr(f"{'n':>3} | {'T(n)':>22} | {'SC(n)':>22} | {'paired':>22} | {'chk'}")
    pr("-" * 78)

    for n in range(1, NMAX + 1):
        t0 = time.time()
        tn = T(n)
        sc = SC(n)
        elapsed = time.time() - t0
        paired = (tn - sc) // 2

        # validation
        chk = ""
        if n in KNOWN_SC:
            chk = "OK" if sc == KNOWN_SC[n] else f"FAIL(sc={KNOWN_SC[n]})"
        if n in KNOWN_T:
            tc = "OK" if tn == KNOWN_T[n] else f"FAIL(T={KNOWN_T[n]})"
            chk = tc + " " + chk

        # sanity: (T - SC) must be even
        if (tn - sc) % 2 != 0:
            chk += " ERR:odd"

        pr(f"{n:3d} | {tn:22d} | {sc:22d} | {paired:22d} | {chk}")
        if elapsed > 0.5:
            pr(f"      ({elapsed:.1f}s)")

    pr()
    pr("Ratios SC(n)/T(n):")
    pr(f"{'n':>3} | {'SC/T':>12} | {'notes'}")
    pr("-" * 40)
    for n in range(1, NMAX + 1):
        tn = T(n)
        sc = SC(n)
        if tn > 0:
            ratio = sc / tn
            note = ""
            if ratio == 1.0:
                note = "all self-converse"
            elif ratio == 0.5:
                note = "exactly half"
            pr(f"{n:3d} | {ratio:12.8f} | {note}")

    pr()
    pr("b-file format (A002785):")
    for n in range(1, NMAX + 1):
        pr(f"{n} {SC(n)}")

    # save output
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', '05-knowledge', 'results',
                           'self_comp_tournaments_s90dl.out')
    with open(outpath, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')
    print(f"\nSaved to {outpath}")

if __name__ == '__main__':
    main()
