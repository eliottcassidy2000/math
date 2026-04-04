#!/usr/bin/env python3
"""
a000568_naive.py — Brute-force computation of OEIS A000568
(Number of non-isomorphic tournaments on n labeled vertices)

This is the NAIVE approach: enumerate all 2^C(n,2) tournaments,
check pairwise isomorphism, count equivalence classes.

Complexity: O(2^{n(n-1)/2} * n!) — feasible only for n <= 7 or so.

This file exists for pedagogical contrast with the optimized version.
"""

from itertools import permutations, combinations
import time
import sys


def all_tournaments(n):
    """Generate all tournaments on n vertices as frozensets of directed arcs.

    A tournament is a complete directed graph: for every pair {u,v},
    exactly one of (u,v) or (v,u) is present.

    There are C(n,2) pairs, so 2^C(n,2) tournaments total.
    """
    pairs = list(combinations(range(n), 2))
    m = len(pairs)  # = n*(n-1)//2

    for bits in range(1 << m):
        arcs = set()
        for k, (u, v) in enumerate(pairs):
            if (bits >> k) & 1:
                arcs.add((u, v))
            else:
                arcs.add((v, u))
        yield frozenset(arcs)


def apply_perm(tournament, perm):
    """Apply vertex permutation to a tournament.

    If sigma is a permutation, the image of arc (u,v) is (sigma(u), sigma(v)).
    """
    return frozenset((perm[u], perm[v]) for (u, v) in tournament)


def canonical_form(tournament, n):
    """Compute canonical form by trying all n! permutations.

    The canonical form is the lexicographically smallest adjacency encoding
    over all permutations of the vertices. Two tournaments are isomorphic
    iff they have the same canonical form.

    This is the simplest correct algorithm — and the slowest.
    """
    best = None
    pairs = list(combinations(range(n), 2))

    for perm in permutations(range(n)):
        # Apply permutation, encode as tuple of bits
        image = apply_perm(tournament, perm)
        encoding = tuple(1 if (u, v) in image else 0 for u, v in pairs)
        if best is None or encoding < best:
            best = encoding
    return best


def count_tournaments_brute(n):
    """Count non-isomorphic tournaments by brute force.

    Enumerate all tournaments, canonicalize each one, count distinct forms.
    """
    if n <= 1:
        return 1

    seen = set()
    total = 0

    for T in all_tournaments(n):
        total += 1
        canon = canonical_form(T, n)
        seen.add(canon)

    return len(seen), total


# ============================================================
# Run
# ============================================================

KNOWN = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

print("OEIS A000568 — BRUTE FORCE ENUMERATION")
print("=" * 60)
print()
print("Method: enumerate all 2^C(n,2) tournaments, canonicalize")
print("        each via all n! permutations, count distinct forms.")
print()
print(f"{'n':>3} | {'T(n)':>8} | {'tournaments':>14} | {'time':>10} | {'check':>6}")
print("-" * 58)

max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 7

for n in range(max_n + 1):
    t0 = time.time()
    if n <= 1:
        result, total = 1, 1
    else:
        result, total = count_tournaments_brute(n)
    elapsed = time.time() - t0

    expected = KNOWN.get(n)
    check = "OK" if expected is None or result == expected else "FAIL"

    print(f"{n:3d} | {result:8d} | {total:14,} | {elapsed:8.2f}s | {check:>6}")

    if elapsed > 120:
        print(f"\n  Stopping: n={n} took {elapsed:.0f}s. The wall is 2^C(n,2) * n!")
        print(f"  n={n+1} would enumerate {2**((n+1)*n//2):,} tournaments")
        print(f"  with {__import__('math').factorial(n+1):,} permutations each.")
        break

print(f"""
COMPLEXITY ANALYSIS:
  n=5:  2^10 * 5!   =      1,024 * 120      =        122,880 operations
  n=6:  2^15 * 6!   =     32,768 * 720      =     23,592,960 operations
  n=7:  2^21 * 7!   =  2,097,152 * 5,040    = 10,569,646,080 operations
  n=8:  2^28 * 8!   = COMPLETELY INFEASIBLE  (2.7 trillion operations)

  The brute force approach hits a WALL at n=7 or n=8.
  To go further, we need a fundamentally different algorithm.
""")
