#!/usr/bin/env python3
"""
prime_numerator_deep_dive.py — kind-pasteur-2026-03-21-S16d

THE PRIME NUMERATOR MYSTERY.

KEY INSIGHT from overnight session: The overlap structure of two
Hamiltonian paths (permutations) decomposes into:

  "agreement" = #{i : tau(i+1) = tau(i) + 1}  (unit ASCENT in tau, rel to identity path)
  "conflict"  = #{i : tau(i+1) = tau(i) - 1}  (unit DESCENT in tau, rel to identity path)

For a FIXED reference path sigma = identity (0->1->2->...->n-1):
  - A unit ascent at position i in tau means tau uses the SAME directed edge (k, k+1)
  - A unit descent means tau uses the REVERSE edge (k+1, k) = conflict

So f(n, a) = n! * #{tau in S_n : exactly a unit ascents AND 0 unit descents}

A "unit ascent" at position i means tau(i+1) - tau(i) = +1.
A "unit descent" at position i means tau(i+1) - tau(i) = -1.

The condition "0 unit descents" means: tau AVOIDS the pattern where
consecutive elements decrease by exactly 1. This is a "forbidden
adjacency" / "succession-free" condition in permutation combinatorics.

THE PLAN:
1. Count permutations by (unit_ascents, unit_descents) for n=3..9
2. Extract f(n,a) = n! * count(a, 0)
3. Compute exact Var(H) and CV^2(H) at n=9
4. Factor numerator -> find the prime
5. Connect to known combinatorial sequences (OEIS)
6. Search for topological interpretation
"""
import sys
from itertools import permutations
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial, gcd, comb
import time

sys.stdout.reconfigure(line_buffering=True)

def count_unit_ascents_descents(perm):
    """Count unit ascents (tau(i+1)-tau(i)=+1) and unit descents (=-1)."""
    asc = 0
    desc = 0
    for i in range(len(perm) - 1):
        diff = perm[i+1] - perm[i]
        if diff == 1:
            asc += 1
        elif diff == -1:
            desc += 1
    return asc, desc

print("=" * 72)
print("  THE PRIME NUMERATOR MYSTERY: DEEP DIVE")
print("  kind-pasteur-2026-03-21-S16d")
print("=" * 72)

# ============================================================
# PART 1: Count permutations by (ascents, descents) for each n
# ============================================================

print("\n" + "=" * 72)
print("  PART 1: Permutations by (unit_ascents, unit_descents)")
print("=" * 72)

# For efficiency at n=9: only need to iterate over n! permutations
# (not n!^2 pairs), since we fix sigma = identity.

asc_desc_table = {}  # n -> {(a, d): count}
f_table = {}  # n -> {a: f(n,a)} where f(n,a) = n! * count(a, 0)

for n in range(3, 10):
    t0 = time.time()
    counts = Counter()

    for perm in permutations(range(n)):
        a, d = count_unit_ascents_descents(perm)
        counts[(a, d)] += 1

    asc_desc_table[n] = counts
    elapsed = time.time() - t0

    # Extract f(n, a) = n! * count(a, d=0)
    f_n = {}
    for (a, d), cnt in counts.items():
        if d == 0:
            f_n[a] = factorial(n) * cnt

    f_table[n] = f_n

    print(f"\n  n={n} ({factorial(n)} perms, {elapsed:.1f}s):")
    print(f"    (asc, desc) distribution:")
    for (a, d) in sorted(counts.keys()):
        print(f"      ({a},{d}): {counts[(a,d)]}", end="")
        if d == 0:
            print(f"  <-- compatible, f({n},{a}) = {f_n[a]}", end="")
        print()

    # Total compatible (d=0)
    total_compat = sum(cnt for (a, d), cnt in counts.items() if d == 0)
    print(f"    Total d=0: {total_compat}/{factorial(n)} = {Fraction(total_compat, factorial(n))}")

# ============================================================
# PART 2: Compute exact Var(H) from f(n, a)
# ============================================================

print("\n" + "=" * 72)
print("  PART 2: Exact Var(H) from f(n, a)")
print("=" * 72)

for n in range(3, 10):
    f_n = f_table[n]

    # E[H^2] = 2^{-2(n-1)} * sum_a f(n,a) * 2^a
    weighted_sum = sum(f_a * (2**a) for a, f_a in f_n.items())
    E_H2 = Fraction(weighted_sum, 4**(n-1))

    E_H = Fraction(factorial(n), 2**(n-1))
    Var_H = E_H2 - E_H**2

    # CV^2(H)
    CV2 = Var_H / E_H**2

    print(f"\n  n={n}:")
    print(f"    E[H] = {E_H} = {float(E_H):.6f}")
    print(f"    E[H^2] = {E_H2} = {float(E_H2):.6f}")
    print(f"    Var(H) = {Var_H}")
    print(f"    CV^2(H) = {CV2}")

    # Factor numerator and denominator of CV^2
    p_num = CV2.numerator
    p_den = CV2.denominator

    def factor(x):
        if x <= 1: return {}
        f = {}
        d = 2
        while d*d <= x:
            while x % d == 0:
                f[d] = f.get(d,0)+1
                x //= d
            d += 1
        if x > 1: f[x] = f.get(x,0)+1
        return f

    print(f"    CV^2 numerator = {p_num}, factors: {factor(p_num)}")
    print(f"    CV^2 denominator = {p_den}, factors: {factor(p_den)}")

    # Compute OCR
    # R^2 = 2(n-2) / [n(n-1) * CV^2]
    R2 = Fraction(2*(n-2), 1) / (n * (n-1) * CV2)
    print(f"    R^2 = {R2} = {float(R2):.12f}")
    print(f"    1-R^2 = {1-R2}")

# ============================================================
# PART 3: The sequence of d=0 counts (no unit descents)
# ============================================================

print("\n" + "=" * 72)
print("  PART 3: Permutations with NO unit descents (OEIS search target)")
print("=" * 72)

print("\n  Sequence: number of permutations of [n] with no unit descents:")
d0_counts = []
for n in range(3, 10):
    counts = asc_desc_table[n]
    total = sum(cnt for (a, d), cnt in counts.items() if d == 0)
    d0_counts.append(total)
    print(f"    n={n}: {total}")

print(f"\n  Sequence: {d0_counts}")
print(f"  (Search OEIS for: {', '.join(map(str, d0_counts))})")

# ============================================================
# PART 4: Refined: permutations with EXACTLY a unit ascents, 0 unit descents
# ============================================================

print("\n" + "=" * 72)
print("  PART 4: Count(a ascents, 0 descents) table")
print("=" * 72)

print(f"\n  {'n':>3} |", end="")
for a in range(10):
    print(f"  a={a:d}  ", end="")
print()
print("  " + "-" * 80)

for n in range(3, 10):
    counts = asc_desc_table[n]
    print(f"  {n:3d} |", end="")
    for a in range(min(n, 10)):
        cnt = counts.get((a, 0), 0)
        print(f"  {cnt:5d}", end="")
    print()

# The a=n-1 column (all unit ascents, 0 descents): only the identity permutation
# The a=0 column: permutations with no unit ascents AND no unit descents

print(f"\n  a=0, d=0 sequence (no consecutive elements differ by 1):")
smarandache = []
for n in range(3, 10):
    cnt = asc_desc_table[n].get((0, 0), 0)
    smarandache.append(cnt)
    print(f"    n={n}: {cnt}")
print(f"  Sequence: {smarandache}")
print(f"  (These are 'Smarandache' or 'discordant' permutations)")

# ============================================================
# PART 5: THE PRIME SEQUENCE
# ============================================================

print("\n" + "=" * 72)
print("  PART 5: THE PRIME SEQUENCE IN CV^2(H)")
print("=" * 72)

primes_found = []
for n in range(3, 10):
    f_n = f_table[n]
    weighted_sum = sum(f_a * (2**a) for a, f_a in f_n.items())
    E_H2 = Fraction(weighted_sum, 4**(n-1))
    E_H = Fraction(factorial(n), 2**(n-1))
    Var_H = E_H2 - E_H**2
    CV2 = Var_H / E_H**2

    p_num = CV2.numerator
    factors = factor(p_num)
    largest_prime = max(factors.keys()) if factors else 1
    primes_found.append(largest_prime)

    print(f"  n={n}: CV^2 = {CV2}, numerator = {p_num}, largest prime = {largest_prime}")

print(f"\n  THE PRIME SEQUENCE: {primes_found}")
print(f"  n:     3  4  5   6   7    8    9")
print(f"  prime: {primes_found}")

# Check if n=9 prime is still 131 or changes
if len(primes_found) >= 7:
    p9 = primes_found[6]  # n=9
    print(f"\n  *** n=9 prime = {p9} ***")
    if p9 == 131:
        print(f"  131 PERSISTS at n=9! The plateau may be real after all.")
    else:
        print(f"  NEW PRIME at n=9! The sequence changes: ...131, 131, {p9}")

# ============================================================
# PART 6: Residual spreading interpretation
# ============================================================

print("\n" + "=" * 72)
print("  PART 6: What the prime MEANS — residual spreading")
print("=" * 72)

print("""
  CV^2(H) = Var(H) / E[H]^2 = p / D

  The numerator p is the "irreducible complexity" of H's distribution.
  It measures how much the variance of Hamiltonian path counts CANNOT
  be reduced by any polynomial combination of simpler invariants.

  ALGEBRAIC INTERPRETATION:
  The permutation group S_n acts on the space of tournaments.
  H is a polynomial function on this space.
  Var(H) = <H - E[H], H - E[H]> in the L^2 inner product.

  The "residual" H - E[H] lives in a representation of S_n.
  The prime p(n) is related to the DIMENSION or STRUCTURE of
  the irreducible representations that appear in this decomposition.

  TOPOLOGICAL INTERPRETATION:
  The compatible pairs (sigma, tau) with a agreements form a
  simplicial complex. The Euler characteristic of this complex
  may equal or be related to the prime p(n).

  Consider: at n=7, the complex of compatible pairs has
  Euler characteristic involving 131. The "hole" structure
  of this complex determines how much variance "leaks through"
  the score prediction.
""")
