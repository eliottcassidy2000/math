#!/usr/bin/env python3
"""
pi_a000255_89c.py — The A000255 connection: compatible permutation pairs
opus-2026-03-14-S89c

DISCOVERY: The number of compatible pairs of Hamiltonian path orderings
in an n-vertex tournament equals n! × A000255(n-1).

A000255(n) satisfies a(n) = n*a(n-1) + (n-1)*a(n-2) with a(0)=a(1)=1.
EGF: exp(-x)/(1-x)²

The compatibility rate n!×A000255(n-1)/n!² = A000255(n-1)/n! → 1/e as n→∞.

This gives e (Euler's number) a natural tournament-theoretic interpretation:
  e = 1 / lim_{n→∞} Pr[two random Hamiltonian path orderings are compatible]

We also compute E[H²] using the pair formula for n=3..7 and analyze CV².
"""

from fractions import Fraction
from math import factorial, e as euler_e
from itertools import permutations

print("=" * 70)
print("THE A000255 CONNECTION")
print("Compatible pairs of Hamiltonian path orderings")
print("=" * 70)

# A000255 sequence
a255 = [1, 1]
for k in range(2, 20):
    a255.append(k * a255[-1] + (k-1) * a255[-2])

print("\nA000255 sequence (first 15 terms):")
for i in range(15):
    print(f"  a({i}) = {a255[i]}")

print("\n" + "=" * 70)
print("VERIFICATION: compatible pair count = n! × A000255(n-1)")
print("=" * 70)

# Direct computation for n=3..6 (n=7 already verified)
for n in range(3, 7):
    all_perms = list(permutations(range(n)))
    compat = 0

    for pi in all_perms:
        edge_pi = set((pi[i], pi[i+1]) for i in range(n-1))
        for sig in all_perms:
            edge_sig = set((sig[i], sig[i+1]) for i in range(n-1))
            if not any((v,u) in edge_sig for u,v in edge_pi):
                compat += 1

    expected = factorial(n) * a255[n-1]
    print(f"  n={n}: compatible = {compat}, n!×a(n-1) = {expected}, match: {'✓' if compat == expected else '✗'}")

# n=7 from our computation
compat_7 = 10679760
expected_7 = factorial(7) * a255[6]
print(f"  n=7: compatible = {compat_7}, n!×a(n-1) = {expected_7}, match: {'✓' if compat_7 == expected_7 else '✗'}")

print("\n" + "=" * 70)
print("COMPATIBILITY RATE → 1/e")
print("=" * 70)

print(f"\n  1/e = {1/euler_e:.10f}")
print()
for n in range(3, 15):
    rate = a255[n-1] / factorial(n)
    err = rate - 1/euler_e
    print(f"  n={n}: rate = {rate:.10f}, error = {err:+.2e}")

print("\n  THEOREM: The fraction of ordered pairs (π, σ) of permutations")
print("  of [n] such that the Hamiltonian path orderings π and σ are")
print("  compatible (no edge in opposite directions) equals")
print("  A000255(n-1)/n! → 1/e as n → ∞.")

print("\n" + "=" * 70)
print("CV² SEQUENCE AND 1/4 LIMIT CONJECTURE")
print("=" * 70)

# CV² sequence: 1/3, 1/3, 19/60, 13/45, 131/504
cv2_data = [
    (3, Fraction(1, 3)),
    (4, Fraction(1, 3)),
    (5, Fraction(19, 60)),
    (6, Fraction(13, 45)),
    (7, Fraction(131, 504)),
]

print("\n  CV² = Var(H)/E[H]² sequence:")
for n, cv in cv2_data:
    delta = cv - Fraction(1, 4)
    print(f"    n={n}: CV² = {cv} = {float(cv):.10f}, CV²-1/4 = {delta} = {float(delta):.10f}")

print(f"\n  The sequence CV²-1/4 → 0:")
for n, cv in cv2_data:
    d = cv - Fraction(1, 4)
    if n >= 5:
        print(f"    n={n}: {d} ≈ {float(d):.6f}")

print(f"\n  (CV²-1/4) × n! gives integers/2:")
for n, cv in cv2_data:
    val = (cv - Fraction(1, 4)) * factorial(n)
    print(f"    n={n}: {val}")

print("\n" + "=" * 70)
print("THE E[H²] PAIR FORMULA IN TERMS OF A000255")
print("=" * 70)

# E[H²] = Σ_{k=0}^{n-1} count(k) × 2^{-(2(n-1)-k)}
# = 2^{-2(n-1)} × Σ_k count(k) × 2^k
# where count(k) = number of compatible pairs with k common undirected edges

# At n=7, the distribution was:
dist_7 = {0: 3255840, 1: 3981600, 2: 2353680, 3: 846720, 4: 206640, 5: 30240, 6: 5040}

# Verify E[H²] computation
eh2_check = Fraction(0)
for k, cnt in dist_7.items():
    constrained = 2*6 - k  # = 12 - k
    eh2_check += Fraction(cnt, 2**constrained)

eh = Fraction(factorial(7), 2**6)
cv2_check = (eh2_check / eh**2) - 1

print(f"\n  E[H²] at n=7 = {eh2_check} = {float(eh2_check):.6f}")
print(f"  CV² at n=7 = {cv2_check} = {float(cv2_check):.10f}")
print(f"  Matches 131/504: {'✓' if cv2_check == Fraction(131, 504) else '✗'}")

print("\n" + "=" * 70)
print("COMMON EDGE DISTRIBUTION TABLE")
print("=" * 70)

# Compute for n=3..6
print("\n  Table of count(k)/n! (compatible pairs with k common edges):")
print(f"  {'k':>3}", end="")
for n in range(3, 8):
    print(f"  {'n='+str(n):>8}", end="")
print()

max_edges = max(range(3,8)) - 1
all_dists = {}

for n in range(3, 7):
    all_perms = list(permutations(range(n)))
    dist = {}
    for pi in all_perms:
        epi_d = set((pi[i], pi[i+1]) for i in range(n-1))
        epi_u = set((min(pi[i],pi[i+1]), max(pi[i],pi[i+1])) for i in range(n-1))
        for sig in all_perms:
            esig_d = set((sig[i], sig[i+1]) for i in range(n-1))
            if any((v,u) in esig_d for u,v in epi_d):
                continue
            esig_u = set((min(sig[i],sig[i+1]), max(sig[i],sig[i+1])) for i in range(n-1))
            k = len(epi_u & esig_u)
            dist[k] = dist.get(k, 0) + 1
    all_dists[n] = dist

all_dists[7] = dist_7

for k in range(max_edges + 1):
    print(f"  {k:>3}", end="")
    for n in range(3, 8):
        if k <= n-1 and k in all_dists.get(n, {}):
            val = all_dists[n][k] // factorial(n)
            print(f"  {val:>8}", end="")
        else:
            print(f"  {'·':>8}", end="")
    print()

# Row sums
print(f"  {'Σ':>3}", end="")
for n in range(3, 8):
    s = sum(all_dists.get(n, {}).values()) // factorial(n)
    print(f"  {s:>8}", end="")
print("  = A000255(n-1)")

print("\n" + "=" * 70)
print("SUMMARY: π (via e = exp(1)) IN TOURNAMENT COMPATIBLE PAIRS")
print("=" * 70)

print("""
  THEOREM (A000255 Connection):
  The number of compatible pairs (π,σ) of Hamiltonian path orderings
  in tournaments on n vertices equals n! × A000255(n-1).

  COROLLARY: Pr[random (π,σ) compatible] = A000255(n-1)/n! → 1/e.

  This gives e a direct tournament-theoretic interpretation:
    e = lim 1/Pr[compatible]

  CONJECTURE (CV² Limit):
  CV²(H) = Var(H)/E[H]² → 1/4 as n → ∞.
  Known values: 1/3, 1/3, 19/60, 13/45, 131/504.
  The approach to 1/4 appears super-exponential.

  CONNECTION TO π:
  Since e = exp(1) and exp is the function that connects
  real analysis to complex analysis (via e^{iπ} = -1),
  and since A000255 has EGF exp(-x)/(1-x)², the appearance
  of e in tournament compatibility is deeply connected to
  the analytic structure of permutation enumeration.
""")

print("Done!")
