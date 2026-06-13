#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Algebraic proof of the key inequality for THM-305.

Need to prove: For every partition λ = (c_1,...,c_k) of odd n into odd parts,
with λ ≠ (n) and k ≥ 2:

    G(λ) - v_2(z_λ) > (k-1)/2

where G(λ) = Σ_{i<j} gcd(c_i, c_j) and z_λ = (Π c_i)(Π a_j!).

STRATEGY: Prove a STRONGER inequality that's easier to analyze:
    G(λ) ≥ k(k-1)/2  (since gcd ≥ 1 for each pair)
    v_2(z_λ) ≤ Σ v_2(a_j!) ≤ Σ (a_j - 1) = k - d  (where d = #distinct parts)

So G - v_2(z) ≥ k(k-1)/2 - k + d = k(k-3)/2 + d.
Need this > (k-1)/2, i.e., k(k-3)/2 + d > (k-1)/2, i.e., k²-3k+2d > k-1,
i.e., k²-4k+2d+1 > 0, i.e., (k-2)² + 2d - 3 > 0.

For d ≥ 2: (k-2)² + 2d - 3 ≥ (k-2)² + 1 > 0. ✓
For d = 1 (all parts equal, c^k = n): (k-2)² + 2 - 3 = (k-2)² - 1.
  k ≥ 4: (k-2)² - 1 ≥ 3 > 0. ✓
  k = 3: (1)² - 1 = 0. TIED! Need sharper bound.
  k = 2: (0)² - 1 = -1 < 0. FAILS! Need sharper bound.

So for d=1 (all equal parts), k ∈ {2, 3}: need a case-specific argument.

Case d=1, k=2: λ = (c, c) with 2c = n. All parts equal to c = n/2.
  G = gcd(c,c) = c. v_2(z) = v_2(c² × 2!) = v_2(c²) + 1 = 2v_2(c) + 1.
  Since n is odd and n = 2c, c = n/2 is NOT an integer. Contradiction!
  So this case is IMPOSSIBLE (n odd, can't partition into two equal parts summing to n).

Case d=1, k=3: λ = (c, c, c) with 3c = n.
  G = 3 × gcd(c,c) = 3c. v_2(z) = v_2(c³ × 3!) = 3v_2(c) + 1.
  Since c is odd (n = 3c, n odd implies c odd), v_2(c) = 0.
  So v_2(z) = 0 + 1 = 1. G - v_2(z) = 3c - 1. (k-1)/2 = 1.
  Need: 3c - 1 > 1, i.e., c ≥ 1. True since c ≥ 1 (actually c ≥ 3 since parts odd, n ≥ 9).
  GAP = 3c - 2 ≥ 7 for n ≥ 9. VERY safe. ✓

So the only dangerous case was k=2, d=1, which is impossible for odd n. ✓

But wait — the original weaker bound also fails for k=2, d≥2 potentially:
For k=2, d=2 (two distinct parts): G = gcd(c_1, c_2) ≥ 1. v_2(z) = 0 (all multiplicities 1).
  Need: gcd(c_1,c_2) > (2-1)/2 = 1/2. Since gcd ≥ 1 > 1/2. ✓

So actually the proof works! Let me formalize this.

For ALL λ ≠ (n) with odd parts summing to odd n, k ≥ 2:

CASE 1: d ≥ 2 (at least 2 distinct part sizes).
  G ≥ C(k,2) = k(k-1)/2 (each gcd ≥ 1).
  v_2(z) ≤ k - d ≤ k - 2.
  G - v_2(z) ≥ k(k-1)/2 - k + 2 = (k² - 3k + 4)/2 = ((k-1)(k-2)+2)/2.
  For k ≥ 2: (k-1)(k-2) ≥ 0, so G - v_2(z) ≥ 1 > (k-1)/2 ... wait, not necessarily.

  For k=2: G - v_2(z) ≥ 1 - 0 = 1. (k-1)/2 = 1/2. 1 > 1/2. ✓
  For k=3: G - v_2(z) ≥ 3 - 1 = 2. (k-1)/2 = 1. 2 > 1. ✓
  For k���4: G - v_2(z) ≥ k(k-1)/2 - k + 2 vs (k-1)/2.
    Need: k(k-1)/2 - k + 2 > (k-1)/2.
    i.e., (k-1)(k-1)/2 - k + 2 + 1/2 > 0... just check:
    k(k-1)/2 - k + 2 - (k-1)/2 = (k-1)(k-1)/2 - k + 2 + 1/2
    = ((k-1)²)/2 - k + 5/2 = (k²-2k+1)/2 - k + 5/2 = (k²-4k+6)/2
    = ((k-2)² + 2)/2 ≥ 1 > 0. ✓

CASE 2: d = 1 (all parts equal: λ = (c^k), kc = n).
  Since n is odd and all c_i = c are odd: k must be odd too (odd × odd = odd).
  k = 1 would give λ = (n), excluded. So k ≥ 3.
  G = C(k,2) × c = k(k-1)c/2.
  v_2(z) = k × v_2(c) + v_2(k!) = 0 + v_2(k!) (since c odd).
  v_2(k!) ≤ k - 1.
  G - v_2(z) ��� k(k-1)c/2 - k + 1.
  For c ≥ 1: ≥ k(k-1)/2 - k + 1 = (k² - 3k + 2)/2 = (k-1)(k-2)/2.
  (k-1)/2 comparison: (k-1)(k-2)/2 > (k-1)/2 iff k-2 > 1 iff k �� 4. ✓

  k = 3: G - v_2(z) ≥ 3c - 1 (computed above). (k-1)/2 = 1. Need 3c - 1 > 1.
  c ≥ 3 (since c odd, kc = n ≥ 9): 3(3) - 1 = 8 > 1. ✓
  c = 1 would give n = 3, λ = (1,1,1): G = 3, v_2(z) = v_2(3!) = 1.
  G - v_2(z) = 2 > 1 = (k-1)/2. ���

So the proof is COMPLETE for all cases. Let me also verify the tightest partition
is always (n-2, 1, 1).
"""

from math import factorial, gcd

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: v += 1; n //= 2
    return v

def odd_partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0:
        yield ()
        return
    start = min(n, max_part)
    if start % 2 == 0: start -= 1
    for k in range(start, 0, -2):
        for rest in odd_partitions(n - k, k):
            yield (k,) + rest

def z_lambda(lam):
    from collections import Counter
    counts = Counter(lam)
    z = 1
    for c in lam:
        z *= c
    for m in counts.values():
        z *= factorial(m)
    return z

def G_lambda(lam):
    k = len(lam)
    return sum(gcd(lam[i], lam[j]) for i in range(k) for j in range(i+1, k))

print("=" * 70)
print("PROOF VERIFICATION: G(λ) - v_2(z_��) > (k-1)/2")
print("Finding the minimum gap for each odd n")
print("=" * 70)
print()

print(f"{'n':>4} {'worst λ':>30} {'G':>5} {'v_2(z)':>7} {'k':>3} {'G-v_2(z)':>9} {'(k-1)/2':>8} {'gap':>5}")
print("-" * 75)

for n in range(3, 30, 2):
    worst = None
    worst_gap = float('inf')
    for lam in odd_partitions(n):
        if lam == (n,): continue
        k = len(lam)
        G = G_lambda(lam)
        vz = v2(z_lambda(lam))
        lhs = G - vz
        rhs = (k - 1) / 2
        gap = lhs - rhs
        if gap < worst_gap:
            worst_gap = gap
            worst = lam
            worst_G = G
            worst_vz = vz
            worst_k = k

    print(f"{n:>4} {str(worst):>30} {worst_G:>5} {worst_vz:>7} {worst_k:>3} {worst_G - worst_vz:>9} {(worst_k-1)/2:>8.1f} {worst_gap:>5.1f}")

print()
print("KEY OBSERVATION: The worst case is ALWAYS λ = (n-2, 1, 1) with gap = 1.")
print()
print("For λ = (n-2, 1, 1):")
print("  G = gcd(n-2,1) + gcd(n-2,1) + gcd(1,1) = 3")
print("  z = (n-2)·1·1·1!·2! = 2(n-2)")
print("  v_2(z) = 1 + v_2(n-2) = 1 (since n-2 odd for odd n)")
print("  k = 3, (k-1)/2 = 1")
print("  Gap = 3 - 1 - 1 = 1 > 0 ✓")
print()
print("ALGEBRAIC PROOF (for all odd n ≥ 3):")
print()
print("Let λ = (c_1,...,c_k) be an all-odd partition of odd n, with k ≥ 2.")
print()
print("CASE 1: d ≥ 2 distinct part sizes.")
print("  G ≥ C(k,2) (each gcd ≥ 1)")
print("  v_2(z) ≤ k - d ≤ k - 2  (since v_2(a!) ≤ a-1 and Σ(a_j-1) = k-d)")
print("  LHS = G - v_2(z) ≥ k(k-1)/2 - k + 2")
print("  RHS = (k-1)/2")
print("  Difference ≥ ((k-2)² + 2)/2 ≥ 1 > 0.  ✓")
print()
print("CASE 2: d = 1, all parts equal: λ = (c^k), kc = n.")
print("  n odd, c odd ⟹ k odd. Since λ ≠ (n): k ≥ 3.")
print("  G = C(k,2)·c = k(k-1)c/2")
print("  v_2(z) = v_2(k!) (since c is odd)")
print("  For k=3: G-v_2(z) = 3c-1, (k-1)/2 = 1. Gap = 3c-2 ≥ 1 (c ≥ 1). ✓")
print("  For k≥5: G-v_2(z) ≥ k(k-1)/2-k+1 = (k-1)(k-2)/2 > (k-1)/2. ✓")
print()
print("QED. The inequality holds for ALL odd n ≥ 3 and all λ ≠ (n). □")
