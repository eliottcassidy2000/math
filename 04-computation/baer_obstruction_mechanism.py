"""
baer_obstruction_mechanism.py
opus-2026-03-14-S71l

Deep investigation of WHY the forbidden values {7, 21, 63} are forbidden,
using the Baer subplane hierarchy and the 7*3^k pattern.

Key findings from S89: H=63 is PERMANENTLY FORBIDDEN at n=7.
Pattern: forbidden = {7*3^0, 7*3^1, 7*3^2} but 7*3^3 = 189 is ACHIEVABLE.

This script investigates:
1. The Baer recursive obstruction mechanism
2. Why the 7*3^k pattern terminates
3. The forbidden value 7 = |2-omega|^2 interpretation
4. The connection to PG(2,4) Baer partition
5. Whether H=63 is forbidden at ALL n, or just n=7
"""

import numpy as np
from math import gcd, factorial, comb
from itertools import permutations, combinations
from collections import Counter

print("=" * 70)
print("BAER OBSTRUCTION MECHANISM")
print("opus-2026-03-14-S71l")
print("=" * 70)

# =====================================================================
# PART 1: THE FORBIDDEN VALUE HIERARCHY
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: THE FORBIDDEN VALUE HIERARCHY {7, 21, 63}")
print("=" * 70)

print("""
  Known forbidden H values (permanently forbidden, all n):
    7  = Phi_3(2)     = 7 * 3^0  (PROVED)
    21 = Phi_3(4)     = 7 * 3^1  (PROVED)

  Newly discovered (S89, exact n=7 enumeration):
    63 = 2^6-1 = 7*9 = 7 * 3^2  (FORBIDDEN at n=7)

  Achievable:
    189 = 7 * 3^3  (achievable at n=7, regular tournament)
    35  = 7 * 5    (achievable at n=7)

  The pattern 7*3^k is forbidden for k = 0, 1, 2.

  CYCLOTOMIC DECOMPOSITION:
    7  = Phi_3(2)
    21 = Phi_3(2) * Phi_6(2) = 7 * 3
    63 = Phi_3(2) * Phi_2(2) * Phi_6(2) = 7 * 3 * 3

  Wait: 63 = Phi_3(2) * Phi_2(2)^2? No.
  63 = 2^6 - 1 = (2^3-1)(2^3+1) = 7 * 9 = 7 * 3^2

  Actually: 2^6 - 1 = prod Phi_d(2) for d | 6, d > 0
  = Phi_1(2) * Phi_2(2) * Phi_3(2) * Phi_6(2)
  = 1 * 3 * 7 * 3 = 63. Yes!

  But Phi_1(2) = 1 (trivial), so:
  63 = Phi_2(2) * Phi_3(2) * Phi_6(2) = 3 * 7 * 3
""")

# Verify
from sympy import factorint, cyclotomic_poly, Symbol
x = Symbol('x')

for val in [7, 21, 63, 189, 35]:
    factors = dict(factorint(val))
    print(f"  {val} = {factors}")

print()

# Cyclotomic values at 2
print("  Phi_d(2) for d = 1..12:")
for d in range(1, 13):
    cp = cyclotomic_poly(d, x)
    val = int(cp.subs(x, 2))
    print(f"    Phi_{d:2d}(2) = {val}")

# =====================================================================
# PART 2: THE RECURSIVE BAER OBSTRUCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: THE RECURSIVE BAER OBSTRUCTION")
print("=" * 70)

print("""
  THEOREM (Baer Partition): PG(2,q^2) decomposes into q^2-q+1 = Phi_6(q)
  disjoint copies of PG(2,q) (Baer subplanes) plus q^2+1 "free points".

  Wait, that's not right. The Baer partition is:
  PG(2,q^2) can be partitioned into Phi_3(q) = q^2+q+1 lines of q^2+1 points...
  No. Let me be precise.

  The BAER PARTITION theorem:
  PG(2,q^2) has q^4+q^2+1 = Phi_3(q^2) points.
  It can be decomposed as a DISJOINT UNION of Phi_6(q) = q^2-q+1
  copies of PG(2,q), ONLY when q is a prime power.

  Actually this isn't standard. What IS standard:
  - PG(2,q^2) contains copies of PG(2,q) as Baer subplanes
  - Through each point of PG(2,q^2), there pass many Baer subplanes
  - The number of Baer subplanes = |PGL(3,q^2)| / |PGL(3,q)| * ...

  For the SPECIFIC case q=2:
  PG(2,4) = 21 points, and 21 = 3*7.
  There are EXACTLY 3 disjoint Baer subplanes (each = PG(2,2) = 7 pts).
  This is because 21/7 = 3.

  For q=4:
  PG(2,16) = 273 points. Phi_6(4) = 13.
  So 273 = 21 * 13 = PG(2,4) * 13.
  There are 13 disjoint copies of PG(2,4)?

  Let me verify: 273/21 = 13 = Phi_6(4). YES!

  THE RECURSIVE STRUCTURE:
  Level 0: 7 points (PG(2,2))
  Level 1: 21 = 7 * 3 points (PG(2,4) = 3 copies of PG(2,2))
  Level 2: 273 = 21 * 13 points (PG(2,16) = 13 copies of PG(2,4))
  Level 3: 65793 = 273 * 241 points (PG(2,256) = 241 copies of PG(2,16))

  THE OBSTRUCTION:
  At level 0: 7 is forbidden (Fano obstruction).
  At level 1: 21 = 3 * 7. If H=21, the HPs must "see" the 3 Baer subplanes.
              Each subplane contributes 7 to the HP count (via restriction?).
              But 7 is forbidden, so the contribution is impossible.
  At level 2: 273 is ACHIEVABLE (proved in S71i).
              The obstruction doesn't propagate to level 2!

  WHY DOES THE OBSTRUCTION STOP AT LEVEL 2?
  Because at level 2, the cofactor 13 is NOT a power of 3.
  The Baer partition of PG(2,16) has 13 copies of PG(2,4),
  and 13 is prime (and NOT 3).
  So the recursive "3 copies of a forbidden value" argument doesn't apply.

  The obstruction requires the cofactor to be a power of 3 (= Phi_6(2)):
  - 21/7 = 3 = Phi_6(2): obstruction applies
  - 63/21 = 3: obstruction applies (63 = 21 * 3)
  - 189/63 = 3: obstruction SHOULD apply, but 189 IS achievable!

  So even 3-powers don't always give obstruction.

  RESOLUTION: The obstruction is about the TOURNAMENT SIZE, not just the value.
  At n=7, max H = 189. The tournament has only 7 vertices and 21 arcs.
  A tournament with H=63 would need a very specific structure,
  and that structure is incompatible with having exactly 7 vertices
  and the Fano plane obstruction.

  At larger n (n >= 9?), H=63 might become achievable!
  This is an OPEN QUESTION: Is H=63 permanently forbidden or only for n <= 7?
""")

# =====================================================================
# PART 3: COMPUTATIONAL CHECK — IS 63 FORBIDDEN AT n=8?
# =====================================================================
print("=" * 70)
print("PART 3: IS H=63 FORBIDDEN AT n=8? (Monte Carlo)")
print("=" * 70)

def compute_H_dp(adj, n):
    """Count Hamiltonian paths via DP (O(2^n * n))."""
    # dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full_mask = (1 << n) - 1
    return sum(dp[full_mask][v] for v in range(n))

n = 8
print(f"  Sampling random tournaments at n={n} to check if H=63 appears...")
import random

h_counts = Counter()
n_samples = 100000
found_63 = False

for trial in range(n_samples):
    # Generate random tournament
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

    h = compute_H_dp(adj, n)
    h_counts[h] += 1
    if h == 63:
        found_63 = True
        if h_counts[63] <= 3:
            print(f"    Found H=63 at trial {trial}!")

    if (trial+1) % 10000 == 0:
        print(f"    {trial+1}/{n_samples} done, {len(h_counts)} distinct H values, H=63 count: {h_counts.get(63, 0)}")

print(f"\n  H=63 found: {found_63} ({h_counts.get(63, 0)} times)")
print(f"  H=7 found: {7 in h_counts} ({h_counts.get(7, 0)} times)")
print(f"  H=21 found: {21 in h_counts} ({h_counts.get(21, 0)} times)")

# Also check other values of interest
for h in [7, 21, 35, 39, 63, 73, 91, 105, 119, 147, 189]:
    print(f"    H={h:3d}: {h_counts.get(h, 0)} times")

# =====================================================================
# PART 4: THE |2-omega|^2 INTERPRETATION
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: FORBIDDEN VALUES AS SQUARED DISTANCES")
print("=" * 70)

omega = np.exp(2j * np.pi / 3)

print("""
  Phi_3(x) = |x - omega|^2 where omega = e^(2*pi*i/3).

  The FORBIDDEN values are |2^k - omega|^2 = Phi_3(2^k):
  k=1: |2 - omega|^2 = 7   (FORBIDDEN)
  k=2: |4 - omega|^2 = 21  (FORBIDDEN)
  k=3: |8 - omega|^2 = 73  (unknown status)
  k=4: |16 - omega|^2 = 273 (ACHIEVABLE)

  But 63 is NOT of the form Phi_3(2^k) = (2^k)^2 + 2^k + 1!
  63 = 2^6 - 1, not (2^k)^2 + 2^k + 1 for any k.

  So the |2^k - omega|^2 pattern captures {7, 21} but NOT {63}.
  The forbidden value 63 has a DIFFERENT origin.

  63 = R_6(2) = (2^6 - 1)/(2 - 1) = sum_{k=0}^{5} 2^k = 111111 in base 2.

  It's a REPUNIT in base 2 (length 6 = the tournament period!).

  The forbidden values so far:
    7  = 111 in base 2 (repunit length 3)
    21 = 10101 in base 2 (NOT a repunit)
    63 = 111111 in base 2 (repunit length 6)

  Hmm, 21 is NOT a repunit in base 2. But:
    21 = 111 in base 4 (repunit length 3 in base 4!)

  So: 7 = 111_2, 21 = 111_4, 63 = 111111_2 = 1111_4?
  63 in base 4: 63 = 3*16 + 3*4 + 3 = 333_4! A REPDIGIT (not repunit).

  Actually:
    7  = 111_2 (repunit base 2)
    21 = 111_4 (repunit base 4)
    63 = 111111_2 (repunit base 2, length 6)
       = 333_4 (repdigit 3 in base 4)
       = 77_8 (repdigit 7 in base 8)
       = 111_8? No, 1*64+1*8+1 = 73 not 63.

  63 = 7*9 = 7 * 3^2.
  In base 8: 63 = 77_8 = 7*8 + 7 = 7*(8+1) = 7*9.
  So 63 = 77_8 = repdigit 7 in base 8!

  And 7 = 7_8 = single digit 7 in base 8.
  And 21 = 25_8 = NOT a repdigit in base 8.

  Different pattern in base 6 (the period base):
    7 = 11_6 (repdigit 1)
    21 = 33_6 (repdigit 3)
    63 = 143_6 (NOT a repdigit)

  So the base-6 repdigit pattern {7, 21} does NOT extend to 63.
""")

# Verify
for h in [7, 21, 63, 189]:
    # Base 2
    b2 = bin(h)[2:]
    # Base 4
    b4 = ""
    t = h
    while t > 0:
        b4 = str(t % 4) + b4
        t //= 4
    # Base 6
    b6 = ""
    t = h
    while t > 0:
        b6 = str(t % 6) + b6
        t //= 6
    # Base 8
    b8 = oct(h)[2:]
    print(f"  {h:3d}: base2={b2}, base4={b4}, base6={b6}, base8={b8}")

# =====================================================================
# PART 5: THE DEEP FACTORIZATION PATTERN
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: DEEP FACTORIZATION PATTERN")
print("=" * 70)

print("""
  Let's look at ALL forbidden values at each n and their factorizations:

  n=3: achievable = {1, 3}. Max = 3. No forbidden (only 2 values possible).
  n=4: achievable = {1, 3}. Max = 3. No forbidden (trivially).
  n=5: achievable = {1, 3, 5, 9, 13}. Max = 13. Forbidden odd in [1,13]: {7, 11}
       Wait: 11 is achievable at n=5? Let me recheck.
""")

# Compute exact H-spectrum for n=5
def all_tournaments(n):
    m = comb(n, 2)
    arcs = [(i,j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(arcs):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

print("  Computing exact H-spectra for n=3,4,5,6...")
for n in range(3, 7):
    h_vals = Counter()
    for adj in all_tournaments(n):
        h = compute_H_dp(adj, n)
        h_vals[h] += 1

    achievable = sorted(h_vals.keys())
    max_h = max(achievable)
    odd_range = list(range(1, max_h+1, 2))
    forbidden = [h for h in odd_range if h not in achievable]

    print(f"\n  n={n}: achievable = {achievable}")
    print(f"    max H = {max_h}")
    print(f"    forbidden odd in [1,{max_h}]: {forbidden}")

    for f in forbidden:
        factors = dict(factorint(f))
        divisible_by_7 = f % 7 == 0
        print(f"      {f} = {factors}, div by 7: {divisible_by_7}")

print("""
  KEY OBSERVATION:
  At n=5: H-spectrum = {1, 3, 5, 9, 13}. Forbidden odd: {7, 11}.
  At n=6: H-spectrum = {1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}.
          Forbidden odd: {7, 21, 35, 39}.

  So 11 is forbidden at n=5 but ACHIEVABLE at n=6!
  35 and 39 are forbidden at n=6 but ACHIEVABLE at n=7!
  Only 7 and 21 remain forbidden at n=7.
  And now 63 joins the permanently-forbidden set at n=7.

  The PERMANENTLY forbidden set (forbidden at ALL n):
    7 = Phi_3(2)     PROVED
    21 = Phi_3(4)    PROVED
    63 = 7*9         STATUS: forbidden at n=7, unknown for n >= 8
""")

print("\n" + "=" * 70)
print("DONE — BAER OBSTRUCTION MECHANISM")
print("=" * 70)
