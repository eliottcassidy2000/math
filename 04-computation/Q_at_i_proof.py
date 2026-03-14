#!/usr/bin/env python3
"""
Prove Q_n(i) = i * N/2 where N = 2^m, m = C(n,2).
opus-2026-03-14-S84

Q_n(q) = Σ_{T} q^{H(T)}
Q_n(i) = Σ_{T} i^{H(T)}

Since H(T) is always odd, i^{H(T)} = i^{2k+1} = i * (-1)^k
So Q_n(i) = i * Σ_T (-1)^{(H(T)-1)/2}

The claim Q_n(i) = i * N/2 means:
Σ_T (-1)^{(H(T)-1)/2} = N/2

Equivalently: #{T : H(T) ≡ 1 mod 4} - #{T : H(T) ≡ 3 mod 4} = N/2

Since H is always odd, H mod 4 is either 1 or 3.
The claim says: #{H ≡ 1} - #{H ≡ 3} = N/2
i.e., #{H ≡ 1} = 3N/4, #{H ≡ 3} = N/4
or equivalently: 75% of tournaments have H ≡ 1 mod 4!

Let's verify this and understand why.
"""

from itertools import permutations
from collections import Counter
import math
import sys

def compute_all_H(n):
    m = n * (n - 1) // 2
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    all_perms = list(permutations(range(n)))
    H_values = []
    for bits in range(N):
        if bits % 10000 == 0 and N > 10000:
            print(f"  n={n}: {bits}/{N}", file=sys.stderr)
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H = sum(1 for p in all_perms if all(adj[p[i]][p[i+1]] == 1 for i in range(n-1)))
        H_values.append(H)
    return H_values

# ============================================================
# Part 1: Verify Q_n(i) = i*N/2
# ============================================================
print("=" * 70)
print("PART 1: VERIFY Q_n(i) = i * N/2")
print("=" * 70)

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    m = n * (n - 1) // 2

    # Compute Q(i)
    val = sum(1j**h for h in H_vals)
    expected = 1j * N / 2

    print(f"\nn={n}: N={N}")
    print(f"  Q(i) = {val.real:.6f} + {val.imag:.6f}i")
    print(f"  i*N/2 = {expected.real:.6f} + {expected.imag:.6f}i")
    print(f"  Match: {abs(val - expected) < 0.001}")

    # Count H mod 4
    mod4 = Counter(h % 4 for h in H_vals)
    print(f"  H mod 4 distribution: {dict(sorted(mod4.items()))}")
    print(f"  #{'{H≡1}'}={mod4.get(1,0)}, #{'{H≡3}'}={mod4.get(3,0)}")
    print(f"  #{'{H≡1}'} - #{'{H≡3}'} = {mod4.get(1,0) - mod4.get(3,0)}")
    print(f"  N/2 = {N//2}")
    print(f"  #{'{H≡1}'}/N = {mod4.get(1,0)/N:.6f} (should be 3/4 = {3/4})")

# ============================================================
# Part 2: Understanding the 3/4 ratio
# ============================================================
print("\n" + "=" * 70)
print("PART 2: WHY 3/4 OF TOURNAMENTS HAVE H ≡ 1 MOD 4")
print("=" * 70)

# H ≡ 1 mod 4 means H = 1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45
# H ≡ 3 mod 4 means H = 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43

# At n=3: H∈{1,3}. H=1 count=6, H=3 count=2.
# 6 ≡ 1 mod 4: count=6 (H=1: 1≡1)
# 3 ≡ 3 mod 4: count=2 (H=3: 3≡3)
# 6-2 = 4 = 8/2 ✓

# KEY INSIGHT: The transitive tournaments all have H=1 ≡ 1 mod 4.
# There are n! transitive tournaments.
# The cyclic tournaments at n=3 have H=3 ≡ 3 mod 4.
# For this to work: we need #{H≡1} = 3N/4

# Let me check: what fraction of each H value class is ≡ 1 mod 4?
for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    dist = Counter(H_vals)

    print(f"\nn={n}:")
    mod1_total = 0
    mod3_total = 0
    for h in sorted(dist.keys()):
        c = dist[h]
        mod = h % 4
        if mod == 1:
            mod1_total += c
            print(f"  H={h:2d} (≡1 mod 4): {c:5d} tours")
        else:
            mod3_total += c
            print(f"  H={h:2d} (≡3 mod 4): {c:5d} tours")

    print(f"  Total ≡1: {mod1_total}, Total ≡3: {mod3_total}")
    print(f"  Ratio: {mod1_total}/{mod3_total} = {mod1_total/mod3_total:.6f}")
    # Prediction: ratio should be 3.0

# ============================================================
# Part 3: Connection to path reversal
# ============================================================
print("\n" + "=" * 70)
print("PART 3: PATH REVERSAL AND MOD-4 STRUCTURE")
print("=" * 70)

# T^op (complement/reversal) has H(T^op) = H(T).
# This pairs tournaments, preserving H.
# Does this pairing respect mod 4?

# YES! Since H(T) = H(T^op), both T and T^op have the same H mod 4.
# The pair (T, T^op) contributes either (+2) to #{≡1} or (+2) to #{≡3}.
# UNLESS T is self-complementary (T = T^op), in which case it contributes +1.

# At n=3: SC tournaments = C_3 and its reverse = same up to relabeling
# Actually, self-complementary means T is isomorphic to T^op.
# But we're counting LABELED tournaments, so T = T^op means
# the adjacency matrix equals its transpose with 0↔1 swap.

# For LABELED: T = T^op iff for all (i,j), arc(i,j) = reverse(arc(j,i))
# which means adj[i][j] = 1 - adj[j][i] for all i≠j
# But adj[i][j] + adj[j][i] = 1, so 1 - adj[j][i] = adj[i][j].
# This is ALWAYS TRUE! Wait, that means T^op flips ALL arcs.
# T^op has adj'[i][j] = adj[j][i] = 1 - adj[i][j].
# So T = T^op means adj[i][j] = adj[j][i] for all i,j.
# But adj[i][j] + adj[j][i] = 1, so adj[i][j] = 1/2 — impossible!
# Therefore NO labeled tournament equals its complement (for n ≥ 2).

# This means T ≠ T^op always, and tournaments pair up as (T, T^op).
# N is even (N = 2^m), so N/2 pairs.
# Each pair (T, T^op) has the same H, so same H mod 4.
# This doesn't directly explain the 3/4 ratio...

# But wait: T^op is obtained by flipping ALL arcs.
# The map T → T^op is a fixed-point-free involution on {0,...,N-1}.
# It acts as bits → (~bits) & ((1<<m)-1) on the bit encoding.
# And H(T) = H(T^op) by Rédei's theorem (path reversal).

print(f"T and T^op always have the same H (path reversal).")
print(f"No labeled tournament equals its complement (n≥2).")
print(f"So tournaments pair up: N/2 pairs with H(T)=H(T^op).")
print(f"This is consistent with 3/4 ratio but doesn't prove it.")

# ============================================================
# Part 4: Proof attempt via Fourier
# ============================================================
print("\n" + "=" * 70)
print("PART 4: FOURIER PROOF OF Q(i) = iN/2")
print("=" * 70)

# Q_n(q) = Σ_T q^{H(T)}
# Since H(T) is always odd: Q_n(q) = q * R_n(q^2)
# where R_n(u) = Σ_T u^{(H(T)-1)/2}

# Q(i) = i * R(i^2) = i * R(-1)

# R(-1) = Σ_T (-1)^{(H(T)-1)/2}
# The claim is R(-1) = N/2.

# Fourier expansion: H(T) = Σ_S Ĥ(S) χ_S(T)
# where χ_S is the character indexed by arc subset S.

# (-1)^{(H(T)-1)/2} = (-1)^{(H-1)/2}

# Let's compute R(-1) directly using the Fourier coefficients.
# H = Ĥ(∅) + Σ_{|S|≥1} Ĥ(S) χ_S
# Ĥ(∅) = mean(H) = n!/2^{n-1}

# (H-1)/2 = (Ĥ(∅) - 1)/2 + (1/2) Σ Ĥ(S) χ_S

# This is getting complex. Let me just verify the pattern and
# try to understand it combinatorially.

# The claim #{H≡1 mod 4}/N = 3/4 can be written as:
# E[(-1)^{(H-1)/2}] = 1/2

# Let's check: (H-1)/2 mod 2 = ?
# H = 1: (H-1)/2 = 0, (-1)^0 = +1
# H = 3: (H-1)/2 = 1, (-1)^1 = -1
# H = 5: (H-1)/2 = 2, (-1)^2 = +1
# H = 7: (H-1)/2 = 3, (-1)^3 = -1
# So (-1)^{(H-1)/2} = 1 iff H ≡ 1 mod 4, = -1 iff H ≡ 3 mod 4

# The generating function approach:
# Σ_T (-1)^{(H-1)/2} = Σ_h #{H=h} * (-1)^{(h-1)/2}

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    dist = Counter(H_vals)

    R_minus1 = sum(c * (-1)**((h-1)//2) for h, c in dist.items())
    print(f"n={n}: R(-1) = {R_minus1}, N/2 = {N//2}, match = {R_minus1 == N//2}")

# ALL MATCH! R(-1) = N/2 for n=3,4,5,6.
# This means: Σ_T (-1)^{(H(T)-1)/2} = 2^{m-1}

# ============================================================
# Part 5: Higher powers of i — Q(i^k) pattern
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Q(ζ_k) FOR k-TH ROOTS OF UNITY")
print("=" * 70)

import cmath

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)

    print(f"\nn={n}:")
    for k in [3, 4, 5, 6, 7, 8]:
        zeta = cmath.exp(2j * cmath.pi / k)
        val = sum(zeta**h for h in H_vals)
        # Check if val is a nice multiple of something
        val_norm = val / N
        print(f"  Q(ζ_{k})/N = {val_norm.real:+.8f} + {val_norm.imag:+.8f}i, |Q|/N = {abs(val)/N:.8f}")

# ============================================================
# Part 6: Q(i) = iN/2 implies mod-4 structure
# ============================================================
print("\n" + "=" * 70)
print("PART 6: COMPLETE MOD-4 AND MOD-8 STRUCTURE")
print("=" * 70)

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    dist = Counter(H_vals)

    # H mod 8 distribution
    mod8 = Counter()
    for h, c in dist.items():
        mod8[h % 8] += c

    print(f"\nn={n}: H mod 8 distribution:")
    for r in sorted(mod8.keys()):
        pct = 100 * mod8[r] / N
        print(f"  H ≡ {r} mod 8: {mod8[r]:6d} ({pct:.2f}%)")

    # Q(ζ_8) = Σ ζ_8^H
    zeta8 = cmath.exp(2j * cmath.pi / 8)
    val8 = sum(c * zeta8**h for h, c in dist.items())
    print(f"  Q(ζ_8) = {val8.real:.4f} + {val8.imag:.4f}i")
    print(f"  |Q(ζ_8)|/N = {abs(val8)/N:.6f}")

# ============================================================
# Part 7: Connection to Rédei's theorem
# ============================================================
print("\n" + "=" * 70)
print("PART 7: RÉDEI PARITY AND MOD-4 STRUCTURE")
print("=" * 70)

# Rédei: H(T) is odd for all tournaments.
# Our finding: #{H≡1 mod 4}/N = 3/4 for all n.
# This is a STRONGER parity result than Rédei!

# Is there a "mod 4 Rédei theorem"?
# Claim: #{H≡1 mod 4} = 3N/4 = 3 * 2^{m-2}

# Let's check if this extends to mod 8:
# From the mod-8 data, what fractions appear?

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    dist = Counter(H_vals)

    mod8 = Counter()
    for h, c in dist.items():
        mod8[h % 8] += c

    # Only odd residues possible: 1, 3, 5, 7
    print(f"\nn={n}: mod-8 fractions:")
    for r in [1, 3, 5, 7]:
        frac = mod8.get(r, 0) / N
        print(f"  #{'{H≡'+str(r)+' mod 8}'}/N = {frac:.6f}")

# ============================================================
# Part 8: Proof sketch
# ============================================================
print("\n" + "=" * 70)
print("PART 8: PROOF SKETCH — WHY Q(i) = iN/2")
print("=" * 70)

print("""
PROOF SKETCH:

Q_n(i) = Σ_T i^{H(T)}

Since H(T) is always odd (Rédei), write H = 2k+1.
Then i^H = i^{2k+1} = i * (i^2)^k = i * (-1)^k.

So Q_n(i) = i * Σ_T (-1)^{(H(T)-1)/2}.

We need to show: Σ_T (-1)^{(H(T)-1)/2} = N/2.

Define f(T) = (-1)^{(H(T)-1)/2}.
This equals +1 if H ≡ 1 mod 4, -1 if H ≡ 3 mod 4.

Consider the GS involution (graph symmetry): T → T^{GS}.
This preserves the tournament and changes H by a known amount.
Actually, GS preserves H (it just relabels).

Consider instead the arc-flip involution on a SINGLE arc e:
T → T^e (flip arc e).
H(T^e) = H(T) + ΔH_e(T).

By the D_e^2 = 2D_e identity (HYP-1242):
ΔH_e can only be ±2 or ±4 at n=4.
In general, ΔH_e is always even (preserves parity, as expected).

But (-1)^{(H-1)/2} depends on H mod 4, and flipping by ±2
changes H mod 4 by ±2, which flips between 1 and 3 mod 4.

So arc flips with ΔH = ±2 swap the sign of f(T).
Arc flips with ΔH = 0 or ±4 preserve the sign.

For Q(i) = iN/2 to hold, we need:
The average of f over all tournaments is 1/2.

This is a FOURIER STATEMENT: f(T) is a function on {0,1}^m.
Its mean is 1/2 iff the Fourier coefficient f̂(∅) = 1/2.

Since f(T) = (-1)^{(H-1)/2} is a NONLINEAR function of H,
and H is a degree-(n-1) function of arcs,
the proof requires understanding the mod-4 structure of H
in Fourier space.

VERIFIED for n=3,4,5,6. Conjecture: holds for ALL n.
""")

# ============================================================
# Part 9: Check: does #{H≡r mod 2^k}/N have clean fractions?
# ============================================================
print("=" * 70)
print("PART 9: MOD-2^k STRUCTURE")
print("=" * 70)

from fractions import Fraction

for n in [3, 4, 5, 6]:
    H_vals = compute_all_H(n)
    N = len(H_vals)
    dist = Counter(H_vals)

    print(f"\nn={n}:")
    for k in [2, 4, 8, 16]:
        mod_k = Counter()
        for h, c in dist.items():
            mod_k[h % k] += c

        # Check for clean fractions
        fracs = {}
        for r in sorted(mod_k.keys()):
            f = Fraction(mod_k[r], N)
            fracs[r] = f

        # Print only odd residues
        odd_residues = [r for r in sorted(fracs.keys()) if r % 2 == 1]
        print(f"  mod {k}: {' '.join(f'{r}:{fracs[r]}' for r in odd_residues)}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — Q(i) = iN/2")
print("=" * 70)
print("""
THEOREM (verified n=3,...,6, conjectured for all n):
  Q_n(i) = i * 2^{m-1}  where m = C(n,2)

Equivalently:
  #{T : H(T) ≡ 1 mod 4} = 3 * 2^{m-2}
  #{T : H(T) ≡ 3 mod 4} = 2^{m-2}

This says 75% of tournaments have H ≡ 1 mod 4.
This is a REFINEMENT of Rédei's parity theorem.

The ratio 3:1 = KEY2:1 is another (2,3) appearance!
75% = 3/4 = KEY2/KEY1² = KEY2/(KEY1+KEY1)

PROOF IDEA: Use the Fourier expansion of H and the structure
of (-1)^{(H-1)/2} as a Boolean function on arc space.
""")
