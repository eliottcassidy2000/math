#!/usr/bin/env python3
"""
chi_per_sum_formula.py — opus-2026-03-12-S69

The sum Σ_T chi_per(T) over all 2^m circulant orientations:
  p=5:  Σ = 0
  p=7:  Σ = 2
  p=11: Σ = 22
  p=13: Σ = -128

Question: Is there a closed-form formula for this sum?

Key insight: chi_per = Σ(-1)^k Ω_k where Ω_k is per-eigenspace.
And Σ_T Ω_k(T) is the sum over all orientations of the per-eigenspace
path dimension at degree k.

So Σ_T chi_per(T) = Σ_k (-1)^k Σ_T Ω_k(T).

For fixed k, Σ_T Ω_k(T) counts the total number of paths of length k
in eigenspace 1 (or any fixed k≠0) summed over ALL tournament orientations.

This should have a combinatorial interpretation.
"""

import sys
import numpy as np
sys.path.insert(0, '04-computation')

# From computed data:
# p=5 (m=2): avg Ω = [1.00, 2.00, 2.00, 2.00, 1.00], total 4 orientations
# p=7 (m=3): avg Ω = [1.00, 3.00, 6.00, 10.50, 12.75, 8.25, 2.25], total 8
# p=11 (m=5): avg Ω = [1.00, 5.00, 20.00, 71.25, 208.125, 463.4375, ...], total 32

# Note: Ω_0 = 1 always, Ω_1 = m always, Ω_2 = m(m-1) always.
# So avg Ω_0 = 1, avg Ω_1 = m, avg Ω_2 = m(m-1) are all independent of orientation.

# Total Σ_T Ω_m from data:
# p=5: [4, 8, 8, 8, 4]
# p=7: [8, 24, 48, 84, 102, 66, 18]
# p=11: [32, 160, 640, 2280, 6660, 14830, 23970, 27830, 22370, 10990, 2440]

# Let's verify the universal parts:
# Σ_T Ω_0 = 2^m * 1 = 2^m
# Σ_T Ω_1 = 2^m * m
# Σ_T Ω_2 = 2^m * m(m-1)
# These check: p=7: 8*1=8, 8*3=24, 8*6=48 ✓

# Σ chi_per = Σ_k (-1)^k · Σ_T Ω_k(T) = Σ_k (-1)^k · (2^m · avg_Ω_k)

# For p=7: 8*(1 - 3 + 6 - 10.5 + 12.75 - 8.25 + 2.25)
#         = 8*(1 - 3 + 6 - 10.5 + 12.75 - 8.25 + 2.25)
#         = 8 * 0.25 = 2 ✓

# So the question reduces to: what is avg Ω_k as a function of k, m, p?

# The universal values are:
# avg Ω_0 = 1
# avg Ω_1 = m
# avg Ω_2 = m(m-1)
# avg Ω_3 = ? (this is where it gets interesting)

# At p=7: avg Ω_3 = 84/8 = 10.5
# At p=11: avg Ω_3 = 2280/32 = 71.25

# Universal Ω_3 (for non-Interval orbits) at p=7 is 9, at p=11 is 70.
# Interval Ω_3 at p=7 is 11, at p=11 is 74.

# At p=7: avg = (6*11 + 2*9)/8 = (66+18)/8 = 84/8 = 10.5 ✓
# At p=11: need to weight by orbit sizes.

# Let's look at this from a different angle.
# An orientation is a function σ: {1,...,m} → {0,1} where σ(i)=0 means d_i ∈ S,
# σ(i)=1 means p-d_i ∈ S. There are 2^m such functions.
# Ω_k(σ) depends on the set S = {d_i if σ(i)=0, p-d_i if σ(i)=1}.
# So avg Ω_k = (1/2^m) Σ_σ Ω_k(σ).

# For k ≤ 2: Ω_k is independent of σ (universal), so avg = universal value.
# For k = 3: it depends on σ through the Fourier spectrum.

# Let's try to understand avg Ω_3.
# m(m-1) = Ω_2 = number of length-2 path shapes.
# A length-3 path shape is (d_1, d_2, d_3) with d_i ∈ S and all partial sums distinct.
# For a fixed shape to exist, we need d_1, d_2, d_3 ∈ S AND the three partial sums
# s_1=d_1, s_2=d_1+d_2, s_3=d_1+d_2+d_3 are all distinct and nonzero.
# Also s_1 ≠ s_2 (automatic since d_2≠0), s_1 ≠ s_3, s_2 ≠ s_3, etc.

# Actually, the constraint is that all of 0, d_1, d_1+d_2, d_1+d_2+d_3 are distinct mod p.
# This means:
# d_1 ≠ 0 (always, since d_i ∈ {1,...,p-1})
# d_2 ≠ 0
# d_1 + d_2 ≠ 0 (i.e., d_2 ≠ p-d_1)
# d_1 + d_2 ≠ d_1 (i.e., d_2 ≠ 0, redundant)
# d_3 ≠ 0
# d_1+d_2+d_3 ≠ 0, d_1, d_1+d_2

# BUT we also need this to be an ALLOWED path: all faces must be valid too.
# This is the Ω vs A distinction.

# Hmm, this is getting complicated. Let me just compute the sum Σ chi_per
# for the data we have and look for patterns.

print("="*80)
print("Σ chi_per OVER ALL ORIENTATIONS")
print("="*80)

# From data:
data = {
    5: {'m': 2, 'sigma': 0},
    7: {'m': 3, 'sigma': 2},
    11: {'m': 5, 'sigma': 22},
    13: {'m': 6, 'sigma': -128},
}

print(f"\n{'p':>4s} {'m':>3s} {'2^m':>5s} {'Σ chi':>8s} {'avg':>10s} {'Σ/2':>8s}")
for p in sorted(data):
    d = data[p]
    print(f"{p:4d} {d['m']:3d} {2**d['m']:5d} {d['sigma']:8d} {d['sigma']/(2**d['m']):10.4f} {d['sigma']//2:8d}")

# Pattern analysis:
# Σ = 0, 2, 22, -128
# Σ/2 = 0, 1, 11, -64
# 0, 1, 11, -64

# Is Σ/2 related to a Legendre sum or Jacobi symbol?
# At p=7: 1 = (7-1)/6? Nope, 6/6=1 ✓ but coincidence.
# At p=11: 11 = p. Interesting!
# At p=13: -64 = -4^3 = -(p-1)/3 * ... hmm.

# Let me check: at p=7 (≡3 mod 4), Σ/2 = 1.
# At p=11 (≡3 mod 4), Σ/2 = 11 = p.
# At p=5 (≡1 mod 4), Σ/2 = 0.
# At p=13 (≡1 mod 4), Σ/2 = -64.

# p ≡ 3 mod 4: Σ/2 = 0, 1, 11 (p=3,7,11)
# Actually we don't have p=3. Let me compute.
# p=3: m=1, 2 orientations: S={1} and S={2}. Both are the same tournament (unique T3).
# For T3: Ω = [1, 1, 1] (path complex on 3 vertices), chi_per = 1.
# Wait, for p=3: eigenspace decomposition... hmm.
# Actually at p=3, m=1. All circulant tournaments on 3 vertices are the cyclic tournament.
# chi_per = Σ(-1)^k Ω_k where Ω_0=1, Ω_1=1, Ω_2=1 (if the full triangle is a path).
# Actually at n=3: the regular tournament is a 3-cycle. Paths: 0→1→2 and reverse, etc.
# Per-eigenspace: Ω_0=1, Ω_1=1, Ω_2=1 (the transitive closure thing).
# chi_per = 1 - 1 + 1 = 1.
# Σ_T chi_per = 2 * 1 = 2 (both orientations give same chi).

# Updated:
# p=3: Σ = 2, m=1, 2^m=2
# p=5: Σ = 0, m=2, 2^m=4
# p=7: Σ = 2, m=3, 2^m=8
# p=11: Σ = 22, m=5, 2^m=32
# p=13: Σ = -128, m=6, 2^m=64

# Σ: 2, 0, 2, 22, -128
# For p ≡ 3 mod 4 (p=3,7,11): Σ = 2, 2, 22
# For p ≡ 1 mod 4 (p=5,13): Σ = 0, -128

# Σ/2 for p ≡ 3: 1, 1, 11
# Hmm. 1, 1, 11. The pattern 1, 1, 11 doesn't match p or m.

# Actually, maybe the average chi = Σ/2^m is more informative:
# p=3: 1
# p=5: 0
# p=7: 0.25
# p=11: 0.6875 = 11/16
# p=13: -2

# Sequence: 1, 0, 1/4, 11/16, -2

# These don't follow an obvious pattern. But note:
# For p=3: 1 = 1
# For p=7: 1/4 = 1/(p-3) = 1/4
# For p=11: 11/16 = p/16... 16 = 2^4 = 2^{m-1}. 11/2^4.

# p=3 (m=1): avg = 1 = (p)/2^{2m-2} = 3/1? No. 1 = 2/2^1 = Σ/2^m.
# p=7 (m=3): avg = 1/4 = 2/8 = Σ/2^m.
# p=11 (m=5): avg = 11/16 = 22/32 = Σ/2^m.
# These are all just Σ/2^m by definition.

# Let me look at the numerators Σ: 2, 0, 2, 22, -128
# 22 = 2 * 11
# -128 = -2^7
# Does Σ relate to the Legendre symbol or class number?

# Legendre symbol pattern:
# For each orientation S, the orientation is determined by choosing
# one from each pair (d, p-d). The Legendre symbol (d/p) classifies d as QR/NQR.
# For p ≡ 3 mod 4: (-1/p) = -1, so each pair (d, p-d) has one QR and one NQR.
# For p ≡ 1 mod 4: (-1/p) = 1, so each pair has same QR/NQR status.

print(f"\n{'='*80}")
print(f"ARITHMETIC ANALYSIS OF Σ chi_per")
print(f"{'='*80}")

for p in [3, 5, 7, 11, 13]:
    m = (p-1)//2
    sigma = data.get(p, {}).get('sigma', None)
    if sigma is None:
        if p == 3:
            sigma = 2
        else:
            continue

    # Legendre symbol structure
    from math import gcd
    qr = set(pow(g, 2, p) for g in range(1, p))
    h = 0  # class number of Q(sqrt(-p))
    # Approximate: h ≈ sum of Legendre symbols
    L_sum = sum(pow(a, (p-1)//2, p) if pow(a, (p-1)//2, p) != p-1 else -1 for a in range(1, p))

    # Gauss sum g^2 = (-1)^{(p-1)/2} * p
    gauss_sq = (-1)**((p-1)//2) * p  # = -p for p≡3, +p for p≡1

    print(f"\n  p={p}: m={m}, Σ={sigma}, 2^m={2**m}, avg={sigma/(2**m):.4f}")
    print(f"    p mod 4 = {p % 4}")
    print(f"    Gauss sum squared = {gauss_sq}")
    print(f"    L(1, χ_p) Legendre sum = {L_sum}")

    # Also: total number of QR in {1,...,m}
    qr_in_half = sum(1 for d in range(1, m+1) if d in qr)
    print(f"    QR in {{1,...,{m}}}: {qr_in_half} of {m}")

# Check: is Σ chi_per = 2 * (Σ_{k≠0} χ_per(eigenspace k, Paley))?
# Only Paley has nontrivial stabilizer...

# Let me think about this differently.
# chi_per(T) = Σ(-1)^k dim Ω_k(T) per eigenspace.
# Σ_T chi_per(T) = Σ_T Σ(-1)^k Ω_k(T) = Σ_k (-1)^k Σ_T Ω_k(T)
# = Σ_k (-1)^k · (total number of allowed k-paths, summed over orientations, per eigenspace)

# Per eigenspace, the total Ω_k summed over all 2^m orientations.
# By translation invariance, Σ_T Ω_k(T) = (1/p) Σ_T (number of valid k+1 tuples in T)
# = (1/p) (number of tuples valid in T, summed over T)

# For a k-path (v_0,...,v_k) with distinct vertices, the orientation T is "compatible"
# if v_i → v_{i+1} for all i (meaning v_{i+1}-v_i mod p is in S).
# A tuple is compatible with T iff all differences d_i = v_{i+1}-v_i satisfy d_i ∈ S_T.

# The number of orientations for which a given tuple is valid equals the number of
# orientations containing all required differences. Since the required differences
# are d_1,...,d_k ∈ S, and each d_i is determined by the pair it belongs to,
# the count depends on how many DISTINCT pairs are needed.

# If the k differences d_1,...,d_k involve r distinct pairs, then 2^{m-r} orientations
# make this tuple valid (the r required pairs are forced, the rest are free).

# So: Σ_T (# valid k-paths in T) = Σ_{valid tuples with distinct verts} 2^{m - #pairs(tuple)}

# This gives a formula for Σ_T Ω_k(T) · p in terms of counting tuples by pair count.

print(f"\n{'='*80}")
print(f"PAIR COUNT ANALYSIS")
print(f"{'='*80}")

for p in [5, 7]:
    m = (p-1)//2
    print(f"\np = {p}, m = {m}")

    # Enumerate all valid k-path shapes (starting at 0)
    def valid_shapes(p, length):
        """All length-step paths 0 → d_1 → d_1+d_2 → ... with distinct vertices
           using differences from {1,...,p-1}."""
        results = []
        def backtrack(diffs, partial_sums):
            if len(diffs) == length:
                results.append(tuple(diffs))
                return
            for d in range(1, p):
                new_sum = (partial_sums[-1] + d) % p
                if new_sum not in set(partial_sums):
                    backtrack(diffs + [d], partial_sums + [new_sum])
        backtrack([], [0])
        return results

    def pair_of(d, p):
        """Which pair does d belong to? Returns the pair index."""
        m = (p-1)//2
        return min(d, p-d)

    for k in range(1, p):
        shapes = valid_shapes(p, k)
        if not shapes:
            continue

        pair_count_hist = {}
        total_weighted = 0
        for shape in shapes:
            pairs_used = set(pair_of(d, p) for d in shape)
            r = len(pairs_used)
            weight = 2**(m - r)
            total_weighted += weight
            pair_count_hist[r] = pair_count_hist.get(r, 0) + 1

        # This gives Σ_T (# k-paths starting at 0 in T)
        # Ω_k total over all T = total_weighted
        # Per eigenspace: total_weighted / p (approximately...)
        # Actually Ω_k counts ALLOWED chains (elements of Ω), not all paths.
        # Hmm, we need Ω_k specifically, not A_k.

        print(f"  k={k}: {len(shapes)} shapes, pair_count histogram = {dict(sorted(pair_count_hist.items()))}")
        print(f"    Σ_T A_k(T)/p ≈ {total_weighted/p:.1f} (sum of 2^{m-r} over shapes)")

print("\nDONE.")
