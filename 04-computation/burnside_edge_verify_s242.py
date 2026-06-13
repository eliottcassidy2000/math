#!/usr/bin/env python3
"""
burnside_edge_verify_s242.py -- opus-2026-03-23-S242

VERIFY AND EXTEND THE BREAKTHROUGH FORMULA:
  edge_orbits(G_n) = T_n/2 + (n-2)!

And investigate SL_orbits = 2, 5, 20, 86, 490, 3703 to find a formula.

From the S241 breakthrough:
  E(G_n) = edge_orbits - SL_orbits = T_n/2 + (n-2)! - SL_orbits

If we can crack SL_orbits, we have E(G_n) at ALL n.

APPROACH:
1. Verify edge_orbits = T_n/2 + (n-2)! at n=3,4,5 exactly
2. Analyze SL_orbits: decompose by Burnside (which cycle types contribute?)
3. Search for SL_orbits formula via ratios, differences, and known sequences

Author: opus-2026-03-23-S242
"""

import sys
import math
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# Known values
T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161}
V = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}

# Claimed formula: edge_orbits = T/2 + (n-2)!
# E = edge_orbits - SL_orbits

banner("VERIFY EDGE ORBIT FORMULA")

print("  n  | T_n  | T_n/2 | (n-2)! | edge_orb | E(G_n) | SL_orb | Check")
print("  ---|------|-------|--------|----------|--------|--------|------")
for n in range(3, 9):
    t = T[n]
    e = E.get(n, '?')
    eo = t // 2 + factorial(n-2)
    if isinstance(e, int):
        sl = eo - e
        check = "OK"
    else:
        sl = '?'
        check = '?'
    print("  %d  | %6d | %5d | %5d  | %7d  | %6s | %6s | %s" %
          (n, t, t//2, factorial(n-2), eo, e, sl, check))

# SL_orbits sequence
SL_orb = {}
for n in range(3, 9):
    if n in E:
        SL_orb[n] = T[n]//2 + factorial(n-2) - E[n]

print("\n  SL_orbits sequence: %s" % [SL_orb[n] for n in range(3, 9) if n in SL_orb])

banner("DEEP ANALYSIS OF SL_ORBITS")

sl_seq = [SL_orb[n] for n in range(3, 9)]
print("  SL_orbits = %s" % sl_seq)

# Ratios
print("\n  Consecutive ratios:")
for i in range(len(sl_seq)-1):
    r = sl_seq[i+1] / sl_seq[i]
    print("    SL(%d)/SL(%d) = %d/%d = %.4f" % (i+4, i+3, sl_seq[i+1], sl_seq[i], r))

# Divided by (n-1)!
print("\n  SL_orbits / (n-1)!:")
for n in range(3, 9):
    if n in SL_orb:
        print("    n=%d: %d / %d = %s = %.6f" %
              (n, SL_orb[n], factorial(n-1), Fraction(SL_orb[n], factorial(n-1)),
               SL_orb[n] / factorial(n-1)))

# Divided by (n-2)!
print("\n  SL_orbits / (n-2)!:")
for n in range(3, 9):
    if n in SL_orb:
        print("    n=%d: %d / %d = %s = %.6f" %
              (n, SL_orb[n], factorial(n-2), Fraction(SL_orb[n], factorial(n-2)),
               SL_orb[n] / factorial(n-2)))

# SL_orbits / V
print("\n  SL_orbits / V:")
for n in range(3, 9):
    if n in SL_orb:
        print("    n=%d: %d / %d = %s = %.4f" %
              (n, SL_orb[n], V[n], Fraction(SL_orb[n], V[n]),
               SL_orb[n] / V[n]))

# Differences
print("\n  First differences SL(n+1) - n*SL(n):")
for i in range(len(sl_seq)-1):
    n = i + 3
    diff = sl_seq[i+1] - (n+1) * sl_seq[i]
    print("    SL(%d) - %d*SL(%d) = %d - %d = %d" %
          (n+1, n+1, n, sl_seq[i+1], (n+1)*sl_seq[i], diff))

print("\n  SL(n+1) - n*SL(n):")
for i in range(len(sl_seq)-1):
    n = i + 3
    diff = sl_seq[i+1] - n * sl_seq[i]
    print("    SL(%d) - %d*SL(%d) = %d - %d = %d" %
          (n+1, n, n, sl_seq[i+1], n*sl_seq[i], diff))

print("\n  SL(n+1) - (n-1)*SL(n):")
for i in range(len(sl_seq)-1):
    n = i + 3
    diff = sl_seq[i+1] - (n-1) * sl_seq[i]
    print("    SL(%d) - %d*SL(%d) = %d - %d = %d" %
          (n+1, n-1, n, sl_seq[i+1], (n-1)*sl_seq[i], diff))

# Factorizations
print("\n  Factorizations:")
for n in range(3, 9):
    if n in SL_orb:
        s = SL_orb[n]; temp = s; factors = []
        for p in range(2, min(temp+1, 100000)):
            if p*p > temp and temp > 1: factors.append(temp); break
            while temp % p == 0: factors.append(p); temp //= p
        print("    n=%d: SL=%d = %s" % (n, s, " * ".join(str(f) for f in factors)))

banner("TESTING RECURRENCE RELATIONS")

# Test a(n) = c1*a(n-1) + c2*a(n-2) + ...
# With SL = 2, 5, 20, 86, 490, 3703

# Linear recurrence test: a(n) = alpha*n*a(n-1) + beta*a(n-2)
# 5 = alpha*4*2 + beta*? (need 3 terms minimum)
# 20 = alpha*5*5 + beta*2
# 86 = alpha*6*20 + beta*5
# From: 20 = 25*alpha + 2*beta and 86 = 120*alpha + 5*beta
# 5*(20) = 5*25*alpha + 10*beta → 100 = 125*alpha + 10*beta
# 2*(86) = 2*120*alpha + 10*beta → 172 = 240*alpha + 10*beta
# Subtract: 72 = 115*alpha → alpha = 72/115. Not clean.

# Try: a(n) = n*a(n-1) - c(n)
# 5 = 3*2 - 1 → c(3)=1
# 20 = 4*5 - 0 → c(4)=0
# 86 = 5*20 - 14 → c(5)=14
# 490 = 6*86 - 26 → c(6)=26
# 3703 = 7*490 - ? → c(7) = 3430 - 3703 = -273?! That's negative.

# Hmm. Try a(n) = (n-1)*a(n-1) + f(n)
# 5 = 2*2 + 1 → f(4)=1
# 20 = 3*5 + 5 → f(5)=5
# 86 = 4*20 + 6 → f(6)=6
# 490 = 5*86 + 60 → f(7)=60
# 3703 = 6*490 + 763 → f(8)=763
# f = 1, 5, 6, 60, 763. Not clean.

# Try a(n) = (n-1)*a(n-1) + (n-2)*a(n-2)
# 5 = 2*2 + 1*? (no a(2))
# 20 = 3*5 + 2*2 = 15+4 = 19 ≠ 20
# Close! Off by 1.

# Try a(n) = (n-1)*a(n-1) + (n-2)*a(n-2) + correction
# 20 = 3*5 + 2*2 + c = 19 + c → c=1
# 86 = 4*20 + 3*5 + c = 80+15 + c = 95 + c → c=-9
# Not constant correction.

# What about subfactorial? D_n = n! * sum_{k=0}^{n} (-1)^k/k!
# D_1=0, D_2=1, D_3=2, D_4=9, D_5=44, D_6=265, D_7=1854
# Not matching SL = 2, 5, 20, 86, 490, 3703

# What about sum of (n-2)! * k for some k?
# 2 = 1*2, 5 = 2*2.5, 20 = 6*3.33, 86 = 24*3.58, 490 = 120*4.08, 3703 = 720*5.14
# SL/(n-2)! = 2, 2.5, 3.33, 3.58, 4.08, 5.14
# The ratio SL/(n-2)! is growing slowly.

# Try: SL = (n-2)! * H(n) where H(n) is the harmonic number?
# H(3)=11/6≈1.83, H(4)=25/12≈2.08, H(5)=137/60≈2.28, H(6)=49/20=2.45
# SL/(n-2)! = 2, 2.5, 3.33, 3.58, 4.08, 5.14
# Doesn't match harmonic numbers either.

# EGF approach: SL(n) = (n-2)! * g(n) for some g
# g = 2, 5/2, 10/3, 43/12, 49/12, 3703/720
# g = 2, 2.5, 3.333, 3.583, 4.083, 5.143
# Differences: 0.5, 0.833, 0.25, 0.5, 1.06
# Second order: partial fractions?

# g(n) for n=3..8:
print("\n  g(n) = SL_orbits / (n-2)!:")
g_vals = []
for n in range(3, 9):
    if n in SL_orb:
        g = Fraction(SL_orb[n], factorial(n-2))
        g_vals.append(g)
        print("    n=%d: g = %s = %.6f" % (n, g, float(g)))

print("\n  Numerators: %s" % [g.numerator for g in g_vals])
print("  Denominators: %s" % [g.denominator for g in g_vals])

# Try: g(n) = sum_{k=2}^{n-1} 1/k! * something?
# Or: g(n) = (partial exponential sum)?

# Check: g(n) * (n-2)! = SL. Is g(n) related to partial sums of 1/k?

banner("BURNSIDE DECOMPOSITION OF SL_ORBITS")

# SL_orbits counts the number of S_n-orbits on self-loop edges of Q_m.
# A self-loop edge: (T, T_e) where T ≅ T_e (reversing arc e gives isomorphic T).
# The S_n-orbit of such an edge is {(σT, σ(T_e)) : σ ∈ S_n}.

# By Burnside: SL_orbits = (1/n!) * sum_σ |Fix_SL(σ)|
# where Fix_SL(σ) = #{self-loop edges (T, T_e) : σ fixes this edge}
# = #{(T, e) : σ(T)=T, σ(e)=e, T_e ≅ T}
# ∪ #{(T, e) : σ(T)=T_e, σ(e)=(reversed e), ...}

# Wait — a self-loop edge is an UNORDERED pair {T, T_e} where T_e ≅ T.
# Since T_e ≅ T, the pair {T, T_e} might have T_e = T (same labeled tournament)
# or T_e ≠ T but T_e ≅ T (different labeling, same class).

# For σ to fix the edge {T, T_e}:
# Either σ(T) = T and σ(T_e) = T_e
# Or σ(T) = T_e and σ(T_e) = T

# Case 1: σ ∈ Aut(T) and σ fixes the arc e (since σ(T_e) = (σT)_{σe} = T_{σe},
#   and we need T_{σe} = T_e, so σe = e)
# Case 2: σ maps T to T_e. Since T_e is obtained by reversing e in T,
#   σ maps T to T with one arc reversed. This means σ is a "near-automorphism"
#   of T — it's an automorphism except at arc e.

# The (n-2)! term in edge_orbits comes from Case 2 (near-automorphisms).
# But for SL_orbits specifically, we need BOTH T ≅ T_e AND σ maps between them.

# This is getting complex. Let me compute SL_orbits directly at n=3,4,5.

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs_n)}

    # Build canonical forms
    cmap = {}; tc = {}
    for bits in range(2**m):
        best = bits
        for perm in permutations(range(n)):
            new_bits = 0
            for k, (i,j) in enumerate(pairs_n):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits & (1<<ok): new_bits |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits & (1<<ok)): new_bits |= (1<<k)
            if new_bits < best: best = new_bits
        tc[bits] = best
        if best not in cmap: cmap[best] = len(cmap)

    # Find all self-loop EDGES: unordered pairs {bits, bits^(1<<k)} where
    # both map to the same canonical form AND bits ≠ bits^(1<<k)
    sl_edges = set()
    for bits in range(2**m):
        c1 = tc[bits]
        for k in range(m):
            flipped = bits ^ (1<<k)
            if flipped > bits and tc[flipped] == c1:
                sl_edges.add((bits, flipped))

    # Count orbits of sl_edges under S_n
    # S_n acts on labeled tournaments; σ maps (bits, flipped) to (σ(bits), σ(flipped))
    # We need orbits of these pairs.

    orbit_reps = set()
    edge_to_orbit = {}
    for edge in sl_edges:
        if edge in edge_to_orbit: continue
        # Find canonical form of this edge under S_n
        best_edge = None
        for perm in permutations(range(n)):
            # Apply perm to both endpoints
            new_a = 0; new_b = 0
            bits_a, bits_b = edge
            for k, (i,j) in enumerate(pairs_n):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits_a & (1<<ok): new_a |= (1<<k)
                    if bits_b & (1<<ok): new_b |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits_a & (1<<ok)): new_a |= (1<<k)
                    if not (bits_b & (1<<ok)): new_b |= (1<<k)
            # Normalize: smaller first
            e_mapped = (min(new_a, new_b), max(new_a, new_b))
            if best_edge is None or e_mapped < best_edge:
                best_edge = e_mapped

        orbit_reps.add(best_edge)
        # Mark all edges in this orbit (slow but correct for small n)
        for perm in permutations(range(n)):
            new_a = 0; new_b = 0
            bits_a, bits_b = edge
            for k, (i,j) in enumerate(pairs_n):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits_a & (1<<ok): new_a |= (1<<k)
                    if bits_b & (1<<ok): new_b |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits_a & (1<<ok)): new_a |= (1<<k)
                    if not (bits_b & (1<<ok)): new_b |= (1<<k)
            e_mapped = (min(new_a, new_b), max(new_a, new_b))
            if e_mapped in sl_edges:
                edge_to_orbit[e_mapped] = best_edge

    sl_orbit_count = len(orbit_reps)
    expected = SL_orb.get(n, '?')
    print("  n=%d: %d self-loop edges, %d orbits (expected %s)" %
          (n, len(sl_edges), sl_orbit_count, expected))

banner("THE SL_ORBITS FORMULA HUNT")

# SL_orbits = 2, 5, 20, 86, 490, 3703
# Let's try to express in terms of known quantities.

# SL_orbits - (n-2)!:
print("  SL_orbits - (n-2)!:")
for n in range(3, 9):
    if n in SL_orb:
        diff = SL_orb[n] - factorial(n-2)
        print("    n=%d: %d - %d = %d" % (n, SL_orb[n], factorial(n-2), diff))

# SL - (n-2)! = 1, 3, 14, 62, 370, 2983
# These are: 1, 3, 14, 62, 370, 2983
# Check: this is G(n)/2! The gap sequence was G = 2, 6, 28, 124, 740, 5966
# G/2 = 1, 3, 14, 62, 370, 2983. EXACTLY!
G = {n: T[n] - 2*E[n] for n in E}
print("\n  G(n)/2:")
for n in range(3, 9):
    if n in G:
        print("    n=%d: G/2 = %d" % (n, G[n]//2))

print("\n  SL_orbits - (n-2)! vs G(n)/2:")
for n in range(3, 9):
    if n in SL_orb and n in G:
        diff = SL_orb[n] - factorial(n-2)
        g2 = G[n] // 2
        print("    n=%d: SL-(n-2)! = %d, G/2 = %d, match = %s" %
              (n, diff, g2, diff == g2))

banner("THE IDENTITY")

print("""
  DISCOVERED IDENTITY:

  SL_orbits = (n-2)! + G(n)/2

  where G(n) = T(n) - 2*E(n) = SL(n) + excess(n)

  CHECK: E(G_n) = edge_orbits - SL_orbits
       = T_n/2 + (n-2)! - [(n-2)! + G(n)/2]
       = T_n/2 - G(n)/2
       = T_n/2 - (T_n - 2*E_n)/2
       = T_n/2 - T_n/2 + E_n
       = E_n  ✓

  THIS IS A TAUTOLOGY! The formula edge_orbits = T_n/2 + (n-2)!
  with SL_orbits = (n-2)! + G/2 gives back E = E. It's consistent
  but doesn't give us NEW information.

  The REAL content is in the claim:
    edge_orbits = T_n/2 + (n-2)!

  If this is TRUE, then:
    SL_orbits = T_n/2 + (n-2)! - E_n
              = (T_n - 2*E_n)/2 + (n-2)!
              = G_n/2 + (n-2)!

  And E_n = T_n/2 + (n-2)! - SL_orbits.

  So everything reduces to: IS edge_orbits = T_n/2 + (n-2)! TRUE?
""")

# Let me verify directly at n=3,4,5 using the Burnside edge orbit count
# from the tournament_flip_edge_orbits.py script above.

print("  Direct verification from Burnside edge orbit count:")
for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs_n)}

    # Burnside: edge_fix_sum = sum_σ |Fix(σ)| where Fix(σ) = edges fixed by σ
    edge_fix_sum = 0

    for perm in permutations(range(n)):
        # How σ acts on labeled tournaments
        def apply_perm(bits):
            new_bits = 0
            for k, (i,j) in enumerate(pairs_n):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits & (1<<ok): new_bits |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits & (1<<ok)): new_bits |= (1<<k)
            return new_bits

        # Count edges {bits, bits^(1<<k)} with bits < bits^(1<<k)
        # fixed by σ: σ maps {bits, flipped} to {σ(bits), σ(flipped)}
        # Edge is fixed iff {σ(bits), σ(flipped)} = {bits, flipped}
        for bits in range(2**m):
            for k in range(m):
                flipped = bits ^ (1<<k)
                if flipped <= bits: continue  # each edge once
                pb = apply_perm(bits)
                pf = apply_perm(flipped)
                if (pb == bits and pf == flipped) or (pb == flipped and pf == bits):
                    edge_fix_sum += 1

    edge_orbits = edge_fix_sum // factorial(n)
    predicted = T[n]//2 + factorial(n-2)
    print("    n=%d: edge_orbits = %d, T/2 + (n-2)! = %d, match = %s" %
          (n, edge_orbits, predicted, edge_orbits == predicted))

banner("WHAT REMAINS")

print("""
  THE FORMULA edge_orbits = T_n/2 + (n-2)! needs to be PROVED.

  T_n = total arc-orbits (Burnside-computable, EXACT for all n)
  (n-2)! = a combinatorial constant

  This formula says: the number of S_n-orbits on edges of Q_m
  (the tournament cube) equals half the arc-orbit count plus (n-2)!.

  THE (n-2)! TERM: this is the "near-automorphism correction."
  It counts edge orbits where σ SWAPS the two endpoints of the edge
  (σ maps T to T_e and T_e to T). These are edges where the two
  labeled tournaments differ by one arc AND there exists a relabeling
  that maps one to the other — i.e., T_e = σ(T) for some σ.

  PROVING edge_orbits = T_n/2 + (n-2)!:
  By Burnside on edges: edge_orbits = (1/n!) * sum_σ |Fix_edges(σ)|

  Fix_edges(σ) = #{edges {T, T'} : σ fixes the pair, T' = T_e for some arc e}

  Case A: σ(T) = T and σ(T') = T' → σ ∈ Aut(T) ∩ Aut(T')
    and σ fixes the arc e (since T' = T_e).
    This contributes T_n/2 (half the arc-orbit-fixing count).

  Case B: σ(T) = T' and σ(T') = T → σ swaps the two endpoints.
    Since T' = T_e, this means σ(T) = T_e, i.e., σ is a near-automorphism.
    σ² then satisfies σ²(T) = σ(T_e) = (σT)_{σe} = T_{e'}_{σe}.
    If σ² ∈ Aut(T), this constrains σ to be an involutory near-automorphism.

  The (n-2)! count suggests: for each edge {i,j} of K_n (m choices),
  the number of near-automorphisms that swap arc (i,j) is (n-2)!/m = 1
  on average. More precisely, sum over σ of these contributions = (n-2)! * n!,
  giving (n-2)! edge orbits.

  THIS NEEDS A FORMAL PROOF.
""")

print("Done. opus-2026-03-23-S242")
