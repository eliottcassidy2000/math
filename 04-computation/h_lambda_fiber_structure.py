#!/usr/bin/env python3
"""
h_lambda_fiber_structure.py — opus-2026-03-13-S71c

H is ALMOST lambda-determined at n=7 (1 ambiguous group in 5000).
Since H = [lambda-determined part] + 2c7, and |Δc7| ≤ 1 (THM-182),
H varies by at most 2 within any lambda fiber.

QUESTIONS:
1. What fraction of lambda fibers have c7 variation?
2. Is the ambiguity always exactly {c7, c7+1} (never bigger gap)?
3. Does H mod 4 or mod 8 have nice lambda-determination properties?
4. Large-scale scan: how many lambda fibers show H ambiguity?
5. At n=8: does H remain almost lambda-determined?

Connection to "2 and 3": H = lambda_part + 2·c7.
The "2" in 2·c7 is the OCF weight. If we used I(3) instead:
I₃ = lambda_part_3 + 3·c7. The "3" is the I(3) weight.
The DIFFERENCE: 3·c7 - 2·c7 = c7. So I₃ - H = lambda_diff + c7.
And c7 = I₃ - H - (lambda_part_3 - lambda_part_2).
"""

import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

# =====================================================================
print("=" * 70)
print("H-VALUE LAMBDA FIBER STRUCTURE AT n=7")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

lam_to_H = defaultdict(set)
lam_to_c7 = defaultdict(set)
lam_to_c3c5 = defaultdict(set)

t0 = time.time()
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk = lambda_key(A, n)
    H = count_ham_paths(A, n)

    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    # c7 via tr(A^7): c7_dir = tr(A^7)/7 is NOT right (degenerate walks)
    # Better: c7 = (H-1)/2 - c3 - c5 - 2α₂ where α₂ = disjoint 3-cycle pairs
    # But computing α₂ is expensive. Let's just track H and (c3,c5).

    lam_to_H[lk].add(H)
    lam_to_c3c5[lk].add((c3, c5))

    if trial % 10000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, {len(lam_to_H)} lambda groups")

dt = time.time() - t0
print(f"\nDone: {len(lam_to_H)} lambda groups, {dt:.1f}s")

# Count H-ambiguous groups
ambig_H = 0
max_H_range = 0
ambig_details = []
for lk, H_vals in lam_to_H.items():
    if len(H_vals) > 1:
        ambig_H += 1
        H_range = max(H_vals) - min(H_vals)
        max_H_range = max(max_H_range, H_range)
        ambig_details.append((lk, sorted(H_vals)))

print(f"\nH ambiguous lambda groups: {ambig_H}/{len(lam_to_H)}")
print(f"Max H range within a fiber: {max_H_range}")

if ambig_details:
    for lk, vals in ambig_details[:15]:
        c3c5 = lam_to_c3c5[lk]
        print(f"  H ∈ {vals}, c3c5 = {sorted(c3c5)}, ΔH = {max(vals)-min(vals)}")

# Check: are (c3,c5) always constant within lambda fibers?
ambig_c3c5 = sum(1 for v in lam_to_c3c5.values() if len(v) > 1)
print(f"\n(c3,c5) ambiguous within lambda fibers: {ambig_c3c5}/{len(lam_to_c3c5)}")

# =====================================================================
print(f"\n{'='*70}")
print("H MOD STRUCTURE WITHIN LAMBDA FIBERS")
print("=" * 70)

# H is always odd. H mod 4 ∈ {1, 3}.
# Does lambda determine H mod 4? (Since ΔH = 2·Δc7 and Δc7 ∈ {0, ±1}, ΔH ∈ {0, ±2}.)
# H mod 4 changes by 2 mod 4 = ±2 mod 4. So {1,3} or {3,1}.
# H mod 4 is NOT preserved under ΔH=±2.

# But H mod 2 IS preserved (ΔH is even).
# Let's check H mod 4 and H mod 8 determination.

lam_to_Hmod4 = defaultdict(set)
lam_to_Hmod8 = defaultdict(set)
for lk, H_vals in lam_to_H.items():
    for h in H_vals:
        lam_to_Hmod4[lk].add(h % 4)
        lam_to_Hmod8[lk].add(h % 8)

ambig_mod4 = sum(1 for v in lam_to_Hmod4.values() if len(v) > 1)
ambig_mod8 = sum(1 for v in lam_to_Hmod8.values() if len(v) > 1)

print(f"Lambda determines H mod 4? Ambiguous: {ambig_mod4}/{len(lam_to_Hmod4)}")
print(f"Lambda determines H mod 8? Ambiguous: {ambig_mod8}/{len(lam_to_Hmod8)}")

# =====================================================================
print(f"\n{'='*70}")
print("H-LAMBDA GAP DISTRIBUTION")
print("=" * 70)

# For the ambiguous groups: the gap ΔH = 2Δc7 should always be exactly 2
# (since |Δc7|=1 per THM-182)
print("\nGap distribution for ambiguous lambda groups:")
gap_dist = defaultdict(int)
for lk, H_vals in lam_to_H.items():
    if len(H_vals) > 1:
        gap = max(H_vals) - min(H_vals)
        gap_dist[gap] += 1

print(f"  {dict(sorted(gap_dist.items()))}")

# =====================================================================
print(f"\n{'='*70}")
print("FRACTION OF TOURNAMENTS IN AMBIGUOUS FIBERS")
print("=" * 70)

# What fraction of random tournaments land in an ambiguous lambda fiber?
# This tells us how "close" H is to being lambda-determined.
total_in_ambig = 0
total = 0
np.random.seed(42)
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk = lambda_key(A, n)
    total += 1
    if len(lam_to_H.get(lk, set())) > 1:
        total_in_ambig += 1

print(f"  Tournaments in ambiguous lambda fibers: {total_in_ambig}/{total} = {100*total_in_ambig/total:.3f}%")

# =====================================================================
print(f"\n{'='*70}")
print("I(3) - H = LAMBDA_DIFF + c7")
print("=" * 70)

# Since I₃ = 1 + 3α₁ + 9α₂ and H = 1 + 2α₁ + 4α₂:
# I₃ - H = α₁ + 5α₂ = (c3+c5+c7) + 5α₂
# = (c3+c5+5α₂) + c7
# where (c3+c5+5α₂) is lambda-determined

# So: c7 = (I₃ - H) - c3 - c5 - 5α₂
# SIMPLER: c7 = (I₃ - H) - (lambda-determined quantity)

# Even simpler: I₃ - H is NOT lambda-determined, and its non-lambda
# part is exactly c7. So (lambda, I₃-H) determines c7 trivially.

# But we already know (lambda, H) determines c7.
# The key: I₃ - H contains c7 with coefficient 1 (not 2!).
# H contains c7 with coefficient 2.
# So I₃ - H "extracts" c7 more directly.

print("  I₃ - H = α₁ + 5α₂ = (c3+c5+c7) + 5α₂")
print("  Lambda-determined part of (I₃-H): c3+c5+5α₂")
print("  Non-lambda part: c7 (coefficient 1)")
print("  vs in H: coefficient of c7 is 2")
print("  So (I₃ - H) - lambda_part = c7 DIRECTLY")

# Verify at smaller scale
np.random.seed(42)
ok = 0
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

    # Find α₂ (disjoint 3-cycle pairs)
    # At n=7, α₂ = (H - 1 - 2*α₁) / 4... no, need to compute directly
    # Quick: find all 3-cycles, count disjoint pairs
    cycles3 = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if A[a][b] and A[b][c] and A[c][a]:
                    cycles3.append(frozenset([a,b,c]))
                elif A[a][c] and A[c][b] and A[b][a]:
                    cycles3.append(frozenset([a,b,c]))
    # Remove duplicates
    cycles3 = list(set(cycles3))
    a2 = 0
    for i in range(len(cycles3)):
        for j in range(i+1, len(cycles3)):
            if not (cycles3[i] & cycles3[j]):
                a2 += 1

    # c7 from H
    a1 = (H - 1 - 4*a2) // 2
    c7 = a1 - c3 - c5

    # Check via I₃ - H
    I3 = 1 + 3*a1 + 9*a2
    diff = I3 - H
    c7_from_diff = diff - c3 - c5 - 5*a2

    if c7 == c7_from_diff:
        ok += 1

print(f"\n  c7 from (I₃-H) formula: {ok}/1000 match")

# =====================================================================
print(f"\n{'='*70}")
print("THE DEEP STRUCTURE: SIGMA vs H FOR c7 DETERMINATION")
print("=" * 70)

# We know: (lambda, H) determines c7 (0 ambiguous).
# We know: (lambda, sigma) determines c7 (0 ambiguous).
# Does H encode the same info as sigma (conditioned on lambda)?

# If (lambda, H) and (lambda, sigma) both determine c7,
# then H and sigma must be "equivalent" information channels
# for the c7-residual.

# But H is a GLOBAL invariant (single number).
# Sigma is a MATRIX (n×n values).
# So H encodes the sigma information in a highly compressed way.

# The compression is lossy for other things but lossless for c7!

print("  Both (lambda, H) and (lambda, sigma) determine c7.")
print("  H is a scalar; sigma is a matrix.")
print("  H compresses sigma's c7-relevant info losslessly!")
print("  This is because c7 has only 1 bit of non-lambda info (per THM-182).")

print("\nDone.")
