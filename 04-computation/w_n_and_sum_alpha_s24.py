#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Derive closed form for W(n) = Σ_T H(T) over labeled tournaments.

W(n) = Σ_T I(Ω(T), 2) = Σ_T (Σ_k α_k(T) 2^k)
     = Σ_k 2^k × Σ_T α_k(T)

where α_k(T) = # independent sets of size k in the conflict graph Ω(T).

Key formula for Σ_T α_1(T):
For each directed odd cycle C of length L in a specific labeled tournament,
the number of labeled tournaments containing C is 2^{C(n,2)-L}.
The number of labeled directed cycles of length L on [n] is C(n,L) × (L-1)!.

So Σ_T α_1(T) = Σ_{L odd, 3≤L≤n} C(n,L) × (L-1)! × 2^{C(n,2)-L}

Similarly for α_2(T): count pairs of vertex-disjoint odd cycles.
"""

from math import comb, factorial
from itertools import combinations
from collections import Counter
import random

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

# ═══════════════════════════════════════════
# PART 1: Compute W(n) = Σ_T H(T) exactly
# ═══════════════════════════════════════════
print("="*60)
print("PART 1: W(n) = Σ_T H(T) (exact)")
print("="*60)

for n in range(1, 8):
    m = n*(n-1)//2
    if 2**m > 5000000:
        break
    total = 0
    for T in all_tournaments(n):
        total += hamiltonian_path_count(T)
    print(f"  n={n}: W(n) = {total}, 2^C(n,2) = {2**m}, W/2^m = {total / 2**m:.6f}")

# ═══════════════════════════════════════════
# PART 2: Theoretical formula for Σ_T α_1(T)
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 2: Σ_T α_1(T) — total odd cycles across all tournaments")
print("="*60)

def sum_alpha1_theoretical(n):
    """Total number of (tournament, odd cycle) pairs = Σ_T α_1(T)."""
    m = n*(n-1)//2
    total = 0
    for L in range(3, n+1, 2):
        # Number of labeled directed L-cycles on [n]
        n_cycles = comb(n, L) * factorial(L-1)
        # Each such cycle exists in 2^{m-L} tournaments
        n_tournaments_with = 2**(m - L)
        contribution = n_cycles * n_tournaments_with
        total += contribution
        print(f"    L={L}: {n_cycles} cycles × 2^{m-L} tournaments = {contribution}")
    return total

for n in range(3, 9):
    print(f"\n  n={n}:")
    theoretical = sum_alpha1_theoretical(n)
    m = n*(n-1)//2
    print(f"    Σ_T α_1(T) = {theoretical}")
    print(f"    = {theoretical} / {2**m} per tournament = {theoretical / 2**m:.6f}")

# Verify against exact computation
print(f"\n  Verification against exact:")
from itertools import permutations as perms

def count_odd_cycles(T):
    """Count ALL directed odd cycles in T."""
    n = len(T)
    count = 0
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            for perm in perms(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(L):
                    if T[path[i]][path[(i+1) % L]] != 1:
                        ok = False
                        break
                if ok:
                    count += 1
                    break
    return count

for n in [3, 4, 5]:
    exact = sum(count_odd_cycles(T) for T in all_tournaments(n))
    theoretical = sum_alpha1_theoretical(n)
    match = "✓" if exact == theoretical else "✗"
    print(f"    n={n}: exact={exact}, formula={theoretical} {match}")

# ═══════════════════════════════════════════
# PART 3: Formula for Σ_T α_2(T)
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 3: Σ_T α_2(T) — pairs of vertex-disjoint odd cycles")
print("="*60)

def sum_alpha2_theoretical(n):
    """Σ_T α_2(T) = total (tournament, pair of disjoint odd cycles) triples."""
    m = n*(n-1)//2
    total = 0
    # For each pair of odd lengths (L1, L2) with L1+L2 ≤ n
    for L1 in range(3, n+1, 2):
        for L2 in range(3, n-L1+1, 2):
            # Number of labeled pairs of disjoint directed cycles
            # Choose L1 vertices, then L2 from remaining
            # But need to account for ordering of the pair
            if L1 < L2:
                # Ordered pairs: 2 orderings
                n_pairs = comb(n, L1) * factorial(L1-1) * comb(n-L1, L2) * factorial(L2-1)
            elif L1 == L2:
                # Same length: divide by 2 for unordered pairs
                n_pairs = comb(n, L1) * factorial(L1-1) * comb(n-L1, L2) * factorial(L2-1) // 2
            else:
                continue  # already counted in L1 < L2 case

            # Each pair exists in 2^{m - L1 - L2} tournaments
            # (all arcs of both cycles fixed, everything else free)
            # Wait — arcs BETWEEN the two cycles' vertex sets are also free
            # The cycles fix L1 + L2 arcs (each cycle's edges)
            # Remaining arcs: m - L1 - L2 ... no, that's wrong.
            # A cycle of length L on vertices S uses L arcs.
            # The C(L,2)-L arcs among S but NOT in the cycle are free.
            # And arcs between S1 and S2 are free (L1*L2 arcs).
            # And arcs from S1∪S2 to rest are free.
            # Total free = m - L1 - L2
            n_tournaments_with = 2**(m - L1 - L2)
            contribution = n_pairs * n_tournaments_with
            total += contribution
            if contribution > 0:
                print(f"    (L1,L2)=({L1},{L2}): {n_pairs} pairs × 2^{m-L1-L2} = {contribution}")
    return total

for n in [3, 4, 5, 6]:
    print(f"\n  n={n}:")
    theoretical = sum_alpha2_theoretical(n)
    m = n*(n-1)//2
    print(f"    Σ_T α_2(T) = {theoretical}")

# Verify α_2 formula
def count_alpha2_exact(T):
    """Count α_2 = number of independent pairs in Ω(T)."""
    n = len(T)
    cycles = []
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            for perm in perms(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(L):
                    if T[path[i]][path[(i+1) % L]] != 1:
                        ok = False
                        break
                if ok:
                    mi = path.index(min(path))
                    canon = tuple(path[(mi + j) % L] for j in range(L))
                    cycles.append(canon)
                    break

    # Count vertex-disjoint pairs
    vsets = [set(c) for c in cycles]
    alpha2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (vsets[i] & vsets[j]):
                alpha2 += 1
    return alpha2

print(f"\n  Verification:")
for n in [3, 4, 5]:
    exact = sum(count_alpha2_exact(T) for T in all_tournaments(n))
    theoretical = sum_alpha2_theoretical(n)
    match = "✓" if exact == theoretical else "✗"
    print(f"    n={n}: exact Σα₂={exact}, formula={theoretical} {match}")

# ═══════════════════════════════════════════
# PART 4: W(n) from the formula
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 4: W(n) = 2^m + 2·Σα₁ + 4·Σα₂ + ...")
print("="*60)

for n in range(3, 9):
    m = n*(n-1)//2
    sa1 = sum_alpha1_theoretical(n)
    sa2 = sum_alpha2_theoretical(n)

    W_approx = 2**m + 2*sa1 + 4*sa2
    print(f"\n  n={n}: 2^m={2**m}, 2·Σα₁={2*sa1}, 4·Σα₂={4*sa2}")
    print(f"    W(n) ≈ {W_approx} (using only α₁, α₂ terms)")

# Compare with exact W(n) for small n
print(f"\nComparison with exact W(n):")
w_exact = {}
for n in range(1, 7):
    m = n*(n-1)//2
    w = sum(hamiltonian_path_count(T) for T in all_tournaments(n))
    w_exact[n] = w

for n in range(3, 7):
    m = n*(n-1)//2
    sa1 = sum_alpha1_theoretical(n)
    sa2 = sum_alpha2_theoretical(n)
    W_formula = 2**m + 2*sa1 + 4*sa2
    print(f"  n={n}: W_exact={w_exact[n]}, W_formula(α₁+α₂)={W_formula}, "
          f"diff={w_exact[n]-W_formula} (= 8·Σα₃ + ...)")

# Compute Σ_T α_3 from the difference
print(f"\n  Extracting Σα₃:")
for n in range(3, 7):
    m = n*(n-1)//2
    sa1 = sum_alpha1_theoretical(n)
    sa2 = sum_alpha2_theoretical(n)
    diff = w_exact[n] - 2**m - 2*sa1 - 4*sa2
    if diff % 8 == 0:
        sa3 = diff // 8
        print(f"  n={n}: Σα₃ = {sa3}")
    else:
        print(f"  n={n}: diff={diff}, NOT divisible by 8 — higher terms matter")

# ═══════════════════════════════════════════
# PART 5: Mean H per tournament
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 5: E[H] = W(n) / 2^{C(n,2)}")
print("="*60)

for n in [1, 2, 3, 4, 5, 6]:
    m = n*(n-1)//2
    if n in w_exact:
        mean_H = w_exact[n] / 2**m
        print(f"  n={n}: E[H] = {w_exact[n]} / {2**m} = {mean_H:.6f}")

# Theoretical: E[H] = 1 + 2·E[α₁] + 4·E[α₂] + ...
# E[α₁] = Σ_{L odd} C(n,L)·(L-1)! / 2^L
print(f"\n  E[α₁] = Σ_L C(n,L)·(L-1)!/2^L:")
for n in range(3, 12):
    ea1 = 0
    for L in range(3, n+1, 2):
        ea1 += comb(n, L) * factorial(L-1) / 2**L
    mean_H_approx = 1 + 2*ea1
    print(f"  n={n}: E[α₁] = {ea1:.4f}, 1+2·E[α₁] = {mean_H_approx:.4f}")

# ═══════════════════════════════════════════
# PART 6: Asymptotic behavior of E[α₁]
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 6: Asymptotic E[α₁]")
print("="*60)

print("  E[α₁] = Σ_{L=3,5,...} C(n,L)·(L-1)!/2^L")
print("  Dominant term at large n: L=3 gives C(n,3)·2!/2^3 = C(n,3)/4")

for n in range(3, 15):
    ea1 = 0
    terms = {}
    for L in range(3, n+1, 2):
        term = comb(n, L) * factorial(L-1) / 2**L
        terms[L] = term
        ea1 += term

    frac_3 = terms.get(3, 0) / ea1 if ea1 > 0 else 0
    print(f"  n={n:2d}: E[α₁]={ea1:12.2f}, L=3 term={terms.get(3,0):10.2f} ({100*frac_3:.1f}%), "
          f"L=5 term={terms.get(5,0):10.2f}")

print(f"\n  E[α₁] for L=3 only: C(n,3)/4 = n(n-1)(n-2)/24")
print(f"  E[α₁] for L=3+5: C(n,3)/4 + C(n,5)·4!/2^5 = C(n,3)/4 + 3C(n,5)/4")

# ═══════════════════════════════════════════
# PART 7: The W(n) sequence — look for patterns
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("PART 7: W(n) patterns")
print("="*60)

# Known: W = 1, 2, 12, ...
# Wait, let me recompute W correctly
print("  Computing W(n) = Σ_T H(T) for n=1..6:")
for n in range(1, 7):
    if n in w_exact:
        print(f"  W({n}) = {w_exact[n]}")

# Ratios
print(f"\n  Ratios W(n)/W(n-1):")
prev = None
for n in range(1, 7):
    if n in w_exact:
        if prev is not None:
            print(f"  W({n})/W({n-1}) = {w_exact[n]/prev:.4f}")
        prev = w_exact[n]

# W(n) / (n! · 2^{C(n-1,2)})
print(f"\n  W(n) / (n! · 2^{{C(n-1,2)}}):")
for n in range(1, 7):
    if n in w_exact:
        denom = factorial(n) * 2**((n-1)*(n-2)//2)
        ratio = w_exact[n] / denom
        print(f"  n={n}: {w_exact[n]} / ({factorial(n)} × {2**((n-1)*(n-2)//2)}) = {ratio:.6f}")

# W(n) / 2^{C(n,2)}
print(f"\n  E[H(T)] = W(n) / 2^C(n,2):")
for n in range(1, 7):
    if n in w_exact:
        mean = w_exact[n] / 2**(n*(n-1)//2)
        print(f"  n={n}: E[H] = {mean:.6f}")

print("\n\nDONE.")
