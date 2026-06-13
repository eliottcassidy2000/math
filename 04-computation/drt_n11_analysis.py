#!/usr/bin/env python3
"""
Analyze DRTs at n=11: cycle counts, H values, Paley maximizer test.

At n=11, there are multiple DRT isomorphism classes (4 according to McKay).
Questions:
1. Do all DRTs at n=11 have the same cycle counts? (Savchenko says yes for c_8)
2. Do all DRTs at n=11 have the same H(T)?
3. Is Paley T_11 the H-maximizer among all DRTs?

Key result from Savchenko: c_m(DR_n) > c_m(RLT_n) for m = 1,2,3 mod 4
This means ALL odd-length cycle counts are higher for DRT than LTT.

We already know H(T_11) = 95095.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import hamiltonian_path_count
from itertools import combinations

def build_paley(p):
    """Build Paley tournament T_p for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

def build_szekeres_whiteman(p):
    """Build Szekeres-Whiteman type DRT using 2-block construction.
    For p = 3 mod 8, the QRs split into two cosets of the 4th power residues.
    We swap one coset to get a different DRT."""
    # For p=11: QR = {1,3,4,5,9}, NQR = {2,6,7,8,10}
    # 4th power residues mod 11: x^4 mod 11 for x=1..10
    # 1^4=1, 2^4=5, 3^4=4, 4^4=4, 5^4=9, 6^4=9, 7^4=4, 8^4=5, 9^4=4, 10^4=1
    # So 4th power residues = {1,4,5,9}
    # Wait, let me be more careful. The multiplicative group mod 11 has order 10.
    # Index-4 subgroup has order 10/gcd(4,10) = 10/2 = 5... that's the QRs themselves.
    # Actually for Szekeres-Whiteman we need p = 3 mod 8.
    # 11 mod 8 = 3, so p=11 qualifies.
    # The construction: split {1,...,p-1} into C0 (QR) and C1 (NQR).
    # Paley uses i->j iff j-i in QR.
    # SW modifies by swapping the role of a specific coset.
    # Let g be a primitive root mod p. For p=11, g=2.
    # QR = {g^0, g^2, g^4, g^6, g^8} = {1, 4, 5, 9, 3}
    # NQR = {g^1, g^3, g^5, g^7, g^9} = {2, 8, 10, 7, 6}
    # The 2-block construction (Szekeres 1969):
    # Define D = {d in Z_p^* : d in S} where S is built from QR but with
    # one coset of index-2 subgroup of QR swapped.
    # Actually, the simplest non-Paley DRT construction at n=11:
    # Use a different difference set.
    pass

def build_drt_from_diff_set(n, diff_set):
    """Build tournament from a difference set D subset of Z_n.
    i->j iff j-i mod n is in D."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in diff_set:
                T[i][j] = 1
    return T

def is_doubly_regular(T, n):
    """Check if T is doubly regular."""
    k = (n-1) // 2
    for v in range(n):
        if sum(T[v]) != k:
            return False
    target = (n - 3) // 4
    for u in range(n):
        for v in range(u+1, n):
            common = sum(1 for w in range(n) if w != u and w != v and T[u][w] and T[v][w])
            if common != target:
                return False
    return True

def count_directed_3cycles(T, n):
    """Count directed 3-cycles."""
    count = 0
    for combo in combinations(range(n), 3):
        a, b, c = combo
        # Check both orientations
        if T[a][b] and T[b][c] and T[c][a]:
            count += 1
        if T[a][c] and T[c][b] and T[b][a]:
            count += 1
    return count

def count_directed_5cycles(T, n):
    """Count directed 5-cycles using DP."""
    total = 0
    for combo in combinations(range(n), 5):
        verts = list(combo)
        v0 = verts[0]
        dp = {}
        dp[(1, 0)] = 1
        for mask in range(1, 1 << 5):
            if not (mask & 1):
                continue
            for vi in range(5):
                if not (mask & (1 << vi)):
                    continue
                c = dp.get((mask, vi), 0)
                if c == 0:
                    continue
                for ui in range(5):
                    if mask & (1 << ui):
                        continue
                    if T[verts[vi]][verts[ui]]:
                        key = (mask | (1 << ui), ui)
                        dp[key] = dp.get(key, 0) + c
        full = (1 << 5) - 1
        for vi in range(1, 5):
            c = dp.get((full, vi), 0)
            if c > 0 and T[verts[vi]][verts[0]]:
                total += c
    return total

# ============================================================
# Build known DRTs at n=11
# ============================================================
print("=" * 70)
print("N=11: DRT ANALYSIS — CYCLE COUNTS AND HAMILTONIAN PATHS")
print("=" * 70)

n = 11

# Paley T_11: QR = {1, 3, 4, 5, 9}
paley_qr = set()
for x in range(1, 11):
    paley_qr.add((x*x) % 11)
print(f"Paley QR mod 11: {sorted(paley_qr)}")

T_paley = build_drt_from_diff_set(11, paley_qr)
assert is_doubly_regular(T_paley, 11), "Paley should be DRT"

# Try to construct non-Paley DRTs at n=11
# At n=11 (p=11, p=3 mod 4), the known DRTs come from:
# 1. Paley (QR difference set)
# 2-4. Non-isomorphic skew Hadamard constructions

# Method: try all possible (11,5,2)-difference sets in Z_11
# A (v,k,lambda)-difference set D has |D|=k, and every nonzero element
# of Z_v appears exactly lambda times as d1-d2 with d1,d2 in D.
# For DRT of order 11: need (11,5,2)-difference set.

print("\nSearching for all (11,5,2)-difference sets in Z_11...")
diff_sets = []
for combo in combinations(range(1, 11), 5):
    D = set(combo)
    # Check: every nonzero element of Z_11 appears exactly 2 times as differences
    diff_counts = {}
    for d1 in D:
        for d2 in D:
            if d1 != d2:
                diff = (d1 - d2) % 11
                diff_counts[diff] = diff_counts.get(diff, 0) + 1
    if all(diff_counts.get(i, 0) == 2 for i in range(1, 11)):
        diff_sets.append(D)

print(f"Found {len(diff_sets)} (11,5,2)-difference sets")

# Group by isomorphism (two diff sets give isomorphic tournaments if one
# is obtained from the other by multiplication by a constant mod 11)
iso_classes = []
used = set()
for i, D in enumerate(diff_sets):
    if i in used:
        continue
    cls = [D]
    used.add(i)
    # Apply multipliers and translations
    for m in range(1, 11):
        D2 = frozenset((m * d) % 11 for d in D)
        for j, D3 in enumerate(diff_sets):
            if j not in used and frozenset(D3) == D2:
                cls.append(D3)
                used.add(j)
    iso_classes.append(cls)

print(f"Isomorphism classes of difference sets: {len(iso_classes)}")

# Build and analyze one tournament per class
print(f"\n{'Class':<8} {'DiffSet':<25} {'DRT?':<6} {'c3':<8} {'c5':<8} {'H':<10}")
print("-" * 70)

for idx, cls in enumerate(iso_classes):
    D = cls[0]
    T = build_drt_from_diff_set(11, D)
    drt = is_doubly_regular(T, 11)
    c3 = count_directed_3cycles(T, 11)
    c5 = count_directed_5cycles(T, 11)

    # Only compute H for DRTs (expensive)
    if drt:
        H = hamiltonian_path_count(T)
        print(f"  {idx:<6} {sorted(D)!s:<25} {'YES':<6} {c3:<8} {c5:<8} {H:<10}")
    else:
        print(f"  {idx:<6} {sorted(D)!s:<25} {'no':<6} {c3:<8} {c5:<8} {'---':<10}")

# Also check: are non-difference-set DRTs possible?
# At n=11 (prime), all DRTs should come from difference sets due to
# the cyclic automorphism group of Z_11.

print(f"\nPaley H(T_11) = {hamiltonian_path_count(T_paley)}")

# Verify Savchenko's claim: c_k is invariant across all DRTs
print("\n--- Savchenko invariance check ---")
print("If all DRTs have the same c3 and c5, cycle counts are DRT invariants.")

print("\nDone.")
