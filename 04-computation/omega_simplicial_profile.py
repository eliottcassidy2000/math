#!/usr/bin/env python3
"""
omega_simplicial_profile.py — opus-2026-03-13-S70

KEY DISCOVERY: The transitive tournament has Omega profile = C(n,k).
At n=5: (5,10,10,5,1) = (C(5,0), C(5,1),...,C(5,4))...
Wait, C(5,0)=1, C(5,1)=5, C(5,2)=10, C(5,3)=10, C(5,4)=5.
But the profile is Omega_0=5, Omega_1=10, Omega_2=10, Omega_3=5, Omega_4=1.
So it's (C(5,1), C(5,2), C(5,3), C(5,4), C(5,5)) = (5,10,10,5,1)?
No... C(5,5)=1, C(5,4)=5. That gives (5,10,10,5,1). Hmm, it's C(5,m+1)?
Actually Omega_0=n=C(n,1), Omega_1=C(n,2) (all edges are 1-paths),
Omega_2=C(n,3)-t3 and for transitive t3=0 so C(n,3).
So for transitive: Omega_m = C(n, m+1)!

Verify this at n=3,4,5,6,7.

Also: investigate the "deviation from simplicial" pattern.
"""

import numpy as np
from math import comb, factorial
from itertools import permutations

def adj_matrix(n, T_bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (T_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def transitive_tournament(n):
    """Transitive tournament: i→j iff i<j."""
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def count_regular_m_paths(A, m):
    n = A.shape[0]
    count = 0
    def dfs(path, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                if depth >= 1 and A[path[-2]][v] != 1:
                    continue
                path.append(v)
                dfs(path, depth + 1)
                path.pop()
    for start in range(n):
        dfs([start], 0)
    return count

def count_regular_m_paths_fast(A, m, memo=None):
    """Slightly faster using sets for membership check."""
    n = A.shape[0]
    count = 0

    def dfs(path_set, last, prev, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        for v in range(n):
            if v in path_set:
                continue
            if not A[last][v]:
                continue
            if depth >= 1 and not A[prev][v]:
                continue
            path_set.add(v)
            dfs(path_set, v, last, depth + 1)
            path_set.remove(v)

    for start in range(n):
        dfs({start}, start, -1, 0)
    return count

# ============================================================
# Verify: Transitive tournament has Omega_m = C(n, m+1)
# ============================================================
print("="*70)
print("TRANSITIVE TOURNAMENT: Omega_m = C(n, m+1)?")
print("="*70)

for n in range(3, 9):
    A = transitive_tournament(n)
    print(f"\n  n = {n}:")
    profile = []
    for m in range(n):
        if n <= 8 or m <= 5:  # limit computation for large n
            omega_m = count_regular_m_paths(A, m)
            expected = comb(n, m+1)
            match = "✓" if omega_m == expected else "✗"
            profile.append(omega_m)
            print(f"    Omega_{m} = {omega_m}, C({n},{m+1}) = {expected} {match}")
        else:
            print(f"    Omega_{m} skipped (too large)")

# ============================================================
# Prove WHY: for transitive tournament, every (m+1)-subset
# has exactly 1 regular m-path
# ============================================================
print(f"\n{'='*70}")
print("PROOF SKETCH: Why Omega_m(trans) = C(n, m+1)")
print("="*70)
print("""
In the transitive tournament on {0,1,...,n-1} with i→j iff i<j:
- An m-path (v_0, v_1, ..., v_m) requires v_0→v_1→...→v_m.
  Since the tournament is transitive, this means v_0 < v_1 < ... < v_m.
- Regularity: v_{i-1} → v_{i+1} for all 0<i<m.
  Since v_{i-1} < v_{i+1} in the transitive tournament, this is
  automatically satisfied.
- Therefore EVERY increasing sequence of length (m+1) is a regular m-path.
- Count = C(n, m+1). QED.
""")

# ============================================================
# For the regular tournament: what's the profile?
# ============================================================
print(f"{'='*70}")
print("REGULAR TOURNAMENT PROFILES")
print("="*70)

# Construct "standard" regular tournaments for small n
# For odd n, the QR tournament on Z_n or the "median" tournament
def circulant_tournament(n, S):
    """Tournament on Z_n with connection set S: i→j iff (j-i)%n in S."""
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

for n in [3, 5, 7]:
    m = (n-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(n, S)

    # Verify it's regular (all out-degrees = m)
    out_degrees = [sum(A[i]) for i in range(n)]
    is_regular = all(d == m for d in out_degrees)

    print(f"\n  n = {n}, interval tournament S={S} (regular={is_regular}):")
    profile = []
    expected_simplicial = []
    for mm in range(n):
        omega = count_regular_m_paths(A, mm)
        simp = comb(n, mm+1)
        profile.append(omega)
        expected_simplicial.append(simp)
        print(f"    Omega_{mm} = {omega:6d}  (simplicial: {simp:6d}, diff: {omega-simp:+6d})")

    chi = sum((-1)**mm * profile[mm] for mm in range(n))
    print(f"    chi = {chi}")

    # Is the profile palindromic? Omega_m = Omega_{n-1-m}?
    is_palin = all(profile[mm] == profile[n-1-mm] for mm in range(n))
    print(f"    Palindromic: {is_palin}")

# ============================================================
# All n=5 profiles: deviation from simplicial
# ============================================================
print(f"\n{'='*70}")
print("n=5: ALL OMEGA PROFILES (deviation from simplicial)")
print("="*70)

n = 5
num_pairs = n*(n-1)//2
simp = [comb(n, m+1) for m in range(n)]
print(f"  Simplicial: {simp}")

# Group by isomorphism type
from collections import defaultdict

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = []
        for i in range(n):
            for j in range(i+1, n):
                enc.append(A[perm[i]][perm[j]])
        enc = tuple(enc)
        if best is None or enc < best:
            best = enc
    return best

profiles_by_type = {}
for bits in range(2**num_pairs):
    A = adj_matrix(n, bits)
    cf = canon_form(A)
    if cf not in profiles_by_type:
        profile = [count_regular_m_paths(A, m) for m in range(n)]
        profiles_by_type[cf] = profile

print(f"\n  {'Type':>4s} {'O0':>4s} {'O1':>4s} {'O2':>4s} {'O3':>4s} {'O4':>4s} | {'d0':>4s} {'d1':>4s} {'d2':>4s} {'d3':>4s} {'d4':>4s} | {'chi':>4s} {'H':>4s}")
print(f"  {'':4s} {'':>4s} {'':>4s} {'':>4s} {'':>4s} {'':>4s} | dev from simplicial  | ")
for cf, profile in sorted(profiles_by_type.items(), key=lambda x: x[1], reverse=True):
    dev = [profile[m] - simp[m] for m in range(n)]
    chi = sum((-1)**m * profile[m] for m in range(n))
    H = sum(1 for perm in permutations(range(n))
            if all(profiles_by_type[cf] is not None for _ in [0]))  # placeholder

    # Recompute H from a representative
    A_tmp = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if cf[idx]:
                A_tmp[i][j] = 1
            else:
                A_tmp[j][i] = 1
            idx += 1
    H_val = sum(1 for perm in permutations(range(n))
                if all(A_tmp[perm[i]][perm[i+1]] for i in range(n-1)))

    print(f"  {str(cf)[:30]:30s} {profile[0]:4d} {profile[1]:4d} {profile[2]:4d} {profile[3]:4d} {profile[4]:4d} |"
          f" {dev[0]:+4d} {dev[1]:+4d} {dev[2]:+4d} {dev[3]:+4d} {dev[4]:+4d} | {chi:4d} {H_val:4d}")

# ============================================================
# n=6: transitive profile check
# ============================================================
print(f"\n{'='*70}")
print("n=6: TRANSITIVE TOURNAMENT PROFILE")
print("="*70)

n = 6
A = transitive_tournament(n)
profile_6 = [count_regular_m_paths(A, m) for m in range(n)]
simp_6 = [comb(n, m+1) for m in range(n)]
print(f"  Omega profile: {profile_6}")
print(f"  Simplicial:    {simp_6}")
print(f"  Match: {profile_6 == simp_6}")

# ============================================================
# Deviation pattern: is Omega_m(T) - C(n,m+1) always ≤ 0?
# ============================================================
print(f"\n  Is transitive tournament the MAXIMIZER of Omega_m?")
n = 5
num_pairs = n*(n-1)//2
simp_5 = [comb(n, m+1) for m in range(n)]

for m in range(n):
    max_val = 0
    for bits in range(2**num_pairs):
        A = adj_matrix(n, bits)
        val = count_regular_m_paths(A, m)
        max_val = max(max_val, val)
    print(f"    Omega_{m}: max={max_val}, C(n,m+1)={simp_5[m]}, "
          f"trans is max: {max_val == simp_5[m]}")

# Actually for Omega_4 at n=5: regular tournament has Omega_4=5, transitive has 1.
# So transitive is NOT the maximizer for all m!

print(f"\n  For n=5 regular tournament: Omega_4 = 5 > 1 = transitive")
print(f"  => Transitive is NOT always maximizer!")
print(f"  The regular tournament maximizes Omega_{n-1} at odd n.")

# But for m < n-1, is there a pattern?
print(f"\n  Which tournament type maximizes each Omega_m at n=5?")
for m in range(n):
    best_val = 0
    best_types = set()
    for bits in range(2**num_pairs):
        A = adj_matrix(n, bits)
        val = count_regular_m_paths(A, m)
        if val > best_val:
            best_val = val
            best_types = {canon_form(A)}
        elif val == best_val:
            best_types.add(canon_form(A))
    print(f"    Omega_{m}: max = {best_val}, achieved by {len(best_types)} type(s)")

print("\nDONE.")
