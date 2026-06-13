#!/usr/bin/env python3
"""
circulant_omega_profiles.py — opus-2026-03-13-S70

Investigate Omega profiles for various circulant tournaments.

KEY OBSERVATION from previous run:
  n=7 interval: Omega = (7, 21, 21, 35, 56, 49, 14)
  - Omega_3 = 35 = C(7,4) = simplicial value
  - Omega_4 = 56 = C(8,5) (??)
  - Omega_5 = 49 = 7^2
  - Omega_6 = 14 = 2*7

Questions:
1. What determines the Omega profile of circulant tournaments?
2. Is there a pattern involving Q_k (Fourier magnitudes)?
3. What is chi for interval tournaments as a function of p?
"""

import numpy as np
from math import comb, gcd
from itertools import permutations

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def count_regular_m_paths(A, m):
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
# Interval tournament profiles for odd n
# ============================================================
print("="*70)
print("INTERVAL TOURNAMENT Omega PROFILES")
print("="*70)

for n in [3, 5, 7, 9]:
    m = (n-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(n, S)

    print(f"\n  n = {n}, S = {{1,...,{m}}}:")
    profile = []
    for mm in range(n):
        omega = count_regular_m_paths(A, mm)
        profile.append(omega)

    simp = [comb(n, mm+1) for mm in range(n)]
    chi = sum((-1)**mm * profile[mm] for mm in range(n))

    print(f"    Omega: {profile}")
    print(f"    Simp:  {simp}")
    print(f"    chi = {chi}")

    # Check for patterns
    # Ratios
    for mm in range(n):
        if simp[mm] > 0:
            ratio = profile[mm] / simp[mm]
            print(f"    Omega_{mm}/C(n,{mm+1}) = {ratio:.4f}")

# ============================================================
# QR tournament profiles (p prime)
# ============================================================
print(f"\n{'='*70}")
print("QR TOURNAMENT Omega PROFILES")
print("="*70)

def legendre(a, p):
    """Legendre symbol."""
    if a % p == 0:
        return 0
    return pow(a, (p-1)//2, p)

for p in [5, 7, 11, 13]:
    QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
    A = circulant_tournament(p, QR)

    print(f"\n  p = {p}, QR = {sorted(QR)}:")
    profile = []
    for mm in range(p):
        if p <= 11 or mm <= 6:
            omega = count_regular_m_paths(A, mm)
            profile.append(omega)
        else:
            profile.append(None)

    print(f"    Omega: {[x for x in profile if x is not None]}")
    chi = sum((-1)**mm * profile[mm] for mm in range(len(profile)) if profile[mm] is not None)
    print(f"    chi (partial) = {chi}")

# ============================================================
# Other circulant sets
# ============================================================
print(f"\n{'='*70}")
print("OTHER CIRCULANT TOURNAMENT PROFILES (n=7)")
print("="*70)

n = 7
for S_tuple in [(1,2,3), (1,3,5), (1,2,4), (2,3,4)]:
    S = set(S_tuple)
    # Check if S ∪ (n-S) = {1,...,n-1}
    complement = {n - s for s in S}
    if S & complement:
        print(f"\n  S={S}: not a valid tournament set (S ∩ -S ≠ ∅)")
        continue
    if S | complement != set(range(1, n)):
        print(f"\n  S={S}: not a valid tournament set (S ∪ -S ≠ {{1,...,{n-1}}})")
        continue

    A = circulant_tournament(n, S)
    print(f"\n  S = {sorted(S)}:")
    profile = []
    for mm in range(n):
        omega = count_regular_m_paths(A, mm)
        profile.append(omega)

    chi = sum((-1)**mm * profile[mm] for mm in range(n))
    print(f"    Omega: {profile}")
    print(f"    chi = {chi}")

# ============================================================
# n=9 interval tournament
# ============================================================
print(f"\n{'='*70}")
print("n=9 INTERVAL TOURNAMENT (S={1,2,3,4})")
print("="*70)

n = 9
m_half = 4
S = set(range(1, m_half+1))
A = circulant_tournament(n, S)

print(f"  Computing Omega profile for n={n}...")
profile = []
for mm in range(n):
    print(f"    Computing Omega_{mm}...", end=" ", flush=True)
    omega = count_regular_m_paths(A, mm)
    profile.append(omega)
    print(f"{omega}")

simp = [comb(n, mm+1) for mm in range(n)]
chi = sum((-1)**mm * profile[mm] for mm in range(n))

print(f"\n  Omega: {profile}")
print(f"  Simp:  {simp}")
print(f"  chi = {chi}")

# Check ratios
for mm in range(n):
    ratio = profile[mm] / simp[mm] if simp[mm] > 0 else float('inf')
    print(f"  Omega_{mm}/C({n},{mm+1}) = {ratio:.4f} (diff={profile[mm]-simp[mm]:+d})")

print("\nDONE.")
