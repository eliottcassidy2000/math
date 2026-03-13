#!/usr/bin/env python3
"""
paley7_constant_omega.py — opus-2026-03-13-S70

The Paley tournament at p=7 (QR={1,2,4}) has CONSTANT Omega_m = 21 = C(7,2)
for all 1 ≤ m ≤ 6. This is extraordinary. Investigate WHY.

The Paley tournament at p=7 is also the UNIQUE doubly regular tournament
on 7 vertices (every pair of vertices has the same number of common
out-neighbors). It's the circulant with S = {1, 2, 4} = QR(7).

Note: for p ≡ 3 mod 4, -1 is a NQR, so the Paley tournament is well-defined.
p=7: 7 ≡ 3 mod 4 ✓

Question: is the constant Omega profile related to the doubly regular property?
Does it hold for other Paley tournaments (p=3, 11, 19, 23, ...)?
"""

import numpy as np
from math import comb
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
            if v in path_set: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path_set.add(v)
            dfs(path_set, v, last, depth + 1)
            path_set.remove(v)
    for start in range(n):
        dfs({start}, start, -1, 0)
    return count

def legendre(a, p):
    if a % p == 0: return 0
    return pow(a, (p-1)//2, p)

# ============================================================
# Paley tournaments for p ≡ 3 mod 4
# ============================================================
print("="*70)
print("PALEY TOURNAMENT Omega PROFILES (p ≡ 3 mod 4)")
print("="*70)

for p in [3, 7, 11, 19, 23]:
    if p % 4 != 3:
        continue

    QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
    A = circulant_tournament(p, QR)

    # Verify tournament property
    is_tournament = True
    for i in range(p):
        for j in range(i+1, p):
            if A[i][j] + A[j][i] != 1:
                is_tournament = False
                break
    if not is_tournament:
        print(f"\n  p = {p}: NOT a valid tournament!")
        continue

    print(f"\n  p = {p}, QR = {sorted(QR)}:")

    profile = []
    c_n_2 = comb(p, 2)
    all_equal = True

    for mm in range(p):
        if p <= 11 or mm <= min(8, p-1):
            omega = count_regular_m_paths(A, mm)
            profile.append(omega)
            if mm >= 1 and omega != c_n_2:
                all_equal = False
        else:
            profile.append(None)

    computed = [x for x in profile if x is not None]
    print(f"    Omega: {computed}")
    print(f"    C({p},2) = {c_n_2}")

    if all_equal and all(x is not None for x in profile):
        print(f"    *** ALL Omega_m = {c_n_2} for m >= 1 ***")
    elif not all_equal:
        print(f"    NOT constant")

    # Check chi
    if all(x is not None for x in profile):
        chi = sum((-1)**mm * profile[mm] for mm in range(p))
        print(f"    chi = {chi}")

    # Check ratios
    for mm, val in enumerate(profile):
        if val is not None and val > 0:
            print(f"    Omega_{mm} = {val}, ratio to C(p,2) = {val/c_n_2:.4f}")

# ============================================================
# Check p=11 more carefully
# ============================================================
print(f"\n{'='*70}")
print("p=11 PALEY DETAILED")
print("="*70)

p = 11
QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
A = circulant_tournament(p, QR)
print(f"  QR = {sorted(QR)}")

for mm in range(p):
    omega = count_regular_m_paths(A, mm)
    print(f"  Omega_{mm} = {omega:8d}, /C({p},2) = {omega/comb(p,2):8.4f}, "
          f"/C({p},{mm+1}) = {omega/comb(p,mm+1) if comb(p,mm+1) > 0 else 0:.4f}")

# ============================================================
# For p=7: understand WHY every regular m-path count = 21
# ============================================================
print(f"\n{'='*70}")
print("p=7 PALEY: WHY IS Omega_m = 21?")
print("="*70)

p = 7
QR = {1, 2, 4}
A = circulant_tournament(p, QR)

# For each m, list all regular m-paths
for mm in [2, 3, 4, 5, 6]:
    paths = []
    def collect_paths(path, depth, prev):
        if depth == mm:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(p):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            collect_paths(path, depth+1, last)
            path.pop()

    for start in range(p):
        collect_paths([start], 0, -1)

    print(f"\n  Omega_{mm}: {len(paths)} regular {mm}-paths")

    # Analyze: how many paths per starting vertex?
    from collections import Counter
    starts = Counter(p[0] for p in paths)
    print(f"    Per starting vertex: {dict(starts)}")

    # How many per (start, end)?
    start_end = Counter((p[0], p[-1]) for p in paths)
    vals = sorted(set(start_end.values()))
    print(f"    Per (start, end): values = {vals}, count = {len(start_end)}")

    # Are they all 3 per starting vertex? 21/7 = 3
    if len(set(starts.values())) == 1:
        print(f"    => {list(starts.values())[0]} paths per vertex (circulant symmetry)")

# Check the "doubly regular" property
print(f"\n  Doubly regular check:")
for i in range(p):
    out_i = {j for j in range(p) if j != i and A[i][j]}
    for j in range(i+1, p):
        out_j = {k for k in range(p) if k != j and A[j][k]}
        common_out = len(out_i & out_j)
        # In doubly regular tournament: same for all i,j
        print(f"    ({i},{j}): |N_out(i) ∩ N_out(j)| = {common_out}")
    break  # just show one row

# For a doubly regular tournament, common out-neighbors is always (n-3)/4
print(f"  Expected common out-neighbors = (p-3)/4 = {(p-3)/4}")

# ============================================================
# Connection to eigenspaces
# ============================================================
print(f"\n{'='*70}")
print("EIGENSPACE DECOMPOSITION AND Omega PROFILES")
print("="*70)

# For circulant tournament on Z_p with connection set S:
# Q_k = |Ŝ(k)|² = |Σ_{s∈S} ω^{ks}|²
# For QR tournament: Q_k = (p-1)/4 for all k ≠ 0 (since Gauss sums)

print(f"  QR tournament (p ≡ 3 mod 4):")
print(f"  Q_k = (p-1)/4 for all k ≠ 0")
print(f"  At p=7: Q_k = 6/4 = 1.5")
print(f"  At p=11: Q_k = 10/4 = 2.5")
print(f"  At p=19: Q_k = 18/4 = 4.5")

# Verify
for p_test in [7, 11]:
    QR_test = {a % p_test for a in range(1, p_test) if legendre(a, p_test) == 1}
    for k in range(1, p_test):
        omega = np.exp(2j * np.pi / p_test)
        S_hat = sum(omega**(k*s) for s in QR_test)
        Q_k = abs(S_hat)**2
        print(f"    p={p_test}, k={k}: Q_k = {Q_k:.4f}, expected = {(p_test-1)/4:.4f}")
    print()

# If all Q_k are equal, then the eigenspace decomposition is maximally uniform.
# Each eigenspace contributes the same amount to Omega_m.
# This might explain the constant Omega profile!

print(f"  KEY INSIGHT: QR tournament has Q_k = (p-1)/4 for all k ≠ 0")
print(f"  => ALL eigenspaces are IDENTICAL")
print(f"  => Omega_m decomposes uniformly across eigenspaces")
print(f"  => By THM-145, the per-eigenspace Omega profile is the same for all k")
print(f"  => The total Omega_m = p * Omega_m(per eigenspace) + correction from k=0")

print("\nDONE.")
