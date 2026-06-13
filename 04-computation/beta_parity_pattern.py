#!/usr/bin/env python3
"""
beta_parity_pattern.py - Investigate even/odd Betti number vanishing pattern

Observation from n=7 data:
- beta_2 = 0 ALWAYS (conjecture)
- beta_4 = 0 for random n=7 (but Paley T_7 has beta_4=6!)
- beta_3 can be nonzero (8.2% of n=7)
- beta_5 can be nonzero

Questions:
1. At n=7: is beta_4=0 for ALL non-Paley? Or just rare?
2. At n=8,9: which beta_p can be nonzero?
3. Is there a pattern: beta_{even} = 0 for p >= 2?
4. Does the vanishing depend on score sequence / regularity?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def paley_tournament(p):
    """Paley tournament on p vertices (p prime, p=3 mod 4)."""
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
                    c3 += 1
    return c3


# ============================================================
# PART 1: n=7 detailed beta survey (larger sample)
# ============================================================
print("=" * 70)
print("n=7: DETAILED BETTI SURVEY (3000 random + Paley)")
print("=" * 70)

n = 7
max_dim = n - 2  # = 5
from collections import Counter, defaultdict

betti_profiles = Counter()
beta4_nonzero = []
beta2_nonzero = 0
t0 = time.time()

for trial in range(3000):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=max_dim)
    # Pad to length max_dim+1
    while len(betti) <= max_dim:
        betti.append(0)
    profile = tuple(int(b) for b in betti)
    betti_profiles[profile] += 1

    if profile[2] != 0:
        beta2_nonzero += 1
    if profile[4] != 0:
        scores = tuple(sorted([sum(row) for row in A]))
        c3 = count_3cycles(A, n)
        beta4_nonzero.append((profile, scores, c3))

elapsed = time.time() - t0
print(f"  3000 random n=7 tournaments ({elapsed:.0f}s)")
print(f"  beta_2 != 0: {beta2_nonzero}")
print(f"  beta_4 != 0: {len(beta4_nonzero)}")

# Add Paley
A_paley = paley_tournament(7)
betti_paley = path_betti_numbers(A_paley, n, max_dim=max_dim)
while len(betti_paley) <= max_dim:
    betti_paley.append(0)
profile_paley = tuple(int(b) for b in betti_paley)
print(f"  Paley T_7: {profile_paley}")

print(f"\n  ALL Betti profiles (n=7):")
for profile, count in sorted(betti_profiles.items(), key=lambda x: -x[1]):
    parity = ""
    for p in range(max_dim + 1):
        if profile[p] != 0:
            parity += f"b{p}={profile[p]} "
    print(f"    {profile}: {count} ({100*count/3000:.1f}%) -- {parity}")

if beta4_nonzero:
    print(f"\n  beta_4 != 0 examples:")
    for profile, scores, c3 in beta4_nonzero[:10]:
        print(f"    beta={profile}, scores={scores}, c3={c3}")


# ============================================================
# PART 2: n=8 beta survey
# ============================================================
print(f"\n{'='*70}")
print("n=8: BETTI SURVEY (500 random)")
print("=" * 70)

n = 8
max_dim = n - 2  # = 6
betti_profiles_8 = Counter()
even_nonzero = defaultdict(int)  # {p: count} for even p >= 2
t0 = time.time()

for trial in range(500):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=max_dim)
    while len(betti) <= max_dim:
        betti.append(0)
    profile = tuple(int(b) for b in betti)
    betti_profiles_8[profile] += 1

    for p in range(2, max_dim + 1, 2):  # even p >= 2
        if profile[p] != 0:
            even_nonzero[p] += 1

elapsed = time.time() - t0
print(f"  500 random n=8 ({elapsed:.0f}s)")

print(f"\n  Even-dimensional nonzero counts:")
for p in sorted(even_nonzero.keys()):
    print(f"    beta_{p} != 0: {even_nonzero[p]}/500")

print(f"\n  ALL Betti profiles (n=8, showing top 15):")
for profile, count in sorted(betti_profiles_8.items(), key=lambda x: -x[1])[:15]:
    parity = ""
    for p in range(max_dim + 1):
        if profile[p] != 0:
            parity += f"b{p}={profile[p]} "
    print(f"    {profile}: {count} ({100*count/500:.1f}%) -- {parity}")

# Check: among beta_3!=0, what's beta_5?
b3_cases = {p: c for p, c in betti_profiles_8.items() if p[3] != 0}
print(f"\n  beta_3 != 0 profiles:")
for profile, count in sorted(b3_cases.items(), key=lambda x: -x[1])[:10]:
    print(f"    {profile}: {count}")


# ============================================================
# PART 3: n=6 exhaustive check of ALL beta_p
# ============================================================
print(f"\n{'='*70}")
print("n=6: EXHAUSTIVE ALL-BETTI CHECK")
print("=" * 70)

n = 6
max_dim = n - 2  # = 4
total = 1 << (n*(n-1)//2)
betti_profiles_6 = Counter()
t0 = time.time()

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=max_dim)
    while len(betti) <= max_dim:
        betti.append(0)
    profile = tuple(int(b) for b in betti)
    betti_profiles_6[profile] += 1

    if (bits+1) % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6 exhaustive ({elapsed:.0f}s):")
for profile, count in sorted(betti_profiles_6.items(), key=lambda x: -x[1]):
    pct = 100*count/total
    parity = ""
    for p in range(max_dim + 1):
        if profile[p] != 0:
            parity += f"b{p}={profile[p]} "
    print(f"  {profile}: {count} ({pct:.1f}%) -- {parity}")

# Key check: ANY even-dimensional nonvanishing at p >= 2?
any_even_2plus = 0
for profile, count in betti_profiles_6.items():
    for p in range(2, max_dim + 1, 2):
        if profile[p] != 0:
            any_even_2plus += count
            break
print(f"\n  ANY beta_p != 0 for even p>=2: {any_even_2plus}/{total}")


# ============================================================
# PART 4: n=5 exhaustive
# ============================================================
print(f"\n{'='*70}")
print("n=5: EXHAUSTIVE ALL-BETTI CHECK")
print("=" * 70)

n = 5
max_dim = n - 2  # = 3
total = 1 << (n*(n-1)//2)
betti_profiles_5 = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=max_dim)
    while len(betti) <= max_dim:
        betti.append(0)
    profile = tuple(int(b) for b in betti)
    betti_profiles_5[profile] += 1

print(f"n=5 exhaustive:")
for profile, count in sorted(betti_profiles_5.items(), key=lambda x: -x[1]):
    pct = 100*count/total
    print(f"  {profile}: {count} ({pct:.2f}%)")


# ============================================================
# PART 5: Euler characteristic formula at n=7
# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC vs BETTI at n=7")
print("=" * 70)

n = 7
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

chi_distribution = Counter()
for trial in range(200):
    A = random_tournament(n)
    # Compute Omega dims
    allowed = {}
    for p in range(-1, n):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega_dims = []
    for p in range(n - 1):
        omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim_omega = omega_p.shape[1] if omega_p.ndim == 2 else 0
        omega_dims.append(dim_omega)

    chi = sum((-1)**p * d for p, d in enumerate(omega_dims))
    betti = path_betti_numbers(A, n, max_dim=n-2)
    while len(betti) <= n-2:
        betti.append(0)
    chi_betti = sum((-1)**p * int(b) for p, b in enumerate(betti))
    b1 = int(betti[1]) if len(betti) > 1 else 0

    chi_distribution[(chi, chi_betti, b1)] += 1

print(f"  (chi_Omega, chi_betti, beta_1) distribution (200 samples):")
for (chi, chi_b, b1), count in sorted(chi_distribution.items(), key=lambda x: -x[1])[:15]:
    print(f"    chi_Omega={chi}, chi_betti={chi_b}, beta_1={b1}: {count}")


print("\n\nDone.")
