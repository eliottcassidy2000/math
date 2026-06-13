#!/usr/bin/env python3
"""
paley_omega_structure.py — opus-2026-03-12-S69

The Paley tournament (QR) at p≡3 mod 4 has FLAT Q spectrum:
  Q_k = (p+1)/4 for all k = 1,...,m.

This maximal symmetry should give a particularly nice Ω profile.

Questions:
1. Is there a closed form for Ω_m of the Paley tournament?
2. Does the Ω generating function relate to the characteristic polynomial?
3. What's the relationship between flat Q and palindromic Ω?
"""

import sys
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        yield S

def zpstar_orbit_rep(p, S):
    best = sorted(S)
    for a in range(2, p):
        S_a = frozenset((a * s) % p for s in S)
        if sorted(S_a) < best:
            best = sorted(S_a)
    return tuple(best)

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

# For p ≡ 3 mod 4: Paley = QR
primes_3mod4 = [7, 11, 19, 23]

print("="*80)
print("PALEY Ω STRUCTURE at p ≡ 3 mod 4")
print("="*80)

for p in primes_3mod4:
    m = (p-1)//2
    QR = frozenset(d for d in range(1,p) if legendre(d,p)==1)

    if len(QR) != m:
        print(f"\np={p}: QR has {len(QR)} elements, need {m}")
        continue

    # Verify it's a valid tournament orientation
    valid = all((d in QR) != ((p-d) in QR) for d in range(1, m+1))
    if not valid:
        print(f"\np={p}: QR is NOT a valid tournament orientation")
        continue

    Q_val = (p+1)/4
    print(f"\np = {p} (m = {m})")
    print(f"  QR = {sorted(QR)}")
    print(f"  Q_k = (p+1)/4 = {Q_val} for all k")

    if p <= 19:
        h = CirculantHomology(n=p, S=QR)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))
        total = sum(omega)

        print(f"  Ω = {omega}")
        print(f"  chi_per = {chi_per}")
        print(f"  total Ω = {total}")

        # Check palindromicity
        n = len(omega)
        is_palin = all(omega[i] == omega[n-1-i] for i in range(n//2))
        print(f"  Palindromic: {is_palin}")

        # Ratios Ω_{k+1} / Ω_k
        ratios = [omega[k+1]/omega[k] if omega[k] > 0 else 'inf' for k in range(len(omega)-1)]
        print(f"  Ratios: [{', '.join(f'{r:.4f}' if isinstance(r, float) else r for r in ratios)}]")

        # Check if Ω_k = C(m,k) * c^k for some c
        # Or Ω_k = f(k) * Q_val^{something}
        print(f"  Ω_k / C(m,k):")
        from math import comb
        for k in range(len(omega)):
            ck = comb(m, k) if k <= m else 0
            if ck > 0:
                print(f"    k={k}: Ω/C({m},{k}) = {omega[k]}/{ck} = {omega[k]/ck:.4f}")

# Also compute Interval Ω for comparison
print(f"\n{'='*80}")
print(f"INTERVAL vs PALEY Ω COMPARISON")
print(f"{'='*80}")

for p in [7, 11]:
    m = (p-1)//2
    S_int = frozenset(range(1, m+1))
    QR = frozenset(d for d in range(1,p) if legendre(d,p)==1)

    h_int = CirculantHomology(n=p, S=S_int)
    h_int._ensure_enumerated(p)
    omega_int = h_int.omega_dims(max_degree=p-1)

    h_pal = CirculantHomology(n=p, S=QR)
    h_pal._ensure_enumerated(p)
    omega_pal = h_pal.omega_dims(max_degree=p-1)

    print(f"\np = {p}:")
    print(f"  Interval Ω = {omega_int}, chi_per = {sum((-1)**k*w for k,w in enumerate(omega_int))}")
    print(f"  Paley   Ω = {omega_pal}, chi_per = {sum((-1)**k*w for k,w in enumerate(omega_pal))}")
    print(f"  Difference: {[omega_pal[k]-omega_int[k] for k in range(p)]}")

# Check: does chi_per = 1 for ALL Paley tournaments at p ≡ 3 mod 4?
print(f"\n{'='*80}")
print(f"chi_per = 1 FOR PALEY AT p ≡ 3 mod 4?")
print(f"{'='*80}")

# Also check: chi_per = 0 for Interval at p ≡ 3 mod 4?
for p in [7, 11, 19]:
    m = (p-1)//2
    S_int = frozenset(range(1, m+1))
    QR = frozenset(d for d in range(1,p) if legendre(d,p)==1)

    # Interval
    h_int = CirculantHomology(n=p, S=S_int)
    h_int._ensure_enumerated(p)
    omega_int = h_int.omega_dims(max_degree=p-1)
    chi_int = sum((-1)**k*w for k,w in enumerate(omega_int))

    # Paley
    h_pal = CirculantHomology(n=p, S=QR)
    h_pal._ensure_enumerated(p)
    omega_pal = h_pal.omega_dims(max_degree=p-1)
    chi_pal = sum((-1)**k*w for k,w in enumerate(omega_pal))

    print(f"  p={p}: Interval chi_per = {chi_int}, Paley chi_per = {chi_pal}")

print("\nDONE.")
