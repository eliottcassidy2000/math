#!/usr/bin/env python3
"""
average_chi_per.py — opus-2026-03-12-S69

Compute the average chi_per over ALL circulant tournament orientations.

Σ_T chi_per(T) where T ranges over all 2^m orientations of C_p.

This is a "character average" — if it has a closed form, it would
connect spectral theory to topology.
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

for p in [5, 7, 11]:
    m = (p-1)//2
    print(f"\np = {p}, m = {m}, 2^m = {2**m}")

    total_chi = 0
    chi_distribution = {}
    count = 0

    for S in gen_orientations(p):
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))
        total_chi += chi_per
        chi_distribution[chi_per] = chi_distribution.get(chi_per, 0) + 1
        count += 1

    print(f"  Σ chi_per = {total_chi}")
    print(f"  Average chi_per = {total_chi}/{count} = {total_chi/count}")
    print(f"  Distribution: {dict(sorted(chi_distribution.items()))}")

    # Also compute Σ Ω_m for each m
    total_omega_by_m = [0] * p
    for S in gen_orientations(p):
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        for mm in range(p):
            total_omega_by_m[mm] += omega[mm]

    print(f"  Σ Ω_m by degree: {total_omega_by_m}")
    avg_omega = [t/count for t in total_omega_by_m]
    print(f"  Avg Ω_m: [{', '.join(f'{a:.2f}' for a in avg_omega)}]")
    chi_from_avg = sum((-1)**mm * avg_omega[mm] for mm in range(p))
    print(f"  Σ(-1)^m <Ω_m> = {chi_from_avg:.4f} (should equal avg chi_per = {total_chi/count:.4f})")

print("\nDONE.")
