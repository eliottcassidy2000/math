#!/usr/bin/env python3
"""
omega_gf_analysis.py — opus-2026-03-12-S69

Analyze the generating function G(t) = Σ Ω_k t^k for various circulant tournaments.

Key question: Does G(t) factor nicely? Is it related to (1+t)^m or other
standard generating functions?

For the TRANSITIVE tournament (on non-circulant n vertices), Ω_k = C(n,k+1).
G(t) = Σ C(n,k+1) t^k = (1+t)^n/t - 1/t = ((1+t)^n - 1)/t.

For circulant tournaments, the per-eigenspace Ω should have a different structure.
"""

import numpy as np
from math import comb, factorial

# Per-eigenspace Ω data (already computed)
omega_data = {
    # p: {orbit_name: [Ω_0, ..., Ω_{p-1}]}
    7: {
        'Interval': [1, 3, 6, 11, 14, 9, 2],
        'Paley':    [1, 3, 6, 9, 9, 6, 3],
    },
    11: {
        'Interval': [1, 5, 20, 74, 224, 522, 926, 1222, 1115, 611, 148],
        'Orbit_A':  [1, 5, 20, 70, 200, 439, 711, 827, 648, 301, 64],
        'Orbit_B':  [1, 5, 20, 70, 201, 430, 620, 596, 384, 151, 26],
        'Paley':    [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30],
    },
    13: {
        'Interval': [1,6,30,140,556,1823,4983,11241,19730,24780,20344,9642,1989],
        'Orbit_1':  [1,6,30,135,517,1622,4129,8436,13292,15078,11316,4988,977],
        'Orbit_2':  [1,6,30,135,524,1675,4176,7707,10295,9681,6081,2298,399],
        'qPaley':   [1,6,30,135,528,1707,4245,7503,9177,7770,4374,1488,237],
        'Orbit_4':  [1,6,30,135,520,1627,3915,6965,9012,8325,5157,1906,321],
        'Orbit_5':  [1,6,30,135,524,1663,4001,6876,8409,7266,4237,1493,238],
    },
}

print("="*80)
print("Ω GENERATING FUNCTION ANALYSIS")
print("="*80)

for p, orbits in omega_data.items():
    m = (p-1)//2
    print(f"\n{'─'*60}")
    print(f"p = {p}, m = {m}")
    print(f"{'─'*60}")

    for name, omega in orbits.items():
        total = sum(omega)
        chi = sum((-1)**k * w for k, w in enumerate(omega))

        print(f"\n  {name}: Ω = {omega}")
        print(f"    Total = {total}, chi_per = {chi}")

        # G(t) at various t values
        for t_val in [1, -1, 2, 0.5]:
            G = sum(w * t_val**k for k, w in enumerate(omega))
            print(f"    G({t_val}) = {G:.4f}", end="")
            # Compare with (1+t)^m
            ref = (1 + t_val)**m
            print(f"  (1+t)^m = {ref:.4f}, ratio = {G/ref:.4f}" if ref != 0 else "")

        # G(1) = total Ω, G(-1) = chi_per
        # G(t)/G_Paley(t) ratio
        if 'Paley' in orbits and name != 'Paley':
            omega_p = orbits['Paley']
            print(f"    Ω/Ω_Paley ratios: [{', '.join(f'{omega[k]/omega_p[k]:.3f}' if omega_p[k] > 0 else '?' for k in range(p))}]")

    # For Paley specifically: check G(t) / (1+ct)^m where c = (p+1)/4
    if 'Paley' in orbits:
        omega_p = orbits['Paley']
        c = (p+1)/4
        print(f"\n  PALEY G(t) vs (1+{c}t)^m:")
        for t_val in [0.1, 0.2, 0.5, 1.0, 2.0]:
            G = sum(w * t_val**k for k, w in enumerate(omega_p))
            ref = (1 + c*t_val)**m
            print(f"    t={t_val}: G = {G:.4f}, (1+{c}t)^m = {ref:.4f}, ratio = {G/ref:.4f}")

    # For Interval: check G(t) / (some function)
    omega_int = orbits.get('Interval', list(orbits.values())[0])
    print(f"\n  INTERVAL G(t) analysis:")
    # The transitive tournament on p vertices has Ω_k = C(p, k+1).
    # But circulant Interval is NOT the transitive tournament.
    # For comparison: C(p, k+1) for k=0,...,p-1
    trans_omega = [comb(p, k+1) for k in range(p)]
    print(f"    Transitive Ω = {trans_omega}")
    print(f"    Interval  Ω = {omega_int}")
    print(f"    Ratio = [{', '.join(f'{omega_int[k]/trans_omega[k]:.4f}' if trans_omega[k] > 0 else '?' for k in range(p))}]")

# Look for the "excess" Ω compared to the Paley base
print(f"\n{'='*80}")
print(f"INTERVAL EXCESS OVER PALEY")
print(f"{'='*80}")

for p in [7, 11]:
    m = (p-1)//2
    omega_int = omega_data[p]['Interval']
    omega_pal = omega_data[p]['Paley']
    excess = [omega_int[k] - omega_pal[k] for k in range(p)]
    print(f"\np = {p}: Excess = {excess}")
    print(f"  Sum of excess = {sum(excess)}")
    print(f"  Alt sum = {sum((-1)**k * e for k, e in enumerate(excess))}")

    # Is the excess a multiple of some simple sequence?
    # Check if excess[k] / C(m,k) has a pattern
    print(f"  Excess / Ω_Paley:")
    for k in range(p):
        if omega_pal[k] > 0:
            print(f"    k={k}: excess/Ω_P = {excess[k]/omega_pal[k]:.6f}")

# Finally: check the Ω profile symmetry
print(f"\n{'='*80}")
print(f"Ω PROFILE SYMMETRY: Ω_k vs Ω_{p-1-k}")
print(f"{'='*80}")

for p in [7, 11, 13]:
    print(f"\np = {p}:")
    for name, omega in omega_data[p].items():
        # Compare Ω_k and Ω_{p-1-k}
        asym = max(abs(omega[k] - omega[p-1-k]) for k in range(p//2))
        sym_ratio = [omega[k]/omega[p-1-k] if omega[p-1-k] > 0 else float('inf') for k in range(p//2)]
        center = omega[p//2] if p % 2 == 1 else None
        print(f"  {name}: max asymmetry = {asym}, center Ω_{p//2} = {center}")

print("\nDONE.")
