#!/usr/bin/env python3
"""
chi_per_sign_mechanism.py — opus-2026-03-12-S69

What determines the SIGN of chi_per for circulant tournaments?

At p=13:
  chi_per = -17: S={1,2,3,5,6,9} (cubic-residue coset)
  chi_per = -8:  S={1,2,3,5,7,9}
  chi_per = -3:  S={1,2,3,4,5,7} (almost contiguous)
  chi_per = +1:  S={1,2,3,4,5,6} (Interval)
  chi_per = +1:  S={1,2,6,8,9,10}
  chi_per = +4:  S={1,2,3,4,6,8}

Hypothesis: The sign is determined by whether the Ω profile is
"right-skewed" (mode at high m → positive chi) or
"left-skewed" (mode at low m → negative chi).

Also test: parity of QR content, additive energy, sumset structure.
"""

import sys
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

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
    return frozenset(best)

p = 13
m = (p-1)//2

# Get all orbit reps
orbits = {}
for S in gen_orientations(p):
    rep = zpstar_orbit_rep(p, S)
    if rep not in orbits:
        orbits[rep] = S

print(f"{'='*80}")
print(f"CHI_PER SIGN MECHANISM at p = {p}")
print(f"{'='*80}")

results = []
for rep, S in orbits.items():
    h = CirculantHomology(n=p, S=S)
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)
    chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))

    # Ω profile analysis
    total_omega = sum(omega)
    mode_m = max(range(len(omega)), key=lambda mm: omega[mm])
    centroid = sum(mm * w for mm, w in enumerate(omega)) / total_omega

    # Skewness
    var_m = sum((mm - centroid)**2 * w for mm, w in enumerate(omega)) / total_omega
    std_m = var_m**0.5
    skew_m = sum((mm - centroid)**3 * w for mm, w in enumerate(omega)) / (total_omega * std_m**3) if std_m > 0 else 0

    # QR content
    qr_in_S = sum(1 for s in S if legendre(s, p) == 1)
    nqr_in_S = m - qr_in_S

    # Legendre product
    leg_prod = 1
    for s in S:
        leg_prod *= legendre(s, p)

    # Additive energy
    E = sum(1 for a in S for b in S for c in S for d in S if (a+b)%p == (c+d)%p)

    # Even/odd Ω balance
    even_omega = sum(omega[mm] for mm in range(0, len(omega), 2))
    odd_omega = sum(omega[mm] for mm in range(1, len(omega), 2))

    # Sumset size
    SS = set((a+b)%p for a in S for b in S)

    results.append({
        'S': sorted(S),
        'chi_per': chi_per,
        'mode_m': mode_m,
        'centroid': centroid,
        'skew': skew_m,
        'total_omega': total_omega,
        'even_omega': even_omega,
        'odd_omega': odd_omega,
        'even_minus_odd': even_omega - odd_omega,
        'qr_count': qr_in_S,
        'leg_prod': leg_prod,
        'E': E,
        'SS_size': len(SS),
        'omega': omega,
    })

results.sort(key=lambda r: r['chi_per'])

print(f"\n  Sorted by chi_per:")
print(f"  {'S':25s} {'chi':>5s} {'mode':>5s} {'centr':>6s} {'skew':>7s} {'Σodd-Σevn':>10s} {'QR':>3s} {'LP':>3s} {'E':>4s} {'|S+S|':>5s}")
for r in results:
    print(f"  {str(r['S']):25s} {r['chi_per']:5d} m={r['mode_m']:1d} {r['centroid']:6.2f} {r['skew']:+7.3f} {-r['even_minus_odd']:+10d} {r['qr_count']:3d} {r['leg_prod']:+3d} {r['E']:4d} {r['SS_size']:5d}")

# Test correlations
chis = [r['chi_per'] for r in results]
centroids = [r['centroid'] for r in results]
skews = [r['skew'] for r in results]
EO_diffs = [r['even_minus_odd'] for r in results]
qr_counts = [r['qr_count'] for r in results]
Es = [r['E'] for r in results]

print(f"\n  Correlations with chi_per:")
for name, vals in [('centroid', centroids), ('skewness', skews),
                    ('Σeven-Σodd', EO_diffs), ('QR count', qr_counts), ('E', Es)]:
    corr = np.corrcoef(chis, vals)[0,1]
    print(f"    corr(chi_per, {name}) = {corr:.4f}")

# KEY INSIGHT: chi_per = Σ(-1)^m Ω_m. The sign depends on whether
# even-indexed Ω's are larger or odd-indexed ones.
print(f"\n{'='*80}")
print(f"EVEN/ODD BALANCE IN Ω")
print(f"{'='*80}")

for r in results:
    print(f"\n  S = {r['S']}, chi_per = {r['chi_per']}")
    for mm in range(p):
        w = r['omega'][mm]
        sign = '+' if mm % 2 == 0 else '-'
        print(f"    m={mm:2d}: Ω = {w:6d}  ({sign}), cumulative chi = {sum((-1)**i * r['omega'][i] for i in range(mm+1)):+6d}")

# The Ω profile shape determines chi_per. Let's look at the
# "effective center of mass" of the alternating-sign weighted Ω:
print(f"\n{'='*80}")
print(f"ALTERNATING WEIGHTED ANALYSIS")
print(f"{'='*80}")

for r in results:
    omega = r['omega']
    # Partial chi_per from each half
    first_half = sum((-1)**mm * omega[mm] for mm in range(m+1))
    second_half = sum((-1)**mm * omega[mm] for mm in range(m+1, p))
    print(f"  S={str(r['S']):25s} chi_per={r['chi_per']:+3d}, first_half(0..{m})={first_half:+6d}, second_half({m+1}..{p-1})={second_half:+6d}")

print("\nDONE.")
