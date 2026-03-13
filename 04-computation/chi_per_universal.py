#!/usr/bin/env python3
"""
chi_per_universal.py — opus-2026-03-12-S69

Compute chi_per = Σ(-1)^m Ω_m / p for all Z_p* orbits at p=7, p=11, p=13.

Key question: can chi_per be predicted from Q-spectrum invariants?

Also investigate: what COMBINATORIAL property of S determines chi_per?
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
    """Get canonical representative of Z_p* orbit of S."""
    best = sorted(S)
    for a in range(2, p):
        S_a = frozenset((a * s) % p for s in S)
        if sorted(S_a) < best:
            best = sorted(S_a)
    return frozenset(best)

def compute_Q(p, S):
    m = (p-1)//2
    Qs = []
    for k in range(1, m+1):
        val = sum(np.exp(2j*np.pi*s*k/p) for s in S)
        Qs.append(abs(val)**2)
    return np.array(Qs)

def difference_set(S, p):
    return set((a-b) % p for a in S for b in S if a != b)

def additive_energy(S, p):
    S_list = sorted(S)
    return sum(1 for a in S_list for b in S_list for c in S_list for d in S_list if (a+b)%p == (c+d)%p)

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n{'='*80}")
    print(f"p = {p}, m = {m}")
    print(f"{'='*80}")

    # Group orientations by Z_p* orbit
    orbits = {}
    for S in gen_orientations(p):
        rep = zpstar_orbit_rep(p, S)
        if rep not in orbits:
            orbits[rep] = []
        orbits[rep].append(S)

    print(f"Number of Z_p* orbits: {len(orbits)}")

    results = []
    for rep, members in sorted(orbits.items(), key=lambda x: sorted(x[0])):
        S = members[0]  # Use first member for computation
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)

        chi = sum((-1)**mm * w for mm, w in enumerate(omega))
        chi_per = chi  # This IS chi_per since omega gives per-eigenspace dims

        Q = compute_Q(p, S)
        Q_sorted = sorted(Q)[::-1]
        ds = difference_set(S, p)
        full_ds = len(ds) == p-1

        # Q-spectrum invariants
        prod_Q = np.prod(Q)
        sum_Q2 = sum(Q**2)
        e2 = sum(Q[i]*Q[j] for i in range(m) for j in range(i+1, m))

        # Additive energy
        E = additive_energy(S, p)

        # Label
        S_paley = frozenset(d for d in range(1, p) if legendre(d, p) == 1)
        S_interval = frozenset(range(1, m+1))
        if S == S_paley:
            label = "Paley"
        elif S == S_interval:
            label = "Interval"
        else:
            label = f"Orbit-{len(results)-1}"

        r = {
            'S': sorted(rep),
            'label': label,
            'orbit_size': len(members),
            'omega': omega,
            'chi': chi,
            'chi_per': chi_per,
            'Q_sorted': Q_sorted,
            'prod_Q': float(prod_Q),
            'sum_Q2': float(sum_Q2),
            'e2': float(e2),
            'full_ds': full_ds,
            'additive_energy': E,
        }
        results.append(r)

    # Print results sorted by chi_per
    results.sort(key=lambda r: r['chi_per'])

    print(f"\n  {'Label':15s} {'S':25s} {'|orb|':>5s} {'chi_per':>7s} {'S-S=Z*':>6s} {'prod_Q':>8s} {'Σ Q²':>8s} {'E':>4s}")
    for r in results:
        print(f"  {r['label']:15s} {str(r['S']):25s} {r['orbit_size']:5d} {r['chi_per']:7d} {'YES' if r['full_ds'] else 'NO':>6s} {r['prod_Q']:8.1f} {r['sum_Q2']:8.1f} {r['additive_energy']:4d}")

    print(f"\n  Ω profiles:")
    for r in results:
        print(f"  {r['label']:15s} chi_per={r['chi_per']:2d}  Ω = {r['omega']}")

    # Key: look at what's DIFFERENT in Ω between chi_per = 0, 1, 2
    print(f"\n  Ω differences (relative to orbit with min chi_per):")
    min_chi_r = results[0]  # min chi_per
    for r in results[1:]:
        diffs = [r['omega'][mm] - min_chi_r['omega'][mm] for mm in range(len(r['omega']))]
        nonzero = [(mm, d) for mm, d in enumerate(diffs) if d != 0]
        print(f"    {r['label']} - {min_chi_r['label']}: {nonzero}")
        print(f"      Sum of diffs: {sum(diffs)}")
        print(f"      Alt sum of diffs: {sum((-1)**mm * d for mm, d in enumerate(diffs))} (= chi_per diff = {r['chi_per'] - min_chi_r['chi_per']})")

    # Q-sorted values
    print(f"\n  Q spectra (sorted, descending):")
    for r in results:
        print(f"    {r['label']:15s} chi={r['chi_per']:2d}  Q = [{', '.join(f'{q:.4f}' for q in r['Q_sorted'])}]")

    # IMPORTANT: check if Q-spectrum uniquely determines the orbit
    q_sets = {}
    for r in results:
        q_key = tuple(round(q, 6) for q in r['Q_sorted'])
        if q_key in q_sets:
            print(f"\n  WARNING: orbits {r['label']} and {q_sets[q_key]} have same Q!")
        q_sets[q_key] = r['label']

print(f"\n{'='*80}")
print("UNIVERSAL PATTERN SEARCH")
print("=" * 80)
print("Looking for a function f(Q_1,...,Q_m) = chi_per...")
print("At each p, chi_per values and Q-invariants listed above.")
print("If chi_per = f(prod_Q, sum_Q², e_2, ...), the formula must work at ALL p.")

print("\nDONE.")
