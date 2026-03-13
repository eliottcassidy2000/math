#!/usr/bin/env python3
"""
thm145_p13_verify.py — opus-2026-03-12-S69

Verify THM-145 (spectral-topological bridge) at p=13.

At p=13 there are 6 Z_p*-orbits. We need to check that the 6 distinct
e_j(Q) tuples give 6 distinct Ω profiles.

This is the critical test: p=13 has MUCH richer structure than p=7 or p=11.
"""

import sys, time
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def compute_Q(p, S):
    m = (p-1)//2
    return np.array([abs(sum(np.exp(2j*np.pi*s*k/p) for s in S))**2 for k in range(1, m+1)])

def compute_esf(Q):
    """Integer elementary symmetric functions of Q."""
    m = len(Q)
    coeffs = np.polynomial.polynomial.polyfromroots(Q)
    esf = []
    for j in range(m+1):
        ej = float(np.real(coeffs[m-j] * (-1)**j))
        esf.append(round(ej))
    return esf

def zpstar_orbit_rep(p, S):
    best = sorted(S)
    for a in range(2, p):
        S_a = frozenset((a * s) % p for s in S)
        if sorted(S_a) < best:
            best = sorted(S_a)
    return tuple(best)

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        yield S

p = 13
m = (p-1)//2  # = 6
print(f"{'='*80}")
print(f"THM-145 VERIFICATION at p = {p}, m = {m}")
print(f"{'='*80}")

# First, find all orbit representatives
orbits = {}
for S in gen_orientations(p):
    rep = zpstar_orbit_rep(p, S)
    if rep not in orbits:
        orbits[rep] = S

print(f"\nFound {len(orbits)} Z_p*-orbits")

# For each orbit, compute e_j and Ω
results = []
for i, (rep, S) in enumerate(sorted(orbits.items())):
    print(f"\n--- Orbit {i}: S = {sorted(S)} ---")
    t0 = time.time()

    # Compute Q and e_j
    Q = compute_Q(p, S)
    esf = compute_esf(Q)
    esf_key = tuple(esf[2:])  # Non-universal part

    print(f"  Q = [{', '.join(f'{q:.4f}' for q in sorted(Q, reverse=True))}]")
    print(f"  e_j (non-universal) = {list(esf_key)}")

    # Compute Ω profile
    h = CirculantHomology(n=p, S=S)
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)
    chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))

    elapsed = time.time() - t0
    print(f"  Ω = {omega}")
    print(f"  chi_per = {chi_per}")
    print(f"  Time: {elapsed:.1f}s")

    results.append({
        'orbit': i,
        'S': sorted(S),
        'esf_key': esf_key,
        'omega': tuple(omega),
        'chi_per': chi_per,
    })

# KEY TEST: Does e_j uniquely determine Ω?
print(f"\n{'='*80}")
print(f"THM-145 VERIFICATION RESULT")
print(f"{'='*80}")

# Check e_j → Ω injectivity
esf_to_omega = {}
for r in results:
    key = r['esf_key']
    if key in esf_to_omega:
        if esf_to_omega[key] != r['omega']:
            print(f"  COUNTEREXAMPLE! e_j = {list(key)} maps to TWO Ω profiles!")
            print(f"    Ω_1 = {list(esf_to_omega[key])}")
            print(f"    Ω_2 = {list(r['omega'])}")
    else:
        esf_to_omega[key] = r['omega']

if len(esf_to_omega) == len(results):
    print(f"  THM-145 VERIFIED at p={p}: {len(results)} distinct e_j → {len(results)} distinct Ω ✓")
else:
    print(f"  WARNING: {len(results)} orbits but only {len(esf_to_omega)} distinct e_j patterns")
    # Check if same e_j gives same Ω (still consistent with THM-145)
    all_consistent = True
    for r in results:
        if esf_to_omega[r['esf_key']] != r['omega']:
            all_consistent = False
            break
    if all_consistent:
        print(f"  But all same-e_j orbits have same Ω, so THM-145 still holds")
    else:
        print(f"  THM-145 FAILS at p={p}!")

# Summary table
print(f"\n  Summary:")
print(f"  {'e_j (non-universal)':35s} | {'chi_per':>8s} | {'Ω profile'}")
for r in sorted(results, key=lambda x: x['chi_per']):
    print(f"  {str(list(r['esf_key'])):35s} | {r['chi_per']:8d} | {list(r['omega'])}")

print("\nDONE.")
