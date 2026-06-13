#!/usr/bin/env python3
"""
omega_from_Q.py — opus-2026-03-12-S69

Can we express Ω_m (or chi_per) as an explicit function of Q_1,...,Q_m?

Approach: The Q_k are the Fourier magnitudes of the tournament.
The Ω_m count valid paths in the tournament.
For circulant tournaments, there should be a relationship.

Strategy:
1. Compute Ω_m and Q_k for ALL orientations at p=7 (8 tournaments)
2. Try to fit Ω_m = f(e_1,...,e_m) where e_j = e_j(Q_1,...,Q_m)
3. If linear in e_j, find the coefficients
4. Test at p=11 and p=13

Note: e_j(Q) are INTEGER for all orientations (HYP-562).
And Ω_m is also an integer. So the relationship, if it exists, should be rational.
"""

import sys, time
import numpy as np
from itertools import product as cartprod
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

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        yield S

# ========== p = 7 ==========
p = 7
m = 3
print(f"{'='*70}")
print(f"p = {p}, m = {m}: Fitting Ω_m = f(e_j)")
print(f"{'='*70}")

data = []
for S in gen_orientations(p):
    h = CirculantHomology(n=p, S=S)
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)
    Q = compute_Q(p, S)
    esf = compute_esf(Q)
    data.append({'S': sorted(S), 'omega': omega, 'esf': esf, 'Q': Q})

print(f"\n  {'S':15s} | {'Ω':35s} | {'e_j':25s}")
for d in data:
    print(f"  {str(d['S']):15s} | {str(d['omega']):35s} | {str(d['esf']):25s}")

# At p=7: e_0=1, e_1=6 always (sum Q = m(p-m)/2 = 6). Only e_2 and e_3 vary.
# So the Q-dependent information is in (e_2, e_3) = (e_2, prod Q).
# We have 8 orientations but 2 Z_p*-orbits, so effectively 2 data points.

# Let's look at degree-by-degree fit: Ω_m = a + b*e_2 + c*e_3
print(f"\n  Degree-by-degree Q-dependence:")
for mm in range(p):
    omegas = [d['omega'][mm] for d in data]
    e2s = [d['esf'][2] for d in data]
    e3s = [d['esf'][3] for d in data]

    # Check if Ω_m depends on (e_2, e_3) only (i.e., same for same esf)
    by_esf = {}
    for d in data:
        key = (d['esf'][2], d['esf'][3])
        if key not in by_esf:
            by_esf[key] = []
        by_esf[key].append(d['omega'][mm])

    all_same = all(len(set(v)) == 1 for v in by_esf.values())
    print(f"  m={mm}: Ω values = {sorted(set(omegas))}, by (e2,e3): {dict((k, v[0]) for k,v in by_esf.items())}, Q-determined: {'YES' if all_same else 'NO'}")

# ========== p = 11 ==========
p = 11
m = 5
print(f"\n{'='*70}")
print(f"p = {p}, m = {m}: Fitting Ω_m = f(e_j)")
print(f"{'='*70}")

data = []
for S in gen_orientations(p):
    h = CirculantHomology(n=p, S=S)
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)
    Q = compute_Q(p, S)
    esf = compute_esf(Q)
    key = tuple(esf[2:])  # Only the non-universal part
    data.append({'S': sorted(S), 'omega': omega, 'esf': esf, 'esf_key': key})

# Group by e_j pattern
by_esf = {}
for d in data:
    key = d['esf_key']
    if key not in by_esf:
        by_esf[key] = []
    by_esf[key].append(d)

print(f"\n  Number of distinct e_j patterns: {len(by_esf)}")
for key, group in sorted(by_esf.items()):
    omega_patterns = set(tuple(d['omega']) for d in group)
    print(f"  e_j = {list(key)}: {len(group)} orientations, {len(omega_patterns)} distinct Ω profiles")
    for omega_pat in omega_patterns:
        chi = sum((-1)**mm * w for mm, w in enumerate(omega_pat))
        print(f"    Ω = {list(omega_pat)}, chi_per = {chi}")

# Key test: does e_j UNIQUELY determine Ω?
print(f"\n  DOES e_j DETERMINE Ω?")
for key, group in sorted(by_esf.items()):
    omega_patterns = set(tuple(d['omega']) for d in group)
    if len(omega_patterns) > 1:
        print(f"  NO! e_j = {list(key)} gives {len(omega_patterns)} different Ω profiles!")
    else:
        print(f"  YES at e_j = {list(key)}")

# Even if e_j doesn't determine Ω exactly, does it determine chi_per?
print(f"\n  DOES e_j DETERMINE chi_per?")
for key, group in sorted(by_esf.items()):
    chis = set(sum((-1)**mm * w for mm, w in enumerate(d['omega'])) for d in group)
    if len(chis) > 1:
        print(f"  NO! e_j = {list(key)} gives chi_per in {chis}!")
    else:
        print(f"  YES at e_j = {list(key)}: chi_per = {chis.pop()}")

# ========== Direct Q → chi_per ==========
print(f"\n{'='*70}")
print(f"DIRECT TEST: Does SORTED Q determine chi_per?")
print(f"{'='*70}")

for p in [7, 11]:
    m = (p-1)//2
    print(f"\n  p = {p}:")
    data = []
    for S in gen_orientations(p):
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        Q = compute_Q(p, S)
        chi = sum((-1)**mm * w for mm, w in enumerate(omega))
        Q_key = tuple(round(q, 4) for q in sorted(Q, reverse=True))
        data.append({'Q_key': Q_key, 'chi': chi})

    by_Q = {}
    for d in data:
        key = d['Q_key']
        if key not in by_Q:
            by_Q[key] = set()
        by_Q[key].add(d['chi'])

    for key, chis in sorted(by_Q.items()):
        status = "UNIQUE" if len(chis) == 1 else f"AMBIGUOUS: {chis}"
        print(f"    Q={list(key)}: chi_per = {status}")

print("\nDONE.")
