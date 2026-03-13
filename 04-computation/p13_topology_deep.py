#!/usr/bin/env python3
"""
p13_topology_deep.py — opus-2026-03-12-S69

Deep analysis of the p=13 topological landscape.

Key mystery: Orbit-2 ({1,2,3,5,6,9}) has chi_per = -17, meaning
total chi = -221 = Σ(-1)^m β_m.

This requires MASSIVE odd-degree Betti numbers. Where does the homology live?

Also: Orbit-4 and Interval both have chi_per = 1 despite different S-S structure.
What do they share?
"""

import sys
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

p = 13
m = (p-1)//2  # = 6

# Define the 6 orbit representatives
orbits = {
    'Interval':   frozenset({1,2,3,4,5,6}),
    'Orbit-0':    frozenset({1,2,3,4,5,7}),  # chi_per = -3
    'Orbit-1':    frozenset({1,2,3,4,6,8}),  # chi_per = 4
    'Orbit-2':    frozenset({1,2,3,5,6,9}),  # chi_per = -17
    'Orbit-3':    frozenset({1,2,3,5,7,9}),  # chi_per = -8
    'Orbit-4':    frozenset({1,2,6,8,9,10}), # chi_per = 1
}

chi_per_known = {
    'Interval': 1, 'Orbit-0': -3, 'Orbit-1': 4,
    'Orbit-2': -17, 'Orbit-3': -8, 'Orbit-4': 1,
}

def compute_Q(p, S):
    m = (p-1)//2
    Qs = []
    for k in range(1, m+1):
        val = sum(np.exp(2j*np.pi*s*k/p) for s in S)
        Qs.append(abs(val)**2)
    return np.array(Qs)

def difference_set(S, p):
    return set((a-b) % p for a in S for b in S if a != b)

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

print(f"{'='*80}")
print(f"p = {p} TOPOLOGICAL LANDSCAPE — DEEP ANALYSIS")
print(f"{'='*80}")

# For each orbit, compute:
# 1. Ω profile (from chi_per_universal.py output)
# 2. Q-spectrum analysis
# 3. Additive structure of S

for name in ['Interval', 'Orbit-0', 'Orbit-1', 'Orbit-2', 'Orbit-3', 'Orbit-4']:
    S = orbits[name]
    chi_per = chi_per_known[name]
    Q = compute_Q(p, S)
    ds = difference_set(S, p)

    print(f"\n{'─'*60}")
    print(f"  {name}: S = {sorted(S)}, chi_per = {chi_per}")
    print(f"  S-S = {sorted(ds)} ({'FULL' if len(ds)==p-1 else f'missing {p-1-len(ds)}'})")
    print(f"  Q = [{', '.join(f'{q:.3f}' for q in sorted(Q, reverse=True))}]")
    print(f"  prod Q = {np.prod(Q):.1f}, sum Q² = {sum(Q**2):.1f}")

    # QR structure
    qr = [s for s in S if legendre(s, p) == 1]
    nqr = [s for s in S if legendre(s, p) == -1]
    print(f"  QR content: {sorted(qr)} ({len(qr)} QR, {len(nqr)} NQR)")

    # Representation function r_S(d) = #{(a,b) : a-b ≡ d, a,b ∈ S}
    rep_func = {}
    for d in range(1, p):
        cnt = sum(1 for a in S for b in S if (a-b) % p == d)
        rep_func[d] = cnt
    rep_vals = sorted(set(rep_func.values()))
    rep_hist = {v: sum(1 for d in range(1,p) if rep_func[d]==v) for v in rep_vals}
    print(f"  Representation r_S(d) histogram: {rep_hist}")
    print(f"  r_S values: {rep_vals}, mean = {sum(rep_func.values())/(p-1):.3f} (should be {m*(m-1)/(p-1):.3f})")

    # Is S a (v,k,lambda)-difference set?
    # That requires r_S(d) = constant for all d
    is_perfect_ds = len(rep_vals) == 1
    print(f"  Perfect difference set: {'YES' if is_perfect_ds else 'NO'}")

    # Additive energy
    E = sum(1 for a in S for b in S for c in S for d in S if (a+b)%p == (c+d)%p)
    print(f"  Additive energy E = {E}")

    # Multiplicative structure: is S a union of cosets?
    # Check if S is closed under multiplication by some g
    generators = []
    for g in range(2, p):
        gS = frozenset((g*s) % p for s in S)
        if gS == S:
            generators.append(g)
    if generators:
        print(f"  S closed under multiplication by: {generators}")
    else:
        print(f"  S not closed under any non-trivial multiplication")

    # Sumset structure
    SS = set((a+b)%p for a in S for b in S)
    print(f"  |S+S| = {len(SS)} (of {p}), 0∈S+S: {'yes' if 0 in SS else 'no'}")

# KEY COMPARISON: What distinguishes chi_per = 1 cases (Interval vs Orbit-4)?
print(f"\n{'='*80}")
print(f"COMPARISON: chi_per = 1 (Interval vs Orbit-4)")
print(f"{'='*80}")

for name in ['Interval', 'Orbit-4']:
    S = orbits[name]
    Q = compute_Q(p, S)
    ds = difference_set(S, p)

    print(f"\n  {name}: S = {sorted(S)}")
    print(f"  Q = [{', '.join(f'{q:.4f}' for q in sorted(Q, reverse=True))}]")
    print(f"  e_j: ", end="")
    esf = []
    for j in range(m+1):
        # Compute e_j using Newton's identities
        if j == 0:
            ej = 1.0
        else:
            # Use polynomial coefficients
            coeffs = np.polynomial.polynomial.polyfromroots(Q)
            # coeffs[k] = (-1)^{m-k} * e_{m-k}
            ej = float(np.real(coeffs[m-j] * (-1)**j))
        esf.append(round(ej))
    print(esf)

# KEY: What makes chi_per = -17?
print(f"\n{'='*80}")
print(f"WHAT MAKES chi_per = -17? (Orbit-2: S={{1,2,3,5,6,9}})")
print(f"{'='*80}")

S = orbits['Orbit-2']
Q = compute_Q(p, S)
print(f"  Q = [{', '.join(f'{q:.6f}' for q in Q)}]")
print(f"  TWO-VALUED Q: three copies of {Q[0]:.6f} and three copies of {Q[3]:.6f}")
print(f"  These equal (p ± √p)/4:")
q_hi = (p + np.sqrt(p))/4
q_lo = (p - np.sqrt(p))/4
print(f"  (p+√p)/4 = {q_hi:.6f}, (p-√p)/4 = {q_lo:.6f}")
print(f"  Actual Q values: {sorted(set(round(q,4) for q in Q))}")

# Check if S is a "half-Paley" — union of QR cosets
QR_13 = {d for d in range(1,p) if legendre(d,p)==1}
NQR_13 = {d for d in range(1,p) if legendre(d,p)==-1}
print(f"  QR_13 = {sorted(QR_13)}")
print(f"  NQR_13 = {sorted(NQR_13)}")
S_sorted = sorted(S)
print(f"  S = {S_sorted}")
print(f"  S ∩ QR = {sorted(S & QR_13)} (size {len(S & QR_13)})")
print(f"  S ∩ NQR = {sorted(S & NQR_13)} (size {len(S & NQR_13)})")

# Is S a (13,6,2)-difference set (like a Paley-type)?
# For (v,k,λ)-DS: r(d) = λ for all d ≠ 0
rep = {}
for d in range(1,p):
    rep[d] = sum(1 for a in S for b in S if (a-b)%p == d)
print(f"  r_S distribution: {sorted(set(rep.values()))} with counts {dict(sorted({v:sum(1 for d in rep if rep[d]==v) for v in set(rep.values())}.items()))}")

# Connection to Paley graph
print(f"\n  At p=13 ≡ 1 mod 4, Paley construction:")
print(f"  Paley graph (undirected): edges where a-b ∈ QR")
print(f"  Conference graph: 13 vertices, regular degree 6")
print(f"  This orbit has Q two-valued → CONFERENCE matrix eigenvalues!")

# What is the Ω profile for Orbit-2?
print(f"\n  Computing Ω for Orbit-2...")
h = CirculantHomology(n=p, S=S)
h._ensure_enumerated(p)
omega = h.omega_dims(max_degree=p-1)
print(f"  Ω = {omega}")
chi_computed = sum((-1)**mm * w for mm, w in enumerate(omega))
print(f"  chi = Σ(-1)^m Ω_m = {chi_computed} (expected {-17})")

# Compare with Interval Ω to understand where the difference comes from
S_int = orbits['Interval']
h_int = CirculantHomology(n=p, S=S_int)
h_int._ensure_enumerated(p)
omega_int = h_int.omega_dims(max_degree=p-1)
print(f"\n  Interval Ω = {omega_int}")
print(f"  Orbit-2 Ω  = {omega}")
print(f"\n  Ω differences (Orbit-2 - Interval):")
for mm in range(p):
    d = omega[mm] - omega_int[mm]
    if d != 0:
        print(f"    m={mm}: Δ = {d:+d}")

# The Ω profile tells us where β must live
print(f"\n  Per-eigenspace Ω for Orbit-2:")
omega_per = [w//p if w%p == 0 else f"{w}/{p}" for w in omega]
print(f"  Ω/p = {omega_per}")
print(f"  chi_per = {chi_computed}")

# From chi = -17 and β ≥ 0: sum of even-degree β - sum of odd-degree β = -17
# So odd-degree Betti DOMINATE by 17 (per eigenspace) or 221 (total)
print(f"\n  Since chi_per = {chi_computed}, odd-degree Betti sum exceeds")
print(f"  even-degree Betti sum by {abs(chi_computed)} per eigenspace")

print("\nDONE.")
