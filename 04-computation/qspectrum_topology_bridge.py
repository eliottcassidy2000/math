#!/usr/bin/env python3
"""
qspectrum_topology_bridge.py — opus-2026-03-12-S69

Bridge between Q-spectrum (Fourier magnitudes) and topological invariants.

Goal: Find what Q-spectrum invariant predicts chi_per (per-eigenspace Euler char).

Known data at p=11:
  Interval: chi_per=0,  Q=[12.34, 0.27, 1.45, 0.35, 0.58]
  Paley:    chi_per=1,  Q=[3, 3, 3, 3, 3]
  Orbit A:  chi_per=2,  Q=[8.74, 1.12, 0.42, 2.49, 2.22]
  Orbit B:  chi_per=0,  Q=[5.99, 4.12, 1.03, 3.61, 0.25]

All have sum Q = 15, but different distributions.
"""

import numpy as np
from math import comb
import mpmath
mpmath.mp.dps = 50

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = set()
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.add(a)
            else:
                S.add(b)
        yield frozenset(S)

def compute_Q(p, S):
    """Compute Q_k = |S_hat(k)|^2 for k=1,...,m."""
    m = (p-1)//2
    Qs = []
    for k in range(1, m+1):
        val = sum(np.exp(2j*np.pi*s*k/p) for s in S)
        Qs.append(abs(val)**2)
    return np.array(Qs)

def compute_esf(Qs):
    """Compute elementary symmetric functions of Q values."""
    m = len(Qs)
    # Use numpy polynomial: prod(t - Q_k) = sum e_j * (-1)^j * t^{m-j}
    coeffs = np.polynomial.polynomial.polyfromroots(Qs)
    # coeffs[j] = e_{m-j} * (-1)^{m-j} ... actually:
    # polyfromroots gives coefficients of t^0, t^1, ..., t^m
    # prod(t - Q_k) = sum_{j=0}^m (-1)^{m-j} e_{m-j} t^j  where e_0 = 1
    # So coeffs[j] = (-1)^{m-j} e_{m-j}
    esf = []
    for j in range(m+1):
        e_j = coeffs[m-j] * (-1)**j if j <= m else 0
        esf.append(float(np.real(e_j)))
    return esf

def zpstar_orbit(p, S):
    """Get the Z_p* orbit of S (all S' = a*S mod p for a in Z_p*)."""
    orbit = set()
    for a in range(1, p):
        S_a = frozenset((a * s) % p for s in S)
        orbit.add(S_a)
    return orbit

def difference_set(S, p):
    """Compute S - S mod p (excluding 0)."""
    diffs = set()
    for a in S:
        for b in S:
            d = (a - b) % p
            if d != 0:
                diffs.add(d)
    return diffs

p = 11
m = (p - 1) // 2
S_paley = frozenset(d for d in range(1, p) if legendre(d, p) == 1)
S_interval = frozenset(range(1, m + 1))

# Known orbit representatives and their Betti/chi data
orbit_reps = {
    'Interval': frozenset({1,2,3,4,5}),
    'Paley': frozenset({1,3,4,5,9}),
    'Orbit_A': frozenset({1,2,3,4,6}),
    'Orbit_B': frozenset({1,6,7,8,9}),
}

# Known topological data
topo_data = {
    'Interval': {'beta': [1,1,0,0,0,0,0,0,0,0,0], 'chi': 0, 'chi_per': 0, 'orbit_size': 10},
    'Paley': {'beta': [1,0,0,0,0,5,15,0,0,0,0], 'chi': 11, 'chi_per': 1, 'orbit_size': 2},
    'Orbit_A': {'beta': [1,0,0,0,0,0,21,0,0,0,0], 'chi': 22, 'chi_per': 2, 'orbit_size': 10},
    'Orbit_B': {'beta': [1,0,0,0,0,14,13,0,0,0,0], 'chi': 0, 'chi_per': 0, 'orbit_size': 10},
}

print("=" * 80)
print("Q-SPECTRUM → TOPOLOGY BRIDGE (p=11)")
print("=" * 80)

# Compute Q-spectrum invariants for each orbit
for name, S in orbit_reps.items():
    Q = compute_Q(p, S)
    esf = compute_esf(Q)
    td = topo_data[name]
    ds = difference_set(S, p)

    # Sort Q for canonical form
    Q_sorted = np.sort(Q)[::-1]

    print(f"\n{'─'*60}")
    print(f"  {name}: S = {sorted(S)}")
    print(f"  chi_per = {td['chi_per']}, chi = {td['chi']}")
    print(f"  beta = {td['beta']}")
    print(f"  S-S = {sorted(ds)} ({'FULL' if len(ds) == p-1 else f'missing {p-1-len(ds)}'})")
    print(f"  Q = [{', '.join(f'{q:.4f}' for q in Q)}]")
    print(f"  Q sorted = [{', '.join(f'{q:.4f}' for q in Q_sorted)}]")

    # Basic invariants
    print(f"\n  Basic invariants:")
    print(f"    sum Q = {sum(Q):.6f} (should be {m*(p-m)/2})")
    print(f"    prod Q = {np.prod(Q):.6f}")
    print(f"    sum Q^2 = {sum(Q**2):.6f}")
    print(f"    sum Q^3 = {sum(Q**3):.6f}")
    print(f"    sum Q^4 = {sum(Q**4):.6f}")
    print(f"    sum 1/Q = {sum(1/Q):.6f} (should be {p-2} for Interval)")

    # ESF
    print(f"\n  Elementary symmetric functions e_j(Q):")
    mv_coeffs = [comb(m+j, 2*j) for j in range(m+1)]
    for j in range(m+1):
        mv = mv_coeffs[j]
        print(f"    e_{j} = {esf[j]:12.4f}  (Morgan-Voyce C({m+j},{2*j}) = {mv})")

    # Power sums (Newton's identities)
    print(f"\n  Power sums p_k = sum Q_i^k:")
    for k in range(1, 6):
        pk = sum(Q**k)
        print(f"    p_{k} = {pk:.6f}")

    # Special invariants
    prod_1pQ = np.prod(1 + Q)
    prod_Qm1 = np.prod(Q - 1)
    print(f"\n  Special products:")
    print(f"    prod(1+Q) = {prod_1pQ:.4f}")
    print(f"    prod(Q-1) = {prod_Qm1:.4f}")
    print(f"    prod(Q)/prod(Q-1) = {np.prod(Q)/prod_Qm1:.4f}")

    # IPR-like
    ipr = sum(Q**2) / sum(Q)**2
    print(f"    IPR = sum(Q^2)/sum(Q)^2 = {ipr:.6f}")

    # Entropy
    Q_norm = Q / sum(Q)
    entropy = -sum(q * np.log(q) if q > 0 else 0 for q in Q_norm)
    print(f"    Entropy H(Q/sum) = {entropy:.6f} (max = {np.log(m):.6f})")

    # Discriminant
    disc = 1.0
    for i in range(m):
        for j in range(i+1, m):
            disc *= (Q[i] - Q[j])**2
    print(f"    disc(Q) = {disc:.4f}")
    if disc > 0:
        print(f"    disc / p^(m-1) = {disc / p**(m-1):.6f}")

# Now look for which invariant separates chi_per
print(f"\n{'='*80}")
print("SEPARATION ANALYSIS: Which invariant predicts chi_per?")
print("=" * 80)

# Collect invariants
inv_data = {}
for name, S in orbit_reps.items():
    Q = compute_Q(p, S)
    td = topo_data[name]
    inv_data[name] = {
        'chi_per': td['chi_per'],
        'sum_Q2': sum(Q**2),
        'sum_Q3': sum(Q**3),
        'sum_Q4': sum(Q**4),
        'prod_Q': np.prod(Q),
        'prod_1pQ': np.prod(1+Q),
        'prod_Qm1': np.prod(Q-1),
        'ipr': sum(Q**2)/sum(Q)**2,
        'entropy': -sum(q*np.log(q) if q > 0 else 0 for q in Q/sum(Q)),
        'esf2': compute_esf(Q)[2],
        'esf3': compute_esf(Q)[3],
        'esf4': compute_esf(Q)[4],
        'sum_1oQ': sum(1/Q),
        'max_Q': max(Q),
        'min_Q': min(Q),
        'range_Q': max(Q) - min(Q),
    }

print(f"\n{'Name':12s} chi_per  sum_Q2    sum_Q3    sum_Q4    prod_Q    prod(1+Q)  prod(Q-1)")
for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    d = inv_data[name]
    print(f"  {name:10s} {d['chi_per']:3d}    {d['sum_Q2']:8.2f}  {d['sum_Q3']:8.2f}  {d['sum_Q4']:10.2f}  {d['prod_Q']:8.4f}  {d['prod_1pQ']:10.2f}  {d['prod_Qm1']:+8.4f}")

print(f"\n{'Name':12s} chi_per  e_2       e_3       e_4       sum_1/Q   entropy   ipr")
for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    d = inv_data[name]
    print(f"  {name:10s} {d['chi_per']:3d}    {d['esf2']:8.2f}  {d['esf3']:8.2f}  {d['esf4']:8.2f}  {d['sum_1oQ']:8.4f}  {d['entropy']:8.6f}  {d['ipr']:.6f}")

# Test specific combinations
print(f"\n{'='*80}")
print("CANDIDATE PREDICTORS FOR chi_per")
print("=" * 80)

for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    Q = compute_Q(p, orbit_reps[name])
    chi_per = topo_data[name]['chi_per']

    # Try: round(sum_Q3 / something)?
    # Try: integer parts of esf?
    # Try: prod_Q relative to Interval

    # Legendre product
    S = orbit_reps[name]
    leg_prod = 1
    for s in S:
        leg_prod *= legendre(s, p)

    # Quadratic residue content
    qr_count = sum(1 for s in S if legendre(s, p) == 1)
    nqr_count = sum(1 for s in S if legendre(s, p) == -1)

    # Additive energy
    E = 0
    S_list = sorted(S)
    for a in S_list:
        for b in S_list:
            for c in S_list:
                for d in S_list:
                    if (a + b) % p == (c + d) % p:
                        E += 1

    print(f"  {name:10s}: chi_per={chi_per}, leg_prod={leg_prod:+d}, QR/NQR={qr_count}/{nqr_count}, E={E}")

    # Number of triples (a,b,c) in S with a+b≡c mod p
    triple_count = 0
    for a in S:
        for b in S:
            if (a + b) % p in S:
                triple_count += 1
    print(f"              additive_triples={triple_count}, sum_Q2={sum(Q**2):.4f}")

# Now test at p=7 as well
print(f"\n{'='*80}")
print("Q-SPECTRUM → TOPOLOGY BRIDGE (p=7)")
print("=" * 80)

p7 = 7
m7 = 3
S_paley7 = frozenset(d for d in range(1, p7) if legendre(d, p7) == 1)
S_interval7 = frozenset(range(1, m7 + 1))

p7_orbits = {
    'Interval': frozenset({1,2,3}),
    'Paley': frozenset({1,2,4}),
}

p7_topo = {
    'Interval': {'chi_per': 0, 'beta': [1,1,0,0,0,0,0]},
    'Paley': {'chi_per': 1, 'beta': [1,0,0,0,6,0,0]},
}

for name, S in p7_orbits.items():
    Q = compute_Q(p7, S)
    td = p7_topo[name]
    ds = difference_set(S, p7)

    S_list = sorted(S)
    E = sum(1 for a in S for b in S for c in S for d in S if (a+b)%p7 == (c+d)%p7)
    triple_count = sum(1 for a in S for b in S if (a+b)%p7 in S)

    qr_count = sum(1 for s in S if legendre(s, p7) == 1)
    leg_prod = 1
    for s in S:
        leg_prod *= legendre(s, p7)

    print(f"\n  {name}: S = {sorted(S)}, chi_per = {td['chi_per']}")
    print(f"    Q = [{', '.join(f'{q:.4f}' for q in Q)}]")
    print(f"    sum Q^2 = {sum(Q**2):.4f}, prod Q = {np.prod(Q):.4f}")
    print(f"    S-S = {sorted(ds)} ({'FULL' if len(ds) == p7-1 else f'missing {p7-1-len(ds)}'})")
    print(f"    leg_prod = {leg_prod}, QR/NQR = {qr_count}/{m7-qr_count}")
    print(f"    additive_energy = {E}, triples = {triple_count}")

# Deeper: look at eigenspace ANOMALY PATTERN vs Q
print(f"\n{'='*80}")
print("ANOMALY PATTERN vs Q-SPECTRUM STRUCTURE (p=11)")
print("=" * 80)

# From THM-144:
anomaly_patterns = {
    'Interval': [(-1, 0, 0, 0, 0, 0)],
    'Paley':    [(-1, 1, -1, 1, -1, -4)],
    'Orbit_A':  [(-1, 1, -1, 1, -1, 1)],
    'Orbit_B':  [(-1, 1, -1, 1, -1, -2)],
}

for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    Q = compute_Q(p, orbit_reps[name])
    chi_per = topo_data[name]['chi_per']
    pattern = anomaly_patterns[name][0]

    # The anomaly alternation stops at m=6 for deep, at m=1 for shallow
    # The LAST anomaly value differs: +1, -4, +1, -2
    # This might encode chi_per!

    last_nonzero = 0
    last_val = 0
    for i, v in enumerate(pattern):
        if v != 0:
            last_nonzero = i + 1  # 1-indexed degree
            last_val = v

    anomaly_sum = sum(pattern)

    print(f"\n  {name}: chi_per = {chi_per}")
    print(f"    Anomaly pattern: Δ = {pattern}")
    print(f"    Anomaly depth = {last_nonzero}, last Δ = {last_val}")
    print(f"    Sum of Δ = {anomaly_sum}")
    print(f"    |Sum of Δ| = {abs(anomaly_sum)}")

# KEY TEST: does |sum_Delta| = chi_per?
print(f"\n  CORRELATION TEST:")
print(f"  Interval: |sum Δ| = {abs(sum(anomaly_patterns['Interval'][0]))}, chi_per = {topo_data['Interval']['chi_per']}")
print(f"  Paley:    |sum Δ| = {abs(sum(anomaly_patterns['Paley'][0]))},  chi_per = {topo_data['Paley']['chi_per']}")
print(f"  Orbit_A:  |sum Δ| = {abs(sum(anomaly_patterns['Orbit_A'][0]))}, chi_per = {topo_data['Orbit_A']['chi_per']}")
print(f"  Orbit_B:  |sum Δ| = {abs(sum(anomaly_patterns['Orbit_B'][0]))}, chi_per = {topo_data['Orbit_B']['chi_per']}")

# TEST 2: sum of ALTERNATING anomaly
print(f"\n  ALTERNATING SUM TEST (Σ (-1)^m Δ_m):")
for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    pattern = anomaly_patterns[name][0]
    alt_sum = sum((-1)**i * v for i, v in enumerate(pattern))
    print(f"    {name}: Σ(-1)^m Δ_m = {alt_sum}, chi_per = {topo_data[name]['chi_per']}")

# TEST 3: weighted anomaly sum
print(f"\n  DEEPER ANOMALY ANALYSIS:")
for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    pattern = anomaly_patterns[name][0]
    chi_per = topo_data[name]['chi_per']

    # chi_per should equal sum over m of (-1)^m * (beta_m^(0) - beta_m^(k≥1))
    # since chi = p * chi_per = sum (-1)^m beta_m
    # and chi = chi^(0) + (p-1) * chi^(k≥1)

    # From the anomaly: Δ_m = r_m^(0) - r_m^(k≥1)
    # beta_m^(k) = Omega_m^(k) - r_m^(k) - r_{m+1}^(k)
    # Since Omega_m^(k) is the SAME for all k (THM-125):
    # beta_m^(0) - beta_m^(k≥1) = -Δ_m + Δ_{m+1}

    print(f"  {name}:")
    for i in range(len(pattern)):
        delta_m = pattern[i]
        delta_m1 = pattern[i+1] if i+1 < len(pattern) else 0
        betti_diff = -delta_m + delta_m1
        print(f"    m={i+1}: Δ_m = {delta_m:+d}, Δ_{i+2} = {delta_m1:+d}, β^(0)-β^(k≥1) = {betti_diff:+d}")

print(f"\n  FORMULA: chi_per = (1/(p-1)) * sum_m (-1)^m beta_m")
print(f"    Since chi = sum (-1)^m beta_m = p * chi_per (from THM-125 + decomposition)")
print(f"    chi = chi^(0) + (p-1) * chi_per_k>=1")
print(f"    chi_per = chi^(0) when all k>=1 contribute equally (which they do!)")
print(f"    So chi_per IS chi^(k=1) = sum_m (-1)^m beta_m^(k=1)")

print("\nDONE.")
