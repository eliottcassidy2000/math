#!/usr/bin/env python3
"""
p13_spectral_analysis.py -- kind-pasteur-2026-03-13-S60

Extend the spectral/overlap analysis to p=13.
At p=13: m=6, 2^6=64 orientations.
Can't compute H directly (too expensive), but CAN compute:
- Eigenvalue spectra (instant)
- sum|lambda|^{2k} (instant)
- c3, c5 cycle counts (fast)
- Possibly c7 (marginal)

Questions:
1. How many H-classes (eigenvalue types)?
2. Does the flat spectrum / H-maximization connection extend?
3. Does the sum_4 separation still work?
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict

p = 13
m = (p - 1) // 2
omega = cmath.exp(2j * cmath.pi / p)

S_qr_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
print(f"p={p}, m={m}, QR={sorted(S_qr_set)}")

def bits_to_S(bits):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))
    return S

def eigenvalues(S, p):
    return [sum(omega**(j*s) for s in S) for j in range(p)]

# Phase 1: Eigenvalue analysis for all 64 orientations (instant)
print(f"\n{'='*70}")
print(f"  EIGENVALUE CLASSIFICATION AT p={p} ({1 << m} orientations)")
print(f"{'='*70}")

spectral_types = defaultdict(list)
for bits in range(1 << m):
    S = bits_to_S(bits)
    eigs = eigenvalues(S, p)
    sum_4 = round(sum(abs(e)**4 for e in eigs).real, 4)
    sum_6 = round(sum(abs(e)**6 for e in eigs).real, 4)
    abs_multiset = tuple(sorted([round(abs(e), 6) for e in eigs], reverse=True))
    key = (sum_4, sum_6)
    spectral_types[key].append(bits)

print(f"  Distinct (sum_4, sum_6) pairs: {len(spectral_types)}")
print(f"  Class sizes: {sorted([len(v) for v in spectral_types.values()], reverse=True)}")

for (s4, s6), bits_list in sorted(spectral_types.items()):
    # Count #QR for each
    qr_counts = []
    for bits in bits_list:
        S = bits_to_S(bits)
        n_qr = sum(1 for s in S if s in S_qr_set)
        qr_counts.append(n_qr)

    # Eigenvalue spectrum for representative
    S = bits_to_S(bits_list[0])
    eigs = eigenvalues(S, p)
    abs_vals = sorted([abs(e) for e in eigs], reverse=True)

    det_val = 1
    for e in eigs:
        det_val *= e
    det_val = round(abs(det_val), 2)

    print(f"\n  sum_4={s4:.1f}, sum_6={s6:.1f}: {len(bits_list)} orientations")
    print(f"    #QR distribution: {sorted(set(qr_counts))} (values: {sorted(qr_counts)[:5]}...)")
    print(f"    |det|={det_val}")
    print(f"    |lam|: {', '.join(f'{x:.3f}' for x in abs_vals)}")

# Phase 2: Paley spectrum analysis
print(f"\n{'='*70}")
print(f"  PALEY SPECTRUM AT p={p}")
print(f"{'='*70}")

S_qr = sorted(S_qr_set)
eigs = eigenvalues(S_qr, p)
abs_vals = sorted([abs(e) for e in eigs], reverse=True)

print(f"  S_QR = {S_qr}")
print(f"  |lambda|: {', '.join(f'{x:.6f}' for x in abs_vals)}")

# Check flatness
non_trivial = abs_vals[1:]
all_equal = all(abs(v - non_trivial[0]) < 1e-6 for v in non_trivial)
print(f"  Flat spectrum? {all_equal}")
if all_equal:
    print(f"  All non-trivial |lambda| = {non_trivial[0]:.6f}")
    print(f"  Expected: sqrt(p)/2 = {math.sqrt(p)/2:.6f}")

# Phase 3: c3 constancy check
print(f"\n{'='*70}")
print(f"  c3 CONSTANCY AT p={p}")
print(f"{'='*70}")

# c3 formula: n(n-1)(n+1)/24 for regular tournaments
c3_formula = p * (p-1) * (p+1) // 24
print(f"  Formula: c3 = p(p-1)(p+1)/24 = {c3_formula}")

# Verify for a few orientations
def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_3cycles(A, p):
    total = 0
    for i in range(p):
        for j in range(i+1, p):
            for k in range(j+1, p):
                if (A[i][j]*A[j][k]*A[k][i]) + (A[i][k]*A[k][j]*A[j][i]) > 0:
                    total += 1
    return total

for bits in [0, 29, 31, 63]:
    S = bits_to_S(bits)
    A = build_adj(p, S)
    c3 = count_3cycles(A, p)
    print(f"  bits={bits}: c3 = {c3} (expected {c3_formula}), match={c3==c3_formula}")

# Phase 4: Quick c5 comparison for a few orientations
print(f"\n{'='*70}")
print(f"  c5 COMPARISON (selected orientations)")
print(f"{'='*70}")

def count_ham_cycles_5(A, verts):
    a, b, c, d, e = verts
    total = 0
    for perm in [(a,b,c,d,e), (a,b,c,e,d), (a,b,d,c,e), (a,b,d,e,c),
                 (a,b,e,c,d), (a,b,e,d,c), (a,c,b,d,e), (a,c,b,e,d),
                 (a,c,d,b,e), (a,c,d,e,b), (a,c,e,b,d), (a,c,e,d,b),
                 (a,d,b,c,e), (a,d,b,e,c), (a,d,c,b,e), (a,d,c,e,b),
                 (a,d,e,b,c), (a,d,e,c,b), (a,e,b,c,d), (a,e,b,d,c),
                 (a,e,c,b,d), (a,e,c,d,b), (a,e,d,b,c), (a,e,d,c,b)]:
        ok = True
        for i in range(5):
            if not A[perm[i]][perm[(i+1)%5]]:
                ok = False
                break
        if ok:
            total += 1
    return total

# Only for the spectral class representatives
reps = {}
for (s4, s6), bits_list in sorted(spectral_types.items()):
    reps[(s4, s6)] = bits_list[0]

print(f"  Computing c5 for {len(reps)} representative orientations...")
for (s4, s6), bits in sorted(reps.items()):
    S = bits_to_S(bits)
    A = build_adj(p, S)

    c5 = 0
    for subset in combinations(range(p), 5):
        c5 += count_ham_cycles_5(A, subset)

    print(f"  sum_4={s4:.1f}: bits={bits}, c5={c5}")

# Is -1 QR at p=13?
print(f"\n  -1 (={p-1}) is {'QR' if (p-1) in S_qr_set else 'NQR'} at p={p}")
print(f"  p mod 4 = {p % 4}")
# p=13: p mod 4 = 1, so -1 IS QR! This means Paley is NOT self-complementary.
# Wait: for p=1 mod 4, -1 is QR, so S_QR = -S_QR, meaning the Paley tournament
# IS self-complementary... no. Tournament is self-complementary if reversing all
# edges gives an isomorphic tournament. The anti-Paley has S_NQR.
# At p=1 mod 4: -1 is QR, so -S_QR = S_QR. This means T_QR and T_QR^op are
# the SAME tournament (sending i -> -i maps a->b to -a->-b, reversing edges
# since -(b-a) = a-b, and -1 being QR means a-b QR iff -(a-b) QR, i.e.,
# b-a QR. So the tournament is SELF-REVERSE, not self-complementary.)
# Actually, self-converse: T and T^op are isomorphic via i -> -i.

print(f"  p=1 mod 4 means Paley T_{p} is self-converse (T ~ T^op)")
print(f"  Anti-Paley has S=NQR = {sorted(set(range(1,p)) - S_qr_set)}")

# At p=13: reversal (complement) of S_QR under the pair structure
paley_bits = None
for bits in range(1 << m):
    if sorted(bits_to_S(bits)) == S_qr:
        paley_bits = bits
        break
anti_bits = (1 << m) - 1 - paley_bits
S_anti = sorted(bits_to_S(anti_bits))
print(f"  Paley bits={paley_bits} (0b{paley_bits:06b})")
print(f"  Anti-Paley bits={anti_bits} (0b{anti_bits:06b}), S={S_anti}")

print("\nDONE.")
