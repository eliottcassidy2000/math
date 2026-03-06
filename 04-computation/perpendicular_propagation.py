#!/usr/bin/env python3
"""
Perpendicular Propagation Analysis
====================================
Test the perpendicular maximizer principle across n=3..8:
- Self-converse (SC) tournaments have T ~ T^op
- SC tournaments sit perpendicular to the transitive<->full diagonal in tiling space
- Hypothesis: SC tournaments always have higher average H(T)
- Hypothesis: The GLOBAL maximizer of H(T) is always SC

Also test:
- How does the SC/NSC H-ratio change with n?
- Is there a pattern in which SC tournaments maximize H?
- Connection to Paley tournaments (always SC)

kind-pasteur-2026-03-05-S16
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, all_tournaments,
                             hamiltonian_path_count, opposite_tournament,
                             random_tournament, find_odd_cycles)
from itertools import permutations
import random

def is_isomorphic(T1, T2):
    """Check if T1 and T2 are isomorphic tournaments."""
    n = len(T1)
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if T1[i][j] != T2[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def is_self_converse(T):
    """Check if T is self-converse (isomorphic to T^op)."""
    Top = opposite_tournament(T)
    return is_isomorphic(T, Top)

def hamming_weight(bits, m):
    """Count set bits."""
    return bin(bits).count('1')

def score_sequence(T):
    """Return sorted score sequence (descending)."""
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return tuple(sorted(scores, reverse=True))

# ── Exhaustive analysis for n=3,4,5 ──
print("=" * 70)
print("PERPENDICULAR PROPAGATION ANALYSIS")
print("=" * 70)

for n in [3, 4, 5]:
    m = n * (n - 1) // 2
    sc_h = []
    nsc_h = []
    max_h = 0
    max_h_sc = False
    max_h_bits = None

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        sc = is_self_converse(T)

        if sc:
            sc_h.append(h)
        else:
            nsc_h.append(h)

        if h > max_h:
            max_h = h
            max_h_sc = sc
            max_h_bits = bits

    avg_sc = sum(sc_h) / len(sc_h) if sc_h else 0
    avg_nsc = sum(nsc_h) / len(nsc_h) if nsc_h else 0
    ratio = avg_sc / avg_nsc if avg_nsc > 0 else float('inf')

    print(f"\nn={n}: {1 << m} tournaments, {len(sc_h)} SC, {len(nsc_h)} NSC")
    print(f"  SC avg H = {avg_sc:.2f}, NSC avg H = {avg_nsc:.2f}, ratio = {ratio:.4f}")
    print(f"  SC max H = {max(sc_h)}, NSC max H = {max(nsc_h) if nsc_h else 'N/A'}")
    print(f"  Global max H = {max_h}, is SC? {max_h_sc}")
    print(f"  SC H distribution: {sorted(set(sc_h))}")
    if nsc_h:
        print(f"  NSC H distribution: {sorted(set(nsc_h))}")

    # Hamming weight analysis
    hw_by_h = {}
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        hw = hamming_weight(bits, m)
        if h not in hw_by_h:
            hw_by_h[h] = []
        hw_by_h[h].append(hw)

    print(f"  Hamming weight by H value:")
    for h_val in sorted(hw_by_h.keys()):
        hws = hw_by_h[h_val]
        avg_hw = sum(hws) / len(hws)
        print(f"    H={h_val}: avg HW={avg_hw:.2f}, range=[{min(hws)},{max(hws)}], count={len(hws)}")

# ── n=6 sampling ──
print(f"\n{'='*70}")
print("n=6: Sampling 2000 random tournaments")
print("=" * 70)

random.seed(42)
sc_h6 = []
nsc_h6 = []
max_h6 = 0
max_h6_sc = False

for _ in range(2000):
    T = random_tournament(6)
    h = hamiltonian_path_count(T)
    sc = is_self_converse(T)
    if sc:
        sc_h6.append(h)
    else:
        nsc_h6.append(h)
    if h > max_h6:
        max_h6 = h
        max_h6_sc = sc

avg_sc6 = sum(sc_h6) / len(sc_h6) if sc_h6 else 0
avg_nsc6 = sum(nsc_h6) / len(nsc_h6) if nsc_h6 else 0
ratio6 = avg_sc6 / avg_nsc6 if avg_nsc6 > 0 else float('inf')

print(f"  {len(sc_h6)} SC, {len(nsc_h6)} NSC")
print(f"  SC avg H = {avg_sc6:.2f}, NSC avg H = {avg_nsc6:.2f}, ratio = {ratio6:.4f}")
print(f"  Sample max H = {max_h6}, is SC? {max_h6_sc}")

# ── n=7 sampling ──
print(f"\n{'='*70}")
print("n=7: Sampling 500 random tournaments")
print("=" * 70)

sc_h7 = []
nsc_h7 = []
max_h7 = 0
max_h7_sc = False

for _ in range(500):
    T = random_tournament(7)
    h = hamiltonian_path_count(T)
    sc = is_self_converse(T)
    if sc:
        sc_h7.append(h)
    else:
        nsc_h7.append(h)
    if h > max_h7:
        max_h7 = h
        max_h7_sc = sc

avg_sc7 = sum(sc_h7) / len(sc_h7) if sc_h7 else 0
avg_nsc7 = sum(nsc_h7) / len(nsc_h7) if nsc_h7 else 0
ratio7 = avg_sc7 / avg_nsc7 if avg_nsc7 > 0 else float('inf')

print(f"  {len(sc_h7)} SC, {len(nsc_h7)} NSC")
print(f"  SC avg H = {avg_sc7:.2f}, NSC avg H = {avg_nsc7:.2f}, ratio = {ratio7:.4f}")
print(f"  Sample max H = {max_h7}, is SC? {max_h7_sc}")

# ── Paley tournaments ──
print(f"\n{'='*70}")
print("PALEY TOURNAMENT ANALYSIS")
print("=" * 70)

def paley_tournament(p):
    """Construct Paley tournament T_p for p ≡ 3 mod 4 prime."""
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j:
                if (j - i) % p in qr:
                    T[i][j] = 1
    return T

for p in [3, 7]:
    T = paley_tournament(p)
    h = hamiltonian_path_count(T)
    sc = is_self_converse(T)
    cycles = find_odd_cycles(T)

    # Count cycles by length
    len_counts = {}
    for c in cycles:
        L = len(c)
        len_counts[L] = len_counts.get(L, 0) + 1

    print(f"\nT_{p}: H={h}, SC={sc}, |Aut|={p*(p-1)//2}")
    print(f"  Odd cycles: {len_counts}")
    print(f"  Total cycles (Omega vertices): {len(cycles)}")

# ── Perpendicular structure in tiling space ──
print(f"\n{'='*70}")
print("PERPENDICULAR STRUCTURE IN TILING SPACE (n=5)")
print("=" * 70)

# For n=5, m=C(4,2)=6 tiles
# Transitive = all 0s (bits=0), Full cycle = all 1s (bits=63)
# "Perpendicular to diagonal" = Hamming distance m/2 from both
n = 5
m = n * (n - 1) // 2
mid = m / 2

print(f"m={m} tiles, midpoint Hamming weight = {mid}")

# For each tiling, compute: Hamming distance from 0 (transitive) and from 2^m-1 (full)
# Perpendicular = HW close to m/2
perp_h = []
diag_h = []

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    hw = hamming_weight(bits, m)
    dist_from_mid = abs(hw - mid)

    if dist_from_mid <= 0.5:  # HW exactly m/2 (if m even) or ±0.5 (if m odd)
        perp_h.append(h)
    elif hw <= 1 or hw >= m - 1:
        diag_h.append(h)

print(f"  Tilings at HW={m//2} or {(m+1)//2} (perpendicular): {len(perp_h)}, avg H = {sum(perp_h)/len(perp_h):.2f}")
print(f"  Tilings at HW<=1 or >={m-1} (near diagonal): {len(diag_h)}, avg H = {sum(diag_h)/len(diag_h):.2f}")
print(f"  Ratio (perp/diag): {sum(perp_h)/len(perp_h) / (sum(diag_h)/len(diag_h)):.4f}")

# More detailed: H as function of Hamming weight
print(f"\n  H vs Hamming weight:")
for hw in range(m + 1):
    h_vals = []
    for bits in range(1 << m):
        if hamming_weight(bits, m) == hw:
            T = tournament_from_bits(n, bits)
            h_vals.append(hamiltonian_path_count(T))
    avg = sum(h_vals) / len(h_vals)
    max_v = max(h_vals)
    print(f"    HW={hw}: count={len(h_vals)}, avg H={avg:.2f}, max H={max_v}")

# ── Self-converse + Hamming weight interaction ──
print(f"\n{'='*70}")
print("SC/NSC x HAMMING WEIGHT INTERACTION (n=5)")
print("=" * 70)

for hw in range(m + 1):
    sc_vals = []
    nsc_vals = []
    for bits in range(1 << m):
        if hamming_weight(bits, m) == hw:
            T = tournament_from_bits(n, bits)
            h = hamiltonian_path_count(T)
            if is_self_converse(T):
                sc_vals.append(h)
            else:
                nsc_vals.append(h)
    sc_avg = sum(sc_vals)/len(sc_vals) if sc_vals else 0
    nsc_avg = sum(nsc_vals)/len(nsc_vals) if nsc_vals else 0
    print(f"  HW={hw}: SC={len(sc_vals)} (avg {sc_avg:.1f}), NSC={len(nsc_vals)} (avg {nsc_avg:.1f})")

print("\nDone.")
