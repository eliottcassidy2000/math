#!/usr/bin/env python3
"""
Deletion-Sum Ratio Theorem Investigation.

Discovery: The H-maximizer MINIMIZES the ratio R(T) = sum_v H(T-v) / H(T).

This means: among all tournaments with H(T) > 0, the maximizer has the
smallest ratio of "descendant paths" to "own paths."

Questions:
1. Is this ALWAYS true (all n)?
2. Does the minimizer of R also maximize H?
3. What is the exact value of R for the maximizer?
4. Is there a clean formula for R in terms of n?

Observed: R_max ≈ n/3 at odd n? (n=5: 5/3=1.667, n=7: 7/3=2.333... no)
Actually n=5: range [1.4, 1.667], n=6: range [1.467, 1.733]

Let me compute exactly for n=3,4,5,6,7.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

def deletion_sum(T):
    return sum(hamiltonian_path_count(delete_vertex(T, v)) for v in range(len(T)))

# ============================================================
# Exact R values for maximizers at each n
# ============================================================
print("=" * 70)
print("DELETION-SUM RATIO FOR MAXIMIZERS")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    # Find maximizers and compute their R values
    r_values = []
    count = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == max_h:
            ds = deletion_sum(T)
            r = ds / h
            r_values.append(r)
            count += 1
            if count >= 10:  # Don't compute all at large n
                break

    if r_values:
        r_min = min(r_values)
        r_max = max(r_values)
        print(f"n={n}: H_max={max_h}, R = {r_min:.6f}" +
              (f" to {r_max:.6f}" if r_max > r_min + 0.0001 else "") +
              f" (sum={int(r_min * max_h)})")

# ============================================================
# Is R minimized by the maximizer? Full check for n=3,4,5
# ============================================================
print(f"\n{'='*70}")
print("IS R MINIMIZED BY H-MAXIMIZER? (n=3,4,5)")
print("=" * 70)

for n in [3, 4, 5]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    all_r = []  # (R, H, bits)
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == 0:
            continue
        ds = deletion_sum(T)
        r = ds / h
        all_r.append((r, h, bits))

    all_r.sort()  # Sort by R
    min_r = all_r[0][0]
    max_r = all_r[-1][0]

    # Which H values achieve minimum R?
    min_r_hs = set(h for r, h, b in all_r if abs(r - min_r) < 0.0001)

    print(f"\nn={n}: R range = [{min_r:.4f}, {max_r:.4f}]")
    print(f"  Min R achieved by H values: {sorted(min_r_hs, reverse=True)}")
    print(f"  Max H = {max_h}, in min R set? {max_h in min_r_hs}")

    # Show the R distribution by H
    from collections import defaultdict
    r_by_h = defaultdict(list)
    for r, h, b in all_r:
        r_by_h[h].append(r)

    print(f"  R by H (top 5 H values):")
    for h in sorted(r_by_h.keys(), reverse=True)[:5]:
        rs = r_by_h[h]
        print(f"    H={h}: R min={min(rs):.4f}, max={max(rs):.4f}, avg={sum(rs)/len(rs):.4f}")

# ============================================================
# Full check at n=6
# ============================================================
print(f"\n{'='*70}")
print("FULL CHECK AT n=6")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
max_h = MAX_H[n]

r_by_h = defaultdict(list)
count = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h == 0:
        continue
    ds = deletion_sum(T)
    r = ds / h
    r_by_h[h].append(r)
    count += 1

print(f"n=6: {count} tournaments with H>0")
print(f"R by H (top 10 H values):")
global_min_r = float('inf')
global_min_h = None
for h in sorted(r_by_h.keys(), reverse=True)[:10]:
    rs = r_by_h[h]
    r_min = min(rs)
    if r_min < global_min_r:
        global_min_r = r_min
        global_min_h = h
    marker = " *** H-MAX" if h == max_h else ""
    marker2 = " <<< R-MIN" if r_min <= global_min_r else ""
    print(f"  H={h}: R min={r_min:.4f}, max={max(rs):.4f}, count={len(rs)}{marker}{marker2}")

print(f"\nGlobal min R = {global_min_r:.4f}, achieved by H = {global_min_h}")
print(f"H-maximizer minimizes R? {global_min_h == max_h}")

# ============================================================
# The exact deletion sums
# ============================================================
print(f"\n{'='*70}")
print("EXACT DELETION SUMS FOR MAXIMIZERS")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    max_h = MAX_H[n]
    found = False
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == max_h:
            ds = deletion_sum(T)
            del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
            print(f"n={n}: H={h}, sum={ds}, R={ds/h:.6f}, del_hs={del_hs}")
            found = True
            break

# ============================================================
# Can we express R in terms of n?
# Look for pattern in R_min values
# ============================================================
print(f"\n{'='*70}")
print("PATTERN IN R_min")
print("=" * 70)

# From above: R values for maximizers
# n=3: H=3, sum=3, R=1
# n=4: H=5, sum=8, R=1.6
# n=5: H=15, sum=?, R=?
# Need to collect these

r_data = {}
for n in range(3, 8):
    m = n * (n - 1) // 2
    max_h = MAX_H[n]
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == max_h:
            ds = deletion_sum(T)
            r_data[n] = (max_h, ds, ds/max_h)
            break

print("n | H_max | sum | R = sum/H | n*H_{n-1}/H_n | (n-1)*(n-2)/n")
for n in sorted(r_data.keys()):
    h, s, r = r_data[n]
    h_prev = MAX_H.get(n-1, 1)
    pred1 = n * h_prev / h
    pred2 = (n-1) * (n-2) / n
    print(f"  {n} | {h:>5} | {s:>5} | {r:.6f} | {pred1:.6f} | {pred2:.6f}")

# Check: is R = n * H_{n-1} / H_n for the maximizer?
print(f"\nIs R = n * H_{{n-1}} / H_n?")
for n in sorted(r_data.keys()):
    h, s, r = r_data[n]
    h_prev = MAX_H.get(n-1, 1)
    pred = n * h_prev / h
    match = abs(r - pred) < 0.001
    print(f"  n={n}: R={r:.6f}, n*H_{{n-1}}/H_n = {pred:.6f}, match={match}")

print("\nDone.")
