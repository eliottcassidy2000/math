#!/usr/bin/env python3
"""
blue_line_n6_s20x.py -- kind-pasteur-2026-03-22-S20x

Extend the blue line skeleton analysis to n=6.
At n=6 there are 32768 tournaments. SC check is expensive (6! perms),
so we sample or use a faster SC test.

Also: the ARBORESCENCE GAP structure -- which arb values are missing?
And the SYMMETRY-BREAKING PHASE TRANSITION: as we continuously deform
from SC to non-SC (by flipping arcs), how do HC and L change?

Author: kind-pasteur-2026-03-22-S20x
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_hc(A, n):
    full = (1 << n) - 1
    dp = defaultdict(int)
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[(full, v)] for v in range(n) if A[v][0])

def count_arborescences(A, n, root=0):
    A_float = A.astype(float)
    D_in = np.diag(A_float.sum(axis=0))
    L = D_in - A_float.T
    indices = [i for i in range(n) if i != root]
    L_minor = L[np.ix_(indices, indices)]
    return int(round(np.linalg.det(L_minor)))

def count_kings(A, n):
    kings = 0
    A2 = A @ A
    reach = A + np.clip(A2, 0, 1)
    np.fill_diagonal(reach, 1)
    for v in range(n):
        if all(reach[v][w] > 0 or v == w for w in range(n)):
            kings += 1
    return kings

def is_sc(A, n):
    """Check if tournament is self-complementary."""
    A_comp = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j:
                A_comp[i][j] = 1 - A[i][j]
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if A[perm[i]][perm[j]] != A_comp[i][j]:
                    match = False
                    break
            if not match: break
        if match:
            return True
    return False

print("=" * 65)
print("  BLUE LINE SKELETON: n=6 EXTENSION")
print("=" * 65)

# ================================================================
# n=6: EXHAUSTIVE (32768 tournaments)
# SC check at n=6 takes 720 perm checks per tournament
# This is feasible but slow. Let's do it.
# ================================================================
n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)  # 15

print(f"\n  Computing n={n} ({2**m} tournaments)...")
print(f"  SC check: {n}! = 720 permutations per tournament")

all_data = []
sc_count = 0
batch = 2**m

for bits in range(batch):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1

    s = A.sum(axis=1).astype(int)
    S2 = int(sum(s*s))
    H = count_hp(A, n)
    HC = count_hc(A, n)
    L = H - n * HC
    c3 = comb(n,3) - (S2 - comb(n,2)) // 2
    arb = count_arborescences(A, n, root=0)
    kings = count_kings(A, n)
    score = tuple(sorted(s))

    # SC check is the bottleneck -- only do it for interesting cases
    # For speed: only check SC for score sequences that CAN be SC
    # SC score sequences satisfy s_i + s_{n-1-i} = n-1
    can_be_sc = all(score[i] + score[n-1-i] == n-1 for i in range(n//2))
    sc_flag = False
    if can_be_sc:
        sc_flag = is_sc(A, n)
    if sc_flag: sc_count += 1

    all_data.append({
        'bits': bits, 'H': H, 'HC': HC, 'L': L, 'c3': c3,
        'S2': S2, 'arb': arb, 'kings': kings,
        'sc': sc_flag, 'score': score, 'can_be_sc': can_be_sc
    })

    if (bits + 1) % 5000 == 0:
        print(f"    processed {bits+1}/{batch}...")

print(f"  Done. {sc_count} SC tournaments out of {batch} total.")

# ================================================================
# 1. H-SPECTRUM: SC vs NON-SC
# ================================================================
print(f"\n  1. H-SPECTRUM: SC vs NON-SC")
print()

sc_data = [d for d in all_data if d['sc']]
nsc_data = [d for d in all_data if not d['sc']]

sc_H_dist = defaultdict(int)
nsc_H_dist = defaultdict(int)
for d in sc_data: sc_H_dist[d['H']] += 1
for d in nsc_data: nsc_H_dist[d['H']] += 1

all_H = sorted(set(d['H'] for d in all_data))
print(f"  {'H':>4s} {'SC#':>6s} {'NSC#':>6s} {'type':>12s}")
for H_val in all_H:
    sc_c = sc_H_dist.get(H_val, 0)
    nsc_c = nsc_H_dist.get(H_val, 0)
    if sc_c > 0 and nsc_c > 0: t = "SHARED"
    elif sc_c > 0: t = "SC-ONLY"
    else: t = "NSC-ONLY"
    print(f"  {H_val:>4d} {sc_c:>6d} {nsc_c:>6d} {t:>12s}")

# ================================================================
# 2. (HC, L) COMPARISON
# ================================================================
print(f"\n  2. (H, HC, L) BY SC STATUS")
print()

sc_HL = defaultdict(int)
nsc_HL = defaultdict(int)
for d in all_data:
    key = (d['H'], d['HC'], d['L'])
    if d['sc']: sc_HL[key] += 1
    else: nsc_HL[key] += 1

all_keys = sorted(set(list(sc_HL.keys()) + list(nsc_HL.keys())))
print(f"  {'H':>4s} {'HC':>4s} {'L':>4s} {'SC#':>6s} {'NSC#':>6s}")
for (H, HC, L) in all_keys:
    sc_c = sc_HL.get((H, HC, L), 0)
    nsc_c = nsc_HL.get((H, HC, L), 0)
    if sc_c > 0 or nsc_c > 0:
        print(f"  {H:>4d} {HC:>4d} {L:>4d} {sc_c:>6d} {nsc_c:>6d}")

# ================================================================
# 3. ARBORESCENCE STRUCTURE
# ================================================================
print(f"\n  3. ARBORESCENCE STRUCTURE")
print()

sc_arbs = sorted(set(d['arb'] for d in sc_data))
nsc_arbs = sorted(set(d['arb'] for d in nsc_data))
all_arbs = sorted(set(d['arb'] for d in all_data))

# Gap analysis
arb_set = set(all_arbs)
arb_max = max(all_arbs)
missing = [a for a in range(arb_max + 1) if a not in arb_set]
print(f"  Total arb values: {len(all_arbs)} distinct (range [0, {arb_max}])")
print(f"  Missing arb values: {missing[:30]}{'...' if len(missing) > 30 else ''}")
print(f"  SC-only arb values: {sorted(set(sc_arbs) - set(nsc_arbs))[:20]}")
print(f"  NSC-only arb values: {sorted(set(nsc_arbs) - set(sc_arbs))[:20]}")

# ================================================================
# 4. KINGS STRUCTURE
# ================================================================
print(f"\n  4. KINGS: SC vs NON-SC")
print()

for k_val in sorted(set(d['kings'] for d in all_data)):
    sc_k = [d for d in sc_data if d['kings'] == k_val]
    nsc_k = [d for d in nsc_data if d['kings'] == k_val]
    sc_avg_H = sum(d['H'] for d in sc_k) / len(sc_k) if sc_k else 0
    nsc_avg_H = sum(d['H'] for d in nsc_k) / len(nsc_k) if nsc_k else 0
    print(f"  kings={k_val}: SC({len(sc_k)}, avg H={sc_avg_H:.1f}), NSC({len(nsc_k)}, avg H={nsc_avg_H:.1f})")

# ================================================================
# 5. SCORE CLASSES ON THE BLUE LINE
# ================================================================
print(f"\n  5. SCORE CLASSES ON THE BLUE LINE")
print()

score_sc = defaultdict(list)
score_nsc = defaultdict(list)
for d in all_data:
    if d['sc']: score_sc[d['score']].append(d)
    else: score_nsc[d['score']].append(d)

for score in sorted(set(d['score'] for d in all_data)):
    sc_list = score_sc.get(score, [])
    nsc_list = score_nsc.get(score, [])
    if sc_list:
        sc_Hs = sorted(set(d['H'] for d in sc_list))
        sc_HCs = sorted(set(d['HC'] for d in sc_list))
        print(f"  {list(score)}: SC({len(sc_list)}) H={sc_Hs}, HC={sc_HCs}")
    if nsc_list and not sc_list:
        nsc_Hs = sorted(set(d['H'] for d in nsc_list))
        print(f"  {list(score)}: NSC-only({len(nsc_list)}) H={nsc_Hs}")

# ================================================================
# 6. SYMMETRY-BREAKING: ARC FLIP FROM SC TO NON-SC
# ================================================================
print(f"\n  6. SYMMETRY-BREAKING PHASE TRANSITION")
print()

# Take a specific SC tournament (the regular one at n=5 embedded) and flip arcs
# Use the regular tournament: 0->1->2->3->4->0, 0->2, 1->3, 2->4, 3->0, 4->1 (the Paley T_5)
# Extended to n=6 with vertex 5

# For n=6, start with an SC tournament and see what happens when we break it
# First find an SC tournament with max H
if sc_data:
    sc_max_H = max(sc_data, key=lambda d: d['H'])
    print(f"  SC max H: H={sc_max_H['H']}, HC={sc_max_H['HC']}, L={sc_max_H['L']}, arb={sc_max_H['arb']}, kings={sc_max_H['kings']}, score={sc_max_H['score']}")

    # Find the bits and reconstruct
    bits = sc_max_H['bits']
    A_sc = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A_sc[i][j] = 1
        else: A_sc[j][i] = 1

    print(f"  SC max tournament adjacency:")
    for i in range(n):
        row = ''.join(str(A_sc[i][j]) for j in range(n))
        print(f"    {row}")

    # Now flip each arc one at a time and see what happens
    print(f"\n  Effect of SINGLE ARC FLIP from SC max:")
    print(f"  {'arc':>8s} {'H_new':>6s} {'dH':>4s} {'HC_new':>6s} {'L_new':>5s} {'arb_new':>7s} {'still_SC':>9s}")
    for i in range(n):
        for j in range(n):
            if i == j or A_sc[i][j] == 0: continue
            A_flip = A_sc.copy()
            A_flip[i][j] = 0
            A_flip[j][i] = 1
            H_new = count_hp(A_flip, n)
            HC_new = count_hc(A_flip, n)
            L_new = H_new - n * HC_new
            arb_new = count_arborescences(A_flip, n, root=0)
            sc_new = is_sc(A_flip, n) if all(
                tuple(sorted(A_flip.sum(axis=1).astype(int)))[k] +
                tuple(sorted(A_flip.sum(axis=1).astype(int)))[n-1-k] == n-1
                for k in range(n//2)) else False
            dH = H_new - sc_max_H['H']
            print(f"  ({i}->{j}): {H_new:>6d} {dH:>+4d} {HC_new:>6d} {L_new:>5d} {arb_new:>7d} {str(sc_new):>9s}")

# ================================================================
# 7. THE GLOBAL SUMMARY
# ================================================================
print(f"\n  7. GLOBAL SUMMARY")
print()

max_H = max(d['H'] for d in all_data)
min_H = min(d['H'] for d in all_data)
sc_max = max(d['H'] for d in sc_data) if sc_data else 0
nsc_max = max(d['H'] for d in nsc_data) if nsc_data else 0
sc_min = min(d['H'] for d in sc_data) if sc_data else 0
nsc_min = min(d['H'] for d in nsc_data) if nsc_data else 0

print(f"  Global H range: [{min_H}, {max_H}]")
print(f"  SC H range: [{sc_min}, {sc_max}]")
print(f"  NSC H range: [{nsc_min}, {nsc_max}]")
print(f"  SC count: {len(sc_data)} ({100*len(sc_data)/len(all_data):.1f}%)")
print(f"  NSC count: {len(nsc_data)} ({100*len(nsc_data)/len(all_data):.1f}%)")
print()

# SC-exclusive and NSC-exclusive H values
sc_H_set = set(d['H'] for d in sc_data)
nsc_H_set = set(d['H'] for d in nsc_data)
print(f"  SC-exclusive H values: {sorted(sc_H_set - nsc_H_set)}")
print(f"  NSC-exclusive H values: {sorted(nsc_H_set - sc_H_set)}")
print(f"  Shared H values: {sorted(sc_H_set & nsc_H_set)}")
print()

# The key: at n=6, does SC still span the full range?
# Does NSC still avoid Hamiltonian cycles?
nsc_hc_vals = sorted(set(d['HC'] for d in nsc_data))
sc_hc_vals = sorted(set(d['HC'] for d in sc_data))
print(f"  SC HC values: {sc_hc_vals}")
print(f"  NSC HC values: {nsc_hc_vals}")
print()

# Forbidden H values
all_H_set = set(d['H'] for d in all_data)
full_range = set(range(min_H, max_H + 1))
forbidden = sorted(full_range - all_H_set)
# Only odd values (Redei)
forbidden_odd = [h for h in forbidden if h % 2 == 1]
print(f"  Forbidden odd H values: {forbidden_odd}")
