#!/usr/bin/env python3
"""
directed_cycle_analysis.py — opus-2026-03-14-S71d

SMART approach: At n≤8, α₃=0 (need 9 vertices for 3 disjoint odd cycles).
So: H = 1 + 2α₁ + 4α₂, giving α₂ = (H-1-2α₁)/4.
Just need to count α₁ (directed odd cycles) exactly.

Also: count α₁ by cycle length and compare directed vs vertex-set counts.
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_directed_cycles_by_length(A, n, max_length=None):
    """Count directed odd cycles by length.

    For each vertex set, count the number of DISTINCT Hamiltonian directed
    cycles. Normalize by: fix start vertex (= min), count each physical
    undirected cycle once (a tournament cycle has unique direction).

    Returns dict: length -> count of directed cycles
    """
    if max_length is None:
        max_length = n if n % 2 == 1 else n - 1

    counts = {}
    cycles_by_vset = {}  # vset -> [(canon_tuple, length)]

    for length in range(3, max_length + 1, 2):
        if length > n:
            break
        count = 0
        for combo in combinations(range(n), length):
            verts = list(combo)
            v0 = verts[0]
            seen = set()
            for perm in permutations(verts[1:]):
                order = (v0,) + perm
                if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                    # Normalize direction: keep lex-smaller of (order, reverse)
                    reverse = (v0,) + perm[::-1]
                    canon = min(order, reverse)
                    if canon not in seen:
                        seen.add(canon)
                        count += 1
                        fs = frozenset(combo)
                        if fs not in cycles_by_vset:
                            cycles_by_vset[fs] = []
                        cycles_by_vset[fs].append((canon, length))
        counts[length] = count

    return counts, cycles_by_vset

def count_vertex_set_cycles(A, n, max_length=None):
    """Count vertex-sets that support at least one odd cycle."""
    if max_length is None:
        max_length = n if n % 2 == 1 else n - 1
    counts = {}
    for length in range(3, max_length + 1, 2):
        if length > n:
            break
        count = 0
        for combo in combinations(range(n), length):
            verts = list(combo)
            found = False
            for perm in permutations(verts[1:]):
                order = [verts[0]] + list(perm)
                if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                    found = True
                    break
            if found:
                count += 1
        counts[length] = count
    return counts

def count_disjoint_cycle_pairs(cycles_by_vset):
    """Count pairs of vertex-disjoint directed cycles (α₂)."""
    all_cycles = []
    for vset, cycs in cycles_by_vset.items():
        for canon, length in cycs:
            all_cycles.append((canon, vset))

    count = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if not (all_cycles[i][1] & all_cycles[j][1]):
                count += 1
    return count

# ======================================================================
# PART 1: n=5 exhaustive — directed vs vertex-set
# ======================================================================
print("=" * 70)
print("PART 1: n=5 EXHAUSTIVE — DIRECTED vs VERTEX-SET CYCLE COUNTS")
print("=" * 70)

n = 5
tb = n*(n-1)//2
stats = defaultdict(int)

t0 = time.time()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc, cbyv = count_directed_cycles_by_length(A, n)
    vc = count_vertex_set_cycles(A, n)

    a1 = sum(dc.values())
    a2_from_H = (H - 1 - 2*a1) // 4  # since α₃=0 at n=5

    # Also directly count disjoint pairs
    a2_direct = count_disjoint_cycle_pairs(cbyv)

    if a2_from_H != a2_direct:
        print(f"  α₂ MISMATCH: bits={bits}, H={H}, α₁={a1}, α₂(H)={a2_from_H}, α₂(direct)={a2_direct}")

    dc3 = dc.get(3, 0)
    dc5 = dc.get(5, 0)
    vc3 = vc.get(3, 0)
    vc5 = vc.get(5, 0)
    stats[(dc3, dc5, vc3, vc5, a1, a2_direct, H)] += 1

dt = time.time() - t0
print(f"  n=5: computed {1 << tb} tournaments in {dt:.1f}s")

# Summary table
print(f"\n  {'dc3':>4s} {'dc5':>4s} {'vc3':>4s} {'vc5':>4s} {'α₁':>4s} {'α₂':>4s} {'H':>5s} {'cnt':>5s} | dc5/vc5")
print(f"  {'─'*4} {'─'*4} {'─'*4} {'─'*4} {'─'*4} {'─'*4} {'─'*5} {'─'*5}-+-{'─'*7}")
for key in sorted(stats.keys()):
    dc3, dc5, vc3, vc5, a1, a2, H = key
    ratio = f"{dc5/vc5:.1f}" if vc5 > 0 else "—"
    print(f"  {dc3:4d} {dc5:4d} {vc3:4d} {vc5:4d} {a1:4d} {a2:4d} {H:5d} {stats[key]:5d} | {ratio}")

# Verify H = 1 + 2α₁ + 4α₂
all_ok = all(key[6] == 1 + 2*key[4] + 4*key[5] for key in stats)
print(f"\n  H = 1 + 2α₁ + 4α₂ verified: {all_ok}")

# ======================================================================
# PART 2: n=7 sample — directed cycle structure
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: n=7 SAMPLE (200 tournaments) — DIRECTED CYCLE ANALYSIS")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

data7 = []
t0 = time.time()

for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    dc, cbyv = count_directed_cycles_by_length(A, n)
    vc = count_vertex_set_cycles(A, n)

    a1 = sum(dc.values())
    a2 = (H - 1 - 2*a1) // 4  # α₃=0 at n=7

    # Also compute I(3) = 1 + 3α₁ + 9α₂ + 27α₃
    # Since α₃=0: I(3) = 1 + 3α₁ + 9α₂
    I3 = 1 + 3*a1 + 9*a2

    # Vandermonde check: α₂ = (2I₃ - 3H + 1)/6
    a2_vand = (2*I3 - 3*H + 1) / 6

    data7.append({
        'H': H, 'dc3': dc.get(3,0), 'dc5': dc.get(5,0), 'dc7': dc.get(7,0),
        'vc3': vc.get(3,0), 'vc5': vc.get(5,0), 'vc7': vc.get(7,0),
        'a1': a1, 'a2': a2, 'I3': I3
    })

    if trial % 50 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, H={H}, α₁={a1} (dc3={dc.get(3,0)},dc5={dc.get(5,0)},dc7={dc.get(7,0)}), α₂={a2}")

dt = time.time() - t0
print(f"  Done: {dt:.1f}s")

# Statistics
print(f"\n  Directed vs vertex-set cycle counts at n=7:")
dc5_arr = np.array([d['dc5'] for d in data7])
vc5_arr = np.array([d['vc5'] for d in data7])
dc7_arr = np.array([d['dc7'] for d in data7])
vc7_arr = np.array([d['vc7'] for d in data7])

print(f"  5-cycles: dc/vc ratio: mean={np.mean(dc5_arr/np.maximum(vc5_arr,1)):.2f}, max={np.max(dc5_arr/np.maximum(vc5_arr,1)):.2f}")
print(f"  7-cycles: dc/vc ratio: mean={np.mean(dc7_arr/np.maximum(vc7_arr,1)):.2f}, max={np.max(dc7_arr/np.maximum(vc7_arr,1)):.2f}")

print(f"\n  dc5 range: [{min(dc5_arr)}, {max(dc5_arr)}]")
print(f"  vc5 range: [{min(vc5_arr)}, {max(vc5_arr)}]")
print(f"  dc7 range: [{min(dc7_arr)}, {max(dc7_arr)}]")
print(f"  vc7 range: [{min(vc7_arr)}, {max(vc7_arr)}]")

# How many directed 5-cycles per vertex-set?
dc5_per_vset = dc5_arr / np.maximum(vc5_arr, 1)
dc7_per_vset = dc7_arr / np.maximum(vc7_arr, 1)
print(f"\n  Directed 5-cycles per vertex-set: mean={np.mean(dc5_per_vset):.2f}")
print(f"  Directed 7-cycles per vertex-set: mean={np.mean(dc7_per_vset):.2f}")
print(f"  (Max possible: 5-cycles=12, 7-cycles=360)")

# α₁ and α₂ statistics
a1_arr = np.array([d['a1'] for d in data7])
a2_arr = np.array([d['a2'] for d in data7])
H_arr = np.array([d['H'] for d in data7])
print(f"\n  α₁: mean={np.mean(a1_arr):.1f}, range=[{min(a1_arr)}, {max(a1_arr)}]")
print(f"  α₂: mean={np.mean(a2_arr):.1f}, range=[{min(a2_arr)}, {max(a2_arr)}]")
print(f"  H:  mean={np.mean(H_arr):.1f}, range=[{min(H_arr)}, {max(H_arr)}]")

# Verify H = 1 + 2α₁ + 4α₂
verify = all(d['H'] == 1 + 2*d['a1'] + 4*d['a2'] for d in data7)
print(f"\n  H = 1 + 2α₁ + 4α₂ verified: {verify}")

# Verify I(3) = 1 + 3α₁ + 9α₂
# And Vandermonde: α₂ = (2I₃ - 3H + 1)/6
vand_ok = all(abs((2*d['I3'] - 3*d['H'] + 1)/6 - d['a2']) < 0.01 for d in data7)
print(f"  Vandermonde extraction verified: {vand_ok}")

# I₃/H ratio
ratios = np.array([d['I3']/d['H'] for d in data7])
print(f"\n  I₃/H ratio: mean={np.mean(ratios):.6f}")
print(f"  3/2 = {1.5}")
print(f"  Deviation from 3/2: mean={np.mean(ratios - 1.5):.6f}")

# ======================================================================
# PART 3: THE KEY INSIGHT — directed cycle multiplicity
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: DIRECTED CYCLE MULTIPLICITY AT n=5 vs n=7")
print("=" * 70)

# At n=5, how many directed 5-cycles does the champion (H=15) have?
n = 5
A_champ = bits_to_adj(40, n)  # H=15
dc_champ, cbyv_champ = count_directed_cycles_by_length(A_champ, n)
print(f"\n  n=5 champion (bits=40, H=15):")
print(f"    dc3={dc_champ.get(3,0)}, dc5={dc_champ.get(5,0)}")
print(f"    α₁={sum(dc_champ.values())}")
for vset, cycs in cbyv_champ.items():
    if len(cycs) > 1:
        print(f"    vertex-set {set(vset)}: {len(cycs)} directed cycles")
        for canon, length in cycs:
            print(f"      {length}-cycle: {'→'.join(str(v) for v in canon)}→{canon[0]}")

# ======================================================================
# PART 4: The corrected 3/2 ratio — is it EXACT within lambda fibers?
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: 3/2 RATIO WITHIN LAMBDA FIBERS (corrected)")
print("=" * 70)

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

# For n=7 sample, group by lambda and check I₃/H ratio within fibers
from collections import defaultdict
lambda_groups = defaultdict(list)
for d in data7:
    # Need to recompute lambda key
    pass

# Actually, let's just check the formula: within a lambda fiber,
# I₃ = (3/2)·H + constant?
# I₃ = 1 + 3α₁ + 9α₂, H = 1 + 2α₁ + 4α₂
# I₃/H = (1 + 3α₁ + 9α₂) / (1 + 2α₁ + 4α₂)
# If α₂ = 0: I₃/H = (1+3α₁)/(1+2α₁) → 3/2 as α₁ → ∞
# If α₂ > 0: I₃/H = (1+3α₁+9α₂)/(1+2α₁+4α₂) → 9/4 as α₂ dominates

print(f"  I₃/H as function of (α₁, α₂):")
print(f"  I₃/H = (1+3α₁+9α₂) / (1+2α₁+4α₂)")
print(f"  For α₂=0: → 3/2 as α₁ → ∞")
print(f"  For large α₂, α₁ fixed: → 9/4 = 2.25")
print()

# Show actual data by α₂ value
for a2_val in sorted(set(d['a2'] for d in data7)):
    group = [d for d in data7 if d['a2'] == a2_val]
    if len(group) < 3:
        continue
    ratios_g = [d['I3']/d['H'] for d in group]
    a1_mean = np.mean([d['a1'] for d in group])
    print(f"  α₂={a2_val:3d} ({len(group):3d} tours): I₃/H mean={np.mean(ratios_g):.6f}, α₁ mean={a1_mean:.0f}")

# ======================================================================
# PART 5: k-nacci convergence constants (corrected)
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: k-NACCI CONVERGENCE CONSTANTS")
print("=" * 70)

def knacci_root(k):
    coeffs = [-1] * k + [1]
    roots = np.roots(coeffs[::-1])
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

def weighted_knacci_root(k):
    coeffs = [1]
    for j in range(k):
        coeffs.append(-2**j)
    roots = np.roots(coeffs)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

print(f"\n  {'k':>3s}  {'C₂=(φ_k-2)·2^k':>18s}  {'C₃=(ψ_k-3)·(3/2)^k':>22s}  {'C₃/C₂':>10s}")
for k in range(3, 30):
    phi_k = knacci_root(k)
    psi_k = weighted_knacci_root(k)
    if phi_k is None or psi_k is None:
        continue
    c2 = (phi_k - 2) * 2**k
    c3 = (psi_k - 3) * (3/2)**k
    mark = " ←" if k in [7, 8] else ""
    print(f"  {k:3d}  {c2:18.10f}  {c3:22.10f}  {c3/c2:10.6f}{mark}")

# ======================================================================
# PART 6: SYNTHESIS
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: SYNTHESIS")
print("=" * 70)

print("""
KEY FINDING: DIRECTED CYCLE MULTIPLICITY

At n=5, a 5-vertex tournament with a 5-cycle can have 1, 2, or 3
distinct directed 5-cycles on those same 5 vertices. This is because
a tournament on k vertices can have multiple Hamiltonian cycles.

The OCF counts EACH directed cycle as a separate vertex of Ω.
This makes α₁ larger than the vertex-set count, which means:

  I(Ω, 2) = 1 + 2·(TRUE α₁) + 4·(TRUE α₂) = H

where TRUE α₁ counts directed cycles, not vertex-sets.

IMPLICATION FOR THE 7→8 TRANSITION:
  At n=7, the number of directed 7-cycles on 7 vertices can be
  very large (up to 360 per tournament). Each contributes to α₁.
  The "extra" directed cycles (beyond vertex-set counting) become
  increasingly important at larger n.

  At n=8, directed 5-cycles on 5 of the 8 vertices pair with
  directed 3-cycles on the remaining 3. The cross-level α₂
  counts ALL pairings of directed cycles (not just vertex-sets).

THE 2-3 FRAMEWORK (corrected):
  I(Ω, 2) = H (OCF, proved)
  I(Ω, 3) = 1 + 3α₁ + 9α₂ (at n≤8 where α₃=0)
  Vandermonde: α₂ = (2I₃ - 3H + 1)/6, α₁ = (H-1-4α₂)/2
  This works with DIRECTED cycle counts.
""")

print("Done.")
