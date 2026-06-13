#!/usr/bin/env python3
"""
swap_sequences_s192.py — New Sequences from the Swap Graph
opus-2026-03-22-S192

FROM S191:
  Within-fiber edge fraction: 1/2, 3/8, 5/16 at n=3,4,5
  Total fiber components: 3, 18, 117 at n=3,4,5

This script:
1. Verifies these sequences and extends to n=6 (exhaustive)
2. Checks within-fiber fraction numerators/denominators for patterns
3. Computes the LINEAR H-FORMULA at n=6
4. Finds new OEIS-checkable sequences from the fiber structure
5. Investigates the (2n-1)/2^n conjecture for within-fiber fraction
"""

from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
from fractions import Fraction
import numpy as np
import time

print("=" * 72)
print("  NEW SEQUENCES FROM THE SWAP GRAPH")
print("  opus-2026-03-22-S192")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def score_seq_labeled(adj, n):
    return tuple(sum(adj[i][j] for j in range(n) if j != i) for i in range(n))

# ==========================================================================
# PART 1: WITHIN-FIBER FRACTION — EXACT COMPUTATION
# ==========================================================================

print("PART 1: WITHIN-FIBER EDGE FRACTION — EXACT")
print("-" * 60)
print()

# Compute for n=3,4,5
fractions_list = []
component_counts = []

for n in range(3, 6):
    m = comb(n, 2)
    t0 = time.time()

    # Classify all tournaments
    b_to_s = {}
    for bits in range(1 << m):
        adj = tournament_adj(n, bits)
        ss = score_seq(adj, n)
        b_to_s[bits] = ss

    within = 0
    total = 0
    for bits in range(1 << m):
        for bit in range(m):
            nbr = bits ^ (1 << bit)
            total += 1
            if b_to_s[bits] == b_to_s[nbr]:
                within += 1

    frac = Fraction(within, total)
    fractions_list.append((n, frac))

    # Components
    s_to_b = defaultdict(list)
    for bits, ss in b_to_s.items():
        s_to_b[ss].append(bits)

    total_comps = 0
    for ss in s_to_b:
        fiber = s_to_b[ss]
        fiber_set = set(fiber)
        remaining = set(fiber)
        while remaining:
            start = remaining.pop()
            queue = deque([start])
            comp = {start}
            while queue:
                curr = queue.popleft()
                for bit in range(m):
                    nbr = curr ^ (1 << bit)
                    if nbr in remaining:
                        remaining.discard(nbr)
                        comp.add(nbr)
                        queue.append(nbr)
            total_comps += 1

    component_counts.append((n, total_comps))
    elapsed = time.time() - t0
    print(f"  n={n}: fraction = {frac} = {float(frac):.6f}, "
          f"components = {total_comps}, time = {elapsed:.1f}s")

# Now n=6 (2^15 = 32768 tournaments)
n = 6
m = comb(n, 2)  # 15
t0 = time.time()
print(f"\n  Computing n={n} ({2**m} tournaments, m={m} arcs)...")

b_to_s = {}
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    ss = score_seq(adj, n)
    b_to_s[bits] = ss

within = 0
total = 0
for bits in range(1 << m):
    for bit in range(m):
        nbr = bits ^ (1 << bit)
        total += 1
        if b_to_s[bits] == b_to_s[nbr]:
            within += 1

frac = Fraction(within, total)
fractions_list.append((n, frac))

# Components at n=6
s_to_b = defaultdict(list)
for bits, ss in b_to_s.items():
    s_to_b[ss].append(bits)

total_comps = 0
fiber_comp_data = {}
for ss in s_to_b:
    fiber = s_to_b[ss]
    fiber_set = set(fiber)
    remaining = set(fiber)
    n_comps = 0
    while remaining:
        start = remaining.pop()
        queue = deque([start])
        comp = {start}
        while queue:
            curr = queue.popleft()
            for bit in range(m):
                nbr = curr ^ (1 << bit)
                if nbr in remaining:
                    remaining.discard(nbr)
                    comp.add(nbr)
                    queue.append(nbr)
        n_comps += 1
    total_comps += n_comps
    fiber_comp_data[ss] = n_comps

component_counts.append((n, total_comps))
elapsed = time.time() - t0
print(f"  n={n}: fraction = {frac} = {float(frac):.6f}, "
      f"components = {total_comps}, time = {elapsed:.1f}s")

print()

# ==========================================================================
# PART 2: PATTERN ANALYSIS
# ==========================================================================

print()
print("PART 2: PATTERN ANALYSIS")
print("-" * 60)
print()

print("  WITHIN-FIBER FRACTIONS:")
for n_val, frac in fractions_list:
    print(f"    n={n_val}: {frac} = {frac.numerator}/{frac.denominator} "
          f"= {float(frac):.6f}")

print()

# Check patterns
nums = [f.numerator for _, f in fractions_list]
dens = [f.denominator for _, f in fractions_list]
print(f"  Numerators:   {nums}")
print(f"  Denominators: {dens}")
print()

# Check if denominators are powers of 2
for n_val, frac in fractions_list:
    d = frac.denominator
    is_pow2 = d & (d-1) == 0
    if is_pow2:
        exp = d.bit_length() - 1
        print(f"    n={n_val}: denom = {d} = 2^{exp}")
    else:
        print(f"    n={n_val}: denom = {d} (NOT a power of 2)")

print()

# Check various conjectures
print("  Checking conjectures:")
for n_val, frac in fractions_list:
    # Is it (2n-1)/2^n?
    conj1 = Fraction(2*n_val - 1, 2**n_val)
    # Is it 1/n?
    conj2 = Fraction(1, n_val)
    # Is it (n-1)/(2*C(n,2))?
    conj3 = Fraction(n_val - 1, 2 * comb(n_val, 2))
    # Is it something/C(n,2)?
    conj4 = frac * comb(n_val, 2)

    print(f"    n={n_val}: frac={float(frac):.6f}, "
          f"(2n-1)/2^n={float(conj1):.6f}, "
          f"1/n={float(conj2):.6f}, "
          f"(n-1)/(2*C(n,2))={float(conj3):.6f}, "
          f"frac*C(n,2)={float(conj4):.6f}")

print()

# COMPONENT COUNTS
print("  TOTAL FIBER COMPONENTS:")
for n_val, tc in component_counts:
    print(f"    n={n_val}: {tc}")

comp_vals = [tc for _, tc in component_counts]
print(f"\n  Values: {comp_vals}")
print(f"  Ratios: {[comp_vals[i+1]/comp_vals[i] for i in range(len(comp_vals)-1)]}")
print()

# ==========================================================================
# PART 3: FIBER-LEVEL DETAIL AT n=6
# ==========================================================================

print()
print("PART 3: FIBER-LEVEL DETAIL AT n=6")
print("-" * 60)
print()

n = 6
m = comb(n, 2)

# Compute H for all tournaments at n=6
print("  Computing H for all n=6 tournaments...")
t0 = time.time()
b_to_h = {}
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    b_to_h[bits] = h

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")
print()

# Score → H values
s_to_h = defaultdict(set)
s_to_hcount = defaultdict(lambda: Counter())
for bits in range(1 << m):
    ss = b_to_s[bits]
    h = b_to_h[bits]
    s_to_h[ss].add(h)
    s_to_hcount[ss][h] += 1

# Show fiber detail
print(f"  {'score':>25s}  {'size':>6s}  {'comps':>5s}  {'H values':>20s}  {'#H':>3s}")
print(f"  {'─'*25}  {'─'*6}  {'─'*5}  {'─'*20}  {'─'*3}")

seqs6 = sorted(s_to_b.keys())
for ss in seqs6:
    fsize = len(s_to_b[ss])
    n_comps = fiber_comp_data[ss]
    h_vals = sorted(s_to_h[ss])
    n_h = len(h_vals)

    h_str = str(h_vals) if len(h_vals) <= 5 else f"[{min(h_vals)}..{max(h_vals)}] ({len(h_vals)} vals)"
    print(f"  {str(ss):>25s}  {fsize:>6d}  {n_comps:>5d}  {h_str:>20s}  {n_h:>3d}")

print()

# Count H-ambiguous
ambiguous = sum(1 for ss in seqs6 if len(s_to_h[ss]) > 1)
print(f"  H-AMBIGUOUS score sequences at n=6: {ambiguous}/{len(seqs6)}")
print()

# ==========================================================================
# PART 4: THE H ≈ H_max - c*S₂ FORMULA AT n=6
# ==========================================================================

print()
print("PART 4: THE H ≈ H_max + b × S₂ FORMULA AT n=6")
print("-" * 60)
print()

n = 6
s_bar = (n - 1) / 2

data_s2 = []
data_h_mean = []
for ss in seqs6:
    s2 = sum(x**2 for x in ss) - n * s_bar**2
    hvals = sorted(s_to_h[ss])
    # Weight by count
    total = sum(s_to_hcount[ss].values())
    mean_h = sum(h * c for h, c in s_to_hcount[ss].items()) / total
    data_s2.append(s2)
    data_h_mean.append(mean_h)

data_s2 = np.array(data_s2)
data_h_mean = np.array(data_h_mean)

# Linear fit
A = np.column_stack([np.ones_like(data_s2), data_s2])
result = np.linalg.lstsq(A, data_h_mean, rcond=None)
a, b = result[0]

r2 = np.corrcoef(data_s2, data_h_mean)[0, 1] ** 2

h_max = max(data_h_mean)
print(f"  LINEAR FIT: H ≈ {a:.2f} + ({b:.3f}) × S₂  (R² = {r2:.4f})")
print(f"  Intercept ≈ {a:.0f}, slope = {b:.4f}")
print()

# Show fit quality
print(f"  {'score':>25s}  {'S₂':>5s}  {'mean H':>8s}  {'pred':>8s}  {'error':>6s}")
for i, ss in enumerate(seqs6):
    h_pred = a + b * data_s2[i]
    err = data_h_mean[i] - h_pred
    print(f"  {str(ss):>25s}  {data_s2[i]:>5.1f}  {data_h_mean[i]:>8.1f}  "
          f"{h_pred:>8.1f}  {err:>+6.1f}")

print()
print(f"  R² = {r2:.4f} at n=6 (vs 0.9792 at n=5, 1.0000 at n=3,4)")
print(f"  The linear approximation degrades slightly with n.")
print()

# ==========================================================================
# PART 5: OEIS CHECK FOR NEW SEQUENCES
# ==========================================================================

print()
print("PART 5: NEW SEQUENCES SUMMARY")
print("-" * 60)
print()

print("  SEQUENCE 1: Within-fiber edge fraction (as fractions)")
for n_val, frac in fractions_list:
    print(f"    n={n_val}: {frac}")

print()
print("  SEQUENCE 2: Total fiber components")
for n_val, tc in component_counts:
    print(f"    n={n_val}: {tc}")

print()

# More sequences from the fiber structure

# Sequence 3: Number of ISOLATED tournaments (in singletons under within-fiber flips)
isolated_counts = []
for n_val in [3, 4, 5]:
    m_val = comb(n_val, 2)
    b_to_s_temp = {}
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        ss = score_seq(adj, n_val)
        b_to_s_temp[bits] = ss

    s_to_b_temp = defaultdict(set)
    for bits, ss in b_to_s_temp.items():
        s_to_b_temp[ss].add(bits)

    isolated = 0
    for bits in range(1 << m_val):
        has_neighbor = False
        for bit in range(m_val):
            nbr = bits ^ (1 << bit)
            if b_to_s_temp[bits] == b_to_s_temp[nbr]:
                has_neighbor = True
                break
        if not has_neighbor:
            isolated += 1
    isolated_counts.append((n_val, isolated))

# n=6
isolated_n6 = 0
for bits in range(1 << comb(6, 2)):
    has_neighbor = False
    for bit in range(comb(6, 2)):
        nbr = bits ^ (1 << bit)
        if b_to_s[bits] == b_to_s[nbr]:
            has_neighbor = True
            break
    if not has_neighbor:
        isolated_n6 += 1
isolated_counts.append((6, isolated_n6))

print("  SEQUENCE 3: Isolated tournaments (no score-preserving flip)")
for n_val, ic in isolated_counts:
    print(f"    n={n_val}: {ic}")

print()

# Sequence 4: Number of connected fibers
connected_fiber_counts = []
for n_val in [3, 4, 5]:
    m_val = comb(n_val, 2)
    b_to_s_temp = {}
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        ss = score_seq(adj, n_val)
        b_to_s_temp[bits] = ss

    s_to_b_temp = defaultdict(list)
    for bits, ss in b_to_s_temp.items():
        s_to_b_temp[ss].append(bits)

    n_connected = 0
    for ss in s_to_b_temp:
        fiber = s_to_b_temp[ss]
        fiber_set = set(fiber)
        remaining = set(fiber)
        n_comps = 0
        while remaining:
            start = remaining.pop()
            queue = deque([start])
            comp = {start}
            while queue:
                curr = queue.popleft()
                for bit in range(m_val):
                    nbr = curr ^ (1 << bit)
                    if nbr in remaining:
                        remaining.discard(nbr)
                        comp.add(nbr)
                        queue.append(nbr)
            n_comps += 1
        if n_comps == 1:
            n_connected += 1

    connected_fiber_counts.append((n_val, n_connected))

# n=6
n_connected_6 = sum(1 for ss in seqs6 if fiber_comp_data[ss] == 1)
connected_fiber_counts.append((6, n_connected_6))

print("  SEQUENCE 4: Number of connected fibers")
for n_val, nc in connected_fiber_counts:
    total = len([1 for _ in (range(3, n_val+1))]) if n_val <= 3 else '?'
    print(f"    n={n_val}: {nc} connected fibers")

print()

# Sequence 5: Max fiber component count
max_comp_counts = []
for n_val in [3, 4, 5]:
    m_val = comb(n_val, 2)
    b_to_s_temp = {}
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        ss = score_seq(adj, n_val)
        b_to_s_temp[bits] = ss

    s_to_b_temp = defaultdict(list)
    for bits, ss in b_to_s_temp.items():
        s_to_b_temp[ss].append(bits)

    max_comps = 0
    for ss in s_to_b_temp:
        fiber = s_to_b_temp[ss]
        fiber_set = set(fiber)
        remaining = set(fiber)
        n_comps = 0
        while remaining:
            start = remaining.pop()
            queue = deque([start])
            comp = {start}
            while queue:
                curr = queue.popleft()
                for bit in range(m_val):
                    nbr = curr ^ (1 << bit)
                    if nbr in remaining:
                        remaining.discard(nbr)
                        comp.add(nbr)
                        queue.append(nbr)
            n_comps += 1
        max_comps = max(max_comps, n_comps)

    max_comp_counts.append((n_val, max_comps))

max_comps_6 = max(fiber_comp_data.values())
max_comp_counts.append((6, max_comps_6))

print("  SEQUENCE 5: Max fiber component count")
for n_val, mc in max_comp_counts:
    print(f"    n={n_val}: {mc}")

print()

# ==========================================================================
# PART 6: WITHIN-FIBER FRACTION DEEPER ANALYSIS
# ==========================================================================

print()
print("PART 6: WITHIN-FIBER FRACTION — ANALYTICAL APPROACH")
print("-" * 60)
print()

print("""  The within-fiber fraction = P(a random flip preserves scores).

  For a uniformly random tournament T on n vertices:
  - Pick a random arc (i,j) with i→j.
  - score[i] = s_i, score[j] = s_j.
  - The flip preserves the sorted score iff |s_i - s_j| ≤ 1.

  So the fraction = E[# arcs connecting vertices with |score diff| ≤ 1] / C(n,2).

  At n=3: every tournament has all scores within ±1 of each other
  (scores are either 0,1,2 or 1,1,1). Fraction = 1/2.

  The fraction EQUALS:
    (1/2^C(n,2)) × Σ_T Σ_{(i,j)} 1[|s_i(T) - s_j(T)| ≤ 1] / C(n,2)
  where the sum is over all labeled tournaments T.
""")

# Compute the fraction differently: for each tournament,
# count what fraction of its arcs connect adjacent-score vertices
for n_val in [3, 4, 5]:
    m_val = comb(n_val, 2)
    total_adj_arcs = 0
    total_arcs = 0

    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        scores_lab = score_seq_labeled(adj, n_val)

        for i in range(n_val):
            for j in range(i+1, n_val):
                total_arcs += 1
                if abs(scores_lab[i] - scores_lab[j]) <= 1:
                    total_adj_arcs += 1

    frac = Fraction(total_adj_arcs, total_arcs)
    print(f"  n={n_val}: P(|score_i - score_j| ≤ 1) = {frac} = {float(frac):.6f}")

print()
print("  NOTE: This is the probability that a random pair of vertices")
print("  in a random tournament has adjacent scores.")
print("  The within-fiber fraction is DIFFERENT (it also requires")
print("  the arc to go from the higher-scored to the lower-scored vertex).")
print()

# ==========================================================================
# PART 7: THE COEFFICIENT b(n) ACROSS n
# ==========================================================================

print()
print("PART 7: THE SLOPE b(n) = (1 - H_max(n)) / S₂_max(n)")
print("-" * 60)
print()

# H_max values (maximum # Hamiltonian paths)
h_max_vals = {3: 3, 4: 5, 5: 15, 6: 81}  # n=6: H_max for PoS(6)

for n_val in [3, 4, 5, 6]:
    s2_max = n_val * (n_val**2 - 1) / 12
    h_max = h_max_vals.get(n_val, '?')
    if isinstance(h_max, int):
        b_exact = (1 - h_max) / s2_max
        print(f"  n={n_val}: H_max = {h_max}, S₂_max = {s2_max:.1f}, "
              f"b = (1-{h_max})/{s2_max:.1f} = {b_exact:.4f}")
    else:
        print(f"  n={n_val}: H_max = {h_max}, S₂_max = {s2_max:.1f}")

print()

# What's H_max at n=6? Let's find it from our data
h_max_6 = max(b_to_h.values())
print(f"  Actual H_max at n=6: {h_max_6}")
print(f"  S₂_max at n=6: {6 * (36 - 1) / 12:.1f}")
b_6 = (1 - h_max_6) / (6 * 35 / 12)
print(f"  b(6) = (1 - {h_max_6}) / {6*35/12:.1f} = {b_6:.4f}")
print()

# H_max sequence: 1, 1, 3, 5, 15, ?
# n=1: 1, n=2: 1, n=3: 3, n=4: 5, n=5: 15, n=6: ?
print(f"  H_max sequence: 1, 1, 3, 5, 15, {h_max_6}")
print(f"  Ratios: {3/1:.1f}, {5/3:.2f}, {15/5:.1f}, {h_max_6/15:.2f}")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: NEW SEQUENCES")
print("=" * 72)
print()

print("  SEQUENCE A: Within-fiber edge fraction")
for n_val, frac in fractions_list:
    print(f"    n={n_val}: {frac}")

print()
print("  SEQUENCE B: Total fiber components")
for n_val, tc in component_counts:
    print(f"    n={n_val}: {tc}")

print()
print("  SEQUENCE C: Isolated tournaments")
for n_val, ic in isolated_counts:
    print(f"    n={n_val}: {ic}")

print()
print("  SEQUENCE D: Connected fibers")
for n_val, nc in connected_fiber_counts:
    print(f"    n={n_val}: {nc}")

print()
print("  SEQUENCE E: Max fiber components")
for n_val, mc in max_comp_counts:
    print(f"    n={n_val}: {mc}")

print()
print("  SEQUENCE F: H_max(n)")
print(f"    3, 5, 15, {h_max_6}")
print(f"    (max # Hamiltonian paths in any n-tournament)")

print()
print("  H-FORMULA QUALITY:")
print(f"    n=3: R² = 1.0000 (exact)")
print(f"    n=4: R² = 1.0000 (exact)")
print(f"    n=5: R² = 0.9792")
print(f"    n=6: R² = {r2:.4f}")

print()
print("opus-2026-03-22-S192")
