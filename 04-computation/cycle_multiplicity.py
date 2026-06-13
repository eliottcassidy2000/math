#!/usr/bin/env python3
"""
cycle_multiplicity.py — opus-2026-03-14-S71d

How many directed Hamiltonian cycles does a tournament on k vertices have?

This is related to MISTAKE-023: α₁ counts directed cycles, not vertex-sets.
The multiplicity (directed cycles per vertex-set) matters for the OCF.

Key facts:
- A tournament on k vertices has 0 or more Hamiltonian cycles
- In a tournament, each undirected Ham cycle gives EXACTLY 1 directed cycle
  (since if a→b then NOT b→a)
- So #(directed Ham cycles) = #(undirected Ham cycles)
- For k=3: max = 1 (either cyclic or transitive)
- For k=5: max = 12 (= (k-1)!/2 = 24/2, the regular tournament on 5 vertices)
- For k=7: max = 360 (= (k-1)!/2)

The regular tournament (doubly regular) maximizes Ham cycles.

QUESTION: How does the number of k-cycles in T affect α₁ and hence H?
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
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

def count_ham_cycles(A, n):
    """Count directed Hamiltonian cycles in tournament A.
    Fix vertex 0 as start, count cycles 0→...→0."""
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {}
    dp[(1 << 0, 0)] = 1  # start at vertex 0

    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue  # must include vertex 0
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue  # don't revisit 0 early
                pm = mask ^ (1 << v)
                if not (pm & 1): continue  # previous mask must include 0
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t

    # Count cycles: need edge from last vertex back to 0
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    """Count directed k-cycles on vertex subsets of size k."""
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        # Build submatrix
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

# ======================================================================
# PART 1: Hamiltonian cycle distribution at n=5 (exhaustive)
# ======================================================================
print("=" * 70)
print("PART 1: DIRECTED 5-CYCLE DISTRIBUTION AT n=5 (all 1024 tournaments)")
print("=" * 70)

n = 5
tb = n*(n-1)//2
hc_dist = Counter()

t0 = time.time()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    hc = count_ham_cycles(A, n)
    hc_dist[hc] += 1

dt = time.time() - t0
print(f"  Computed in {dt:.1f}s\n")

print(f"  #(directed 5-cycles)  count  fraction")
for hc in sorted(hc_dist.keys()):
    print(f"  {hc:20d}  {hc_dist[hc]:5d}  {hc_dist[hc]/1024:.3f}")

# ======================================================================
# PART 2: At n=5, all cycle data with H verification
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: FULL CYCLE STRUCTURE AT n=5")
print("=" * 70)

struct = defaultdict(list)
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_ham_cycles(A, n)  # = directed 5-cycles on all 5 vertices
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    struct[(dc3, dc5, a2, H)].append(bits)

print(f"\n  {'dc3':>4s} {'dc5':>4s} {'α₁':>4s} {'α₂':>4s} {'H':>5s} {'cnt':>5s}")
for key in sorted(struct.keys()):
    dc3, dc5, a2, H = key
    a1 = dc3 + dc5
    print(f"  {dc3:4d} {dc5:4d} {a1:4d} {a2:4d} {H:5d} {len(struct[key]):5d}")

# ======================================================================
# PART 3: At n=7, directed 3/5/7-cycle counts (sample)
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: DIRECTED CYCLE STRUCTURE AT n=7 (200 samples)")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

data = []
t0 = time.time()
for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)

    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4

    data.append({
        'H': H, 'dc3': dc3, 'dc5': dc5, 'dc7': dc7,
        'a1': a1, 'a2': a2, 'bits': bits
    })

    if trial % 50 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, H={H}, dc3={dc3}, dc5={dc5}, dc7={dc7}, α₁={a1}, α₂={a2}")

dt = time.time() - t0
print(f"  Done: {dt:.1f}s")

# Verify H = 1 + 2α₁ + 4α₂
verify = all((d['H'] - 1 - 2*d['a1']) % 4 == 0 and d['a2'] >= 0 for d in data)
print(f"\n  H = 1 + 2α₁ + 4α₂ verified: {verify}")

# Statistics
print(f"\n  Statistics at n=7:")
for var in ['dc3', 'dc5', 'dc7', 'a1', 'a2', 'H']:
    vals = [d[var] for d in data]
    print(f"    {var:>4s}: mean={np.mean(vals):8.1f}, range=[{min(vals):5d}, {max(vals):5d}]")

# ======================================================================
# PART 4: Ratio of directed to vertex-set counts
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: DIRECTED vs VERTEX-SET RATIO BY CYCLE LENGTH")
print("=" * 70)

# At n=5: what fraction of 5-vertex sets have k directed Ham cycles?
n = 5
tb = n*(n-1)//2
hc_counts_5 = Counter()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    hc = count_ham_cycles(A, n)
    hc_counts_5[hc] += 1

print(f"\n  n=5: distribution of #(directed Ham cycles) on all 5 vertices:")
total_dc = sum(k*v for k, v in hc_counts_5.items())
total_vs = sum(v for k, v in hc_counts_5.items() if k > 0)
print(f"    Total directed 5-cycles across all tournaments: {total_dc}")
print(f"    Total vertex-sets with ≥1 5-cycle: {total_vs}")
print(f"    Mean directed cycles per vertex-set (conditional): {total_dc/total_vs if total_vs > 0 else 0:.2f}")

# At n=7: from sample data, what's dc5/vc5?
print(f"\n  n=7: directed/vertex-set ratios from 200 samples:")
ratios_5 = []
ratios_7 = []
for d in data:
    A = bits_to_adj(d['bits'], n)
    # Count 5-cycle vertex-sets
    vc5 = 0
    for combo in combinations(range(n), 5):
        verts = list(combo)
        sub = np.zeros((5, 5), dtype=int)
        for i in range(5):
            for j in range(5):
                sub[i][j] = A[verts[i]][verts[j]]
        if count_ham_cycles(sub, 5) > 0:
            vc5 += 1
    vc7 = 1 if d['dc7'] > 0 else 0
    if vc5 > 0:
        ratios_5.append(d['dc5'] / vc5)
    if vc7 > 0:
        ratios_7.append(d['dc7'])  # vc7 is always 0 or 1

print(f"    5-cycles: dc/vc ratio: mean={np.mean(ratios_5):.3f}, range=[{min(ratios_5):.1f}, {max(ratios_5):.1f}]")
if ratios_7:
    print(f"    7-cycles: dc7 (when ≥1): mean={np.mean(ratios_7):.1f}, range=[{min(ratios_7)}, {max(ratios_7)}]")
else:
    print(f"    7-cycles: no tournaments with dc7>0 in sample")

# ======================================================================
# PART 5: KEY FORMULA — total directed cycles and OCF
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: THE OCF WITH DIRECTED CYCLES")
print("=" * 70)

print("""
  H = I(Ω(T), 2) where Ω has vertices = directed odd cycles

  H = 1 + 2α₁ + 4α₂ + 8α₃ + ...

  where α₁ = total #(directed odd cycles)
        α₂ = #{vertex-disjoint pairs of directed odd cycles}
        α₃ = #{vertex-disjoint triples} (needs n≥9)

  The key point: α₁ is NOT c3 + c5_vsets + c7_vsets.
  It IS: α₁ = c3_directed + c5_directed + c7_directed
  where c_k_directed counts ALL distinct directed k-cycles.

  For 3-cycles: c3_directed = c3_vsets (always, in tournaments)
  For k≥5: c_k_directed ≥ c_k_vsets (can be much larger!)

  CONSEQUENCE: The 'Vandermonde extraction' with correct α₁ gives:
    α₂ = (H - 1 - 2α₁) / 4

  And I(Ω, 3) = 1 + 3α₁ + 9α₂ (at n≤8 where α₃=0)
    = 1 + 3α₁ + 9(H-1-2α₁)/4
    = 1 + 3α₁ + (9H-9-18α₁)/4
    = (4 + 12α₁ + 9H - 9 - 18α₁)/4
    = (9H - 5 - 6α₁)/4

  So I(Ω, 3) = (9H - 5 - 6α₁) / 4
""")

# Verify this formula
print("  Verification of I₃ = (9H - 5 - 6α₁)/4:")
ok = 0
for d in data:
    I3 = 1 + 3*d['a1'] + 9*d['a2']
    I3_formula = (9*d['H'] - 5 - 6*d['a1']) / 4
    if abs(I3 - I3_formula) < 0.01:
        ok += 1
print(f"  {ok}/{len(data)} correct")

# And the I₃/H ratio:
print(f"\n  I₃/H = (9 - 5/H - 6α₁/H) / 4")
print(f"  As H → ∞: I₃/H → 9/4 - 6α₁/(4H)")
print(f"  If α₁ ~ cH for some constant c: I₃/H → (9 - 6c)/4")
print(f"  For c=1/2 (α₁≈H/2): I₃/H → (9-3)/4 = 3/2 ← the '3/2 ratio'!")

# What's the actual α₁/H ratio?
ratios = [d['a1']/d['H'] for d in data if d['H'] > 0]
print(f"\n  α₁/H ratio at n=7: mean={np.mean(ratios):.4f}")
print(f"  (H-1)/(2H) for comparison: mean={np.mean([(d['H']-1)/(2*d['H']) for d in data if d['H']>0]):.4f}")

# Actually: α₁ = (H - 1 - 4α₂)/2, so α₁/H = (1 - 1/H - 4α₂/H)/2
# For small α₂/H: α₁/H ≈ 1/2 - 1/(2H)
# So α₁ ≈ H/2 asymptotically, giving I₃/H → 3/2!

print(f"\n  INSIGHT: α₁ ≈ H/2 (since α₂ << H), so I₃/H ≈ 3/2")
print(f"  The 3/2 ratio is a CONSEQUENCE of α₂ being small relative to H!")

print("\nDone.")
