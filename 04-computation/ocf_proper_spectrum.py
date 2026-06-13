#!/usr/bin/env python3
"""
ocf_proper_spectrum.py — opus-2026-03-14-S71d

Correct I(Ω, x) computation using DIRECTED odd cycles (not vertex sets).
MISTAKE-023: vertex-set counting undercounts α₁ when multiple directed
cycles exist on the same vertex set.

Key questions:
1. Re-verify the 3/2 ratio with correct α values
2. Re-verify Vandermonde extraction
3. What does α₁(directed) vs α₁(vertex-set) look like?
4. k-nacci convergence analysis (no cycle issues here)
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

def enumerate_directed_odd_cycles(A, n, max_length=None):
    """Enumerate all directed odd cycles as (normalized tuple, vertex frozenset) pairs.

    Normalization: start at minimum vertex, choose direction so neighbor
    of min vertex with smaller index comes second (if tied, arbitrary).
    Each physical directed cycle is counted exactly once.
    """
    if max_length is None:
        max_length = n if n % 2 == 1 else n - 1

    cycles = []
    for length in range(3, max_length + 1, 2):
        if length > n:
            break
        for combo in combinations(range(n), length):
            verts = list(combo)
            v0 = verts[0]  # minimum vertex (combo is sorted)
            # Find ALL directed cycles on this vertex set
            for perm in permutations(verts[1:]):
                order = (v0,) + perm
                if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                    # This is a valid directed cycle
                    # Normalize: start at v0 (already done), choose direction
                    # so that the cycle is lexicographically smaller
                    # Forward: order
                    # Reverse would be: v0, perm[-1], perm[-2], ..., perm[0]
                    reverse = (v0,) + perm[::-1]
                    # Keep the lexicographically smaller one
                    canon = min(order, reverse)
                    cycles.append((canon, frozenset(combo)))

    # Remove duplicates (same canonical form)
    seen = set()
    unique = []
    for canon, vset in cycles:
        if canon not in seen:
            seen.add(canon)
            unique.append((canon, vset))

    return unique

def compute_I_omega(A, n, x_vals, max_length=None):
    """Compute I(Ω(T), x) using proper directed cycle enumeration."""
    cycles = enumerate_directed_odd_cycles(A, n, max_length)
    nc = len(cycles)

    # Build conflict graph: two cycles conflict iff vertex sets intersect
    vsets = [c[1] for c in cycles]
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                adj[i][j] = True
                adj[j][i] = True

    # Count independent sets by size
    alpha = [0] * (nc + 1)
    for mask in range(1 << nc):
        bits_list = [i for i in range(nc) if mask & (1 << i)]
        indep = True
        for a in range(len(bits_list)):
            for b in range(a+1, len(bits_list)):
                if adj[bits_list[a]][bits_list[b]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(bits_list)] += 1

    results = {}
    for x in x_vals:
        val = sum(alpha[k] * x**k for k in range(nc + 1))
        results[x] = val

    return results, alpha, nc

# ======================================================================
# PART 1: Verify I(Ω, 2) = H at n=5 (exhaustive)
# ======================================================================
print("=" * 70)
print("PART 1: VERIFY I(Ω, 2) = H WITH DIRECTED CYCLES (n=5 exhaustive)")
print("=" * 70)

n = 5
tb = n*(n-1)//2
mismatches = 0
alpha_diff = 0  # cases where directed ≠ vertex-set count

t0 = time.time()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    ix, alpha, nc = compute_I_omega(A, n, [2])
    if ix[2] != H:
        mismatches += 1
        if mismatches <= 3:
            print(f"  MISMATCH bits={bits}: H={H}, I(2)={ix[2]}")

dt = time.time() - t0
print(f"  n=5: mismatches = {mismatches}/{1 << tb}, time = {dt:.1f}s")

# ======================================================================
# PART 2: Corrected Vandermonde at n=5 (exhaustive)
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: CORRECTED VANDERMONDE EXTRACTION (n=5)")
print("=" * 70)

n = 5
tb = n*(n-1)//2
data5 = defaultdict(list)

for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    ix, alpha, nc = compute_I_omega(A, n, [2, 3, 6])
    a1 = alpha[1] if len(alpha) > 1 else 0
    a2 = alpha[2] if len(alpha) > 2 else 0
    data5[(a1, a2)].append({'H': H, 'I2': ix[2], 'I3': ix[3], 'I6': ix[6], 'bits': bits})

print(f"\n  {'(α₁,α₂)':>10s} | {'cnt':>4s} | {'H':>5s} {'I(2)':>5s} {'I(3)':>5s} {'I(6)':>5s} | {'H=I(2)':>6s}")
print(f"  {'-'*10}-+-{'-'*4}-+-{'-'*23}-+-{'-'*6}")
for key in sorted(data5.keys()):
    group = data5[key]
    ex = group[0]
    match = "✓" if ex['H'] == ex['I2'] else "✗"
    # Verify I(2) = 1 + 2a1 + 4a2
    a1, a2 = key
    expected = 1 + 2*a1 + 4*a2
    print(f"  ({a1:2d},{a2:2d})    | {len(group):4d} | {ex['H']:5d} {ex['I2']:5d} {ex['I3']:5d} {ex['I6']:5d} | {match:>6s}  (1+2·{a1}+4·{a2}={expected})")

# Vandermonde: from I(2) and I(3), recover α₁ and α₂
print(f"\n  Vandermonde extraction: I(2)=1+2α₁+4α₂, I(3)=1+3α₁+9α₂")
print(f"  α₂ = (2·I(3) - 3·I(2) + 1) / 6")
print(f"  α₁ = (I(2) - 1 - 4α₂) / 2")
vandermonde_ok = 0
vandermonde_fail = 0
for key, group in data5.items():
    a1, a2 = key
    ex = group[0]
    a2_est = (2*ex['I3'] - 3*ex['I2'] + 1) / 6
    a1_est = (ex['I2'] - 1 - 4*a2_est) / 2
    if abs(a1_est - a1) < 0.01 and abs(a2_est - a2) < 0.01:
        vandermonde_ok += len(group)
    else:
        vandermonde_fail += len(group)
        print(f"  FAIL: α=({a1},{a2}), recovered=({a1_est:.1f},{a2_est:.1f})")
print(f"  Vandermonde: {vandermonde_ok}/{vandermonde_ok+vandermonde_fail} correct")

# ======================================================================
# PART 3: n=7 directed cycle analysis (sample)
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: n=7 DIRECTED CYCLE SPECTRUM (100 samples)")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

data7 = []
t0 = time.time()
for trial in range(100):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Count directed cycles by length
    cycles = enumerate_directed_odd_cycles(A, n)
    c3 = sum(1 for c, vs in cycles if len(c) == 3)
    c5 = sum(1 for c, vs in cycles if len(c) == 5)
    c7 = sum(1 for c, vs in cycles if len(c) == 7)
    c5_vsets = len(set(vs for c, vs in cycles if len(c) == 5))
    c7_vsets = len(set(vs for c, vs in cycles if len(c) == 7))

    # Compute I(Ω, x)
    ix, alpha, nc = compute_I_omega(A, n, [2, 3, 6])

    data7.append({
        'H': H, 'I2': ix[2], 'I3': ix[3], 'I6': ix[6],
        'c3': c3, 'c5': c5, 'c7': c7,
        'c5_vsets': c5_vsets, 'c7_vsets': c7_vsets,
        'a1': alpha[1] if len(alpha) > 1 else 0,
        'a2': alpha[2] if len(alpha) > 2 else 0,
        'a3': alpha[3] if len(alpha) > 3 else 0,
        'nc': nc
    })

    if trial % 20 == 0:
        dt = time.time() - t0
        match = "✓" if ix[2] == H else "✗"
        print(f"  trial {trial}: {dt:.1f}s, H={H}, I(2)={ix[2]} {match}, α₁={alpha[1] if len(alpha)>1 else 0}, c3={c3}, c5={c5}({c5_vsets}vs), c7={c7}({c7_vsets}vs)")

dt = time.time() - t0
print(f"  Done: {dt:.1f}s")

# Verify I(2) = H
mismatches = sum(1 for d in data7 if d['I2'] != d['H'])
print(f"\n  I(Ω,2) = H? Mismatches: {mismatches}/{len(data7)}")

# Show the directed vs vertex-set discrepancy
print(f"\n  Directed vs vertex-set cycle counts:")
print(f"  5-cycles: directed/vset ratio: mean={np.mean([d['c5']/d['c5_vsets'] if d['c5_vsets']>0 else 1 for d in data7]):.2f}")
print(f"  7-cycles: directed/vset ratio: mean={np.mean([d['c7']/d['c7_vsets'] if d['c7_vsets']>0 else 1 for d in data7]):.2f}")

# ======================================================================
# PART 4: Re-check Vandermonde and 3/2 ratio at n=7
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: VANDERMONDE AND 3/2 RATIO AT n=7 (corrected)")
print("=" * 70)

# Vandermonde extraction
vand_ok = 0
for d in data7:
    a2_est = (2*d['I3'] - 3*d['I2'] + 1) / 6
    a1_est = (d['I2'] - 1 - 4*a2_est) / 2
    if abs(a1_est - d['a1']) < 0.01 and abs(a2_est - d['a2']) < 0.01:
        vand_ok += 1

print(f"  Vandermonde (I(2),I(3)) → (α₁,α₂): {vand_ok}/{len(data7)} correct")

# I₃/I₂ ratio
ratios = [d['I3']/d['I2'] for d in data7 if d['I2'] > 0]
print(f"  I(3)/I(2) ratio: mean={np.mean(ratios):.6f}, std={np.std(ratios):.6f}")
print(f"  3/2 = {1.5}")

# Does α₁(directed) differ meaningfully from c3+c5_vsets+c7_vsets?
print(f"\n  α₁(directed) vs c3+c5_vsets+c7_vsets:")
for d in data7[:5]:
    vset_total = d['c3'] + d['c5_vsets'] + d['c7_vsets']
    dir_total = d['a1']
    extra = dir_total - vset_total
    print(f"    α₁={dir_total}, vset={vset_total} (c3={d['c3']}, c5vs={d['c5_vsets']}, c7vs={d['c7_vsets']}), extra directed={extra}")

# ======================================================================
# PART 5: k-nacci convergence (no cycle issues)
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: k-NACCI CONVERGENCE RATES AND CONSTANTS")
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

print(f"\n  Convergence constants:")
print(f"  φ_k - 2 ≈ C₂ · (1/2)^k")
print(f"  ψ_k - 3 ≈ C₃ · (2/3)^k")
print()
print(f"  {'k':>3s}  {'C₂ = (φ_k-2)·2^k':>18s}  {'C₃ = (ψ_k-3)·(3/2)^k':>22s}  {'ratio C₃/C₂':>12s}")
print(f"  {'─'*3}  {'─'*18}  {'─'*22}  {'─'*12}")

for k in range(3, 25):
    phi_k = knacci_root(k)
    psi_k = weighted_knacci_root(k)
    if phi_k is None or psi_k is None:
        continue
    c2 = (phi_k - 2) * 2**k
    c3 = (psi_k - 3) * (3/2)**k
    ratio = c3/c2 if abs(c2) > 1e-15 else float('nan')
    mark = " ←" if k in [7, 8] else ""
    print(f"  {k:3d}  {c2:18.10f}  {c3:22.10f}  {ratio:12.6f}{mark}")

# What are the limiting constants?
phi_30 = knacci_root(30)
psi_30 = weighted_knacci_root(30)
if phi_30 and psi_30:
    C2_limit = (phi_30 - 2) * 2**30
    C3_limit = (psi_30 - 3) * (3/2)**30
    print(f"\n  Limiting constants (from k=30):")
    print(f"    C₂ ≈ {C2_limit:.10f}")
    print(f"    C₃ ≈ {C3_limit:.10f}")
    print(f"    C₃/C₂ ≈ {C3_limit/C2_limit:.10f}")

# ======================================================================
# PART 6: SYNTHESIS — HOW 7 AND 8 CHANGE THE 2-3 FRAMEWORK
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: SYNTHESIS — THE 7→8 TRANSITION IN THE 2-3 FRAMEWORK")
print("=" * 70)

print("""
CORRECTED PICTURE (after MISTAKE-023):

The conflict graph Ω(T) uses DIRECTED odd cycles as vertices.
A 5-vertex set can have up to 12 directed 5-cycles.
A 7-vertex set can have up to 360 directed 7-cycles.

So α₁ (total number of vertices of Ω) is MUCH larger than just
counting vertex-sets with at least one cycle.

This means the independence polynomial I(Ω, x) has potentially
very high degree, and the Vandermonde framework I(2), I(3), I(6)
extracts α₁ and α₂ from a RICH structure.

AT n=7:
  - 3-cycle vertex-sets: up to C(7,3)=35 (each supports ≤1 directed cycle)
  - 5-cycle vertex-sets: up to C(7,5)=21 (each supports up to 12 directed cycles)
  - 7-cycle vertex-sets: 1 (supports up to 360 directed cycles)
  - α₁ can be very large (up to 35+252+360=647 theoretical max)
  - α₂ includes disjoint pairs of directed cycles

AT n=8:
  - Same as n=7 but with C(8,k) more vertex-sets
  - PLUS (3,5) cross-level pairs: disjoint 3-cycle and 5-cycle (directed!)
  - 8-vertex set: up to (8-1)!/2 = 2520 directed 8-cycles (EVEN, not counted)
  - Cross-level α₂^{35}: pairs of directed cycles from different levels

THE 2-3 KEYS:
  I(Ω, 2) = H (OCF, proved)
  I(Ω, 3) gives a related but richer count
  Together: Vandermonde extraction of α₁, α₂ (at n≤7 where α₃=0)
  At n≥9: need I(Ω, 6) or higher for α₃

  The k-nacci root approaches 2 at rate 1/2 per step (binary refining)
  The weighted root approaches 3 at rate 2/3 per step (ternary refining)
  These are the CONVERGENCE RATES of the recurrence towers.
""")

print("Done.")
