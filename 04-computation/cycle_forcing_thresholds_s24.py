#!/usr/bin/env python3
"""
cycle_forcing_thresholds_s24.py — opus-2026-04-05-S24

PRECISE QUESTION: For a tournament on n vertices with exactly c_p directed
p-cycles, what is the MINIMUM number of q-cycles guaranteed to exist?

Equivalently: what is the forcing threshold — the smallest c_p value at
which c_q > 0 is guaranteed for ALL tournaments?

Compute this for all pairs (p, q) with p < q, both odd, across n = 5..9.

KEY METHODOLOGICAL POINT: We must be careful about what holds at each n
separately. A pattern at n=5,6 may BREAK at n=7. So we compute each n
independently and only claim patterns when verified across multiple n.

CYCLE COUNTING:
- c₃: exact via score formula c₃ = C(n,3) - Σ C(s_v, 2)
- c₅: exact via DP on all C(n,5) five-subsets: c₅(sub) = (H(sub)-1)/2 - c₃(sub)
- c₇: exact via DP on all C(n,7) seven-subsets (need α₂ correction)
- c₉: direct enumeration on 9-subsets (only at n=9)
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import time
import sys
import math

sys.stdout.reconfigure(encoding='utf-8')

# ================================================================
# Fast tournament primitives
# ================================================================

def bits_to_adj(bits, n):
    """Convert bit-packed tournament to adjacency list-of-lists."""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament_bits(n):
    m = n * (n - 1) // 2
    return random.randint(0, (1 << m) - 1)

def score_c3(A, n):
    """c₃ = C(n,3) - Σ C(s_v, 2). O(n) time."""
    total = n * (n-1) * (n-2) // 6
    for v in range(n):
        s = sum(A[v])
        total -= s * (s-1) // 2
    return total

def dp_ham_paths(A, n):
    """Count Hamiltonian paths via bitmask DP. Works for small n."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))

def extract_sub(A, verts):
    """Extract subtournament on given vertices."""
    k = len(verts)
    S = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            S[i][j] = A[verts[i]][verts[j]]
    return S

def c5_via_subtournaments(A, n):
    """Count directed 5-cycles by summing over all 5-subsets.

    On 5 vertices: H = 1 + 2*(c₃ + c₅) because α₂ = 0 (two disjoint
    3-cycles need 6 vertices). So c₅(sub) = (H-1)/2 - c₃(sub).
    Each 5-cycle lives in exactly one 5-subset, so total c₅ = sum over subsets.
    """
    total_c5 = 0
    for verts in combinations(range(n), 5):
        S = extract_sub(A, verts)
        H = dp_ham_paths(S, 5)
        c3_sub = score_c3(S, 5)
        c5_sub = (H - 1) // 2 - c3_sub
        total_c5 += c5_sub
    return total_c5

def c7_via_subtournaments(A, n):
    """Count directed 7-cycles by summing over all 7-subsets.

    On 7 vertices: H = 1 + 2*(c₃+c₅+c₇) + 4*α₂
    where α₂ = number of pairs of disjoint 3-cycles (6 vertices, possible).
    A 3+5 pair needs 8 > 7 vertices, so α₂ = disjoint 3-cycle pairs only.

    c₇(sub) = (H - 1 - 4*α₂)/2 - c₃(sub) - c₅(sub)
    """
    total_c7 = 0
    for verts in combinations(range(n), 7):
        S = extract_sub(A, verts)
        H = dp_ham_paths(S, 7)
        c3_sub = score_c3(S, 7)
        c5_sub = c5_via_subtournaments(S, 7)

        # Count disjoint 3-cycle pairs in S
        # A pair is two vertex-disjoint 3-cycles
        # Find all 3-cycles, then count disjoint pairs
        cycles_3 = []
        for triple in combinations(range(7), 3):
            a, b, c = triple
            if S[a][b] and S[b][c] and S[c][a]:
                cycles_3.append(frozenset(triple))
            elif S[a][c] and S[c][b] and S[b][a]:
                cycles_3.append(frozenset(triple))

        alpha2 = 0
        for i in range(len(cycles_3)):
            for j in range(i+1, len(cycles_3)):
                if len(cycles_3[i] & cycles_3[j]) == 0:
                    alpha2 += 1

        c7_sub = (H - 1 - 4 * alpha2) // 2 - c3_sub - c5_sub
        total_c7 += c7_sub
    return total_c7

def count_kcycles_direct(A, n, k):
    """Direct count of directed k-cycles. Fix smallest vertex, try (k-1)! orders."""
    count = 0
    for verts in combinations(range(n), k):
        start = verts[0]
        rest = list(verts[1:])
        for perm in permutations(rest):
            cycle = (start,) + perm
            valid = True
            for i in range(k):
                if A[cycle[i]][cycle[(i+1) % k]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
    return count

# ================================================================
# Main computation
# ================================================================

def compute_all_cycle_counts(A, n, max_k=None):
    """Compute c₃, c₅, c₇, ... for tournament A on n vertices."""
    if max_k is None:
        max_k = n if n % 2 == 1 else n - 1

    result = {}
    result[3] = score_c3(A, n)

    if max_k >= 5 and n >= 5:
        result[5] = c5_via_subtournaments(A, n)

    if max_k >= 7 and n >= 7:
        if n == 7:
            result[7] = c7_via_subtournaments(A, n)
        elif n <= 9:
            result[7] = count_kcycles_direct(A, n, 7)

    if max_k >= 9 and n >= 9:
        result[9] = count_kcycles_direct(A, n, 9)

    return result

def analyze_forcing(n, exhaustive=True, num_samples=50000):
    """For tournament size n, compute the full (c_p, c_q) joint distribution."""
    m = n * (n - 1) // 2
    total = 1 << m
    max_k = n if n % 2 == 1 else n - 1
    odd_ks = [k for k in range(3, max_k + 1, 2)]

    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m}, cycle lengths: {odd_ks}")
    if exhaustive:
        print(f"EXHAUSTIVE: {total} tournaments")
    else:
        print(f"SAMPLED: {num_samples} tournaments")
    print(f"{'='*70}")

    # For each pair (p, q), track: c_p -> {c_q: count}
    pair_data = {}
    for p in odd_ks:
        for q in odd_ks:
            if q > p:
                pair_data[(p, q)] = defaultdict(lambda: defaultdict(int))

    # Also track individual distributions
    cycle_dist = {k: Counter() for k in odd_ks}

    # Main loop
    t0 = time.time()
    count = 0

    if exhaustive:
        iterator = range(total)
    else:
        iterator = (random_tournament_bits(n) for _ in range(num_samples))

    for bits in iterator:
        A = bits_to_adj(bits, n)
        cc = compute_all_cycle_counts(A, n, max_k)

        for k in odd_ks:
            if k in cc:
                cycle_dist[k][cc[k]] += 1

        for (p, q), data in pair_data.items():
            if p in cc and q in cc:
                data[cc[p]][cc[q]] += 1

        count += 1
        if count % 10000 == 0:
            elapsed = time.time() - t0
            rate = count / elapsed
            if exhaustive:
                eta = (total - count) / rate
                print(f"  {count}/{total} ({count/total*100:.1f}%), {rate:.0f}/s, ETA {eta:.0f}s", end='\r')
            else:
                print(f"  {count}/{num_samples}, {rate:.0f}/s", end='\r')

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s ({count/elapsed:.0f} tournaments/s)")

    # Analysis
    print(f"\n  Individual cycle count distributions:")
    for k in odd_ks:
        if k in cycle_dist and cycle_dist[k]:
            vals = sorted(cycle_dist[k].keys())
            print(f"    c_{k}: min={vals[0]}, max={vals[-1]}, distinct values={len(vals)}")
            if len(vals) <= 20:
                print(f"         values: {vals}")

    print(f"\n  FORCING THRESHOLDS:")
    print(f"  (p -> q): min c_p that GUARANTEES c_q > 0")
    print(f"  {'pair':>10} | {'threshold':>9} | {'max_c_p':>7} | {'ratio':>6} | details")
    print(f"  {'-'*10}-|-{'-'*9}-|-{'-'*7}-|-{'-'*6}-|--------")

    results = {}

    for (p, q), data in sorted(pair_data.items()):
        if not data:
            continue

        # For each c_p value, find min(c_q)
        min_cq_by_cp = {}
        max_cq_by_cp = {}
        for cp_val in sorted(data.keys()):
            cq_vals = data[cp_val]
            min_cq = min(cq_vals.keys())
            max_cq = max(cq_vals.keys())
            min_cq_by_cp[cp_val] = min_cq
            max_cq_by_cp[cp_val] = max_cq

        # Forcing threshold: smallest cp such that min(cq) > 0
        threshold = None
        max_cp = max(data.keys())
        for cp_val in sorted(data.keys()):
            if min_cq_by_cp[cp_val] > 0:
                threshold = cp_val
                break

        # Also find: max cp with cq=0 possible
        max_cp_with_cq0 = None
        for cp_val in sorted(data.keys(), reverse=True):
            if min_cq_by_cp[cp_val] == 0:
                max_cp_with_cq0 = cp_val
                break

        ratio = threshold / max_cp if threshold and max_cp > 0 else None
        ratio_str = f"{ratio:.3f}" if ratio else "N/A"

        # Compact detail: show min_cq for each cp near threshold
        detail_parts = []
        for cp_val in sorted(data.keys()):
            if threshold is not None and abs(cp_val - threshold) <= 3:
                detail_parts.append(f"c{p}={cp_val}:min(c{q})={min_cq_by_cp[cp_val]}")

        threshold_str = str(threshold) if threshold else "never"
        print(f"  c{p}->c{q}: {threshold_str:>8s} | {max_cp:>7} | {ratio_str:>6} | {'; '.join(detail_parts[:4])}")

        results[(p, q)] = {
            'threshold': threshold,
            'max_cp': max_cp,
            'max_cp_with_cq0': max_cp_with_cq0,
            'min_cq_by_cp': dict(min_cq_by_cp),
            'max_cq_by_cp': dict(max_cq_by_cp),
        }

    # Detailed tables
    for (p, q), data in sorted(pair_data.items()):
        if not data:
            continue
        r = results[(p, q)]
        print(f"\n  Full table: c_{p} -> c_{q} at n={n}")
        print(f"  {'c_'+str(p):>6} | {'min(c_'+str(q)+')':>10} | {'max(c_'+str(q)+')':>10} | {'#tournaments':>12}")
        print(f"  {'-'*6}-|-{'-'*10}-|-{'-'*10}-|-{'-'*12}")
        for cp_val in sorted(data.keys()):
            total_at_cp = sum(data[cp_val].values())
            min_cq = r['min_cq_by_cp'][cp_val]
            max_cq = r['max_cq_by_cp'][cp_val]
            marker = " *" if min_cq > 0 and (cp_val == r['threshold']) else ""
            print(f"  {cp_val:>6} | {min_cq:>10} | {max_cq:>10} | {total_at_cp:>12}{marker}")

    return results

# ================================================================
# Run for each n
# ================================================================

all_results = {}

# n=5: exhaustive (1024 tournaments, trivial)
all_results[5] = analyze_forcing(5, exhaustive=True)

# n=6: exhaustive (32768 tournaments, ~30s)
all_results[6] = analyze_forcing(6, exhaustive=True)

# n=7: exhaustive (2M tournaments) — c₅ via 21 subtournament DPs
# This is the computationally intensive one. If too slow, fallback to sampling.
print("\n\nAttempting n=7 exhaustive (2,097,152 tournaments)...")
print("If this takes >10 minutes, consider sampling instead.\n")
try:
    all_results[7] = analyze_forcing(7, exhaustive=True)
except KeyboardInterrupt:
    print("\n  Interrupted! Falling back to sampling for n=7.")
    all_results[7] = analyze_forcing(7, exhaustive=False, num_samples=100000)

# n=8: sampled
all_results[8] = analyze_forcing(8, exhaustive=False, num_samples=30000)

# n=9: sampled (small, c₉ counting is expensive)
all_results[9] = analyze_forcing(9, exhaustive=False, num_samples=5000)

# ================================================================
# Cross-n comparison: the recursive pattern
# ================================================================

print("\n\n" + "=" * 70)
print("CROSS-n COMPARISON: FORCING THRESHOLDS")
print("=" * 70)

# Table: threshold(n, p->q)
pairs = [(3, 5), (3, 7), (3, 9), (5, 7), (5, 9), (7, 9)]
print(f"\n  {'n':>3} | {'max c₃':>6} | ", end="")
for (p, q) in pairs:
    print(f"{'c'+str(p)+'->c'+str(q):>8} | ", end="")
print()
print(f"  {'-'*3}-|-{'-'*6}-|-", end="")
for _ in pairs:
    print(f"{'-'*8}-|-", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    max_c3 = max(r.get((3, 5), {}).get('min_cq_by_cp', {0: 0}).keys()) if (3, 5) in r else "?"
    print(f"  {n:>3} | {max_c3:>6} | ", end="")
    for (p, q) in pairs:
        if (p, q) in r:
            thresh = r[(p, q)]['threshold']
            max_cp = r[(p, q)]['max_cp']
            if thresh is not None:
                print(f"{thresh:>4}/{max_cp:<3} | ", end="")
            else:
                print(f"{'never':>8} | ", end="")
        else:
            print(f"{'N/A':>8} | ", end="")
    print()

# Also show thresholds as fractions of max
print(f"\n  Thresholds as fraction of max c_p:")
print(f"  {'n':>3} | ", end="")
for (p, q) in pairs:
    print(f"{'c'+str(p)+'->c'+str(q):>8} | ", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"  {n:>3} | ", end="")
    for (p, q) in pairs:
        if (p, q) in r:
            thresh = r[(p, q)]['threshold']
            max_cp = r[(p, q)]['max_cp']
            if thresh is not None and max_cp > 0:
                print(f"{thresh/max_cp:>8.3f} | ", end="")
            else:
                print(f"{'N/A':>8} | ", end="")
        else:
            print(f"{'N/A':>8} | ", end="")
    print()

# Show max c_p with c_q=0 (the "escape boundary")
print(f"\n  Max c_p with c_q = 0 still possible:")
print(f"  {'n':>3} | ", end="")
for (p, q) in pairs:
    print(f"{'c'+str(p)+'->c'+str(q):>8} | ", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"  {n:>3} | ", end="")
    for (p, q) in pairs:
        if (p, q) in r:
            max_0 = r[(p, q)]['max_cp_with_cq0']
            max_cp = r[(p, q)]['max_cp']
            if max_0 is not None:
                print(f"{max_0:>4}/{max_cp:<3} | ", end="")
            else:
                print(f"{'0/'+str(max_cp):>8} | ", end="")
        else:
            print(f"{'N/A':>8} | ", end="")
    print()
