#!/usr/bin/env python3
"""
cycle_forcing_v2_s24.py — opus-2026-04-05-S24

Fast cycle forcing threshold computation.

For each n and pair (p,q), find: min c_p that GUARANTEES c_q > 0.

Uses direct cycle enumeration (not DP) for speed.
Exhaustive at n=5,6. Sampled at n=7,8,9.

The key question at each n:
  "How many 3-cycles can you have while avoiding ALL 5-cycles?"
  "How many 3-cycles force a 7-cycle?"
  "How many 5-cycles force a 7-cycle?"
  etc.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

# Precompute permutations for direct cycle checking
PERMS = {}
for k in [3, 5, 7, 9]:
    PERMS[k] = list(permutations(range(k - 1)))

def bits_to_adj(bits, n):
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

def score_c3(A, n):
    """c₃ via score formula. O(n)."""
    total = n * (n-1) * (n-2) // 6
    for v in range(n):
        s = 0
        for j in range(n):
            s += A[v][j]
        total -= s * (s-1) // 2
    return total

def count_kcycles(A, n, k):
    """Count directed k-cycles. Fix smallest vertex in each k-subset."""
    count = 0
    perms = PERMS[k]
    for combo in combinations(range(n), k):
        start = combo[0]
        rest = combo[1:]
        for p in perms:
            # Check cycle: start -> rest[p[0]] -> rest[p[1]] -> ... -> rest[p[k-2]] -> start
            valid = True
            prev = start
            for idx in p:
                nxt = rest[idx]
                if not A[prev][nxt]:
                    valid = False
                    break
                prev = nxt
            if valid and A[prev][start]:
                count += 1
    return count

def has_kcycle(A, n, k):
    """Check if ANY directed k-cycle exists. Early termination."""
    perms = PERMS[k]
    for combo in combinations(range(n), k):
        start = combo[0]
        rest = combo[1:]
        for p in perms:
            prev = start
            valid = True
            for idx in p:
                nxt = rest[idx]
                if not A[prev][nxt]:
                    valid = False
                    break
                prev = nxt
            if valid and A[prev][start]:
                return True
    return False

def analyze_n(n, exhaustive=True, num_samples=100000):
    """Analyze cycle forcing thresholds at tournament size n."""
    m = n * (n - 1) // 2
    total = 1 << m
    max_k = n if n % 2 == 1 else n - 1
    odd_ks = [k for k in range(3, min(max_k + 1, 10), 2)]  # cap at 9

    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m}, cycle lengths checked: {odd_ks}")
    if exhaustive:
        print(f"EXHAUSTIVE: {total} tournaments")
    else:
        print(f"SAMPLED: {num_samples} tournaments")
    print(f"{'='*70}")

    # For each pair (p,q), track joint distribution
    # pair_data[(p,q)][cp_val] = Counter of cq values
    pair_data = {}
    for i, p in enumerate(odd_ks):
        for q in odd_ks[i+1:]:
            pair_data[(p, q)] = defaultdict(Counter)

    cycle_dist = {k: Counter() for k in odd_ks}

    t0 = time.time()
    processed = 0

    if exhaustive:
        source = range(total)
        target = total
    else:
        source = (random.randint(0, total - 1) for _ in range(num_samples))
        target = num_samples

    for bits in source:
        A = bits_to_adj(bits, n)

        # Compute all cycle counts
        cc = {}
        cc[3] = score_c3(A, n)
        for k in odd_ks:
            if k > 3:
                cc[k] = count_kcycles(A, n, k)

        for k in odd_ks:
            cycle_dist[k][cc[k]] += 1

        for (p, q) in pair_data:
            pair_data[(p, q)][cc[p]][cc[q]] += 1

        processed += 1
        if processed % 5000 == 0:
            elapsed = time.time() - t0
            rate = processed / elapsed
            eta = (target - processed) / rate if rate > 0 else 0
            print(f"  {processed}/{target} ({processed/target*100:.1f}%), {rate:.0f}/s, ETA {eta:.0f}s     ", end='\r')
            sys.stdout.flush()

    elapsed = time.time() - t0
    print(f"\n  Completed: {processed} tournaments in {elapsed:.1f}s ({processed/elapsed:.0f}/s)")

    # Report results
    print(f"\n  Cycle count ranges:")
    for k in odd_ks:
        vals = sorted(cycle_dist[k].keys())
        print(f"    c_{k}: {vals[0]} .. {vals[-1]} ({len(vals)} distinct values)")

    results = {}

    for (p, q), data in sorted(pair_data.items()):
        if not data:
            continue

        min_cq_by_cp = {}
        max_cq_by_cp = {}
        count_by_cp = {}
        for cp_val in sorted(data.keys()):
            cq_counter = data[cp_val]
            min_cq_by_cp[cp_val] = min(cq_counter.keys())
            max_cq_by_cp[cp_val] = max(cq_counter.keys())
            count_by_cp[cp_val] = sum(cq_counter.values())

        # Threshold: first cp where min(cq) > 0
        threshold = None
        max_cp = max(data.keys())
        for cp_val in sorted(data.keys()):
            if min_cq_by_cp[cp_val] > 0:
                threshold = cp_val
                break

        max_cp_with_cq0 = None
        for cp_val in sorted(data.keys(), reverse=True):
            if min_cq_by_cp[cp_val] == 0:
                max_cp_with_cq0 = cp_val
                break

        results[(p, q)] = {
            'threshold': threshold,
            'max_cp': max_cp,
            'max_cp_with_cq0': max_cp_with_cq0,
            'min_cq_by_cp': dict(min_cq_by_cp),
            'max_cq_by_cp': dict(max_cq_by_cp),
            'count_by_cp': dict(count_by_cp),
        }

    # Print forcing threshold summary
    print(f"\n  FORCING THRESHOLDS at n={n}:")
    for (p, q), r in sorted(results.items()):
        thresh = r['threshold']
        max_cp = r['max_cp']
        last_zero = r['max_cp_with_cq0']
        if thresh is not None:
            print(f"    c_{p} >= {thresh} => c_{q} > 0  (max c_{p} = {max_cp}, last c_{p} with c_{q}=0: {last_zero})")
        else:
            print(f"    c_{p} NEVER forces c_{q} > 0  (c_{q}=0 possible even at max c_{p}={max_cp})")

    # Print detailed tables
    for (p, q), r in sorted(results.items()):
        data = pair_data[(p, q)]
        print(f"\n  c_{p} -> c_{q} table at n={n}:")
        print(f"    {'c_'+str(p):>5} | {'min':>5} | {'max':>5} | {'count':>8} | note")
        print(f"    {'-'*5}-|-{'-'*5}-|-{'-'*5}-|-{'-'*8}-|-----")
        for cp_val in sorted(data.keys()):
            min_cq = r['min_cq_by_cp'][cp_val]
            max_cq = r['max_cq_by_cp'][cp_val]
            cnt = r['count_by_cp'][cp_val]
            note = ""
            if cp_val == r['threshold']:
                note = " <-- THRESHOLD (c_q always > 0)"
            elif min_cq == 0 and cp_val == r['max_cp_with_cq0']:
                note = " <-- last c_p with c_q=0 possible"
            print(f"    {cp_val:>5} | {min_cq:>5} | {max_cq:>5} | {cnt:>8} | {note}")

    return results


# ================================================================
# Run computations
# ================================================================

all_results = {}

# n=5: exhaustive (instant)
all_results[5] = analyze_n(5, exhaustive=True)

# n=6: exhaustive (~30s)
all_results[6] = analyze_n(6, exhaustive=True)

# n=7: sampled (exhaustive too slow for c₅, c₇)
all_results[7] = analyze_n(7, exhaustive=False, num_samples=200000)

# n=8: sampled
all_results[8] = analyze_n(8, exhaustive=False, num_samples=50000)

# n=9: sampled (c₉ counting is very expensive, use fewer samples)
all_results[9] = analyze_n(9, exhaustive=False, num_samples=5000)


# ================================================================
# Cross-n comparison
# ================================================================

print("\n\n" + "=" * 70)
print("CROSS-n SYNTHESIS")
print("=" * 70)

pairs_of_interest = [(3, 5), (3, 7), (3, 9), (5, 7), (5, 9), (7, 9)]

# Table 1: Forcing thresholds
print(f"\n  Table 1: Forcing threshold (min c_p guaranteeing c_q > 0)")
print(f"  {'':>4}", end="")
for (p, q) in pairs_of_interest:
    print(f" | {'c'+str(p)+'->c'+str(q):>9}", end="")
print()
print(f"  {'n':>4}", end="")
for _ in pairs_of_interest:
    print(f" | {'thr/max':>9}", end="")
print()
print(f"  {'----':>4}", end="")
for _ in pairs_of_interest:
    print(f"-|-{'-'*9}", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"  {n:>4}", end="")
    for (p, q) in pairs_of_interest:
        if (p, q) in r:
            thresh = r[(p, q)]['threshold']
            max_cp = r[(p, q)]['max_cp']
            if thresh is not None:
                print(f" | {thresh:>4}/{max_cp:<4}", end="")
            else:
                print(f" |  nev/{max_cp:<4}", end="")
        else:
            print(f" | {'N/A':>9}", end="")
    print()

# Table 2: Max c_p with c_q = 0
print(f"\n  Table 2: Max c_p where c_q = 0 is still possible")
print(f"  {'':>4}", end="")
for (p, q) in pairs_of_interest:
    print(f" | {'c'+str(p)+'->c'+str(q):>9}", end="")
print()
print(f"  {'n':>4}", end="")
for _ in pairs_of_interest:
    print(f" | {'max0/max':>9}", end="")
print()
print(f"  {'----':>4}", end="")
for _ in pairs_of_interest:
    print(f"-|-{'-'*9}", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"  {n:>4}", end="")
    for (p, q) in pairs_of_interest:
        if (p, q) in r:
            max0 = r[(p, q)].get('max_cp_with_cq0')
            max_cp = r[(p, q)]['max_cp']
            if max0 is not None:
                print(f" | {max0:>4}/{max_cp:<4}", end="")
            else:
                print(f" |    0/{max_cp:<4}", end="")
        else:
            print(f" | {'N/A':>9}", end="")
    print()

# Table 3: Threshold as fraction of max
print(f"\n  Table 3: Threshold / max_c_p (fraction)")
print(f"  {'':>4}", end="")
for (p, q) in pairs_of_interest:
    print(f" | {'c'+str(p)+'->c'+str(q):>9}", end="")
print()
print(f"  {'----':>4}", end="")
for _ in pairs_of_interest:
    print(f"-|-{'-'*9}", end="")
print()

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"  {n:>4}", end="")
    for (p, q) in pairs_of_interest:
        if (p, q) in r:
            thresh = r[(p, q)]['threshold']
            max_cp = r[(p, q)]['max_cp']
            if thresh is not None and max_cp > 0:
                frac = thresh / max_cp
                print(f" | {frac:>9.4f}", end="")
            else:
                print(f" | {'N/A':>9}", end="")
        else:
            print(f" | {'N/A':>9}", end="")
    print()

# Max values of each cycle count
print(f"\n  Table 4: Maximum cycle counts by n")
print(f"  {'n':>4} | {'max c₃':>7} | {'max c₅':>7} | {'max c₇':>7} | {'max c₉':>7} | C(n,3)")
print(f"  {'----':>4}-|-{'-'*7}-|-{'-'*7}-|-{'-'*7}-|-{'-'*7}-|-------")
for n in sorted(all_results.keys()):
    r = all_results[n]
    parts = [f"  {n:>4}"]
    for k in [3, 5, 7, 9]:
        found = False
        for (p, q) in r:
            if p == k:
                parts.append(f" {r[(p,q)]['max_cp']:>7}")
                found = True
                break
        if not found:
            # Check if k is the larger cycle
            for (p, q) in r:
                if q == k:
                    max_cq = max(r[(p,q)]['max_cq_by_cp'].values())
                    parts.append(f" {max_cq:>7}")
                    found = True
                    break
        if not found:
            parts.append(f" {'N/A':>7}")
    cn3 = n*(n-1)*(n-2)//6
    parts.append(f" {cn3:>7}")
    print(" |".join(parts))
