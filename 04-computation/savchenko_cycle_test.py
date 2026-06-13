#!/usr/bin/env python3
"""
Test Savchenko's observation that c_k(DRT) is invariant across all DRTs of same order.

We already know from THM-028 that at n=7 (only one DRT: Paley), all cycle
counts are determined. But the question is:
  - Is c_k the same for ALL regular tournaments, or only DRTs?
  - How do cycle counts vary across isomorphism classes?

Also: verify the DRT phase-transition observation using our n=7 data.

At n=7, locally transitive tournament (LTT) = tournament where out-set and
in-set of every vertex are transitive. Compare cycle counts with DRT.

kind-pasteur-2026-03-06-S19
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict
from itertools import combinations

def build_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

def is_transitive(adj, verts):
    """Check if subtournament on verts is transitive (acyclic)."""
    n = len(verts)
    scores = []
    for i in verts:
        s = sum(adj[i][j] for j in verts if j != i)
        scores.append(s)
    scores.sort()
    return scores == list(range(n))

def is_locally_transitive(T, n):
    """Check if T is locally transitive (out-set and in-set of each vertex are transitive)."""
    for v in range(n):
        out_set = [u for u in range(n) if u != v and T[v][u]]
        in_set = [u for u in range(n) if u != v and T[u][v]]
        if not is_transitive(T, out_set):
            return False
        if not is_transitive(T, in_set):
            return False
    return True

def is_doubly_regular(T, n):
    """Check if T is doubly regular (score all (n-1)/2, co-score balanced)."""
    k = (n-1) // 2
    for v in range(n):
        if sum(T[v]) != k:
            return False
    # Check: for every pair (u,v), number of common out-neighbors is (n-3)/4
    target = (n - 3) // 4
    for u in range(n):
        for v in range(u+1, n):
            common = sum(1 for w in range(n) if w != u and w != v and T[u][w] and T[v][w])
            if common != target:
                return False
    return True

def count_directed_cycles_by_length(T, n):
    """Count directed Hamiltonian cycles by length using DP."""
    counts = {}
    for k in range(3, n+1, 2):
        total_dir = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
            v0 = verts[0]
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            num_cycles = 0
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    num_cycles += c
            total_dir += num_cycles
        counts[k] = total_dir
    return counts

# ============================================================
# Classify all n=7 tournaments by type
# ============================================================
print("=" * 70)
print("N=7: CYCLE COUNTS BY TOURNAMENT TYPE (DRT vs LTT vs OTHER REGULAR)")
print("=" * 70)

n = 7
m = n*(n-1)//2

drt_data = []
ltt_data = []
other_reg = []
count = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = sorted(sum(T[i]) for i in range(n))
    if scores != [3]*7:
        continue
    count += 1

    h = hamiltonian_path_count(T)
    dc = count_directed_cycles_by_length(T, n)

    is_drt = is_doubly_regular(T, n)
    is_lt = is_locally_transitive(T, n)

    entry = {'bits': bits, 'h': h, 'dc': dc, 'drt': is_drt, 'lt': is_lt}

    if is_drt:
        drt_data.append(entry)
    elif is_lt:
        ltt_data.append(entry)
    else:
        other_reg.append(entry)

    if count % 200 == 0:
        print(f"  {count} regular tournaments found...", flush=True)

print(f"\nTotal regular: {count}")
print(f"  DRT: {len(drt_data)}")
print(f"  Locally transitive: {len(ltt_data)}")
print(f"  Other regular: {len(other_reg)}")

# Report cycle counts for each type
print(f"\n--- DRT (Doubly Regular) ---")
if drt_data:
    dc_vals = set()
    for d in drt_data:
        dc_vals.add(tuple(sorted(d['dc'].items())))
    for dc_val in dc_vals:
        example = [d for d in drt_data if tuple(sorted(d['dc'].items())) == dc_val][0]
        n_tours = sum(1 for d in drt_data if tuple(sorted(d['dc'].items())) == dc_val)
        print(f"  {n_tours} tournaments, H={example['h']}, cycle counts: {dict(dc_val)}")

print(f"\n--- Locally Transitive ---")
if ltt_data:
    dc_vals = set()
    for d in ltt_data:
        dc_vals.add(tuple(sorted(d['dc'].items())))
    for dc_val in dc_vals:
        example = [d for d in ltt_data if tuple(sorted(d['dc'].items())) == dc_val][0]
        n_tours = sum(1 for d in ltt_data if tuple(sorted(d['dc'].items())) == dc_val)
        print(f"  {n_tours} tournaments, H={example['h']}, cycle counts: {dict(dc_val)}")

print(f"\n--- Other Regular ---")
if other_reg:
    dc_vals = set()
    for d in other_reg:
        dc_vals.add(tuple(sorted(d['dc'].items())))
    for dc_val in dc_vals:
        example = [d for d in other_reg if tuple(sorted(d['dc'].items())) == dc_val][0]
        n_tours = sum(1 for d in other_reg if tuple(sorted(d['dc'].items())) == dc_val)
        print(f"  {n_tours} tournaments, H={example['h']}, cycle counts: {dict(dc_val)}")

# Summary: is cycle count a DRT invariant?
print(f"\n--- Summary ---")
all_reg = drt_data + ltt_data + other_reg
dc_by_h = defaultdict(set)
for d in all_reg:
    dc_by_h[d['h']].add(tuple(sorted(d['dc'].items())))
for h in sorted(dc_by_h.keys(), reverse=True):
    print(f"  H={h}: {len(dc_by_h[h])} distinct cycle-count vectors")

print("\nDone.")
