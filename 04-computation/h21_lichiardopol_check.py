#!/usr/bin/env python3
"""
Check: at n=9, cycle-rich tournaments with min outdegree >= 5 vs < 5.

Lichiardopol (proved by Bang-Jensen et al. for l=3):
  min outdegree >= (3-1)*3 - 1 = 5 => 3 disjoint 3-cycles.

So if min_outdeg >= 5: Part C blocks H=21.
If min_outdeg <= 4: need safe deletion argument.

This script checks:
1. What fraction of cycle-rich n=9 tournaments have min_outdeg >= 5?
2. For min_outdeg <= 4: does mm >= 3 ever hold anyway?
3. For min_outdeg <= 4 and mm <= 2: does safe deletion always work?

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from math import comb
import random

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def find_3cycle_sets(adj, n):
    cycle_sets = []
    for vs in combinations(range(n), 3):
        a, b, c = vs
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            cycle_sets.append(frozenset(vs))
    return cycle_sets

def max_matching_3cycles(cycle_sets):
    nc = len(cycle_sets)
    if nc == 0:
        return 0
    for a in range(nc):
        for b in range(a+1, nc):
            if cycle_sets[a] & cycle_sets[b]:
                continue
            for c in range(b+1, nc):
                if (not (cycle_sets[c] & cycle_sets[a])) and \
                   (not (cycle_sets[c] & cycle_sets[b])):
                    return 3
    for a in range(nc):
        for b in range(a+1, nc):
            if not (cycle_sets[a] & cycle_sets[b]):
                return 2
    return 1 if nc > 0 else 0

def check_deletion_cycle_rich(adj, n, v):
    """Check if T-v is cycle-rich (every remaining vertex in 3-cycle, no src/sink)."""
    remaining = [w for w in range(n) if w != v]
    # Check source/sink
    for w in remaining:
        out_deg = sum(1 for u in remaining if u != w and adj[w][u])
        if out_deg == 0 or out_deg == len(remaining) - 1:
            return False
    # Check every vertex in 3-cycle
    for w in remaining:
        found = False
        out_w = [j for j in remaining if j != w and adj[w][j]]
        in_w = [j for j in remaining if j != w and adj[j][w]]
        for u in out_w:
            for x in in_w:
                if u != x and adj[u][x]:
                    found = True
                    break
            if found:
                break
        if not found:
            return False
    return True

def main():
    n = 9
    random.seed(42)

    print(f"=== n={n}: Lichiardopol threshold check ===")
    print(f"Lichiardopol for k=3 disjoint 3-cycles: min outdeg >= 5")

    stats = {
        'total_cr': 0,
        'min_od_ge5': 0,
        'min_od_le4': 0,
        'min_od_le4_mm3': 0,
        'min_od_le4_mm2': 0,
        'min_od_le4_mm1': 0,
        'min_od_le4_safe_del': 0,
        'min_od_le4_no_safe': 0,
    }

    min_od_dist = {}

    for trial in range(2000000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        scores = [sum(adj[i]) for i in range(n)]
        if 0 in scores or (n-1) in scores:
            continue

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 15:
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        stats['total_cr'] += 1
        min_od = min(scores)
        min_od_dist[min_od] = min_od_dist.get(min_od, 0) + 1

        if min_od >= 5:
            stats['min_od_ge5'] += 1
        else:
            stats['min_od_le4'] += 1

            c3_sets = find_3cycle_sets(adj, n)
            mm = max_matching_3cycles(c3_sets)

            if mm >= 3:
                stats['min_od_le4_mm3'] += 1
            elif mm == 2:
                stats['min_od_le4_mm2'] += 1
                # Check safe deletion
                has_safe = False
                for v in range(n):
                    if check_deletion_cycle_rich(adj, n, v):
                        has_safe = True
                        break
                if has_safe:
                    stats['min_od_le4_safe_del'] += 1
                else:
                    stats['min_od_le4_no_safe'] += 1
                    print(f"  !!! NO SAFE DELETION at trial {trial}!")
                    print(f"      scores={tuple(sorted(scores))}, mm={mm}, t3={t3}")
            else:
                stats['min_od_le4_mm1'] += 1
                # mm=1: check safe deletion
                has_safe = False
                for v in range(n):
                    if check_deletion_cycle_rich(adj, n, v):
                        has_safe = True
                        break
                if has_safe:
                    stats['min_od_le4_safe_del'] += 1
                else:
                    stats['min_od_le4_no_safe'] += 1
                    print(f"  !!! MM=1 NO SAFE DELETION at trial {trial}!")

        if stats['total_cr'] % 2000 == 0:
            print(f"  {stats['total_cr']} cycle-rich: od>=5={stats['min_od_ge5']}, "
                  f"od<=4={stats['min_od_le4']} (mm3={stats['min_od_le4_mm3']}, "
                  f"mm2={stats['min_od_le4_mm2']}, mm1={stats['min_od_le4_mm1']}, "
                  f"safe={stats['min_od_le4_safe_del']}, nosafe={stats['min_od_le4_no_safe']})")

    print(f"\n=== RESULTS ===")
    print(f"Total cycle-rich: {stats['total_cr']}")
    print(f"Min outdeg >= 5 (Lichiardopol): {stats['min_od_ge5']} "
          f"({100*stats['min_od_ge5']/max(1,stats['total_cr']):.1f}%)")
    print(f"Min outdeg <= 4: {stats['min_od_le4']} "
          f"({100*stats['min_od_le4']/max(1,stats['total_cr']):.1f}%)")
    print(f"  Of which mm >= 3: {stats['min_od_le4_mm3']}")
    print(f"  Of which mm = 2: {stats['min_od_le4_mm2']}")
    print(f"  Of which mm = 1: {stats['min_od_le4_mm1']}")
    print(f"  Safe deletion found: {stats['min_od_le4_safe_del']}")
    print(f"  NO safe deletion: {stats['min_od_le4_no_safe']}")
    print(f"\nMin outdeg distribution: {dict(sorted(min_od_dist.items()))}")

    if stats['min_od_le4_no_safe'] == 0:
        print(f"\n*** ALL cycle-rich with min_od <= 4 have safe deletion! ***")
        print(f"*** Combined with Lichiardopol (min_od >= 5 => 3 disjoint 3-cycles): ***")
        print(f"*** DICHOTOMY HOLDS at n={n} ***")


if __name__ == "__main__":
    main()
