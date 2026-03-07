#!/usr/bin/env python3
"""
Source/sink avoidance analysis for the dichotomy proof.

When we delete v from R (a P-source in the poisoning graph), we need
T-v to have no source or sink.

New source: w has score 1 in T and w -> v (score becomes 0 in T-v).
New sink: w has score n-2 in T and v -> w (score stays n-2, which is max in T-v).

This script analyzes:
1. How many P-sources exist in R?
2. How many are "blocked" by score-1 or score-(n-2) constraints?
3. Is there always an unblocked P-source?
4. If not, can we delete from A or B instead?

Also check: at LARGER n (10, 11), does the source/sink issue ever block ALL options?

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

def find_disjoint_pair(cycle_sets):
    for i in range(len(cycle_sets)):
        for j in range(i+1, len(cycle_sets)):
            if not (cycle_sets[i] & cycle_sets[j]):
                return cycle_sets[i], cycle_sets[j]
    return None, None

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

def get_vertex_3cycles(cycle_sets, v):
    return [c for c in cycle_sets if v in c]

def compute_poisoning_graph(cycle_sets, R, AB):
    edges = {}
    S = set()
    for w in R:
        w_cycles = get_vertex_3cycles(cycle_sets, w)
        if not w_cycles:
            continue
        has_ab_cycle = any(c <= ({w} | AB) for c in w_cycles)
        if has_ab_cycle:
            S.add(w)
            continue
        common_r = None
        for c in w_cycles:
            r_in_c = (c - {w}) & set(R)
            if common_r is None:
                common_r = r_in_c
            else:
                common_r = common_r & r_in_c
        if len(common_r) == 1:
            target = list(common_r)[0]
            edges[w] = target
        else:
            S.add(w)
    return edges, S

def check_deletion_cycle_rich(adj, n, v):
    remaining = [w for w in range(n) if w != v]
    for w in remaining:
        out_deg = sum(1 for u in remaining if u != w and adj[w][u])
        if out_deg == 0 or out_deg == len(remaining) - 1:
            return False
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
    for n in [9, 10]:
        random.seed(42)
        max_trials = 2000000 if n == 9 else 500000

        print(f"\n=== n={n}: Source/sink avoidance analysis ===")

        stats = {
            'mm2_total': 0,
            'p_sources_count': {},  # distribution of number of P-sources
            'blocked_all': 0,  # all P-sources blocked by source/sink
            'unblocked_exists': 0,
            'safe_from_AB': 0,
            'total_safe': 0,
            'total_unsafe': 0,
        }

        for trial in range(max_trials):
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
            if t3 > 20:
                continue

            all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
            if not all_in:
                continue

            c3_sets = find_3cycle_sets(adj, n)
            mm = max_matching_3cycles(c3_sets)
            if mm != 2:
                continue

            stats['mm2_total'] += 1

            A, B = find_disjoint_pair(c3_sets)
            AB = A | B
            R = [v for v in range(n) if v not in AB]

            edges, S = compute_poisoning_graph(c3_sets, R, AB)

            # Find P-sources (in-degree 0 in P)
            targets = set(edges.values())
            p_sources = [v for v in R if v not in targets]

            ns = len(p_sources)
            stats['p_sources_count'][ns] = stats['p_sources_count'].get(ns, 0) + 1

            # Check which P-sources can be safely deleted (no source/sink)
            safe_sources = []
            for v in p_sources:
                if check_deletion_cycle_rich(adj, n, v):
                    safe_sources.append(v)

            if safe_sources:
                stats['unblocked_exists'] += 1
                stats['total_safe'] += 1
            else:
                stats['blocked_all'] += 1
                # Try A or B
                ab_safe = False
                for v in list(AB):
                    if check_deletion_cycle_rich(adj, n, v):
                        ab_safe = True
                        break
                if ab_safe:
                    stats['safe_from_AB'] += 1
                    stats['total_safe'] += 1
                else:
                    stats['total_unsafe'] += 1
                    print(f"  !!! TRULY UNSAFE: trial={trial}, t3={t3}")
                    print(f"      scores={tuple(sorted(scores))}")
                    print(f"      P-sources={p_sources}, targets={targets}")

            if stats['mm2_total'] % 5000 == 0:
                print(f"  n={n}: {stats['mm2_total']} mm=2, "
                      f"safe={stats['total_safe']}, "
                      f"blocked_all={stats['blocked_all']}, "
                      f"safe_AB={stats['safe_from_AB']}, "
                      f"unsafe={stats['total_unsafe']}")

        print(f"\n=== n={n} RESULTS ===")
        print(f"Total mm=2: {stats['mm2_total']}")
        print(f"P-source count dist: {dict(sorted(stats['p_sources_count'].items()))}")
        print(f"Unblocked P-source exists: {stats['unblocked_exists']} "
              f"({100*stats['unblocked_exists']/max(1,stats['mm2_total']):.1f}%)")
        print(f"All P-sources blocked: {stats['blocked_all']}")
        print(f"  Of which safe from AB: {stats['safe_from_AB']}")
        print(f"TOTAL UNSAFE: {stats['total_unsafe']}")


if __name__ == "__main__":
    main()
