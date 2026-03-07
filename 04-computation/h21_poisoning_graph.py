#!/usr/bin/env python3
"""
Poisoning graph analysis for mm=2 cycle-rich tournaments.

For cycle-rich T on n vertices with max 3-cycle matching = 2:
- A, B = two disjoint 3-cycles
- R = V \ (A union B), |R| >= 3
- Poisoning graph P on R: w->v iff ALL w's 3-cycles go through v

Key structural claim: out-degree(w) <= 1 in P for all w.
(If w->v and w->v', then {w,v,v'} is a 3-cycle in R, giving 3-matching with A,B.)

This script checks:
1. Does the all-permutation case (every R vertex has out-degree 1) occur?
2. How often is S (vertices with out-degree 0, having 3-cycle in {w}∪A∪B) non-empty?
3. When S = emptyset, can we safely delete from A or B?

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
    """All 3-cycle sets containing vertex v."""
    return [c for c in cycle_sets if v in c]

def compute_poisoning_graph(cycle_sets, R, AB):
    """
    Compute poisoning graph on R.
    w -> v iff ALL of w's 3-cycles contain v.
    """
    edges = {}  # w -> v
    S = set()  # vertices with out-degree 0 (have 3-cycle in {w}∪AB)

    for w in R:
        w_cycles = get_vertex_3cycles(cycle_sets, w)
        if not w_cycles:
            continue

        # Check if w has a 3-cycle within {w}∪AB (not using any R vertex besides w)
        has_ab_cycle = False
        for c in w_cycles:
            if c <= ({w} | AB):
                has_ab_cycle = True
                break

        if has_ab_cycle:
            S.add(w)
            continue

        # All w's 3-cycles use some R vertex besides w
        # Find the common R vertex (if unique)
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
        elif len(common_r) == 0:
            # w's 3-cycles go through different R vertices
            S.add(w)  # No single dependency

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
    n = 9
    random.seed(42)

    print(f"=== n={n}: Poisoning graph analysis for mm=2 ===")

    stats = {
        'mm2_total': 0,
        'S_nonempty': 0,
        'S_empty_permutation': 0,
        'S_empty_other': 0,
        'safe_from_R': 0,
        'safe_from_AB': 0,
        'no_safe': 0,
    }

    outdeg_dist = {}

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

        c3_sets = find_3cycle_sets(adj, n)
        mm = max_matching_3cycles(c3_sets)

        if mm != 2:
            continue

        stats['mm2_total'] += 1

        A, B = find_disjoint_pair(c3_sets)
        if A is None:
            continue

        AB = A | B
        R = [v for v in range(n) if v not in AB]

        edges, S = compute_poisoning_graph(c3_sets, R, AB)

        # Record out-degree distribution
        for w in R:
            od = 1 if w in edges else 0
            outdeg_dist[od] = outdeg_dist.get(od, 0) + 1

        if S:
            stats['S_nonempty'] += 1
        else:
            # All R vertices have out-degree 1 (permutation)
            if len(edges) == len(R):
                stats['S_empty_permutation'] += 1
            else:
                stats['S_empty_other'] += 1

        # Check if safe deletion exists from R
        safe_r = False
        for v in R:
            if check_deletion_cycle_rich(adj, n, v):
                safe_r = True
                break

        if safe_r:
            stats['safe_from_R'] += 1
        else:
            # Try A and B
            safe_ab = False
            for v in list(AB):
                if check_deletion_cycle_rich(adj, n, v):
                    safe_ab = True
                    break
            if safe_ab:
                stats['safe_from_AB'] += 1
            else:
                stats['no_safe'] += 1
                print(f"  !!! NO SAFE DELETION: trial={trial}, t3={t3}")
                print(f"      scores={tuple(sorted(scores))}")
                print(f"      A={A}, B={B}, R={R}")
                print(f"      Poisoning edges: {edges}")
                print(f"      S (out-deg 0): {S}")

        if stats['mm2_total'] <= 5:
            print(f"  Example #{stats['mm2_total']}:")
            print(f"    A={A}, B={B}, R={R}")
            print(f"    Poisoning edges: {edges}")
            print(f"    S (have 3-cycle in w+AB): {S}")
            print(f"    t3={t3}, scores={tuple(sorted(scores))}")

        if stats['mm2_total'] % 5000 == 0:
            print(f"  {stats['mm2_total']} mm=2: S_nonempty={stats['S_nonempty']}, "
                  f"perm={stats['S_empty_permutation']}, "
                  f"safe_R={stats['safe_from_R']}, safe_AB={stats['safe_from_AB']}, "
                  f"nosafe={stats['no_safe']}")

    print(f"\n=== RESULTS ===")
    print(f"Total mm=2 cycle-rich: {stats['mm2_total']}")
    print(f"S non-empty (some R vertex has 3-cycle in w+AB): {stats['S_nonempty']} "
          f"({100*stats['S_nonempty']/max(1,stats['mm2_total']):.1f}%)")
    print(f"S empty (all-permutation): {stats['S_empty_permutation']} "
          f"({100*stats['S_empty_permutation']/max(1,stats['mm2_total']):.1f}%)")
    print(f"Safe deletion from R: {stats['safe_from_R']} "
          f"({100*stats['safe_from_R']/max(1,stats['mm2_total']):.1f}%)")
    print(f"Safe deletion from AB only: {stats['safe_from_AB']} "
          f"({100*stats['safe_from_AB']/max(1,stats['mm2_total']):.1f}%)")
    print(f"NO safe deletion: {stats['no_safe']}")
    print(f"\nR vertex out-degree dist: {dict(sorted(outdeg_dist.items()))}")

    if stats['no_safe'] == 0:
        print(f"\n*** ALL mm=2 cycle-rich have safe deletion! ***")


if __name__ == "__main__":
    main()
