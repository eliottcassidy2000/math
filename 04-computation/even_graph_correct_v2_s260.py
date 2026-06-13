#!/usr/bin/env python3
"""
even_graph_correct_v2_s260.py — Correct cycle-space projection via brute force
opus-2026-03-23-S260

Uses brute-force search over 2^(n-1) cut-space elements to find the unique
decomposition v = v_cut + v_cycle where v_cycle has all even degrees.

This avoids the MM^T singularity issue at even n.
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import comb, factorial
import subprocess
import sys

def get_tournaments(n):
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]
    tournaments = []
    for line in lines:
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if line[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
        tournaments.append(adj)
    return tournaments, P

def build_cut_basis(n, P):
    """Build cut space basis: delta(i) for i=0,...,n-2."""
    m = len(P)
    deltas = []
    for i in range(n-1):
        d = [0]*m
        for k, (a,b) in enumerate(P):
            if a == i or b == i:
                d[k] = 1
        deltas.append(d)
    return deltas

def cycle_projection_bruteforce(vec, n, P, deltas):
    """Find the unique cycle-space projection of vec."""
    m = len(P)
    # Try all 2^(n-1) elements of cut space
    for mask in range(1 << (n-1)):
        v_cut = [0]*m
        for i in range(n-1):
            if mask & (1 << i):
                for k in range(m):
                    v_cut[k] ^= deltas[i][k]

        v_cycle = tuple((vec[k] ^ v_cut[k]) for k in range(m))

        # Check all degrees even
        ok = True
        for v in range(n):
            deg = 0
            for k, (a,b) in enumerate(P):
                if (a == v or b == v) and v_cycle[k]:
                    deg += 1
            if deg % 2 != 0:
                ok = False
                break

        if ok:
            return v_cycle, tuple(v_cut)

    raise ValueError("No valid cycle projection found!")

def canonical_undirected(adj_bits, n, P):
    """Canonical form of an undirected graph given as upper-triangle bits."""
    # Convert to adjacency matrix
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(P):
        if adj_bits[k]:
            adj[i][j] = 1
            adj[j][i] = 1

    degrees = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    groups = defaultdict(list)
    for v in range(n):
        groups[degrees[v]].append(v)
    sorted_groups = [groups[d] for d in sorted(groups.keys())]

    best = None
    count = 0
    def gen_perms(gs):
        if not gs:
            yield []
            return
        for p in permutations(gs[0]):
            for rest in gen_perms(gs[1:]):
                yield list(p) + rest

    for perm in gen_perms(sorted_groups):
        sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or sig < best:
            best = sig
        count += 1
        if count > 500000:
            break
    return best

def canonical_tournament(adj, n):
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    groups = defaultdict(list)
    for v in range(n):
        groups[scores[v]].append(v)
    sorted_groups = [groups[s] for s in sorted(groups.keys())]

    best = None
    count = 0
    def gen_perms(gs):
        if not gs:
            yield []
            return
        for p in permutations(gs[0]):
            for rest in gen_perms(gs[1:]):
                yield list(p) + rest

    for perm in gen_perms(sorted_groups):
        sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or sig < best:
            best = sig
        count += 1
        if count > 500000:
            break
    return best

def hamiltonian_paths(adj, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u), u)] = dp.get((S | (1 << u), u), 0) + val
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def analyze(n):
    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"{'#'*70}")

    tournaments, P = get_tournaments(n)
    m = comb(n, 2)
    deltas = build_cut_basis(n, P)
    cycle_dim = m - (n - 1)

    print(f"  m = {m}, cut_dim = {n-1}, cycle_dim = {cycle_dim}")
    print(f"  Tournament classes: {len(tournaments)}")

    # Project each tournament to cycle space
    data = []
    cycle_class_groups = defaultdict(list)  # even_canon -> list of tidx

    for tidx, adj in enumerate(tournaments):
        vec = tuple(adj[i][j] for (i,j) in P)
        v_cycle, v_cut = cycle_projection_bruteforce(vec, n, P, deltas)

        h = hamiltonian_paths(adj, n)
        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

        # SC check
        comp = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j: comp[i][j] = 1 - adj[i][j]
        comp_canon = canonical_tournament(comp, n)
        tourn_canon = canonical_tournament(adj, n)
        is_sc = (tourn_canon == comp_canon)

        # Canonical even graph
        even_canon = canonical_undirected(v_cycle, n, P)

        # Complement tournament's cycle projection
        comp_vec = tuple(comp[i][j] for (i,j) in P)
        comp_cycle, _ = cycle_projection_bruteforce(comp_vec, n, P, deltas)
        comp_even_canon = canonical_undirected(comp_cycle, n, P)

        # Even graph degree sequence
        degs = []
        for v in range(n):
            d = sum(1 for k, (a,b) in enumerate(P) if (a==v or b==v) and v_cycle[k])
            degs.append(d)
        degs = tuple(degs)

        data.append({
            'tidx': tidx, 'v_cycle': v_cycle, 'H': h,
            'scores': scores, 'is_sc': is_sc,
            'even_canon': even_canon, 'comp_even_canon': comp_even_canon,
            'degs': degs, 'n_edges': sum(v_cycle),
        })
        cycle_class_groups[even_canon].append(tidx)

    n_even_classes = len(cycle_class_groups)
    print(f"  Even graph iso classes: {n_even_classes}")
    ratio = len(tournaments) / n_even_classes
    print(f"  V_tournament/V_even = {len(tournaments)}/{n_even_classes} = {ratio:.3f}")

    # Fiber structure
    fiber_sizes = sorted(len(tids) for tids in cycle_class_groups.values())
    print(f"  Fiber sizes: {fiber_sizes}")

    # Complement interaction
    same_class_under_comp = 0
    diff_class_under_comp = 0
    for d in data:
        if d['even_canon'] == d['comp_even_canon']:
            same_class_under_comp += 1
        else:
            diff_class_under_comp += 1
    print(f"\n  COMPLEMENT: same even class: {same_class_under_comp}/{len(data)}, "
          f"different: {diff_class_under_comp}/{len(data)}")

    # H variance decomposition
    import numpy as np
    all_H = [d['H'] for d in data]
    mean_H = np.mean(all_H)
    total_var = np.var(all_H)

    within_ss = 0
    between_ss = 0
    for canon, tids in cycle_class_groups.items():
        h_vals = [data[t]['H'] for t in tids]
        fiber_mean = np.mean(h_vals)
        within_ss += sum((h - fiber_mean)**2 for h in h_vals)
        between_ss += len(h_vals) * (fiber_mean - mean_H)**2

    N = len(data)
    print(f"\n  H VARIANCE:")
    print(f"  Total: {total_var:.1f}")
    if total_var > 0:
        print(f"  Within-fiber (score-driven):   {100*within_ss/(N*total_var):.1f}%")
        print(f"  Between-fiber (cycle-driven):  {100*between_ss/(N*total_var):.1f}%")

    # Edge projection
    tourn_canons = {}
    for d in data:
        tc = canonical_tournament(tournaments[d['tidx']], n)
        tourn_canons[tc] = d['tidx']

    tourn_edges = set()
    same_even = 0
    diff_even = 0

    for d in data:
        adj = tournaments[d['tidx']]
        for p_idx in range(m):
            i, j = P[p_idx]
            new_adj = [row[:] for row in adj]
            new_adj[i][j] = 1 - adj[i][j]
            new_adj[j][i] = 1 - adj[j][i]
            nc = canonical_tournament(new_adj, n)
            if nc in tourn_canons:
                nb = tourn_canons[nc]
                if nb != d['tidx']:
                    a, b = min(d['tidx'], nb), max(d['tidx'], nb)
                    if (a,b) not in tourn_edges:
                        tourn_edges.add((a,b))
                        ea = data[a]['even_canon']
                        eb = data[b]['even_canon']
                        if ea == eb:
                            same_even += 1
                        else:
                            diff_even += 1

    total_te = same_even + diff_even
    print(f"\n  METAGRAPH EDGES:")
    print(f"  Total: {len(tourn_edges)}")
    if total_te > 0:
        print(f"  Same even class: {same_even} ({100*same_even/total_te:.1f}%)")
        print(f"  Diff even class: {diff_even} ({100*diff_even/total_te:.1f}%)")

    # Fiber details
    print(f"\n  TOP FIBERS:")
    for canon, tids in sorted(cycle_class_groups.items(), key=lambda x: -len(x[1]))[:5]:
        h_vals = [data[t]['H'] for t in tids]
        sc_count = sum(1 for t in tids if data[t]['is_sc'])
        scores_set = set(data[t]['scores'] for t in tids)
        ne = data[tids[0]]['n_edges']
        dg = data[tids[0]]['degs']
        print(f"    Size {len(tids)}: {ne} edges, degs={dg}, "
              f"H=[{min(h_vals)}..{max(h_vals)}], {sc_count} SC, {len(scores_set)} scores")

    return {
        'n': n, 'V_t': len(tournaments), 'V_e': n_even_classes,
        'ratio': ratio, 'cycle_dim': cycle_dim,
        'within': within_ss/(N*total_var) if total_var > 0 else 0,
        'between': between_ss/(N*total_var) if total_var > 0 else 0,
        'comp_same': same_class_under_comp,
        'same_even_pct': 100*same_even/total_te if total_te > 0 else 0,
    }

if __name__ == '__main__':
    print("=" * 70)
    print("  CORRECT EVEN GRAPH PROJECTION (brute force)")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    results = []
    for n in range(3, 8):
        try:
            r = analyze(n)
            results.append(r)
        except Exception as e:
            print(f"  ERROR at n={n}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("  SUMMARY TABLE")
    print("=" * 70)
    print(f"{'n':>3} {'V_t':>5} {'V_e':>5} {'ratio':>7} {'cyc':>4} "
          f"{'within':>7} {'between':>8} {'comp_same':>10} {'same_e':>7}")
    for r in results:
        print(f"{r['n']:>3} {r['V_t']:>5} {r['V_e']:>5} "
              f"{r['ratio']:>7.3f} {r['cycle_dim']:>4} "
              f"{100*r['within']:>6.1f}% {100*r['between']:>7.1f}% "
              f"{r['comp_same']:>10} {r['same_even_pct']:>6.1f}%")

    print("\nDONE.")
