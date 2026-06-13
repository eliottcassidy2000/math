#!/usr/bin/env python3
"""
even_graph_fiber_analysis_s260.py — Deeper analysis of the fiber structure
opus-2026-03-23-S260

Key questions:
1. What determines fiber size? Is it related to |Aut| of the even graph?
2. Does the same-even fraction stabilize? What is the asymptotic value?
3. How does H distribute within fibers? Do fibers concentrate H?
4. What's the cut/cycle information ratio at larger n?
5. Does the cycle-space projection respect the SC/NS partition?
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import comb, factorial, gcd
import subprocess
import numpy as np

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

def build_cut_space_basis(n, P):
    """Build RREF of cut space of K_n."""
    m = len(P)
    basis = []
    for v in range(n-1):
        vec = [0]*m
        for k, (i,j) in enumerate(P):
            if i == v or j == v:
                vec[k] = 1
        basis.append(vec)

    # Gaussian elimination
    pivot_cols = []
    row = 0
    for col in range(m):
        found = -1
        for r in range(row, len(basis)):
            if basis[r][col]:
                found = r
                break
        if found == -1:
            continue
        basis[row], basis[found] = basis[found], basis[row]
        for r in range(len(basis)):
            if r != row and basis[r][col]:
                basis[r] = [(basis[r][k] + basis[row][k]) % 2 for k in range(m)]
        pivot_cols.append(col)
        row += 1

    return basis[:row], pivot_cols

def project_to_cycle_space(vec, basis, pivot_cols, m):
    """Project a tournament vector onto cycle space."""
    v_list = list(vec)
    v_cut = [0]*m
    for r_idx, pc in enumerate(pivot_cols):
        if v_list[pc]:
            v_cut = [(v_cut[k] + basis[r_idx][k]) % 2 for k in range(m)]
            v_list = [(v_list[k] + basis[r_idx][k]) % 2 for k in range(m)]
    v_cycle = tuple((vec[k] + v_cut[k]) % 2 for k in range(m))
    return v_cycle, tuple(v_cut)

def canonical_form(adj, n, directed=False):
    """Canonical form of a graph."""
    scores = [sum(adj[i][j] for j in range(n) if (not directed or i != j)) for i in range(n)]
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
        if directed:
            sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        else:
            sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or sig < best:
            best = sig
        count += 1
        if count > 200000:
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

def aut_size(adj, n, directed=False):
    """Count automorphisms."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if directed:
                    if i == j: continue
                    if adj[i][j] != adj[perm[i]][perm[j]]:
                        ok = False
                        break
                else:
                    if i >= j: continue
                    if adj[i][j] != adj[perm[i]][perm[j]]:
                        ok = False
                        break
            if not ok: break
        if ok: count += 1
    return count

def analyze(n):
    print(f"\n{'='*70}")
    print(f"  n = {n}: FIBER ANALYSIS")
    print(f"{'='*70}")

    tournaments, P = get_tournaments(n)
    m = comb(n, 2)
    basis, pivots = build_cut_space_basis(n, P)
    cycle_dim = m - len(pivots)

    print(f"  m = {m}, cut dim = {len(pivots)}, cycle dim = {cycle_dim}")
    print(f"  Tournament classes: {len(tournaments)}")

    # Project each tournament to cycle space
    data = []
    cycle_groups = defaultdict(list)

    for tidx, adj in enumerate(tournaments):
        vec = tuple(adj[i][j] for (i,j) in P)
        v_cycle, v_cut = project_to_cycle_space(vec, basis, pivots, m)
        h = hamiltonian_paths(adj, n)
        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

        # Tournament complement
        comp = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j: comp[i][j] = 1 - adj[i][j]
        comp_canon = canonical_form(comp, n, directed=True)
        tourn_canon = canonical_form(adj, n, directed=True)
        is_sc = (tourn_canon == comp_canon)

        # Cycle projection complement
        comp_cycle = tuple(1 - x for x in v_cycle)  # bitwise NOT of cycle projection

        d = {
            'tidx': tidx, 'vec': vec, 'v_cycle': v_cycle, 'v_cut': v_cut,
            'H': h, 'scores': scores, 'is_sc': is_sc,
            'comp_cycle': comp_cycle,
        }
        data.append(d)
        cycle_groups[v_cycle].append(tidx)

    # Labeled even graph → canonical even graph
    even_canon_map = {}
    for v_cycle in cycle_groups:
        even_adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if v_cycle[k]:
                even_adj[i][j] = 1
                even_adj[j][i] = 1
        canon = canonical_form(even_adj, n)
        even_canon_map[v_cycle] = canon

    # Group by canonical even graph
    even_class_groups = defaultdict(list)  # canon → list of tournament indices
    for v_cycle, tids in cycle_groups.items():
        canon = even_canon_map[v_cycle]
        even_class_groups[canon].extend(tids)

    n_even_classes = len(even_class_groups)
    print(f"  Even graph classes (from projection): {n_even_classes}")

    # Fiber analysis
    fiber_sizes = sorted([len(tids) for tids in even_class_groups.values()])
    print(f"  Fiber size distribution: {fiber_sizes}")
    print(f"  Max fiber: {max(fiber_sizes)}, Mean fiber: {len(tournaments)/n_even_classes:.2f}")

    # H statistics within fibers
    print(f"\n  H STATISTICS WITHIN FIBERS:")
    total_H_var_within = 0
    total_H_var_between = 0
    all_H = [data[t]['H'] for t in range(len(tournaments))]
    mean_H = np.mean(all_H)
    total_H_var = np.var(all_H)

    for canon, tids in even_class_groups.items():
        h_vals = [data[t]['H'] for t in tids]
        fiber_mean = np.mean(h_vals)
        total_H_var_within += sum((h - fiber_mean)**2 for h in h_vals)
        total_H_var_between += len(h_vals) * (fiber_mean - mean_H)**2

    total_N = len(tournaments)
    within_fraction = total_H_var_within / (total_N * total_H_var) if total_H_var > 0 else 0
    between_fraction = total_H_var_between / (total_N * total_H_var) if total_H_var > 0 else 0

    print(f"  Total H variance: {total_H_var:.1f}")
    print(f"  Within-fiber fraction: {within_fraction:.4f} ({100*within_fraction:.1f}%)")
    print(f"  Between-fiber fraction: {between_fraction:.4f} ({100*between_fraction:.1f}%)")
    print(f"  (within = H variance explained by scores; between = H variance explained by cycle structure)")

    # Score diversity within fibers
    print(f"\n  SCORE DIVERSITY WITHIN FIBERS:")
    max_score_div = 0
    for canon, tids in sorted(even_class_groups.items(), key=lambda x: -len(x[1])):
        if len(tids) > 1:
            scores_set = set(data[t]['scores'] for t in tids)
            sc_count = sum(1 for t in tids if data[t]['is_sc'])
            max_score_div = max(max_score_div, len(scores_set))
            if len(tids) >= 5:
                print(f"    Fiber size {len(tids)}: {len(scores_set)} distinct scores, "
                      f"{sc_count} SC, H range [{min(data[t]['H'] for t in tids)}, "
                      f"{max(data[t]['H'] for t in tids)}]")
    print(f"  Max score diversity in any fiber: {max_score_div}")

    # SC enrichment in fibers
    n_sc = sum(1 for d in data if d['is_sc'])
    sc_rate = n_sc / len(data)
    print(f"\n  SC ENRICHMENT:")
    print(f"  Overall SC rate: {n_sc}/{len(data)} = {sc_rate:.3f}")

    singleton_sc = sum(1 for canon, tids in even_class_groups.items()
                       if len(tids) == 1 and data[tids[0]]['is_sc'])
    singleton_total = sum(1 for tids in even_class_groups.values() if len(tids) == 1)
    if singleton_total > 0:
        print(f"  SC rate in singleton fibers: {singleton_sc}/{singleton_total} = "
              f"{singleton_sc/singleton_total:.3f}")

    large_fiber_sc = sum(1 for canon, tids in even_class_groups.items()
                         if len(tids) > 3
                         for t in tids if data[t]['is_sc'])
    large_fiber_total = sum(len(tids) for tids in even_class_groups.values()
                            if len(tids) > 3)
    if large_fiber_total > 0:
        print(f"  SC rate in large fibers (>3): {large_fiber_sc}/{large_fiber_total} = "
              f"{large_fiber_sc/large_fiber_total:.3f}")

    # Complement of cycle projection
    print(f"\n  COMPLEMENT INTERACTION:")
    # Does complementing a tournament complement its cycle projection?
    comp_matches = 0
    comp_mismatches = 0
    for d in data:
        comp_adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j: comp_adj[i][j] = 1 - tournaments[d['tidx']][i][j]
        comp_vec = tuple(comp_adj[i][j] for (i,j) in P)
        comp_v_cycle, _ = project_to_cycle_space(comp_vec, basis, pivots, m)

        # Is comp_v_cycle == complement of d['v_cycle']?
        bit_comp = tuple(1 - x for x in d['v_cycle'])
        if comp_v_cycle == bit_comp:
            comp_matches += 1
        else:
            comp_mismatches += 1

    print(f"  complement(T) maps to complement(cycle_proj(T)): "
          f"{comp_matches}/{len(data)} ({100*comp_matches/len(data):.1f}%)")
    print(f"  Mismatches: {comp_mismatches}")

    # Labeled even graph orbit sizes
    print(f"\n  ORBIT STRUCTURE:")
    # For each canonical even graph, count how many labeled even graphs map to it
    labeled_even = defaultdict(int)
    for v_cycle in cycle_groups:
        canon = even_canon_map[v_cycle]
        labeled_even[canon] += len(cycle_groups[v_cycle])

    # Wait, that's labeled tournaments per canonical even graph.
    # What we want is labeled EVEN GRAPHS per canonical even graph.
    labeled_even_graphs = defaultdict(int)
    for v_cycle in cycle_groups:
        canon = even_canon_map[v_cycle]
        labeled_even_graphs[canon] += 1  # each v_cycle is a labeled even graph

    orbit_sizes = sorted(labeled_even_graphs.values())
    print(f"  Even graph orbit sizes: {orbit_sizes}")
    print(f"  Sum of orbit sizes: {sum(orbit_sizes)} (should be {2**cycle_dim} labeled even graphs)")

    # Tournament orbit sizes per even class (total labeled tournaments in fiber)
    tourn_per_even = []
    for canon, tids in even_class_groups.items():
        total_labeled = sum(factorial(n) // aut_size(tournaments[t], n, directed=True) for t in tids) if n <= 6 else None
        tourn_per_even.append((len(tids), labeled_even_graphs.get(canon, 0), total_labeled))

    if n <= 6:
        print(f"\n  FIBER DETAIL (tournament_classes, labeled_even_graphs, labeled_tournaments):")
        for fc, le, lt in sorted(tourn_per_even, key=lambda x: -x[0]):
            if fc > 1:
                print(f"    {fc} classes, {le} labeled even graphs, {lt} labeled tournaments")

    return {
        'n': n, 'V_tourn': len(tournaments), 'V_even': n_even_classes,
        'cycle_dim': cycle_dim,
        'within_var': within_fraction, 'between_var': between_fraction,
        'fiber_sizes': fiber_sizes,
        'comp_matches': comp_matches, 'total': len(data),
    }

if __name__ == '__main__':
    print("=" * 70)
    print("  EVEN GRAPH FIBER ANALYSIS")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    results = []
    for n in range(3, 8):
        r = analyze(n)
        results.append(r)

    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"{'n':>3} {'V_t':>5} {'V_e':>5} {'ratio':>6} {'cycle_dim':>9} {'H_within':>9} {'H_between':>10} {'comp%':>6}")
    for r in results:
        print(f"{r['n']:>3} {r['V_tourn']:>5} {r['V_even']:>5} "
              f"{r['V_tourn']/r['V_even']:.3f} {r['cycle_dim']:>9} "
              f"{100*r['within_var']:>8.1f}% {100*r['between_var']:>9.1f}% "
              f"{100*r['comp_matches']/r['total']:>5.1f}%")

    print("\nDONE.")
