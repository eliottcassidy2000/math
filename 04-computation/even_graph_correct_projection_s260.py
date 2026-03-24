#!/usr/bin/env python3
"""
even_graph_correct_projection_s260.py — CORRECT GF(2) cycle-space projection
opus-2026-03-23-S260

The previous projection used RREF complement, which is NOT the orthogonal
complement (cycle space). This version uses the incidence matrix approach:

    v = v_cut + v_cycle where v_cut ∈ cut space, v_cycle ∈ cycle space

    Cut space = row space of incidence matrix M
    Cycle space = null space of M (= graphs with all even degrees)

    Projection: solve MM^T s = Mv, then v_cut = M^T s, v_cycle = v + v_cut.
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import comb, factorial
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

def build_incidence_matrix(n, P):
    """Build vertex-edge incidence matrix M (n x m) over GF(2)."""
    m = len(P)
    M = [[0]*m for _ in range(n)]
    for k, (i,j) in enumerate(P):
        M[i][k] = 1
        M[j][k] = 1
    return M

def solve_gf2(A, b, n_vars):
    """Solve Ax = b over GF(2). Returns one solution or None."""
    n_eq = len(A)
    # Augmented matrix [A | b]
    aug = [row[:] + [b[i]] for i, row in enumerate(A)]
    pivot_cols = []
    row_idx = 0
    for col in range(n_vars):
        found = -1
        for r in range(row_idx, n_eq):
            if aug[r][col]:
                found = r
                break
        if found == -1:
            continue
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        for r in range(n_eq):
            if r != row_idx and aug[r][col]:
                aug[r] = [(aug[r][k] + aug[row_idx][k]) % 2 for k in range(n_vars + 1)]
        pivot_cols.append(col)
        row_idx += 1

    # Check consistency
    for r in range(row_idx, n_eq):
        if aug[r][n_vars]:
            return None  # inconsistent

    # Read solution (free variables = 0)
    x = [0] * n_vars
    for ri, col in zip(range(len(pivot_cols)), pivot_cols):
        x[col] = aug[ri][n_vars]
    return x

def correct_cycle_projection(vec, n, P, M):
    """Correctly project tournament vector to cycle space over GF(2)."""
    m = len(P)

    # Compute syndrome Mv = degree parities
    Mv = [sum(M[i][k] * vec[k] for k in range(m)) % 2 for i in range(n)]

    # If all even already, v is in cycle space
    if all(x == 0 for x in Mv):
        return tuple(vec), tuple(0 for _ in range(m))

    # Build MM^T (n x n) over GF(2)
    MMT = [[sum(M[i][k] * M[j][k] for k in range(m)) % 2 for j in range(n)] for i in range(n)]

    # Solve MM^T s = Mv over GF(2)
    s = solve_gf2(MMT, Mv, n)
    if s is None:
        # Shouldn't happen for connected graphs, but try removing last equation
        s = solve_gf2(MMT[:n-1], Mv[:n-1], n)

    # v_cut = M^T s (mod 2)
    v_cut = tuple(sum(M[i][k] * s[i] for i in range(n)) % 2 for k in range(m))

    # v_cycle = v + v_cut (mod 2)
    v_cycle = tuple((vec[k] + v_cut[k]) % 2 for k in range(m))

    # Verify: all degrees even in v_cycle
    for v in range(n):
        deg = sum(v_cycle[k] for k, (i,j) in enumerate(P) if i == v or j == v)
        assert deg % 2 == 0, f"Vertex {v} has odd degree {deg} in cycle projection!"

    return v_cycle, v_cut

def canonical_undirected(adj, n):
    """Canonical form of an undirected graph under S_n."""
    scores = [sum(adj[i][j] for j in range(n) if j != i) for i in range(n)]
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
        sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or sig < best:
            best = sig
        count += 1
        if count > 200000:
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

def analyze(n):
    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"{'#'*70}")

    tournaments, P = get_tournaments(n)
    m = comb(n, 2)
    M = build_incidence_matrix(n, P)
    cycle_dim = m - (n - 1)  # cycle space dimension

    print(f"  m = {m}, cut_dim = {n-1}, cycle_dim = {cycle_dim}")
    print(f"  Tournament classes: {len(tournaments)}")
    print(f"  Labeled even graphs: 2^{cycle_dim} = {2**cycle_dim}")

    # Project each tournament
    data = []
    cycle_vectors = set()

    for tidx, adj in enumerate(tournaments):
        vec = tuple(adj[i][j] for (i,j) in P)
        v_cycle, v_cut = correct_cycle_projection(vec, n, P, M)

        h = hamiltonian_paths(adj, n)
        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

        # Complement
        comp = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j: comp[i][j] = 1 - adj[i][j]
        comp_canon = canonical_tournament(comp, n)
        tourn_canon = canonical_tournament(adj, n)
        is_sc = (tourn_canon == comp_canon)

        # Even graph as undirected graph
        even_adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if v_cycle[k]:
                even_adj[i][j] = 1
                even_adj[j][i] = 1
        even_canon = canonical_undirected(even_adj, n)

        # Check complement's cycle projection
        comp_vec = tuple(comp[i][j] for (i,j) in P)
        comp_cycle, _ = correct_cycle_projection(comp_vec, n, P, M)

        # Is complement cycle = bitwise complement of cycle?
        bit_comp = tuple(1 - x for x in v_cycle)
        comp_cycle_is_bitcomp = (comp_cycle == bit_comp)

        # Number of edges in even graph
        n_edges = sum(v_cycle)
        degs = tuple(sum(even_adj[i][j] for j in range(n) if j != i) for i in range(n))

        cycle_vectors.add(v_cycle)

        data.append({
            'tidx': tidx, 'v_cycle': v_cycle, 'v_cut': v_cut,
            'H': h, 'scores': scores, 'is_sc': is_sc,
            'even_canon': even_canon, 'n_edges': n_edges, 'degs': degs,
            'comp_cycle': comp_cycle, 'comp_is_bitcomp': comp_cycle_is_bitcomp,
        })

    print(f"  Distinct labeled even graphs from projections: {len(cycle_vectors)}")

    # Group by canonical even graph
    even_class_groups = defaultdict(list)
    for d in data:
        even_class_groups[d['even_canon']].append(d['tidx'])

    n_even_classes = len(even_class_groups)
    print(f"  Even graph iso classes: {n_even_classes}")
    print(f"  Ratio V_tournament / V_even = {len(tournaments)} / {n_even_classes} = "
          f"{len(tournaments)/n_even_classes:.3f}")

    # Verify all degrees are even
    all_even = all(all(d % 2 == 0 for d in dat['degs']) for dat in data)
    print(f"  All projections have even degrees: {all_even}")

    # Fiber structure
    fiber_sizes = sorted(len(tids) for tids in even_class_groups.values())
    print(f"  Fiber sizes: {fiber_sizes}")

    # Complement interaction
    comp_bitcomp_count = sum(1 for d in data if d['comp_is_bitcomp'])
    print(f"\n  COMPLEMENT INTERACTION:")
    print(f"  comp(T)_cycle == ~T_cycle: {comp_bitcomp_count}/{len(data)}")

    # Check: does complement tournament have SAME even graph class?
    comp_same = 0
    comp_diff = 0
    for d in data:
        comp_adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j: comp_adj[i][j] = 1 - tournaments[d['tidx']][i][j]
        comp_even_adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if d['comp_cycle'][k]:
                comp_even_adj[i][j] = 1
                comp_even_adj[j][i] = 1
        comp_even_canon = canonical_undirected(comp_even_adj, n)
        if comp_even_canon == d['even_canon']:
            comp_same += 1
        else:
            comp_diff += 1

    print(f"  comp(T) same even class: {comp_same}/{len(data)}")
    print(f"  comp(T) different even class: {comp_diff}/{len(data)}")

    # H variance decomposition
    all_H = [d['H'] for d in data]
    mean_H = np.mean(all_H)
    total_var = np.var(all_H)

    within_ss = 0
    between_ss = 0
    for canon, tids in even_class_groups.items():
        h_vals = [data[t]['H'] for t in tids]
        fiber_mean = np.mean(h_vals)
        within_ss += sum((h - fiber_mean)**2 for h in h_vals)
        between_ss += len(h_vals) * (fiber_mean - mean_H)**2

    N = len(data)
    print(f"\n  H VARIANCE DECOMPOSITION (cycle projection):")
    print(f"  Total var(H) = {total_var:.1f}")
    if total_var > 0:
        print(f"  Within-fiber (score variation):  {100*within_ss/(N*total_var):.1f}%")
        print(f"  Between-fiber (cycle variation): {100*between_ss/(N*total_var):.1f}%")

    # Fiber details for large fibers
    print(f"\n  LARGEST FIBERS:")
    for canon, tids in sorted(even_class_groups.items(), key=lambda x: -len(x[1]))[:8]:
        h_vals = [data[t]['H'] for t in tids]
        sc_count = sum(1 for t in tids if data[t]['is_sc'])
        scores_set = set(data[t]['scores'] for t in tids)
        n_edges = data[tids[0]]['n_edges']
        degs = data[tids[0]]['degs']
        print(f"    Size {len(tids)}: {n_edges} edges, degs={degs}, "
              f"H=[{min(h_vals)}..{max(h_vals)}], {sc_count} SC, "
              f"{len(scores_set)} score types")

    # Edge projection: tournament metagraph → even graph class partition
    tourn_to_canon = {}
    for d in data:
        tc = canonical_tournament(tournaments[d['tidx']], n)
        tourn_to_canon[tc] = d['tidx']

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
            if nc in tourn_to_canon:
                nb = tourn_to_canon[nc]
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

    total_edges = same_even + diff_even
    print(f"\n  METAGRAPH EDGE PROJECTION:")
    print(f"  Tournament edges: {len(tourn_edges)}")
    if total_edges > 0:
        print(f"  Same even class: {same_even} ({100*same_even/total_edges:.1f}%)")
        print(f"  Different even class: {diff_even} ({100*diff_even/total_edges:.1f}%)")

    return {
        'n': n, 'V_t': len(tournaments), 'V_e': n_even_classes,
        'cycle_dim': cycle_dim,
        'within_pct': 100*within_ss/(N*total_var) if total_var > 0 else 0,
        'between_pct': 100*between_ss/(N*total_var) if total_var > 0 else 0,
        'comp_same': comp_same, 'comp_diff': comp_diff,
        'same_even': same_even, 'diff_even': diff_even,
    }

if __name__ == '__main__':
    print("=" * 70)
    print("  CORRECT EVEN GRAPH CYCLE-SPACE PROJECTION")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    results = []
    for n in range(3, 8):
        try:
            r = analyze(n)
            results.append(r)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"{'n':>3} {'V_t':>5} {'V_e':>5} {'ratio':>6} {'cyc':>4} "
          f"{'within%':>8} {'between%':>9} {'comp_same':>10} {'same_e%':>8}")
    for r in results:
        total = r['same_even'] + r['diff_even']
        se_pct = 100*r['same_even']/total if total > 0 else 0
        print(f"{r['n']:>3} {r['V_t']:>5} {r['V_e']:>5} "
              f"{r['V_t']/r['V_e']:.3f} {r['cycle_dim']:>4} "
              f"{r['within_pct']:>7.1f}% {r['between_pct']:>8.1f}% "
              f"{r['comp_same']:>10} {se_pct:>7.1f}%")

    print("\nDONE.")
