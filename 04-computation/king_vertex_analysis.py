#!/usr/bin/env python3
"""
king_vertex_analysis.py
opus-2026-05-27-S1

Rigorous investigation of king vertices, strong connectivity,
court-rivals structure, and their relationship to the OCF (H = I(Omega, 2)).

Topics covered:
  1. Strong connectivity characterization via the "cut theorem"
  2. H-increment lower bound from court-rivals structure
  3. Apex tile effect on H
  4. King position in Hamiltonian paths
  5. Regular tournaments and universal king property
  6. Score sequence and the Q-P axis
"""

import itertools
from collections import defaultdict

# ============================================================
# Core tournament utilities
# ============================================================

def all_tournaments(n):
    """Generate all tournaments on {0,...,n-1} as adjacency matrices."""
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(arcs)):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def outdegree(T, v):
    return sum(T[v])

def score_seq(T):
    n = len(T)
    return sorted(outdegree(T, v) for v in range(n))

def is_strongly_connected(T):
    """Check strong connectivity via BFS."""
    n = len(T)
    def reachable(src):
        vis = {src}
        q = [src]
        while q:
            u = q.pop()
            for w in range(n):
                if T[u][w] and w not in vis:
                    vis.add(w)
                    q.append(w)
        return vis
    return len(reachable(0)) == n and all(
        0 in reachable(v) for v in range(n)
    )

def hamilton_paths(T):
    """Count directed Hamiltonian paths (brute force)."""
    n = len(T)
    count = 0
    for perm in itertools.permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not T[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            count += 1
    return count

def odd_cycles(T):
    """Return all directed odd cycles (as frozensets of vertices)."""
    n = len(T)
    cycles = set()
    def dfs(start, v, path, pathset):
        for w in range(n):
            if T[v][w]:
                if w == start and len(path) >= 3 and len(path) % 2 == 1:
                    cycles.add(frozenset(path))
                elif w not in pathset:
                    pathset.add(w)
                    path.append(w)
                    dfs(start, w, path, pathset)
                    path.pop()
                    pathset.remove(w)
    for s in range(n):
        dfs(s, s, [s], {s})
    return list(cycles)

def conflict_graph(cycles):
    """Build conflict graph: two cycles adjacent iff they share a vertex."""
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def independence_poly_at_2(cycles):
    """Compute I(Omega(T), 2) directly."""
    if not cycles:
        return 1
    adj = conflict_graph(cycles)
    m = len(cycles)
    total = 0
    for mask in range(1 << m):
        subset = [i for i in range(m) if mask & (1 << i)]
        # Check if independent set
        ok = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if adj[subset[a]][subset[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**len(subset)
    return total

# ============================================================
# Tiling model utilities
# ============================================================

def tiling_to_tournament(n, tile_bits):
    """
    Construct tournament from tiling.
    Base path: n-1 -> n-2 -> ... -> 0 (0-indexed, base path k -> k-1).
    Tiles: (x, y) with x >= y+2, ordered as the staircase enumeration.
    tile_bits[k] = 0 means x->y (default), 1 means y->x (flipped).
    Returns adjacency matrix T where T[i][j]=1 means i->j.
    """
    tiles = []
    for y in range(n-2):
        for x in range(n-1, y+1, -1):
            tiles.append((x, y))
    assert len(tiles) == (n-1)*(n-2)//2

    T = [[0]*n for _ in range(n)]
    # Base path arcs
    for k in range(1, n):
        T[k][k-1] = 1
    # Tile arcs
    for idx, (x, y) in enumerate(tiles):
        if tile_bits[idx] == 0:
            T[x][y] = 1
        else:
            T[y][x] = 1
    return T

def apex_tile_index(n):
    """Index of the apex tile (n-1, 0) in 0-indexed vertices."""
    # Apex tile connects vertex n-1 (top) to vertex 0 (bottom)
    # Enumeration: for y in range(n-2): for x in range(n-1, y+1, -1)
    # (n-1, 0) appears when y=0, x=n-1 => first tile in y=0 row
    return 0  # (n-1, 0) is the first tile when y=0, x=n-1

def all_tilings(n):
    """Generate all 2^m tilings for n vertices."""
    m = (n-1)*(n-2)//2
    for bits in range(1 << m):
        tile_bits = [(bits >> k) & 1 for k in range(m)]
        yield tile_bits

# ============================================================
# THEOREM 1: Strong Connectivity Cut Theorem
# ============================================================
# In the tiling model (0-indexed vertices 0,...,n-1, base path k->k-1):
# The cut at position k separates {k,...,n-1} (upper) from {0,...,k-1} (lower).
# Tiles crossing cut k: tile (x,y) with x >= k and y < k (0-indexed).
# A cut is "satisfied" if at least one crossing tile is flipped (y->x direction).
# T is strongly connected iff every cut k (k=1,...,n-1) is satisfied.

def tiles_for_cut(n, k):
    """Return indices of tiles crossing cut k (upper={k,...,n-1}, lower={0,...,k-1})."""
    tiles = []
    idx = 0
    for y in range(n-2):
        for x in range(n-1, y+1, -1):
            if x >= k and y < k:
                tiles.append(idx)
            idx += 1
    return tiles

def is_strongly_connected_by_cuts(n, tile_bits):
    """
    Check strong connectivity using the cut theorem.
    Returns: (is_strong, list of failed cuts)
    """
    failed_cuts = []
    for k in range(1, n):
        cut_tiles = tiles_for_cut(n, k)
        # Check if any tile in this cut is "upward" (bit=1, y->x direction)
        satisfied = any(tile_bits[i] == 1 for i in cut_tiles)
        if not satisfied:
            failed_cuts.append(k)
    return len(failed_cuts) == 0, failed_cuts

def verify_cut_theorem(n):
    """Verify that the cut theorem correctly characterizes strong connectivity."""
    correct = 0
    total = 0
    errors = []
    for tile_bits in all_tilings(n):
        T = tiling_to_tournament(n, tile_bits)
        sc_direct = is_strongly_connected(T)
        sc_cuts, failed = is_strongly_connected_by_cuts(n, tile_bits)
        if sc_direct != sc_cuts:
            errors.append((tile_bits, sc_direct, sc_cuts, failed))
        else:
            correct += 1
        total += 1
    return correct, total, errors

# ============================================================
# THEOREM 2: Court-Rivals and H-Increment Lower Bound
# ============================================================
# For max-degree vertex Q: H(T) - H(T-Q) >= 2 * |N^-(Q)|
# where N^-(Q) = vertices that beat Q.
# This follows from Claim A + King theorem (every rival in >=1 3-cycle with Q).

def tournament_minus_v(T, v):
    """Remove vertex v from tournament (relabel vertices)."""
    n = len(T)
    vertices = [u for u in range(n) if u != v]
    m = len(vertices)
    Tv = [[0]*m for _ in range(m)]
    for i, u in enumerate(vertices):
        for j, w in enumerate(vertices):
            Tv[i][j] = T[u][w]
    return Tv

def verify_king_h_bound(n):
    """
    For each tournament: H(T) - H(T-Q) >= 2 * |N^-(Q)|
    where Q is the unique max-degree vertex (if not unique, use any max-degree vertex).
    Also check: H(T) - H(T-Q) >= 2 * min_over_max_deg_vertices |N^-(v)|
    """
    results = []
    for T in all_tournaments(n):
        # Find max-degree vertices
        degs = [outdegree(T, v) for v in range(n)]
        max_deg = max(degs)
        Q_candidates = [v for v in range(n) if degs[v] == max_deg]

        H_T = hamilton_paths(T)

        for Q in Q_candidates:
            Tv = tournament_minus_v(T, Q)
            H_Tv = hamilton_paths(Tv)

            # Rivals: vertices that beat Q
            rivals = [v for v in range(n) if v != Q and T[v][Q]]
            n_rivals = len(rivals)

            delta_H = H_T - H_Tv
            bound = 2 * n_rivals

            if delta_H < bound:
                results.append({
                    'n': n,
                    'T': T,
                    'Q': Q,
                    'deg_Q': max_deg,
                    'n_rivals': n_rivals,
                    'H_T': H_T,
                    'H_Tv': H_Tv,
                    'delta_H': delta_H,
                    'bound': bound,
                    'violation': True
                })
            else:
                results.append({
                    'n': n,
                    'Q': Q,
                    'deg_Q': max_deg,
                    'n_rivals': n_rivals,
                    'H_T': H_T,
                    'delta_H': delta_H,
                    'bound': bound,
                    'violation': False,
                    'excess': delta_H - bound
                })
    return results

# ============================================================
# THEOREM 3: Apex Tile Monotonicity
# ============================================================
# Conjecture: Flipping the apex tile from default (n-1 -> 0) to flipped (0 -> n-1)
# always increases or preserves H.

def apex_effect_on_H(n):
    """
    For each tiling with apex bit=0, compute H difference when apex is flipped to 1.
    Returns list of (tile_bits_no_apex, H_0, H_1, delta).
    """
    m = (n-1)*(n-2)//2
    apex_idx = apex_tile_index(n)
    results = []

    for tile_bits in all_tilings(n):
        if tile_bits[apex_idx] == 0:  # only consider "apex-default" tilings
            T0 = tiling_to_tournament(n, tile_bits)

            # Flip apex
            flipped = list(tile_bits)
            flipped[apex_idx] = 1
            T1 = tiling_to_tournament(n, flipped)

            H0 = hamilton_paths(T0)
            H1 = hamilton_paths(T1)
            delta = H1 - H0

            results.append({
                'tile_bits': tuple(tile_bits),
                'H0': H0,
                'H1': H1,
                'delta': delta,
                'sc0': is_strongly_connected(T0),
                'sc1': is_strongly_connected(T1),
                'score0': score_seq(T0),
                'score1': score_seq(T1),
            })
    return results

# ============================================================
# THEOREM 4: King Position in Hamiltonian Paths
# ============================================================
# For max-degree vertex Q: does Q appear more often in certain positions of HPs?
# Compute position distribution of Q in all HPs.

def king_position_distribution(T):
    """
    For each HP, record position of the max-degree vertex Q.
    Returns: (degree_of_Q, {position: count})
    """
    n = len(T)
    degs = [outdegree(T, v) for v in range(n)]
    max_deg = max(degs)
    Q = degs.index(max_deg)  # take first max-degree vertex

    pos_count = defaultdict(int)
    for perm in itertools.permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not T[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            pos = perm.index(Q)
            pos_count[pos] += 1

    return max_deg, dict(pos_count)

# ============================================================
# THEOREM 5: Universal King Property for Regular Tournaments
# ============================================================
# For regular tournaments (all degrees equal), every vertex is a king.
# Verify computationally and check if this implies the uniform position matrix.

def is_regular(T):
    n = len(T)
    degs = [outdegree(T, v) for v in range(n)]
    return len(set(degs)) == 1

def is_king(T, v):
    """Check if v can reach all other vertices in <=2 steps."""
    n = len(T)
    reachable = set()
    reachable.add(v)
    # 1-step
    for w in range(n):
        if T[v][w]:
            reachable.add(w)
    # 2-step
    for w in list(reachable):
        if w != v:
            for u in range(n):
                if T[w][u]:
                    reachable.add(u)
    return len(reachable) == n

def all_kings(T):
    """Return list of king vertices."""
    n = len(T)
    return [v for v in range(n) if is_king(T, v)]

def verify_regular_king_property(n):
    """
    For all tournaments: check if regular <=> every vertex is a king.
    Also check weaker: regular => every vertex is a king.
    """
    results = {'regular_all_king': 0, 'regular_not_all_king': 0,
               'notregular_all_king': 0, 'notregular_not_all_king': 0}

    for T in all_tournaments(n):
        reg = is_regular(T)
        kings = all_kings(T)
        all_king = len(kings) == n

        if reg and all_king:
            results['regular_all_king'] += 1
        elif reg and not all_king:
            results['regular_not_all_king'] += 1
        elif not reg and all_king:
            results['notregular_all_king'] += 1
        else:
            results['notregular_not_all_king'] += 1

    return results

# ============================================================
# THEOREM 6: Score Sequence and the Q-P Gap
# ============================================================
# The "Q-P gap" = d^+(Q) - d^+(P) = max_degree - min_degree.
# Conjecture: H(T) is monotone decreasing in the Q-P gap (on average).
# Also: compute correlation between gap and H.

def qp_gap(T):
    n = len(T)
    degs = [outdegree(T, v) for v in range(n)]
    return max(degs) - min(degs)

def analyze_qp_gap(n):
    """For each tournament, compute (Q-P gap, H, strong connectivity)."""
    data = []
    for T in all_tournaments(n):
        gap = qp_gap(T)
        H = hamilton_paths(T)
        sc = is_strongly_connected(T)
        ss = tuple(sorted(outdegree(T, v) for v in range(n)))
        data.append({'gap': gap, 'H': H, 'sc': sc, 'score': ss})
    return data

# ============================================================
# THEOREM 7: Court-Rivals Triangle Density Bounds
# ============================================================
# c3(Q) = #{3-cycles through Q} >= |N^-(Q)|
# Also: is there an upper bound? c3(Q) <= |N^+(Q)| * |N^-(Q)|

def count_3cycles_through(T, v):
    """Count directed 3-cycles passing through vertex v."""
    n = len(T)
    count = 0
    for a in range(n):
        if a == v or not T[v][a]:
            continue
        for b in range(n):
            if b == v or b == a or not T[b][v] or not T[a][b]:
                continue
            count += 1
    return count

def analyze_court_rivals_triangles(n):
    """
    For max-degree vertex Q: verify c3(Q) >= |N^-(Q)| and c3(Q) <= |N^+(Q)|*|N^-(Q)|.
    """
    results = []
    for T in all_tournaments(n):
        degs = [outdegree(T, v) for v in range(n)]
        max_deg = max(degs)

        for Q in [v for v in range(n) if degs[v] == max_deg]:
            court = [v for v in range(n) if v != Q and T[Q][v]]
            rivals = [v for v in range(n) if v != Q and T[v][Q]]
            c3 = count_3cycles_through(T, Q)

            lower_ok = c3 >= len(rivals)
            upper_ok = c3 <= len(court) * len(rivals)

            results.append({
                'court_size': len(court),
                'rivals_size': len(rivals),
                'c3_Q': c3,
                'lower_ok': lower_ok,
                'upper_ok': upper_ok,
                'lower_tight': (c3 == len(rivals)),
                'upper_tight': (c3 == len(court) * len(rivals))
            })
    return results

# ============================================================
# RUN ALL ANALYSES
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("KING VERTEX ANALYSIS — opus-2026-05-27-S1")
    print("=" * 70)

    # ---- THEOREM 1: Strong Connectivity Cut Theorem ----
    print("\n" + "="*60)
    print("THEOREM 1: Strong Connectivity Cut Theorem")
    print("="*60)
    print("Verifying: T is strongly connected <=> every cut k has an upward tile")
    for n in range(3, 8):
        correct, total, errors = verify_cut_theorem(n)
        print(f"  n={n}: {correct}/{total} correct, {len(errors)} errors")
        if errors:
            print(f"    ERRORS: {errors[:3]}")

    # Compute fraction of tilings that are strongly connected
    print("\nFraction of tilings that are strongly connected:")
    for n in range(3, 8):
        m = (n-1)*(n-2)//2
        total = 2**m
        sc_count = 0
        for tile_bits in all_tilings(n):
            sc, _ = is_strongly_connected_by_cuts(n, tile_bits)
            if sc:
                sc_count += 1
        print(f"  n={n}: {sc_count}/{total} = {sc_count/total:.4f}")

    # ---- THEOREM 2: King H-increment bound ----
    print("\n" + "="*60)
    print("THEOREM 2: H-increment Lower Bound from King Theorem")
    print("="*60)
    print("Verifying: H(T) - H(T-Q) >= 2 * |N^-(Q)|")
    for n in range(3, 7):
        results = verify_king_h_bound(n)
        violations = [r for r in results if r['violation']]
        tight = [r for r in results if not r['violation'] and r['excess'] == 0]
        print(f"  n={n}: {len(results)} cases, {len(violations)} violations, {len(tight)} tight cases")
        if violations:
            print(f"    VIOLATION: {violations[0]}")
        else:
            # Show excess distribution
            excesses = [r['excess'] for r in results if not r['violation']]
            if excesses:
                min_e = min(excesses)
                max_e = max(excesses)
                avg_e = sum(excesses) / len(excesses)
                print(f"    Excess stats: min={min_e}, max={max_e}, avg={avg_e:.2f}")

    # ---- THEOREM 3: Apex tile effect ----
    print("\n" + "="*60)
    print("THEOREM 3: Apex Tile Effect on H")
    print("="*60)
    print("Flipping apex tile (n-1 -> 0) to (0 -> n-1): does H always increase?")
    for n in range(3, 7):
        results = apex_effect_on_H(n)
        deltas = [r['delta'] for r in results]
        neg = [r for r in results if r['delta'] < 0]
        zero = [r for r in results if r['delta'] == 0]
        pos = [r for r in results if r['delta'] > 0]
        print(f"  n={n}: {len(results)} cases, delta<0: {len(neg)}, delta=0: {len(zero)}, delta>0: {len(pos)}")
        if neg:
            print(f"    DECREASES: {neg[:2]}")
        else:
            print(f"    delta stats: min={min(deltas)}, max={max(deltas)}, avg={sum(deltas)/len(deltas):.3f}")

        # Strong connectivity relationship
        sc_before = sum(1 for r in results if r['sc0'])
        sc_after = sum(1 for r in results if r['sc1'])
        print(f"    SC before apex flip: {sc_before}/{len(results)}, after: {sc_after}/{len(results)}")

    # ---- THEOREM 4: King position distribution ----
    print("\n" + "="*60)
    print("THEOREM 4: King Position in Hamiltonian Paths")
    print("="*60)
    for n in range(3, 6):
        pos_by_deg_excess = defaultdict(lambda: defaultdict(int))
        h_by_deg = defaultdict(list)
        for T in all_tournaments(n):
            degs = [outdegree(T, v) for v in range(n)]
            max_deg = max(degs)
            deg_excess = max_deg - (n-1)//2  # excess above average
            max_deg_Q, pos_dist = king_position_distribution(T)
            total_h = sum(pos_dist.values())
            if total_h > 0:
                frac_first = pos_dist.get(0, 0) / total_h
                frac_last = pos_dist.get(n-1, 0) / total_h
                h_by_deg[max_deg].append((total_h, frac_first, frac_last))

        print(f"\n  n={n}: Position of Q in Hamiltonian paths (grouped by max_degree)")
        for deg in sorted(h_by_deg.keys()):
            entries = h_by_deg[deg]
            avg_first = sum(e[1] for e in entries) / len(entries)
            avg_last = sum(e[2] for e in entries) / len(entries)
            print(f"    deg_Q={deg}: {len(entries)} tourn, avg_frac_first={avg_first:.3f}, avg_frac_last={avg_last:.3f}")

    # ---- THEOREM 5: Regular tournaments and universal kings ----
    print("\n" + "="*60)
    print("THEOREM 5: Regular Tournaments and Universal King Property")
    print("="*60)
    for n in range(3, 8):
        results = verify_regular_king_property(n)
        print(f"  n={n}: reg+all_king={results['regular_all_king']}, "
              f"reg+not_all_king={results['regular_not_all_king']}, "
              f"notregular+all_king={results['notregular_all_king']}, "
              f"notregular+not_all_king={results['notregular_not_all_king']}")

    # ---- THEOREM 6: Q-P gap vs H ----
    print("\n" + "="*60)
    print("THEOREM 6: Q-P Gap and H Distribution")
    print("="*60)
    for n in range(3, 7):
        data = analyze_qp_gap(n)
        by_gap = defaultdict(list)
        for d in data:
            by_gap[d['gap']].append(d['H'])

        print(f"\n  n={n}: Q-P gap distribution")
        for gap in sorted(by_gap.keys()):
            hs = by_gap[gap]
            print(f"    gap={gap}: count={len(hs)}, min_H={min(hs)}, max_H={max(hs)}, avg_H={sum(hs)/len(hs):.2f}")

        # Correlation between gap and H
        all_gaps = [d['gap'] for d in data]
        all_hs = [d['H'] for d in data]
        n_pts = len(data)
        mean_gap = sum(all_gaps)/n_pts
        mean_h = sum(all_hs)/n_pts
        cov = sum((g-mean_gap)*(h-mean_h) for g,h in zip(all_gaps,all_hs)) / n_pts
        var_gap = sum((g-mean_gap)**2 for g in all_gaps) / n_pts
        var_h = sum((h-mean_h)**2 for h in all_hs) / n_pts
        if var_gap > 0 and var_h > 0:
            corr = cov / (var_gap**0.5 * var_h**0.5)
            print(f"    Pearson corr(gap, H) = {corr:.4f}")

    # ---- THEOREM 7: Court-rivals triangle bounds ----
    print("\n" + "="*60)
    print("THEOREM 7: Court-Rivals 3-Cycle Density")
    print("="*60)
    print("Verifying: c3(Q) >= |N^-(Q)| and c3(Q) <= |N^+(Q)| * |N^-(Q)|")
    for n in range(3, 7):
        results = analyze_court_rivals_triangles(n)
        lower_fails = [r for r in results if not r['lower_ok']]
        upper_fails = [r for r in results if not r['upper_ok']]
        lower_tight = [r for r in results if r['lower_tight']]
        upper_tight = [r for r in results if r['upper_tight']]
        print(f"  n={n}: {len(results)} (T,Q) pairs")
        print(f"    Lower bound c3(Q)>=|rivals|: {len(lower_fails)} fails, {len(lower_tight)} tight")
        print(f"    Upper bound c3(Q)<=|court|*|rivals|: {len(upper_fails)} fails, {len(upper_tight)} tight")

        # Also: what is the ratio c3(Q) / |rivals| distribution?
        ratios = [r['c3_Q'] / r['rivals_size'] for r in results if r['rivals_size'] > 0]
        if ratios:
            print(f"    c3(Q)/|rivals|: min={min(ratios):.2f}, max={max(ratios):.2f}, avg={sum(ratios)/len(ratios):.2f}")

        # Average c3(Q) vs average |court|*|rivals|
        avg_c3 = sum(r['c3_Q'] for r in results) / len(results)
        avg_prod = sum(r['court_size']*r['rivals_size'] for r in results) / len(results)
        print(f"    avg c3(Q) = {avg_c3:.2f}, avg |court|*|rivals| = {avg_prod:.2f}")

    # ---- BONUS: Q-P tightness analysis ----
    print("\n" + "="*60)
    print("BONUS: When is H(T) - H(T-Q) = 2 * |N^-(Q)| EXACTLY? (tight bound)")
    print("="*60)
    for n in range(3, 6):
        tight_cases = []
        for T in all_tournaments(n):
            degs = [outdegree(T, v) for v in range(n)]
            max_deg = max(degs)
            Q_candidates = [v for v in range(n) if degs[v] == max_deg]
            H_T = hamilton_paths(T)

            for Q in Q_candidates:
                Tv = tournament_minus_v(T, Q)
                H_Tv = hamilton_paths(Tv)
                rivals = [v for v in range(n) if v != Q and T[v][Q]]
                delta = H_T - H_Tv
                bound = 2 * len(rivals)

                if delta == bound:
                    sc = is_strongly_connected(T)
                    tight_cases.append({
                        'score': tuple(sorted(degs)),
                        'sc': sc,
                        'n_rivals': len(rivals),
                        'H_T': H_T,
                        'deg_Q': max_deg
                    })

        print(f"\n  n={n}: {len(tight_cases)} tight cases (delta=2*|rivals|)")
        # Distribution by whether strongly connected
        sc_tight = [c for c in tight_cases if c['sc']]
        nsc_tight = [c for c in tight_cases if not c['sc']]
        print(f"    SC: {len(sc_tight)}, not-SC: {len(nsc_tight)}")
        # Distribution by score sequence
        by_score = defaultdict(int)
        for c in tight_cases:
            by_score[c['score']] += 1
        for ss, cnt in sorted(by_score.items()):
            print(f"    score {ss}: {cnt} times")
