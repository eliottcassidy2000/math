import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
drt_contraction_tournament.py
kind-pasteur-2026-03-07-S39b

For DRT (doubly regular tournaments), every pair (u,v) has exactly
(n-3)/4 common out-neighbors. The contraction T/e is a tournament
iff u and v have IDENTICAL out-neighborhoods.

For DRTs: out-degree = (n-1)/2, common out-nbrs = (n-3)/4.
So they agree on (n-3)/4 out-nbrs and disagree on... let's compute.

If u->v: then among the n-2 other vertices:
  - (n-3)/4 are common out-nbrs of u and v
  - How many does u beat that v loses to?
    u has (n-1)/2 out-nbrs total. One is v itself. So (n-1)/2 - 1 = (n-3)/2
    out-nbrs among the other n-2 vertices.
    v has (n-1)/2 out-nbrs. One might be u (but u->v so v does NOT beat u).
    So v has (n-1)/2 out-nbrs among the other n-2 vertices.

    Common out: (n-3)/4
    u-only out (u beats, v loses): (n-3)/2 - (n-3)/4 = (n-3)/4
    v-only out (v beats, u loses): (n-1)/2 - (n-3)/4 = (2(n-1) - (n-3))/4 = (n+1)/4
    Both lose: n-2 - (n-3)/4 - (n-3)/4 - (n+1)/4 = n-2 - (3n-5)/4 = (n-3)/4

So for DRTs with u->v:
  Common out-nbrs: (n-3)/4
  u-only out: (n-3)/4
  v-only out: (n+1)/4
  Both-in: (n-3)/4

For T/e to be a tournament: need u-only out = 0 AND v-only out = 0.
That means (n-3)/4 = 0 AND (n+1)/4 = 0.
(n-3)/4 = 0 => n = 3.
So for DRTs with n > 3: T/e is NEVER a tournament!

This is expected since the contraction criterion requires identical
out-neighborhoods, which for DRTs means (n-3)/4 = 0.

Let's verify this and explore what happens with NEAR-DRT or general
regular tournaments.

Also: for NON-regular tournaments, some edges may give tournament
contractions. Which ones? How does this relate to H(T)?
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict


def contract_edge(T, u, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    m = len(verts)
    D = [[0]*m for _ in range(m)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            if i == u:
                D[i_idx][j_idx] = T[v][j]
            elif j == u:
                D[i_idx][j_idx] = T[i][u]
            else:
                D[i_idx][j_idx] = T[i][j]
    return D


def is_tournament(adj, n):
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


def count_tc_edges(T):
    """Count edges whose contraction gives a tournament."""
    n = len(T)
    count = 0
    total = 0
    for u in range(n):
        for v in range(n):
            if u == v or not T[u][v]:
                continue
            # Check: do u and v have same out-nbrs on V\{u,v}?
            same = all(T[u][x] == T[v][x] for x in range(n) if x != u and x != v)
            if same:
                count += 1
            total += 1
    return count, total


# ============================================================
# Analysis
# ============================================================
print("=" * 70)
print("TOURNAMENT-CONTRACTION EDGES: which edges give T/e = tournament?")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total_tours = 1 << m

    if n > 6:
        # Sample
        import random
        random.seed(42)
        sample = [random.randint(0, total_tours - 1) for _ in range(5000)]
        sample_type = "sampled (5000)"
    else:
        sample = range(total_tours)
        sample_type = "exhaustive"

    tc_by_score = defaultdict(list)

    for bits in sample:
        T = tournament_from_bits(n, bits)
        scores = tuple(sorted(sum(T[i]) for i in range(n)))
        tc_count, total_edges = count_tc_edges(T)
        tc_by_score[scores].append(tc_count)

    print(f"\nn={n} ({sample_type}):")
    print(f"  {'Score':<25} {'Mean TC edges':>15} {'Min':>5} {'Max':>5} {'Total arcs':>12}")
    for scores in sorted(tc_by_score.keys()):
        vals = tc_by_score[scores]
        mean_tc = sum(vals) / len(vals)
        total_arcs = n * (n - 1) // 2 * 2  # directed edges
        print(f"  {str(scores):<25} {mean_tc:>15.2f} {min(vals):>5} {max(vals):>5} {total_arcs:>12}")

# ============================================================
# H-value correlation with number of TC edges
# ============================================================
print("\n" + "=" * 70)
print("CORRELATION: TC edges vs H(T)")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    tc_h_pairs = []

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        tc_count, _ = count_tc_edges(T)
        H = hamiltonian_path_count(T)
        tc_h_pairs.append((tc_count, H))

    # Compute correlation
    tc_vals = [p[0] for p in tc_h_pairs]
    h_vals = [p[1] for p in tc_h_pairs]
    n_pts = len(tc_h_pairs)
    mean_tc = sum(tc_vals) / n_pts
    mean_h = sum(h_vals) / n_pts
    cov = sum((t - mean_tc) * (h - mean_h) for t, h in tc_h_pairs) / n_pts
    var_tc = sum((t - mean_tc)**2 for t in tc_vals) / n_pts
    var_h = sum((h - mean_h)**2 for h in h_vals) / n_pts

    if var_tc > 0 and var_h > 0:
        corr = cov / (var_tc**0.5 * var_h**0.5)
    else:
        corr = 0

    print(f"\n  n={n}: corr(TC_edges, H) = {corr:.4f}")
    # Group by TC count
    tc_groups = defaultdict(list)
    for tc, h in tc_h_pairs:
        tc_groups[tc].append(h)
    for tc in sorted(tc_groups.keys()):
        hs = tc_groups[tc]
        print(f"    TC={tc}: {len(hs)} tournaments, mean H={sum(hs)/len(hs):.1f}, "
              f"min={min(hs)}, max={max(hs)}")

# ============================================================
# Transitive tournament: which edges are TC?
# ============================================================
print("\n" + "=" * 70)
print("TRANSITIVE TOURNAMENT: TC edges")
print("=" * 70)

for n in range(3, 8):
    # Transitive: T[i][j] = 1 iff i < j
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1

    tc_edges = []
    for u in range(n):
        for v in range(n):
            if u == v or not T[u][v]:
                continue
            same = all(T[u][x] == T[v][x] for x in range(n) if x != u and x != v)
            if same:
                tc_edges.append((u, v))

    print(f"  n={n}: {len(tc_edges)} TC edges: {tc_edges}")
    # For transitive: u < v always. u beats {u+1,...,n-1}, v beats {v+1,...,n-1}.
    # Same out-nbrs on V\{u,v}: for x != u,v, T[u][x] = T[v][x]
    # T[u][x] = 1 iff x > u; T[v][x] = 1 iff x > v (where v > u).
    # For x < u: both 0 (agree)
    # For u < x < v: T[u][x] = 1, T[v][x] = 0 (disagree unless no such x)
    # For x > v: both 1 (agree)
    # So agree iff NO x with u < x < v, i.e., v = u+1 (consecutive).

print("\nDone.")
