#!/usr/bin/env python3
r"""beta2_cocycle_algebra.py - Algebraic structure forcing b1 <= 1

KEY QUESTION: In Z_1 (dim = (n-1)(n-2)/2), TT boundaries give linear
constraints. Why do they always have corank <= 1?

APPROACH: Analyze the TT constraint matrix structure.

For a 1-cycle z in Z_1:
  - TT cocycle condition: z(a,b) + z(b,c) = z(a,c) for each TT (a,b,c)
  - 1-cycle condition: sum_in z(u,v) = sum_out z(v,w) for each vertex v

The TT cocycle condition means z is a "coboundary" on transitive edges.
On 3-cycle edges, z has one free parameter per 3-cycle (the cycle sum).

Question: How many independent cycle sums can there be?

If all cycle sums are forced equal by the structure, then dim <= 1.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_3cycles(A, n):
    """Return list of directed 3-cycles (i,j,k) with i->j->k->i, i<j<k (unordered)."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i, k, j))  # oriented: i->k->j->i
    return cycles


def analyze_cocycle_structure(A, n):
    """Analyze TT cocycle space structure."""
    # Build edge list
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)

    # D1 matrix
    D1 = np.zeros((n, num_edges))
    for idx, (a, b) in enumerate(edges):
        D1[b, idx] += 1
        D1[a, idx] -= 1

    # Z_1 basis
    U, S, Vt = np.linalg.svd(D1, full_matrices=True)
    rk_d1 = int(sum(s > 1e-8 for s in S))
    Z1_basis = Vt[rk_d1:]  # rows are basis vectors of Z_1
    dim_Z1 = Z1_basis.shape[0]

    # TT triples
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not A[b][c]:
                    continue
                if A[a][c]:
                    tt_triples.append((a, b, c))

    # 3-cycles
    cycles_3 = get_3cycles(A, n)
    c3 = len(cycles_3)

    # TT boundary matrix in Z_1
    if not tt_triples:
        return {
            'dim_Z1': dim_Z1, 'num_tt': 0, 'c3': c3,
            'rk_tt': 0, 'b1': dim_Z1
        }

    M_TT = np.zeros((num_edges, len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        M_TT[edge_idx[(b, c)], j] += 1
        M_TT[edge_idx[(a, c)], j] -= 1
        M_TT[edge_idx[(a, b)], j] += 1

    M_TT_Z1 = Z1_basis @ M_TT
    sv = np.linalg.svd(M_TT_Z1, compute_uv=False)
    rk_tt = int(sum(s > 1e-8 for s in sv))
    b1 = dim_Z1 - rk_tt

    # If b1 > 0, find the cocycle
    cocycle = None
    if b1 > 0:
        U2, S2, Vt2 = np.linalg.svd(M_TT_Z1.T, full_matrices=True)
        # Null space of M_TT_Z1.T = cocycle space in Z_1 coords
        null_vecs = Vt2[rk_tt:]  # shape (b1, dim_Z1)
        # Convert to edge coordinates
        cocycle_edges = null_vecs @ Z1_basis  # shape (b1, num_edges)
        cocycle = cocycle_edges[0]  # first cocycle

    # Identify 3-cycle edges
    cycle_edges = set()
    for (i, j, k) in cycles_3:
        # oriented cycle: find the actual edges
        if A[i][j]:
            cycle_edges.add((i, j))
        if A[j][k]:
            cycle_edges.add((j, k))
        if A[k][i]:
            cycle_edges.add((k, i))
        if A[j][i]:
            cycle_edges.add((j, i))
        if A[k][j]:
            cycle_edges.add((k, j))
        if A[i][k]:
            cycle_edges.add((i, k))

    # Non-cycle edges (edges not in any 3-cycle)
    non_cycle_edges = [e for e in edges if e not in cycle_edges]

    # Edge classification: for each edge, how many TT triples use it?
    tt_per_edge = Counter()
    for (a, b, c) in tt_triples:
        tt_per_edge[(a, b)] += 1
        tt_per_edge[(b, c)] += 1
        tt_per_edge[(a, c)] += 1

    # Edges NOT in any TT triple
    unconstrained_edges = [e for e in edges if tt_per_edge[e] == 0]

    return {
        'dim_Z1': dim_Z1, 'num_tt': len(tt_triples), 'c3': c3,
        'rk_tt': rk_tt, 'b1': b1,
        'num_cycle_edges': len(cycle_edges),
        'num_non_cycle_edges': len(non_cycle_edges),
        'num_unconstrained_edges': len(unconstrained_edges),
        'unconstrained_edges': unconstrained_edges,
        'non_cycle_edges': non_cycle_edges,
        'cocycle': cocycle,
        'edges': edges, 'edge_idx': edge_idx,
        'cycles_3': cycles_3
    }


# ============================================================
# Part 1: Edge classification at n=4,5
# ============================================================
print("=" * 70)
print("EDGE CLASSIFICATION")
print("=" * 70)

for n in [4, 5]:
    total = 2 ** (n * (n - 1) // 2)
    results = []

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        r = analyze_cocycle_structure(A, n)
        results.append(r)

    print(f"\nn={n}: {total} tournaments")

    # Summarize by (b1, c3)
    by_b1_c3 = defaultdict(list)
    for r in results:
        by_b1_c3[(r['b1'], r['c3'])].append(r)

    for (b1, c3), rs in sorted(by_b1_c3.items()):
        nc_dist = Counter(r['num_non_cycle_edges'] for r in rs)
        uc_dist = Counter(r['num_unconstrained_edges'] for r in rs)
        print(f"  b1={b1}, c3={c3}: {len(rs)} tournaments")
        print(f"    Non-cycle edges: {dict(sorted(nc_dist.items()))}")
        print(f"    Unconstrained edges: {dict(sorted(uc_dist.items()))}")


# ============================================================
# Part 2: Why is corank <= 1? Dimension counting
# ============================================================
print(f"\n{'=' * 70}")
print("DIMENSION COUNTING: WHY CORANK <= 1")
print("=" * 70)

for n in [4, 5]:
    total = 2 ** (n * (n - 1) // 2)
    print(f"\nn={n}:")
    print(f"  dim(Z_1) = {(n-1)*(n-2)//2}")
    print(f"  num_edges = {n*(n-1)//2}")

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        r = analyze_cocycle_structure(A, n)
        if r['b1'] == 1:
            # Print detailed info for ONE example
            print(f"\n  EXAMPLE b1=1: bits={bits}, c3={r['c3']}")
            print(f"    dim_Z1={r['dim_Z1']}, #TT={r['num_tt']}, rk_TT={r['rk_tt']}")
            print(f"    Non-cycle edges: {r['non_cycle_edges']}")
            print(f"    Unconstrained edges: {r['unconstrained_edges']}")

            # Compute: which edges are in 3-cycles?
            cycle_edge_list = []
            for (i, j, k) in r['cycles_3']:
                ce = []
                if A[i][j]: ce.append((i, j))
                if A[j][k]: ce.append((j, k))
                if A[k][i]: ce.append((k, i))
                if A[j][i]: ce.append((j, i))
                if A[k][j]: ce.append((k, j))
                if A[i][k]: ce.append((i, k))
                cycle_edge_list.append(((i, j, k), ce))
            print(f"    3-cycles and their edges:")
            for cyc, ces in cycle_edge_list:
                print(f"      {cyc}: {ces}")

            # Cocycle values
            if r['cocycle'] is not None:
                z = r['cocycle']
                norm = max(abs(z))
                z_norm = z / norm
                print(f"    Cocycle (normalized):")
                for e_idx, e in enumerate(r['edges']):
                    if abs(z_norm[e_idx]) > 1e-8:
                        print(f"      {e}: {z_norm[e_idx]:.4f}")

                # Cycle sums
                for (i, j, k) in r['cycles_3']:
                    # Find directed edges of cycle
                    if A[i][j] and A[j][k] and A[k][i]:
                        s = z_norm[r['edge_idx'][(i,j)]] + z_norm[r['edge_idx'][(j,k)]] + z_norm[r['edge_idx'][(k,i)]]
                    else:
                        s = z_norm[r['edge_idx'][(i,k)]] + z_norm[r['edge_idx'][(k,j)]] + z_norm[r['edge_idx'][(j,i)]]
                    print(f"    Cycle sum {(i,j,k)}: {s:.4f}")

            break  # Just one example


# ============================================================
# Part 3: Graph of 3-cycles (cycle-adjacency)
# ============================================================
print(f"\n{'=' * 70}")
print("3-CYCLE ADJACENCY GRAPH")
print("=" * 70)

for n in [5, 6]:
    total_bits = 2 ** (n * (n - 1) // 2)
    max_check = total_bits if n <= 5 else 1000
    t0 = time.time()

    cycle_connectivity = Counter()
    all_connected = True

    for trial in range(max_check):
        if n <= 5:
            bits = trial
        else:
            bits = random.randint(0, total_bits - 1)

        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        cycles = get_3cycles(A, n)
        c3 = len(cycles)
        if c3 == 0:
            continue

        # Build adjacency: two 3-cycles adjacent if they share an edge
        adj = [[False]*c3 for _ in range(c3)]
        cycle_edge_sets = []
        for (i, j, k) in cycles:
            es = set()
            if A[i][j]: es.add((i,j))
            if A[j][k]: es.add((j,k))
            if A[k][i]: es.add((k,i))
            if A[j][i]: es.add((j,i))
            if A[k][j]: es.add((k,j))
            if A[i][k]: es.add((i,k))
            cycle_edge_sets.append(es)

        for p in range(c3):
            for q in range(p+1, c3):
                if cycle_edge_sets[p] & cycle_edge_sets[q]:
                    adj[p][q] = adj[q][p] = True

        # Check connectivity
        visited = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in range(c3):
                if adj[u][v] and v not in visited:
                    visited.add(v)
                    stack.append(v)

        connected = (len(visited) == c3)
        cycle_connectivity[(c3, connected)] += 1
        if not connected:
            all_connected = False

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s, {max_check} tournaments):")
    print(f"  Cycle graph connectivity: {dict(sorted(cycle_connectivity.items()))}")
    print(f"  ALL cycle graphs connected: {all_connected}")


# ============================================================
# Part 4: Key insight - TT cocycle = cycle sum function
# ============================================================
print(f"\n{'=' * 70}")
print("TT-COCYCLE = CYCLE-SUM FUNCTION")
print("=" * 70)

# THESIS: When b1=1, the unique TT-cocycle z assigns a constant value
# to each 3-cycle sum. This is because:
# 1. z restricted to TT edges is a coboundary (potential difference)
# 2. The 3-cycle sums are the only "free" parameters
# 3. The cycle graph is connected (edge-sharing), so equal cycle sums
#    propagate through the cycle graph
# 4. Therefore at most 1 DOF: the common cycle sum value

# Let's verify: if two 3-cycles share an edge, do their sums MUST be equal?
# Consider cycles C1 = (a,b,c) with a->b->c->a and C2 sharing edge (a,b).
# C2 must have a->b and a third vertex d.
# If C2 = (a,b,d) with a->b->d->a, then:
#   sum(C1) = z(a,b) + z(b,c) + z(c,a)
#   sum(C2) = z(a,b) + z(b,d) + z(d,a)
#
# For these to be equal: z(b,c) + z(c,a) = z(b,d) + z(d,a)
# i.e., z(b,c) - z(b,d) = z(d,a) - z(c,a)
#
# Now, look at the 4-vertex subtournament on {a,b,c,d}.
# We know: a->b, c->a (from C1), d->a (from C2).
# The edge between c and d: either c->d or d->c.
#
# Case 1: c->d. Then (b,c,d) forms a TT if b->c,c->d, and b->d (?)
#   We need to check what edges exist.
#   From C1: a->b, b->c, c->a
#   From C2: a->b, b->d, d->a
#   Edge c-d: c->d
#   Edge b->c: yes, b->d: yes
#   TT (b,c,d) requires c->d AND b->d: have both! So z(b,c)+z(c,d)=z(b,d)
#   Also check (c,d,a): c->d, d->a, and c->a: yes! TT. z(c,d)+z(d,a)=z(c,a)
#   From TT(b,c,d): z(b,c)-z(b,d) = -z(c,d)
#   From TT(c,d,a): z(d,a)-z(c,a) = -z(c,d)
#   So z(b,c)-z(b,d) = z(d,a)-z(c,a). EQUAL CYCLE SUMS!
#
# Case 2: d->c. Then (b,d,c) if b->d, d->c, b->c: TT! z(b,d)+z(d,c)=z(b,c)
#   Also (d,c,a): d->c, c->a, d->a: TT! z(d,c)+z(c,a)=z(d,a)
#   From TT(b,d,c): z(b,c)-z(b,d) = z(d,c)
#   From TT(d,c,a): z(d,a)-z(c,a) = z(d,c)
#   Again equal!

print("""
THEOREM: If two directed 3-cycles C1 and C2 share a directed edge,
then for any TT-cocycle z: sum_C1(z) = sum_C2(z).

PROOF: Let C1 = (a->b->c->a) and C2 = (a->b->d->a) share edge a->b.
Then c->a and d->a. For edge c-d:

Case 1: c->d.
  TT(b,c,d): b->c->d, b->d => z(b,c)+z(c,d) = z(b,d)
  TT(c,d,a): c->d->a, c->a => z(c,d)+z(d,a) = z(c,a)
  => z(b,c)-z(b,d) = -z(c,d) = z(d,a)-z(c,a)
  => sum_C1 - sum_C2 = [z(b,c)+z(c,a)] - [z(b,d)+z(d,a)]
                      = [z(b,c)-z(b,d)] - [z(d,a)-z(c,a)] = 0

Case 2: d->c.
  TT(b,d,c): b->d->c, b->c => z(b,d)+z(d,c) = z(b,c)
  TT(d,c,a): d->c->a, d->a => z(d,c)+z(c,a) = z(d,a)
  => z(b,c)-z(b,d) = z(d,c) = z(d,a)-z(c,a)
  => sum_C1 - sum_C2 = 0

COROLLARY: If the 3-cycle graph (adjacency by shared edge) is connected,
then ALL cycle sums are equal, and the cocycle space is at most 1-dim.

COROLLARY: b1(T) <= (number of connected components of 3-cycle graph).
When cycle graph is connected: b1 <= 1.
When no 3-cycles: b1 = 0 (tournament is transitive, acyclic => beta_2=0).
""")

# Verify the theorem computationally
print("VERIFICATION: Shared-edge cycles have equal sums")
n = 5
count_verified = 0
count_total = 0
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    r = analyze_cocycle_structure(A, n)
    if r['b1'] == 0 or r['cocycle'] is None:
        continue

    z = r['cocycle']
    cycles = r['cycles_3']
    edge_idx = r['edge_idx']

    # Compute cycle sums
    cycle_sums = []
    for (i, j, k) in cycles:
        if A[i][j] and A[j][k] and A[k][i]:
            s = z[edge_idx[(i,j)]] + z[edge_idx[(j,k)]] + z[edge_idx[(k,i)]]
        else:
            s = z[edge_idx[(i,k)]] + z[edge_idx[(k,j)]] + z[edge_idx[(j,i)]]
        cycle_sums.append(s)

    # Check: all cycle sums equal?
    if cycle_sums:
        ref = cycle_sums[0]
        all_equal = all(abs(s - ref) < 1e-6 for s in cycle_sums)
        count_total += 1
        if all_equal:
            count_verified += 1
        else:
            print(f"  COUNTEREXAMPLE at bits={bits}: sums={cycle_sums}")

print(f"  n=5: {count_verified}/{count_total} b1=1 tournaments have all cycle sums equal")


# ============================================================
# Part 5: When is the 3-cycle graph connected?
# ============================================================
print(f"\n{'=' * 70}")
print("3-CYCLE GRAPH CONNECTIVITY vs b1")
print("=" * 70)

for n in [4, 5]:
    total = 2 ** (n*(n-1)//2)
    components_vs_b1 = Counter()

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        cycles = get_3cycles(A, n)
        c3 = len(cycles)

        if c3 == 0:
            num_components = 0
        else:
            # Build edge-sharing adjacency
            cycle_edge_sets = []
            for (i, j, k) in cycles:
                es = set()
                if A[i][j]: es.add((i,j))
                if A[j][k]: es.add((j,k))
                if A[k][i]: es.add((k,i))
                if A[j][i]: es.add((j,i))
                if A[k][j]: es.add((k,j))
                if A[i][k]: es.add((i,k))
                cycle_edge_sets.append(es)

            adj = [[False]*c3 for _ in range(c3)]
            for p in range(c3):
                for q in range(p+1, c3):
                    if cycle_edge_sets[p] & cycle_edge_sets[q]:
                        adj[p][q] = adj[q][p] = True

            # Count components
            visited = set()
            num_components = 0
            for start in range(c3):
                if start in visited:
                    continue
                num_components += 1
                stack = [start]
                visited.add(start)
                while stack:
                    u = stack.pop()
                    for v in range(c3):
                        if adj[u][v] and v not in visited:
                            visited.add(v)
                            stack.append(v)

        r = analyze_cocycle_structure(A, n)
        components_vs_b1[(c3, num_components, r['b1'])] += 1

    print(f"\nn={n}: (c3, #components, b1) => count")
    for key, cnt in sorted(components_vs_b1.items()):
        print(f"  c3={key[0]}, components={key[1]}, b1={key[2]}: {cnt}")


# ============================================================
# Part 6: At n=6 sampled
# ============================================================
print(f"\n{'=' * 70}")
print("n=6 SAMPLED: CYCLE GRAPH CONNECTIVITY")
print("=" * 70)

random.seed(42)
n = 6
trials = 500
components_vs_b1 = Counter()

for trial in range(trials):
    if trial < 2**(n*(n-1)//2) and trial < 500:
        A = [[0] * n for _ in range(n)]
        idx = 0
        bits = random.randint(0, 2**(n*(n-1)//2)-1)
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
    else:
        A = random_tournament(n)

    cycles = get_3cycles(A, n)
    c3 = len(cycles)

    if c3 == 0:
        num_components = 0
    else:
        cycle_edge_sets = []
        for (i, j, k) in cycles:
            es = set()
            if A[i][j]: es.add((i,j))
            if A[j][k]: es.add((j,k))
            if A[k][i]: es.add((k,i))
            if A[j][i]: es.add((j,i))
            if A[k][j]: es.add((k,j))
            if A[i][k]: es.add((i,k))
            cycle_edge_sets.append(es)

        adj = [[False]*c3 for _ in range(c3)]
        for p in range(c3):
            for q in range(p+1, c3):
                if cycle_edge_sets[p] & cycle_edge_sets[q]:
                    adj[p][q] = adj[q][p] = True

        visited = set()
        num_components = 0
        for start in range(c3):
            if start in visited:
                continue
            num_components += 1
            stack = [start]
            visited.add(start)
            while stack:
                u = stack.pop()
                for v in range(c3):
                    if adj[u][v] and v not in visited:
                        visited.add(v)
                        stack.append(v)

    r = analyze_cocycle_structure(A, n)
    components_vs_b1[(num_components, r['b1'])] += 1

print(f"n=6: (#components, b1) => count ({trials} trials)")
for key, cnt in sorted(components_vs_b1.items()):
    print(f"  components={key[0]}, b1={key[1]}: {cnt}")


# ============================================================
# Part 7: Can we prove connectivity?
# ============================================================
print(f"\n{'=' * 70}")
print("CONNECTIVITY ANALYSIS")
print("=" * 70)

# For any tournament with c3 >= 1:
# Two 3-cycles share a directed edge iff they share 2 vertices
# and the common edge is oriented the same way in both cycles.
# But since tournaments have exactly one orientation per vertex pair,
# sharing 2 vertices means sharing exactly 1 directed edge.

# So 3-cycle graph connectivity = vertex-pair-sharing connectivity of cycles
# = the "cycle graph" where cycles adjacent iff they share 2 of 3 vertices.

# When is this NOT connected?
# Need two sets of 3-cycles with no shared vertex pair.
# At n=4: only 2 possible 3-cycles (disjoint vertex sets impossible on 4 vertices)
# So any two 3-cycles at n=4 share at least 2 vertices => connected.

# At n=5 with c3=5: every pair of 3 vertices from {0,1,2,3,4} shares 1 vertex.
# So any 3-cycle on {a,b,c} and another on {d,e,f} share at least 1 vertex...
# Actually C(5,3)=10 but 3-sets can share 0,1,2 vertices.
# {0,1,2} and {3,4,?}: need 3 vertices from remaining, but only 2 left.
# So at n=5 any two 3-element subsets share at least 1 vertex.
# But sharing 1 vertex doesn't mean sharing a directed edge...

# Let me check: at n=5, when is the cycle graph disconnected?
n = 5
disc_count = 0
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    cycles = get_3cycles(A, n)
    c3 = len(cycles)
    if c3 <= 1:
        continue

    cycle_edge_sets = []
    for (i, j, k) in cycles:
        es = set()
        if A[i][j]: es.add((i,j))
        if A[j][k]: es.add((j,k))
        if A[k][i]: es.add((k,i))
        if A[j][i]: es.add((j,i))
        if A[k][j]: es.add((k,j))
        if A[i][k]: es.add((i,k))
        cycle_edge_sets.append(es)

    adj_list = [[] for _ in range(c3)]
    for p in range(c3):
        for q in range(p+1, c3):
            if cycle_edge_sets[p] & cycle_edge_sets[q]:
                adj_list[p].append(q)
                adj_list[q].append(p)

    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj_list[u]:
            if v not in visited:
                visited.add(v)
                stack.append(v)

    if len(visited) < c3:
        disc_count += 1
        if disc_count <= 3:
            print(f"  DISCONNECTED at n=5, bits={bits}, c3={c3}")
            r = analyze_cocycle_structure(A, n)
            print(f"    b1={r['b1']}")

print(f"  n=5: {disc_count} tournaments with disconnected cycle graph (c3>=2)")


# At larger n
print(f"\n  Larger n (sampled):")
for n in [7, 8, 9]:
    random.seed(42)
    disc = 0
    trials = 500
    for trial in range(trials):
        A = random_tournament(n)
        cycles = get_3cycles(A, n)
        c3 = len(cycles)
        if c3 <= 1:
            continue

        cycle_edge_sets = []
        for (i, j, k) in cycles:
            es = set()
            for x, y in [(i,j),(j,k),(k,i),(j,i),(k,j),(i,k)]:
                if A[x][y]:
                    es.add((x,y))
            cycle_edge_sets.append(es)

        adj_list = [[] for _ in range(c3)]
        for p in range(c3):
            for q in range(p+1, c3):
                if cycle_edge_sets[p] & cycle_edge_sets[q]:
                    adj_list[p].append(q)
                    adj_list[q].append(p)

        visited = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in adj_list[u]:
                if v not in visited:
                    visited.add(v)
                    stack.append(v)

        if len(visited) < c3:
            disc += 1

    print(f"  n={n}: {disc}/{trials} disconnected")


print("\n\nDone.")
