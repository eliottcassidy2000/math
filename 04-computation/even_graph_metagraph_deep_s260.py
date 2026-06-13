#!/usr/bin/env python3
"""
even_graph_metagraph_deep_s260.py — Deep analysis of even graphs through the metagraph lens
opus-2026-03-23-S260

Key questions:
1. Even graph flip graph G_n^{even}: what are its properties vs tournament G_n?
2. The fiber structure: how many tournaments project to each even graph class?
3. The "even shadow" of the merged metagraph: what does complement pairing look like?
4. The cut/cycle decomposition: is each tournament class determined by (score_class, even_class)?
5. Cross-space bipartite graph: tournament classes ↔ even graph classes
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb, factorial
import subprocess
import sys

def get_tournaments(n):
    """Get all tournament iso class reps via nauty's gentourng."""
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

def adj_to_edge_set(adj, n):
    """Convert adjacency matrix to frozenset of undirected edges (for even graph)."""
    edges = set()
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]: edges.add((i,j))
    return frozenset(edges)

def tournament_to_gf2_vector(adj, n, P):
    """Tournament adjacency matrix -> GF(2) bit vector (upper triangle)."""
    return tuple(adj[i][j] for (i,j) in P)

def gf2_to_undirected(vec, n, P):
    """GF(2) vector -> undirected graph (even graph = cycle space projection)."""
    # For tournament T, the underlying undirected graph is K_n (all edges present).
    # The EVEN GRAPH associated to T is its cycle-space projection.
    # Over GF(2): cut space spanned by delta(v) for v in V.
    # delta(v) = characteristic vector of edges incident to v.
    # Cycle space = orthogonal complement of cut space.
    # Project T's bit vector onto cycle space.

    m = len(P)
    # Build cut space basis: delta(v) for v = 0,...,n-2 (n-1 vectors, but rank n-1)
    cut_basis = []
    for v in range(n-1):  # n-1 generators suffice
        vec_v = [0]*m
        for k, (i,j) in enumerate(P):
            if i == v or j == v:
                vec_v[k] = 1
        cut_basis.append(vec_v)

    # Gaussian elimination to get row echelon form
    basis = [list(b) for b in cut_basis]
    pivots = []
    for col_target in range(len(basis)):
        # Find pivot
        found = False
        for r in range(len(pivots), len(basis)):
            if basis[r][col_target if col_target < m else 0]:
                pass  # continue search
        # Actually do proper GE

    # Simpler: project directly. The cycle-space projection of v is v - proj_cut(v).
    # In GF(2), projection onto cut space: express v in terms of cut basis.
    # But it's easier to just compute the cycle space component.

    # Method: the cycle space projection = v XOR (sum of cut basis vectors that make v orthogonal to cut space)
    # Actually in GF(2): v = v_cut + v_cycle. We need v_cycle = v + v_cut.
    # v_cut is the unique element of cut space closest to v in Hamming distance.

    # For tournaments: the cut-space component is determined by the score sequence.
    # Specifically: for edge (i,j) with i<j, the cut component bit equals:
    #   sum_{v != i,j} delta(v)_{(i,j)} * c_v mod 2 where c_v are the coefficients
    # This is getting complicated. Let me use a direct approach.

    # Direct GF(2) linear algebra:
    # Build the full cut space (as a set of vectors)
    # Then for each tournament vector, find its cycle-space projection

    # Build cut space via GF(2) row reduction
    mat = [list(b) for b in cut_basis]
    pivot_cols = []
    row = 0
    for col in range(m):
        # Find a row with a 1 in this column at or below current row
        found = -1
        for r in range(row, len(mat)):
            if mat[r][col]:
                found = r
                break
        if found == -1:
            continue
        # Swap
        mat[row], mat[found] = mat[found], mat[row]
        # Eliminate
        for r in range(len(mat)):
            if r != row and mat[r][col]:
                mat[r] = [(mat[r][k] + mat[row][k]) % 2 for k in range(m)]
        pivot_cols.append(col)
        row += 1

    # Now mat[:row] is RREF of cut space. rank = row = n-1.
    # Project v onto cut space: for each basis vector, if v has a 1 in the pivot column, XOR.
    v_list = list(vec)
    v_cut = [0]*m
    for r_idx in range(row):
        pc = pivot_cols[r_idx]
        if v_list[pc]:
            v_cut = [(v_cut[k] + mat[r_idx][k]) % 2 for k in range(m)]
            v_list = [(v_list[k] + mat[r_idx][k]) % 2 for k in range(m)]

    # v_cycle = v XOR v_cut
    v_cycle = tuple((vec[k] + v_cut[k]) % 2 for k in range(m))
    return v_cycle, tuple(v_cut)

def canonical_graph(adj_or_edges, n):
    """Canonical form of an undirected graph under S_n."""
    if isinstance(adj_or_edges, (list, tuple)) and len(adj_or_edges) == n:
        # It's an adjacency matrix
        adj = adj_or_edges
    else:
        # It's a frozenset of edges
        adj = [[0]*n for _ in range(n)]
        for (i,j) in adj_or_edges:
            adj[i][j] = 1
            adj[j][i] = 1

    best = None
    # For small n, enumerate all permutations
    if n <= 6:
        for perm in permutations(range(n)):
            sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
            if best is None or sig < best:
                best = sig
        return best
    else:
        # For n=7,8: use score-based grouping to reduce search
        scores = sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n))
        score_key = tuple(scores)

        # Group vertices by degree
        deg = [sum(adj[i][j] for j in range(n) if j != i) for i in range(n)]
        groups = defaultdict(list)
        for v in range(n):
            groups[deg[v]].append(v)
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
            if count > 100000:  # safety limit
                break
        return best

def canonical_tournament(adj, n):
    """Canonical form of a tournament."""
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
        if count > 100000:
            break
    return best

def hamiltonian_paths(adj, n):
    """Count Hamiltonian paths in tournament."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u), u)] = dp.get((S | (1 << u), u), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def analyze_n(n, verbose=True):
    """Full analysis at given n."""
    if verbose:
        print(f"\n{'#'*70}")
        print(f"  n = {n}")
        print(f"{'#'*70}")

    tournaments, P = get_tournaments(n)
    m = comb(n, 2)

    if verbose:
        print(f"  Tournament iso classes: {len(tournaments)}")

    # For each tournament, compute:
    # 1. Its GF(2) vector
    # 2. Its cycle-space projection (even graph component)
    # 3. Its cut-space projection (score component)
    # 4. The canonical form of the even graph
    # 5. H(T)

    t_data = []
    even_classes = {}  # canonical_even -> class_id
    even_id = 0

    for tidx, adj in enumerate(tournaments):
        vec = tournament_to_gf2_vector(adj, n, P)
        v_cycle, v_cut = gf2_to_undirected(vec, n, P)

        # Convert cycle projection to undirected graph
        even_edges = set()
        for k, (i,j) in enumerate(P):
            if v_cycle[k]:
                even_edges.add((i,j))

        # Canonical form of the even graph
        even_adj = [[0]*n for _ in range(n)]
        for (i,j) in even_edges:
            even_adj[i][j] = 1
            even_adj[j][i] = 1

        # Check even: all degrees even
        degs = [sum(even_adj[i][j] for j in range(n) if j != i) for i in range(n)]
        is_even = all(d % 2 == 0 for d in degs)

        canon_even = canonical_graph(even_adj, n)
        if canon_even not in even_classes:
            even_classes[canon_even] = even_id
            even_id += 1

        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
        h = hamiltonian_paths(adj, n)

        # Complement tournament
        comp_adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    comp_adj[i][j] = 1 - adj[i][j]
        comp_canon = canonical_tournament(comp_adj, n)
        tourn_canon = canonical_tournament(adj, n)
        is_sc = (tourn_canon == comp_canon)

        t_data.append({
            'tidx': tidx,
            'vec': vec,
            'v_cycle': v_cycle,
            'v_cut': v_cut,
            'even_class': even_classes[canon_even],
            'scores': scores,
            'H': h,
            'is_even': is_even,
            'is_sc': is_sc,
            'canon_even': canon_even,
            'canon_tourn': tourn_canon,
            'n_even_edges': sum(v_cycle),
            'degs': tuple(degs),
        })

    n_even_classes = len(even_classes)

    if verbose:
        print(f"  Even graph iso classes (from tournament projections): {n_even_classes}")
        print(f"  Ratio V_tournament / V_even = {len(tournaments)} / {n_even_classes} = {len(tournaments)/n_even_classes:.3f}")

    # Fiber structure: how many tournaments project to each even graph class?
    fiber_sizes = Counter()
    for td in t_data:
        fiber_sizes[td['even_class']] += 1

    if verbose:
        print(f"\n  FIBER STRUCTURE (tournament classes per even graph class):")
        for eid in sorted(fiber_sizes.keys()):
            tids = [td['tidx'] for td in t_data if td['even_class'] == eid]
            h_vals = [t_data[t]['H'] for t in tids]
            scores_vals = [t_data[t]['scores'] for t in tids]
            sc_vals = [t_data[t]['is_sc'] for t in tids]
            n_edges = t_data[tids[0]]['n_even_edges']
            print(f"    Even class {eid} ({n_edges} edges, degs={t_data[tids[0]]['degs']}): "
                  f"{len(tids)} tournaments, H={h_vals}, SC={sc_vals}")
            if len(set(scores_vals)) > 1:
                print(f"      Score diversity: {len(set(scores_vals))} distinct scores")

    # Is the map tournament_class -> even_class injective?
    is_injective = all(v == 1 for v in fiber_sizes.values())
    if verbose:
        print(f"\n  Is tournament -> even_class INJECTIVE? {is_injective}")
        print(f"  Fiber size distribution: {sorted(fiber_sizes.values())}")

    # Cross-space bipartite graph: tournament classes ↔ even graph classes
    # Build the even graph flip graph (edge flips in even graphs)
    # An even graph edge flip = toggle one edge of K_n in the cycle space projection
    # But wait: a single edge toggle doesn't preserve even-ness.
    # To stay in the cycle space, we need to flip a CYCLE (fundamental circuit).
    # The natural operation is: flip all edges of a triangle (3-edge cycle).
    # Or: the even graph flip graph uses symmetric difference with small even graphs.

    # Actually, the even graph metagraph as defined by kind-pasteur uses
    # single edge flips among ALL graphs, then restricts to even-even pairs.
    # But a single edge flip of an even graph gives an ODD graph (parity changes).
    # So the EE sub-metagraph uses paths of length 2 through odd graphs.

    # Let me instead compute the TOURNAMENT METAGRAPH edges and see how they
    # project onto even graph classes.

    if verbose:
        print(f"\n  TOURNAMENT METAGRAPH AND EVEN GRAPH PROJECTION:")

    # Build tournament metagraph
    tourn_edges = set()
    tourn_to_canon = {}
    for td in t_data:
        tourn_to_canon[td['canon_tourn']] = td['tidx']

    edge_projections = defaultdict(int)  # (even_a, even_b) -> count
    same_even_edges = 0
    diff_even_edges = 0

    for td in t_data:
        adj = tournaments[td['tidx']]
        for p_idx in range(m):
            i, j = P[p_idx]
            # Flip arc (i,j) <-> (j,i)
            new_adj = [row[:] for row in adj]
            new_adj[i][j] = 1 - adj[i][j]
            new_adj[j][i] = 1 - adj[j][i]

            new_canon = canonical_tournament(new_adj, n)
            if new_canon in tourn_to_canon:
                nb = tourn_to_canon[new_canon]
                if nb != td['tidx']:
                    a, b = min(td['tidx'], nb), max(td['tidx'], nb)
                    tourn_edges.add((a, b))

                    ea = t_data[a]['even_class']
                    eb = t_data[b]['even_class']
                    if ea == eb:
                        same_even_edges += 1
                    else:
                        diff_even_edges += 1
                        edge_projections[(min(ea,eb), max(ea,eb))] += 1

    if verbose:
        print(f"  Tournament metagraph edges: {len(tourn_edges)}")
        print(f"  Edges where both endpoints have SAME even class: {same_even_edges//2}")
        print(f"  Edges where endpoints have DIFFERENT even class: {diff_even_edges//2}")
        print(f"  Fraction same-even: {same_even_edges/(same_even_edges+diff_even_edges):.3f}"
              if (same_even_edges+diff_even_edges) > 0 else "")

    # The "even graph shadow" of the tournament metagraph
    # = project tournament metagraph edges onto even graph classes
    even_shadow_edges = set()
    for (ea, eb), cnt in edge_projections.items():
        if ea != eb:
            even_shadow_edges.add((ea, eb))

    if verbose:
        print(f"\n  EVEN GRAPH SHADOW (projected metagraph):")
        print(f"    Vertices: {n_even_classes}")
        print(f"    Edges: {len(even_shadow_edges)}")
        if n_even_classes > 0:
            print(f"    Edge density: {2*len(even_shadow_edges)/(n_even_classes*(n_even_classes-1)):.3f}"
                  if n_even_classes > 1 else "    (trivial)")

    # SC structure in even graph classes
    sc_even_classes = set()
    for td in t_data:
        if td['is_sc']:
            sc_even_classes.add(td['even_class'])

    if verbose:
        print(f"\n  SC TOURNAMENTS AND EVEN GRAPH CLASSES:")
        print(f"    SC tournament classes: {sum(1 for td in t_data if td['is_sc'])}")
        print(f"    Even classes containing SC tournaments: {len(sc_even_classes)}")
        print(f"    Even classes with ONLY SC tournaments: ", end="")
        only_sc = sum(1 for eid in sc_even_classes
                      if all(t_data[t]['is_sc'] for t in range(len(t_data))
                             if t_data[t]['even_class'] == eid))
        print(only_sc)

    # How does arc flip change the cycle-space projection?
    if verbose:
        print(f"\n  ARC FLIP EFFECT ON CYCLE SPACE:")

    hamming_dists_cycle = []
    hamming_dists_cut = []
    for (a, b) in list(tourn_edges)[:50]:  # sample
        va = t_data[a]['v_cycle']
        vb = t_data[b]['v_cycle']
        hd_c = sum(x != y for x, y in zip(va, vb))
        hamming_dists_cycle.append(hd_c)

        va_cut = t_data[a]['v_cut']
        vb_cut = t_data[b]['v_cut']
        hd_k = sum(x != y for x, y in zip(va_cut, vb_cut))
        hamming_dists_cut.append(hd_k)

    if verbose and hamming_dists_cycle:
        print(f"    Cycle-space Hamming distance per arc flip: {sorted(Counter(hamming_dists_cycle).items())}")
        print(f"    Cut-space Hamming distance per arc flip: {sorted(Counter(hamming_dists_cut).items())}")

    # The merged metagraph view
    if verbose:
        print(f"\n  MERGED METAGRAPH VIEW:")

    # Merge complement pairs
    merge_map = {}
    mid = 0
    for td in t_data:
        if td['tidx'] in merge_map:
            continue
        # Find complement
        comp_adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    comp_adj[i][j] = 1 - tournaments[td['tidx']][i][j]
        comp_canon = canonical_tournament(comp_adj, n)
        comp_idx = tourn_to_canon.get(comp_canon, -1)

        if td['is_sc']:
            merge_map[td['tidx']] = mid
            mid += 1
        else:
            merge_map[td['tidx']] = mid
            if comp_idx >= 0 and comp_idx != td['tidx']:
                merge_map[comp_idx] = mid
            mid += 1

    # Merged edges
    merged_edges = set()
    for (a, b) in tourn_edges:
        ma, mb = merge_map[a], merge_map[b]
        if ma != mb:
            merged_edges.add((min(ma, mb), max(ma, mb)))

    # Classify merged edges by SC status
    merged_info = {}
    for mm in range(mid):
        tids = [td['tidx'] for td in t_data if merge_map[td['tidx']] == mm]
        is_sc = t_data[tids[0]]['is_sc']
        even_cls = t_data[tids[0]]['even_class']
        h_val = t_data[tids[0]]['H']
        merged_info[mm] = {'sc': is_sc, 'even_class': even_cls, 'H': h_val}

    spine = sum(1 for (a,b) in merged_edges if merged_info[a]['sc'] and merged_info[b]['sc'])
    ribs = sum(1 for (a,b) in merged_edges if merged_info[a]['sc'] != merged_info[b]['sc'])
    sea = sum(1 for (a,b) in merged_edges if not merged_info[a]['sc'] and not merged_info[b]['sc'])

    if verbose:
        print(f"    V_merged: {mid}")
        print(f"    E_merged: {len(merged_edges)}")
        print(f"    Spine (SC-SC): {spine}")
        print(f"    Ribs (SC-NS): {ribs}")
        print(f"    Sea (NS-NS): {sea}")

    # Even graph class distribution in merged metagraph
    even_in_merged = defaultdict(list)
    for mm in range(mid):
        even_in_merged[merged_info[mm]['even_class']].append(mm)

    if verbose:
        print(f"\n  EVEN GRAPH CLASSES IN MERGED METAGRAPH:")
        for eid in sorted(even_in_merged.keys()):
            mids = even_in_merged[eid]
            sc_count = sum(1 for mm in mids if merged_info[mm]['sc'])
            print(f"    Even class {eid}: {len(mids)} merged nodes ({sc_count} SC)")

    # The key question: does the even graph class partition respect the metagraph structure?
    # How many merged edges connect nodes in the SAME even class vs DIFFERENT?
    same_even_merged = 0
    diff_even_merged = 0
    for (a, b) in merged_edges:
        if merged_info[a]['even_class'] == merged_info[b]['even_class']:
            same_even_merged += 1
        else:
            diff_even_merged += 1

    if verbose:
        total_me = same_even_merged + diff_even_merged
        print(f"\n  METAGRAPH EDGES BY EVEN CLASS:")
        print(f"    Same even class: {same_even_merged} ({100*same_even_merged/total_me:.1f}%)" if total_me > 0 else "")
        print(f"    Different even class: {diff_even_merged} ({100*diff_even_merged/total_me:.1f}%)" if total_me > 0 else "")

    return {
        'n': n,
        'V_tourn': len(tournaments),
        'V_even': n_even_classes,
        'E_tourn': len(tourn_edges),
        'E_shadow': len(even_shadow_edges),
        'fiber_sizes': sorted(fiber_sizes.values()),
        'injective': is_injective,
        'V_merged': mid,
        'E_merged': len(merged_edges),
        'spine': spine,
        'ribs': ribs,
        'sea': sea,
        'same_even_merged': same_even_merged,
        'diff_even_merged': diff_even_merged,
    }

if __name__ == '__main__':
    print("=" * 70)
    print("  EVEN GRAPHS THROUGH THE MERGED METAGRAPH LENS")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    results = []
    for n in range(3, 8):  # n=3..7 (n=8 would be slow with this approach)
        try:
            r = analyze_n(n)
            results.append(r)
        except Exception as e:
            print(f"\n  ERROR at n={n}: {e}")
            import traceback
            traceback.print_exc()

    # Summary table
    print("\n" + "=" * 70)
    print("  SUMMARY TABLE")
    print("=" * 70)
    print(f"{'n':>3} {'V_t':>6} {'V_e':>6} {'ratio':>6} {'E_t':>6} {'E_sha':>6} {'inj':>5} {'V_m':>6} {'E_m':>6} {'same%':>6}")
    for r in results:
        total = r['same_even_merged'] + r['diff_even_merged']
        pct = 100*r['same_even_merged']/total if total > 0 else 0
        print(f"{r['n']:>3} {r['V_tourn']:>6} {r['V_even']:>6} {r['V_tourn']/r['V_even']:.3f} "
              f"{r['E_tourn']:>6} {r['E_shadow']:>6} {str(r['injective']):>5} "
              f"{r['V_merged']:>6} {r['E_merged']:>6} {pct:>5.1f}%")

    print(f"\n  Fiber sizes per n:")
    for r in results:
        print(f"    n={r['n']}: {r['fiber_sizes']}")

    print("\nDONE.")
