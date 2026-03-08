import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
deletion_contraction_test.py
Test the Mitrovic deletion-contraction identity for tournaments.

For a tournament T on n vertices with directed edge e=(u,v):
  - T\\e  = digraph with edge e removed (still n vertices)
  - T/e  = contraction: merge u,v into w (n-1 vertices)
            edges: for x != u,v:  (x,w) iff (x,u),  (w,x) iff (v,x)
            edges among others: unchanged

We test whether H(T) = H(T\\e) + H(T/e) or some other linear relation,
where H(D) = number of Hamiltonian paths in digraph D.

RESULT: H(T) = H(T\\e) + H(T/e) holds UNIVERSALLY (100% of all edges,
all tournaments, at n=4 and n=5). The contraction convention that works
is: w inherits IN-edges from u (the tail), OUT-edges from v (the head).

Uses bitmask DP for Hamiltonian path counting.
Exhaustive at n=4 (64 tournaments) and n=5 (1024 tournaments).
"""

from itertools import permutations
from collections import defaultdict


def tournament_from_bits(n, bits):
    """
    Encode a tournament on n vertices as an adjacency matrix.
    bits encodes the upper triangle: for i<j, bit index = i*(2n-i-1)//2 + (j-i-1).
    If bit is 1: edge i->j; if 0: edge j->i.
    """
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj


def ham_paths_dp(adj, n):
    """
    Count Hamiltonian paths in a digraph on n vertices using bitmask DP.
    dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v.
    """
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def get_edges(adj, n):
    """Return list of directed edges (u,v) where u->v."""
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and adj[i][j]:
                edges.append((i, j))
    return edges


def delete_edge(adj, n, u, v):
    """Return adjacency matrix of T\\e where e=(u,v). Same n vertices, edge removed."""
    new_adj = [row[:] for row in adj]
    new_adj[u][v] = False
    return new_adj


def contract_edge(adj, n, u, v):
    """
    Return adjacency matrix of T/e where e=(u,v), on n-1 vertices.
    Merge u and v into a single vertex w.
    Convention: w inherits IN-edges from u, OUT-edges from v.
      For x != u,v:
        (x, w) iff (x, u) in T
        (w, x) iff (v, x) in T
    Edges among other vertices: unchanged.

    Vertex relabeling: vertices {0,...,n-1} minus {u,v} plus {w}.
    We relabel: keep vertices < min(u,v) as is, shift others down, w gets index of min(u,v).
    """
    # Map old vertices to new indices
    # w will take the position of the smaller of u,v
    w_new = 0  # new index for merged vertex
    old_to_new = {}
    new_n = n - 1

    # Vertices other than u and v, sorted
    others = sorted([x for x in range(n) if x != u and x != v])

    # Assign w_new = 0, then others get 1, 2, ...
    # Actually let's just do a clean mapping
    new_idx = 0
    # Put w first
    w_pos = new_idx
    new_idx += 1
    other_map = {}
    for x in others:
        other_map[x] = new_idx
        new_idx += 1

    new_adj = [[False]*new_n for _ in range(new_n)]

    # Edges among others: unchanged
    for x in others:
        for y in others:
            if adj[x][y]:
                new_adj[other_map[x]][other_map[y]] = True

    # Edges involving w:
    for x in others:
        # (x, w) iff (x, u) in T
        if adj[x][u]:
            new_adj[other_map[x]][w_pos] = True
        # (w, x) iff (v, x) in T
        if adj[v][x]:
            new_adj[w_pos][other_map[x]] = True

    return new_adj, new_n


def test_deletion_contraction(n, verbose=False):
    """Test the deletion-contraction relation for all tournaments on n vertices."""
    num_edges_total = n * (n - 1) // 2
    num_tournaments = 1 << num_edges_total

    # Track relation: H(T) vs H(T\e) + H(T/e)
    relation_counts = defaultdict(int)
    total_tests = 0
    sum_relation_holds = 0

    # Track differences
    diffs = defaultdict(int)

    # For discovering other relations
    all_triples = []  # (H_T, H_del, H_con) for analysis

    print(f"\n{'='*70}")
    print(f"Testing n={n}: {num_tournaments} tournaments, {n*(n-1)//2} edges each")
    print(f"{'='*70}")

    for bits in range(num_tournaments):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)

        edges = get_edges(adj, n)

        for (u, v) in edges:
            # Only test each undirected edge once (the directed version present in T)
            # edges already gives directed edges, so we test each one

            # Deletion
            del_adj = delete_edge(adj, n, u, v)
            H_del = ham_paths_dp(del_adj, n)

            # Contraction
            con_adj, con_n = contract_edge(adj, n, u, v)
            H_con = ham_paths_dp(con_adj, con_n)

            total_tests += 1

            diff = H_T - H_del - H_con
            diffs[diff] += 1

            if H_T == H_del + H_con:
                sum_relation_holds += 1

            all_triples.append((H_T, H_del, H_con, n, bits, u, v))

            if verbose and bits < 4:
                print(f"  T=bits({bits:0{num_edges_total}b}), e=({u},{v}): "
                      f"H(T)={H_T}, H(T\\e)={H_del}, H(T/e)={H_con}, "
                      f"diff={diff}")

    print(f"\nTotal edge tests: {total_tests}")
    print(f"H(T) = H(T\\e) + H(T/e) holds: {sum_relation_holds}/{total_tests} "
          f"({100*sum_relation_holds/total_tests:.1f}%)")

    print(f"\nDifference distribution: H(T) - H(T\\e) - H(T/e)")
    for d in sorted(diffs.keys()):
        print(f"  diff = {d:+d}: {diffs[d]} cases ({100*diffs[d]/total_tests:.1f}%)")

    # Check other linear relations: H(T) = a*H(T\e) + b*H(T/e)?
    # Try to find a,b by least squares or exact check
    print(f"\n--- Checking other linear relations ---")

    # Check H(T) = H(T\e) + 2*H(T/e)
    count_rel2 = sum(1 for (ht, hd, hc, *_) in all_triples if ht == hd + 2*hc)
    print(f"H(T) = H(T\\e) + 2*H(T/e): {count_rel2}/{total_tests}")

    # Check H(T) = H(T\e) - H(T/e)
    count_rel3 = sum(1 for (ht, hd, hc, *_) in all_triples if ht == hd - hc)
    print(f"H(T) = H(T\\e) - H(T/e): {count_rel3}/{total_tests}")

    # Check 2*H(T) = 2*H(T\e) + H(T/e)
    count_rel4 = sum(1 for (ht, hd, hc, *_) in all_triples if 2*ht == 2*hd + hc)
    print(f"2*H(T) = 2*H(T\\e) + H(T/e): {count_rel4}/{total_tests}")

    # Check H(T) + H(T/e) = H(T\e) + 2*H(T/e), i.e., H(T) = H(T\e) + H(T/e) again
    # Try ratio analysis
    print(f"\n--- Ratio analysis ---")
    ratios_del = defaultdict(int)
    ratios_con = defaultdict(int)
    for (ht, hd, hc, *_) in all_triples:
        if hd > 0:
            r = ht / hd
            ratios_del[f"{r:.4f}"] += 1
        if hc > 0:
            r = ht / hc
            ratios_con[f"{r:.4f}"] += 1

    # Show top 5 most common ratios
    print("Top H(T)/H(T\\e) ratios:")
    for r, c in sorted(ratios_del.items(), key=lambda x: -x[1])[:5]:
        print(f"  ratio={r}: {c} cases")

    print("Top H(T)/H(T/e) ratios:")
    for r, c in sorted(ratios_con.items(), key=lambda x: -x[1])[:5]:
        print(f"  ratio={r}: {c} cases")

    # New idea: maybe H(T) = H(T\e) + H(T/e) works only for specific edge types?
    # Check by edge "type" (u->v where score(u) and score(v) matter)
    print(f"\n--- Analysis by edge direction (u->v, score of u and v) ---")
    score_analysis = defaultdict(lambda: defaultdict(int))
    for (ht, hd, hc, nn, bb, uu, vv) in all_triples:
        a = tournament_from_bits(nn, bb)
        su = sum(a[uu])
        sv = sum(a[vv])
        diff = ht - hd - hc
        score_analysis[(su, sv)][diff] += 1

    for (su, sv) in sorted(score_analysis.keys()):
        d = score_analysis[(su, sv)]
        total_here = sum(d.values())
        zeros = d.get(0, 0)
        print(f"  scores=({su},{sv}): {zeros}/{total_here} exact, "
              f"diffs={dict(sorted(d.items()))}")

    # CRUCIAL: Check the REVERSE contraction convention
    # Maybe (x,w) iff (x,v) and (w,x) iff (u,x)?
    print(f"\n--- Testing REVERSE contraction convention ---")
    print(f"  (x,w) iff (x,v), (w,x) iff (u,x)")
    rev_holds = 0
    rev_diffs = defaultdict(int)
    for bits in range(num_tournaments):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)
        edges = get_edges(adj, n)
        for (u, v) in edges:
            # Reverse contraction: w inherits IN from v, OUT from u
            con_adj_rev, con_n_rev = contract_edge_reverse(adj, n, u, v)
            H_con_rev = ham_paths_dp(con_adj_rev, con_n_rev)
            del_adj = delete_edge(adj, n, u, v)
            H_del = ham_paths_dp(del_adj, n)
            d = H_T - H_del - H_con_rev
            rev_diffs[d] += 1
            if d == 0:
                rev_holds += 1
    print(f"  H(T) = H(T\\e) + H(T/e_rev): {rev_holds}/{total_tests}")
    print(f"  Diff distribution:")
    for d in sorted(rev_diffs.keys()):
        print(f"    diff = {d:+d}: {rev_diffs[d]} cases")

    # Also check subtraction: H(T) = H(T\e) - H(T/e_rev)
    sub_holds = sum(1 for d in rev_diffs if d + 2*rev_diffs.get(d,0) for _ in [0])  # dummy
    # Actually let's just check
    sub_rev_diffs = defaultdict(int)
    for bits in range(num_tournaments):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)
        edges = get_edges(adj, n)
        for (u, v) in edges:
            con_adj_rev, con_n_rev = contract_edge_reverse(adj, n, u, v)
            H_con_rev = ham_paths_dp(con_adj_rev, con_n_rev)
            del_adj = delete_edge(adj, n, u, v)
            H_del = ham_paths_dp(del_adj, n)
            d = H_T - H_del + H_con_rev
            sub_rev_diffs[d] += 1
    print(f"\n  H(T) = H(T\\e) - H(T/e_rev):")
    for d in sorted(sub_rev_diffs.keys()):
        print(f"    diff = {d:+d}: {sub_rev_diffs[d]} cases")

    return all_triples


def contract_edge_reverse(adj, n, u, v):
    """
    Reverse contraction: w inherits IN-edges from v, OUT-edges from u.
    (x,w) iff (x,v), (w,x) iff (u,x)
    """
    others = sorted([x for x in range(n) if x != u and x != v])
    new_n = n - 1
    w_pos = 0
    other_map = {}
    idx = 1
    for x in others:
        other_map[x] = idx
        idx += 1

    new_adj = [[False]*new_n for _ in range(new_n)]
    for x in others:
        for y in others:
            if adj[x][y]:
                new_adj[other_map[x]][other_map[y]] = True
    for x in others:
        if adj[x][v]:
            new_adj[other_map[x]][w_pos] = True
        if adj[u][x]:
            new_adj[w_pos][other_map[x]] = True

    return new_adj, new_n


def deep_analysis(all_triples_n4, all_triples_n5):
    """Look for any consistent relation across both n=4 and n=5."""
    print(f"\n{'='*70}")
    print("DEEP ANALYSIS: searching for universal relation")
    print(f"{'='*70}")

    # Try: H(T) = a * H(T\e) + b * H(T/e) for rational a, b
    # Use two equations to solve for a, b, then check all others

    for label, triples in [("n=4", all_triples_n4), ("n=5", all_triples_n5)]:
        print(f"\n--- {label}: Trying H(T) = a*H(T\\e) + b*H(T/e) ---")

        # Find two linearly independent equations
        t0 = triples[0]
        ht0, hd0, hc0 = t0[0], t0[1], t0[2]

        for i in range(1, len(triples)):
            ti = triples[i]
            ht1, hd1, hc1 = ti[0], ti[1], ti[2]
            det = hd0 * hc1 - hd1 * hc0
            if det != 0:
                a = (ht0 * hc1 - ht1 * hc0) / det
                b = (hd0 * ht1 - hd1 * ht0) / det
                # Check all
                fails = 0
                for (ht, hd, hc, *_) in triples:
                    if abs(ht - a*hd - b*hc) > 1e-9:
                        fails += 1
                if fails == 0:
                    print(f"  UNIVERSAL: H(T) = {a}*H(T\\e) + {b}*H(T/e)")
                else:
                    print(f"  Candidate a={a:.6f}, b={b:.6f} from eqs 0,{i}: {fails} failures")
                break

    # Try polynomial: H(T) = H(T\e) + c*H(T/e) for various c
    print(f"\n--- Sweep: H(T) = H(T\\e) + c*H(T/e) ---")
    for c_num in range(-5, 6):
        for c_den in [1, 2, 3]:
            c = c_num / c_den
            ok4 = sum(1 for (ht, hd, hc, *_) in all_triples_n4
                      if abs(ht - hd - c*hc) < 1e-9)
            ok5 = sum(1 for (ht, hd, hc, *_) in all_triples_n5
                      if abs(ht - hd - c*hc) < 1e-9)
            if ok4 > len(all_triples_n4) * 0.3 or ok5 > len(all_triples_n5) * 0.3:
                print(f"  c={c:+.4f}: n=4 {ok4}/{len(all_triples_n4)}, "
                      f"n=5 {ok5}/{len(all_triples_n5)}")

    # Summary of the additive identity
    print(f"\n--- Summary: H(T) = H(T\\e) + H(T/e) ---")
    ok4 = sum(1 for (ht, hd, hc, *_) in all_triples_n4 if ht == hd + hc)
    ok5 = sum(1 for (ht, hd, hc, *_) in all_triples_n5 if ht == hd + hc)
    print(f"  n=4: {ok4}/{len(all_triples_n4)}")
    print(f"  n=5: {ok5}/{len(all_triples_n5)}")


def print_sample_details(n, num_samples=3):
    """Print detailed examples for understanding."""
    print(f"\n{'='*70}")
    print(f"DETAILED EXAMPLES at n={n}")
    print(f"{'='*70}")

    num_edges_total = n * (n - 1) // 2

    for bits in range(min(num_samples, 1 << num_edges_total)):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)

        print(f"\nTournament bits={bits:0{num_edges_total}b}, H(T)={H_T}")
        print(f"  Adjacency: ", end="")
        for i in range(n):
            row = ""
            for j in range(n):
                row += "1" if adj[i][j] else "0"
            print(row, end=" ")
        print()

        edges = get_edges(adj, n)
        for (u, v) in edges:
            del_adj = delete_edge(adj, n, u, v)
            H_del = ham_paths_dp(del_adj, n)
            con_adj, con_n = contract_edge(adj, n, u, v)
            H_con = ham_paths_dp(con_adj, con_n)
            con_adj_r, con_n_r = contract_edge_reverse(adj, n, u, v)
            H_con_r = ham_paths_dp(con_adj_r, con_n_r)

            diff1 = H_T - H_del - H_con
            diff2 = H_T - H_del - H_con_r
            print(f"  e=({u}->{v}): H(T\\e)={H_del}, H(T/e)={H_con}, H(T/e_rev)={H_con_r}  |  "
                  f"H-Hdel-Hcon={diff1:+d}, H-Hdel-Hcon_rev={diff2:+d}")


if __name__ == "__main__":
    # Detailed examples first
    print_sample_details(4, num_samples=3)

    # Exhaustive tests
    triples_4 = test_deletion_contraction(4, verbose=True)
    triples_5 = test_deletion_contraction(5, verbose=False)

    # Deep analysis
    deep_analysis(triples_4, triples_5)

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")
