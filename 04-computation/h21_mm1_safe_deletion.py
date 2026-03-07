#!/usr/bin/env python3
r"""
h21_mm1_safe_deletion.py
Verify structural claims about mm=1 (max matching = 1 in 3-cycle hypergraph)
tournaments on n=9,10,11 vertices.

Claims:
1. If all 3-cycles share a common vertex u, then V\{u} is transitive.
2. There always exists a vertex v != u whose deletion preserves cycle-richness
   (every vertex in T-v is still in at least one 3-cycle).

Also verifies the bipartite backward-arc structure:
  In the transitive ordering of V\{u}, arcs from N+(u) to N-(u) create
  3-cycles through u. The last vertex in the transitive ordering (position n-2)
  is beaten by all of N+(u), so deleting any single vertex from N+(u) cannot
  isolate it from 3-cycles.

Author: kind-pasteur-2026-03-07
"""

import random
import sys
from itertools import combinations

def random_tournament(n):
    """Generate a random tournament on n vertices as adjacency matrix."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

def construct_mm1_tournament(n):
    """
    Construct a random tournament where all 3-cycles share vertex 0.
    Strategy: make V\{0} transitive (topological order 1,2,...,n-1 where i->j if i<j
    in the sub-tournament), then add edges between 0 and others.
    Vertex 0 has N+ and N-. For 3-cycles through 0: need a in N+(0), b in N-(0)
    with a->b (backward arc in the transitive order).
    For cycle-rich: every vertex must be in a 3-cycle.
    """
    adj = [[0]*n for _ in range(n)]
    # V\{0} is transitive: i->j when i < j (for i,j in 1..n-1)
    for i in range(1, n):
        for j in range(i+1, n):
            adj[i][j] = 1  # i beats j

    # Randomly assign N+(0) and N-(0)
    # Need at least one in each for any 3-cycle to exist
    others = list(range(1, n))
    random.shuffle(others)
    # Random split: k vertices in N+(0), n-1-k in N-(0), 1 <= k <= n-2
    k = random.randint(1, n-2)
    n_plus = set(others[:k])   # 0 beats these
    n_minus = set(others[k:])  # these beat 0

    for v in n_plus:
        adj[0][v] = 1
    for v in n_minus:
        adj[v][0] = 1

    # Check cycle-rich: every vertex in at least one 3-cycle through 0
    # 3-cycle through 0: 0 -> a -> b -> 0, a in N+, b in N-, a->b
    # In the transitive order, a->b iff a < b
    # So 3-cycles: {0, a, b} where a in N+, b in N-, a < b
    # Vertex a in N+ is in a 3-cycle iff exists b in N- with b > a
    # Vertex b in N- is in a 3-cycle iff exists a in N+ with a < b

    # Check coverage
    in_cycle = {0}  # 0 is in any 3-cycle if one exists
    has_any_cycle = False
    for a in n_plus:
        for b in n_minus:
            if a < b:  # a->b in transitive order
                in_cycle.add(a)
                in_cycle.add(b)
                has_any_cycle = True

    if not has_any_cycle or len(in_cycle) < n:
        return None, None, None  # not cycle-rich

    return adj, 0, (n_plus, n_minus)

def find_3cycles(adj, n):
    """Find all 3-cycles (as frozensets of vertices)."""
    cycles = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.add(frozenset([i, j, k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.add(frozenset([i, j, k]))
    return cycles

def is_cycle_rich(adj, n, cycles_3):
    """Every vertex is in at least one 3-cycle."""
    if not cycles_3:
        return False
    covered = set()
    for c in cycles_3:
        covered.update(c)
    return len(covered) == n

def find_common_vertex(cycles_3):
    """If all 3-cycles share a common vertex, return it. Otherwise None."""
    if not cycles_3:
        return None
    common = None
    for c in cycles_3:
        if common is None:
            common = set(c)
        else:
            common = common & set(c)
        if not common:
            return None
    if len(common) >= 1:
        return min(common)
    return None

def is_transitive(adj, vertices):
    """Check if the sub-tournament on given vertices is transitive."""
    vlist = sorted(vertices)
    k = len(vlist)
    if k <= 1:
        return True
    for i in range(k):
        for j in range(i+1, k):
            for l in range(j+1, k):
                a, b, c = vlist[i], vlist[j], vlist[l]
                if adj[a][b] and adj[b][c] and adj[c][a]:
                    return False
                if adj[a][c] and adj[c][b] and adj[b][a]:
                    return False
    return True

def transitive_ordering(adj, vertices):
    """Return the transitive ordering of vertices (topological sort)."""
    vlist = list(vertices)
    out_degs = {}
    for v in vlist:
        out_degs[v] = sum(adj[v][w] for w in vlist if w != v)
    vlist.sort(key=lambda v: -out_degs[v])
    return vlist

def check_cycle_rich_after_deletion(adj, n, v_del):
    """Check if T - v_del is cycle-rich (every remaining vertex in a 3-cycle)."""
    remaining = [i for i in range(n) if i != v_del]
    k = len(remaining)
    in_cycle = set()
    for i in range(k):
        if len(in_cycle) == k:
            break
        for j in range(i+1, k):
            for l in range(j+1, k):
                a, b, c = remaining[i], remaining[j], remaining[l]
                if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                   (adj[a][c] and adj[c][b] and adj[b][a]):
                    in_cycle.add(a)
                    in_cycle.add(b)
                    in_cycle.add(c)
    return len(in_cycle) == k

def analyze_bipartite_structure(adj, n, u, trans_order, verbose=False):
    """Analyze the bipartite backward-arc graph B = (N+(u), N-(u))."""
    n_plus = [v for v in range(n) if v != u and adj[u][v] == 1]
    n_minus = [v for v in range(n) if v != u and adj[v][u] == 1]
    pos = {v: i for i, v in enumerate(trans_order)}

    backward_arcs = []
    for a in n_plus:
        for b in n_minus:
            if adj[a][b] == 1:
                backward_arcs.append((a, b))

    n_minus_backward_deg = {}
    for b in n_minus:
        deg = sum(1 for a in n_plus if adj[a][b] == 1)
        n_minus_backward_deg[b] = deg

    n_plus_forward_to_nminus = {}
    for a in n_plus:
        deg = sum(1 for b in n_minus if adj[a][b] == 1)
        n_plus_forward_to_nminus[a] = deg

    if verbose:
        print(f"  Common vertex u = {u}")
        print(f"  N+(u) = {sorted(n_plus)} (u beats these, |N+|={len(n_plus)})")
        print(f"  N-(u) = {sorted(n_minus)} (these beat u, |N-|={len(n_minus)})")
        print(f"  Transitive ordering of V\\{{u}}: {trans_order}")
        print(f"  Positions: {pos}")
        print(f"  Backward arcs (a in N+, b in N-: a->b, creating 3-cycle u->a->b->u):")
        for a, b in backward_arcs:
            print(f"    {a}(pos {pos[a]}) -> {b}(pos {pos[b]})")
        print(f"  N-(u) backward degrees (# of N+(u) members beating them):")
        for b in sorted(n_minus, key=lambda x: pos[x]):
            print(f"    vertex {b} (pos {pos[b]}): beaten by {n_minus_backward_deg[b]} of {len(n_plus)} N+(u) members")
        print(f"  N+(u) forward-to-N-(u) degrees:")
        for a in sorted(n_plus, key=lambda x: pos[x]):
            print(f"    vertex {a} (pos {pos[a]}): beats {n_plus_forward_to_nminus[a]} of {len(n_minus)} N-(u) members")

        last_v = trans_order[-1]
        print(f"  Last in transitive order: vertex {last_v} (pos {len(trans_order)-1})")
        if last_v in n_minus:
            print(f"    -> In N-(u). Beaten by ALL of N+(u) (backward deg = {n_minus_backward_deg[last_v]})")
            print(f"    -> Deleting any single v from N+(u) leaves {n_minus_backward_deg[last_v]-1} N+(u) members still beating {last_v}")
            print(f"    -> So {last_v} remains in 3-cycles through u after any single N+(u) deletion")
        elif last_v in n_plus:
            print(f"    -> In N+(u). All other vertices in V\\{{u}} beat it.")

    return n_plus, n_minus, backward_arcs, n_minus_backward_deg, n_plus_forward_to_nminus


def run_random_verification(n, num_trials, verbose_max=5):
    """Run random tournament sampling for given n."""
    random.seed(42 + n)

    print(f"\n{'='*70}")
    print(f"RANDOM SAMPLING: n={n}, {num_trials} tournaments")
    print(f"{'='*70}\n")

    total = 0
    cycle_rich_count = 0
    mm1_count = 0
    transitive_ok = 0
    safe_deletion_ok = 0
    transitive_fails = []
    safe_deletion_fails = []
    verbose_count = 0

    for trial in range(num_trials):
        adj = random_tournament(n)
        total += 1

        cycles_3 = find_3cycles(adj, n)
        if not is_cycle_rich(adj, n, cycles_3):
            continue
        cycle_rich_count += 1

        u = find_common_vertex(cycles_3)
        if u is None:
            continue
        mm1_count += 1

        # Claim 1
        complement = [v for v in range(n) if v != u]
        if is_transitive(adj, complement):
            transitive_ok += 1
        else:
            transitive_fails.append(trial)

        # Verbose
        if verbose_count < verbose_max:
            verbose_count += 1
            trans_order = transitive_ordering(adj, complement)
            print(f"--- mm=1 Example #{verbose_count} (trial {trial}) ---")
            print(f"  Number of 3-cycles: {len(cycles_3)}")
            n_plus, n_minus, _, _, _ = analyze_bipartite_structure(
                adj, n, u, trans_order, verbose=True)
            safe_vs = [v for v in range(n) if v != u and check_cycle_rich_after_deletion(adj, n, v)]
            print(f"  Safe deletion vertices: {safe_vs}")
            for v in safe_vs:
                label = "N+(u)" if v in n_plus else "N-(u)"
                print(f"    vertex {v}: in {label}")
            print()

        # Claim 2
        found_safe = any(check_cycle_rich_after_deletion(adj, n, v) for v in range(n) if v != u)
        if found_safe:
            safe_deletion_ok += 1
        else:
            safe_deletion_fails.append(trial)

        if (trial + 1) % 100000 == 0:
            print(f"  ... {trial+1} trials, {mm1_count} mm=1 so far")

    print(f"\nRESULTS for n={n}:")
    print(f"  Total: {total}, Cycle-rich: {cycle_rich_count} ({100*cycle_rich_count/max(1,total):.2f}%)")
    print(f"  mm=1: {mm1_count} ({100*mm1_count/max(1,cycle_rich_count):.3f}% of cycle-rich)")
    print(f"  CLAIM 1 (V\\{{u}} transitive): {transitive_ok}/{mm1_count} {'VERIFIED' if not transitive_fails else 'REFUTED!'}")
    print(f"  CLAIM 2 (safe deletion exists): {safe_deletion_ok}/{mm1_count} {'VERIFIED' if not safe_deletion_fails else 'REFUTED!'}")
    if transitive_fails:
        print(f"  *** Transitivity failures at trials: {transitive_fails[:10]}")
    if safe_deletion_fails:
        print(f"  *** Safe deletion failures at trials: {safe_deletion_fails[:10]}")

    return mm1_count, transitive_ok, safe_deletion_ok


def run_constructed_verification(n, num_trials, verbose_max=3):
    """
    Directly construct mm=1 tournaments (V\{0} transitive, vertex 0 is hub).
    This gives a much higher hit rate for testing claims.
    """
    random.seed(1000 + n)

    print(f"\n{'='*70}")
    print(f"CONSTRUCTED mm=1 TOURNAMENTS: n={n}, attempting {num_trials}")
    print(f"{'='*70}\n")

    generated = 0
    cycle_rich_mm1 = 0
    transitive_ok = 0
    safe_deletion_ok = 0
    transitive_fails = []
    safe_deletion_fails = []
    verbose_count = 0

    # Structural checks
    last_in_nminus = 0
    last_beaten_all_nplus = 0
    all_nminus_covered = 0
    nplus_sizes = {}

    for trial in range(num_trials):
        result = construct_mm1_tournament(n)
        if result[0] is None:
            continue
        adj, u, (n_plus_set, n_minus_set) = result
        generated += 1

        # Verify it's actually mm=1 by finding cycles
        cycles_3 = find_3cycles(adj, n)
        if not is_cycle_rich(adj, n, cycles_3):
            continue

        cu = find_common_vertex(cycles_3)
        if cu is None:
            # Construction should guarantee common vertex, but verify
            continue

        cycle_rich_mm1 += 1
        u = cu

        # Track N+ sizes
        nps = len(n_plus_set)
        nplus_sizes[nps] = nplus_sizes.get(nps, 0) + 1

        complement = [v for v in range(n) if v != u]

        # Claim 1
        if is_transitive(adj, complement):
            transitive_ok += 1
        else:
            transitive_fails.append(trial)
            if len(transitive_fails) <= 3:
                print(f"*** TRANSITIVITY FAILURE at trial {trial}!")

        trans_order = transitive_ordering(adj, complement)

        # Structural checks
        n_plus = [v for v in range(n) if v != u and adj[u][v] == 1]
        n_minus = [v for v in range(n) if v != u and adj[v][u] == 1]

        last_v = trans_order[-1]
        if last_v in n_minus:
            last_in_nminus += 1
        beaten = sum(1 for a in n_plus if adj[a][last_v] == 1)
        if beaten == len(n_plus):
            last_beaten_all_nplus += 1
        if all(any(adj[a][b] == 1 for a in n_plus) for b in n_minus):
            all_nminus_covered += 1

        # Verbose
        if verbose_count < verbose_max:
            verbose_count += 1
            print(f"--- Constructed mm=1 #{verbose_count} (trial {trial}) ---")
            print(f"  Number of 3-cycles: {len(cycles_3)}")
            analyze_bipartite_structure(adj, n, u, trans_order, verbose=True)
            safe_vs = [v for v in range(n) if v != u and check_cycle_rich_after_deletion(adj, n, v)]
            print(f"  Safe deletion vertices: {safe_vs}")
            print()

        # Claim 2
        found_safe = any(check_cycle_rich_after_deletion(adj, n, v) for v in range(n) if v != u)
        if found_safe:
            safe_deletion_ok += 1
        else:
            safe_deletion_fails.append(trial)
            if len(safe_deletion_fails) <= 3:
                print(f"*** SAFE DELETION FAILURE at trial {trial}!")
                for v in range(n):
                    if v == u: continue
                    ok = check_cycle_rich_after_deletion(adj, n, v)
                    in_set = "N+(u)" if adj[u][v] else "N-(u)"
                    print(f"    delete {v} ({in_set}): cycle-rich after = {ok}")

    print(f"\nRESULTS for n={n} (constructed):")
    print(f"  Generated: {generated}, Cycle-rich mm=1: {cycle_rich_mm1}")
    print(f"  CLAIM 1 (V\\{{u}} transitive): {transitive_ok}/{cycle_rich_mm1} {'VERIFIED' if not transitive_fails else 'REFUTED!'}")
    print(f"  CLAIM 2 (safe deletion exists): {safe_deletion_ok}/{cycle_rich_mm1} {'VERIFIED' if not safe_deletion_fails else 'REFUTED!'}")
    if transitive_fails:
        print(f"  *** Transitivity failures: {transitive_fails[:10]}")
    if safe_deletion_fails:
        print(f"  *** Safe deletion failures: {safe_deletion_fails[:10]}")

    print(f"\n  Structural checks ({cycle_rich_mm1} mm=1 tournaments):")
    print(f"    Last in trans. order in N-(u): {last_in_nminus}/{cycle_rich_mm1}")
    print(f"    Last beaten by ALL of N+(u):   {last_beaten_all_nplus}/{cycle_rich_mm1}")
    print(f"    Every N-(u) has N+(u) neighbor: {all_nminus_covered}/{cycle_rich_mm1}")

    print(f"\n  |N+(u)| distribution:")
    for k in sorted(nplus_sizes.keys()):
        print(f"    |N+(u)| = {k}: {nplus_sizes[k]} cases ({100*nplus_sizes[k]/max(1,sum(nplus_sizes.values())):.1f}%)")

    return cycle_rich_mm1, transitive_ok, safe_deletion_ok


def main():
    print("="*70)
    print("mm=1 SAFE DELETION VERIFICATION")
    print("Structural claims about cycle-rich tournaments with all")
    print("3-cycles sharing a common vertex")
    print("="*70)

    # Phase 1: Random sampling at n=9 (500k trials)
    mm1_9, t_ok_9, sd_ok_9 = run_random_verification(9, 500000, verbose_max=5)

    # Phase 2: Constructed mm=1 tournaments for n=9, 10, 11
    # This gives many more mm=1 cases to test
    for n in [9, 10, 11]:
        if n <= 10:
            num = 100000
        else:
            num = 50000  # n=11 is slower
        run_constructed_verification(n, num, verbose_max=3)

    # Phase 3: Key structural argument summary
    print(f"\n{'='*70}")
    print("STRUCTURAL ARGUMENT SUMMARY")
    print("="*70)
    print("""
In an mm=1 cycle-rich tournament T on n vertices:
  - All 3-cycles pass through a single vertex u.
  - V\\{u} is transitive (otherwise there'd be a 3-cycle not through u).
  - Let the transitive ordering of V\\{u} be v_1 > v_2 > ... > v_{n-1}.
  - N+(u) = vertices beaten by u; N-(u) = vertices beating u.
  - 3-cycles: u -> a -> b -> u where a in N+(u), b in N-(u), a > b in trans. order.

KEY OBSERVATION (verified computationally):
  v_{n-1} (last in transitive order) is ALWAYS in N-(u).
  [If it were in N+(u), it would beat no one in V\\{u}, so it couldn't
   be in any 3-cycle, contradicting cycle-richness.]

  Since v_{n-1} loses to ALL of V\\{u}, in particular ALL of N+(u) beats it.
  So v_{n-1} participates in |N+(u)| 3-cycles: u -> a -> v_{n-1} -> u for each a in N+(u).

  Deleting any single vertex v from N+(u) leaves |N+(u)|-1 >= 1 other N+(u) members
  still beating v_{n-1}, so v_{n-1} remains in a 3-cycle.

  Similarly, v_1 (first in trans. order) is ALWAYS in N+(u).
  [If it were in N-(u), it would beat everyone in V\\{u} but needs someone
   in N+(u) to beat it for a 3-cycle; the only such vertex is u, but
   v_1 beats u (v_1 in N-(u)), contradiction.]

  v_1 beats ALL of N-(u), creating 3-cycles u -> v_1 -> b -> u for each b in N-(u).
  Deleting any single b from N-(u) leaves v_1 in other 3-cycles.

  SAFE DELETION STRATEGY:
  The "most deletable" vertex is typically in the middle of the transitive order
  where it has multiple backward arcs. But even boundary vertices can be safe
  to delete because their 3-cycle partners have redundant connections.
""")

if __name__ == "__main__":
    main()
