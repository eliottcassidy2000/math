"""
Key insight from beta2_four_bad_contradiction.py:
- Middle vertex of bad triple has 0 uniquely-dominated cycles
- Deleting middle vertex SPLITS the 3-cycle graph from 1 component to 2 components
- One of those 2 components is all-free (giving b1=1)

This means: the middle vertex's cycles are the ONLY connection between a free
subgraph and a dominated subgraph in the cycle graph.

For 4 bad vertices v1->v2->v3->v4 (transitive):
- v1 (top): 1 uniquely-dominated from above, deletion keeps 1 component
- v4 (bot): 1 uniquely-dominated from below, deletion keeps 1 component
- v2, v3 (bridges): 0 uniquely-dominated, deletion splits into 2 components

Question: Can TWO vertices both be cycle-graph bridges simultaneously?

Test: For each bad vertex, compute:
1. Whether the cycle graph changes component count upon deletion
2. The MECHANISM: which cycles through v connect which parts?
3. Can two "bridge" vertices coexist?
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def get_dom(A, n, cyc):
    doms = []
    a, b, c = cyc
    for d in range(n):
        if d in {a,b,c}: continue
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def cycle_graph_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]: continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]: continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]: stack.append(u)
        components.append(comp)
    return components

def is_dominated(A, n, cyc):
    return len(get_dom(A, n, cyc)) > 0

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    return sum(1 for comp in comps if all(not is_dominated(A, n, cycles[ci]) for ci in comp))

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def is_sc(A, n):
    for start in range(n):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited: continue
            visited.add(v)
            for u in range(n):
                if A[v][u] and u not in visited:
                    stack.append(u)
        if len(visited) < n:
            return False
    return True

def main():
    print("=" * 70)
    print("BRIDGE VERTEX ANALYSIS: Component splitting mechanism")
    print("=" * 70)

    # PART 1: Classify bad vertices by component-splitting behavior
    print("\nPART 1: Component change upon bad vertex deletion")
    print("-" * 50)

    rng = np.random.RandomState(42)
    for n in [6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        # For each bad vertex: track (unique_dom_count, delta_components)
        bad_type = defaultdict(int)
        # Track the component delta for each position in the transitive chain
        pos_delta = defaultdict(lambda: defaultdict(int))

        n_with_bad = 0
        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue

            bad = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.append(v)

            if len(bad) == 0:
                continue
            n_with_bad += 1

            cycles = find_3cycles(A, n)
            comps = cycle_graph_components(cycles)
            n_comps = len(comps)

            # For each bad vertex: compute unique dom count and component delta
            for v in bad:
                # unique dom count
                n_unique = sum(1 for cyc in cycles
                              if v not in set(cyc) and get_dom(A, n, cyc) == [v])

                # Component count of cycle graph restricted to cycles not through v
                restricted = [cyc for cyc in cycles if v not in set(cyc)]
                if restricted:
                    r_comps = cycle_graph_components(restricted)
                    n_r_comps = len(r_comps)
                else:
                    n_r_comps = 0

                delta = n_r_comps - n_comps
                bad_type[(n_unique, delta)] += 1

            # Classify positions in transitive ordering
            if len(bad) >= 2:
                # Sort bad by score or try all orderings
                for perm in permutations(bad):
                    is_transitive = True
                    for i in range(len(perm)):
                        for j in range(i+1, len(perm)):
                            if not A[perm[i]][perm[j]]:
                                is_transitive = False
                                break
                        if not is_transitive:
                            break
                    if is_transitive:
                        ordered = perm
                        break

                for pos_idx, v in enumerate(ordered):
                    n_unique = sum(1 for cyc in cycles
                                  if v not in set(cyc) and get_dom(A, n, cyc) == [v])
                    restricted = [cyc for cyc in cycles if v not in set(cyc)]
                    if restricted:
                        r_comps = cycle_graph_components(restricted)
                        n_r_comps = len(r_comps)
                    else:
                        n_r_comps = 0
                    delta = n_r_comps - n_comps

                    # Position: 0=top, -1=bot, middle=anything else
                    if pos_idx == 0:
                        pos_label = "top"
                    elif pos_idx == len(ordered) - 1:
                        pos_label = "bot"
                    else:
                        pos_label = "mid"

                    pos_delta[pos_label][(n_unique, delta)] += 1

        print(f"\nn={n} ({method}): {n_with_bad} tournaments with bad vertices")
        print(f"  Bad vertex types (unique_dom, delta_components):")
        for key, cnt in sorted(bad_type.items()):
            print(f"    unique_dom={key[0]}, delta_comp={key[1]}: {cnt}")

        if pos_delta:
            print(f"  Position breakdown:")
            for pos in ['top', 'mid', 'bot']:
                if pos in pos_delta:
                    print(f"    {pos}:")
                    for key, cnt in sorted(pos_delta[pos].items()):
                        print(f"      unique_dom={key[0]}, delta_comp={key[1]}: {cnt}")

    # PART 2: For tournaments with 2 bad vertices - what's the structure?
    print("\n" + "=" * 70)
    print("PART 2: Two bad vertices - are both 'bridges' ever?")
    print("-" * 50)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            rng2 = np.random.RandomState(123)
            nsamp = 5000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        two_bad_both_bridge = 0
        two_bad_one_bridge = 0
        two_bad_no_bridge = 0
        three_bad_bridges = defaultdict(int)
        total_two_bad = 0
        total_three_bad = 0

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng2)

            if compute_b1(A, n) != 0:
                continue

            bad = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.append(v)

            if len(bad) < 2:
                continue

            cycles = find_3cycles(A, n)
            comps = cycle_graph_components(cycles)
            n_comps = len(comps)

            # Classify each bad vertex as bridge (delta > 0) or not
            bridge_count = 0
            for v in bad:
                restricted = [cyc for cyc in cycles if v not in set(cyc)]
                if restricted:
                    r_comps = cycle_graph_components(restricted)
                    delta = len(r_comps) - n_comps
                else:
                    delta = -n_comps
                if delta > 0:
                    bridge_count += 1

            if len(bad) == 2:
                total_two_bad += 1
                if bridge_count == 2:
                    two_bad_both_bridge += 1
                elif bridge_count == 1:
                    two_bad_one_bridge += 1
                else:
                    two_bad_no_bridge += 1
            elif len(bad) == 3:
                total_three_bad += 1
                three_bad_bridges[bridge_count] += 1

        print(f"\nn={n} ({method}):")
        print(f"  2-bad tournaments: {total_two_bad}")
        if total_two_bad > 0:
            print(f"    both bridge: {two_bad_both_bridge}")
            print(f"    one bridge: {two_bad_one_bridge}")
            print(f"    no bridge: {two_bad_no_bridge}")
        print(f"  3-bad tournaments: {total_three_bad}")
        for bc, cnt in sorted(three_bad_bridges.items()):
            print(f"    {bc} bridges: {cnt}")

    # PART 3: Deeper analysis of the bridge mechanism
    # When mid is deleted and cycle graph splits, what are the two parts?
    print("\n" + "=" * 70)
    print("PART 3: Anatomy of the cycle graph split")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    split_anatomy = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        if not is_sc(A, n):
            continue
        kappa1 = any(not is_sc(delete_vertex(A, n, v), n-1) for v in range(n))
        if kappa1:
            continue

        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        if len(bad) != 3:
            continue

        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)

        # Remove cycles through mid
        restricted = [cyc for cyc in cycles if mid not in set(cyc)]
        r_comps = cycle_graph_components(restricted)

        # Classify each component: all-free or has-dominated
        comp_types = []
        for comp in r_comps:
            all_free = all(not is_dominated(A, n, restricted[ci]) for ci in comp)
            verts = set()
            for ci in comp:
                verts.update(restricted[ci])
            comp_types.append(('free' if all_free else 'mixed', len(comp), sorted(verts)))

        # Sort by type
        comp_types.sort()
        key = tuple((t, nc) for t, nc, v in comp_types)
        split_anatomy[key] += 1

        # Print first few
        if bits < 500:
            print(f"  bits={bits}: mid={mid}, components after removing mid's cycles:")
            for t, nc, v in comp_types:
                print(f"    {t}: {nc} cycles, verts={v}")

    print(f"\nSplit anatomy distribution:")
    for key, cnt in sorted(split_anatomy.items(), key=lambda x: -x[1]):
        print(f"  {key}: {cnt}")

    # PART 4: For EACH tournament with 3 bad, what's the detailed freed-set picture?
    # Key: For the bridge (mid) vertex, its freed set = free cycles not through mid.
    # These must form a connected all-free component in T\mid.
    # In T, these free cycles exist but are connected to dominated cycles via mid's cycles.
    # Removing mid's cycles severs this connection.
    #
    # This means: every dominated cycle in T that's connected to a free cycle
    # must be connected THROUGH a cycle that goes through mid.
    #
    # For a 4th bad vertex to also be a bridge, it would need a DIFFERENT set of
    # connecting cycles. But if mid's cycles are the ONLY connection, removing a
    # different vertex's cycles can't also sever the connection.
    print("\n" + "=" * 70)
    print("PART 4: Are mid's cycles the ONLY bridge between free and dominated parts?")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        if not is_sc(A, n):
            continue
        kappa1 = any(not is_sc(delete_vertex(A, n, v), n-1) for v in range(n))
        if kappa1:
            continue

        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        if len(bad) != 3:
            continue

        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)

        # Identify free and dominated cycles
        free_idxs = [i for i, cyc in enumerate(cycles) if not is_dominated(A, n, cyc)]
        dom_idxs = [i for i, cyc in enumerate(cycles) if is_dominated(A, n, cyc)]

        # Build cycle adjacency
        nc = len(cycles)
        adj = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if shared_directed_edge(cycles[i], cycles[j]):
                    adj[i].add(j)
                    adj[j].add(i)

        # Find all edges between free and dominated cycles
        bridge_edges = []
        for fi in free_idxs:
            for di in dom_idxs:
                if di in adj[fi]:
                    bridge_edges.append((fi, di))

        # How many of these bridge edges involve mid?
        mid_bridges = [(fi, di) for fi, di in bridge_edges
                       if mid in set(cycles[fi]) or mid in set(cycles[di])]
        other_bridges = [(fi, di) for fi, di in bridge_edges
                        if mid not in set(cycles[fi]) and mid not in set(cycles[di])]

        if bits < 300:
            print(f"  bits={bits}: mid={mid}, #bridge_edges={len(bridge_edges)}, "
                  f"mid_bridges={len(mid_bridges)}, other_bridges={len(other_bridges)}")

        if other_bridges:
            print(f"  *** OTHER BRIDGES EXIST at bits={bits}! ***")
            break

    print("\n  (If no '*** OTHER BRIDGES' message, mid's cycles are the ONLY free-dom bridge)")

    # PART 5: At n=7 - verify the pattern holds
    print("\n" + "=" * 70)
    print("PART 5: n=7 verification of bridge structure")
    print("-" * 50)

    rng = np.random.RandomState(42)
    n = 7
    nsamp = 3000
    n_three_bad = 0
    bridge_count_dist = defaultdict(int)
    mid_only_bridge_count = 0
    other_bridge_exists = 0

    for _ in range(nsamp):
        A = random_tournament(n, rng)
        if compute_b1(A, n) != 0:
            continue

        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)

        if len(bad) != 3:
            continue
        n_three_bad += 1

        # Find transitive order
        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)
        comps = cycle_graph_components(cycles)
        n_comps = len(comps)

        # Count bridges
        bridges = 0
        for v in bad:
            restricted = [cyc for cyc in cycles if v not in set(cyc)]
            if restricted:
                r_comps = cycle_graph_components(restricted)
                if len(r_comps) > n_comps:
                    bridges += 1
        bridge_count_dist[bridges] += 1

        # Check if mid's cycles are the only free-dom bridge
        free_idxs = [i for i, cyc in enumerate(cycles) if not is_dominated(A, n, cyc)]
        dom_idxs = [i for i, cyc in enumerate(cycles) if is_dominated(A, n, cyc)]

        nc = len(cycles)
        adj_set = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if shared_directed_edge(cycles[i], cycles[j]):
                    adj_set[i].add(j)
                    adj_set[j].add(i)

        bridge_edges = [(fi, di) for fi in free_idxs for di in dom_idxs if di in adj_set[fi]]
        other_br = [(fi, di) for fi, di in bridge_edges
                    if mid not in set(cycles[fi]) and mid not in set(cycles[di])]

        if other_br:
            other_bridge_exists += 1
        else:
            mid_only_bridge_count += 1

    print(f"n=7: {n_three_bad} tournaments with 3 bad vertices")
    print(f"  Bridge count distribution: {dict(bridge_count_dist)}")
    print(f"  Mid is only free-dom bridge: {mid_only_bridge_count}/{n_three_bad}")
    print(f"  Other bridges exist: {other_bridge_exists}/{n_three_bad}")

if __name__ == '__main__':
    main()
