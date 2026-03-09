"""
Two-mechanism proof attempt for good-vertex existence.

Mechanism 1: v uniquely dominates some cycle C (dom(C) = {v}).
  Removing v frees C, potentially creating a free component.

Mechanism 2: v is in ALL dominated cycles of some component K.
  Removing v removes all dominated members from K, making it all-free.

For v to be bad, at least one mechanism must apply.
For v to be GOOD, NEITHER mechanism applies.

Claim: there always exists v with neither mechanism.

Test: for each tournament, find vertices safe from BOTH mechanisms.
"""
import numpy as np
from itertools import combinations
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

def main():
    rng = np.random.RandomState(42)

    print("=" * 70)
    print("TWO-MECHANISM PROOF: good vertex = safe from both mechanisms")
    print("=" * 70)

    for n in [5, 6, 7, 8, 9, 10]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000, 10: 1000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        safe_exists = 0
        safe_is_good = 0  # among safe vertices, how many are actually good?
        safe_is_bad = 0   # safe vertex is bad (would disprove the approach)
        no_safe = 0       # no vertex is safe from both mechanisms
        n_with_bad = 0
        total_tested = 0

        # More refined: mechanism 2 needs component awareness
        m2_bad_sizes = defaultdict(int)

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            b1 = compute_b1(A, n)
            if b1 != 0:
                continue
            total_tested += 1

            cycles = find_3cycles(A, n)
            if not cycles:
                continue
            comps = cycle_graph_components(cycles)

            # Mechanism 1: v uniquely dominates some cycle not through v
            m1_bad = set()  # vertices with unique_dom > 0
            for v in range(n):
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    if get_dom(A, n, cyc) == [v]:
                        m1_bad.add(v)
                        break

            # Mechanism 2: v is in ALL dominated cycles of some component
            m2_bad = set()
            for comp in comps:
                # Find dominated cycles in this component
                dom_in_comp = []
                for ci in comp:
                    if is_dominated(A, n, cycles[ci]):
                        dom_in_comp.append(ci)

                if not dom_in_comp:
                    continue  # all-free component (shouldn't exist if b1=0... unless it's there)

                # Intersection of vertex sets of dominated cycles
                inter = set(cycles[dom_in_comp[0]])
                for ci in dom_in_comp[1:]:
                    inter &= set(cycles[ci])

                m2_bad.update(inter)

            # Safe vertices: not in m1_bad and not in m2_bad
            safe = set(range(n)) - m1_bad - m2_bad

            # Compute actual bad vertices
            actual_bad = set()
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    actual_bad.add(v)

            if actual_bad:
                n_with_bad += 1

            if safe:
                safe_exists += 1
                # Check if safe vertices are actually good
                for v in safe:
                    if v in actual_bad:
                        safe_is_bad += 1
                    else:
                        safe_is_good += 1
            else:
                no_safe += 1

            m2_bad_sizes[len(m2_bad)] += 1

        print(f"\nn={n} ({method}): {total_tested} tournaments with b1=0")
        print(f"  Safe vertex exists: {safe_exists}/{total_tested} ({100*safe_exists/max(total_tested,1):.1f}%)")
        print(f"  No safe vertex: {no_safe}/{total_tested}")
        print(f"  Safe and actually GOOD: {safe_is_good}")
        print(f"  Safe but actually BAD: {safe_is_bad} {'*** DISPROVES APPROACH ***' if safe_is_bad > 0 else '(0 = approach valid so far)'}")
        print(f"  Mechanism 2 bad set sizes: {dict(m2_bad_sizes)}")

    # PART 2: More detailed analysis of safe-but-bad cases
    print("\n" + "=" * 70)
    print("PART 2: Analyzing cases where 'safe' vertex is bad")
    print("-" * 50)

    rng = np.random.RandomState(42)
    for n in [7, 8]:
        nsamp = {7: 10000, 8: 5000}[n]
        bad_safe_count = 0

        for _ in range(nsamp):
            A = random_tournament(n, rng)
            if compute_b1(A, n) != 0:
                continue

            cycles = find_3cycles(A, n)
            if not cycles:
                continue
            comps = cycle_graph_components(cycles)

            m1_bad = set()
            for v in range(n):
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    if get_dom(A, n, cyc) == [v]:
                        m1_bad.add(v)
                        break

            m2_bad = set()
            for comp in comps:
                dom_in_comp = [ci for ci in comp if is_dominated(A, n, cycles[ci])]
                if not dom_in_comp:
                    continue
                inter = set(cycles[dom_in_comp[0]])
                for ci in dom_in_comp[1:]:
                    inter &= set(cycles[ci])
                m2_bad.update(inter)

            safe = set(range(n)) - m1_bad - m2_bad

            for v in safe:
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_safe_count += 1
                    print(f"\n  n={n}: SAFE vertex {v} is BAD!")
                    print(f"    m1_bad={m1_bad}, m2_bad={m2_bad}, safe={safe}")

                    # Diagnose: how does removing v create b1=1?
                    # Check which cycles go through v
                    thru_v = [cyc for cyc in cycles if v in set(cyc)]
                    not_thru_v = [cyc for cyc in cycles if v not in set(cyc)]
                    print(f"    #cycles through v: {len(thru_v)}")

                    # Check free cycles not through v
                    free_not_v = [cyc for cyc in not_thru_v if not is_dominated(A, n, cyc)]
                    dom_not_v = [cyc for cyc in not_thru_v if is_dominated(A, n, cyc)]
                    print(f"    free not through v: {len(free_not_v)}")
                    print(f"    dominated not through v: {len(dom_not_v)}")

                    # Check which dominated cycles not through v lose ALL dominators
                    for cyc in dom_not_v:
                        doms = get_dom(A, n, cyc)
                        remaining_doms = [d for d in doms if d != v]
                        if not remaining_doms:
                            print(f"    CYCLE {cyc} loses all dominators!")
                            print(f"      Original doms: {doms}")

                    # Check component structure
                    r_comps = cycle_graph_components(not_thru_v)
                    for ci, comp in enumerate(r_comps):
                        all_free = all(not is_dominated(A, n, not_thru_v[j]) for j in comp)
                        verts = set()
                        for j in comp:
                            verts.update(not_thru_v[j])
                        if all_free:
                            print(f"    RESTRICTED FREE COMP: {len(comp)} cycles, verts={sorted(verts)}")

                    if bad_safe_count >= 5:
                        break
            if bad_safe_count >= 5:
                break

        if bad_safe_count == 0:
            print(f"\n  n={n}: 0 safe-but-bad vertices found!")

    # PART 3: Refined mechanism 2 - track per-component vertex intersection
    # Maybe the mechanism is more nuanced: not just "in ALL dom cycles"
    # but "in all dom cycles that are the only connection between free and dom parts"
    print("\n" + "=" * 70)
    print("PART 3: Component-internal bridge analysis")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        if not cycles:
            continue
        comps = cycle_graph_components(cycles)

        actual_bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                actual_bad.add(v)

        if len(actual_bad) != 3:
            continue

        # For each component: analyze internal bridge structure
        for comp in comps:
            dom_cycles = [ci for ci in comp if is_dominated(A, n, cycles[ci])]
            free_cycles = [ci for ci in comp if not is_dominated(A, n, cycles[ci])]

            if not dom_cycles or not free_cycles:
                continue

            # Find vertices that appear in ALL dominated cycles of this component
            inter = set(cycles[dom_cycles[0]])
            for ci in dom_cycles[1:]:
                inter &= set(cycles[ci])

            # Also: vertices that appear in dominated cycles that BRIDGE to free cycles
            nc = len(cycles)
            adj_list = [set() for _ in range(nc)]
            for i in range(nc):
                for j in range(i+1, nc):
                    if shared_directed_edge(cycles[i], cycles[j]):
                        adj_list[i].add(j)
                        adj_list[j].add(i)

            # Bridge dominated cycles: those adjacent to at least one free cycle
            bridge_dom = []
            for di in dom_cycles:
                if any(fi in adj_list[di] for fi in free_cycles):
                    bridge_dom.append(di)

            # Vertices in ALL bridge dominated cycles
            if bridge_dom:
                bridge_inter = set(cycles[bridge_dom[0]])
                for ci in bridge_dom[1:]:
                    bridge_inter &= set(cycles[ci])
            else:
                bridge_inter = set()

            if bits < 500 and len(actual_bad) == 3:
                print(f"  bits={bits}: #dom={len(dom_cycles)}, #free={len(free_cycles)}, "
                      f"#bridge_dom={len(bridge_dom)}")
                print(f"    inter_all_dom={sorted(inter)}, "
                      f"inter_bridge_dom={sorted(bridge_inter)}, "
                      f"bad={sorted(actual_bad)}")
                break  # Just first component

if __name__ == '__main__':
    main()
