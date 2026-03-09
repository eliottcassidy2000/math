"""
Attempt to prove |BAD| <= 3 (and hence good vertex exists for n >= 4).

Key structural facts:
1. Bad set is always TRANSITIVE (HYP-287)
2. Domination direction: top=above only, bot=below only, mid=bridge (HYP-289)
3. f(v) >= n-3 necessary for bad (HYP-288)
4. b1 <= 1 (THM-107)

Strategy: Show that 4 bad vertices create a contradiction.

If v1 -> v2 -> v3 -> v4 (transitive chain of length 4), then:
- v1 uniquely dominates from above: exists cycle C with dom(C) = {v1}, v1 beats all of C
- v4 uniquely dominates from below: exists cycle C with dom(C) = {v4}, C beats all of v4
- v2, v3 are bridges: no unique domination

For v_i to be bad: the freed cycles in T\v_i must span all n-1 remaining vertices,
specifically including all other bad vertices.

Key question: for v2 (a bridge), how does it "free" cycles? It must be that:
- v2 is NOT a unique dominator of any cycle (above=0, below=0)
- But deleting v2 must create a free component spanning n-1 vertices

This means: all cycles freed by deleting v2 are ALREADY free (dom=0),
and these free cycles don't contain v2 but DO span n-1 vertices.
Wait - if they're already free, b1(T) > 0. Contradiction with b1(T) = 0!

Unless: cycles with dom = {v2} DO exist but v2 dominates them from MIXED direction
(some from above, some from below). Let me re-check the data.

Actually HYP-289 says mid=(0,0) means 0 uniquely-above AND 0 uniquely-below.
So the middle vertex in a 3-bad triple has NO uniquely-dominated cycles at all.
How can it be bad?

The freed set for v_mid = {free cycles not through v_mid} union {cycles with dom={v_mid}}.
If dom={v_mid} is empty (no unique domination), then freed set = free cycles not through v_mid.
For b1(T) = 0, there are no free components, meaning every connected component of the
3-cycle graph has at least one dominated cycle.

But wait: there CAN be free cycles in T even when b1=0! A "free component" means a
connected component ALL of whose cycles are free. If free cycles exist but are all
connected (via shared edges) to at least one dominated cycle, then b1=0.

When we delete v_mid: the dominated cycles that had dom={v_mid} become free.
But HYP-289 says mid has 0 such cycles! So deleting v_mid only changes cycles
CONTAINING v_mid (which are removed, not freed).

Wait, I need to be more careful. When we delete v_mid:
- Cycles not involving v_mid remain, with their domination status possibly changed
  (if some dominator was v_mid, but v_mid is removed - wait, v_mid IS removed,
   so it can't dominate anymore)

Actually: a cycle C not containing v_mid was dominated by v_mid if v_mid was a dominator.
After deleting v_mid, C loses v_mid as a dominator. If v_mid was the ONLY dominator,
C becomes free in T\v_mid.

dom(C) = {v_mid} in T means C has unique dominator v_mid.
HYP-289 says mid has 0 such cycles! (both above=0 and below=0)

So the freed set for v_mid = {free cycles not through v_mid} + {nothing new}
= {free cycles not through v_mid}.

For v_mid to be bad, b1(T\v_mid) = 1 while b1(T) = 0.
The free cycles not through v_mid must form a connected spanning component
(spanning all n-1 vertices except v_mid) in T\v_mid.

But these same free cycles exist in T! And in T, they're connected to dominated
cycles too (otherwise b1(T) >= 1). So the free component in T\v_mid must arise
from: some component in T that wasn't all-free becomes all-free after removing
cycles through v_mid.

This means: v_mid was in some cycle C that was the ONLY dominated cycle in its
connected component. Removing v_mid removes C, making the rest all-free.

Let me verify this computationally and then see if it creates a contradiction
when combined with 4 bad vertices.
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
    print("PROOF: |BAD| <= 3 via domination structure analysis")
    print("=" * 70)

    # PART 1: For middle vertex in bad triple, HOW does it become bad?
    # If mid has no uniquely-dominated cycles, what mechanism creates b1(T\mid) = 1?
    print("\nPART 1: Mechanism for middle vertex badness")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    mid_mechanism = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        if not is_sc(A, n):
            continue
        # Check kappa >= 2
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

        # Find transitive order
        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)

        # For the middle vertex: analyze its role
        # Count cycles through mid
        mid_cycles = [cyc for cyc in cycles if mid in set(cyc)]
        mid_not_cycles = [cyc for cyc in cycles if mid not in set(cyc)]

        # Which of mid_not_cycles have dom = {mid}?
        uniquely_mid_dom = []
        for cyc in mid_not_cycles:
            doms = get_dom(A, n, cyc)
            if doms == [mid]:
                uniquely_mid_dom.append(cyc)

        # Free cycles not through mid
        free_not_mid = []
        for cyc in mid_not_cycles:
            if not is_dominated(A, n, cyc):
                free_not_mid.append(cyc)

        # Cycles through mid that are the ONLY dominated cycle in their component
        # (removing mid would remove them, possibly freeing a component)
        comps = cycle_graph_components(cycles)
        cycle_to_comp = {}
        for ci, comp in enumerate(comps):
            for idx in comp:
                cycle_to_comp[idx] = ci

        # For each component: count dominated cycles that go through mid
        comp_analysis = defaultdict(lambda: {'total': 0, 'dom': 0, 'dom_thru_mid': 0, 'free': 0})
        for i, cyc in enumerate(cycles):
            ci = cycle_to_comp[i]
            comp_analysis[ci]['total'] += 1
            if is_dominated(A, n, cyc):
                comp_analysis[ci]['dom'] += 1
                if mid in set(cyc):
                    comp_analysis[ci]['dom_thru_mid'] += 1
            else:
                comp_analysis[ci]['free'] += 1

        # If a component has ALL dominated cycles going through mid,
        # then deleting mid makes it all-free
        mechanism = "unknown"
        for ci, info in comp_analysis.items():
            if info['free'] > 0 and info['dom'] > 0 and info['dom'] == info['dom_thru_mid']:
                mechanism = "dom_removal"
                break

        mid_mechanism[mechanism] += 1
        mid_mechanism[f"unique_mid_dom={len(uniquely_mid_dom)}"] += 1
        mid_mechanism[f"free_not_mid={len(free_not_mid)}"] += 1

    print("Middle vertex mechanism distribution:")
    for k, v in sorted(mid_mechanism.items()):
        print(f"  {k}: {v}")

    # PART 2: Key counting argument for 4-bad impossibility
    print("\n" + "=" * 70)
    print("PART 2: Why 4 bad vertices are impossible")
    print("-" * 50)

    # For each tournament with 3 bad: what is the FREED set structure?
    # For each bad vertex v: which cycles become newly free in T\v?
    # (cycles with dom={v} that don't contain v)

    for n in [6, 7]:
        print(f"\n--- n={n} ---")
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
        else:
            rng = np.random.RandomState(42)
            sample = range(3000)

        # Track: for each bad vertex, count of cycles that participate in
        # each "freed" mechanism
        freed_stats = defaultdict(int)
        dom_overlap = defaultdict(int)

        count = 0
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

            if len(bad) < 2:
                continue
            count += 1

            cycles = find_3cycles(A, n)

            # For each bad vertex: which cycles have dom = {v}?
            dom_sets = {}  # v -> list of cycles uniquely dominated by v
            for v in bad:
                dom_sets[v] = []
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    doms = get_dom(A, n, cyc)
                    if doms == [v]:
                        dom_sets[v].append(cyc)

            # Do any two bad vertices share a uniquely-dominated cycle?
            # (That would mean the cycle has dom = {v1, v2}, not {v1} or {v2})
            # Actually no - if dom(C) = {v} means v is the ONLY dominator.
            # Two different bad vertices can uniquely dominate DIFFERENT cycles.

            # Key: cycle C with dom(C) = {v1} is freed only by deleting v1.
            # Cycle C with dom(C) = {v2} is freed only by deleting v2.
            # No cycle can be in both freed sets (since unique dom means exactly one).

            # But: cycles with dom=empty (free cycles) appear in ALL freed sets.

            freed_stats[f"n_bad={len(bad)}"] += 1

            # Count: for each pair of bad vertices, how many free cycles
            # contain NEITHER?
            free_cycles = [cyc for cyc in cycles if len(get_dom(A, n, cyc)) == 0]

            if len(bad) >= 2:
                for i in range(len(bad)):
                    for j in range(i+1, len(bad)):
                        vi, vj = bad[i], bad[j]
                        shared_free = sum(1 for cyc in free_cycles
                                         if vi not in set(cyc) and vj not in set(cyc))
                        dom_overlap[f"shared_free_pair"] += shared_free
                        dom_overlap[f"shared_free_pair_count"] += 1

            # For bad vertex with NO unique domination (bridge):
            # its freed set = free cycles not through it
            for v in bad:
                n_unique = len(dom_sets[v])
                if n_unique == 0:
                    freed_stats["bridge_bad"] += 1
                    # How many free cycles not through v?
                    free_not_v = sum(1 for cyc in free_cycles if v not in set(cyc))
                    freed_stats[f"bridge_free_not_v={free_not_v}"] += 1

        print(f"  Analyzed {count} tournaments with >= 2 bad")
        for k, v in sorted(freed_stats.items()):
            print(f"  {k}: {v}")
        if dom_overlap:
            avg = dom_overlap["shared_free_pair"] / max(dom_overlap["shared_free_pair_count"], 1)
            print(f"  avg shared free cycles per bad pair: {avg:.2f}")

    # PART 3: The DEFINITIVE counting argument
    # For each vertex v, define:
    #   R(v) = set of vertices reachable by freed cycles from v
    # If v is bad: R(v) must be all of V\{v}.
    # The freed cycles for v are: {C : v not in C, dom(C) subset {v}}
    # = {free cycles not through v} + {cycles with dom={v}}
    #
    # For 4 bad vertices v1,v2,v3,v4:
    # Each freed set must span n-1 vertices.
    # But the cycles with dom={v_i} are disjoint across different i.
    # The free cycles are shared.
    #
    # Total "freed cycle-vertex incidences":
    # = (n-3) * 3 * |free cycles| + 3 * sum_i |dom={v_i} cycles|
    # Must be >= 4 * (n-1) (since each of 4 freed sets must span n-1 vertices)
    #
    # But actually need CONNECTIVITY too, not just coverage.

    print("\n" + "=" * 70)
    print("PART 3: Counting freed-cycle vertex incidences")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)

        if len(bad) != 3:
            continue

        cycles = find_3cycles(A, n)

        # Classify all cycles
        free_cycs = []
        unique_dom_cycs = defaultdict(list)  # v -> cycles with dom={v}
        multi_dom_cycs = []

        for cyc in cycles:
            doms = get_dom(A, n, cyc)
            if len(doms) == 0:
                free_cycs.append(cyc)
            elif len(doms) == 1:
                unique_dom_cycs[doms[0]].append(cyc)
            else:
                multi_dom_cycs.append(cyc)

        # For each bad vertex: verify freed set spans n-1 vertices
        for v in bad:
            freed = [cyc for cyc in free_cycs if v not in set(cyc)]
            freed += unique_dom_cycs.get(v, [])
            spanned = set()
            for cyc in freed:
                spanned.update(cyc)
            assert len(spanned) >= n - 1, f"v={v} freed set spans only {len(spanned)} vertices"

        # Print first example
        if bits < 500:
            print(f"\nbits={bits}, bad={bad}")
            print(f"  #free={len(free_cycs)}, #unique_dom per vertex: {[(v,len(unique_dom_cycs[v])) for v in range(n)]}")
            print(f"  #multi_dom={len(multi_dom_cycs)}")

            for v in bad:
                freed = [cyc for cyc in free_cycs if v not in set(cyc)]
                freed += unique_dom_cycs.get(v, [])
                spanned = set()
                for cyc in freed:
                    spanned.update(cyc)
                print(f"  bad v={v}: #freed={len(freed)}, span={sorted(spanned)}")

    # PART 4: Can freed cycles for ALL 4 hypothetical bad vertices span n-1?
    # Key constraint: each cycle has 3 vertices.
    # Free cycles contribute to n-3 freed sets each.
    # Unique-dom cycles contribute to 1 freed set each.
    #
    # For B bad vertices each needing n-1 vertex coverage:
    # Total vertex-incidences needed >= B * (n-1)
    #
    # From free cycles: each contributes 3 to n-3 freed sets = 3*(n-3) incidences total
    # But only to freed sets of v not in cycle.
    # If v IS in the cycle, that cycle doesn't help v.
    #
    # From unique-dom cycles: each contributes 3 incidences to exactly 1 freed set.
    #
    # So: 3*(n-3)*|F| + 3*sum|D_v| >= B*(n-1)
    # where F = free cycles, D_v = unique-dom by v, sum over bad v only.
    #
    # But sum|D_v| over ALL v = total unique-dom cycles.
    # Only the bad vertices' D_v matter.
    #
    # At n=6: we need B*(5) = 5B coverage.
    # Free cycles: each gives 3*3=9 incidences spread over 3 freed sets.
    # Unique-dom: 3 each to 1 freed set.
    #
    # This doesn't immediately give a contradiction. Need tighter analysis.

    # PART 5: Bridge vertex deep analysis
    print("\n" + "=" * 70)
    print("PART 5: Bridge vertex — how does deletion create b1=1?")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)
    bridge_analysis = defaultdict(int)

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

        # Find transitive order
        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)

        # For EACH bad vertex: analyze the 3-cycle graph of T\v
        for label, v in [('top', top), ('mid', mid), ('bot', bot)]:
            # In T: which cycles go through v?
            thru_v = [cyc for cyc in cycles if v in set(cyc)]
            not_thru_v = [cyc for cyc in cycles if v not in set(cyc)]

            # Among not_thru_v: how many are free/dominated?
            free_not_v = [cyc for cyc in not_thru_v if not is_dominated(A, n, cyc)]
            dom_not_v = [cyc for cyc in not_thru_v if is_dominated(A, n, cyc)]

            # Among dom_not_v: how many have dom = {v}?
            newly_freed = [cyc for cyc in dom_not_v if get_dom(A, n, cyc) == [v]]

            # Among thru_v: were any of them the ONLY dominated cycle
            # connecting to a free subgraph?
            bridge_analysis[f"{label}_thru_v"] += len(thru_v)
            bridge_analysis[f"{label}_free_not_v"] += len(free_not_v)
            bridge_analysis[f"{label}_newly_freed"] += len(newly_freed)
            bridge_analysis[f"{label}_dom_not_v"] += len(dom_not_v)

    n_tours = sum(1 for k, v in bridge_analysis.items() if k == "top_thru_v")
    if n_tours > 0:
        print(f"Averages over {bridge_analysis.get('top_thru_v',0)//1 if n_tours==0 else n_tours} tournaments:")
        # Actually compute properly
        for label in ['top', 'mid', 'bot']:
            total_t = bridge_analysis.get(f'{label}_thru_v', 0)
            # Count is accumulated over all tours, divide by count of tours
            # We need to track count separately

    # Print raw totals
    print("Raw totals (divide by #tours for averages):")
    for k in sorted(bridge_analysis.keys()):
        print(f"  {k}: {bridge_analysis[k]}")

    # PART 6: For each bad vertex, track which component gains a free subcomponent
    print("\n" + "=" * 70)
    print("PART 6: Component structure change upon vertex deletion")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)
    comp_change = defaultdict(int)
    n_analyzed = 0

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
        n_analyzed += 1

        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)
        comps = cycle_graph_components(cycles)
        n_comps_T = len(comps)

        for label, v in [('top', top), ('mid', mid), ('bot', bot)]:
            B = delete_vertex(A, n, v)
            keep = [i for i in range(n) if i != v]

            # Cycles in T\v = cycles of T not through v (with relabeled vertices)
            # Plus possibly new cycles? No - T\v is a subtournament, cycles are subset.
            cyc_B = find_3cycles(B, n-1)
            comps_B = cycle_graph_components(cyc_B)
            n_comps_B = len(comps_B)

            # How many free components in T\v?
            free_comps_B = sum(1 for comp in comps_B
                              if all(not is_dominated(B, n-1, cyc_B[ci]) for ci in comp))

            comp_change[f"{label}_n_comps_T"] += n_comps_T
            comp_change[f"{label}_n_comps_Tv"] += n_comps_B
            comp_change[f"{label}_free_comps_Tv"] += free_comps_B

    print(f"Analyzed {n_analyzed} SC kappa>=2 tournaments with 3 bad vertices")
    for k in sorted(comp_change.keys()):
        avg = comp_change[k] / n_analyzed if n_analyzed > 0 else 0
        print(f"  {k}: total={comp_change[k]}, avg={avg:.2f}")

    # PART 7: The crucial test - can we find ANY tournament at ANY n with 4 bad vertices?
    print("\n" + "=" * 70)
    print("PART 7: Exhaustive/sampled search for 4+ bad vertices")
    print("-" * 50)

    rng = np.random.RandomState(999)
    for n in [5, 6, 7, 8, 9]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        max_bad = 0
        bad_dist = defaultdict(int)

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue

            bad_count = 0
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_count += 1

            bad_dist[bad_count] += 1
            max_bad = max(max_bad, bad_count)

        print(f"\n  n={n} ({method}): max bad = {max_bad}")
        for bc in sorted(bad_dist.keys()):
            print(f"    {bc} bad: {bad_dist[bc]}")

if __name__ == '__main__':
    main()
