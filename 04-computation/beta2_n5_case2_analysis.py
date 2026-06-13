"""
beta2_n5_case2_analysis.py — opus-2026-04-04-S1
Analyze the exact structure of n=5 Case 2 tournaments
(b1=0, SC, kappa>=2, free cycles exist) to find an algebraic proof
that a good vertex always exists.

The goal: close the ONLY remaining gap in the beta_2=0 proof (THM-108/109).
"""
import itertools
import numpy as np
from collections import defaultdict

def tournament_from_bits(n, bits):
    """Create adjacency matrix from bit encoding."""
    T = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def find_3cycles(T, n):
    """Find all directed 3-cycles."""
    cycles = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                if T[a][b] and T[b][c] and T[c][a]:
                    triple = tuple(sorted([a,b,c]))
                    if triple not in [tuple(sorted(x)) for x in cycles]:
                        cycles.append((a,b,c))
    # Deduplicate by vertex set
    seen = set()
    unique = []
    for c in cycles:
        key = tuple(sorted(c))
        if key not in seen:
            seen.add(key)
            unique.append(key)
    return unique

def is_dominated(T, n, cycle_verts):
    """Check if a 3-cycle is dominated (has external vertex that beats all or loses to all)."""
    a, b, c = cycle_verts
    external = [v for v in range(n) if v not in (a,b,c)]
    dominators = []
    for d in external:
        beats_all = T[d][a] and T[d][b] and T[d][c]
        loses_all = T[a][d] and T[b][d] and T[c][d]
        if beats_all or loses_all:
            dominators.append(d)
    return dominators

def share_directed_edge(T, c1, c2):
    """Check if two 3-cycles (given as sorted vertex tuples) share a directed edge in T."""
    edges1 = set()
    for i, a in enumerate(c1):
        for j, b in enumerate(c1):
            if a != b and T[a][b]:
                edges1.add((a,b))
    edges2 = set()
    for i, a in enumerate(c2):
        for j, b in enumerate(c2):
            if a != b and T[a][b]:
                edges2.add((a,b))
    return len(edges1 & edges2) > 0

def is_strongly_connected(T, n):
    """Check if tournament is strongly connected."""
    # BFS from 0
    reachable = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if u not in reachable and T[v][u]:
                reachable.add(u)
                queue.append(u)
    if len(reachable) < n:
        return False
    # BFS in reverse
    reachable = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if u not in reachable and T[u][v]:
                reachable.add(u)
                queue.append(u)
    return len(reachable) == n

def vertex_connectivity(T, n):
    """Compute vertex connectivity kappa(T)."""
    if not is_strongly_connected(T, n):
        return 0
    # Try removing each subset of size k
    for k in range(1, n):
        for subset in itertools.combinations(range(n), k):
            remaining = [v for v in range(n) if v not in subset]
            if len(remaining) <= 1:
                continue
            # Check if T restricted to remaining is SC
            n_rem = len(remaining)
            T_rem = np.zeros((n_rem, n_rem), dtype=int)
            for i, u in enumerate(remaining):
                for j, v in enumerate(remaining):
                    T_rem[i][j] = T[u][v]
            if not is_strongly_connected(T_rem, n_rem):
                return k
    return n - 1

def compute_b1(T, n):
    """Compute b1(T) = number of free components in cycle adjacency graph."""
    cycles = find_3cycles(T, n)
    if not cycles:
        return 0

    # Build adjacency of cycle graph (share directed edge)
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if share_directed_edge(T, cycles[i], cycles[j]):
                adj[i][j] = adj[j][i] = True

    # Find connected components
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        queue = [start]
        visited[start] = True
        while queue:
            v = queue.pop(0)
            comp.append(v)
            for u in range(nc):
                if not visited[u] and adj[v][u]:
                    visited[u] = True
                    queue.append(u)
        components.append(comp)

    # Count free components
    free_comps = 0
    for comp in components:
        all_free = all(len(is_dominated(T, n, cycles[i])) == 0 for i in comp)
        if all_free:
            free_comps += 1

    return free_comps

def analyze_n5():
    n = 5
    num_arcs = n*(n-1)//2  # 10

    case2_tours = []  # b1=0, SC, kappa>=2, has free cycles

    for bits in range(1 << num_arcs):
        T = tournament_from_bits(n, bits)

        b1 = compute_b1(T, n)
        if b1 != 0:
            continue

        if not is_strongly_connected(T, n):
            continue

        kappa = vertex_connectivity(T, n)
        if kappa < 2:
            continue

        cycles = find_3cycles(T, n)
        has_free = any(len(is_dominated(T, n, c)) == 0 for c in cycles)
        if not has_free:
            continue

        case2_tours.append((bits, T.copy(), cycles))

    print(f"n=5 Case 2 tournaments: {len(case2_tours)} / 1024")

    # For each Case 2 tournament, analyze the free-dom structure
    for idx, (bits, T, cycles) in enumerate(case2_tours[:20]):  # first 20 for inspection
        free_cycles = [c for c in cycles if len(is_dominated(T, n, c)) == 0]
        dom_cycles = [c for c in cycles if len(is_dominated(T, n, c)) > 0]

        # Find good vertices
        good = []
        bad = []
        for v in range(n):
            T_del = np.delete(np.delete(T, v, 0), v, 1)
            remaining = [u for u in range(n) if u != v]
            # Remap T_del
            n_del = n - 1
            T_v = np.zeros((n_del, n_del), dtype=int)
            for i, u in enumerate(remaining):
                for j, w in enumerate(remaining):
                    T_v[i][j] = T[u][w]
            b1_v = compute_b1(T_v, n_del)
            if b1_v <= 0:
                good.append(v)
            else:
                bad.append(v)

        scores = [sum(T[v]) for v in range(n)]
        c3_per_v = []
        for v in range(n):
            c3v = sum(1 for c in cycles if v in c)
            c3_per_v.append(c3v)

        if idx < 5 or bad:
            print(f"\n  Tour #{idx}: bits={bits}, scores={scores}, #cycles={len(cycles)}, "
                  f"free={len(free_cycles)}, dom={len(dom_cycles)}")
            print(f"    c3/vertex: {c3_per_v}")
            print(f"    good={good}, bad={bad}")
            print(f"    free cycles: {free_cycles}")
            print(f"    dom cycles: {[(c, is_dominated(T,n,c)) for c in dom_cycles]}")

    # Count bad vertex statistics across ALL Case 2 tournaments
    bad_counts = defaultdict(int)
    bad_score_patterns = defaultdict(int)
    total_bad = 0
    max_bad = 0

    for bits, T, cycles in case2_tours:
        free_cycles = [c for c in cycles if len(is_dominated(T, n, c)) == 0]

        bad_verts = []
        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            n_del = n - 1
            T_v = np.zeros((n_del, n_del), dtype=int)
            for i, u in enumerate(remaining):
                for j, w in enumerate(remaining):
                    T_v[i][j] = T[u][w]
            b1_v = compute_b1(T_v, n_del)
            if b1_v > 0:
                bad_verts.append(v)

        nb = len(bad_verts)
        bad_counts[nb] += 1
        total_bad += nb
        max_bad = max(max_bad, nb)

        if nb > 0:
            scores = [sum(T[v]) for v in range(n)]
            bad_scores = tuple(sorted([scores[v] for v in bad_verts]))
            bad_score_patterns[bad_scores] += 1

    print(f"\n\n=== STATISTICS ACROSS {len(case2_tours)} CASE 2 TOURNAMENTS ===")
    print(f"Bad vertex count distribution:")
    for k in sorted(bad_counts.keys()):
        print(f"  {k} bad: {bad_counts[k]} tournaments")
    print(f"Max bad vertices: {max_bad}")
    print(f"Total bad vertex occurrences: {total_bad}")

    if bad_score_patterns:
        print(f"\nBad vertex score patterns:")
        for pat, count in sorted(bad_score_patterns.items()):
            print(f"  scores {pat}: {count} tournaments")

    # KEY ANALYSIS: For the worst case (max bad), look at what makes good vertices good
    print("\n\n=== DETAILED ANALYSIS OF WORST-CASE TOURNAMENTS ===")
    for bits, T, cycles in case2_tours:
        free_cycles = [c for c in cycles if len(is_dominated(T, n, c)) == 0]
        dom_cycles = [(c, is_dominated(T, n, c)) for c in cycles if len(is_dominated(T, n, c)) > 0]

        bad_verts = []
        good_verts = []
        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            T_v = np.zeros((n-1, n-1), dtype=int)
            for i, u in enumerate(remaining):
                for j, w in enumerate(remaining):
                    T_v[i][j] = T[u][w]
            b1_v = compute_b1(T_v, n-1)
            if b1_v > 0:
                bad_verts.append(v)
            else:
                good_verts.append(v)

        if len(bad_verts) == max_bad and max_bad > 0:
            scores = [sum(T[v]) for v in range(n)]
            c3_per_v = [sum(1 for c in cycles if v in c) for v in range(n)]

            # For each vertex, analyze freed/remaining-dom structure
            print(f"\n  Tour bits={bits}, scores={scores}")
            print(f"  c3/vertex: {c3_per_v}")
            print(f"  bad={bad_verts}, good={good_verts}")
            print(f"  free: {free_cycles}")
            print(f"  dom: {dom_cycles}")

            for v in range(n):
                freed = [c for c in cycles if v not in c and
                         (len(is_dominated(T,n,c)) == 0 or is_dominated(T,n,c) == [v])]
                rem_dom = [c for c in cycles if v not in c and
                          any(d != v for d in is_dominated(T,n,c))]
                freed_verts = set()
                for c in freed:
                    freed_verts.update(c)
                freed_verts.discard(v)

                has_iso_edge = False
                for f in freed:
                    for r in rem_dom:
                        if share_directed_edge(T, f, r):
                            has_iso_edge = True
                            break
                    if has_iso_edge:
                        break

                tag = "BAD" if v in bad_verts else "GOOD"
                print(f"    v={v} [{tag}]: score={scores[v]}, c3={c3_per_v[v]}, "
                      f"#freed={len(freed)}, freed_span={len(freed_verts)}, "
                      f"#rem_dom={len(rem_dom)}, iso_edge={has_iso_edge}")

            # Only show first 3 worst-case examples
            if len([1 for b,t,c in case2_tours[:case2_tours.index((bits,T,cycles))+1]
                    if any(compute_b1(np.zeros((1,1)),1) >= 0 for _ in [1])]) > 3:
                break

    # FINAL: look for algebraic pattern
    print("\n\n=== PATTERN SEARCH ===")
    # Hypothesis: at n=5, every Case 2 tournament has at most 2 bad vertices
    # (so at least 3 good), and the good vertices are characterized by...

    # For each pair (free cycle F, adjacent dom cycle D), count how many vertices
    # are rescued by different structural properties
    rescue_mechanisms = defaultdict(int)
    for bits, T, cycles in case2_tours:
        free_cycles = [c for c in cycles if len(is_dominated(T, n, c)) == 0]
        dom_cycles_with_doms = [(c, is_dominated(T, n, c)) for c in cycles if len(is_dominated(T, n, c)) > 0]

        # Find a free-dom adjacent pair
        for fc in free_cycles:
            for dc, doms in dom_cycles_with_doms:
                if share_directed_edge(T, fc, dc):
                    shared_verts = set(fc) & set(dc)
                    all_verts = set(fc) | set(dc)
                    dominator_set = set(doms)
                    outside = set(range(n)) - all_verts - dominator_set

                    rescue_mechanisms[('shared', len(shared_verts))] += 1
                    rescue_mechanisms[('all_cycle_verts', len(all_verts))] += 1
                    rescue_mechanisms[('outside_verts', len(outside))] += 1
                    rescue_mechanisms[('doms_in_cycle', len(dominator_set & all_verts))] += 1
                    rescue_mechanisms[('multi_dom', len(doms) > 1)] += 1

    for key, count in sorted(rescue_mechanisms.items()):
        print(f"  {key}: {count}")

if __name__ == "__main__":
    analyze_n5()
